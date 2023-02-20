"""General page routes."""
from flask import Blueprint, request, make_response, send_file

from flask import render_template
from lithium.backend.models.substance import Substance
from lithium.backend.models.fingerprints import Fingerprints
from lithium.app import db
from rdkit.Chem import AllChem
from sqlalchemy import select, func
import io
import json
from sqlalchemy.sql.expression import cast
from matplotlib import colors
from rdkit.Chem import MolFromSmiles, rdFMCS, MolFromSmarts
from rdkit.Chem.Draw import MolToImage
from ctypes import ArgumentError
from lithium import celery

from sqlalchemy.sql import text
from celery import group
import psycopg2


# Blueprint Configuration
similarity_bp = Blueprint(
    "similarity_bp", __name__, template_folder="templates", static_folder="static"
)

@similarity_bp.route("/search/similarity/", methods=["POST", "GET"])
@similarity_bp.route("/search/similarity.<file_type>", methods=["POST"])
def similarity_search(file_type=None, smiles=None, page=None, asyncSearch=False, limit=1000):
    if not request.args.get('smiles') and not request.form.get('smiles'):
        return {"error":"No smiles submitted"}, 400
    if request.args.get('smiles'):
        smiles = request.args.get('smiles')
        
    if request.args.get('limit'):
        limit = request.args.get('limit')
        
    data = request.form

    if not smiles:
        smiles = str(data['smiles'])
        
    mol = MolFromSmiles(smiles)
    if not mol:
        return None

    mol = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=512)

    if not mol:
        return None

    res = group(get_similar_fingerprints.s(i, smiles, limit)
                for i in range(0, 3)).apply_async()
    
    res.save()
    
    if asyncSearch:
        return res.id
    else:
        data = res.get()
       
        res = []
        for i in data:
            res.extend(i)
        
        if limit:
            res = sorted(res, key=lambda x: x['similarity'], reverse=True)
            res = res[:int(limit)]

        return json.dumps(res)


@celery.task
def get_similar_fingerprints(i, smiles, limit=None):
    db.session.execute('set rdkit.tanimoto_threshold=0.4')

    mol = MolFromSmiles(smiles)
    mol = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=512)

    get_similar = Fingerprints.query.filter(Fingerprints.ecfp4.tanimoto_sml(cast(mol, Fingerprints.ecfp4.type)))\
        .filter(Fingerprints.id == Substance.id)\
        .with_entities(
            Fingerprints.id,
            func.tanimoto_sml(Fingerprints.ecfp4, cast(mol, Fingerprints.ecfp4.type)).label('similarity'), 
            func.mol_to_smiles(Substance.mol).label('smiles')
        ).order_by(func.tanimoto_sml(Fingerprints.ecfp4, cast(mol, Fingerprints.ecfp4.type)).desc())\
        .filter(Fingerprints.id >= i*10000000).filter(Fingerprints.id < (i+1)*10000000)    # use sqlalchemy to query postgres
        
    if limit:
        res = list(get_similar.limit(limit).all())
    else:
        res = list(get_similar.all())
    res = [{'id': x.id, 'similarity': x.similarity, 'smiles': x.smiles}
           for x in res]

    return list(res)
