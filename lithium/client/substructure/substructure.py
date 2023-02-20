"""General page routes."""
from flask import Blueprint, request, make_response, send_file
from flask import current_app as app
from flask import render_template
from lithium.backend.models.substance import Substance
from lithium.backend.models.fingerprints import Fingerprints

from rdkit.Chem import AllChem
from sqlalchemy import select, func
import io
import json
import pandas as pd
from matplotlib import colors
from rdkit.Chem import MolFromSmiles, rdFMCS, MolFromSmarts
from rdkit.Chem.Draw import MolToImage
from lithium import celery
from celery import group

# Blueprint Configuration
substructure_bp = Blueprint(
    "substructure_bp", __name__, template_folder="templates", static_folder="static"
)

@substructure_bp.route("/search/substructure", methods=["POST", "GET"])
@substructure_bp.route("/search/substructure.<file_type>", methods=["POST"])
def substructure_search(file_type=None, smiles=None, page=None, asyncSearch=False, limit=1000):
    if not request.args.get('smiles') and not request.form.get('smiles'):
        return {"error":"No smiles submitted"}, 400
    
    if request.args.get('smiles'):
        smiles = request.args.get('smiles')
    
    if request.args.get('limit'):
        limit = request.args.get('limit')
        limit = int(limit/3)
        
    data = request.form

    if not smiles:
        smiles = str(data['smiles'])

    mol = MolFromSmiles(smiles)
    if not mol:
        return None
    
    res = group(get_substructures.s(x,smiles, limit)
                for x in range(0, 3)).apply_async()
    res.save()
    
    if asyncSearch:
        return res.id
    
    data = res.get()
    
    res = []
    for x in data:
        res.extend(x)

    # if page:
    #     res = get_subtructures.paginate(page=page, per_page=100).items
    # else:
    #     res = get_subtructures.all()

    
    return json.dumps(res)
    


@celery.task
def get_substructures(i,smiles, limit=None):
    mol = MolFromSmiles(smiles)

    res = Substance.query.filter(Substance.mol.hassubstruct(mol))\
        .with_entities(Substance.id, func.mol_to_smiles(Substance.mol).label('smiles')).filter(Substance.id >= i*10000000).filter(Substance.id < (i+1)*10000000)
        
    
    if limit:
        res = res.limit(limit)
            
    res = list(res.all())

    res = [{'id': x[0], 'smiles': x[1]} for x in res]
    return list(res)
