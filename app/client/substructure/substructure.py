"""General page routes."""
from flask import Blueprint, request, make_response, send_file
from flask import current_app as app
from flask import render_template
from app.backend.models.substance import Substance
from app.backend.models.fingerprints import Fingerprints
from app import db
from rdkit.Chem import AllChem
from sqlalchemy import select, func
import io, json
import pandas as pd
from matplotlib import colors
from rdkit.Chem import MolFromSmiles, rdFMCS, MolFromSmarts
from rdkit.Chem.Draw import MolToImage
# Blueprint Configuration
substructure_bp = Blueprint(
    "substructure_bp", __name__, template_folder="templates", static_folder="static"
)

@substructure_bp.route("/substructure", methods=["GET"])
def substructure(message = None):
    page = request.args.get('page', 1, type=int)
    
    if request.args.get('smiles'):
        res = substructure_search(smiles=request.args.get('smiles'))
        if not res:
            res = []
            message = "Error: invalid smiles submitted: %s" % request.args.get('smiles')        
        
        return render_template(
            "substructure.jinja2",
            title="Substructure Search",
            results=res,
            message = message,
            smiles = request.args.get('smiles'),
        )
    else:
        return render_template(
            "substructure.jinja2",
            title="Substructure Search",
            template="home-template",
        )
    
@substructure_bp.route("/substructure/search", methods=["POST"])    
@substructure_bp.route("/substructure/search.<file_type>", methods=["POST"])    
def substructure_search(file_type=None, smiles=None, page=None):
    data = request.form 
    
    if not smiles: 
        smiles = str(data['smiles'])
        
    mol = MolFromSmiles(smiles)
    if not mol:
        return None
    
    
    get_subtructures = Substance.query.filter(Substance.mol.hassubstruct(mol))\
            .with_entities(Substance.id, func.mol_to_smiles(Substance.mol).label('smiles')).limit(100)
    
    if page:
        res = get_subtructures.paginate(page=page, per_page=100).items
    else:
        res = get_subtructures.all()

    if file_type == 'json':
        return json.dumps([{'id': x.id, 'smiles': x.smiles} for x in res])
    else:
        return res
    
