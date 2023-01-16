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
from sqlalchemy.sql.expression import cast
from matplotlib import colors
from rdkit.Chem import MolFromSmiles, rdFMCS, MolFromSmarts
from rdkit.Chem.Draw import MolToImage
# Blueprint Configuration
similarity_bp = Blueprint(
    "similarity_bp", __name__, template_folder="templates", static_folder="static"
)

@similarity_bp.route("/similarity", methods=["GET"])
def similarity(message = None):    
    page = request.args.get('page', 1, type=int)
    if request.args.get('smiles'):
        try:
            res = similarity_search(smiles=request.args.get('smiles'))
        except:
            res = []
            message = "Error: invalid smiles submitted: %s" % request.args.get('smiles')
        return render_template(
            "similarity.jinja2",
            title="Similarity Search",
            results=res,
            message = message,
            smiles = request.args.get('smiles'),
        )
    else:
        return render_template(
            "similarity.jinja2",
            title="Similarity Search",
            template="home-template",
        )
    
    

@similarity_bp.route("/similarity/search", methods=["POST"])    
@similarity_bp.route("/similarity/search.<file_type>", methods=["POST"])    
def similarity_search(file_type=None, smiles=None, page=None):
    data = request.form 
    
    if not smiles: 
        smiles = str(data['smiles'])
        
    mol = MolFromSmiles(smiles)
    
    mol= AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=512)
    
    if not mol:
        return "Error: invalid smiles submitted: %s" % smiles
    
    db.session.execute('set rdkit.tanimoto_threshold=0.4')
    
    get_similar = Fingerprints.query.filter(Fingerprints.ecfp4.tanimoto_sml(cast(mol, Fingerprints.ecfp4.type)))\
        .filter(Fingerprints.id == Substance.id)\
        .with_entities( 
            Fingerprints.id, \
            func.tanimoto_sml(Fingerprints.ecfp4, cast(mol, Fingerprints.ecfp4.type)).label('similarity')\
            , func.mol_to_smiles(Substance.mol).label('smiles')\
        ).order_by(func.tanimoto_sml(Fingerprints.ecfp4, cast(mol, Fingerprints.ecfp4.type)).desc())
        
    if page:
        get_similar = get_similar.paginate(page=page, per_page=100)
        res = get_similar.items
    else:
        res = get_similar.all()

   
    if file_type == 'json':
        return json.dumps([{'id': x.id, 'similarity': x.similarity, 'smiles':x.smiles} for x in res])
    else:
        return res
    

    
