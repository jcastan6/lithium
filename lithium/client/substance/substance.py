"""General page routes."""
from flask import Blueprint, request, make_response, send_file
from flask import current_app as app

from lithium.backend.models.substance import Substance
from sqlalchemy import select, func
import io
import json
import pandas as pd
from lithium import celery


# Blueprint Configuration
substance_bp = Blueprint(
    "substance_bp", __name__,
)

@substance_bp.route("/substance", methods=["POST", "GET"])
@substance_bp.route("/search/substructure.<file_type>", methods=["POST"])
def substance(file_type=None, smiles=None, page=None, per_page=50, limit=1000):

    if request.args.get('page'):
        page = int(request.args.get('page'))
    
    if request.args.get('per_page'):
        per_page = int(request.args.get('per_page'))
        
    res = get_substances.s(page, per_page).apply_async() 
    res = res.get()
    print(res)
    return json.dumps(res)

@substance_bp.route("/substance/total", methods=["POST", "GET"])
def substance_total(file_type=None, smiles=None, page=None, per_page=50, limit=1000):
    res = get_total.s().apply_async()
    res = res.get()
    return json.dumps(res)

@celery.task
def get_substances(page, per_page):
    if page:
        res = Substance.query.with_entities(Substance.id, func.mol_to_smiles(Substance.mol).label('smiles'))\
        .order_by(Substance.id)\
        .paginate(page=page, per_page=per_page).items
    else:
        res = []
    
    res = [{"id":x[0], "smiles":x[1]} for x in res]
  
    print(res)
    return res
    
    
@celery.task
def get_total():
    #get total number of substances
    res = Substance.query.with_entities(func.count(Substance.id)).first()
    return res[0]
    