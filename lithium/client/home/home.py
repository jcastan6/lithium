"""General page routes."""
from flask import Blueprint, request, make_response, send_file
from flask import current_app as app
from flask import render_template

from lithium import celery

# global definitions 

home_bp = Blueprint(
    "home_bp", __name__,
)


@home_bp.route("/asyncResult/<id>", methods=["GET"])
def get_result(id):
    res = celery.GroupResult.restore(id)

    results = []
    
    for i in res.children:
        
        if i.ready():
            results.extend(i.get())
    return results
        