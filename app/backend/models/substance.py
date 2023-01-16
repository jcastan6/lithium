from app import db
from app.backend.models.rdkit_postgresql.types import Mol

class Substance(db.Model):
    id = db.Column('id', db.BigInteger, primary_key=True, nullable=False)    
    mol = db.Column('mol', Mol , nullable=False)

    
