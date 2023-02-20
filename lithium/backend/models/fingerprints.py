from lithium.app import db
from lithium.backend.models.rdkit_postgresql.types import Bfp


class Fingerprints(db.Model):
    id = db.Column('id', db.BigInteger, primary_key=True, nullable=False)
    ecfp4 = db.Column('ecfp4', Bfp, nullable=False)
    xor_data = db.Column('xor_data', Bfp)
    num_bits = db.Column('num_bits', db.Integer)
    
    # set up many to many relationship with substance
    substance = db.relationship('Substance', secondary='fingerprints_substance', backref=db.backref('fingerprints', lazy='dynamic'))
