from lithium.app import db
from lithium.backend.models.rdkit_postgresql.types import Mol
from sqlalchemy import Index


class Substance(db.Model):
    id = db.Column('id', db.BigInteger, primary_key=True, nullable=False)
    mol = db.Column('mol', Mol, nullable=False)

    def as_dict(self):
       return {c.name: getattr(self, c.name) for c in self.__table__.columns}