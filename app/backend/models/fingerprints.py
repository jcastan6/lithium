from app import db
from app.backend.models.rdkit_postgresql.types import Bfp

class Fingerprints(db.Model):
    id = db.Column('id', db.BigInteger, primary_key=True, nullable=False)
    ecfp4 = db.Column('ecfp4', Bfp, nullable=False)
    
# create 3 partitions for the table, for min and max values of the id column
# db.session.execute('create range substance_id_range (bigint, bigint) \
#     partition by range (id) (start (1) end (1000000) every (1000000), \
#     start (1000000) end (2000000) every (1000000), \
#     start (2000000) end (3000000) every (1000000))')
