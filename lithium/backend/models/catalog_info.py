from lithium.app import db

class catalog_info(db.Model):
    cat_id = db.Column('cat_id', db.BigInteger, nullable=False)

    
    