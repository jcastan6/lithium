from lithium.app import db

class catalog_substance(db.Model):
    cat_content_fk = db.relationship('catalog_content', backref=db.backref('catalog_substance', lazy=True))
    sub_id_fk = db.relationship('substance', backref=db.backref('catalog_substance', lazy=True))
    cat_sub_id = db.Column('cat_sub_itm_id', db.BigInteger, primary_key=True, nullable=False)
    
    
    

    
    