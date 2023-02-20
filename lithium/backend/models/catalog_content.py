from lithium.app import db

class catalog_content(db.Model):
    cat_content_id = db.Column('cat_content_id', db.BigInteger, primary_key=True, nullable=False)
    cat_id_fk = db.relationship('catalog_into', backref=db.backref('catalog_content', lazy=True), nullable=False)
    supplier_code = db.Column('supplier_code', db.String(255), nullable=False)
    depleted = db.Column('depleted', db.Boolean, nullable=False)
    pack_size_mg = db.Column('pack_size_mg', db.Integer, nullable=False)
    price_usd = db.Column('price_usd', db.Float, nullable=False)
    price_eur = db.Column('price_eur', db.Float, nullable=False)
    
    def as_dict(self):
         return {c.name: getattr(self, c.name) for c in self.__table__.columns}

    
    
