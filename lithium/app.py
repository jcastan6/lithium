from flask import Flask
from flask_assets import Environment
from flask_bootstrap  import Bootstrap5
from flask_sqlalchemy import SQLAlchemy
from config import Config
import urllib
from sqlalchemy.pool import NullPool
from flask_cors import CORS

db = SQLAlchemy()

app = Flask(__name__, instance_relative_config=False)


def init_app():
    app.config.from_object(Config)
    
    assets = Environment()
    assets.init_app(app)
    
    bootstrap = Bootstrap5(app)
    
    
    with app.app_context():
        # Import parts of our application
    
       
        from .client.home import home
        from .client.similarity import similarity
        from .client.substructure import substructure
        from .client.substance import substance
        # Register Blueprints
        app.register_blueprint(home.home_bp)
        app.register_blueprint(similarity.similarity_bp)
        app.register_blueprint(substructure.substructure_bp)
        app.register_blueprint(substance.substance_bp)

        app.config["SQLALCHEMY_DATABASE_URI"] = "postgresql:///molecules"
        app.config['SQLALCHEMY_ENGINE_OPTIONS'] = {
           'poolclass': NullPool,
        }
        CORS(app)
        db.init_app(app)
        db.create_all()
        app.jinja_env.filters['parse_url'] = lambda u: urllib.parse.quote(u)
        return app