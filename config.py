
"""Class-based Flask app configuration."""
from os import environ, path

from dotenv import load_dotenv

basedir = path.abspath(path.dirname(__file__))
load_dotenv(path.join(basedir, ".env"))


class Config:
    """Configuration from environment variables."""

    SECRET_KEY = environ.get("SECRET_KEY")
    FLASK_ENV = environ.get("FLASK_DEBUG")
    FLASK_APP = "lithium.py"

    # Flask-Assets
    LESS_BIN = environ.get("LESS_BIN")
    ASSETS_DEBUG = True
    LESS_RUN_IN_DEBUG = True

    # Static Assets
    STATIC_FOLDER = "static"
    TEMPLATES_FOLDER = "templates"
    COMPRESSOR_DEBUG = True
    
    CELERY_RESULT_BACKEND = environ.get("CELERY_RESULT_BACKEND")
    CELERY_BROKER_URL = environ.get("CELERY_BROKER_URL")
    SQLALCHEMY_DATABASE_URI = environ.get("SQLALCHEMY_DATABASE_URI")
    CELERY_IMPORTS = ("app.client.similarity.similarity")
    CELERY_ACCEPT_CONTENT = ['json']
    CELERY_RESULT_SERIALIZER = 'json'
    CELERY_TASK_SERIALIZER = 'json'