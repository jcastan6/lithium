from lithium import celery
from lithium.app import init_app
from lithium.celery_worker import init_celery

app = init_app()


init_celery(app, celery)
