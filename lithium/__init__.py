from celery import Celery


celery = Celery(__name__, 
                broker='redis://localhost:6379/0',
                backend='pyamqp://guest@localhost//',
                include=['lithium.client.similarity.similarity']
                )
