# from celery import group
# from celery import Task
# from lithium import celery


# def parallelize_query(queries, per_query_limit=0, total_limit=0):
#     results = []

#     # if there is a query limit, run queries in parallel and yield results until limit is reached
#     print(queries)
#     group_results = group(run_query.s(query, per_query_limit)
#                           for query in queries)()
#     print(group_results)
#     return group_results


# @celery.task
# def run_query(query, per_query_limit=0):
#     # if there is a limit, yield results until limit is reached
#     return query.all()
