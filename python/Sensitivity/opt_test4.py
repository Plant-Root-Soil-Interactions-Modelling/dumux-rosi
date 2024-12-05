import math
from hyperopt import fmin, tpe, hp
from hyperopt.mongoexp import MongoTrials

trials = MongoTrials('mongo://localhost:1234/foo_db/jobs', exp_key = 'exp1')

best = fmin(math.sin, hp.uniform('x', -2, 2), trials = trials, algo = tpe.suggest, max_evals = 100)

#
# import time
# from hyperopt import fmin, tpe, hp, Trials
# from hyperopt.mongoexp import MongoTrials
# import random
#
# # MongoDB connection string
# uri = 'mongo://localhost:27017/hyperopt'  # Adjust to your MongoDB URI
# exp_name = 'experiment'
#
#
# # Function to optimize (for example, minimizing a simple function)
# def objective(params):
#     x = params['x']
#     # Example objective function: a simple quadratic function
#     return {'loss': x ** 2, 'status': 'ok'}
#
#
# # Define the search space
# space = {
#     'x': hp.uniform('x', -10, 10)
# }
#
# # Set up MongoTrials for parallel optimization
# trials = MongoTrials(uri + '/' + exp_name, exp_key = exp_name)
#
# # Perform the optimization with parallel trials
# best = fmin(
#     fn = objective,  # The objective function to minimize
#     space = space,  # The search space for hyperparameters
#     algo = tpe.suggest,  # The algorithm to use (Tree-structured Parzen Estimator)
#     max_evals = 100,  # Number of evaluations
#     trials = trials,  # Store trials in MongoDB
#     verbose = 1  # Print progress to console
# )
#
# # Output the best result
# print("Best hyperparameters found: ", best)
