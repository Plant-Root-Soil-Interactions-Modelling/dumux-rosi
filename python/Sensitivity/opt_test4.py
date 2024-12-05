import time
from hyperopt import fmin, tpe, hp, Trials
from hyperopt.mongoexp import MongoTrials
import random

# MongoDB connection string
MONGO_URI = 'mongodb://localhost:27017/hyperopt'  # Adjust to your MongoDB URI
EXPERIMENT_NAME = 'my_experiment'


# Function to optimize (for example, minimizing a simple function)
def objective(params):
    x = params['x']
    # Example objective function: a simple quadratic function
    return {'loss': x ** 2, 'status': 'ok'}


# Define the search space
space = {
    'x': hp.uniform('x', -10, 10)
}

# Set up MongoTrials for parallel optimization
trials = MongoTrials(MONGO_URI + '/' + EXPERIMENT_NAME, exp_key = EXPERIMENT_NAME)

# Perform the optimization with parallel trials
best = fmin(
    fn = objective,  # The objective function to minimize
    space = space,  # The search space for hyperparameters
    algo = tpe.suggest,  # The algorithm to use (Tree-structured Parzen Estimator)
    max_evals = 100,  # Number of evaluations
    trials = trials,  # Store trials in MongoDB
    verbose = 1  # Print progress to console
)

# Output the best result
print("Best hyperparameters found: ", best)
