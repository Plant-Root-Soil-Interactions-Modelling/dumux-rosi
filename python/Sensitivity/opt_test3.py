# Import HyperOpt Library
from hyperopt import tpe, hp, fmin, Trials
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


def objective(params):
    x, y = params['x'], params['y']
    # print("*")
    return np.sin(np.sqrt(x ** 2 + y ** 2))

# x = np.linspace(-6, 6, 30)
# y = np.linspace(-6, 6, 30)
# x, y = np.meshgrid(x, y)
#
# z = objective({'x': x, 'y': y})
#
# fig = plt.figure()
# ax = plt.axes(projection = '3d')
# ax.plot_surface(x, y, z, cmap = cm.coolwarm)
#
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z')
#
# plt.show()


space = {
    'x': hp.uniform('x', -6, 6),
    'y': hp.uniform('y', -6, 6)
}

trials = Trials()

best = fmin(
    fn = objective,  # Objective Function to optimize
    space = space,  # Hyperparameter's Search Space
    algo = tpe.suggest,  # Optimization algorithm (representative TPE)
    max_evals = 1000,  # Number of optimization attempts
    trials = trials
)

print(best)
print(objective(best))
print()
print(trials)

