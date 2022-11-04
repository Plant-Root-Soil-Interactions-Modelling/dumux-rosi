import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm


# defining a model
def model(x, a, b):
    return a * np.exp(-b * x)


# defining the x vector and the real value of some parameters
x_vector = np.arange(100)
a_real, b_real = 1, 0.05

# some toy data with multiplicative uncertainty
y_vector = model(x_vector, a_real, b_real) * (1 + norm.rvs(scale = 0.08, size = 100))

# fit the parameters, equal weighting on all data points
params, cov = curve_fit(model, x_vector, y_vector)
print(params)
print(cov)

# fit the parameters, weighting each data point by its inverse value
params, cov = curve_fit(model, x_vector, y_vector,
                        sigma = 1 / y_vector, absolute_sigma = False)
print(params)
print(cov)

# with absolute_sigma=False:
# # multiplicative transformations of y_data don't matter
params, cov = curve_fit(model, x_vector, y_vector,
                        sigma = 100 / y_vector, absolute_sigma = False)
print(params)
print(cov)

# but absolute_sigma=True:
# # multiplicative transformations of sigma carry through to pcov
params, cov = curve_fit(model, x_vector, y_vector,
                        sigma = 100 / y_vector, absolute_sigma = True)
print(params)
print(cov)

