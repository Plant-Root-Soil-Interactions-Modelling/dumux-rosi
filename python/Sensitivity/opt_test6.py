from bayes_opt import BayesianOptimization
from bayes_opt.logger import JSONLogger
from bayes_opt.event import Events


def black_box_function(x, y):
    print("*******************************************")
    return -x ** 2 - (y - 1) ** 2 + 1


def example1():

    pbounds = {'x': (2, 4), 'y': (-3, 3)}

    optimizer = BayesianOptimization(
        f = black_box_function,
        pbounds = pbounds,
        verbose = 2,  # verbose = 1 prints only when a maximum is observed, verbose = 0 is silent
        random_state = 1,
    )

    logger = JSONLogger(path = "./logs.log")
    optimizer.subscribe(Events.OPTIMIZATION_STEP, logger)

    optimizer.maximize(
        init_points = 2,
        n_iter = 3,
    )

    print(optimizer.max)


def example2():

    pbounds = {'x': (2, 4), 'y': (-3, 3)}

    optimizer = BayesianOptimization(
        f = black_box_function,
        pbounds = pbounds,
        verbose = 2,  # verbose = 1 prints only when a maximum is observed, verbose = 0 is silent
        random_state = 1,
    )

    next_point_to_probe = optimizer.suggest()
    print("Next point to probe is:", next_point_to_probe)

    print(optimizer.suggest())
    print(optimizer.suggest())
    print(optimizer.suggest())
    print(optimizer.suggest())


if __name__ == "__main__":
    example2()
