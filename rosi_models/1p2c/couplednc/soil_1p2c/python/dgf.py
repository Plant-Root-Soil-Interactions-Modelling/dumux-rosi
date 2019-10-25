import numpy as np


def createDGF_1D(filename, N, depth, top, bot, domainId, vertex3d = False):

    z_ = np.linspace(0, -depth, N)
    initial = np.linspace(top, bot, N)  # per node
    initialC = np.linspace(top, bot, N - 1)  # per cell
    id = range(0, N)

    file = open(filename, "w")

    file.write("DGF\nVertex\n")
    file.write('parameters 2\n')  # initial data, domain index
    for i in range(0, N):
        if vertex3d:
            file.write('{:g} {:g} {:g} {:g} {:g}\n'.format(np.zeros(N,), np.zeros(N,), z_[i], initialC[i], domainId[i]))
        else:
            file.write('{:g} {:g} {:g}\n'.format(z_[i], initial[i], domainId[i]))

    file.write('#\nSimplex\n')
    file.write('parameters 2\n')  # initial data, domain index
    for i in range(0, N - 1):
        file.write('{:g} {:g} {:g} {:g}\n'.format(id[i], id[i + 1], initialC[i], domainId[i]));

    file.write('#\nBOUNDARYSEGMENTS\n')  # not used...
    file.write('2 0\n')
    file.write('3 {:g}\n'.format(N - 1))  # vertex id, index starts with 0
    file.write('#\nBOUNDARYDOMAIN\ndefault 1\n#\n')
    file.close()


if __name__ == "__main__":

    # Benchmark 1, for Figure 2abc
    domain_b1 = np.hstack((np.ones(50,), 2 * np.ones(151,)))
    createDGF_1D("../grids/b1.dgf", 201, 2., -200, -200., domain_b1)

    # Benchmark 2, for Figure 3
    domain_b2 = np.ones(55,)
    createDGF_1D("../grids/b2.dgf", 55, .54, -54, 0, domain_b2)

    # Benchmark 3, for Figure 4abc
    domain_b3 = np.ones(201,)
    createDGF_1D("../grids/b3.dgf", 201, 2., -400, -400., domain_b3)

    # Benchmark 4, for Figure 5abcd
    createDGF_1D("../grids/b4.dgf", 201, 2., -200, -200., domain_b3)
    createDGF_1D("../grids/b4sand.dgf", 201, 2., -40, -40., domain_b3)

    # Benchmark 4, for Figure 5abcd
    domain_b4hr = np.ones(801,)
    createDGF_1D("../grids/b4hr.dgf", 801, 2., -200, -200., domain_b4hr)
    createDGF_1D("../grids/b4hr_sand.dgf", 801, 2., -40, -40., domain_b4hr)

    print("its done.")
