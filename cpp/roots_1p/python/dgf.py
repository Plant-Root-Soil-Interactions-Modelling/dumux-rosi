import numpy as np

'''
Creates dgf files for the root system benchmarks 1 (and old benchmark 2 and 3), and C11
'''


def createDGF_1Droots(filename, nodes, seg, params = np.zeros((0, 0))):
    file = open(filename, "w")  # write file

    nop = params.shape[0]  # number of parameters
    file.write("DGF\n")
    file.write('Vertex\n')
    for i in range(0, len(nodes)):
        file.write('{:g} {:g} {:g} \n'.format(nodes[i, 0], nodes[i, 1], nodes[i, 2]))

    file.write('#\n');
    file.write('Simplex\n');
    if nop > 0:
        file.write('parameters {:d}\n'.format(nop))
    for i in range(0, len(seg)):
        file.write('{:g} {:g}'.format(seg[i, 0], seg[i, 1]))
        for j in range(0, nop):
            file.write(' {:g}'.format(params[j, i]))
        file.write(' \n')

    # not used...
    file.write('#\nBOUNDARYSEGMENTS\n2 0\n')
    file.write('3 {:g}\n'.format(len(seg)))
    file.write('#\nBOUNDARYDOMAIN\ndefault 1\n#\n')

    file.close()


if __name__ == "__main__":

    # create grid for benchmark 1
    L = 0.5  # length of single straight root (m)
    nnz = 100
    nodes = np.zeros((nnz, 3))
    seg = np.zeros(((nnz - 1), 2), dtype = int)
    for i in range(1, nnz):
        seg[i - 1, 0] = i - 1
        seg[i - 1, 1] = i
        nodes[i, :] = [0., 0., -i * L / (nnz - 1)]
#  0 order, 1 brnID, 2 surf [cm2], 3 length [cm], 4 radius [cm],
#  5 kz [cm4 hPa-1 d-1],  6 kr [cm hPa-1 d-1],  7 emergence time [d],
#  8 subType, 9 organType
    z_ = np.zeros((1, nnz - 1))
    o_ = np.ones((1, nnz - 1))
    a_ = o_ * 0.2;  # 0.2 cm
    ct_ = np.linspace(0, 25, nnz - 1)  # days (assuming 2 cm /day)
    kr_ = 0.;  # cm/hPa/day, = 2.e-9 m/Pa/s;
    kz_ = 0.;  # cm^4/hPa/day, = 5.e-13 m/Pa/s;
    dx_ = o_ * (L / (nnz - 1))
    params = np.vstack((z_, o_, 2 * np.pi * np.multiply(dx_, a_), dx_, a_, z_, z_, ct_, o_, 2 * o_))
    createDGF_1Droots("../../../grids/singleroot.dgf", nodes, seg, params)
    nodes2 = np.transpose(np.vstack((nodes[:, 2], nodes[:, 0], nodes[:, 1] - 1.e-4)))  # 1e.-4 to make it easier with the BC
    nodes2[0, 2] = 0  # collar BC
    createDGF_1Droots("../../../grids/singlerootH.dgf", nodes2, seg, params)

    # create grid for benchmark 2
    nodes = [ [0.00, 0.00, -2.9], [0.00, 0.00, -3.00], [-0.00, -0.01, -3.48], [-0.85, 0.48, -3.71], [-1.69, 0.99, -3.90], [-2.58, 1.32, -4.21], [-3.48, 1.67, -4.49], [-4.38, 2.00, -4.77], [-5.24, 2.40, -5.09], [-6.08, 2.82, -5.42], [-6.93, 3.27, -5.69], [-6.96, 3.29, -5.70], [-0.00, 0.01, -3.97], [0.20, -0.95, -4.20], [0.43, -1.88, -4.49], [0.65, -2.81, -4.77], [0.84, -3.75, -5.06], [1.04, -4.70, -5.31], [1.27, -5.64, -5.54], [1.43, -6.58, -5.84], [1.48, -6.91, -5.94], [-0.01, 0.03, -4.45], [0.75, 0.68, -4.48], [1.52, 1.32, -4.50], [2.30, 1.94, -4.46], [3.07, 2.58, -4.41], [3.88, 3.16, -4.46], [4.73, 3.69, -4.50], [5.34, 4.05, -4.53], [-0.03, 0.06, -4.97], [-0.73, 0.63, -5.40], [-1.46, 1.20, -5.79], [-2.15, 1.80, -6.18], [-2.76, 2.48, -6.59], [-3.16, 3.16, -7.21], [-3.63, 3.90, -7.64], [-0.06, 0.07, -5.42], [0.07, 1.04, -5.61], [0.23, 2.00, -5.82], [0.46, 2.93, -6.11], [0.72, 3.85, -6.40], [0.99, 4.77, -6.69], [1.06, 5.11, -6.78], [-0.10, 0.08, -5.96], [0.35, 0.90, -6.31], [0.80, 1.72, -6.65], [1.23, 2.55, -7.01], [1.62, 3.41, -7.34], [1.70, 3.58, -7.41], [-0.14, 0.10, -6.46], [0.05, -0.81, -6.83], [0.29, -1.72, -7.18], [0.43, -2.62, -7.59], [0.50, -3.12, -7.84], [-0.19, 0.12, -7.02], [0.32, 0.98, -7.03], [0.82, 1.84, -6.95], [1.07, 2.23, -6.88], [-0.23, 0.17, -7.56], [0.68, -0.18, -7.80], [1.12, -0.35, -7.93], [-0.24, 0.21, -8.14], [-0.01, 0.36, -8.22], [-0.25, 0.24, -8.69], [-0.24, 0.25, -9.25], [-0.24, 0.26, -9.71], [-0.26, 0.26, -10.09], [-0.28, 0.25, -10.57], [-0.26, 0.24, -11.05], [-0.26, 0.21, -11.58], [-0.25, 0.19, -12.06], [-0.25, 0.17, -12.55], [-0.23, 0.15, -13.00], [-0.23, 0.12, -13.46], [-0.22, 0.12, -13.99], [-0.21, 0.15, -14.54], [-0.19, 0.20, -15.07], [-0.17, 0.24, -15.60], [-0.17, 0.31, -16.17], [-0.15, 0.36, -16.64], [-0.12, 0.43, -17.18], [-0.11, 0.48, -17.70], [-0.07, 0.52, -18.23], [-0.06, 0.53, -18.53] ]
    seg = [ [0, 1], [1, 2], [2, 12], [12, 21], [21, 29], [29, 36], [36, 43], [43, 49], [49, 54], [54, 58], [58, 61], [61, 63], [63, 64], [64, 65], [65, 66], [66, 67], [67, 68], [68, 69], [69, 70], [70, 71], [71, 72], [72, 73], [73, 74], [74, 75], [75, 76], [76, 77], [77, 78], [78, 79], [79, 80], [80, 81], [81, 82], [82, 83], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10], [10, 11], [12, 13], [13, 14], [14, 15], [15, 16], [16, 17], [17, 18], [18, 19], [19, 20], [21, 22], [22, 23], [23, 24], [24, 25], [25, 26], [26, 27], [27, 28], [29, 30], [30, 31], [31, 32], [32, 33], [33, 34], [34, 35], [36, 37], [37, 38], [38, 39], [39, 40], [40, 41], [41, 42], [43, 44], [44, 45], [45, 46], [46, 47], [47, 48], [49, 50], [50, 51], [51, 52], [52, 53], [54, 55], [55, 56], [56, 57], [58, 59], [59, 60], [61, 62] ]
    age = [ 8, 7.76, 7.52, 7.29, 7.03, 6.80, 6.53, 6.28, 6.00, 5.73, 5.43, 5.16, 4.87, 4.64, 4.44, 4.20, 3.95, 3.68, 3.43, 3.17, 2.94, 2.70, 2.42, 2.14, 1.85, 1.57, 1.27, 1.02, 0.73, 0.45, 0.16, -0.00, 2.29, 2.02, 1.74, 1.43, 1.11, 0.77, 0.41, 0.01, 0.00, 2.04, 1.77, 1.48, 1.18, 0.85, 0.51, 0.14, 0.00, 1.78, 1.51, 1.22, 0.92, 0.60, 0.26, 0.00, 1.51, 1.24, 0.96, 0.66, 0.34, 0.00, 1.28, 1.01, 0.73, 0.43, 0.12, 0.00, 0.97, 0.69, 0.39, 0.07, 0.00, 0.72, 0.45, 0.17, 0.00, 0.41, 0.14, 0.00, 0.13, 0.00, 0.00 ]
    types = [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 ]
    a_tap = 0.2  # tap root radius (cm)
    a_lateral = 0.1  # lateral root radius (cm)
    kz0, kz1 = 4.32e-5, 8.64e-4  # cm^4 / hPa / day
    kr0, kr1 = 1.4688e-4, 1.728e-5  # cm / hPa / day
    kz = lambda age: kz1  # kz0 * (age <= 3) + kz1 * (age > 3)
    kr = lambda age: kr1  #  * (age <= 3) + kr1 * (age > 3)
    a = lambda t: (a_tap * (t == 1) + a_lateral * (t == 2))
    a_ = list(map(a, types))
    kr_ = list(map(kr, age))  # m / (Pa s)
    kz_ = list(map(kz, age))  # m^4 / (Pa s)

    nodes = np.array(nodes) * 1.e-2  # cm->m, and convert from list to numpy array
    seg = np.array(seg)  # convert from list to numpy array
    order = np.array(types) - 1;
    params = np.vstack((order, a_, age, kr_, kz_))
    createDGF_1Droots("../../../grids/rootsystem_b2.dgf", nodes, seg, params)

    kz = lambda age: kz0 * (age <= 3) + kz1 * (age > 3)
    kr = lambda age: kr1 * (age <= 3) + kr1 * (age > 3)
    kr_ = list(map(kr, age))  # m / (Pa s)
    kz_ = list(map(kz, age))  # m^4 / (Pa s)
    params = np.vstack((order, a_, age, kr_, kz_))
    createDGF_1Droots("../../../grids/rootsystem_b3.dgf", nodes, seg, params)

    # create grid for benchmark 1
    L = 0.01  #  1 cm
    nnz = 3
    nodes = np.zeros((nnz, 3))
    seg = np.zeros(((nnz - 1), 2), dtype = int)
    for i in range(1, nnz):
        seg[i - 1, 0] = i - 1
        seg[i - 1, 1] = i
        nodes[i, :] = [0., 0., -i * L / (nnz - 1)]
    order = np.zeros((1, nnz))
    a_ = np.ones((1, nnz)) * 0.02;  # 0.2 cm
    age = np.zeros((1, nnz))
    kr_ = np.ones((1, nnz)) * 1.;  # cm/hPa/day
    kz_ = np.ones((1, nnz)) * 1.;  # cm^4/hPa/day
    params = np.vstack((order, a_, age, kr_, kz_))
    createDGF_1Droots("../../../grids/singleC11.dgf", nodes, seg, params)

    print("its done.")

