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

    createDGF_1Droots("../grids/singleroot.dgf", nodes, seg, params)

    print("its done.")

