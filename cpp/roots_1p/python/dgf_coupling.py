import numpy as np

'''
Creates a single root dgf files for dumux-rootgrowth module
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

    # create single root for dumux-rootgrowth
    L = 0.2  # length of single straight root (m)
    nnz = 100  # 200 causes problems..., 50 maybe too...
    nodes = np.zeros((nnz, 3))
    seg = np.zeros(((nnz - 1), 2), dtype = int)
    for i in range(1, nnz):
        seg[i - 1, 0] = i - 1
        seg[i - 1, 1] = i
        nodes[i, :] = [0., 0., -i * L / (nnz - 1)]
    order = np.zeros((1, nnz))
    a_ = np.ones((1, nnz)) * 0.1;  # 0.2 cm
    age = np.zeros((1, nnz))
    kr_ = np.ones((1, nnz)) * 1.8e-5  # cm hPa-1 d-1 #
    kz_ = np.ones((1, nnz)) * 4.3;  #  cm4 hPa-1 d-1  # not used by dumux root-growth
    params = np.vstack((order, order, age, age, a_, age, age, age))  # ONLY put 0 and radius (a_)
    createDGF_1Droots("../../../grids/singlerootC.dgf", nodes, seg, params)

    print("its done.")

