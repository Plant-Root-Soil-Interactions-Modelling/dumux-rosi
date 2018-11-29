import numpy as np


def createDGF_1D(filename, N, depth, top, bot, domainId, vertex3d = False):

    # prepare arrays
    z_ = np.linspace(0, -depth, N)
    initial = np.linspace(top, bot, N)  # per node
    initialC = np.linspace(top, bot, N - 1)  # per cell
    id = range(0, N)

    # write file
    file = open(filename, "w")

    file.write("DGF\n")
    file.write('Vertex\n')
    file.write('parameters 2\n');
    for i in range(0, N):
        if vertex3d:
            file.write('{:g} {:g} {:g} {:g} {:g}\n'.format(0., 0., z_[i], initialC[i], domainId[i]))
        else:
            file.write('{:g} {:g} {:g}\n'.format(z_[i], initial[i], domainId[i]))

    file.write('#\n');
    file.write('Simplex\n');
    file.write('parameters 2\n');
    for i in range(0, N - 1):
        file.write('{:g} {:g} {:g} {:g}\n'.format(id[i], id[i + 1], initialC[i], domainId[i]));

    # not used...
    file.write('#\n')
    file.write('BOUNDARYSEGMENTS\n')  # how do i get the boundary segments into DUMUX ?
    file.write('2 0\n')
    file.write('3 {:g}\n'.format(N - 1))  # vertex id, but index starts with 0
    file.write('#\n');
    file.write('BOUNDARYDOMAIN\n')
    file.write('default 1\n');
    file.write('#\n')

    file.close()


def createDGF_1Droots(filename, nodes, seg):

    file = open(filename, "w")  # write file

    file.write("DGF\n")
    file.write('Vertex\n')
    # file.write('parameters 2\n');
    for i in range(0, len(nodes)):
        file.write('{:g} {:g} {:g} \n'.format(nodes[i, 0], nodes[i, 1], nodes[i, 2]))

    file.write('#\n');
    file.write('Simplex\n');
    # file.write('parameters 2\n');
    for i in range(0, len(seg)):
        file.write('{:g} {:g} \n'.format(seg[i, 0], seg[i, 1]));

    # not used...
    file.write('#\n')
    file.write('BOUNDARYSEGMENTS\n')  # how do i get the boundary segments into DUMUX ?
    file.write('2 0\n')
    file.write('3 {:g}\n'.format(len(seg)))  # vertex id, but index starts with 0
    file.write('#\n');
    file.write('BOUNDARYDOMAIN\n')
    file.write('default 1\n');
    file.write('#\n')

    file.close()


if __name__ == "__main__":

    L = 0.5  # length of single straight root (m)

    # create grid
    nnz = 100
    nodes = np.zeros((nnz, 3))
    seg = np.zeros(((nnz - 1), 2), dtype = int)
    c = 0
    for i in range(1, nnz):
        seg[c, 0] = i - 1
        seg[c, 1] = i
        c += 1
        nodes[i, :] = [0., 0., -i * L / (nnz - 1)]

    sn = len(seg)
    nn = len(nodes)

    createDGF_1Droots("../grids/singleroot.dgf", nodes, seg)

    print("its done.")

