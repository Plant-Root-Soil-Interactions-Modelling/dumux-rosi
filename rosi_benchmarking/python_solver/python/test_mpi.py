from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


def magic_formula(i :int):
    return 2 * i


N = 1000
if rank == 0:
    data = []
    c = 0
    chunk_size = round(N / size - 0.5)
    for i in range(0, size - 1):
        data.append(list(range(c, c + chunk_size)))
        c += chunk_size
    data.append(list(range(c, N)))
else:
    data = None
data = comm.scatter(data, root = 0)
print(data)

res = []
for d in data:
    res.append(magic_formula(d))

res = comm.gather(res, root = 0)
if rank == 0:
    print(res)

