import numpy as np
#import cppimport
#dummy_work1 = cppimport.imp("dummy_work1")
import dummy_work1

from mpi4py import MPI
comm = MPI.COMM_WORLD
comm.Barrier()
time0 = MPI.Wtime()
size = comm.Get_size()
rank = comm.Get_rank()

N = 2 ** 13
N_part = N // size

if rank == 0:
    np.random.seed(123)
    array1 = np.random.rand(N, N)
    array2 = np.random.rand(N, N)
else:
    array1 = None
    array2 = None
    
array1_part = np.empty((N_part, N))
array2_part = np.empty((N_part, N))

comm.Scatter(array1, array1_part, root=0)
comm.Scatter(array2, array2_part, root=0)

prod = np.empty(1)
prod_part = np.array(dummy_work1.scalar_prod(array1_part, array2_part))

comm.Reduce(prod_part, prod, root=0)

if rank == 0:
    print(prod[0])
    print("Finished in", MPI.Wtime() - time0, "sec.")
