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

np.random.seed(rank)

N = 2 ** 13
N_part = N // size
    
array1_part = np.random.rand(N_part, N)
array2_part = np.random.rand(N_part, N)

prod = np.empty(1)
prod_part = np.array(dummy_work1.scalar_prod(array1_part, array2_part))

comm.Reduce(prod_part, prod, root=0)

if rank == 0:
    print(prod[0])
    print("Finished in", MPI.Wtime() - time0, "sec.")
