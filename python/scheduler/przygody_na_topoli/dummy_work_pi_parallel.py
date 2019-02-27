import numpy as np
#import cppimport
#dummy_work_pi = cppimport.imp("dummy_work_pi")
import dummy_work_pi

from mpi4py import MPI
comm = MPI.COMM_WORLD
comm.Barrier()
time0 = MPI.Wtime()
size = comm.Get_size()
rank = comm.Get_rank()

pi = np.empty(1)
pi_part = np.array(dummy_work_pi.omp_pi(1000000000, rank, size))

comm.Reduce(pi_part, pi, root=0)

if rank == 0:
    print(pi[0])
    print("Finished in", MPI.Wtime() - time0, "sec.")
    