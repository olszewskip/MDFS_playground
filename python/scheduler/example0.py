from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

a_size = 4

recvdata = np.empty(a_size,dtype=np.float64)
senddata = None
counts = None
dspls = None
if rank == 0:
    senddata = np.arange(100, dtype=numpy.float64)
    counts=(1,2,3)
    dspls=(4,3,10)
comm.Scatterv([senddata, counts, dspls, MPI.DOUBLE], recvdata,root=0)

print(rank, recvdata)