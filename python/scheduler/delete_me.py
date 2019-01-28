import numpy as np
from scipy.stats import rankdata
from pandas import read_csv
from itertools import chain, starmap
# np.random.seed(123)

from mpi4py import MPI
comm = MPI.COMM_WORLD
comm.Barrier()
time0 = MPI.Wtime()
size = comm.Get_size()
rank = comm.Get_rank()
print(rank, size)