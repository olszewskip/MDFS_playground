import numpy as np
from scipy.stats import rankdata
from pandas import read_csv
from itertools import chain, starmap

from mpi4py import MPI
comm = MPI.COMM_WORLD
comm.Barrier()
time0 = MPI.Wtime()
size = comm.Get_size()
rank = comm.Get_rank()

k = 3
window = 10
divisions = 10
range_ = 0.1
seed = 123

print(rank, "starting")

def discretize(seq, divisions=divisions, range_=range_, seed=seed):
    '''
    >>> discretize([3, 4, 1, 8, 13, 8], divisions=4, range_=0, seed=123) = array([1, 1, 0, 2, 3, 2])
    where
    ranks = [2., 3., 1., 4.5, 6., 4.5]
    tresholds = [1.5,  3.,  4.5]
    '''
    np.random.seed(seed)
    ranks = rankdata(seq)
    tresholds = (range_ * (np.random.random(divisions - 1) - .5) + np.arange(1, divisions)) / divisions * len(seq)
    discrete_seq = np.zeros(len(seq), dtype='int64')
    for treshold in tresholds:
        discrete_seq[ranks > treshold] += 1
    return discrete_seq

discretize_vec = np.vectorize(discretize, signature='(n)->(n)', excluded=['divisions', 'range_', 'seed'])

if rank == 0:
    print(rank, "attempting to read data")
    file = "madelon_tiny.csv"
    data = discretize_vec(read_csv(file, dtype='float64', header=None).values.T[:-2])
    data_shape = data.shape
else:
    data_shape = None    

dim0, dim1 = comm.bcast(data_shape, root=0)
M = (dim0 - 1) // window + 1

print(rank, "dims:", dim0, dim1, "tile_index_span:", M)

if rank != 0:
    data = np.empty((dim0, dim1), dtype='int64')

comm.Bcast([data, MPI.INT], root=0)
if rank == 0:
    print("Elapsed:", MPI.Wtime() - time0, "sec")
    
    
def tiles_generator(k, M, skip_diag=True):
    '''
    Python-generator.
    E.g. output for k=2 and skippig the diagonal elements:
    {0,1}, {0,2}, ..., {0, M-1}, {1,2}, ..., {1,M-1}, ..., {M-2, M-2}
    '''
    sd = int(skip_diag)

    def embedd(*indeces):
        for i in range(indeces[-1] + sd, M):
            yield (*indeces, i)

    indeces = ((idx,) for idx in range(M))
    for _ in range(k - 1):
        indeces = chain.from_iterable(starmap(embedd, indeces))
        
    return indeces

def dummy_job(indeces):
    '''
    Job function. A dummy version.
    Output: {indexA: (number, list-of-indeces), indexB: ..., ...}
    '''
    results = {}
    for index in indeces:
        other_indeces = list(set(indeces) - set([index]))
        result = np.sum(data[other_indeces]) - np.sum(data[index])
        results[index] = (result, other_indeces)
    return results



if rank == 0:
    final_results = {index: (0, None) for index in range(dim0)}
    def record(results):
        '''
        Accepts output of the job function and updates the final_results dict
        '''
        for index, score in results.items():
            if score[0] > final_results[index][0]:
                final_results[index] = score
    
    print(rank, "entering the for loop")
    status = MPI.Status()
    for tile in tiles_generator(k, dim0):
        results = comm.recv(status=status)
        comm.isend(tile, dest=status.source)
        record(results)
    
    print("Elapsed:", MPI.Wtime() - time0, "sec")
    print(rank, "starts terminating workers")
    for _ in range(size - 1):
        last_results = comm.recv(status=status)
        record(last_results)
        comm.isend(None, dest=status.source)
    
    print(rank, "says goodbye")
    print("final_results:", final_results)
        
    
else:
    print(rank, "attempting to send blank_results")
    comm.send({0: (0,)}, dest=0)
    print(rank, "entering the while loop")
    while True:
        tile = comm.recv(source = 0)
        if tile:
            results = dummy_job(tile)
        else:
            print(rank, "says goodbye")
            break
        comm.isend(results, dest=0)
