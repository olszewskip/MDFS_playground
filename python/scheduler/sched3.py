import csv
import numpy as np
from scipy.stats import rankdata
#from pandas import read_csv, DataFrame
from itertools import product, chain, starmap, combinations, combinations_with_replacement
import pickle

from mpi4py import MPI
comm = MPI.COMM_WORLD
comm.Barrier()
time0 = MPI.Wtime()
size = comm.Get_size()
rank = comm.Get_rank()

k = 3
window = 10
divisions = 10
range_ = 0.00
seed = 123

# 1. Function definitions

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
    discrete_seq = np.zeros(len(seq), dtype='float64')
    for treshold in tresholds:
        discrete_seq[ranks > treshold] += 1
    return discrete_seq

discretize_vec = np.vectorize(discretize, signature='(n)->(n)', excluded=['divisions', 'range_', 'seed'])

# Read the data in each rank

file = "madelon_tiny.csv"
#data = read_csv(file, dtype='float64', header=None).values.T[:-1]

data = []
with open(file) as csvfile:
    reader = csv.reader(csvfile, delimiter=',',
                        quoting=csv.QUOTE_NONNUMERIC)
    for row in reader:
        data.append(row)
data = np.array(data, dtype='float64').T[:-1]

data[:-1] = discretize_vec(data[:-1])
data = data.astype('int64')


dim0, dim1 = data[:-1].shape
M = (dim0 - 1) // window + 1
border_cols = range( (M-1) * window, dim0)

n_classes = len(np.unique(data[-1]))

if rank == 0:
    print(rank, "dims:", dim0, dim1, ", tile_index_span:", M, ", number of classes:", n_classes, " window:", window)
    print(rank, "Elapsed:", MPI.Wtime() - time0, "sec")

# 2. Function definitions

def tiles_generator(k=k, M=M):
    '''
    Python-generator.
    E.g. output for k=2:
    {0,0}, {0,1}, ..., {0, M-1}, {1,1}, ..., {1,M-1}, ..., {M-1, M-1}
    Go with combinations(range(M), k) to exclude diagonal tuples
    '''        
    return combinations_with_replacement(range(M), k)    


def jobs_generator(tile, window=window, M=M, border_cols=border_cols):
    '''
    Map tile into sequence of fundamental-tiles
    (i.e. elements of the cartesian product of data columns)
    '''
    index_counts = {index: tile.count(index) for index in tile}
    index_to_cols = lambda index: range(index * window, (index + 1) * window) if index != (M - 1) else border_cols 
    cols_tile = (combinations(index_to_cols(index), count) for (index, count) in index_counts.items())
    return (list(chain.from_iterable(col_indeces)) for col_indeces in product(*cols_tile))


def dummy_work(indeces):
    '''
    Work-function. A dummy version.
    Output: {indexA: (number, list-of-indeces), indexB: ..., ...}
    '''
    results = {}
    for index in indeces:
        other_indeces = list(set(indeces) - set([index]))
        result = np.sum(data[other_indeces]) - np.sum(data[index])
        results[index] = (result, other_indeces)
    return results


def neg_H(p):
    return p * np.log(p + 1E-5)

def neg_H_cond(matrix):
    return np.sum(neg_H(matrix)) - np.sum(neg_H(np.sum(matrix, axis=0)))

def work(indeces, n_classes=n_classes, divisions=divisions, k=k):
    '''
    Work-function.
    Output: {indexA: (number, list-of-indeces), indexB: ..., ...}
    indeces -> list
    '''
    results = {}
    
    contingency_m = np.zeros([n_classes] + [divisions] * k)
    for c_index in data[[-1] + indeces].T:
        contingency_m[tuple(c_index)] += 1
        
    for i, index in enumerate(indeces):
        result = neg_H_cond(contingency_m) - neg_H_cond(np.sum(contingency_m, axis=i+1))
        results[index] = (result, indeces)
    return results


def record(results, records):
    '''
    results, records -> dicts
    Accepts output of the work-function and updates the dict that accumulates global results
    '''
    for index, score in results.items():
        if index not in records or score[0] > records[index][0]:
            records[index] = score

# Work loop

if rank == 0:
    final_results = {}
    current_assignements = {rank: (0, None) for rank in range(1,size)}
    
    print(rank, "entering the for loop")
    status = MPI.Status()
    for tile in tiles_generator(k, M):
        tile_results = comm.recv(status=status)
        record(tile_results, final_results)
        comm.isend(tile, dest=status.source)
        job_count = current_assignements[status.source][0]
        current_assignements[status.source] = (job_count + 1, tile)
        
        #print(rank, "currently:", current_assignements)
    
    print(rank, "Work queue is empty")
    for _ in range(size - 1):
        tile_results = comm.recv(status=status)
        record(tile_results, final_results)
        comm.isend(None, dest=status.source)

    # Save the results to a file
    with open("delete_me3.pkl", "wb") as file:
        pickle.dump(final_results, file)

    # DataFrame(final_results).T.rename(columns={0: 'IG_max', 1: 'tuple'}).to_pickle("delete_me3.pkl")

    print(rank, "says goodbye")
    # print("Final_results:", final_results)
    print(rank, "Elapsed:", MPI.Wtime() - time0, "sec")
    
else:
    comm.send({}, dest=0)
    print(rank, "entering the while loop")
    while True:
        tile = comm.recv(source = 0)
        if tile:
            tile_results = {}
            for job in jobs_generator(tile):
                results = work(job)
                record(results, tile_results)
            comm.isend(tile_results, dest=0)
        else:
            print(rank, "says goodbye")
            break