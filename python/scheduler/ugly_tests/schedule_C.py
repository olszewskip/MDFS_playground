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
window = 5
divisions = 1
range_ = 0.0
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
    ranks = rankdata(seq, method='ordinal') # method='ordinal'/'average' ?
    
    random_blocks = np.cumsum(range_ * (2 * np.random.random(divisions + 1) - 1) + np.ones(divisions + 1))
    tresholds = random_blocks[:-1] / random_blocks[-1] * len(seq)
    
    discrete_seq = np.zeros(len(seq), dtype='float64')
    for treshold in tresholds:
        discrete_seq[ranks > treshold] += 1
    return discrete_seq

discretize_vec = np.vectorize(discretize, signature='(n)->(n)', excluded=['divisions', 'range_', 'seed'])

# 2. Read the data in each rank

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
labels, counts = np.unique(data[-1], return_counts=True)
label_counts = {int(label): label_count for (label, label_count) in zip(labels, counts)}
min_count = np.min(counts)

dim0, dim1 = data[:-1].shape
M = (dim0 - 1) // window + 1
border_cols = range( (M-1) * window, dim0)


if rank == 0:
    print(rank, "dims:", dim0, dim1, ", tile_index_span:", M, ", labels:", label_counts, " window:", window)
    print(rank, "Elapsed:", MPI.Wtime() - time0, "sec")

# 3. More function definitions

def tiles_generator(k=k, M=M):
    '''
    Python-generator.
    E.g. output for k=2:
    {0,0}, {0,1}, ..., {0, M-1}, {1,1}, ..., {1,M-1}, ..., {M-1, M-1}
    Go with combinations(range(M), k) to exclude diagonal tuples
    '''        
    return combinations_with_replacement(range(M), k)    


def jobs_generator(tile):
    '''
    Map tile into sequence of fundamental-tiles
    (i.e. elements of the cartesian product of data columns)
    '''
    index_counts = {index: tile.count(index) for index in tile}
    index_to_cols = lambda index: range(index * window, (index + 1) * window) if index != (M - 1) else border_cols 
    cols_tile = (combinations(index_to_cols(index), count) for (index, count) in index_counts.items())
    return (list(chain.from_iterable(col_indeces)) for col_indeces in product(*cols_tile))


def neg_H(p):
    return p * np.log2(p)

def neg_H_cond(matrix):
    return np.sum(neg_H(matrix)) - np.sum(neg_H(np.sum(matrix, axis=-1)))

xi = 1e-5
def work(indeces):
    '''
    Work-function.
    Output: {indexA: (number, list-of-indeces), indexB: ..., ...}
    indeces -> list
    '''
    
    # contingency-matrix: begin with pseudo-counts
    contingency_m = np.ones([divisions + 1] * k + [len(labels)], dtype='float64')
    for label, count in label_counts.items():
        contingency_m[..., label] *= xi * (count / min_count)
        
    # contingency-matrix: normal counts
    for c_index in data[indeces + [-1]].T:
        contingency_m[tuple(c_index)] += 1

    results = {}
    for i, index in enumerate(indeces):
        result = neg_H_cond(contingency_m) - neg_H_cond(np.sum(contingency_m, axis=i))
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

# 4 Work loop

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
    with open("final_results_C.pkl", "wb") as file:
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
                # print(rank, results)
                record(results, tile_results)
            comm.isend(tile_results, dest=0)
        else:
            print(rank, "says goodbye")
            break