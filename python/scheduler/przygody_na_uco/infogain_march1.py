import csv
import numpy as np
from scipy.stats import rankdata
from itertools import product, chain, starmap, combinations, combinations_with_replacement
import pickle

import wrap_discretize
import infogain

from mpi4py import MPI
comm = MPI.COMM_WORLD
comm.Barrier()
time0 = MPI.Wtime()
size = comm.Get_size()
rank = comm.Get_rank()

k = 3
window = 5
divisions = 31
range_ = 0.0
seed = 123
    
wrap_discretize_switch = True

# 1. Read in the data

file = "madelon_tiny.csv"
input_ = []
with open(file) as csvfile:
    reader = csv.reader(csvfile, delimiter=',',
                        quoting=csv.QUOTE_NONNUMERIC)
    for row in reader:
        input_.append(row)
        
input_ = np.ascontiguousarray(np.array(input_, dtype='float64').T[:-1])
data = np.zeros_like(input_, dtype='uint8')

# Discretize the data

if wrap_discretize_switch:

    for col_idx in range(input_.shape[0] - 1):
        wrap_discretize.discretize(
            seed = 123,
            discretization_index = 0,
            feature_id = col_idx,  # ?
            divisions = divisions,
            object_count = input_.shape[1],
            py_in_data = input_[col_idx],
            py_sorted_in_data = np.sort(input_[col_idx]),
            py_out_data = data[col_idx],
            range_ = range_
        )
else:
    
    def discretize(seq, divisions=divisions, range_=range_, seed=seed):
        '''
        >>> discretize([3, 4, 1, 8, 13, 8], divisions=4, range_=0, seed=123) = array([1, 1, 0, 2, 3, 2])
        where
        ranks = [2., 3., 1., 4.5, 6., 4.5]
        tresholds = [1.5,  3.,  4.5]
        '''
        np.random.seed(seed)
        ranks = rankdata(seq, method='average') # method='ordinal'/'average' ?

        random_blocks = np.cumsum(range_ * (2 * np.random.random(divisions + 1) - 1) + np.ones(divisions + 1))
        tresholds = random_blocks[:-1] / random_blocks[-1] * len(seq)

        discrete_seq = np.zeros(len(seq), dtype='uint8')
        for treshold in tresholds:
            discrete_seq[ranks > treshold] += 1
        return discrete_seq

    discretize_vec = np.vectorize(discretize, signature='(n)->(n)', excluded=['divisions', 'range_', 'seed'])
    data[:-1] = discretize_vec(input_[:-1])

    
data[-1] = input_[-1:].astype('uint8')

labels, counts = np.unique(data[-1], return_counts=True)
n_classes = len(labels)

xi = 0.25
pseudo_counts = xi * counts / np.min(counts)

dim0, dim1 = data[:-1].shape
M = (dim0 - 1) // window + 1
border_cols = range( (M-1) * window, dim0)

if rank == 0:
    print(rank, "dim0 =", dim0, "; dim1 =", dim1)

# 3. More function definitions

def tile_generator(k=k, M=M):
    '''
    Python-generator.
    E.g. output for k=2:
    {0,0}, {0,1}, ..., {0, M-1}, {1,1}, ..., {1,M-1}, ..., {M-1, M-1}
    Go with combinations(range(M), k) to exclude diagonal tuples
    '''        
    return combinations_with_replacement(range(M), k)    

def tuple_generator(tile, window=window, border_cols=border_cols):
    '''
    Map tile into sequence of k-tuples, i.e. fundamental-tiles,
    i.e. elements of the cartesian product of the data-columns
    '''
    index_counts = {index: tile.count(index) for index in tile}
    index_to_cols = lambda index: range(index * window, (index + 1) * window) if index != (M - 1) else border_cols 
    cols_tile = (combinations(index_to_cols(index), count) for (index, count) in index_counts.items())
    return (list(chain.from_iterable(col_indeces)) for col_indeces in product(*cols_tile))

def neg_H(p):
    return p * np.log2(p)

def neg_H_cond(matrix):
    return np.sum(neg_H(matrix)) - np.sum(neg_H(np.sum(matrix, axis=-1)))

def slow_work(tuple_, divisions=divisions, n_classes=n_classes, pseudo_counts=pseudo_counts):
    '''
    tuple_ -> list # dammit...
    Work-function.
    Output: tuple of Information Gains implicitly corresponding to column-indeces in the tuple_
    '''
    # contingency-matrix: begin with pseudo-counts
    contingency_m = np.empty([divisions + 1] * k + [n_classes], dtype='float64')
    for label, pseudo_count in enumerate(pseudo_counts):
        contingency_m[..., label] = pseudo_count
    
    # contingency-matrix: normal counts
    for c_index in data[tuple_ + [-1]].T:
        contingency_m[tuple(c_index)] += 1
    
    IGs = tuple(neg_H_cond(contingency_m) - neg_H_cond(np.sum(contingency_m, axis=i)) for i in range(len(tuple_)))
    return IGs


def record_tuple(tuple_, IGs, records):
    for column, IG in zip(tuple_, IGs):
        if column not in records or IG > records[column][0]:
            records[column] = (IG, tuple_)
            
def record_tile(tile_results, records):
    for column, (IG, tuple_) in tile_results.items():
        if column not in records or IG > records[column][0]:
            records[column] = (IG, tuple_)

# 4 Work loop

if rank == 0:
    final_results = {}
    current_assignements = {rank: (0, None) for rank in range(1,size)}
    
    print(rank, "entering the for loop")
    status = MPI.Status()
    for tile in tile_generator():
        tile_results = comm.recv(status=status)
        record_tile(tile_results, final_results)
        comm.isend(tile, dest=status.source)
        job_count = current_assignements[status.source][0]
        current_assignements[status.source] = (job_count + 1, tile)
        
        #print(rank, "currently:", current_assignements)
    
    print(rank, "Work queue is empty")
    for _ in range(size - 1):
        tile_results = comm.recv(status=status)
        record_tile(tile_results, final_results)
        comm.isend(None, dest=status.source)

    # Save the results to a file
    with open("march1_mpi_results.pkl", "wb") as file:
        pickle.dump(final_results, file)

    print(rank, "says goodbye")
    print(rank, "Elapsed:", MPI.Wtime() - time0, "sec")
    
else:
    comm.send({}, dest=0)
    print(rank, "entering the while loop")
    while True:
        tile = comm.recv(source = 0)
        try:
            tile_results = {}
            for tuple_ in tuple_generator(tile):
                
                #IGs = slow_work(tuple_)
                IGs = infogain.work_3a(dim1, divisions, data[tuple_[0]], data[tuple_[1]], data[tuple_[2]], n_classes, pseudo_counts, data[-1])
                #IGs = infogain.work_3b(dim1, divisions, data[tuple_[0]], data[tuple_[1]], data[tuple_[2]], n_classes, pseudo_counts, data[-1])
                #IGs = infogain.work_3c(dim1, divisions, data[tuple_[0]], data[tuple_[1]], data[tuple_[2]], n_classes, pseudo_counts, data[-1])
                record_tuple(tuple_, IGs, tile_results)
            comm.isend(tile_results, dest=0)
        except:
            print(rank, "says goodbye")
            break
