import csv
import numpy as np
from scipy.stats import rankdata
from pandas import read_csv, DataFrame
from itertools import product, chain, starmap, combinations, combinations_with_replacement
from time import time
import pickle

import wrap_discretize
import fast

time0 = time()

k = 3
divisions = 1
range_ = 0.00
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
            discretization_index = 4,
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

# 3. More function definitions

def tuple_generator(k=k, dim0=dim0):
    '''
    Python-generator.
    E.g. output for k=2:
    {0,1}, {0,2}, ..., {0, dim0-1}, {1,2}, ..., {1,dim0-1}, ..., {dim0-2, dim0-1}
    Go with combinations_with_replacement() to include diagonal tuples like {0, 0}
    '''        
    return combinations(range(dim0), k)    

def neg_H(p):
    return p * np.log2(p)

def neg_H_cond(matrix):
    return np.sum(neg_H(matrix)) - np.sum(neg_H(np.sum(matrix, axis=-1)))

def slow_work(indeces):
    '''
    indeces -> tuple
    Work-function.
    Output: tuple of Information Gains implicitly corresponding to the indeces
    '''
    # contingency-matrix: begin with pseudo-counts
    contingency_m = np.empty([divisions + 1] * k + [len(labels)], dtype='float64')
    for label, pseudo_count in enumerate(pseudo_counts):
        contingency_m[..., label] = pseudo_count
    
    # contingency-matrix: normal counts
    for c_index in data[list(indeces) + [-1]].T:
        contingency_m[tuple(c_index)] += 1
    
    results = []
    for i in range(len(indeces)):
        result = neg_H_cond(contingency_m) - neg_H_cond(np.sum(contingency_m, axis=i))
        results.append(result)
    return tuple(results)


def record(tuple_, IGs, records):
    for column, IG in zip(tuple_, IGs):
        if column not in records or IG > records[column][0]:
            records[column] = (IG, tuple_)

final_results = {}

for tuple_ in tuple_generator():
    #IGs = slow_work(tuple_)
    #IGs = fast.work_3a(dim1, divisions, data[tuple_[0]], data[tuple_[1]], data[tuple_[2]], n_classes, pseudo_counts, data[-1])
    #IGs = fast.work_3b(dim1, divisions, data[tuple_[0]], data[tuple_[1]], data[tuple_[2]], n_classes, pseudo_counts, data[-1])
    #IGs = fast.work_3c(dim1, divisions, data[tuple_[0]], data[tuple_[1]], data[tuple_[2]], n_classes, pseudo_counts, data[-1])
    IGs = fast.work_3_old(dim1, divisions, data[tuple_[0]], data[tuple_[1]], data[tuple_[2]], n_classes, pseudo_counts, data[-1])
    record(tuple_, IGs, final_results)

# result
print("Finished in", time() - time0, "sec.")

with open("march1_n1_results.pkl", "wb") as file:
    pickle.dump(final_results, file)