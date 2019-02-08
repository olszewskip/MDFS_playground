import csv
import numpy as np
from scipy.stats import rankdata
from pandas import read_csv, DataFrame
from itertools import product, chain, starmap, combinations, combinations_with_replacement
from time import time

time0 = time()

k = 3
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

# Read and broadcast the data

file = "madelon_tiny.csv"
data = read_csv(file, dtype='float64', header=None).values.T[:-1]
data[:-1] = discretize_vec(data[:-1])
data = data.astype('int64') 

dim0, dim1 = data[:-1].shape 
n_classes = len(np.unique(data[-1]))

print("dims:", dim0, dim1, ", number of classes:", n_classes, ", k=", k)

# 2. Function definitions

def jobs_generator(k=k, dim0=dim0):
    '''
    Python-generator.
    E.g. output for k=2:
    {0,1}, {0,2}, ..., {0, dim0-1}, {1,2}, ..., {1,dim0-1}, ..., {dim0-2, dim0-1}
    Go with combinations(range(M), k) to exclude diagonal tuples
    '''        
    return combinations(range(dim0), k)    


def info(p):
    return p * np.log(p + 1E-5)

def neg_H_cond(matrix):
    return np.sum(info(matrix)) - np.sum(info(np.sum(matrix, axis=0)))

def work(indeces, n_classes=n_classes, divisions=divisions, k=k):
    '''
    Work-function.
    Output: {indexA: (number, list-of-indeces), indexB: ..., ...}
    '''
    results = {}
    
    contingency_m = np.zeros([n_classes] + [divisions] * k)
    for c_index in data[[-1] + list(indeces)].T:
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


final_results = {}

for job in jobs_generator():
    results = work(job)
    record(results, final_results)

# result
DataFrame(final_results).T.rename(columns={0: 'IG_max', 1: 'tuple'}).to_pickle("delete_me2.pkl")
print("Finished in", time()-time0, "sec.")