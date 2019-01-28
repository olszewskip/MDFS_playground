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

#tag_idle = 100
tag_results = 101
tag_terminate = 102


divisions = 10
range_ = 0.1
seed = 123
k = 2
window = 10
# self.M = (N-1) // self.w + 1

print(rank, "stage 1")

def discretize(seq, divisions=divisions, range_=range_, seed=seed):
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
    file = "madelon.csv"
    data = discretize_vec(read_csv(file, dtype='float64', header=None).values.T[:-2])
    data_shape = data.shape
else:
    data_shape = None    

dim0, dim1 = comm.bcast(data_shape, root=0)
print(rank, "dims:", dim0, dim1)

if rank != 0:
    data = np.empty((dim0, dim1), dtype='int64')

comm.Bcast([data, MPI.INT], root=0)
if rank == 0:
    print(MPI.Wtime() - time0)
    
def tiles_generator(k, M, skip_diag=True):
    sd = int(skip_diag)

    def embedd(*indeces):
        for i in range(indeces[-1] + sd, M):
            yield (*indeces, i)

    indeces = ((idx,) for idx in range(M))
    for _ in range(k - 1):
        indeces = chain.from_iterable(starmap(embedd, indeces))
        
    return indeces

def dummy_job(indeces):
    results = {}
    for index in indeces:
        other_indeces = tuple(set(indeces) - set([index]))
        result = np.mean(data[other_indeces]) - np.mean(data[index])
        results[index] = (result, other_indeces)
    return results


local_tiles = tiles_generator(k, dim0)
status = MPI.Status()
print(rank, "stage 2")

if rank == 0:
    final_results = {index: (0, None) for index in range(dim0)}
    for _ in range(size - 1):
        next(local_tiles)
    workers = set(range(1,size))
    print(rank, "entering while loop")
#     for worker in workers:
#         results = comm.recv(source = MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
#         print(rank, "received from", status.source)
    while workers:
        results = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        print(rank, "received", status.tag, "from", status.source)
#         for index, score in results.items():
#             if score[0] > final_results[index][0]:
#                 final_results[index] = score
        comm.send(None, dest=status.source, tag=tag_terminate)
        workers = workers - set([status.source])
    print(rank, "says goodbye")
        
    
else:
    for _ in range(rank - 1):
        next(local_tiles)
    first_tile = next(local_tiles)
    first_results = dummy_job(first_tile)
    print(rank, "attempting to send first results")
    comm.send(first_results, dest=0, tag=tag_results)
    print(rank, "entering while loop")
    while True:
        tile = comm.recv(source = 0, tag=MPI.ANY_TAG, status=status)
        if status.tag == tag_terminate:
            print(rank, "says goodbye")
            break
        else:
            results = dummy_job(tile)
        comm.send(results, dest=0, tag=tag_results)



# sample = np.zeros(5, dtype='float64')

# if rank == 1:
#     sample = np.ones(5, dtype='float64')
#     comm.Send([sample, MPI.DOUBLE], dest=0, tag=tag_idle)
# elif rank == 0:
#     #for _ in range(1, size):
#     status = MPI.Status()
#     comm.Probe(MPI.ANY_SOURCE, MPI.ANY_TAG, status=status)
#     print(status.tag)
#     print(status.source)
#     print(status.Get_elements(MPI.DOUBLE))
#     comm.Recv(sample, source=status.source)


# if rank == 0:
#     status = MPI.Status()
#     task_queue = range( (dim0 - 1) // window + 1)
#     while(task_queue)
# else:
#     comm.send(None, dest=0, tag=tag_idle)
    
#     request = comm.Irecv(MPI.ANY_SOURCE, 
#     packet = data[task*window: (task+1)*window]
#     comm.Send(data, )


# if rank != 0:
#     data = np.empty(data_shape, dtype=dtype)

# if rank == 0:
#     colsA = [0, 1]
#     colsB = [2, 3]
#     column_bunch = data[:, colsA]
#     comm.Send(column_bunch, dest=1, tag=0)
#     .
# elif rank == 1:
#     column_bunch = np.empty((dim0, 2), dtype=dtype)
#     comm.Recv(column_bunch, source=0, tag=0)

# print(rank, column_bunch)