# file -> numpy.loadtxt
# See: https://stackoverflow.com/questions/18259393/numpy-loading-csv-too-slow-compared-to-matlab

import time
import numpy

M = numpy.loadtxt("madelonX16.csv", delimiter=',')

runs = 100_000
cols = [9, 99, 199]

tick = time.time()
for _ in range(runs):
    column_bunch = M[:, cols].copy()
    del column_bunch
    
tock = time.time()

print(f"{runs} runs finished in: {tock-tick} sec.")
