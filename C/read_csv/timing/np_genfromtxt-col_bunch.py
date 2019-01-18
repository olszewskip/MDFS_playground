# file -> numpy.genfromtxt
# See: https://stackoverflow.com/questions/18259393/numpy-loading-csv-too-slow-compared-to-matlab

import time
import numpy

runs = 100
cols = [9, 99, 199]

tick = time.time()
for _ in range(runs):
    M = numpy.genfromtxt("madelon.csv", delimiter=',')
    column_bunch = M[:, cols].copy()
    del M
    del column_bunch

tock = time.time()

print(f"{runs} runs finished in: {tock-tick} sec.")
