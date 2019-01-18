# file -> pandas.read_csv -> numpy.array

import time
import numpy
from pandas import read_csv

runs = 100
# cols = [9, 99, 199]

tick = time.time()
for _ in range(runs):
    M = read_csv("madelonX9.csv").values
    #column_bunch = M[:, cols].copy()
    del M
    #del column_bunch

tock = time.time()

print(f"{runs} runs finished in: {tock-tick} sec.")
