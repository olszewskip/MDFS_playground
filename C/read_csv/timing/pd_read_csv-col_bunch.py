# file -> pandas.read_csv -> numpy.array

import time
import numpy
from pandas import read_csv

M = read_csv("madelonX4.csv").values

runs = 100_000
cols = [9, 99, 199]

tick = time.time()
for _ in range(runs):
    column_bunch = M[:, cols].copy()
    del column_bunch

tock = time.time()

print(f"{runs} runs finished in: {tock-tick} sec.")
