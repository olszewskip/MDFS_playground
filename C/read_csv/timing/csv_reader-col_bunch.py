# file -> list of lists -> numpy.array

import time
import csv
import numpy as np

runs = 100_000
cols = [9, 99, 199]

M = []
with open('madelonX16.csv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',',
                        quoting=csv.QUOTE_NONNUMERIC)
    for row in reader:
        M.append(row)
M = np.array(M)

tick = time.time()

for _ in range(runs):
    column_bunch = M[:, cols].copy()
    del column_bunch

tock = time.time()

print(f"{runs} runs finished in: {tock-tick} sec.")
