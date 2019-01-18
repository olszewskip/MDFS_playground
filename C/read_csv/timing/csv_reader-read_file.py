# file -> list of lists -> numpy.array

import time
import csv
import numpy as np

runs = 100
cols = [9, 99, 199]

tick = time.time()

for _ in range(runs):
    M = []
    with open('madelon.csv', newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
        for row in reader:
            M.append(row)
    M = np.array(M)
    # column_bunch = M[:, cols].copy()
    del M
    # del column_bunch

tock = time.time()

print(f"{runs} runs finished in: {tock-tick} sec.")
