# import time
# import numpy

# runs = 100
# cols = [9, 99, 199]

# tick = time.time()
# for _ in range(runs):
#     M = numpy.loadtxt("madelon.csv", delimiter=',')
#     column_bunch = M[:, cols].copy()
# tock = time.time()

# print(f"{runs} runs finished in: {tock-tick}s")

import time
import numpy
import pandas

runs = 100
cols = [9, 99, 199]

tick = time.time()
for _ in range(runs):
    M = pandas.read_csv("madelonX4.csv")
    column_bunch = M.iloc[:, cols].copy()
tock = time.time()

print(f"{runs} runs finished in: {tock-tick}s")
