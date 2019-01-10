import time
import numpy

runs = 10
cols = [200, 500]

tick = time.time()
for _ in range(runs):
    M = numpy.loadtxt("madelon.csv", delimiter=',')
    column_bunch = M[:, cols].copy()
tock = time.time()

print(f"Average time of {runs} runs: {round((tock-tick)/runs, 3)}s")
