import numpy as np
import time
import dummy_work1

time_0 = time.time()

N = 2 ** 13

np.random.seed(0)
#array1 = np.random.rand(N, N)
#array2 = np.random.rand(N, N)

#prod = dummy_work1.scalar_prod(array1, array2)
prod = dummy_work1.parallel_pi(1000000000)

print("Finished in", time.time() - time_0, "sec.")
print(prod)