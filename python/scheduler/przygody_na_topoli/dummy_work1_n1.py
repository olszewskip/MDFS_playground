import numpy as np
import time
time_0 = time.time()

#import cppimport
#dummy_work1 = cppimport.imp("dummy_work1")
import dummy_work1

np.random.seed(0)

N = 2 ** 12
    
array1 = np.random.rand(N, N)
array2 = np.random.rand(N, N)

prod = dummy_work1.scalar_prod(array1, array2)

print(prod)
print("Finished in", time.time() - time_0, "sec.")

