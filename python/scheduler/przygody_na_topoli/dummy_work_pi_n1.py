import numpy as np
#import cppimport
#dummy_work_pi = cppimport.imp("dummy_work_pi")
import dummy_work_pi

import time

time_0 = time.time()
pi = dummy_work_pi.omp_pi(1000000000)
print(pi)
print("Finished in", time.time() - time_0)
