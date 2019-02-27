import time
import my_pi

time_0 = time.time()
pi = my_pi.parallel(1000000000)
print(pi)
print("Finished in", time.time() - time_0)