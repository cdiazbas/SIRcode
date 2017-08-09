import os
import time

# Current time
time0 = time.time()

# Removing last compiled version
os.system('rm *.o')
os.system('rm *.a')
os.system('rm *.x')

# Execution time
time1 = time.time()

# Elapsed time
print('==> tc: {0:2.2f}s'.format(time1-time0))