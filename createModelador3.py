import os
import time

# Current dir:
currdir = os.getcwd()

# Current time
time0 = time.time()

# Compiling sir.x
os.chdir(currdir+'/SIR2015/')
os.system('make fc=gfortran modelador3.x')

# Copying sir.x inside test folder
os.system('cp modelador3.x ../.')

# Cleaning the last version
os.system('python clean.py')

# Execution time
time1 = time.time()

# Elapsed time
print('==> tc: {0:2.2f}s'.format(time1-time0))