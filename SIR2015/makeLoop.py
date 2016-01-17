import os
import time
import numpy as np

# Time compilation array
tcArray = []

# Time execution array
teArray = []

# n ite
nite = 3

for i in range(nite):
	# Current time
	time0 = time.time()

	# Remove last compiled version
	os.system('rm *.o')
	os.system('rm *.a')
	os.system('rm *.x')

	# Compile sir.x
	os.system('make fc=gfortran')

	# Copy sir.x inside test folder
	os.system('cp sir.x ../test_CARLOS/.')

	# Execution time
	time1 = time.time()

	# Execute sir.x
	os.chdir('../test_CARLOS')
	os.system('echo sir.trol | ./sir.x')

	# Elapsed time
	print('==> tc: {0:2.2f}s'.format(time1-time0))
	print('==> te: {0:2.2f}s'.format(time.time()-time1))

	# Save the results
	tcArray.append(time1-time0)
	teArray.append(time.time()-time1)

	os.chdir('../SIR2015')

# Array conversion
tcArray = np.array(tcArray)
teArray = np.array(teArray)


print('TC ==> m: {0:2.2f} std: {1:2.2f}'.format(np.mean(tcArray),np.std(tcArray)))
print('TE ==> m: {0:2.2f} std: {1:2.2f}'.format(np.mean(teArray),np.std(teArray)))


# -O4 -ffixed-line-length-none -fno-automatic
# 521 KB
# TC ==> m: 17.11 std: 1.34
# TE ==> m: 2.95 std: 0.24

# -O4 -ffixed-line-length-none -fno-automatic -ffast-math
# 534 KB
# TC ==> m: 16.72 std: 0.38
# TE ==> m: 2.21 std: 0.03

# -O4 -ffixed-line-length-none -fno-automatic -ffast-math -finline-functions
# 534 KB
# TC ==> m: 16.39 std: 0.28
# TE ==> m: 2.21 std: 0.03

# -O4 -ffixed-line-length-none -fno-automatic -ffast-math -funsafe-math-optimizations -ffinite-math-only -fno-trapping-math
# 534 KB
# TC ==> m: 16.26 std: 0.15
# TE ==> m: 2.20 std: 0.03

# -Ofast -ffixed-line-length-none -fno-automatic
# 530 KB
# TC ==> m: 16.69 std: 0.44
# TE ==> m: 2.37 std: 0.04

# ++ -funroll-loops
# 755 KB
# TC ==> m: 22.31 std: 0.34
# TE ==> m: 2.19 std: 0.03

# ++ -fstrength-reduce
# 755 KB
# TC ==> m: 22.29 std: 0.31
# TE ==> m: 2.19 std: 0.03

# -O4 -ffixed-line-length-none -fno-automatic -ffast-math -fstrength-reduce -fexpensive-optimizations
# 534 KB
# TC ==> m: 16.96 std: 0.31
# TE ==> m: 2.24 std: 0.04

# -O4 -ffixed-line-length-none -fno-automatic -ffast-math -funroll-loops
# 755 KB
# TC ==> m: 22.55 std: 0.31
# TE ==> m: 2.20 std: 0.03 ==> Mejora del 25 porciento

# 2015 ONLY
# 433 KB
# TC ==> m: 12.68 std: 0.35
# TE ==> m: 3.48 std: 0.10

# 2015 -ffast-math -funroll-loops
# 659 KB
# TC ==> m: 17.91 std: 0.25
# TE ==> m: 2.73 std: 0.02 ==> Mejora del 20 porciento

