import sys
import numpy as np

data = np.genfromtxt("data/table.dat")
data = data.transpose()

for line in data:
	for value in line:
		print("%.2lf " % (value)),
	print("")
