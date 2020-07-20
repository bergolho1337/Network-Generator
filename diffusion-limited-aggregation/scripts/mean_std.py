import sys
import numpy as np

filename = sys.argv[1]
data = np.genfromtxt(filename)

mean = np.mean(data)
std = np.std(data)
min_value = np.min(data)
max_value = np.max(data)

print("Size = %u" % (len(data)))
print("%g $\pm$ %g" % (mean,std))
print("Min = %g" % (min_value))
print("Max = %g" % (max_value))
