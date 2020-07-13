import sys
import numpy as np

filename = sys.argv[1]
data = np.genfromtxt(filename)

mean = np.mean(data)
std = np.std(data)

print("%.2lf $\pm$ %.2lf" % (mean,std))
