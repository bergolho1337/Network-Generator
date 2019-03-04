import numpy as np
from pylab import *

# Gaussian parameters
sigma = 0.1
xp = 0.25

# Domain
x = np.linspace(0,2,200)

# Gaussian function
y = (1.0/(sigma*np.sqrt(2.0*np.pi)))*np.exp(-0.5*((x-xp)/sigma)**2.0)

# Plotting the data
plot(x,y)

# Output to the user
show()
