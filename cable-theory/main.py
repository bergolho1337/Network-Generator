import numpy as np
from matplotlib import pyplot as plt

Gi = 7.9

def Gm (d):
    return 0.005*d + 0.2

def lambda_m (Rm,Ri,d):
    return ( (Rm*d)/(4.0*Ri) )**(0.5)

def Ri (d):
    D = d * 1.0e-04
    return (D**3 * Gm(d) * Gi)**(-0.5) * (2.0/3.14159)

d = np.linspace(0,300,100)
G_m = Gm(d)
R_i = Ri(d)

#plt.plot(d,G_m)
#plt.ylim([0,2])
plt.plot(d,R_i)
#plt.ylim([0,800])
plt.show()