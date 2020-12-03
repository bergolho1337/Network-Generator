import numpy as np
from matplotlib import pyplot as plt

Gi = 7.9
Cf = 3.4
tauf = 0.1

#Gi = 4.8
#Cf = 3.4
#tauf = 0.1

def Gm (d):
    return 0.005*d + 0.2

def lambda_m (Rm,Ri,d):
    return ( (Rm*d)/(4.0*Ri) )**(0.5)

def Ri (d):
    D = d * 1.0e-04
    return (D**3 * Gm(d) * Gi)**(-0.5) * (2.0/3.14159)

def v (d):
    return ( (Gi*d)/(4.0*Cf*tauf) )**(0.5) * 0.1

d = np.linspace(0,300,100)
s = v(d)

plt.plot(d,s,color='black',linewidth=3.0)
plt.title("Cable equation - Purkinje fiber")
plt.xlabel(r"d ($\mu m$)",fontsize=15)
plt.ylabel(r"v ($m/s$)",fontsize=15)
plt.show()
#plt.savefig("cable_equation.svg")
