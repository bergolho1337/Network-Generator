import math

d = 0.0041153
r = d / 2.0
i_max = 100
l = 1.0
rho = 1.724e-08

R = (rho * l) / (math.pi * r * r)

print(R)

V_in = 1.0
V_out = V_in - 0.01296

print(V_out)
