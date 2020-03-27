import matplotlib as plt

#r_p = 1.0
r_p = 0.00102269
gamma = 3.0
factor = 0.5**(1.0/gamma)

print("gamma = %g" % gamma)
print("factor = %g" % factor)
for i in range(50):
	print("r_p = %g" % r_p)
	r_f = factor * r_p
	r_p = r_f
