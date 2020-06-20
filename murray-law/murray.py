import matplotlib.pyplot as plt

def eval_murray (d_start,gamma,nlevel):
	d = []
	cur_d = d_start
	factor = 0.5**(1.0/gamma)
	print("alpha = %g -- gamma = %g" % (factor,gamma))

	for i in range(nlevel):
		d.append(cur_d)
		cur_d = cur_d * factor
	
	return d

gammas = range(10,20)
nlevel = 20
d_start = 200.0

#for gamma in gammas:
	#level = range(nlevel)
	#d = eval_murray(d_start,gamma,nlevel)
	#for value in d:
	#	print("\t%g" % value)
	#print

level = range(nlevel)
d = eval_murray(d_start,19,nlevel)
plt.title(r"Murray law: $r_{1}^{\gamma} = r_{21}^{\gamma} + r_{22}^{\gamma}$")
plt.xlim([0,20])
plt.ylabel(r"$d (\mu m)$",fontsize=15)
plt.xlabel(r"level",fontsize=15)
plt.grid()
plt.plot(level,d,marker='o',label=r'$\gamma$ = 19')
plt.legend(loc=0)
#plt.show()
plt.savefig("purkinje-murray.pdf")
