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
nlevel = 50
d_start = 200.0
gamma = 99
gamma_lv = 19
gamma_rv = 99

#for gamma in gammas:
	#level = range(nlevel)
	#d = eval_murray(d_start,gamma,nlevel)
	#for value in d:
	#	print("\t%g" % value)
	#print

level = range(nlevel)
d_lv = eval_murray(d_start,gamma_lv,nlevel)
d_rv = eval_murray(d_start,gamma_rv,nlevel)
plt.title(r"Murray's law: $d_{1}^{\gamma} = d_{21}^{\gamma} + d_{22}^{\gamma}$")
plt.ylabel(r"$d (\mu m)$",fontsize=15)
plt.xlabel(r"level",fontsize=15)
#plt.grid()
#plt.plot(level,d_lv,marker='o',label=r'$\gamma$ = %.1lf' % (gamma_lv))
plt.plot(level,d_lv,label=r'$\gamma$ = %.0lf' % (gamma_lv),linewidth=3.0,color="black")
#plt.plot(level,d_rv,label=r'$\gamma$ = %.0lf' % (gamma_rv),linewidth=3.0)
plt.legend(loc=0)
#plt.show()
plt.savefig("murray_law.svg")
