L = 0.0287941
L_2 = L/2
L_4 = L/4

print("%g %g %g" % (-L,0,0))
print("%g %g %g" % (0,0,0))
print("%g %g %g" % (L_2,L_2,0))
print("%g %g %g" % (L_2,-L_2,0))
print("%g %g %g" % (L_2,L_2+L_4,0))
print("%g %g %g" % (L_2+L_4,L_2,0))
print("%g %g %g" % (L_2+L_4,-L_2,0))
print("%g %g %g" % (L_2,-L_2-L_4,0))