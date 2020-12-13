# ===========================================================================================================
# This script calibrates the diameter of the Purkinje cells in order to adjust to the reference value
# ===========================================================================================================

import numpy as np

Gi = 7.9
Cf = 3.4
tauf = 0.1

def calc_diameter (v):
    return (v*v*4.0*Cf*tauf) / (Gi) * 100.0

def calc_velocity (d):
    return ((Gi*d)/(4.0*Cf*tauf) )**(0.5) * 0.1

def calc_proportion (d,ref_lat,aprox_lat):
    return d*aprox_lat/ref_lat

ref_velocity = 1.9
ref_diameter = calc_diameter(ref_velocity)
print("Reference diameter = %g" % (ref_diameter))

ref_lat = 67.25
aprox_lat = 57.59
new_d = calc_proportion(ref_diameter,ref_lat,aprox_lat)
print("New diameter = %g" % (new_d))

new_v = calc_velocity(new_d)
print("New velocity = %g" % (new_v))

print("CV(%g) = %g" % (100,calc_velocity(100)))
print("CV(%g) = %g" % (200,calc_velocity(200)))
