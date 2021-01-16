import numpy as np

seeds = [1562002891,1562002894,1562005513,1562005553,1562006177,1562007596,1562008172,1562008424,1562009134,1562009769,1562046115,1562013988,1562042299,1562005512,1562003066,1562009768,1562044566,1562008423,1562036996,1562020974]
all_segments = []
all_angles = []
for seed in seeds:
	# Segment length
	filename = "/home/berg/Github/Network-Generator/cco-project-3d-purkinje/outputs/test_LV__seed:%d/segments_length.dat" % (seed)
	segments = np.genfromtxt(filename)
	for line in segments:
		length = float(line)
		all_segments.append(length)
	mean_length = np.mean(segments)
	std_length = np.std(segments)
	
	# Bifurcation angles
	filename = "/home/berg/Github/Network-Generator/cco-project-3d-purkinje/outputs/test_LV__seed:%d/bifurcation_angle.dat" % (seed)
	angles = np.genfromtxt(filename)
	for line in angles:
		angle = float(line)
		all_angles.append(angle)
	mean_angle = np.mean(angles)
	std_angle = np.std(angles)
	print("%.2lf +/- %.2lf\t%.2lf +/- %.2lf\t1299" % (mean_length,std_length,mean_angle,std_angle))

all_mean_segment_length	= np.mean(all_segments)
all_std_segment_length = np.std(all_segments)
all_mean_bifurcation_angle = np.mean(all_angles)
all_std_bifurcation_angle = np.std(all_angles)

#print("%.2lf +/- %.2lf" % (all_mean_segment_length,all_std_segment_length))
#print("%.2lf +/- %.2lf" % (all_mean_bifurcation_angle,all_std_bifurcation_angle))

