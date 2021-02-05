import numpy as np

#seeds = [1562002891,1562002894,1562005513,1562005553,1562006177,1562007596,1562008172,1562008424,1562009134,1562009769,1562046115,1562013988,1562042299,1562005512,1562003066,1562009768,1562044566,1562008423,1562036996,1562020974]
#seeds = [1562002891,1562002894,1562005513,1562005553,1562006177,1562007596,1562008172,1562008424,1562009134,1562009769]
#seeds = [1562002891,1562002894,1562005513,1562005553,1562006177,1562008172]
#seeds = [1562002891,1562002894,1562005513,1562005553,1562008172] # Patient_Specific (Experiment 1)
seeds = [1562002894,1562005513,1562005553,1562006177,1562008172] # Patient_Specific (Experiment 2)
#seeds = [1562002891,1562002894,1562005513,1562007596,1562008172] # Patient_Specific (Experiment 3)
folder_name = "/home/berg/Github/Network-Generator/cco-project-3d-purkinje/outputs/Patient-Specific-Experiment/Homogenous-Cloud-Min.Length-Min.Error"
all_segments = []
all_angles = []
num_segments = 706
for seed in seeds:
	# Segment length
	filename = "%s/test_LV__seed:%d_homogenous_min:length+error/segments_length.dat" % (folder_name,seed)
	segments = np.genfromtxt(filename)
	for line in segments:
		length = float(line)
		all_segments.append(length)
	mean_length = np.mean(segments)
	std_length = np.std(segments)
	
	# Bifurcation angles
	filename = "%s/test_LV__seed:%d_homogenous_min:length+error/bifurcation_angle.dat" % (folder_name,seed)
	angles = np.genfromtxt(filename)
	for line in angles:
		angle = float(line)
		all_angles.append(angle)
	mean_angle = np.mean(angles)
	std_angle = np.std(angles)
	print("%.2lf +/- %.2lf\t%.2lf +/- %.2lf\t%d" % (mean_length,std_length,mean_angle,std_angle,num_segments))
	#print("Min.Length+Error & %d & %.2lf $\\pm$ %.2lf & %.2lf $\\pm$ %.2lf & 1299 \\\\" % (seed,mean_length,std_length,mean_angle,std_angle))

all_mean_segment_length	= np.mean(all_segments)
all_std_segment_length = np.std(all_segments)
all_mean_bifurcation_angle = np.mean(all_angles)
all_std_bifurcation_angle = np.std(all_angles)

print("%.2lf $\\pm$ %.2lf" % (all_mean_segment_length,all_std_segment_length))
print("%.2lf $\\pm$ %.2lf" % (all_mean_bifurcation_angle,all_std_bifurcation_angle))

