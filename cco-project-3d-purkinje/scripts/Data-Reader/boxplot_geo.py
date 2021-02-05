import numpy as np
import matplotlib.pyplot as plt

def read_all_data (foldername,basename,seeds):
    all_segments_LV = []
    all_segments_RV = []  
    all_angles_LV = []
    all_angles_RV = []  

    for seed in seeds:
        # Segment length - LV
        filename = "%s/test_LV__seed:%d_%s/segments_length.dat" % (foldername,seed,basename)
        segments = np.genfromtxt(filename)
        for line in segments:
            length = float(line)
            all_segments_LV.append(length)
        # Segment length - RV
        filename = "%s/test_RV__seed:%d_%s/segments_length.dat" % (foldername,seed,basename)
        segments = np.genfromtxt(filename)
        for line in segments:
            length = float(line)
            all_segments_RV.append(length)
        
        # Bifurcation angles - LV
        filename = "%s/test_LV__seed:%d_%s/bifurcation_angle.dat" % (foldername,seed,basename)
        angles = np.genfromtxt(filename)
        for line in angles:
            angle = float(line)
            all_angles_LV.append(angle)
        # Bifurcation angles - RV
        filename = "%s/test_RV__seed:%d_%s/bifurcation_angle.dat" % (foldername,seed,basename)
        angles = np.genfromtxt(filename)
        for line in angles:
            angle = float(line)
            all_angles_RV.append(angle)

    return all_segments_LV, all_segments_RV, all_angles_LV, all_angles_RV

# ===================================================================================================================================================
# Boxplot configuration
medianprops = dict(linestyle='-', linewidth=2.5, color='firebrick')
meanpointprops = dict(marker='D', markeredgecolor='black',
                      markerfacecolor='firebrick')

# Reference values
ref_segment_length_LV = 1.04
ref_angle_LV = 57.23
ref_segment_length_RV = 1.70
ref_angle_RV = 68.13
ref_num_segment_LV = 1349
ref_num_segment_RV = 563
num_segments_LV = 1299
num_segments_RV = 706
segment_length_ref_LV = [ref_segment_length_LV,ref_segment_length_LV,ref_segment_length_LV]
angle_ref_LV = [ref_angle_LV,ref_angle_LV,ref_angle_LV]
segment_length_ref_RV = [ref_segment_length_RV,ref_segment_length_RV,ref_segment_length_RV]
angle_ref_RV = [ref_angle_RV,ref_angle_RV,ref_angle_RV]
all_ref_num_segments_LV = [ref_num_segment_LV,ref_num_segment_LV,ref_num_segment_LV]
all_ref_num_segments_RV = [ref_num_segment_RV,ref_num_segment_RV,ref_num_segment_RV]
all_num_segments_LV = [num_segments_LV,num_segments_LV,num_segments_LV]
all_num_segments_RV = [num_segments_RV,num_segments_RV,num_segments_RV]

# Experiment 1 data
seeds = [1562002891,1562002894,1562005513,1562005553,1562008172] # Patient_Specific (Experiment 1)
folder_name = "/home/berg/Github/Network-Generator/cco-project-3d-purkinje/outputs/Patient-Specific-Experiment/Homogenous-Cloud-Min.Length"
base_name = "homogenous_min:length"
exp1_all_segments_LV, exp1_all_segments_RV, exp1_all_angles_LV, exp1_all_angles_RV = read_all_data(folder_name,base_name,seeds)

# Experiment 2 data
seeds = [1562002894,1562005513,1562005553,1562006177,1562008172] # Patient_Specific (Experiment 2)
folder_name = "/home/berg/Github/Network-Generator/cco-project-3d-purkinje/outputs/Patient-Specific-Experiment/Homogenous-Cloud-Min.Length-Min.Error"
base_name = "homogenous_min:length+error"
exp2_all_segments_LV, exp2_all_segments_RV, exp2_all_angles_LV, exp2_all_angles_RV = read_all_data(folder_name,base_name,seeds)

# Experiment 3 data
seeds = [1562002891,1562002894,1562005513,1562007596,1562008172] # Patient_Specific (Experiment 3)
folder_name = "/home/berg/Github/Network-Generator/cco-project-3d-purkinje/outputs/Patient-Specific-Experiment/Heterogenous-Cloud-with-Round-Robin/Inactives:Min.Length_Actives:Min.Error"
base_name = "min:length+error"
exp3_all_segments_LV, exp3_all_segments_RV, exp3_all_angles_LV, exp3_all_angles_RV = read_all_data(folder_name,base_name,seeds)

# LV
#segment_length_box_plot_data = [ exp1_all_segments_LV, exp2_all_segments_LV, exp3_all_segments_LV ]
#angle_box_plot_data = [ exp1_all_angles_LV, exp2_all_angles_LV, exp3_all_angles_LV ]
# RV
segment_length_box_plot_data = [ exp1_all_segments_RV, exp2_all_segments_RV, exp3_all_segments_RV ]
angle_box_plot_data = [ exp1_all_angles_RV, exp2_all_angles_RV, exp3_all_angles_RV ]

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))
# Segment length
axes[0].boxplot(segment_length_box_plot_data, medianprops=medianprops, meanprops=meanpointprops, showmeans=False, showfliers=False)
axes[0].scatter([1,2,3],segment_length_ref_RV,s=200,c="red",marker='D')
axes[0].set_title("Segment length")
axes[0].set_ylabel("Size (mm)")

# Bifrucation angle
axes[1].boxplot(angle_box_plot_data, medianprops=medianprops, meanprops=meanpointprops, showmeans=False, showfliers=False)
axes[1].scatter([1,2,3],angle_ref_RV,s=200,c="red",marker='D')
axes[1].set_title("Bifurcation angle")
axes[1].set_ylabel(r"Angle ($^{\circ}$)")

# Number of segments
axes[2].scatter([1,2,3],all_ref_num_segments_RV,s=200,c="red",marker='D')
axes[2].scatter([1,2,3],all_num_segments_RV,s=200,c="blue",marker='D')
axes[2].set_title("Number of segments")
axes[2].set_ylabel(r"#")
#axes[2].set_ylim([1250,1400]) # LV
axes[2].set_ylim([500,750]) # RV

plt.setp(axes, xticks=[y + 1 for y in range(3)],
         xticklabels=['Exp1', 'Exp2', 'Exp3'])
fig.suptitle("Right Ventricle - Geometric Results",fontsize=20)
fig.subplots_adjust(hspace=0.3)
plt.savefig("geometric_results_RV.pdf")
#plt.show()