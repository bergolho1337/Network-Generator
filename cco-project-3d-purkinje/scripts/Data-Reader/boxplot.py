import numpy as np
import matplotlib.pyplot as plt

# Boxplot configuration
medianprops = dict(linestyle='-', linewidth=2.5, color='firebrick')
meanpointprops = dict(marker='D', markeredgecolor='black',
                      markerfacecolor='firebrick')

# Data column indexes
MINLAT = 0
MAXLAT = 1
MAXERROR = 2
RMSE = 3
RRMSE = 4
EPS2 = 5
EPS5 = 6

# Reference values
#ref_min_lat = 16.23 # LV
#ref_max_lat = 50.37 # LV
ref_min_lat = 18.58 # RV
ref_max_lat = 67.25 # RV

# Reading data
data_exp1 = np.genfromtxt("data/electric_data_RV_exp1.dat")
data_exp2 = np.genfromtxt("data/electric_data_RV_exp2.dat")
data_exp3 = np.genfromtxt("data/electric_data_RV_exp3.dat")

# Min. LAT
min_lat_box_plot_data = [ data_exp1[:,0], data_exp2[:,0], data_exp3[:,0] ]
min_lat_ref_value = [ref_min_lat, ref_min_lat, ref_min_lat]
# Max. LAT
max_lat_box_plot_data = [ data_exp1[:,1], data_exp2[:,1], data_exp3[:,1] ]
max_lat_ref_value = [ref_max_lat, ref_max_lat, ref_max_lat]
# Max. Error
max_error_box_plot_data = [ data_exp1[:,2], data_exp2[:,2], data_exp3[:,2] ]
# RMSE
rmse_box_plot_data = [ data_exp1[:,3], data_exp2[:,3], data_exp3[:,3] ]
# RRMSE
rrmse_box_plot_data = [ data_exp1[:,4], data_exp2[:,4], data_exp3[:,4] ]
# eps < 2ms
eps2_box_plot_data = [ data_exp1[:,5], data_exp2[:,5], data_exp3[:,5] ]
# eps < 5ms
eps5_box_plot_data = [ data_exp1[:,6], data_exp2[:,6], data_exp3[:,6] ]

fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(12, 12))
# Min. LAT
axes[0, 0].boxplot(min_lat_box_plot_data, medianprops=medianprops, meanprops=meanpointprops, showmeans=False)
axes[0, 0].scatter([1,2,3],min_lat_ref_value,s=200,c="red",marker='D')
axes[0, 0].set_title("Min.LAT")
axes[0, 0].set_ylabel("LAT (ms)")
#axes[0, 0].set_ylim([14,17])

# Max. LAT
axes[0, 1].boxplot(max_lat_box_plot_data, medianprops=medianprops, meanprops=meanpointprops, showmeans=False)
axes[0, 1].scatter([1,2,3],max_lat_ref_value,s=200,c="red",marker='D')
axes[0, 1].set_title("Max.LAT")
axes[0, 1].set_ylabel("LAT (ms)")
#axes[0, 1].set_ylim([40,60])

# Max. Error
axes[0, 2].boxplot(max_error_box_plot_data, medianprops=medianprops, meanprops=meanpointprops, showmeans=False)
axes[0, 2].set_title("Max.Error")
axes[0, 2].set_ylabel("LAT (ms)")

# RMSE
axes[1, 0].boxplot(rmse_box_plot_data, medianprops=medianprops, meanprops=meanpointprops, showmeans=False)
axes[1, 0].set_title("RMSE")
axes[1, 0].set_ylabel("LAT (ms)")

# RRMSE
axes[1, 1].boxplot(rrmse_box_plot_data, medianprops=medianprops, meanprops=meanpointprops, showmeans=False)
axes[1, 1].set_title("RRMSE")
axes[1, 1].set_ylabel("%")

# eps < 2ms
axes[1, 2].boxplot(eps2_box_plot_data, medianprops=medianprops, meanprops=meanpointprops, showmeans=False)
axes[1, 2].set_title(r"$\epsilon$ < 2ms")
axes[1, 2].set_ylabel("%")

# eps < 5ms
axes[2, 0].boxplot(eps5_box_plot_data, medianprops=medianprops, meanprops=meanpointprops, showmeans=False)
axes[2, 0].set_title(r"$\epsilon$ < 5ms")
axes[2, 0].set_ylabel("%")

axes[2, 1].set_visible(False)
axes[2, 2].set_visible(False)

plt.setp(axes, xticks=[y + 1 for y in range(3)],
         xticklabels=['Exp1', 'Exp2', 'Exp3'])
fig.suptitle("Right Ventricle - Electric Results",fontsize=20)
fig.subplots_adjust(hspace=0.3)
plt.savefig("electric_results_RV.pdf")
#plt.show()