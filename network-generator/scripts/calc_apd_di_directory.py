import os
import sys
import numpy as np
from glob import glob
import matplotlib.pyplot as plt

# Build a tuple where: T{i} = (t_i,dvdt_i)
def compute_derivatives (t,vm):
    n = len(vm)-1
    dvdt = {}
    for i in range(n):
        dvdt[t[i]] = vm[i+1]-vm[i] 
    return dvdt

def get_cell_id_from_filename (filename):
    aux = filename.split('-')
    aux2 = aux[1].split('.')
    cell_id = aux2[0]
    return int(cell_id)

def calc_apd (sv):
    # Copy the state vector data to local variables
    t = sv[:,0]
    vm = sv[:,1]
    
    # Calculate APD resting value
    min_value = np.amin(vm)
    max_value = np.amax(vm)
    APD_90_value = abs(max_value - min_value)*0.1 + min_value

    # Calculate the derivatives for each interval and sort then in decreasing order
    dvdt = compute_derivatives(t,vm)
    max_dvdt = sorted(dvdt.iteritems(), key=lambda (k,v): (v,k), reverse=True)

    #for i in range(10):
    #    print max_dvdt[i]

    # Get the peak times of each action potential
    #default: peak_times = max_dvdt[0:2]
    peak_times = max_dvdt[0:3]

    # Sort by time first
    peak_times.sort()

    # Sort by derivative next
    #for i in range(0,3):
    #	for j in range(i+1,3):
    #		if ( abs(peak_times[i][1]-peak_times[j][1]) > 10.0 ):
    #			aux = peak_times[i]
    #			peak_times[i] = peak_times[j]
    #			peak_times[j] = aux
    peak_times = peak_times[0:3]

    #print("Peak times")
    #print peak_times
    
    apds = []
    for i in range(len(peak_times)-1):
        start_index = 0
        end_index = 0
        t1 = peak_times[i][0]
        t2 = peak_times[i][0]

        # Find the time where this stimulus occur in the data
        for j in range(len(t)):
            diff = abs(t[j] - t1)
            if (diff < 1.0e-08):
                start_index = j
        # Until the transmembrane potential reaches the APD_90_value we iterate over the data
        # When the difference between the current potential and the APD_value is minimum and when
        # we do not pass the time of the next stimulus we can save the end time of the APD
        min_diff = max_value
        min_t = t1
        for j in range(start_index,len(t)):
            # If the current time pass the start time of the next stimulus we break the loop
            if (i < len(peak_times)-1 and t[j]+5.0 > peak_times[i+1][0]):
                #min_t = t[j]
                #end_index = j
                break
            # Until then we find the minimum difference between the APD value and the potential
            diff = abs(APD_90_value - vm[j])
            # Save the time when the difference is minimum
            if (diff < min_diff):
                min_diff = diff
                end_index = j
                min_t = t[j]
        t2 = min_t
        
        apd = t2-t1
        apds.append(apd)
        
	'''
	print("------------------------------")
        print("APD %d" % i)
        print("t1 = %.10lf" % t1)
        print("t2 = %.10lf" % t2)
        print("start_index = %d" % start_index)
        print("end_index = %d" % end_index)
        print("diff = %.10lf" % min_diff)
        print("value = %.10lf" % (apd))
        print("------------------------------")
	'''
        
    return apds[0], apds[1]

def calc_di (even_apd,odd_apd,bcl):
	even_di = bcl - even_apd
	odd_di = bcl - odd_apd

	'''
	print("------------------------------")
        print("APD 1 = %.10lf" % even_apd)
	print("DI 1 = %.10lf" % even_di)
	print("BCL = %.10lf" % bcl)
        print("------------------------------")
	print("------------------------------")
        print("APD 2 = %.10lf" % odd_apd)
	print("DI 2 = %.10lf" % odd_di)
	print("BCL = %.10lf" % bcl)
        print("------------------------------")
	'''

	return even_di, odd_di

def calc_apd_di (sv,bcl):
	even_apd, odd_apd = calc_apd(sv)
	even_di, odd_di = calc_di(even_apd,odd_apd,bcl)
	
	return even_apd, odd_apd, even_di, odd_di

def sort_apds_by_cell_id (cell_ids,even_apds,odd_apds,even_dis,odd_dis):
    n = len(cell_ids)
    for i in range(n):
        for j in range(n):
            if (cell_ids[j] < cell_ids[i]):
                cell_ids[j], cell_ids[i] = cell_ids[i], cell_ids[j]
                even_apds[j],even_apds[i] = even_apds[i],even_apds[j]
                odd_apds[j],odd_apds[i] = odd_apds[i],odd_apds[j]
		even_dis[j],even_dis[i] = even_dis[i],even_dis[j]
                odd_dis[j],odd_dis[i] = odd_dis[i],odd_dis[j]
    return cell_ids, even_apds, odd_apds, even_dis, odd_dis

def plot_apd_over_cable (cell_ids,even_apds,odd_apds,dir_name,bcl):
    plt.grid()
    plt.plot(cell_ids,even_apds,label="even",c="blue",marker='s')
    plt.plot(cell_ids,odd_apds,label="odd",c="red",marker='o')
    plt.xlabel("cell id",fontsize=15)
    plt.ylabel("APD",fontsize=15)
    plt.ylim([50,250])
    plt.title("BCL = %d" % bcl,fontsize=14)
    plt.legend(loc=0,fontsize=14)
    plt.savefig(dir_name + "/cable_APD_BCL%dms" % bcl + ".pdf")
    #plt.show()

def main():

    if (len(sys.argv) != 3):
	print("******************************************************************************************")        
	print("Usage:> %s <directory_name> <BCL>" % (sys.argv[0]))
	print("\t<BCL> = Basic cycle length (period)")
	print("\t<directory_name> = Path to the directory for the data")
	print("******************************************************************************************")
        sys.exit(1)
    else:
	dir_name = sys.argv[1]        
	bcl = float(sys.argv[2])	

	# Capture the data from 
        #data = np.genfromtxt(input_filename)

	# List to store the APD of each measured cell over the cable
        # Even = First AP
        # Odd = Second AP
        even_apds = []
        odd_apds = []
	even_dis = []
	odd_dis = []
        cell_ids = []

        # List of files contaning the measured AP of the marked cells
        files = [f for f in glob(dir_name + '/*.dat')]
		
	# For each file calculate the Even and Odd APD and store then in a list
        for f in files:		
		cell_id = get_cell_id_from_filename(f)
            	data = np.genfromtxt(f)

            	print("**********************************************************")
            	print("Cell id = %s" % cell_id)

		even_apd, odd_apd = calc_apd(data)
		even_di, odd_di = calc_di(even_apd,odd_apd,bcl)
		
		cell_ids.append(cell_id)
            	even_apds.append(even_apd)     
            	odd_apds.append(odd_apd)
		even_dis.append(even_di)     
            	odd_dis.append(odd_di)
	
	# Sort the APD by ascending order of the cell id
	cell_ids, even_apds, odd_apds, even_dis, odd_dis = sort_apds_by_cell_id(cell_ids,even_apds,odd_apds,even_dis,odd_dis)
	
	# Plot the APDs over the cable
	plot_apd_over_cable(cell_ids,even_apds,odd_apds,dir_name,bcl)
            
if __name__ == "__main__":
    main()
