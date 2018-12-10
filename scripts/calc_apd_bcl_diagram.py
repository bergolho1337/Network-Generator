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
            
def get_cell_id_from_filename (filename):
    aux = filename.split('-')
    aux2 = aux[1].split('.')
    cell_id = aux2[0]
    return int(cell_id)

def plot_apd_over_cable (cell_ids,even_apds,odd_apds,dir_name,output_filename):
    plt.grid()
    plt.plot(cell_ids,even_apds,label="even",c="blue",marker='s')
    plt.plot(cell_ids,odd_apds,label="odd",c="red",marker='o')
    plt.xlabel("cell id",fontsize=15)
    plt.ylabel("APD",fontsize=15)
    plt.ylim([50,250])
    plt.title("%s" % output_filename,fontsize=14)
    plt.legend(loc=0,fontsize=14)
    plt.savefig(dir_name+output_filename+".pdf")
    #plt.show()

def swap_values (a,b):
    aux = a
    a = b
    b = aux

def sort_apds_by_cell_id (cell_ids,even_apds,odd_apds):
    n = len(cell_ids)
    for i in range(n):
        for j in range(n):
            if (cell_ids[j] < cell_ids[i]):
                cell_ids[j], cell_ids[i] = cell_ids[i], cell_ids[j]
                even_apds[j],even_apds[i] = even_apds[i],even_apds[j]
                odd_apds[j],odd_apds[i] = odd_apds[i],odd_apds[j]
    return cell_ids, even_apds, odd_apds

def plot_bcl_apd (even_apds,odd_apds,bcls,cell_id):
    plt.clf()
    plt.grid()
    plt.scatter(bcls,even_apds,label="even",c="red",marker='o')
    plt.scatter(bcls,odd_apds,label="odd",c="blue",marker='o')
    plt.xlabel("BCL",fontsize=15)
    plt.ylabel("APD",fontsize=15)
    plt.title("BCL x APD (cell %d)" % cell_id,fontsize=14)
    plt.legend(loc=0,fontsize=14)
    #plt.show()
    plt.savefig("BCL-APD/bcl-apd-cell-%d.pdf" % cell_id)

def main():

    if (len(sys.argv) != 1):
	print("===============================================================")
        print("Usage:> %s" % (sys.argv[0]))
	print("===============================================================")
        sys.exit(1)
    else:

        #cell_ids = [100,400,800,1200,1600,2000]
        #cell_ids = [71,142,213,284,355,426]
        cell_ids = [0,100,200,300,400]

        for cell_id in cell_ids:
            print("[!] Working on cell number %d ..." % cell_id)

            even_apds = []
            odd_apds = []
            bcls = []
            out_filename = "BCL-APD/cell-"+str(cell_id)+".dat"
            file = open(out_filename,"w")
            
            for bcl in range(100,300,5):
                dir_name = str(bcl) + "ms"
                file_name = "sv-" + str(cell_id) + ".dat"
                
                data = np.genfromtxt(dir_name+"/"+file_name) 
                even_apd, odd_apd = calc_apd(data)

                file.write("%d %.10lf %10lf\n" % (bcl,even_apd,odd_apd))

                even_apds.append(even_apd)
                odd_apds.append(odd_apd)
                bcls.append(bcl)
            
            file.close()

            plot_bcl_apd(even_apds,odd_apds,bcls,cell_id)

if __name__ == "__main__":
    main()
