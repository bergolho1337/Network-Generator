import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def plot_data (data,xname,yname):
    plt.grid()
    plt.plot(data[:,0],data[:,1],label=yname,c="black",linewidth=3.0)
    plt.xlabel(xname,fontsize=15)
    plt.ylabel(yname,fontsize=15)
    plt.title("Data",fontsize=14)
    #plt.xlim([2000,3000])
    plt.legend(loc=2,fontsize=14)
    plt.show()
    #plt.savefig("ap.pdf")

def main():

    if (len(sys.argv) != 3):
        print("==========================================================")
        print("Usage:> %s <x-axis-name> <y-axis-name>" % (sys.argv[0]))
        print("==========================================================")
        sys.exit(1)
    else:
        input_file_1 = "../output/velocity_squid_axon.txt"
        input_file_2 = "../output/velocity_lobster_axon.txt"
        input_file_3 = "../output/velocity_crab_axon.txt"
        input_file_4 = "../output/velocity_earthworm_axon.txt"
        input_file_5 = "../output/velocity_marine_worm_axon.txt"
        input_file_6 = "../output/velocity_mammalian_cardiac.txt"
        input_file_7 = "../output/velocity_barnacle_muscle.txt"
        xname = sys.argv[1]
        yname = sys.argv[2]

        data_1 = np.genfromtxt(input_file_1)
        data_2 = np.genfromtxt(input_file_2)
        data_3 = np.genfromtxt(input_file_3)
        data_4 = np.genfromtxt(input_file_4)
        data_5 = np.genfromtxt(input_file_5)
        data_6 = np.genfromtxt(input_file_6)
        data_7 = np.genfromtxt(input_file_7)

        plt.plot(data_1[:,0],data_1[:,1],label="Squid",linewidth=3.0)
        plt.plot(data_2[:,0],data_2[:,1],label="Lobster",linewidth=3.0)
        plt.plot(data_3[:,0],data_3[:,1],label="Crab",linewidth=3.0)
        plt.plot(data_4[:,0],data_4[:,1],label="Earthworm",linewidth=3.0)
        plt.plot(data_5[:,0],data_5[:,1],label="Marine worm",linewidth=3.0)
        plt.plot(data_6[:,0],data_6[:,1],label="Mammalian cardiac",linewidth=3.0)
        plt.plot(data_7[:,0],data_7[:,1],label="Barnacle muscle",linewidth=3.0)

        plt.grid()
        plt.xlim([0,600])
        plt.ylim([0,30])
        plt.xlabel(xname,fontsize=15)
        plt.ylabel(yname,fontsize=15)
        plt.title("Data",fontsize=14)
        plt.legend(loc=0,fontsize=8)
        #plt.show()
        plt.savefig("../figures/propagation_velocity_excitable_cells.pdf")
        #plot_data(data,xaxis_name,yaxis_name)



if __name__ == "__main__":
    main()
