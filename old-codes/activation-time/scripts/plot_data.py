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

    if (len(sys.argv) != 4):
        print("==========================================================")
        print("Usage:> %s <datafile> <x-axis-name> <y-axis-name>" % (sys.argv[0]))
        print("==========================================================")
        sys.exit(1)
    else:
        input_file = sys.argv[1]
        xaxis_name = sys.argv[2]
        yaxis_name = sys.argv[3]

        data = np.genfromtxt(input_file)

        plot_data(data,xaxis_name,yaxis_name)


if __name__ == "__main__":
    main()
