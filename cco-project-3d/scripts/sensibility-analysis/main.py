import sys
import chaospy as cp
from cco_sections import write_cco_config_file
from co_sections import write_co_config_file
from co_activation_time_sections import write_co_activation_time_config_file
from monoalg3d_sections import write_monoalg_config_file

def main():
	
    seeds = [1562046115,1562013988,1562042299,1562005513,1562009134,1562009769,1562044567,1562008424,1562036996,1562020974,1562049673,1562023024,1562025823,1562007596,1562028813,1562005553,1562017900,1562034865,1562047507,1562011558,1562051289,1562006177,1562002894,1562002891,1562018879,1562021648,1562035947,1562008172,1562043344,1562021993]
    #seeds = [1562046115,1562013988,1562042299,1562005513,1562009134,1562009769]

    #rand_offsets = [2,3,4,5,6,7,8,9,10]
    rand_offsets = [2]

    start_radius = 0.1

    for seed in seeds:
        for rand_offset in rand_offsets:
            #write_cco_config_file(seed,rand_offset)
            #write_monoalg_config_file(seed,rand_offset,1)
            #write_co_config_file(seed,rand_offset)
            #write_monoalg_config_file(seed,rand_offset,2)
            write_co_activation_time_config_file(seed,rand_offset)
            write_monoalg_config_file(seed,rand_offset,3)

if __name__ == "__main__":
	main()


'''
# ChaosPy code

min_rand_offset = 1
max_rand_offset = 10

rand_offset = cp.Uniform(min_rand_offset,max_rand_offset)
#b = cp.Uniform(0,5)

#distribution = cp.J(a,b)
distribution = cp.J(rand_offset)

samples = distribution.sample(10,"L")

print samples.T.astype(int)
'''
