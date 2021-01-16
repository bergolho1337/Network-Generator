seeds = [1562002891,1562002894,1562005513,1562005553,1562006177,1562007596,1562008172,1562008424,1562009134,1562009769,1562046115,1562013988,1562042299,1562005512,1562003066,1562009768,1562044566,1562008423,1562036996,1562020974]

for seed in seeds:
	filename = "/home/berg/Github/Network-Generator/cco-project-3d-purkinje/outputs/test_LV__seed:%d/network_info.txt" % (seed)
	file = open(filename,"r")
	for line in file:
		if (line.find("[INFO]") == -1):
			tokens = line.split()
			min_lat, max_lat, max_error, rmse, rrmse, eps2, eps5 = float(tokens[0]), float(tokens[1]), float(tokens[2]), float(tokens[3]), float(tokens[4]), float(tokens[5]), float(tokens[6])
			print("%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf" % (min_lat,max_lat,max_error,rmse,rrmse,eps2,eps5))
