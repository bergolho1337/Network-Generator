import sys
import numpy as np

def print_data (data):
    n, m = data.shape
    for i in range(n):
        for j in range(m):
            print("%.2lf " % data[i][j]),
        print("")

def get_closest_point (x1, y1, z1, data):
	closest_id = 0
	closest_dist = sys.float_info.max
	for i in range(len(data)):
		x2, y2, z2 = data[i][0], data[i][1], data[i][2]
		dist = np.sqrt( (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2 )
		if (dist < closest_dist):
			closest_dist = dist
			closest_id = i
	return closest_id, closest_dist
		

def calculate_errors (ref, aprox):
	a = []
	b = []
	errors = []
	errors_2 = []
	ids = []
	x, y, z, lat = ref[:,0], ref[:,1], ref[:,2], ref[:,3]
	for i in range(len(ref)):
		x1, y1, z1, lat1 = x[i], y[i], z[i], lat[i]
		a.append(lat1)
		k, dist = get_closest_point(x1,y1,z1,aprox)
		x2, y2, z2, lat2 = aprox[k][0], aprox[k][1], aprox[k][2], aprox[k][3]
		error = abs(lat1-lat2)
		errors.append(error)
		errors_2.append(lat1-lat2)
		b.append(lat2)
		ids.append(k)
	max_error, rmse, rrmse = calculate_statistics(a,b,errors)
	return max_error, rmse, rrmse*100.0, ids
        
def calculate_statistics (ref,aprox,errors):
    ref_min_lat, ref_max_lat = np.min(ref), np.max(ref)
    aprox_min_lat, aprox_max_lat = np.min(aprox), np.max(aprox)
    sum_num = 0.0
    sum_den = 0.0
    max_error = sys.float_info.min
    for i in range(len(errors)):
        if (errors[i] > max_error):
            max_error = errors[i]
        sum_num = sum_num + errors[i]**2
        sum_den = sum_den + ref[i]**2
    l2_norm = np.sqrt(sum_den)
    rmse = np.sqrt(sum_num/len(ref))
    rrmse = np.sqrt(sum_num/sum_den)
    print("Min.LAT Reference = %g" % (ref_min_lat))
    print("Max.LAT Reference = %g" % (ref_max_lat))
    print("Min.LAT Aproximation = %g" % (aprox_min_lat))
    print("Max.LAT Aproximation = %g" % (aprox_max_lat))
    return max_error, rmse, rrmse

def get_active_pmjs (term, ids):
	pmjs = []
	for i in range(len(ids)):
		pmjs.append(term[ids[i]])
	return np.array(pmjs)
	
analitical_data = np.genfromtxt("data/Model_Zero/analitical_ELIZABETH_reference.dat")
monodomain_reference_pk_data = np.genfromtxt("data/Model_Zero/monodomain_purkinje_only_ELIZABETH_reference.dat")
monodomain_aproximation_pk_data = np.genfromtxt("data/Model_Zero/monodomain_purkinje_only_data_ELIZABETH_best_cco.dat")

# Analitical x Monodomain (Reference)
max_error, rmse, rrmse, ids = calculate_errors(analitical_data,monodomain_reference_pk_data)
print("Max. Error = %g ms" % max_error)
print("RMSE = %g ms" % rmse)
print("RRMSE = %g %%" % rrmse)

# Monodomain x Monodomain (Reference x Aproximation)
monodomain_reference_pk_data_pmjs = get_active_pmjs(monodomain_reference_pk_data,ids)
max_error, rmse, rrmse, ids = calculate_errors(monodomain_reference_pk_data_pmjs,monodomain_aproximation_pk_data)
print("Max. Error = %g ms" % max_error)
print("RMSE = %g ms" % rmse)
print("RRMSE = %g %%" % rrmse)
