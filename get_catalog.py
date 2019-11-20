import numpy as np
import matplotlib.pyplot as plt
import os, sys


# define some useful functions
def get_ra(hour, minutes, seconds):
	"""Return the RA in degrees."""
	return hour * 15. + minutes * 0.25 + (15./3600.)*seconds

def get_dec(deg, arcmin, arcsec):
	"""Return the dec in degrees."""
	return deg + (1./60.)*arcmin + (1./3600.)*arcsec


def get_angular_dist(ra1, dec1, ra2, dec2):
	"""Find the angular distance between points (ra1, dec1) and (ra2, dec2), given in degrees.
	Returns a tuple containing the angular distance in radians, degrees, and arcminutes."""
	ra1 = np.radians(ra1)
	ra2 = np.radians(ra2)
	dec1 = np.radians(dec1)
	dec2 = np.radians(dec2)
	cosx = np.cos((np.pi/2.)-dec1)*np.cos((np.pi/2.)-dec2) + np.sin((np.pi/2.)-dec1)*np.sin((np.pi/2.)-dec2)*np.cos(ra1-ra2)
	x = np.arccos(cosx)
	return x, np.degrees(x), np.degrees(x)*60.

def get_distances(ra, dec, coords):
	"""Given the coordinates (ra, dec) and a list coords of coordinate values [(ra1,dec1), (ra2,dec2), ...],
	returns a list of the distance from the point ra, dec to each point in the coordinate list. The inputs are expected to be in degrees,
	and the output is also in degrees."""
	distances = []
	for coord in coords:
		distances.append(get_angular_dist(ra, dec, coord[0], coord[1])[1])
	return np.array(distances)


# offset found from using the known magnitudes of ref stars
offset = {'B': 19.6759, 'V': 20.2690}
offset_err = {'B': 0.0001, 'V': 0.0002}
# extinction found using NASA extragalactic database
ext = {'B': 0.863, 'V': 0.653}

fmax = 50000. # max counts per pixel

filters  = ['B', 'V']
stars = ['atlas', 'hd23778', 'hd23479', 'alcyone', 'celaeno', 'merope', 'hd23511',  'electra']

# get a list of catalog filenames
cat_fnames = os.listdir(os.path.join(os.getcwd(), 'source_extractor_output'))
cat_files = {star: {'B': [], 'V': []} for star in stars}
for fname in cat_fnames:
	for star in stars:
		if star in fname:
			if 'B' in fname:
				cat_files[star]['B'].append(fname)
			elif 'V' in fname:
				cat_files[star]['V'].append(fname)




# TODO: use large aperture files to get magnitudes for the very bright/large stars, and make sure those are the ONLY values used for those stars

data = {} 
# loop through each star/frame
for star in stars:
	# loop through one filter at a time
	for filt in filters:
		# loop through each catalog file for that star and filter
		for cat_fname in cat_files[star][filt]:
			if (star in ['atlas', 'celaeno', 'alcyone', 'merope']) and ('0.1s' in cat_fname):
				# these exposures are no good (at least for now); source extractor is just detecting noise
				# TODO: see if some of the exposures are still useable
				continue
			exp_time = cat_fname.split('_')[2] # keep track of the exposure time for the image used to produce catalog file
			# we need to make sure we only take measurements for v and b mag at same exp time

			# load data from cat file
			ras, decs, fluxs, flux_errs, mags, mag_errs, background, max_flux = np.loadtxt('source_extractor_output/'+cat_fname, unpack=True)
			for i,ra in enumerate(ras):
				coords_i = (ras[i], decs[i]) # get coordinates at position i in catalog
				mag_i = mags[i] + offset[filt] - ext[filt]
				err_i = np.sqrt(mag_errs[i]**2 + offset_err[filt]**2)

				# check if the star is saturated; if so, ignore it, and if not, add its magnitude to the list
				if (max_flux[i]+background[i] > fmax) or (mag_i >= 95.) or (err_i >= 95.):
					pass
				else:
					#print cat_fname, ras[i], decs[i], mag_i, err_i
					# check if we already have these coordinates
					coords = np.array(data.keys()) # list of coordinates we already have
					if len(coords) > 0: # if this is False, the the data dictionary is empty, and we need to make it
						dists = get_distances(ras[i], decs[i], coords) # list of distances between current point and this list of coords
						loc = np.where(dists == np.min(dists)) # find point with the minimum distance
						if len(loc[0]) == 0:
							# quick and lazy bug fix; if there's a problem, just throw this data point away
							continue
						min_dist = dists[loc] # minimum distance
						guess_coords = coords[loc] # guess of the coordinates 
						#print min_dist[0], guess_coords[0][0], guess_coords[0][1]	
						if min_dist[0] < 0.001:
							data[(guess_coords[0][0], guess_coords[0][1])][filt]['mag'].append(mag_i) # add the magnitude and correct it
							data[(guess_coords[0][0], guess_coords[0][1])][filt]['err'].append(err_i)
							data[(guess_coords[0][0], guess_coords[0][1])][filt]['exp_time'].append(exp_time)
						else:
							data[coords_i] = {filt: {'mag': [], 'err': [], 'exp_time': []} for filt in filters}
							data[coords_i][filt]['mag'].append(mag_i)
							data[coords_i][filt]['err'].append(err_i)
							data[coords_i][filt]['exp_time'].append(exp_time)
					else: # if we don't have this star yet, begin a list of its measureed magnitudes
						data[coords_i] = {filt: {'mag': [], 'err': [], 'exp_time': []} for filt in filters}
						data[coords_i][filt]['mag'].append(mag_i)
						data[coords_i][filt]['err'].append(err_i)
						data[coords_i][filt]['exp_time'].append(exp_time)


coords = data.keys()
mags = {}
for coord in coords:
	# need to match up mags obtained from diff exp times - we can't use data for one filter at a given exposure time if we don't have the data for the other filter
	
	if (len(data[coord]['V']['mag']) > 0) and (len(data[coord]['B']['mag']) > 0):
		# temporary lists to hold magnitudes for the star
		tmpB = []
		tmpB_errs = []
		tmpV = []
		tmpV_errs = []
		# make sure we have at least one V and B magnitude for each exp time
		exp_times = set(data[coord]['B']['exp_time'] + data[coord]['V']['exp_time'])
		for exp_time in exp_times:
			if (exp_time in data[coord]['B']['exp_time']) and (exp_time in data[coord]['V']['exp_time']):
				# this means we have data for B and V using that exposure time
				idxsB = np.where(np.array(data[coord]['B']['exp_time']) == exp_time)
				idxsV = np.where(np.array(data[coord]['V']['exp_time']) == exp_time)
				magsB = np.array(data[coord]['B']['mag'])[idxsB]
				magsV = np.array(data[coord]['V']['mag'])[idxsV]
				errsB = np.array(data[coord]['B']['err'])[idxsB]
				errsV = np.array(data[coord]['V']['err'])[idxsV]
				for i in range(len(magsB)):
					tmpB.append(magsB[i])
					tmpB_errs.append(errsB[i])
				for i in range(len(magsV)):
					tmpV.append(magsV[i])
					tmpV_errs.append(errsV[i])
		if (len(tmpB) > 0) and (len(tmpV)>0):
			mags[coord] = {filt: {} for filt in filters}
			wtsB = 1. / np.array(tmpB_errs)**2
			wtsV = 1. / np.array(tmpV_errs)**2
			if len(tmpB) > 1:
				mags[coord]['B']['mag'] = np.average(np.array(tmpB), weights=wtsB)
				mags[coord]['B']['err'] = np.sqrt(1. / np.sum(wtsB))
			else:
				mags[coord]['B']['mag'] = tmpB[0]
				mags[coord]['B']['err'] = tmpB_errs[0]
			if len(tmpV) > 1:
				mags[coord]['V']['mag'] = np.average(np.array(tmpV), weights=wtsV)
				mags[coord]['V']['err'] = np.sqrt(1. / np.sum(wtsV))
			else:
				mags[coord]['V']['mag'] = tmpV[0]
				mags[coord]['V']['err'] = tmpV_errs[0]


# save the magnitudes in a text file
Bmags = []
Vmags = []
ra_vals = []
dec_vals = []
for coord in mags.keys():
	Bmags.append(mags[coord]['B']['mag'])
	Vmags.append(mags[coord]['V']['mag'])
	ra_vals.append(coord[0])
	dec_vals.append(coord[1])
B = np.array(Bmags)
V = np.array(Vmags)
np.savetxt('cat.txt', np.column_stack((ra_vals, dec_vals, B, V)), header="ra dec B V")


# plots for testing:

# plot every point we detected
ra, dec = zip(*data.keys())
plt.figure(figsize=(10,10))
plt.plot(ra, dec, 'x')
plt.grid()
plt.savefig('test.png')

# plot every point we're using
ra, dec = zip(*mags.keys())
plt.figure(figsize=(10,10))
plt.plot(ra, dec, 'x')
plt.grid()
plt.savefig('test2.png')

# CMD
print "plot: {} points".format(len(mags.keys()))
plt.figure(figsize=(10,10))
plt.plot(B-V, V, 'o')
plt.grid()
plt.xlabel("B-V")
plt.ylabel("V")
plt.savefig('cmd.png')
