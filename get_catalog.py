import numpy as np
import matplotlib.pyplot as plt
import os, sys

# use stars with known magnitudes to calibrate our data

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

'''
B: offset = 19.675894639 +/- 0.000109913049358
V: offset = 20.2689758958 +/- 0.000218948842336
offset = known - measured
so want to add offset to mags
'''

offset = {'B': 19.6759, 'V': 20.2690}
offset_err = {'B': 0.0001, 'V': 0.0002}

fmax = 50000. # max counts per pixel

filters  = ['B', 'V']
stars = ['atlas', 'hd23778', 'hd23479', 'alcyone', 'celaeno', 'merope', 'hd23511',  'electra']
cat_fnames = os.listdir(os.path.join(os.getcwd(), 'source_extractor_output'))

cat_files = {star: {'B': [], 'V': []} for star in stars}
for fname in cat_fnames:
	for star in stars:
		if star in fname:
			if 'B' in fname:
				cat_files[star]['B'].append(fname)
			elif 'V' in fname:
				cat_files[star]['V'].append(fname)




# want data[star] = {ra: ra, dec: dec, mag_V: magv, mag_B: magb}
#coords = {star: {'B': [], 'V': []} for star in stars}

#data = {star: {} for star in stars}

data = {}

for star in stars[1:2]: # for now just testing one star
	for filt in filters:
		for cat_fname in cat_files[star][filt]:
			exp_time = cat_fname.split('_')[2] # keep track of the exposure time for the image used to produce catalog file
			# we need to make sure we only take measurements for v and b mag at same exp time

			# load data from cat file
			ras, decs, fluxs, flux_errs, mags, mag_errs, background, max_flux = np.loadtxt('source_extractor_output/'+cat_fname, unpack=True)
			for i,ra in enumerate(ras):
				coords_i = (ras[i], decs[i]) # get coordinates at position i in catalog
				mag_i = mags[i] + offset[filt]
				err_i = np.sqrt(mag_errs[i]**2 + offset_err[filt]**2)

				# check if the star is saturated; if so, ignore it, and if not, add its magnitude to the list
				if (max_flux[i]+background[i] > fmax) or (mag_i >= 95.) or (err_i >= 95.):
					pass
				else:
					print cat_fname, ras[i], decs[i], mag_i, err_i
					# check if we already have these coordinates
					coords = np.array(data.keys()) # list of coordinates we already have
					if len(coords) > 0:
						dists = get_distances(ras[i], decs[i], coords) # list of distances between current point and this list of coords
						loc = np.where(dists == np.min(dists)) # find point with the minimum distance
						print loc
						min_dist = dists[loc] # minimum distance
						guess_coords = coords[loc] # guess of the coordinates 
						print min_dist[0], guess_coords[0][0], guess_coords[0][1]	
						if min_dist[0] < 0.001:
							data[(guess_coords[0][0], guess_coords[0][1])][filt]['mag'].append(mag_i) # add the magnitude and correct it
							data[(guess_coords[0][0], guess_coords[0][1])][filt]['err'].append(err_i)
							data[(guess_coords[0][0], guess_coords[0][1])][filt]['exp_time'].append(exp_time)
						else:
							data[coords_i] = {filt: {'mag': [], 'err': [], 'exp_time': []} for filt in filters}
							data[coords_i][filt]['mag'].append(mag_i)
							data[coords_i][filt]['err'].append(err_i)
							cata[coords_i][filt]['exp_time'].append(exp_time)
					else: # if we don't have this star yet, begin a list of its measureed magnitudes
						data[coords_i] = {filt: {'mag': [], 'err': [], 'exp_time': []} for filt in filters}
						data[coords_i][filt]['mag'].append(mag_i)
						data[coords_i][filt]['err'].append(err_i)
						cata[coords_i][filt]['exp_time'].append(exp_time)


#print data.keys()	
coords = data.keys()
for coord in coords:
	# need to match up mags obtained from diff catalogs
	# ex - if we get v mags from cat 0, 1, and 3 at 10s exp time, and get B mags from cat 1, 3, 4 at 10 s exp,
	# then we want to only keep mags from cat 1 and 3?
	# or, that doesn't really matter - just make sure same exp time (then automatically same aperture)
	
	if (len(data[coord]['V']['mag']) > 0) and (len(data[coord]['B']['mag']) > 0):
		print coord, len(data[coord]['V']['mag']), len(data[coord]['B']['mag'])

		
ra, dec = zip(*data.keys())
plt.figure(figsize=(10,10))
plt.scatter(ra, dec)
plt.savefig('test.png')

