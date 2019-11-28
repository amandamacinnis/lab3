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

'''
B: offset = 19.061 +/- 0.071
V: offset = 19.882 +/- 0.050
'''

# offset found from using the known magnitudes of ref stars
#offset = {'B': 19.6701, 'V': 20.4728}
#offset_err = {'B': 0.0001, 'V': 0.0002}
#offset = {'B': 20.648, 'V': 21.219}
#offset_err = {'B': 0.034, 'V': 0.014}
offset = {'B': 19.061, 'V': 19.882}
offset_err = {'B': 0.071, 'V': 0.050}
# extinction found using NASA extragalactic database
ext = {'B': 0.863, 'V': 0.653}

fmax = 60000. # max counts per pixel

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





data = {} 
# loop through each star/frame
for star in stars:
	# loop through one filter at a time
	for filt in filters:
		# loop through each catalog file for that star and filter
		for cat_fname in cat_files[star][filt]:
			if (cat_fname=='atlas_V_0.1s_01.cat') or (cat_fname=='atlas_V_0.1s_03.cat'):# or ('large' in cat_fname):
				continue # something wrong with these specific files

			# load data from cat file
			ras, decs, fluxs, flux_errs, mags, mag_errs, background, max_flux = np.loadtxt('source_extractor_output/'+cat_fname, unpack=True)
			print "loading {} ... {} entries in catalog".format(cat_fname,len(ras))

			for i,ra in enumerate(ras):
				coords_i = (ras[i], decs[i]) # get coordinates at position i in catalog
				mag_i = mags[i] + offset[filt] - ext[filt] # get the magnitude and correct it
				err_i = np.sqrt(mag_errs[i]**2 + offset_err[filt]**2)

				# check if the star is saturated; if so, ignore it, and if not, add its magnitude to the list
				if (max_flux[i]+background[i] > fmax) or (mag_i >= 95.) or (err_i >= 95.):
					pass
				else:
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
						if min_dist[0] < 0.01: # if it's within 30 arcseconds of a coordinate we've already found
							data[(guess_coords[0][0], guess_coords[0][1])][filt]['mag'].append(mag_i) # add the magnitude 
							data[(guess_coords[0][0], guess_coords[0][1])][filt]['err'].append(err_i)
						else:   # if we don't have this star yet, begin a list of its measureed magnitudes
							data[coords_i] = {filt: {'mag': [], 'err': []} for filt in filters}
							data[coords_i][filt]['mag'].append(mag_i)
							data[coords_i][filt]['err'].append(err_i)
					else: # begin the list of coordinates
						data[coords_i] = {filt: {'mag': [], 'err': []} for filt in filters}
						data[coords_i][filt]['mag'].append(mag_i)
						data[coords_i][filt]['err'].append(err_i)


coords = data.keys()

max_err = 1.  # maximum error for magnitude to include it as a data point

Bmags = []
Berrs = []
Vmags = []
Verrs = []
BV = []
BV_err = []
ra_vals = []
dec_vals = []

for coord in coords:
	# need to match up mags obtained from diff exp times - we can't use data for one filter at a given exposure time if we don't have the data for the other filter
	
	if (len(data[coord]['V']['mag']) > 0) and (len(data[coord]['B']['mag']) > 0):
		# temporary lists to hold magnitudes for the star
		tmpB = data[coord]['B']['mag'][:]
		tmpB_errs = data[coord]['B']['err'][:]
		tmpV = data[coord]['V']['mag'][:]
		tmpV_errs = data[coord]['V']['err'][:]
		mags[coord] = {filt: {} for filt in filters}
		wtsB = 1. / np.array(tmpB_errs)**2
		wtsV = 1. / np.array(tmpV_errs)**2
		mag_B = None
		mag_V = None
		# only save the data point if the errors aren't huge (this was checked manually- only very dim stars will be removed)
		if len(tmpB) > 1:
			mag_B = np.average(np.array(tmpB), weights=wtsB)
			mag_B_err = np.sqrt(1. / np.sum(wtsB))
		else:
			mag_B = tmpB[0]
			mag_B_err = tmpB_errs[0]
		if len(tmpV) > 1:
			mag_V = np.average(np.array(tmpV), weights=wtsV)
			mag_V_err = np.sqrt(1. / np.sum(wtsV))
		else:
			mag_V = tmpV[0]
			mag_V_err = tmpV_errs[0]
		if (mag_B_err <= max_err) and (mag_V_err <= max_err):
			Bmags.append(mag_B)
			Berrs.append(mag_B_err)
			Vmags.append(mag_V)
			Verrs.append(mag_V_err)
			BV.append(mag_B - mag_V)
			BV_err.append(np.sqrt(mag_B_err**2 + mag_V_err**2))
			ra_vals.append(coord[0])
			dec_vals.append(coord[1])


# save the magnitudes in a text file
B = np.array(Bmags)
Berrs = np.array(Berrs)
V = np.array(Vmags)
Verrs = np.array(Verrs)
BV = np.array(BV)
BV_err = np.array(BV_err)
np.savetxt('cat_v7.txt', np.column_stack((ra_vals, dec_vals, B, Berrs, V, Verrs, BV, BV_err)), header="ra, dec, B, B err, V, V err, B-V, B-V err")

print "{} points out of {} detected stars".format(len(V), len(data.keys()))

# plots for testing:

# plot every point we detected
ra, dec = zip(*data.keys())
plt.figure(figsize=(10,10))
plt.plot(ra, dec, 'x')
plt.grid()
plt.savefig('test_v4.png')

# plot every point we're using
ra, dec = zip(*mags.keys())
plt.figure(figsize=(10,10))
plt.plot(ra, dec, 'x')
plt.grid()
plt.savefig('test2_v4.png')

# CMD
plt.figure(figsize=(10,10))
plt.plot(B-V, V, 'o')
plt.grid()
plt.gca().invert_yaxis()
plt.xlabel("B-V")
plt.ylabel("V")
plt.savefig('cmd_v4.png')

plt.figure(figsize=(10,10))
plt.errorbar(BV, V, xerr=BV_err, yerr=Verrs, fmt='.', lw=0.9, capsize=2, alpha=0.85)
plt.grid()
plt.gca().invert_yaxis()
plt.xlabel("B-V")
plt.ylabel("V")
plt.savefig('cmd_err_v4.png')
