import numpy as np
import os

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



filters  = ['B', 'V']
cat_fnames = os.listdir(os.path.join(os.getcwd(), 'source_extractor_output'))

cat_files = {'B': [], 'V': []}
for fname in cat_fnames:
	if ('B' in fname):
		cat_files['B'].append(fname)
	elif ('V' in fname):
		cat_files['V'].append(fname)


# note the magnitudes are found from the https://www.aavso.org/ finding charts
# create a dict of reference stars to use for calibration
ref_stars = {}
ref_stars['merope'] = {'ra': get_ra(3., 46., 19.55), 
	'dec': get_dec(23., 56., 53.3), 
	'mag_V': 4.180, 'mag_V_err': 0.100, 'mag_B': 4.120, 'mag_B_err': 0.141,   # website lists source as 'BSC'
	'cat_files': [('merope', 0.1, False), ('merope', 1, True)]} # specify catalog files to be used 

ref_stars['hd23778_ref1'] = {'ra': get_ra(3., 48., 30.07),
	'dec': get_dec(24, 20., 44.5),
	'mag_V': 6.964, 'mag_V_err': 0.100, 'mag_B': 7.023, 'mag_B_err': 0.141, # WBVR
	'cat_files': [('hd23778', 1, False)]} 

ref_stars['hd23778_ref2'] = {'ra': get_ra(3., 49., 25.98),
	'dec': get_dec(24., 14., 52.5),
	'mag_V': 7.993, 'mag_V_err': 0.100, 'mag_B': 8.095, 'mag_B_err': 0.141, # WBVR
	'cat_files': [('hd23778', 1, False), ('atlas', 1, False)]}



# this is horrible, but it gets the job done ... for each star, get only the catalog filenames we want to use
for star in ref_stars.keys(): # loop through stars
	ref_stars[star]['fnames'] = {}
	for filt in filters: # loop through filters
		ref_stars[star]['fnames'][filt] = [] # list of catalog filenames for each star in each filter
		for cat_info in ref_stars[star]['cat_files']: # for each catalog file we've specified
			star_name, exp_time, use_large_ap = cat_info # get exp time and aperture for catalog
			exp_str = '_{}s'.format(exp_time) # used to identify correct filenames
			for fname in cat_files[filt]: # loop through each filename and add correct ones to list
				if (star_name in fname) and (exp_str in fname): 
					if use_large_ap and ('large' in fname):
						ref_stars[star]['fnames'][filt].append(fname)
					elif (not use_large_ap) and ('large' not in fname):
						ref_stars[star]['fnames'][filt].append(fname)


# look at the measured magnitude of each
for star in ref_stars.keys():  # loop through the stars
	ra = ref_stars[star]['ra']
	dec = ref_stars[star]['dec']
	# loop through each catalog for the given  filter
	for filt in filters: # loop through filters
		measured_mags = []
		measured_mag_errs = []
		for fname in ref_stars[star]['fnames'][filt]:  # loop through catalogs
			# load the catalog
			ras, decs, fluxs, flux_errs, mags, mag_errs, background, max_flux = np.loadtxt('source_extractor_output/'+fname, unpack=True)
			# get a list of angular distances (in arcmin) between where we think the star is, and the objects detected by source extractor
 			_, _, distances = get_angular_dist(ra, dec, ras, decs)
			# make sure we have some reasonable star matching the coordinates we gave
			if np.min(distances) <= 1.:
				# the star we want will (hopefully) have the minimum angular distance 
				loc = np.where(distances == np.min(distances))
				print star, fname, max_flux[loc]+background[loc], mags[loc]
				# add the values of the magnitude and its error to a list if the magnitude is detected and the star isn't saturated
				if (max_flux[loc]+background[loc] <= 60000.) and (mags[loc] < 99.) and (mag_errs[loc] < 99.):
					measured_mags.append(mags[loc])
					measured_mag_errs.append(mag_errs[loc])
		# take the average of the measurements
		measured_mags = np.array(measured_mags)
		measured_mag_errs = np.array(measured_mag_errs)
		wts = 1. / measured_mag_errs**2
		print "\n", star, filt
		for m, merr in zip(measured_mags, measured_mag_errs):
			print "{} +/- {}".format(m,merr)
		if len(measured_mags) > 1:
			mag = np.average(measured_mags, weights=wts)
			#mag_err = np.sqrt(1. / np.sum(wts))
			mag_err = np.std(measured_mags) / float(len(measured_mags)) 
		elif len(measured_mags) == 1:
			mag = measured_mags[0][0]
			mag_err = measured_mag_errs[0][0]
		else:
			print "no magnitudes found for {} in {} filter.".format(star, filt)
			mag = None
			mag_err = None
		if mag is not None:
			# save the info
			ref_stars[star]['measured_mag_{}'.format(filt)] = mag
			ref_stars[star]['measured_mag_{}_err'.format(filt)] = mag_err
			ref_stars[star]['mag_{}_offset'.format(filt)] = ref_stars[star]['mag_{}'.format(filt)] - mag
			ref_stars[star]['mag_{}_offset_err'.format(filt)] = np.sqrt(mag_err**2 + ref_stars[star]['mag_{}_err'.format(filt)]**2)
			print "{} {} known mag = {} +/- {}, measured mag = {:1.3f} +/- {:1.3f}, diff = {:1.3f} +/- {:1.4f}, standard dev of measurements = {:1.3f}, npts = {}.".format(star, filt, ref_stars[star]['mag_{}'.format(filt)], ref_stars[star]['mag_{}_err'.format(filt)], mag, mag_err, ref_stars[star]['mag_{}_offset'.format(filt)], ref_stars[star]['mag_{}_offset_err'.format(filt)], np.std(measured_mags), len(measured_mags))

# calculate the zero offset
mag_diffs = {'V': {'diff': [], 'err': []}, 'B': {'diff': [], 'err': []}}
for star in ref_stars.keys():
	for filt in filters:
		if 'mag_{}_offset'.format(filt) in ref_stars[star].keys():
			mag_diffs[filt]['diff'].append(ref_stars[star]['mag_{}_offset'.format(filt)])
			mag_diffs[filt]['err'].append(ref_stars[star]['mag_{}_offset_err'.format(filt)])

# take the average for each filter
print "\n"
mag_offset = {}
for filt in filters:
	wts = 1. / np.array(mag_diffs[filt]['err'])**2
	mag_offset[filt] = {}
	mag_offset[filt]['offset'] = np.average(np.array(mag_diffs[filt]['diff']), weights=wts)
	mag_offset[filt]['err'] = np.sqrt(1. / np.sum(wts))
	print "{}: offset = {:1.3f} +/- {:1.3f}".format(filt, mag_offset[filt]['offset'], mag_offset[filt]['err']) 

# just print out some info
for star in ref_stars.keys():
	print "\n",star
	star_info = ref_stars[star]
	print "known mags: B = {} +/- {}, V = {} +/- {}".format(star_info['mag_B'], star_info['mag_B_err'], star_info['mag_V'], star_info['mag_V_err'])
	B = ref_stars[star]['measured_mag_{}'.format('B')] + mag_offset['B']['offset']
	Berr = np.sqrt(ref_stars[star]['measured_mag_{}_err'.format('B')]**2 + star_info['mag_B_err']**2)
	V = ref_stars[star]['measured_mag_{}'.format('V')] + mag_offset['V']['offset']
	Verr = np.sqrt(ref_stars[star]['measured_mag_{}_err'.format('V')]**2 + star_info['mag_V_err']**2)
	print "corrected mag: B = {:1.3f} +/- {:1.3f}, V = {:1.3f} +/- {:1.3f}".format(B,Berr,V,Verr)
	print "fractional difference: B: {:1.3f}, V: {:1.3f}".format((B - star_info['mag_B'])/star_info['mag_B'], (V - star_info['mag_V'])/star_info['mag_V'])
	print "known B-V = {:1.3f} +/- {:1.3f}, measured = {:1.3f}, corrected = {:1.3f} (corrected for reddening = {:1.3f} +/- {:1.3f}; fdiff = {:1.3f}; diff/err = {:1.3f})".format(star_info['mag_B'] - star_info['mag_V'],  np.sqrt(star_info['mag_B_err']**2 + star_info['mag_V_err']**2), ref_stars[star]['measured_mag_{}'.format('B')] - ref_stars[star]['measured_mag_{}'.format('V')], B-V, (B-0.863)-(V-0.653), np.sqrt(Berr**2 + Verr**2), ((B-V) - (star_info['mag_B'] - star_info['mag_V']))/(star_info['mag_B'] - star_info['mag_V']), ((B-V) - (star_info['mag_B'] - star_info['mag_V'])) / np.sqrt(np.sqrt(Berr**2 + Verr**2)**2 + np.sqrt(star_info['mag_B_err']**2 + star_info['mag_V_err']**2)**2))
