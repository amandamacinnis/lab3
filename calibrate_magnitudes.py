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
	if 'B' in fname:
		cat_files['B'].append(fname)
	elif 'V' in fname:
		cat_files['V'].append(fname)

# note the magnitudes are found from the https://www.aavso.org/ finding charts
# create a dict of reference stars to use for calibration
ref_stars = {}
ref_stars['merope'] = {'ra': get_ra(3., 46., 19.5738), 
	'dec': get_dec(23., 56., 54.0812), 
	'mag_V': 4.180, 'mag_B': 4.120}       # website lists source as 'BSC' 
ref_stars['atlas'] = {'ra': get_ra(3., 49., 9.57),
	'dec': get_dec(24., 3., 12.3),
	'mag_V': 3.630, 'mag_B': 3.540} # BSC
ref_stars['hd23778'] = {'ra': get_ra(3., 48., 34.78),
	'dec': get_dec(24., 10., 53.1),
	'mag_V': 9.063, 'mag_B': 9.519} # B from Tycho-2, V from ASAS
ref_stars['hd23778_ref1'] = {'ra': get_ra(3., 48., 30.07),
	'dec': get_dec(24, 20., 44.5),
	'mag_V': 6.964, 'mag_B': 7.023} # WBVR 
ref_stars['hd23778_ref2'] = {'ra': get_ra(3., 49., 25.98),
	'dec': get_dec(24., 14., 52.5),
	'mag_V': 7.993, 'mag_B': 8.095} # WBVR
ref_stars['hd23511_ref1'] = {'ra': get_ra(3., 46., 19.55),
	'dec': get_dec(23., 56., 53.3),
	'mag_V': 4.180, 'mag_B': 4.120} # BSC
ref_stars['hd23511_ref2'] = {'ra': get_ra(3., 46., 34.26),
	'dec': get_dec(24., 8., 17.3),
	'mag_V': 11.745, 'mag_B': 12.438} # APASS
ref_stars['celaeno_ref1'] = {'ra': get_ra(3., 45., 12.51),
	'dec': get_dec(24., 28., 1.9),
	'mag_V': 4.3, 'mag_B': 4.19} # BSC
#ref_stars['hd23479_ref1'] = {'ra': get_ra(3., 45., 49.57),
#	'dec': get_dec(24., 22., 3.4),
#	'mag_V': 3.870, 'mag_B': 3.800} # BSC




# look at the measured magnitude of each
for star in ref_stars.keys():
	ra = ref_stars[star]['ra']
	dec = ref_stars[star]['dec']
	# loop through each catalog for the given  filter
	for filt in filters:
		measured_mags = []
		measured_mag_errs = []
		for fname in cat_files[filt]:
			# load the catalog
			ras, decs, fluxs, flux_errs, mags, mag_errs, background, max_flux = np.loadtxt('source_extractor_output/'+fname, unpack=True)
			# get a list of angular distances (in arcmin) between where we think the star is, and the objects detected by source extractor
 			_, _, distances = get_angular_dist(ra, dec, ras, decs)
			# make sure we have some reasonable star matching the coordinates we gave
			if np.min(distances) <= 1.:
				# the star we want will (hopefully) have the minimum angular distance 
				loc = np.where(distances == np.min(distances))

				# add the values of the magnitude and its error to a list if the magnitude is detected and the star isn't saturated
				if (max_flux[loc]+background[loc] <= 50000.) and (mags[loc] < 99.) and (mag_errs[loc] < 99.):
					measured_mags.append(mags[loc])
					measured_mag_errs.append(mag_errs[loc])
		# take the average of the measurements
		measured_mags = np.array(measured_mags)
		measured_mag_errs = np.array(measured_mag_errs)
		wts = 1. / measured_mag_errs**2
		if len(measured_mags) > 1.:
			mag = np.average(measured_mags, weights=wts)
			mag_err = np.sqrt(1. / np.sum(wts))
			ref_stars[star]['measured_mag_{}'.format(filt)] = mag
			ref_stars[star]['measured_mag_{}_err'.format(filt)] = mag_err
			ref_stars[star]['mag_{}_offset'.format(filt)] = ref_stars[star]['mag_{}'.format(filt)] - mag
			print "{} {} known mag = {}, measured mag = {:1.3f} +/- {:1.3f}, diff = {:1.3f}, standard dev = {:1.3f}.".format(star, filt, ref_stars[star]['mag_{}'.format(filt)], mag, mag_err, ref_stars[star]['mag_{}'.format(filt)] - mag, np.std(measured_mags))

# calculate the zero offset
mag_diffs = {'V': {'diff': [], 'err': []}, 'B': {'diff': [], 'err': []}}
for star in ref_stars.keys():
	for filt in filters:
		if 'mag_{}_offset'.format(filt) in ref_stars[star].keys():
			mag_diffs[filt]['diff'].append(ref_stars[star]['mag_{}_offset'.format(filt)])
			mag_diffs[filt]['err'].append(ref_stars[star]['measured_mag_{}_err'.format(filt)])

# take the average for each filter
mag_offset = {}
for filt in filters:
	wts = 1. / np.array(mag_diffs[filt]['err'])**2
	mag_offset[filt] = {}
	mag_offset[filt]['offset'] = np.average(np.array(mag_diffs[filt]['diff']), weights=wts)
	mag_offset[filt]['err'] = np.sqrt(1. / np.sum(wts))
	print "{}: offset = {} +/- {}".format(filt, mag_offset[filt]['offset'], mag_offset[filt]['err']) 
