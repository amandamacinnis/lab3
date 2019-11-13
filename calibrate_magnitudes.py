import numpy as np

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


# TODO: re-run source extractor and check the max flux above background level to make sure star is not saturated

# note the magnitudes are found from the https://www.aavso.org/ finding charts
# create a dict of reference stars to use for calibration
ref_stars = {}
ref_stars['merope'] = {'ra': get_ra(3., 46., 19.5738), 
	'dec': get_dec(23., 56., 54.0812), 
	'mag_V': 4.180, 'mag_B': 4.120,       # website lists source as 'BSC' 
	'cat_fname': 'source_extractor_output/merope_{}_1s_{:02d}.cat', 
	'cat_aper': 11.} # radius of aperture in pixels
ref_stars['atlas'] = {'ra': get_ra(3., 49., 9.57),
	'dec': get_dec(24., 3., 12.3),
	'mag_V': 3.630, 'mag_B': 3.540,
	'cat_fname': 'source_extractor_output/atlas_{}_1s_{:02d}.cat',
	'cat_aper': 12.}



# look at the measured magnitude of each
for star in ref_stars.keys():
	ra = ref_stars[star]['ra']
	dec = ref_stars[star]['dec']
	# loop through each catalog for the given exposure time and filter
	for filt in filters:
		measured_mags = []
		measured_mag_errs = []
		for i in range(5):
			# load the catalog
			ras, decs, fluxs, flux_errs, mags, mag_errs = np.loadtxt(ref_stars[star]['cat_fname'].format(filt, i), unpack=True)
			# get a list of angular distances (in arcmin) between where we think the star is, and the objects detected by source extractor
 			_, _, distances = get_angular_dist(ra, dec, ras, decs)
			# make sure we have some reasonable star matching the coordinates we gave
			if np.min(distances) > 1.:
				print "WARNING: Could not find a star within 1 arcmin of the given coordinates."

			# the star we want will (hopefully) have the minimum angular distance 
			loc = np.where(distances == np.min(distances))

			# make sure the star is not oversaturated
			avg_count_per_px = fluxs[loc] / (np.pi * ref_stars[star]['cat_aper']**2)
			if avg_count_per_px > 60000:
				print "WARNING: average counts per pixel is {:d}".format(avg_count_per_px)

			# add the values of the magnitude and its error to a list
			measured_mags.append(mags[loc])
			measured_mag_errs.append(mag_errs[loc])
			#print i, ra, dec
			#print loc, ras[loc], decs[loc], np.min(distances)
		# take the average of the measurements
		measured_mags = np.array(measured_mags)
		measured_mag_errs = np.array(measured_mag_errs)
		wts = 1. / measured_mag_errs**2
		mag = np.average(measured_mags, weights=wts)
		mag_err = np.sqrt(1. / np.sum(wts))
		ref_stars[star]['measured_mag_{}'.format(filt)] = mag
		ref_stars[star]['measured_mag_{}_err'.format(filt)] = mag_err
		ref_stars[star]['mag_{}_offset'] = ref_stars[star]['mag_{}'.format(filt)] - mag
		print star, filt, ref_stars[star]['mag_{}'.format(filt)], mag, mag_err, ref_stars[star]['mag_{}'.format(filt)] - mag
