#Pythons math package with all the constants and functions
import math as m
#Library for python plotting
import matplotlib.pyplot as plt
#Library to manipulate images, used to smooth with gaussian kernel
from scipy import ndimage

from conversion import convert_resolution_adv
from plotting import condon
from fitting import fct_lsq, fct_odr, fit_lsq, fit_odr, fit

###Flatten 2D arrays in to a long list; this is easier to plot
def flatten( data ):
	data.flatten()
	data = [item for sublist in data for item in sublist]
	return data

###Function that calculates the fit parameter a (slope) as a function of gaussian kernel width sigma
def fct_gauss(sigma, data_s, pixels_l, pixels_h, cutoff, config, opt, PRINTALL, CALIB_ERR, case):
	#Compute new spectral indices (may have changed due to different cutting)
	conv_alpha = []
	for i in range(len(pixels_l)):
		conv_alpha.append( m.log10( m.fabs(pixels_l[i] /pixels_h[i]) ) / m.log10( config.getfloat('values','freq_low') / config.getfloat('values','freq_high') ))
	#Convolve with gaussian 
	conv_s2d = convolve_gauss( data_s, config['values'], sigma, opt, PRINTALL )
	conv_pix2d = convert_resolution_adv( conv_s2d, config['values'] )
	conv_pix = flatten( conv_pix2d )
	conv_pix_cut = []
	conv_pix_cut_err = []
	conv_pix_l_cut = []
	conv_pix_l_cut_err = []
	conv_pix_h_cut = []
	conv_pix_h_cut_err = []
	conv_alpha_cut = []
	#applying the same 3 sigma cut as for the non-covolved data
	for i in range(len(pixels_l)):
		if(pixels_l[i] > 3. * cutoff['low'] ):
			if(pixels_h[i] > 3. * cutoff['high'] ): 
				if(conv_pix[i] > 3. * cutoff['sfr'] ):
					if( conv_alpha[i] <= config.getfloat( 'boundaries','low' ) ):
						conv_alpha_cut.append( conv_alpha[i] )
						conv_pix_cut.append( conv_pix[i] )
						conv_pix_cut_err.append( calc_error( conv_pix[i], cutoff['sfr'], CALIB_ERR ) )
						conv_pix_l_cut.append( pixels_l[i] )
						conv_pix_l_cut_err.append( calc_error( pixels_l[i], cutoff['low'], CALIB_ERR ) )
						conv_pix_h_cut.append( pixels_h[i] )
						conv_pix_h_cut_err.append( calc_error( pixels_h[i], cutoff['high'], CALIB_ERR ) )
		#Now apply the condon relation to the radio map set via parameter 'opt' and return values
		###CHANGE TO INCLUDE ERRORS!
	if(opt == 'high'):
		conv_pix_h_cut = condon( conv_pix_h_cut, config.getfloat('values','FWHM'), config.getfloat('values','freq_high') )
		conv_pix_h_cut_err = condon( conv_pix_h_cut_err, config.getfloat('values','FWHM'), config.getfloat('values','freq_high') )
		a_h , _ , b_h, _, _ = fit(conv_pix_cut, conv_pix_h_cut, val_x_err=conv_pix_cut_err, val_y_err=conv_pix_h_cut_err, case=case)
		return conv_pix_cut, conv_pix_cut_err, conv_pix_h_cut, conv_pix_h_cut_err, conv_alpha_cut, a_h, b_h
	elif(opt == 'low'):
		conv_pix_l_cut = condon( conv_pix_l_cut, config.getfloat('values','FWHM'), config.getfloat('values','freq_low') )
		conv_pix_l_cut_err = condon( conv_pix_l_cut_err, config.getfloat('values','FWHM'), config.getfloat('values','freq_low') )
		a_l , _ , b_l, _, _ = fit(conv_pix_cut, conv_pix_l_cut, val_x_err=conv_pix_cut_err, val_y_err=conv_pix_l_cut_err, case=case)
		return conv_pix_cut, conv_pix_cut_err, conv_pix_l_cut, conv_pix_l_cut_err, conv_alpha_cut, a_l, b_l

###Gaussian Kernel function for optimization purposes, because optimize searches for zeros by default
def fct_gauss_fit(sigma, data_s, pixels_l, pixels_h, cutoff, config, opt, PRINTALL, CALIB_ERR, case ):
	_,_,_,_,_,x,_ = fct_gauss(sigma, data_s, pixels_l, pixels_h, cutoff, config, opt, PRINTALL, CALIB_ERR, case)
	return x-1.

###Convolution of the 2d image with gaussian kernel
###CHECK THE UNITS OF SIGMA HERE!
def convolve_gauss( data , cfg, sigma_in, opt , PRINTALL):
	sigma = convert_kpc2px(sigma_in, cfg)
	if(PRINTALL == True):
		print 'current sigma:','%0.3f' %  convert_px2kpc(sigma,cfg) , 'in kpc',  '%0.3f' % sigma, 'in px'
	#Convolution is identical with gaussian_filter, needs standard-deviation/pixel!
	res = ndimage.filters.gaussian_filter(data, sigma)
	if( PRINTALL == True ):
		#Now plot the map (compare to fits file if something doesnt work)
		plt.imshow(res, cmap='gray')
		plt.colorbar()
		plt.ylim(0,1024)
		fname_image = 'convolution_'+opt+'.png'
		plt.savefig( fname_image )
		plt.clf()
	return res

### Unit conversion from pixel distance to distance in kiloparsec
def convert_px2kpc(px, cfg):
	kpc_per_arcsec = cfg.getfloat('distance') * m.tan( 2 * m.pi / (360. * 3600.) ) / 1000.
	return float( px / cfg.getfloat('pixel_per_arcsec') * kpc_per_arcsec )
### Unit conversion from distance in kiloparsec to pixel distance
def convert_kpc2px(kpc, cfg):
	kpc_per_arcsec = cfg.getfloat('distance') * m.tan( 2 * m.pi / (360.* 3600.) ) / 1000.
	return float(kpc * cfg.getfloat('pixel_per_arcsec') / kpc_per_arcsec )


def calc_error(val, noise, calibration):
	error = m.sqrt( m.pow(noise,2.) + m.pow(val*calibration,2.) )
	return error

