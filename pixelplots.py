#System tools (file in and output ...)
import sys

#Ability to call external programs, shell commands etc
import os
from subprocess import call

#Pythons math package with all the constants and functions
import math as m

#Library includes a zero-finder method (to find the best value for the kernel), also includes methods fit data/curves
from scipy import optimize

#For reading config files, based on key = value structure,
import configparser

#Supressing runtime warning (mainly for maximum iterations of optimize)
import warnings
warnings.filterwarnings("ignore")

#If this is set to false, only minimal console output and no diagnostic images are produced
#only the final plots and datafiles are written
#can be set to True as a second command line argument
PRINTALL = True

#Global variable to set the fitting mechanism, default is lsq, which accepts y-errors and uses Levenberg-Marquart
#alternatively use option odr to also include x-errors
FIT_METHOD = 'lsq'

#################
### FUNCTIONS ###
#################

#For more information on what those functions do, check their respective files for documentation
from plotting import condon, fct_f, fct_result, map_sqrt, plot
from conversion import ceil, floor, convert1200, convert1200_adv, conv_px_per_box, calculate_rms
from convolution import convolve_gauss, fct_gauss, fct_gauss_fit, convert_kpc2px, convert_px2kpc, flatten
from data import print_conv_data, print_data, read_fits
from fitting import fct_lsq, fct_odr, fit_lsq, fit_odr, fit

###########################
### PROGRAM STARTS HERE ###
###########################
###Read in config file from command line argument, exit if it cant be found
if(len(sys.argv) < 2):
	print 'Need filname of config file'
	sys.exit(-1)
tmp_name = str(sys.argv[1])
if(len(sys.argv)==3):
	if (sys.argv[2] == '1'):
		PRINTALL = True
else:
	print 'Too many command line arguments! Quitting...'
	quit()
print 'Full output is set to: ', PRINTALL, ' (change through command line)'
print 'Currently using ', FIT_METHOD, ' as fitting method (change directly in the code).'

#Read config file, object config works almost like a dictoinary, print everything
config = configparser.ConfigParser()
config.read(tmp_name)
if( PRINTALL == True ):
	print config.sections()
	for key in config['names']: print key, ' = ',config['names'][key] 
	for key in config['values']: print key, ' = ',config['values'][key] 

###change working directory to data-path (specified in the ini file)
os.chdir(config['names']['fullpath'])
print  'Current working directory is:', os.getcwd()

###Create png versions of the full fits files, including the (cutting) box using ds9
###also includes the 'rmsbox' to calculate the sigmas for the 3 sigma cutoff
boxsize = 	str(config.getint('values','n_boxes') * conv_px_per_box( config['values'] ) )
cmd = {}
center_x = config.get('values','center_x')
center_y = config.get('values','center_y')
for fname in (['low','high','sfr']):
	oname = config.get('names',fname).rstrip('.fits') + '_ds9box.png'
	cmd[fname] =	('ds9 '+ config.get('names',fname) +
			' -regions system image '+
			' -regions command "box '+
			config.get(fname+'_cutoff_box','center_x') +' '+
			config.get(fname+'_cutoff_box','center_y') +' '+
			config.get(fname+'_cutoff_box','size_x') +' '+
			config.get(fname+'_cutoff_box','size_y') +
			' # color=red" '+
			'-regions command "box '+
			str(center_x) +' '+
			str(center_y) +' '+
			str(boxsize)+' '+
			str(boxsize)+
			' # color=yellow" '+
			' -zoom to fit -scale mode 90  -saveimage png '+
			oname +' -quit')
	if( PRINTALL==True ):
		os.system(cmd[fname])
#Save the DS9 command as string to print to results file later
for key in cmd:
	cmd[key] = cmd[key].rstrip('-quit')

###read in the three maps, 2 radio maps with different frequencys and the SFR map from galex/spitzer
data_l = read_fits(config.get('names','low'), PRINTALL)
data_h = read_fits(config.get('names','high'), PRINTALL)
data_s = read_fits(config.get('names','sfr'), PRINTALL)

###Call the convert function 3 times to convert each map to 1200pc^2 pixels
pixels2d_l = convert1200_adv( data_l, config['values'] )
pixels2d_h = convert1200_adv( data_h, config['values'] )
pixels2d_s = convert1200_adv( data_s, config['values'] )

###Flatten the data to 1D lists
pixels_l = flatten ( pixels2d_l )
pixels_h = flatten ( pixels2d_h )
pixels_s = flatten ( pixels2d_s )

###Apply the cuts to the data based on the signal
###create three empty lists to store the 3 sigma cut data
pix_l_cut = []
pix_h_cut = []
pix_s_cut = []

sigma_high = calculate_rms( data_h, config['high_cutoff_box'])
sigma_low = calculate_rms( data_l, config['low_cutoff_box'])
sigma_sfr = calculate_rms( data_s, config['sfr_cutoff_box'])

sigma = {'low' : sigma_low, 'high' : sigma_high, 'sfr' : sigma_sfr }
#3 sigma cutoff for all datasets
for i in range(len(pixels_l)):
	if(pixels_l[i] > 3. * sigma_low ):
		if(pixels_h[i] > 3. * sigma_high): 
			if(pixels_s[i] > 5. * sigma_sfr ): 
				pix_l_cut.append(pixels_l[i])
				pix_h_cut.append(pixels_h[i])
				pix_s_cut.append(pixels_s[i])

print 'RMS for the 3 different maps, used as sigma for 3 sigma cutoff:\n', sigma
###Define alpha array to store the spectral index (assuming a simple power law)
alpha = []
for i in range(len(pix_l_cut)):
	alpha.append( m.log10( pix_l_cut[i] /pix_h_cut[i] ) / m.log10( config.getfloat('values','freq_low') / config.getfloat('values','freq_high') )) 

#Some diagnostic output of how many points were cut
print 'Number of points in sample:', len( pixels_l ) 
print 'Number of points after 3 sigma cutoff:', len( pix_l_cut )

###Print all pixel data to file (after the 3 sigma cut and before conversion to SFR)
mean = print_data( config, pix_l_cut, pix_h_cut, pix_s_cut, alpha )
mean_old = 0
tmp = flatten( data_l )
for i in range(len(tmp)):
	mean_old += tmp[i]
mean = 0
for i in range(len(pixels_l)):
	mean += pixels_l[i]


mean *=  conv_px_per_box( config['values'] )**2
print 'Total total sum of LOFAR image:\t', '%0.3f' % mean_old
print 'Total sum in the cut LOFAR map is:', '%0.3f' % mean

###Now convert lofar and swrt to SFR sufrace density using condon-relation
pix_l_cut = condon( pix_l_cut, config.getfloat('values','FWHM'), config.getfloat('values','freq_low') )
pix_h_cut = condon( pix_h_cut, config.getfloat('values','FWHM'), config.getfloat('values','freq_high') )

### Diganosis only: Create square rooted map of the 1200pc data (to check if the cutting worked as intended)
if( PRINTALL == True ):
	map_sqrt(pixels2d_l, config, 'low')
	map_sqrt(pixels2d_h, config, 'high')
	map_sqrt(pixels2d_s, config, 'sfr')


###Fitting for both datasets
#Fitting low /LOFAR
a_l , a_l_err , b_l, b_l_err, chi_l = fit(	pix_s_cut, 
											pix_l_cut,
											output=True,
											case=FIT_METHOD)
#Fitting high / WSRT
a_h , a_h_err , b_h, b_h_err, chi_h = fit(	pix_s_cut,
											pix_h_cut,
											output=True,
											case=FIT_METHOD)

############
### This concludes the creation of the 'simple' pixel plots (without the actual plotting)
### Below the data for the convolved pixel plots is calculated
############

###Finding optimal gaussian kernel for both radio maps
### Calculating pixel values and fits based on the optimal kernel

###low frequency radio map
print 'Finding optimal gaussian kernel for LOFAR (low) data. This may take a moment...'
optimal_sigma_l = optimize.fsolve(fct_gauss_fit, config.getfloat('values','sigma_conv'), args=(data_s, pixels_l, pixels_h, sigma , config, 'low', PRINTALL), maxfev = 10 )

conv_pix_cut_low, conv_pix_l_cut, conv_alpha_l, a_smooth_l, b_smooth_l = fct_gauss(optimal_sigma_l[0], data_s, pixels_l, pixels_h, sigma , config, 'low', PRINTALL )

fit(conv_pix_cut_low, conv_pix_l_cut, output=True)

#Diffusion length is the FWHM/2 of the final image, FWHM are square added first
optimal_sigma_l[0] = m.sqrt( m.pow(2.3548 *optimal_sigma_l[0],2) + m.pow(1.2,2) )/2.
print 'Final value for Diffusion length:\t','%0.3f' % optimal_sigma_l[0], 'kpc'

###high frequency radio map
print 'Finding optimal gaussian kernel for WSRT (high) data. This may take a moment...'
optimal_sigma_h = optimize.fsolve(fct_gauss_fit, config.getfloat('values','sigma_conv'), args=(data_s, pixels_l, pixels_h, sigma , config, 'high', PRINTALL ), maxfev = 10 )

conv_pix_cut_high, conv_pix_h_cut, conv_alpha_h, a_smooth_h, b_smooth_h = fct_gauss(optimal_sigma_h[0], data_s, pixels_l, pixels_h, sigma, config, 'high', PRINTALL )

fit(conv_pix_cut_high, conv_pix_h_cut, output=True)

#Diffusion length is the FWHM/2 of the final image, FWHM are square added first
optimal_sigma_h[0] = m.sqrt( m.pow(2.3548 *optimal_sigma_h[0],2) + m.pow(1.2,2) )/2.
print 'Final value for Diffusion length:\t', '%0.3f' % optimal_sigma_h[0], 'kpc'

############
###Finally creating the pixel plots for all the data 
############
print 'Creating final images and writing results to file ...'
'''if(PRINTALL==True):
	try:
		os.system('gnuplot fit_test.gp')
	except:
		pass
'''
#Plotting low / LOFAR
plot(pix_s_cut, pix_l_cut, alpha, a_l, b_l, config['names'],config['boundaries'] , 'low' )
#Plotting high / WSRT
plot(pix_s_cut, pix_h_cut, alpha, a_h, b_h, config['names'],config['boundaries'], 'high' )
#Plotting convolved data (low freq.)
plot(conv_pix_cut_low, conv_pix_l_cut, conv_alpha_l, a_smooth_l, b_smooth_l, config['names'],config['boundaries'], 'conv_low' , optimal_sigma_l[0])
#Plotting convolved data (high freq.)
plot(conv_pix_cut_high, conv_pix_h_cut, conv_alpha_h, a_smooth_h, b_smooth_h, config['names'],config['boundaries'], 'conv_high' , optimal_sigma_h[0])

###Print all convolved pixel data to file (after the 3 sigma cut and conversion to SFR)
mean = print_conv_data( config,conv_pix_cut_low, conv_pix_l_cut , conv_alpha_l ,'_conv_low' )
mean = print_conv_data( config, conv_pix_cut_high, conv_pix_h_cut, conv_alpha_h, '_conv_high' )

### Write the final results (fits and diffusion lengths) to file
# this makes use of the configparser library again 
tmp = config.get('names','galaxy').split(' ')
results_ini = tmp[0]+'_'+tmp[1]+'_results.ini'
res_out = configparser.ConfigParser()

res_out['name'] = {'fname': config.get('names','galaxy'),
					'ds9cmd_low' : cmd['low'],
					'ds9cmd_high' : cmd['high'],
					'ds9cmd_sfr' : cmd['sfr']}

res_out['Low_freq_fit'] =	{'#Fit results for pixel plots, with std errors\n'
							'a': str(a_l),
							'a_err': str(a_l_err),
							'b': str(b_l_err),
							'b_err': str(b_l_err),
							'chi_sqr': str(chi_l)}

res_out['High_freq_fit'] =	{'a': str(a_h),
							'a_err': str(a_h_err),
							'b': str(b_h_err),
							'b_err': str(b_h_err),
							'chi_sqr': str(chi_h)}

res_out['Conv_results'] =	{'#Diffusion length from gaussian kernel in kpc \n'
							'sigma_l': str(optimal_sigma_l[0]),
							'sigma_h': str(optimal_sigma_h[0])}

with open(results_ini, 'w') as configfile:
	res_out.write(configfile)

print "Finished!"



