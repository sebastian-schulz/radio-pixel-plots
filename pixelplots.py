#System tools (file in and output ...)
import sys

#Ability to call external programs, shell commands etc
import os
from subprocess import call
#Mighty numerical library of ptyhon
import numpy as np

#Pythons math package with all the constants and functions
import math as m

#Library for plotting and images
import matplotlib
import matplotlib.pyplot as plt

#Library to understand the fits format
from astropy.io import fits

#Library to perform linear regression (with x and y errors!)
from scipy import odr
#Library to manipulate images, used to smooth with gaussian kernel
from scipy import ndimage
#Library includes a zero-finder method (to find the best value for the kernel)
from scipy import optimize

#For LateX fonts and symbols in plots
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['T1']})# T1 is the LaTex standard font, also used by scr familiy documents as default
rc('text', usetex=True)

#For reading config files, based on key = value structure,
import configparser

#Supressing runtime warning (mainly for maximum iterations of optimize)
import warnings
warnings.filterwarnings("ignore")

#If this is set to false, only minimal console output and no diagnostic images are produced
#only the final plots and datafiles are written
#can be set to True as a second command line argument
PRINTALL = False

#################
### FUNCTIONS ###
#################

###Read image from fits file and store it in 2d numpy array
def read_fits( fname ):
	#Reading in the image file and checking size as well as shape
	data = fits.getdata( fname )
	#if it has the shape (1,1, X, Y), SLICE it to a 2D array (x,y)
	try:
		data = data[0][0][:][:]
	except:
		pass
	if (PRINTALL == True):
		#Check for size and shape of the array
		#print 'Data type and shape of', fname ,':',type(data), data.shape 
		#Now plot all the map (compare to fits file if something doesnt work)
		plt.imshow(data, cmap='gray')
		plt.colorbar()
		plt.ylim(0,1024)
		fname_image = fname.rstrip('.fits') + '.png'
		plt.savefig( fname_image )
		plt.clf()
	return data

###Flatten 2D arrays in to a long list; this is easier to plot
def flatten( data ):
	data.flatten()
	data = [item for sublist in data for item in sublist]
	return data

###Converts 2d array (image) to pixels of 1.2kpc length and cuts it to the given box size; basically produces a smaller and lower-resolution version
def convert1200( data, cfg ):
	#calculate starting positions (bottom left corner) in pixels
	s_x = cfg.getint('center_x') - cfg.getint('n_boxes')/2 * cfg.getint('pixels_per_box')
	s_y = cfg.getint('center_y') - cfg.getint('n_boxes')/2 * cfg.getint('pixels_per_box')
	#Creating empty array to store newly created image
	pixels = np.zeros(shape = ( cfg.getint('n_boxes'), cfg.getint('n_boxes') ) )

	#Loop over the 2d array. This makes use of cutting decimals by int-cutting
	#saves some computation time
	for x in range( cfg.getint('n_boxes') * cfg.getint('pixels_per_box') ):
		for y in range( cfg.getint('n_boxes') * cfg.getint('pixels_per_box') ):
			pixels[ x/cfg.getint('pixels_per_box'), y/cfg.getint('pixels_per_box') ] += data[ x+s_x][ y+s_y ]
	#Return average values (total signal divided by number of old pixels per new pixel)
	return pixels / (cfg.getfloat('pixels_per_box'))**2
	

#PASSED BUGTET (July 2018), gives the same results as DS9 with the same boxes
def convert1200_adv( data, cfg ):
	''' Idea: allow for pixels_per_box to be float!
assumption 1 : middle pixel for p and p' have the same position (center)
calculate everything in p distances: (n_boxes should be odd!)
bottom left x = center_x - px_per_box * n_boxes/2
same for y
loop over all p'
	loop over p inside p'
		go to bottom left corner
		find pixel p fully inside
		same for top right corner
		-> calculate sum of all pixels fully inside
		now borders:
		4 corners
		2 edges horizontal
		2 edges vertical
'''
	#position of center of bottom left new pixel in old pixels
	#for some reason x and y axis are in the wrong order in the array
	img_bot_left_y = cfg.getint('center_x')-1 - (cfg.getint('n_boxes')/2. -0.5) * cfg.getfloat('pixels_per_box')
	img_bot_left_x = cfg.getint('center_y')-1 - (cfg.getint('n_boxes')/2.-0.5) * cfg.getfloat('pixels_per_box')
	#new pixel array
	N=0
	pixels = np.zeros(shape = ( cfg.getint('n_boxes'), cfg.getint('n_boxes') ) )
	for x in range( cfg.getint('n_boxes') ):
		for y in range( cfg.getint('n_boxes') ):
			#position of current pixel in old pixels
			img_current_x = img_bot_left_x + x * cfg.getfloat('pixels_per_box')
			img_current_y = img_bot_left_y + y * cfg.getfloat('pixels_per_box')
			#positions of the corners of the new pixels in old pixels
			pix_x_min = img_current_x - cfg.getfloat('pixels_per_box')/2.
			pix_x_max = img_current_x + cfg.getfloat('pixels_per_box')/2.
			pix_y_min = img_current_y - cfg.getfloat('pixels_per_box')/2.
			pix_y_max = img_current_y + cfg.getfloat('pixels_per_box')/2.
			#loop over the old pixels fully inside the current new pixel
			for p_x in range( ceil( pix_x_min ), floor( pix_x_max ), 1):
				for p_y in range( ceil( pix_y_min ), floor( pix_y_max ), 1):
					pixels[x,y] += data[ p_x, p_y ]
					N+=1
			#now loop over the edges, without the corners
			factor_max_x = m.fabs(floor(pix_x_max)-pix_x_max)
			factor_min_x = m.fabs(ceil(pix_x_min)-pix_x_min)
			factor_max_y = m.fabs(floor(pix_y_max)-pix_y_max)
			factor_min_y = m.fabs(ceil(pix_y_min)-pix_y_min)
			#if (factor_max_x > 1 or factor_max_y > 1 or factor_min_x > 1 or factor_min_y > 1):
			#	print ' ERROR!'
			#print 'x and y factors \t', factor_max_x, '\t', factor_max_y
			for p_y in range( ceil( pix_y_min ), floor( pix_y_max ), 1):
				N+= factor_min_x + factor_max_x
				pixels[x,y] += data[floor( pix_x_min ) ,p_y ] * factor_min_x
				pixels[x,y] += data[ceil( pix_x_max ) ,p_y ] * factor_max_x
			for p_x in range( ceil( pix_x_min ), floor( pix_x_max ), 1):
				N+= factor_min_y + factor_max_y
				pixels[x,y] += data[ p_x , floor( pix_y_min ) ] * factor_min_y
				pixels[x,y] += data[ p_x, ceil( pix_y_max ) ] * factor_max_y

			#finally the four corner pixels:
			pixels[x,y] += data[floor(pix_x_min), floor(pix_y_min)] *factor_min_x * factor_min_y
			pixels[x,y] += data[floor(pix_x_min), ceil(pix_y_max)] * factor_min_x *factor_max_y
			pixels[x,y] += data[ceil(pix_x_max), ceil(pix_y_max)] * factor_max_x * factor_max_y
			pixels[x,y] += data[ceil(pix_x_max), floor(pix_y_min)] * factor_max_x * factor_min_y
			N+= factor_min_x * factor_min_y + factor_min_x *factor_max_y +factor_max_x * factor_max_y + factor_max_x * factor_min_y
	#if(PRINTALL==True ):
	#print 'Total number of original pixels used in conv:\t', N
	return pixels / (cfg.getfloat('pixels_per_box'))**2


###Prints a map of the 1200pc pixel data (with sqrt(values) for better visibility)
### only used to debug/compare with volker
def map_sqrt(pixels , cfg, galaxyname ):
	tmp = np.zeros(shape = ( len(pixels), len(pixels)) )
	for i in range(len(pixels)):
		for j in range(len(pixels)):
			tmp[i,j] = m.sqrt(m.fabs(pixels[i,j]))
	#Plot the modified data as a greyscale map (compare with original)
	plt.imshow( tmp, cmap='gray')
	plt.colorbar()
	plt.ylim(0, cfg.getint('values','n_boxes') - 1 )
	image_name = config.get('names', galaxyname )
	image_name = image_name.rstrip('.fits')
	image_name += '_sqrtmap.png'
	#And save the plot to file:
	plt.savefig( image_name )
	plt.clf() #clears the plot for further plotting

### Calculates SFR from radio data using condons relation
### x in Jy/beam, fwhm in arcsec, freq in MHz
def condon( x, fwhm, freq ):
	for i in range(len(x)):
		x[i] = x[i]*3.31e3*(freq/1400.)**(.8) *(fwhm)**(-2.)
	return x

###Fitting function with x and y errors see Python Scipy ODR
### fits a *x +b =y
### needs x and y values as list or array
def fit_odr (val_x, val_y, val_x2d_err=None, val_y2d_err=None, output=False):
	val_log_x = []
	val_log_y = []
	for i in range(len(val_x)):
		if (val_x[i] < 0 ):
			val_log_x.append(0.)
		else:
			val_log_x.append( m.log10( val_x[i] ) )
		if( val_y[i] < 0 ):
			val_log_y.append(0.)
		else:
			val_log_y.append( m.log10( val_y[i] ) )
	
	myodr = odr.odr(fct_odr, [1.,0.], val_log_y, val_log_x, full_output=1)
	outodr = odr.Output(myodr)
	if(output == True):
		print '########## FIT RESULTS ###########'
		print 'a =\t', '%0.3f' % outodr.beta[0], '+/-\t', '%0.3f' % outodr.sd_beta[0]
		print 'b =\t','%0.3f' %  outodr.beta[1], '+/-\t', '%0.3f' % outodr.sd_beta[1]
		print 'Chi squared = \t','%0.3f' %  outodr.sum_square
		#print outodr.sum_square_delta, '\t', outodr.sum_square_eps
	return outodr.beta[0], outodr.sd_beta[0], outodr.beta[1] , outodr.sd_beta[1], outodr.sum_square

###Linear fitting function for ODR
def fct_odr(B, val_x):
	return B[0]*val_x+B[1]

### y=x function
def fct_f (x):
	return x

###Function with Volkers (gnuplot) best fit values JUST FOR TESTING NGC5194
### needs 'case' parameter = low LOFAR or =high WSRT
def fct_volker(x, case):
	if(case=='low'):
		return x**0.51 *10**-0.69
	elif(case=='high'):
		return x**0.76 *10**-0.08
	else:
		print 'Value error for gnuplot fit!'
		return 0.

###Exponential function to correctly plot the linear (log-log) fit
def fct_result(x, a, b):
	return x**a * 10**b

###Plotting function, also sorts data according to spectral index 
def plot(val_x, val_y, alpha, a, b, cfg, boundaries, case, sigma=None ):
	#Define empty arrys/lists to store portions of data based on different spectral indices
	val_flat_x = []
	val_flat_y = []
	val_med_x = []
	val_med_y = []
	val_steep_x = []
	val_steep_y = []
	val_outlier_x = []
	val_outlier_y = []
	#Sort the data into different ranges of alpha
	for i in range(len(val_x)):
		if( float(boundaries['low']) > alpha[i] > float(boundaries['med']) ):
			val_flat_x.append( val_x[i] )
			val_flat_y.append( val_y[i] )
		elif( float(boundaries['med']) > alpha[i] > float(boundaries['high']) ):
			val_med_x.append( val_x[i] )
			val_med_y.append( val_y[i] )
		elif( float(boundaries['low']) > alpha[i] ):
			val_steep_x.append( val_x[i] )
			val_steep_y.append( val_y[i] )
		else:
			val_outlier_x.append[i]
			val_outlier_y.append[i]
	#Create the plot: 
	t = np.linspace(1e-4,1e-1) #datapoints in the plotting range
	plt.clf() #clear pervious plots
	#double-logarithmic plots with advanced options (color, marker, ...)
	plt.loglog( val_steep_x , val_steep_y, marker='^', linestyle='None', color='b')
	plt.loglog( val_med_x , val_med_y , marker='o', linestyle='None', color='g')
	plt.loglog( val_flat_x, val_flat_y , marker='v', linestyle='None', color='r')
	plt.loglog( t, fct_f(t) ,linestyle='--')
	plt.loglog(t, fct_result(t, a, b) , linestyle='-')
	plt.loglog( val_outlier_x, val_outlier_y , marker='s', linestyle='None', color='tab:gray')
	#if(case == 'high' or case == 'low'):
	#	plt.loglog(t, fct_volker(t, case) , linestyle='-.')
	#Create a legend, labels ... be careful with latex symbols ...
	plt.legend([boundaries['low']+r'$\,> {\alpha}\,$', boundaries['low']+r'$\,> {\alpha}>\,$' +boundaries['med'] ,boundaries['med'] + r'$\,> {\alpha} >\,$' +boundaries['med'] , r'{Condon}', r'{Python fit}', r'{Outlier}']) #r'{Gnuplot fit}'
	plt.grid(True)
	plt.ylabel(r'{SFR surface density based on Radio data}')
	plt.xlabel(r'{SFR surface density based on GALEX/Spitzer}')
	plt.xlim(1e-4,1e-1)
	plt.ylim(1e-4,1e-1)
	#Save plots as pdf
	if(case == 'radio_high' or case == 'radio_low'):
		outfile = cfg.get(case).rstrip('.fits') + '_pixel.pdf'
	else:
		outfile = cfg.get(case) + '_pixel.pdf'
	
	if(case == 'conv_low' or case == 'conv_high'):
		title = cfg.get('galaxy') + ' with Gaussian kernel, $l = $'+'%0.2f' % sigma
	else:
		title = cfg.get('galaxy')
	plt.title( title )
	plt.savefig( outfile )
	plt.clf() #clean for further plotting
	return None

###Convolution of the 2d image with gaussian kernel
###CHECK THE UNITS OF SIGMA HERE!
def convolve_gauss( data , cfg, sigma_in, opt ):
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

###Function that calculates the fit parameter a (slope) as a function of gaussian kernel width sigma
def fct_gauss(sigma, data_s, pixels_l, pixels_h, config, opt):
	#Convolve with gaussian 
	conv_s2d = convolve_gauss( data_s, config['values'], sigma, opt )
	conv_pix2d = convert1200_adv( conv_s2d, config['values'] )
	conv_pix = flatten( conv_pix2d )
	conv_pix_cut = []
	conv_pix_l_cut = []
	conv_pix_h_cut = []
	#applying the same 3 sigma cut as for the non-covolved data
	for i in range(len(pixels_l)):
		if(pixels_l[i] > 3* config.getfloat('values','sigma_low')):
			if(pixels_h[i] > 3.* config.getfloat('values','sigma_high')): 
				if(conv_pix[i] > 1.*config.getfloat('values','sigma_diff')): 
					conv_pix_cut.append( conv_pix[i] )
					conv_pix_l_cut.append( pixels_l[i] )
					conv_pix_h_cut.append( pixels_h[i] )
	#Compute new spectral indices (may have changed due to different cutting)
	conv_alpha = []
	for i in range(len(conv_pix_cut)):
		conv_alpha.append( m.log10( conv_pix_l_cut[i] /conv_pix_h_cut[i] ) / m.log10( config.getfloat('values','freq_low') / config.getfloat('values','freq_high') ))
		#Now apply the condon relation to the radio map set via parameter 'opt' and return values
	if(opt == 'high'):
		conv_pix_h_cut = condon( conv_pix_h_cut, config.getfloat('values','FWHM'), config.getfloat('values','freq_high') )
		a_h , _ , b_h, _, _ = fit_odr(conv_pix_cut, conv_pix_h_cut)
		return conv_pix_cut, conv_pix_h_cut, conv_alpha, a_h, b_h
	elif(opt == 'low'):
		conv_pix_l_cut = condon( conv_pix_l_cut, config.getfloat('values','FWHM'), config.getfloat('values','freq_low') )
		a_l , _ , b_l, _, _ = fit_odr(conv_pix_cut, conv_pix_l_cut)
		return conv_pix_cut, conv_pix_l_cut, conv_alpha, a_l, b_l

###Gaussian Kernel function for optimization purposes, because optimize searches for zeros by default
def fct_gauss_fit(sigma, data_s, pixels_l, pixels_h, config, opt ):
	_,_,_,x,_ = fct_gauss(sigma, data_s, pixels_l, pixels_h, config, opt)
	return x-1.

###Prints data arrays to file as tab separated values
def print_data( cfg, data_low, data_high, data_hyb, alpha):
	mean = 0.
	tmp = config.get('names','galaxy').split(' ')
	dataname = tmp[0]+'_'+tmp[1]+'_pixels.dat'
	f_tmp = open(dataname, 'w')
	f_tmp.write('#low freq. radio \t high freq. radio \t hybrid SFR \t spectral index \n' )
	for i in range(len(data_low)):
		f_tmp.write(str(data_low[i])) #LOFAR / low frequency
		f_tmp.write('\t')
		f_tmp.write(str(data_high[i])) #WSRT / high frequency
		f_tmp.write('\t')
		f_tmp.write(str(data_hyb[i])) #GALEX/Spitzer hybrid sfr
		f_tmp.write('\t')
		f_tmp.write(str(alpha[i])) # spectral index
		f_tmp.write('\n')
		mean += data_low[i]
	f_tmp.close()
	return float(mean)

def print_conv_data( cfg, data_radio, data_hyb, alpha, case ):
	mean = 0.
	tmp = config.get('names','galaxy').split(' ')
	dataname = tmp[0]+'_'+tmp[1]+'_pixels'+case +'.dat'
	f_tmp = open(dataname, 'w')
	f_tmp.write('#radio SFR \t conv hybrid SFR \t spectral index \n' )
	for i in range(len(data_radio)):
		f_tmp.write(str(data_radio[i])) 
		f_tmp.write('\t')
		f_tmp.write(str(data_hyb[i])) 
		f_tmp.write('\t')
		f_tmp.write(str(alpha[i])) 
		f_tmp.write('\n')
		mean += data_radio[i]
	f_tmp.close()
	return float(mean)

### Unit conversion from pixel distance to distance in kiloparsec
def convert_px2kpc(px, cfg):
	kpc_per_arcsec = cfg.getfloat('distance')*m.tan(2*m.pi / (360.* 3600.)) / 1000.
	return float(px / cfg.getfloat('pixel_per_arcsec') * kpc_per_arcsec )
### Unit conversion from distance in kiloparsec to pixel distance
def convert_kpc2px(kpc, cfg):
	kpc_per_arcsec = cfg.getfloat('distance')*m.tan(2*m.pi / (360.* 3600.)) / 1000.
	return float(kpc * cfg.getfloat('pixel_per_arcsec') / kpc_per_arcsec )

def ceil(x):
	return int(m.ceil(x))

def floor(x):
	return int(m.floor(x))

###########################
### PROGRAM STARTS HERE ###
#TODO: 
#		Format number output?!
###########################
###Read in config file from command line argument, exit if it cant be found
if(len(sys.argv) < 2):
	print 'Need filname of config file'
	sys.exit(-1)
tmp_name = str(sys.argv[1])
if(len(sys.argv)>2):
	if (sys.argv[2] == '1'):
		PRINTALL = True
print 'Full output is set to: ', PRINTALL, ' (change through command line)'

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

cmd = ''

boxsize = 	(str(config.getint('values','n_boxes') *
			config.getfloat('values','pixels_per_box')))
center_x = config.get('values','center_x')
center_y = config.get('values','center_y')
for fname in (['radio_low','radio_high','hybrid']):
	oname = config.get('names',fname).rstrip('.fits') + '_ds9box.png'
	cmd =	('ds9 '+ config.get('names',fname) +' -regions system image  -regions command "box '+
			str(center_x) +' '+ str(center_y) +' '+str(boxsize)+' '+str(boxsize)+
			' # color=yellow" -zoom to fit -scale mode 99.5  -saveimage png '+
			oname +' -quit')
	if( PRINTALL==True ):
		os.system(cmd)
		
cmd = cmd.rstrip('-quit')


###read in the three maps, 2 radio maps with different frequencys and the SFR map from galex/spitzer
data_l = read_fits(config.get('names','radio_low'))
data_h = read_fits(config.get('names','radio_high'))
data_s = read_fits(config.get('names','hybrid'))

###Call the convert function 3 times to convert each map to 1200pc^2 pixels
pixels2d_l = convert1200_adv( data_l, config['values'] )
pixels2d_h = convert1200_adv( data_h, config['values'] )
pixels2d_s = convert1200_adv( data_s, config['values'] )

###Flatten the data to 1D lists
pixels_l = flatten ( pixels2d_l )
pixels_h = flatten ( pixels2d_h )
pixels_s = flatten ( pixels2d_s )

###Apply the cuts to the data based on the signal
###create three empty lists to store the 3sigma cut data
pix_l_cut = []
pix_h_cut = []
pix_s_cut = []

#3 sigma cutoff for all datasets
for i in range(len(pixels_l)):
	if(pixels_l[i] > 3* config.getfloat('values','sigma_low')):
		if(pixels_h[i] > 3.* config.getfloat('values','sigma_high')): 
			if(pixels_s[i] > 3.* config.getfloat('values','sigma_hybrid')): 
				pix_l_cut.append(pixels_l[i])
				pix_h_cut.append(pixels_h[i])
				pix_s_cut.append(pixels_s[i])

###Define alpha array to store the spectral index (assuming a simple power law)
alpha = []
for i in range(len(pix_l_cut)):
	alpha.append( m.log10( pix_l_cut[i] /pix_h_cut[i] ) / m.log10( config.getfloat('values','freq_low') / config.getfloat('values','freq_high') )) 

#Some diagnostic output of how many points were cut
print 'Number of points in sample:', len( pixels_l ) 
print 'Number of points after 3sigma cutoff:', len( pix_l_cut )

###Print all pixel data to file (after the 3 sigma cut and before conversion to SFR)
mean = print_data( config['names'], pix_l_cut, pix_h_cut, pix_s_cut, alpha )
mean_old = 0
tmp = flatten( data_l )
for i in range(len(tmp)):
	mean_old += tmp[i]
mean = 0
for i in range(len(pixels_l)):
	mean += pixels_l[i]


mean *=  config.getfloat('values','pixels_per_box')**2
print 'Total total sum of LOFAR image:\t', '%0.3f' % mean_old
print 'Total sum in the cut LOFAR map is:', '%0.3f' % mean

#beamsize = 1.133 * config.getfloat('values','FWHM')**2 * config.getfloat('values','pixel_per_arcsec')**2
#flux = mean / beamsize 

#print 'Total flux in the cut LOFAR map is:', flux

###Now convert lofar and swrt to SFR sufrace density using condon-relation
pix_l_cut = condon( pix_l_cut, config.getfloat('values','FWHM'), config.getfloat('values','freq_low') )
pix_h_cut = condon( pix_h_cut, config.getfloat('values','FWHM'), config.getfloat('values','freq_high') )

### Diganosis only: Create square rooted map of the 1200pc data (to check if the cutting worked as intended)
if( PRINTALL == True ):
	map_sqrt(pixels2d_l, config, 'radio_low')
	map_sqrt(pixels2d_h, config, 'radio_high')
	map_sqrt(pixels2d_s, config, 'hybrid')


###Fitting for both datasets
#Fitting low /LOFAR
a_l , a_l_err , b_l, b_l_err, chi_l = fit_odr(pix_s_cut, pix_l_cut, output=True)
#Fitting high / WSRT
a_h , a_h_err , b_h, b_h_err, chi_h = fit_odr(pix_s_cut, pix_h_cut, output=True)

############
### This concludes the creation of the 'simple' pixel plots (without the actual plotting)
### Below the data for the convolved pixel plots is calculated
############

###Finding optimal gaussian kernel for both radio maps
### Calculating pixel values and fits based on the optimal kernel
#low frequency radio map
print 'Finding optimal gaussian kernel for LOFAR (low) data. This may take a moment...'
optimal_sigma_l = optimize.fsolve(fct_gauss_fit, config.getfloat('values','sigma_conv'), args=(data_s, pixels_l, pixels_h, config, 'low' ), maxfev = 10 )

conv_pix_cut_low, conv_pix_l_cut, conv_alpha_l, a_smooth_l, b_smooth_l = fct_gauss(optimal_sigma_l[0], data_s, pixels_l, pixels_h, config, 'low' )

fit_odr(conv_pix_cut_low, conv_pix_l_cut, output=True)

#Diffusion length is the FWHM/2 of the final image, FWHM are square added first
optimal_sigma_l[0] = m.sqrt( m.pow(2.3548 *optimal_sigma_l[0],2) + m.pow(1.2,2) )/2.
print 'Final value for Diffusion length:\t','%0.3f' % optimal_sigma_l[0], 'kpc'

#high frequency radio map
print 'Finding optimal gaussian kernel for WSRT (high) data. This may take a moment...'
optimal_sigma_h = optimize.fsolve(fct_gauss_fit, config.getfloat('values','sigma_conv'), args=(data_s, pixels_l, pixels_h, config, 'high' ), maxfev = 10 )

conv_pix_cut_high, conv_pix_h_cut, conv_alpha_h, a_smooth_h, b_smooth_h = fct_gauss(optimal_sigma_h[0], data_s, pixels_l, pixels_h, config, 'high' )

fit_odr(conv_pix_cut_high, conv_pix_h_cut, output=True)

#Diffusion length is the FWHM/2 of the final image, FWHM are square added first
optimal_sigma_h[0] = m.sqrt( m.pow(2.3548 *optimal_sigma_h[0],2) + m.pow(1.2,2) )/2.
print 'Final value for Diffusion length:\t', '%0.3f' % optimal_sigma_h[0], 'kpc'

############
###Finally creating the pixel plots for all the data 
############
if(PRINTALL==True):
	try:
		os.system('gnuplot fit_test.gp')
	except:
		pass

#Plotting low / LOFAR
plot(pix_s_cut, pix_l_cut, alpha, a_l, b_l, config['names'],config['boundaries'] , 'radio_low' )
#Plotting high / WSRT
plot(pix_s_cut, pix_h_cut, alpha, a_h, b_h, config['names'],config['boundaries'], 'radio_high' )
#Plotting convolved data (low freq.)
plot(conv_pix_cut_low, conv_pix_l_cut, conv_alpha_l, a_smooth_l, b_smooth_l, config['names'],config['boundaries'], 'conv_low' , optimal_sigma_l[0])
#Plotting convolved data (high freq.)
plot(conv_pix_cut_high, conv_pix_h_cut, conv_alpha_h, a_smooth_h, b_smooth_h, config['names'],config['boundaries'], 'conv_high' , optimal_sigma_h[0])

###Print all convolved pixel data to file (after the 3 sigma cut and conversion to SFR)
mean = print_conv_data( config['names'],conv_pix_cut_low, conv_pix_l_cut , conv_alpha_l ,'_conv_low' )
mean = print_conv_data( config['names'], conv_pix_cut_high, conv_pix_h_cut, conv_alpha_h, '_conv_high' )

### Write the final results (fits and diffusion lengths) to file
tmp = config.get('names','galaxy').split(' ')
results_ini = tmp[0]+'_'+tmp[1]+'_results.ini'
res_out = configparser.ConfigParser()


res_out['name'] = {'fname': config.get('names','galaxy'),
					'ds9cmd' : cmd}

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
