#Library for python plotting
import matplotlib.pyplot as plt
#Mighty numerical library of ptyhon
import numpy as np
#Pythons math package with all the constants and functions
import math as m

#For LateX fonts and symbols in plots
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['T1']})# T1 is the LaTex standard font, also used by scr familiy documents as default
rc('text', usetex=True)

###Prints a map of the 1200pc pixel data (with sqrt(values) for better visibility)
### used for debug purposes
def map_sqrt(pixels , cfg, galaxyname ):
	tmp = np.zeros(shape = ( len(pixels), len(pixels)) )
	for i in range(len(pixels)):
		for j in range(len(pixels)):
			tmp[i,j] = m.sqrt(m.fabs(pixels[i,j])) 
	#Plot the modified data as a greyscale map (compare with original)
	plt.imshow( tmp, cmap='gray')# clim=(edges[9]) )
	plt.colorbar()
	plt.ylim(0, cfg.getint('values','n_boxes') - 1 )
	image_name = cfg.get('names', galaxyname )
	image_name = image_name.rstrip('.fits')
	image_name += '_sqrtmap.png'
	#And save the plot to file:
	plt.savefig( image_name )
	plt.clf() #clears the plot for further plotting

### Calculates SFR from radio data using condons relation
### x in Jy/beam, fwhm in arcsec, freq in MHz
def condon( x, fwhm, freq ):
	for i in range(len(x)):
		x[i] = x[i] * 3.31e3 * (freq / 1400.)**(.8) * (fwhm)**(-2.)
	return x

### y=x function
def fct_f (x):
	return x

###Exponential function to correctly plot the linear (log-log) fit
def fct_result(x, a, b):
	return x**a * 10**b

###Plotting function, also sorts data according to spectral index 
def plot(val_x, val_y, alpha, a, b, cfg, case, sigma=None, x_err=None, y_err=None ):
	#Define empty arrys/lists to store portions of data based on different spectral indices
	#add in support for plots with error-bars both in x and y, maybe both?
	val_flat_x = []
	val_flat_x_err = []
	val_flat_y = []
	val_flat_y_err = []
	val_med_x = []
	val_med_x_err = []
	val_med_y = []
	val_med_y_err = []
	val_steep_x = []
	val_steep_x_err = []
	val_steep_y = []
	val_steep_y_err = []
	val_outlier_x = []
	val_outlier_x_err = []
	val_outlier_y = []
	val_outlier_y_err = []
	#Sort the data into different ranges of alpha
	for i in range(len(val_x)):
		if( cfg.getfloat('boundaries','low') > alpha[i] > cfg.getfloat('boundaries','med') ):
			val_flat_x.append( val_x[i] )
			val_flat_x_err.append( x_err[i] )
			val_flat_y.append( val_y[i] )
			val_flat_y_err.append( y_err[i] )
		elif( cfg.getfloat('boundaries','med') > alpha[i] > cfg.getfloat('boundaries','high') ):
			val_med_x.append( val_x[i] )
			val_med_x_err.append( x_err[i] )
			val_med_y.append( val_y[i] )
			val_med_y_err.append( y_err[i] )
		elif( cfg.getfloat('boundaries','high') > alpha[i] ):
			val_steep_x.append( val_x[i] )
			val_steep_x_err.append( x_err[i] )
			val_steep_y.append( val_y[i] )
			val_steep_y_err.append( y_err[i] )
		else:
			val_outlier_x.append( val_x[i] )
			val_outlier_x_err.append( x_err[i] )
			val_outlier_y.append( val_y[i] )
			val_outlier_y_err.append( y_err[i] )
	
	#Create the plot: 
	t = np.linspace(1e-5,1e-1) #datapoints in the plotting range
	plt.clf() #clear pervious plots
	
	#defining labels
	l1 = cfg.get('boundaries','high') + r'$\,>{\alpha}$'
	l2 = cfg.get('boundaries','med') + r'$\,> {\alpha}>\,$' + cfg.get('boundaries','high')
	l3 = cfg.get('boundaries','low') + r'$\,> {\alpha}>\,$' + cfg.get('boundaries','med')
	l6 = r'{Outliers}'
	l4 = '{Condon}'
	l5 = r'{Least Square Fit}'

	ax = plt.subplot(111)
	#double-logarithmic plots with advanced options (color, marker, ...)
	ax.errorbar( val_steep_x , val_steep_y, xerr=val_steep_x_err, yerr=val_steep_y_err, marker='.', linestyle='None', color='b', label=l1)
	ax.errorbar( val_med_x , val_med_y, xerr=val_med_x_err, yerr=val_med_x_err, marker='.', linestyle='None', color='g', label=l2)
	ax.errorbar( val_flat_x, val_flat_y, xerr=val_flat_x_err, yerr= val_flat_y_err, marker='.', linestyle='None', color='r', label=l3)
	ax.plot( t, fct_f(t) ,linestyle='--', label=l4)
	ax.plot(t, fct_result(t, a, b) , linestyle='-', label=l5)
	ax.errorbar( val_outlier_x, val_outlier_y, xerr=val_outlier_x_err, yerr=val_outlier_y_err, marker='.', linestyle='None', color='tab:gray', label=l6)
	ax.legend()

	ax.set_xscale("log", nonposx='clip')
	ax.set_yscale("log", nonposy='clip')

	ax.grid(True)
	ax.set_ylabel(r'{$\left(\Sigma_{SFR}\right) _{RC}$ in M$_{Sun}$ yr$^{-1}$ kpc$^{-2}$}')
	ax.set_xlabel(r'{$\left(\Sigma_{SFR}\right) _{hyb}$ in M$_{Sun}$ yr$^{-1}$ kpc$^{-2}$}')
	#plt.xlim(1e-4,1e-1)
	#plt.ylim(1e-4,1e-1)
	#Save plots as pdf
	if(case == 'high' or case == 'low'):
		outfile = cfg.get('names',case).rstrip('.fits') + '_pixel.pdf'
	else:
		outfile = cfg.get('names',case) + '_pixel.pdf'
		tmp = cfg.get('names','galaxy').split(' ')
		outfile = tmp[0]+'_'+tmp[1]+'_'+ cfg.get('names',case) +'_pixel_conv.pdf'
	
	if(case == 'conv_low'):
		title = cfg.get('names','galaxy') + r' with Gaussian kernel, $l_{CRE} = $'+'%0.2f' % sigma + r' kpc , @' + cfg.get('values','freq_low') + r' MHz'
	elif(case == 'conv_high'):
		title = cfg.get('names','galaxy') + r' with Gaussian kernel, $l_{CRE} = $'+'%0.2f' % sigma + r' kpc, @' + cfg.get('values','freq_high') + r' MHz'
	elif(case == 'low'):
		title = cfg.get('names','galaxy') + r', @' + cfg.get('values','freq_low') + r' MHz'
	elif(case == 'high'):
		title = cfg.get('names','galaxy') + r', @' + cfg.get('values','freq_high') + r' MHz'
	
	ax.set_title( title )
	plt.savefig( outfile )
	plt.clf() #clean for further plotting
	return None




