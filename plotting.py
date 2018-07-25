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
### only used to debug/compare with volker
def map_sqrt(pixels , cfg, galaxyname ):
	tmp = np.zeros(shape = ( len(pixels), len(pixels)) )
	for i in range(len(pixels)):
		for j in range(len(pixels)):
			tmp[i,j] = m.sqrt(m.fabs(pixels[i,j])) 
	#Plot the modified data as a greyscale map (compare with original)
	plt.imshow( tmp, cmap='gray') #clim=(0 , 0.05)
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
		elif( float(boundaries['high']) > alpha[i] ):
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
	plt.legend([boundaries['low']+r'$\,> {\alpha}\,$', boundaries['low']+r'$\,> {\alpha}>\,$' +boundaries['med'] ,  r'${\alpha} >\,$' + boundaries['high'] , r'{Condon}', r'{Least Square Fit}']) #r'{Gnuplot fit}', r'{Outlier}'
	plt.grid(True)
	plt.ylabel(r'{SFR surface density based on Radio data}')
	plt.xlabel(r'{SFR surface density based on GALEX/Spitzer}')
	plt.xlim(1e-4,1e-1)
	plt.ylim(1e-4,1e-1)
	#Save plots as pdf
	if(case == 'high' or case == 'low'):
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




