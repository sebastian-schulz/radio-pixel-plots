#Library for python plotting
import matplotlib.pyplot as plt
#Library to understand the fits format
from astropy.io import fits
#Mighty numerical library of ptyhon
import numpy as np


###Read image from fits file and store it in 2d numpy array
def read_fits( fname , PRINTALL):
	#Reading in the image file and checking size as well as shape
	data = fits.getdata( fname )
	#if it has the shape (1,1, X, Y), SLICE it to a 2D array (x,y)
	try:
		data = data[0][0][:][:]
	except:
		pass
	if (PRINTALL == True):
		#Check for size and shape of the array
<<<<<<< HEAD
		#print('Data type and shape of', fname ,':',type(data), data.shape )
=======
		#print( 'Data type and shape of', fname ,':',type(data), data.shape ) 
>>>>>>> 39e6c05434ac1f792f6c59921ecae3e02b27653c
		#Now plot all the map (compare to fits file if something doesnt work)
		
		#first create a histogram with 1e6 bins to set the z-scale to 90%
		values, edges = np.histogram(data, bins = 1000000, range=(0, data.max()))
		total = 0
		for l in range(len(values)):
			total += values[l]
		i = 0
		cmax = 0
		while True:
			cmax +=values[i]
			if( cmax > total * .99 ):
				break
			i += 1
		#print(str(i)+'\n'+str(edges[i]))
		#Create the plot with boundaries between 0 and 99% of the histogram
		plt.imshow(data, cmap='gray', clim=( 0 , edges[i] ) )
		plt.colorbar()
		plt.ylim(0,1024)
		fname_image = fname.rstrip('.fits') + '.png'
		plt.savefig( fname_image )
		plt.clf()
	return data

###Prints data arrays to file as tab separated values
def print_data( cfg, data_low, data_low_err, data_high, data_high_err, data_hyb, data_hyb_err, alpha):
	mean = 0.
	tmp = cfg.get('names','galaxy').split(' ')
	dataname = tmp[0]+'_'+tmp[1]+'_pixel.dat'
	f_tmp = open(dataname, 'w')
	f_tmp.write('#low freq. radio \t error \t high freq. radio \t error \t hybrid SFR \t error \t spectral index \n' )
	for i in range(len(data_low)):
		f_tmp.write(str(data_low[i])) #LOFAR / low frequency
		f_tmp.write('\t')
		f_tmp.write(str(data_low_err[i])) #errors
		f_tmp.write('\t')
		f_tmp.write(str(data_high[i])) #WSRT / high frequency
		f_tmp.write('\t')
		f_tmp.write(str(data_high_err[i])) #errors
		f_tmp.write('\t')
		f_tmp.write(str(data_hyb[i])) #GALEX/Spitzer hybrid sfr
		f_tmp.write('\t')
		f_tmp.write(str(data_hyb_err[i])) #errors
		f_tmp.write('\t')
		f_tmp.write(str(alpha[i])) # spectral index
		f_tmp.write('\n')
		mean += data_low[i]
	f_tmp.close()
	return float(mean)

###Similar to print data, but for the convolved data
def print_conv_data( cfg, data_radio, data_radio_err, data_hyb, data_hyb_err, alpha, case ):
	mean = 0.
	tmp = cfg.get('names','galaxy').split(' ')
	dataname = tmp[0]+'_'+tmp[1]+'_'+ cfg.get('names',case) +'_pixel_conv.dat'
	f_tmp = open(dataname, 'w')
	f_tmp.write('#radio SFR \t error \t conv hybrid SFR \t error \t spectral index \n' )
	for i in range(len(data_radio)):
		f_tmp.write(str(data_radio[i])) 
		f_tmp.write('\t')
		f_tmp.write(str(data_radio_err[i])) 
		f_tmp.write('\t')
		f_tmp.write(str(data_hyb[i])) 
		f_tmp.write('\t')
		f_tmp.write(str(data_hyb_err[i])) 
		f_tmp.write('\t')
		f_tmp.write(str(alpha[i])) 
		f_tmp.write('\n')
		mean += data_radio[i]
	f_tmp.close()
	return float(mean)


