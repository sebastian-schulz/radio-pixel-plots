#Library for python plotting
import matplotlib.pyplot as plt
#Library to understand the fits format
from astropy.io import fits

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
		#print 'Data type and shape of', fname ,':',type(data), data.shape 
		#Now plot all the map (compare to fits file if something doesnt work)
		plt.imshow(data, cmap='gray')
		plt.colorbar()
		plt.ylim(0,1024)
		fname_image = fname.rstrip('.fits') + '.png'
		plt.savefig( fname_image )
		plt.clf()
	return data

###Prints data arrays to file as tab separated values
def print_data( cfg, data_low, data_high, data_hyb, alpha):
	mean = 0.
	tmp = cfg.get('names','galaxy').split(' ')
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

###Similar to print data, but for the convolved data
def print_conv_data( cfg, data_radio, data_hyb, alpha, case ):
	mean = 0.
	tmp = cfg.get('names','galaxy').split(' ')
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




