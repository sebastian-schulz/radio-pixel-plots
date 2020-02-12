#Pythons math package with all the constants and functions
import math as m
#Library includes a zero-finder method (to find the best value for the kernel), also includes methods fit data/curves
from scipy import optimize
#Mighty numerical library of ptyhon
import numpy as np
#Library to perform linear regression (with x and y errors!)
from scipy import odr

### Main fitting function, calls one of the 2 possible ones based on case
def fit(val_x, val_y, val_x_err=None, val_y_err=None, output=False, case='None'):
	if(case=='odr'):
		return fit_odr (val_x, val_y, val_x_err=val_x_err, val_y_err=val_y_err, outp=output)
	else:
<<<<<<< HEAD
		print('Error in fitting function call!')
=======
		print('Error in fitting function call! (probably tried to call an depreceated fitting method)')
>>>>>>> 39e6c05434ac1f792f6c59921ecae3e02b27653c
		quit()

###Fitting function with x and y errors see Python Scipy ODR
### fits a *x +b =y
### needs x and y values as list or array
def fit_odr (val_x, val_y, val_x_err=None, val_y_err=None, outp=False):
	#First sort the data into different arrays based on their spectral index
	val_log_x = []
	val_log_y = []
	val_log_x_err = []
	val_log_y_err = []
	for i in range(len(val_x)):
		if (val_x[i] < 0 ):
			val_log_x.append(0.)
		else:
			val_log_x.append( m.log10( val_x[i] ) )
		if( val_y[i] < 0 ):
			val_log_y.append(0.)
		else:
			val_log_y.append( m.log10( val_y[i] ) )

	#Define the odr fitting function and then start the fit
	linear = odr.Model(fct_odr)
	
	if( val_x_err and val_y_err ):
		for i in range(len(val_x_err)):
			val_log_x_err.append( m.log10( val_x_err[i] ) )
			val_log_y_err.append( m.log10( val_y_err[i] ) )
		mydata = odr.RealData( val_log_x, val_log_y, sx=val_log_x_err ,sy=val_log_y_err)
	else:
		print('No errors defined!')
		mydata = odr.Data( val_log_x, val_log_y )
	
	myodr = odr.ODR(mydata, linear, beta0=[1., 0.])
	outodr = myodr.run()

	#Print results to screen and then return all the relevant data (values, std errors and chi squared)
	if(outp == True):
		print('########## FIT RESULTS ###########')
		print('a =\t', '%0.3f' % outodr.beta[0], '+/-\t', '%0.3f' % outodr.sd_beta[0])
		print('b =\t','%0.3f' %  outodr.beta[1], '+/-\t', '%0.3f' % outodr.sd_beta[1])
		print('Chi squared = \t','%0.3f' %  outodr.sum_square)
	return outodr.beta[0], outodr.sd_beta[0], outodr.beta[1] , outodr.sd_beta[1], outodr.sum_square

### 2nd Fitting function, this one cannot deal with x-errors, uses LM algorithm
def fit_lsq (val_x, val_y, val_y_err=None, outp=False):
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
	#Fit with the LM-algorithm from scipy.optimize
	popt, pcov = optimize.curve_fit( fct_lsq, val_log_x, val_log_y, absolute_sigma=False, method='lm')
	perr = np.sqrt(np.diag(pcov))
	#The chisquared needs to be calculated, as the curve_fit does not have an option to do so automatically
	chisq = 0
	for i in range(len(val_log_x)):
		chisq += m.pow( fct_lsq(val_log_x[i], popt[0], popt[1]) - val_log_y[i] ,2)
	
	if(outp == True):
		print('########## FIT RESULTS ###########')
		print('a =\t', '%0.3f' % popt[0], '+/-\t', '%0.3f' % perr[0])
		print('b =\t','%0.3f' %  popt[1], '+/-\t', '%0.3f' % perr[1])
		print('Chi squared = \t','%0.3f' % chisq)
	return popt[0], perr[0], popt[1] , perr[1], chisq

###Linear model for ODR-fit
def fct_odr(B, val_x):
	return B[0]*val_x+B[1]





