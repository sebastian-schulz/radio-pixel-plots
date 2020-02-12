#System tools (file in and output ...)
import sys

#Ability to call external programs, shell commands etc
import os
from subprocess import call

#Pythons math package with all the constants and functions
import math as m
import numpy as np
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

#The fitting function supported other fitting methods earlier, this is no longer the case! (too much work maintaining the different fitting methods). DO NOT CHANGE
FIT_METHOD = 'odr'

#rel. calibration error for both radio and sfr maps
CALIB_ERR = 0.05

#################
### FUNCTIONS ###
#################

#For more information on what those functions do, check their respective files for documentation
from plotting import condon, fct_f, fct_result, map_sqrt, plot
from conversion import ceil, floor, convert_resolution_adv, conv_px_per_box, calculate_rms
from convolution import convolve_gauss, fct_gauss, fct_gauss_fit, convert_kpc2px, convert_px2kpc, flatten, calc_error
from data import print_conv_data, print_data, read_fits
from fitting import fct_odr, fit_lsq, fit_odr, fit #fct_lsq

###########################
### PROGRAM STARTS HERE ###
###########################

###TODO###
# Change convolution and regular fit not to include outliers (outside alpha range)
# BUT still include it in plots (as outliers)
##########

###Read in config file from command line argument, exit if it cant be found
if(len(sys.argv) == 1):
    print('Need filname of config file')
    sys.exit(-1)
elif(len(sys.argv)==2):
    tmp_name = str(sys.argv[1])
elif(len(sys.argv)==3):
    tmp_name = str(sys.argv[1])
    if (sys.argv[2] == '0'):
        PRINTALL = False
else:
    print('Too many command line arguments! Quitting...')
    sys.exit(-1)
print('Full output is set to: ', PRINTALL, ' (change through command line)')

#Read config file, object config works almost like a dictoinary, print everything
config = configparser.ConfigParser()
config.read(tmp_name)
if( PRINTALL == True ):
    print(config.sections())
    for key in config['names']: print(key, ' = ',config['names'][key] )
    for key in config['values']: print(key, ' = ',config['values'][key])

###change working directory to data-path (specified in the ini file)
os.chdir(config['names']['fullpath'])
print('Current working directory is:', os.getcwd())

###Create png versions of the full fits files, including the (cutting) box using ds9
###also includes the 'rmsbox' to calculate the sigmas for the 3 sigma cutoff
boxsize = str(config.getint('values','n_boxes') * conv_px_per_box( config['values'] ) )
cmd = {}
center_x = config.get('values','center_x')
center_y = config.get('values','center_y')
for fname in (['low','high','sfr']):
    oname = config.get('names',fname).rstrip('.fits') + '_ds9box.png'
    cmd[fname] =    ('ds9 '+ config.get('names',fname) +
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
pixels2d_l = convert_resolution_adv( data_l, config['values'] )
pixels2d_h = convert_resolution_adv( data_h, config['values'] )
pixels2d_s = convert_resolution_adv( data_s, config['values'] )

###Flatten the data to 1D lists
pixels_l = flatten ( pixels2d_l )
pixels_h = flatten ( pixels2d_h )
pixels_s = flatten ( pixels2d_s )

###Define alpha array to store the spectral index (assuming a simple power law)
alpha_tmp = []
for i in range(len(pixels_l)):
    alpha_tmp.append( m.log10( m.fabs(pixels_l[i] /pixels_h[i]) ) / m.log10( config.getfloat('values','freq_low') / config.getfloat('values','freq_high') )) 

###Apply the cuts to the data based on the signal
###create three empty lists to store the 3 sigma cut data
pix_l_cut = []
pix_l_cut_err = []
pix_h_cut = []
pix_h_cut_err = []
pix_s_cut = []
pix_s_cut_err = []
alpha_cut = []

pix_l_fit = []
pix_l_fit_err = []
pix_h_fit = []
pix_h_fit_err = []
pix_s_fit = []
pix_s_fit_err = []
alpha_fit = []
#calculate rms from the boxes set in the config file
sigma_high = calculate_rms( data_h, config['high_cutoff_box'])
sigma_low = calculate_rms( data_l, config['low_cutoff_box'])
sigma_sfr = calculate_rms( data_s, config['sfr_cutoff_box'])

#Use these lines if you want to change the cutoff directly *quick and dirty*
#sigma_high /= 10.
#sigma_low /= 10.
#sigma_sfr /= 10.



sigma = {'low' : sigma_low, 'high' : sigma_high, 'sfr' : sigma_sfr }
#3 sigma cutoff for all datasets
# for each dataset, there is a corresponding set containing the errors
# name is always the same with _err attached
for i in range(len(pixels_l)):
    if(pixels_l[i] > 3. * sigma_low ):
        if(pixels_h[i] > 3. * sigma_high): 
            if(pixels_s[i] > 3. * sigma_sfr ):
                    alpha_cut.append(alpha_tmp[i])
                    pix_l_cut.append(pixels_l[i])
                    pix_l_cut_err.append( calc_error( pixels_l[i],sigma_low ,CALIB_ERR ) )
                    pix_h_cut.append(pixels_h[i])
                    pix_h_cut_err.append( calc_error( pixels_h[i],sigma_high ,CALIB_ERR ) )
                    pix_s_cut.append(pixels_s[i])
                    pix_s_cut_err.append( calc_error( pixels_s[i],sigma_sfr ,CALIB_ERR ) )
                    if(alpha_tmp[i] < config.getfloat('boundaries','low') ):
                        alpha_fit.append(alpha_tmp[i])
                        pix_l_fit.append(pixels_l[i])
                        pix_l_fit_err.append( calc_error( pixels_l[i],sigma_low ,CALIB_ERR ) )
                        pix_h_fit.append(pixels_h[i])
                        pix_l_cut_err.append( calc_error( pixels_h[i],sigma_high ,CALIB_ERR ) )
                        pix_s_fit.append(pixels_s[i])
                        pix_s_fit_err.append( calc_error( pixels_s[i],sigma_sfr ,CALIB_ERR ) )

print('RMS for the 3 different maps, used as sigma for 3 sigma cutoff:\n', sigma)

#Some diagnostic output of how many points were cut
print('Number of points in sample:', len( pixels_l ))
print('Number of points after 3 sigma cutoff:', len( pix_l_cut ))

###Print all pixel data to file (after the 3 sigma cut and before conversion to SFR)
mean = print_data( config, pix_l_cut, pix_l_cut_err, pix_h_cut, pix_h_cut_err, pix_s_cut, pix_s_cut_err, alpha_cut )
mean_old = 0
tmp = flatten( data_l )
for i in range(len(tmp)):
    mean_old += m.fabs(tmp[i])
mean = 0
for i in range(len(pixels_l)):
    mean += pixels_l[i]


mean *=  conv_px_per_box( config['values'] )**2
print('Total total sum of lower freq. image:\t', '%0.3f' % mean_old)
print('Total sum in the cut lower freq. map is:', '%0.3f' % mean)

###Now convert lofar and swrt to SFR sufrace density using condon-relation
pix_l_cut = condon( pix_l_cut,
                   config.getfloat('values','FWHM'),
                   config.getfloat('values','freq_low') )

pix_h_cut = condon( pix_h_cut,
                   config.getfloat('values','FWHM'),
                   config.getfloat('values','freq_high') )

pix_l_fit = condon( pix_l_fit, config.getfloat('values','FWHM'),
                   config.getfloat('values','freq_low') )

pix_h_fit = condon( pix_h_fit, config.getfloat('values','FWHM'),
                   config.getfloat('values','freq_high') )

#The same conversion is performed on the errors
pix_l_cut_err = condon( pix_l_cut_err,
                       config.getfloat('values','FWHM'),
                       config.getfloat('values','freq_low') )

pix_h_cut_err = condon( pix_h_cut_err,
                       config.getfloat('values','FWHM'),
                       config.getfloat('values','freq_high') )

pix_l_fit_err = condon( pix_l_fit_err,
                       config.getfloat('values','FWHM'),
                       config.getfloat('values','freq_low') )

pix_h_fit_err = condon( pix_h_fit_err,
                       config.getfloat('values','FWHM'),
                       config.getfloat('values','freq_high') )


### Diganosis only: Create square rooted map of the 1200pc data (to check if the cutting worked as intended)
if( PRINTALL == True ):
    map_sqrt(pixels2d_l, config, 'low')
    map_sqrt(pixels2d_h, config, 'high')
    map_sqrt(pixels2d_s, config, 'sfr')

###Fitting for both datasets
#Fitting low /LOFAR
a_l , a_l_err , b_l, b_l_err, chi_l = fit(    pix_s_fit, 
                                            pix_l_fit,
                                            val_x_err=pix_s_fit_err,
                                            val_y_err=pix_l_fit_err,
                                            output=True,
                                            case=FIT_METHOD)
#Fitting high / WSRT
a_h , a_h_err , b_h, b_h_err, chi_h = fit(    pix_s_fit,
                                            pix_h_fit,
                                            val_x_err=pix_s_fit_err,
                                            val_y_err=pix_h_fit_err,
                                            output=True,
                                            case=FIT_METHOD)

############
### This concludes the creation of the 'simple' pixel plots (without the actual plotting)
### Below the data for the convolved pixel plots is calculated
############

###Finding optimal gaussian kernel for both radio maps
### Calculating pixel values and fits based on the optimal kernel

###low frequency radio map
print('Finding optimal gaussian kernel for lower freqency data. This may take a moment...')
optimal_sigma_l = optimize.fsolve(fct_gauss_fit,
                                  config.getfloat('values','sigma_conv'),
                                  args=(data_s, pixels_l, pixels_h, sigma , config, 'low', PRINTALL, CALIB_ERR, FIT_METHOD),
                                  maxfev = 20 )

conv_pix_cut_low,\
conv_pix_cut_low_err,\
conv_pix_l_cut,\
conv_pix_l_cut_err,\
conv_alpha_l,\
a_smooth_l,\
b_smooth_l = fct_gauss( optimal_sigma_l[0],
                         data_s, pixels_l,
                         pixels_h,
                         sigma,
                         config,
                         'low',
                         PRINTALL,
                         CALIB_ERR,
                         FIT_METHOD )

_, a_l_conv_err, _, _, _ = fit(conv_pix_cut_low,
                               conv_pix_l_cut,
                               val_x_err=conv_pix_cut_low_err,
                               val_y_err=conv_pix_l_cut_err,
                               output=True,
                               case=FIT_METHOD)

#Diffusion length is the FWHM/2 of the final image, FWHM are square added first
optimal_sigma_l[0] = m.sqrt( m.pow(2.3548 *optimal_sigma_l[0],2) + m.pow(1.2,2) )/2.
print('Final value for Diffusion length:\t','%0.3f' % optimal_sigma_l[0], 'kpc')

###high frequency radio map
print('Finding optimal gaussian kernel for higher frequency data. This may take a moment...')
optimal_sigma_h = optimize.fsolve(fct_gauss_fit,
                                  config.getfloat('values','sigma_conv'),
                                  args=(data_s, pixels_l, pixels_h, sigma , config, 'high', PRINTALL, CALIB_ERR, FIT_METHOD ),
                                  maxfev = 20 )

conv_pix_cut_high,\
conv_pix_cut_high_err,\
conv_pix_h_cut, conv_pix_h_cut_err,\
conv_alpha_h,\
a_smooth_h,\
b_smooth_h = fct_gauss(optimal_sigma_h[0],
                       data_s, pixels_l,
                       pixels_h,
                       sigma,
                       config,
                       'high',
                       PRINTALL,
                       CALIB_ERR,
                       FIT_METHOD )

_, a_h_conv_err, _, _, _ = fit(conv_pix_cut_high,
                               conv_pix_h_cut,
                               val_x_err=conv_pix_cut_high_err,
                               val_y_err=conv_pix_h_cut_err,
                               output=True,
                               case=FIT_METHOD)

#Diffusion length is the FWHM/2 of the final image, FWHM are square added first
optimal_sigma_h[0] = m.sqrt( m.pow(2.3548 *optimal_sigma_h[0],2) + m.pow(1.2,2) )/2.
print('Final value for Diffusion length:\t', '%0.3f' % optimal_sigma_h[0], 'kpc')

############
###Finally creating the pixel plots for all the data 
############
print('Creating final images and writing results to file ...')

#Plotting low / LOFAR
plot(pix_s_cut,
     pix_l_cut,
     alpha_cut,
     a_l,
     b_l,
     config,
     'low',
     x_err=pix_s_cut_err,
     y_err=pix_l_cut_err )
#Plotting high / WSRT
plot(pix_s_cut,
     pix_h_cut,
     alpha_cut,
     a_h,
     b_h,
     config,
     'high',
     x_err=pix_s_cut_err,
     y_err=pix_h_cut_err )
#Plotting convolved data (low freq.)
plot(conv_pix_cut_low,
     conv_pix_l_cut,
     conv_alpha_l,
     a_smooth_l,
     b_smooth_l,
     config,
     'conv_low',
     optimal_sigma_l[0],
     x_err=conv_pix_cut_low_err,
     y_err=conv_pix_l_cut_err)
#Plotting convolved data (high freq.)
plot(conv_pix_cut_high,
     conv_pix_h_cut,
     conv_alpha_h,
     a_smooth_h,
     b_smooth_h,
     config,
     'conv_high',
     optimal_sigma_h[0],
     x_err=conv_pix_cut_high_err,
     y_err=conv_pix_h_cut_err)

###Print all convolved pixel data to file (after the 3 sigma cut and conversion to SFR)
mean = print_conv_data( config,
                       conv_pix_cut_low,
                       conv_pix_cut_low_err,
                       conv_pix_l_cut,
                       conv_pix_l_cut_err,
                       conv_alpha_l,
                       'conv_low')

mean = print_conv_data( config,
                       conv_pix_cut_high,
                       conv_pix_cut_high_err,
                       conv_pix_h_cut,
                       conv_pix_h_cut_err,
                       conv_alpha_h,
                       'conv_high')

### Write the final results (fits and diffusion lengths) to file
# this makes use of the configparser library again 
tmp = config.get('names','galaxy').split(' ')
results_ini = tmp[0]+'_'+tmp[1]+'_results.ini'
res_out = configparser.ConfigParser()

res_out['name'] = {'fname': config.get('names','galaxy'),
                    'ds9cmd_low' : cmd['low'],
                    'ds9cmd_high' : cmd['high'],
                    'ds9cmd_sfr' : cmd['sfr']}

res_out['Low_freq_fit'] =    {'#Fit results for pixel plots, with std errors\n'
                            'a': str(a_l),
                            'a_err': str(a_l_err),
                            'b': str(b_l_err),
                            'b_err': str(b_l_err),
                            'chi_sqr': str(chi_l)}

res_out['High_freq_fit'] =    {'a': str(a_h),
                            'a_err': str(a_h_err),
                            'b': str(b_h_err),
                            'b_err': str(b_h_err),
                            'chi_sqr': str(chi_h)}

res_out['Conv_results'] =    {'#Diffusion length from gaussian kernel in kpc \n'
                            'sigma_l': str(optimal_sigma_l[0]),
                            'a_err_l': str(a_l_conv_err),
                            'a_err_h': str(a_h_conv_err),
                            'sigma_h': str(optimal_sigma_h[0])}
                            
res_out['rms_sigma'] =         {'#The RMS values from calculated from the rms boxes\n'
                            'rms_l': str(sigma_low),
                            'rms_h': str(sigma_high),
                            'rms_s': str(sigma_sfr)}

res_out['frequencies'] =         {'#Frequencies for the two radio maps\n'
                            'freq_low': config.get('values','freq_low'),
                            'freq_high': config.get('values','freq_high')}


with open(results_ini, 'w') as configfile:
    res_out.write(configfile)

tmp = config.get('names','galaxy').split(' ')
tmp2 = tmp[0]+'_'+tmp[1]+'_plots_combined.pdf'

os.system('rm *plots_combined.pdf')
os.system('pdfunite *pdf '+ tmp2 )
print("Finished!")





