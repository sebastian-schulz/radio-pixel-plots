# Mighty numerical library of ptyhon
import numpy as np
# Pythons math package with all the constants and functions
import math as m
# Library to perform linear regression (with x and y errors!)
from scipy import odr


# PASSED BUGTEST (July 2018), gives the same results as DS9 with the same boxes
def convert_resolution_adv(data, x_center, y_center, n_boxes, px_per_box):
    # Idea: allow for pixels_per_box to be float!
    # assumption 1 : middle pixel for p and p' have the same position (center)
    # calculate everything in p distances: (n_boxes should be odd!)
    # bottom left x = center_x - px_per_box * n_boxes/2
    # same for y
    # loop over all p'
    # 	loop over p inside p'
    # 		go to bottom left corner,
    # 		find pixel p fully inside,
    # 		same for top right corner,
    # 		-> calculate sum of all pixels fully inside,
    # 		now borders:
    # 		4 corners
    # 		2 edges horizontal
    # 		2 edges vertical

    # position of center of bottom left new pixel in old pixels
    # for some reason x and y axis are in the wrong order in the array
    img_bot_left_y = x_center - 1 - (n_boxes / 2. - 0.5) * px_per_box
    img_bot_left_x = y_center - 1 - (n_boxes / 2. - 0.5) * px_per_box
    # new pixel array
    n = 0
    n_boxes = int(n_boxes)
    pixels = np.zeros((n_boxes, n_boxes))

    for x in range(n_boxes):
        for y in range(n_boxes):
            # position of current pixel in old pixels
            img_current_x = img_bot_left_x + x * px_per_box
            img_current_y = img_bot_left_y + y * px_per_box
            # positions of the corners of the new pixels in old pixels
            pix_x_min = img_current_x - px_per_box / 2.
            pix_x_max = img_current_x + px_per_box / 2.
            pix_y_min = img_current_y - px_per_box / 2.
            pix_y_max = img_current_y + px_per_box / 2.
            # loop over the old pixels fully inside the current new pixel
            for p_x in range(ceil(pix_x_min), floor(pix_x_max), 1):
                for p_y in range(ceil(pix_y_min), floor(pix_y_max), 1):
                    pixels[x, y] += data[p_x, p_y]
                    n += 1
            # now loop over the edges, without the corners
            factor_max_x = m.fabs(floor(pix_x_max) - pix_x_max)
            factor_min_x = m.fabs(ceil(pix_x_min) - pix_x_min)
            factor_max_y = m.fabs(floor(pix_y_max) - pix_y_max)
            factor_min_y = m.fabs(ceil(pix_y_min) - pix_y_min)
            # if (factor_max_x > 1 or factor_max_y > 1 or factor_min_x > 1 or factor_min_y > 1):
            # print( 'ERROR!')
            # print('x and y factors \t', factor_max_x, '\t', factor_max_y)
            for p_y in range(ceil(pix_y_min), floor(pix_y_max), 1):
                n += factor_min_x + factor_max_x
                pixels[x, y] += data[floor(pix_x_min), p_y] * factor_min_x
                pixels[x, y] += data[ceil(pix_x_max), p_y] * factor_max_x
            for p_x in range(ceil(pix_x_min), floor(pix_x_max), 1):
                n += factor_min_y + factor_max_y
                pixels[x, y] += data[p_x, floor(pix_y_min)] * factor_min_y
                pixels[x, y] += data[p_x, ceil(pix_y_max)] * factor_max_y

            # finally the four corner pixels:
            pixels[x, y] += data[floor(pix_x_min), floor(pix_y_min)] * factor_min_x * factor_min_y
            pixels[x, y] += data[floor(pix_x_min), ceil(pix_y_max)] * factor_min_x * factor_max_y
            pixels[x, y] += data[ceil(pix_x_max), ceil(pix_y_max)] * factor_max_x * factor_max_y
            pixels[x, y] += data[ceil(pix_x_max), floor(pix_y_min)] * factor_max_x * factor_min_y
            n += factor_min_x * factor_min_y + \
                 factor_min_x * factor_max_y + \
                 factor_max_x * factor_max_y + \
                 factor_max_x * factor_min_y
    return pixels / px_per_box ** 2


# Custom defined ceiling and floor functions with int as return type
def ceil(x):
    return int(m.ceil(x))


def floor(x):
    return int(m.floor(x))


def calculate_rms(data, x_center, y_center, x_size, y_size):
    rms = 0.
    # calculate borders of the rms box in pixel, makes use of int cutting
    # therefor the box may be 1 pixel length off in any direction
    # but that only has a small impact
    y_min = y_center - y_size / 2.
    y_max = y_center + y_size / 2.
    x_min = x_center - x_size / 2.
    x_max = x_center + x_size / 2.
    n = int(x_size * y_size)
    for i in range(int(x_min), int(x_max)):
        for j in range(int(y_min), int(y_max)):
            rms += m.pow(data[j][i], 2) #IMPORTANT: MUST BE data[Y][X]
    rms = m.sqrt(rms / n)
    # print('ROOT MEAN SQUARE OF BOX: ', rms)
    return rms


# Calculate the error for each value based on background noise and
# a calibration error that depends on the magnitude of the value
def calc_rms_error(val, noise, calibration):
    error = m.sqrt(m.pow(noise, 2.) + m.pow(val * calibration, 2.))
    return error


# Calculates SFR from radio data using condons relation
# x in Jy/beam, fwhm in arcsec, freq in MHz
def condon(x, fwhm, freq):
    for i in range(len(x)):
        x[i] = x[i] * 3.31e3 * (freq / 1400.) ** .8 * fwhm ** (-2.)
    return x


def fit_odr(val_x, val_y, val_x_err=None, val_y_err=None):
    # First sort the data into different arrays based on their spectral index
    val_log_x = []
    val_log_y = []
    val_log_x_err = []
    val_log_y_err = []
    for i in range(len(val_x)):
        if val_x[i] < 0:
            val_log_x.append(0.)
        else:
            val_log_x.append(m.log10(val_x[i]))
        if val_y[i] < 0:
            val_log_y.append(0.)
        else:
            val_log_y.append(m.log10(val_y[i]))

    # Define the odr fitting function and then start the fit
    linear = odr.Model(fct_odr)

    if val_x_err and val_y_err:
        for i in range(len(val_x_err)):
            val_log_x_err.append(m.log10(val_x_err[i]))
            val_log_y_err.append(m.log10(val_y_err[i]))
        mydata = odr.RealData(val_log_x, val_log_y, sx=val_log_x_err, sy=val_log_y_err)
    else:
        print('No errors defined!')
        mydata = odr.Data(val_log_x, val_log_y)

    myodr = odr.ODR(mydata, linear, beta0=[1., 0.])
    outodr = myodr.run()

    # Print results to screen and then return all the relevant data (values, std errors and chi squared)
    print('########## FIT RESULTS ###########')
    print('a =\t', '%0.3f' % outodr.beta[0], '+/-\t', '%0.3f' % outodr.sd_beta[0])
    print('b =\t', '%0.3f' % outodr.beta[1], '+/-\t', '%0.3f' % outodr.sd_beta[1])
    print('Chi squared = \t', '%0.3f' % outodr.sum_square)

    return [outodr.beta[0], outodr.sd_beta[0]], [outodr.beta[1], outodr.sd_beta[1]], outodr.sum_square


# Linear model for ODR-fit
def fct_odr(b, val_x):
    return b[0] * val_x + b[1]


def convert_px2kpc(px, distance, px_per_arcsec):
    kpc_per_arcsec = distance * m.tan(2 * m.pi / (360. * 3600.)) / 1000.
    return float(px / px_per_arcsec * kpc_per_arcsec)


def convert_kpc2px(kpc, distance, pixel_per_arcsec):
    kpc_per_arcsec = distance * m.tan(2 * m.pi / (360. * 3600.)) / 1000.
    return float(kpc * pixel_per_arcsec / kpc_per_arcsec)


# y=x function
def fct_f(x):
    return x


# Exponential function to correctly plot the linear (log-log) fit
def fct_result(x, a, b):
    return x ** a * 10 ** b
