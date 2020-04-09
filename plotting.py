# Library for python plotting
import matplotlib.pyplot as plt
# Mighty numerical library of ptyhon
import numpy as np
import configparser

from calc_functions import fct_f, fct_result
#For LateX fonts and symbols in plots
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['T1']})# T1 is the LaTex standard font, also used by scr familiy documents as default
rc('text', usetex=True)


class Plotting:
    def __init__(self, pp, conv):
        self.pp = pp
        self.conv = conv

    def print_all(self):
        self.__print_data(self.pp)
        self.__print_conv_data_low(self.pp, self.conv)
        self.__print_conv_data_high(self.pp, self.conv)
        self.__print_results(self.pp, self.conv)

    def plot_all(self):
        self.__pplot(self.pp,
                     self.pp.pixel_sfr_cut[0],
                     self.pp.pixel_low_cut[2],
                     self.pp.alpha_cut,
                     self.pp.fit_low_a[0],
                     self.pp.fit_low_b[0],
                    'low',
                     x_err=self.pp.pixel_sfr_cut[1],
                     y_err=self.pp.pixel_low_cut[3])
        self.__pplot(self.pp,
                     self.pp.pixel_sfr_cut[0],
                     self.pp.pixel_high_cut[2],
                     self.pp.alpha_cut,
                     self.pp.fit_high_a[0],
                     self.pp.fit_high_b[0],
                    'high',
                     x_err=self.pp.pixel_sfr_cut[1],
                     y_err=self.pp.pixel_high_cut[3])

        self.__pplot(self.pp,
                     self.conv.pixel_sfr_conv_cut_low[0],
                     self.conv.pixel_low_conv_cut_low[2],
                     self.conv.alpha_conv_cut_low,
                     self.conv.fit_low_conv_a[0],
                     self.conv.fit_low_conv_b[0],
                    'conv_low',
                     sigma= self.conv.optimal_sigma_low[0],
                     x_err=self.conv.pixel_sfr_conv_cut_low[1],
                     y_err=self.conv.pixel_low_conv_cut_low[3])
        self.__pplot(self.pp,
                     self.conv.pixel_sfr_conv_cut_high[0],
                     self.conv.pixel_high_conv_cut_high[2],
                     self.conv.alpha_conv_cut_high,
                     self.conv.fit_high_conv_a[0],
                     self.conv.fit_high_conv_b[0],
                    'conv_high',
                     sigma=self.conv.optimal_sigma_high[0],
                     x_err=self.conv.pixel_sfr_conv_cut_high[1],
                     y_err=self.conv.pixel_high_conv_cut_high[3])

    # Prints data arrays to file as tab separated values
    def __print_data(self, pp):
        # mean = 0.
        tmp = pp.config.get('names', 'galaxy').split(' ')
        dataname = tmp[0] + '_' + tmp[1] + '_pixel.dat'
        f_tmp = open(dataname, 'w')
        print(dataname)
        f_tmp.write(
            '#low freq. radio \t error \t high freq. radio \t error \t hybrid SFR \t error \t spectral index \n')
        for i in range(len(pp.alpha_cut)):
            f_tmp.write(str(pp.pixel_low_cut[0][i]))  # LOFAR / low frequency
            f_tmp.write('\t')
            f_tmp.write(str(pp.pixel_low_cut[1][i]))  # errors
            f_tmp.write('\t')
            f_tmp.write(str(pp.pixel_high_cut[0][i]))  # WSRT / high frequency
            f_tmp.write('\t')
            f_tmp.write(str(pp.pixel_high_cut[1][i]))  # errors
            f_tmp.write('\t')
            f_tmp.write(str(pp.pixel_sfr_cut[0][i]))  # GALEX/Spitzer hybrid sfr
            f_tmp.write('\t')
            f_tmp.write(str(pp.pixel_sfr_cut[1][i]))  # errors
            f_tmp.write('\t')
            f_tmp.write(str(pp.alpha_cut[i]))  # spectral index
            f_tmp.write('\n')
            # mean += data_low[i]
        f_tmp.close()
        # return float(mean)

    # Similar to print data, but for the convolved data
    def __print_conv_data_high(self, pp, conv):
        # mean = 0.
        tmp = pp.config.get('names', 'galaxy').split(' ')
        dataname = tmp[0] + '_' + tmp[1] + 'high' + '_pixel_conv.dat'
        f_tmp = open(dataname, 'w')
        print(dataname)
        f_tmp.write('#radio SFR \t error \t conv hybrid SFR \t error \t spectral index \n')
        for i in range(len(conv.alpha_conv_cut_high)):
            f_tmp.write(str(conv.pixel_high_conv_cut_high[0][i]))
            f_tmp.write('\t')
            f_tmp.write(str(conv.pixel_high_conv_cut_high[1][i]))
            f_tmp.write('\t')
            f_tmp.write(str(conv.pixel_sfr_conv_cut_high[0][i]))
            f_tmp.write('\t')
            f_tmp.write(str(conv.pixel_sfr_conv_cut_high[0][i]))
            f_tmp.write('\t')
            f_tmp.write(str(conv.alpha_conv_cut_high[i]))
            f_tmp.write('\n')
            # mean += data_radio[i]
        f_tmp.close()

    def __print_conv_data_low(self, pp, conv):
        # mean = 0.
        tmp = pp.config.get('names', 'galaxy').split(' ')
        dataname = tmp[0] + '_' + tmp[1] + 'low' + '_pixel_conv.dat'
        f_tmp = open(dataname, 'w')
        print(dataname)
        f_tmp.write('#radio SFR \t error \t conv hybrid SFR \t error \t spectral index \n')
        for i in range(len(conv.alpha_conv_cut_low)):
            f_tmp.write(str(conv.pixel_low_conv_cut_low[0][i]))
            f_tmp.write('\t')
            f_tmp.write(str(conv.pixel_low_conv_cut_low[1][i]))
            f_tmp.write('\t')
            f_tmp.write(str(conv.pixel_sfr_conv_cut_low[0][i]))
            f_tmp.write('\t')
            f_tmp.write(str(conv.pixel_sfr_conv_cut_low[0][i]))
            f_tmp.write('\t')
            f_tmp.write(str(conv.alpha_conv_cut_low[i]))
            f_tmp.write('\n')
            # mean += data_radio[i]
        f_tmp.close()
        # return float(mean)

    # Plotting function, also sorts data according to spectral index
    def __pplot(self, pp, val_x, val_y, alpha, a, b, case, sigma=None, x_err=None, y_err=None):
        # Define empty arrys/lists to store portions of data based on different spectral indices
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
        # Sort the data into different ranges of alpha
        for i in range(len(val_x)):
            if pp.config.getfloat('boundaries', 'low') > alpha[i] > pp.config.getfloat('boundaries', 'med'):
                val_flat_x.append(val_x[i])
                val_flat_x_err.append(x_err[i])
                val_flat_y.append(val_y[i])
                val_flat_y_err.append(y_err[i])
            elif pp.config.getfloat('boundaries', 'med') > alpha[i] > pp.config.getfloat('boundaries', 'high'):
                val_med_x.append(val_x[i])
                val_med_x_err.append(x_err[i])
                val_med_y.append(val_y[i])
                val_med_y_err.append(y_err[i])
            elif pp.config.getfloat('boundaries', 'high') > alpha[i]:
                val_steep_x.append(val_x[i])
                val_steep_x_err.append(x_err[i])
                val_steep_y.append(val_y[i])
                val_steep_y_err.append(y_err[i])
            else:
                val_outlier_x.append(val_x[i])
                val_outlier_x_err.append(x_err[i])
                val_outlier_y.append(val_y[i])
                val_outlier_y_err.append(y_err[i])

        # Create the plot:
        t = np.linspace(1e-5, 1e-1)  # datapoints in the plotting range
        plt.clf()  # clear pervious plots

        # defining labels
        l1 = pp.config.get('boundaries', 'high') + r'$\,>{\alpha}$'
        l2 = pp.config.get('boundaries', 'med') + r'$\,> {\alpha}>\,$' + pp.config.get('boundaries', 'high')
        l3 = pp.config.get('boundaries', 'low') + r'$\,> {\alpha}>\,$' + pp.config.get('boundaries', 'med')
        l6 = r'{Outliers}'
        l4 = '{Condon}'
        l5 = r'{Least Square Fit}'

        ax = plt.subplot(111)
        # double-logarithmic plots with advanced options (color, marker, ...)
        ax.errorbar(val_steep_x, val_steep_y, xerr=val_steep_x_err, yerr=val_steep_y_err, marker='.', linestyle='None',
                    color='b', label=l1)
        ax.errorbar(val_med_x, val_med_y, xerr=val_med_x_err, yerr=val_med_x_err, marker='.', linestyle='None',
                    color='g', label=l2)
        ax.errorbar(val_flat_x, val_flat_y, xerr=val_flat_x_err, yerr=val_flat_y_err, marker='.', linestyle='None',
                    color='r', label=l3)
        ax.plot(t, fct_f(t), linestyle='--', label=l4)
        ax.plot(t, fct_result(t, a, b), linestyle='-', label=l5)
        ax.errorbar(val_outlier_x, val_outlier_y, xerr=val_outlier_x_err, yerr=val_outlier_y_err, marker='.',
                    linestyle='None', color='tab:gray', label=l6)

        # Create a legend, labels ... be careful with latex symbols ...

        ax.legend()

        ax.set_xscale("log", nonposx='clip')
        ax.set_yscale("log", nonposy='clip')

        ax.grid(True)
        ax.set_ylabel(r'{$\left(\Sigma_{SFR}\right) _{RC}$ in M$_{Sun}$ yr$^{-1}$ kpc$^{-2}$}')
        ax.set_xlabel(r'{$\left(\Sigma_{SFR}\right) _{hyb}$ in M$_{Sun}$ yr$^{-1}$ kpc$^{-2}$}')
        # plt.xlim(1e-4,1e-1)
        # plt.ylim(1e-4,1e-1)
        # Save plots as pdf
        if (case == 'high' or case == 'low'):
            outfile = pp.config.get('names', case).rstrip('.fits') + '_pixel.pdf'
        else:
            outfile = pp.config.get('names', case) + '_pixel.pdf'
            tmp = pp.config.get('names', 'galaxy').split(' ')
            outfile = tmp[0] + '_' + tmp[1] + '_' + pp.config.get('names', case) + '_pixel_conv.pdf'

        if (case == 'conv_low'):
            title = pp.config.get('names',
                            'galaxy') + r' with Gaussian kernel, $l_{CRE} = $' + '%0.2f' % sigma + r' kpc , @' + pp.config.get(
                'values', 'freq_low') + r' MHz'
        elif (case == 'conv_high'):
            title = pp.config.get('names',
                            'galaxy') + r' with Gaussian kernel, $l_{CRE} = $' + '%0.2f' % sigma + r' kpc, @' + pp.config.get(
                'values', 'freq_high') + r' MHz'
        elif (case == 'low'):
            title = pp.config.get('names', 'galaxy') + r', @' + pp.config.get('values', 'freq_low') + r' MHz'
        elif (case == 'high'):
            title = pp.config.get('names', 'galaxy') + r', @' + pp.config.get('values', 'freq_high') + r' MHz'

        ax.set_title(title)
        plt.savefig(outfile)
        plt.clf()  # clean for further plotting

    def __print_results(self,pp, conv):
        ### Write the final results (fits and diffusion lengths) to file
        # this makes use of the configparser library again
        tmp = pp.config.get('names', 'galaxy').split(' ')
        results_ini = tmp[0] + '_' + tmp[1] + '_results.ini'
        res_out = configparser.ConfigParser()

        res_out['name'] = {'fname': pp.config.get('names', 'galaxy')}
                           # 'ds9cmd_low': cmd['low'],
                           # 'ds9cmd_high': cmd['high'],
                           # 'ds9cmd_sfr': cmd['sfr']}

        res_out['Low_freq_fit'] = {'#Fit results for pixel plots, with std errors\n'
                                   'a': str(pp.fit_low_a[0]),
                                   'a_err': str(pp.fit_low_a[1]),
                                   'b': str(pp.fit_low_b[0]),
                                   'b_err': str(pp.fit_low_b[1]),
                                   'chi_sqr': str(pp.fit_low_chi)}

        res_out['High_freq_fit'] = {'a': str(pp.fit_high_a[0]),
                                   'a_err': str(pp.fit_high_a[1]),
                                   'b': str(pp.fit_high_b[0]),
                                   'b_err': str(pp.fit_high_b[1]),
                                   'chi_sqr': str(pp.fit_high_chi)}

        res_out['Conv_results'] = {'#Diffusion length from gaussian kernel in kpc \n'
                                   'sigma_l': str(conv.optimal_sigma_low[0]),
                                   'a_err_l': str(conv.fit_low_conv_a[1]),
                                   'sigma_h': str(conv.optimal_sigma_high[0]),
                                   'a_err_h': str(conv.fit_high_conv_a[1])}

        res_out['rms_sigma'] = {'#The RMS values from calculated from the rms boxes\n'
                                'rms_l': str(pp.sigma['low']),
                                'rms_h': str(pp.sigma['high']),
                                'rms_s': str(pp.sigma['sfr'])}

        res_out['frequencies'] = {'#Frequencies for the two radio maps\n'
                                  'freq_low': pp.config.get('values', 'freq_low'),
                                  'freq_high': pp.config.get('values', 'freq_high')}

        with open(results_ini, 'w') as configfile:
            res_out.write(configfile)


