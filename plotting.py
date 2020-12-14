# Library for python plotting
import configparser
import matplotlib.pyplot as plt
# Mighty numerical library of ptyhon
import numpy as np
# For LateX fonts and symbols in plots
from matplotlib import rc, rcParams
import os
from calc_functions import fct_f, fct_result


hi=3
wi=4.5
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['T1'],'size': 8})# T1 is the LaTex standard font, also used by scr familiy documents as default
rc('text', usetex=True)
rcParams.update({'figure.autolayout': True})

class Plotting:
    def __init__(self, pp, conv):
        self.pp = pp
        self.conv = conv

    def png_all(self):
        self.map_plt(self.pp.data_high, 'high')
        self.map_plt(self.pp.data_low, 'low')
        self.map_plt(self.pp.data_sfr, 'sfr')
        self.__ds9_map(self.pp)
        self.col_map_plt(self.pp.data_sfr, 'sfr_colormap.png')
        self.col_map_plt(self.conv.data_sfr_conv_low, 'low_convolved_colormap.png')
        self.col_map_plt(self.conv.data_sfr_conv_high, 'high_convolved_colormap.png')
        self.__sfr_1200pc_plt(self.pp.pixel_sfr_2d, 'sfr')

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
                     self.pp.fit_low_a,
                     self.pp.fit_low_b,
                    'low',
                     x_err=self.pp.pixel_sfr_cut[1],
                     y_err=self.pp.pixel_low_cut[3])
        self.__pplot(self.pp,
                     self.pp.pixel_sfr_cut[0],
                     self.pp.pixel_high_cut[2],
                     self.pp.alpha_cut,
                     self.pp.fit_high_a,
                     self.pp.fit_high_b,
                    'high',
                     x_err=self.pp.pixel_sfr_cut[1],
                     y_err=self.pp.pixel_high_cut[3])

        self.__pplot(self.pp,
                     self.conv.pixel_sfr_conv_cut_low[0],
                     self.conv.pixel_low_conv_cut_low[2],
                     self.conv.alpha_conv_cut_low,
                     self.conv.fit_low_conv_a,
                     self.conv.fit_low_conv_b,
                    'conv_low',
                     sigma= self.conv.optimal_l_low[0],
                     x_err=self.conv.pixel_sfr_conv_cut_low[1],
                     y_err=self.conv.pixel_low_conv_cut_low[3])
        self.__pplot(self.pp,
                     self.conv.pixel_sfr_conv_cut_high[0],
                     self.conv.pixel_high_conv_cut_high[2],
                     self.conv.alpha_conv_cut_high,
                     self.conv.fit_high_conv_a,
                     self.conv.fit_high_conv_b,
                    'conv_high',
                     sigma=self.conv.optimal_l_high[0],
                     x_err=self.conv.pixel_sfr_conv_cut_high[1],
                     y_err=self.conv.pixel_high_conv_cut_high[3])

    # TEST
    def map_plt(self, data, case):
        values, edges = np.histogram(data, bins=1000000, range=(0, data.max()))
        cmax = 0
        total = 0
        i = 0
        for l in range(len(values)):
            total += values[l]
        while True:
            cmax +=values[i]
            if cmax > total * .99:
                break
            i += 1
        #tmp = np.zeros(shape=(data.shape))
        #for i in range(len(data)):
        #    for j in range(len(data)):
        #        tmp[i, j] = m.sqrt(m.fabs(data[i, j]))
            # Plot the modified data as a greyscale map (compare with original)
        plt.imshow(data, cmap='gray', clim=(0, edges[i]))  # clim=(edges[9]) )
        plt.colorbar()
        plt.ylim(0, len(data) - 1)
#        image_name = cfg.get('names', galaxyname)
#        image_name = image_name.rstrip('.fits')
        tmp = self.pp.config.get('names', 'galaxy').split(' ')
        image_name = 'n' + tmp[1] + '_' + case + '_099_map.png'
        # And save the plot to file:
        plt.savefig(image_name, dpi=500)
        plt.clf()  # clears the plot for further plotting


    def col_res_plt(self, data, fname):
        fig, ax = plt.subplots(1, 1)
        pcm = ax.pcolormesh(data, cmap='RdBu_r', vmax=np.max(data), vmin=0)
        fig.colorbar(pcm, ax=ax, extend='both')
        plt.savefig(fname, dpi=500)
        plt.clf()

    def col_map_plt(self, data, fname):
        fig, ax = plt.subplots(1, 1)
        pcm = ax.pcolormesh(data, cmap='RdBu_r', vmax=np.max(data), vmin=0)
        fig.colorbar(pcm, ax=ax, extend='both')
        plt.savefig(fname, dpi=500)
        plt.clf()

    def __sfr_1200pc_plt(self, data, case):
           plt.imshow(data, cmap='gray')  # clim=(edges[9]) )
           plt.colorbar()
           plt.ylim(0, len(data) - 1)
           tmp = self.pp.config.get('names', 'galaxy').split(' ')
           image_name = 'n' + tmp[1] + '_' + case + '_1200pc_map.png'
           # And save the plot to file:
           plt.savefig(image_name, dpi=500)
           plt.clf()  # clears the plot for further plotting

    # Prints data arrays to file as tab separated values
    def __print_data(self, pp):
        # mean = 0.
        tmp = pp.config.get('names', 'galaxy').split(' ')
        dataname = tmp[0] + '_' + tmp[1] + '_pixel.dat'
        f_tmp = open(dataname, 'w')
        #print(dataname)
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
        dataname = tmp[0] + '_' + tmp[1] + '_' +'high' + '_pixel_conv.dat'
        f_tmp = open(dataname, 'w')
        #print(dataname)
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
        dataname = tmp[0] + '_' + tmp[1] + '_' + 'low' + '_pixel_conv.dat'
        f_tmp = open(dataname, 'w')
        #print(dataname)
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
        t = np.linspace(1e-10, 1e-0)  # datapoints in the plotting range
        plt.clf()  # clear pervious plots

        # defining labels
        l1 = pp.config.get('boundaries', 'high') + r'$\,>{\alpha}$'
        l2 = pp.config.get('boundaries', 'med') + r'$\,> {\alpha}>\,$' + pp.config.get('boundaries', 'high')
        l3 = pp.config.get('boundaries', 'low') + r'$\,> {\alpha}>\,$' + pp.config.get('boundaries', 'med')
        l6 = r'{Outliers}'
        l4 = r'{Condon--relation}'
        l5 = r'ODR fit $a=\,$'+ f'{a[0]:.2f}'+ r'$\,\pm\,$' + f'{a[1]:.2f}'

        point_size = 3
        error_line_width = 0.45
        plt.gcf().set_size_inches(wi,hi)
        ax = plt.subplot(111)
        # double-logarithmic plots with advanced options (color, marker, ...)
        ax.errorbar(val_steep_x, val_steep_y, xerr=val_steep_x_err, yerr=val_steep_y_err, marker='^', ms=point_size, elinewidth=error_line_width, linestyle='None',
                    color='b', label=l1)
        ax.errorbar(val_med_x, val_med_y, xerr=val_med_x_err, yerr=val_med_x_err, marker='o', ms=point_size, elinewidth=error_line_width, linestyle='None',
                    color='g', label=l2)
        ax.errorbar(val_flat_x, val_flat_y, xerr=val_flat_x_err, yerr=val_flat_y_err, marker='v', ms=point_size, elinewidth=error_line_width, linestyle='None',
                    color='r', label=l3)
        ax.plot(t, fct_f(t), linestyle='--', label=l4)
        ax.plot(t, fct_result(t, a[0], b[0]), linestyle='-', label=l5)
        ax.errorbar(val_outlier_x, val_outlier_y, xerr=val_outlier_x_err, yerr=val_outlier_y_err, marker='.',
                    linestyle='None', color='tab:gray', label=l6)

        # Create a legend, labels ... be careful with latex symbols ...

        ax.legend()

        ax.set_xscale("log", nonposx='clip')
        ax.set_yscale("log", nonposy='clip')

        ax.grid(True)
        ax.set_ylabel(r'{$\left(\Sigma_{\mbox{SFR}}\right) _{\mbox{RC}}$ in M$_{\odot}$ yr$^{-1}$ kpc$^{-2}$}')
        ax.set_xlabel(r'{$\left(\Sigma_{\mbox{SFR}}\right) _{\mbox{hyb}}$ in M$_{\odot}$ yr$^{-1}$ kpc$^{-2}$}')
        scale=2.5
        ax.set_xlim(min(val_x)/scale, max(val_x)*scale)
        ax.set_ylim(min(val_y)/scale, max(val_y)*scale)
        # Save plots as pdf
        if (case == 'high' or case == 'low'):
            outfile = pp.config.get('names', case).rstrip('.fits') + '_pixel.pdf'
        else:
            outfile = pp.config.get('names', case) + '_pixel.pdf'
            tmp = pp.config.get('names', 'galaxy').split(' ')
            outfile = tmp[0] + '_' + tmp[1] + '_' + pp.config.get('names', case) + '_pixel_conv.pdf'

        if (case == 'conv_low'):
            title = pp.config.get('names',
                            'galaxy') +r' @' +\
                            pp.config.get('values', 'freq_low') +\
                            r' MHz' + r' convolved with a Gaussian ($l_{\mbox{CRE}} = $' + '%0.2f' % sigma + r' kpc)'
        elif (case == 'conv_high'):
            title = pp.config.get('names',
                            'galaxy') + r' @' + \
                            pp.config.get('values', 'freq_high') + \
                            r' MHz' + r' convolved with a Gaussian ($l_{\mbox{CRE}} = $' + '%0.2f' % sigma + r' kpc)'
        elif (case == 'low'):
            title = pp.config.get('names', 'galaxy') + r' @' + pp.config.get('values', 'freq_low') + r' MHz'
        elif (case == 'high'):
            title = pp.config.get('names', 'galaxy') + r' @' + pp.config.get('values', 'freq_high') + r' MHz'

        ax.set_title(title)
        plt.savefig(outfile)
        plt.clf()  # clean for further plotting

    def __print_results(self,pp, conv):
        ### Write the final results (fits and diffusion lengths) to file
        # this makes use of the configparser library again
        tmp = pp.config.get('names', 'galaxy').split(' ')
        results_ini = tmp[0] + '_' + tmp[1] + '_results.ini'
        res_out = configparser.ConfigParser()
        res_out.read(results_ini)
        if not res_out.has_section('name'):
            res_out.add_section('name')
        res_out['name'] = {'fname': pp.config.get('names', 'galaxy')}
#                            'ds9cmd_low': cmd['low'],
#                           'ds9cmd_high': cmd['high'],
#                            'ds9cmd_sfr': cmd['sfr']}
        if not res_out.has_section('Low_freq_fit'):
            res_out.add_section('Low_freq_fit')
        res_out['Low_freq_fit'] = {'#Fit results for pixel plots, with std errors\n'
                                   'a': str(pp.fit_low_a[0]),
                                   'a_err': str(pp.fit_low_a[1]),
                                   'b': str(pp.fit_low_b[0]),
                                   'b_err': str(pp.fit_low_b[1]),
                                   'chi_sqr': str(pp.fit_low_chi)}
        if not res_out.has_section('High_freq_fit'):
            res_out.add_section('High_freq_fit')
        res_out['High_freq_fit'] = {'a': str(pp.fit_high_a[0]),
                                   'a_err': str(pp.fit_high_a[1]),
                                   'b': str(pp.fit_high_b[0]),
                                   'b_err': str(pp.fit_high_b[1]),
                                   'chi_sqr': str(pp.fit_high_chi)}
        if not res_out.has_section('Conv_results'):
            res_out.add_section('Conv_results')
        res_out['Conv_results'] = {'#Diffusion length from gaussian kernel in kpc \n'
                                   'sigma_l': str(conv.optimal_l_low[0]),
                                   'a_err_l': str(conv.fit_low_conv_a[1]),
                                   'sigma_h': str(conv.optimal_l_high[0]),
                                   'a_err_h': str(conv.fit_high_conv_a[1])}
        if not res_out.has_section('rms_sigma'):
            res_out.add_section('rms_sigma')
        res_out['rms_sigma'] = {'#The RMS values from calculated from the rms boxes\n'
                                'rms_l': str(pp.sigma['low']),
                                'rms_h': str(pp.sigma['high']),
                                'rms_s': str(pp.sigma['sfr'])}
        if not res_out.has_section('frequencies'):
            res_out.add_section('frequencies')
        res_out['frequencies'] = {'#Frequencies for the two radio maps\n'
                                  'freq_low': pp.config.get('values', 'freq_low'),
                                  'freq_high': pp.config.get('values', 'freq_high')}
        if not res_out.has_section('number_of_points'):
            res_out.add_section('number_of_points')
        res_out['number_of_points'] = {'#Error clipping comparison\n'
                                  'pixelplot': pp.no_points_cut,
                                  'convolution': conv.no_cut_points}

        with open(results_ini, 'w+') as configfile:
            res_out.write(configfile)

    def __ds9_map(self, pp):
        for fname in (['low', 'high', 'sfr']):
            oname = pp.config.get('names', fname).rstrip('.fits') + '_ds9box.png'
            cmd = ('ds9 ' + pp.config.get('names', fname) +
                          ' -regions system image ' +
                          ' -regions command "box ' +
                          pp.config.get(fname + '_cutoff_box', 'center_x') + ' ' +
                          pp.config.get(fname + '_cutoff_box', 'center_y') + ' ' +
                          pp.config.get(fname + '_cutoff_box', 'size_x') + ' ' +
                          pp.config.get(fname + '_cutoff_box', 'size_y') +
                          ' # color=red" ' +
                          '-regions command "box ' +
                          pp.config.get('values', 'center_x') + ' ' +
                          pp.config.get('values', 'center_y') + ' ' +
                          str(pp.px_per_box * pp.config.getint('values', 'n_boxes')) + ' ' +
                          str(pp.px_per_box * pp.config.getint('values', 'n_boxes')) +
                          ' # color=yellow" ' +
                          ' -zoom to fit -scale mode 95  -saveimage png ' +
                          oname + ' -quit')
            os.system(cmd)
            #print(cmd.rstrip(' -quit'))
        # Save the DS9 command as string to print to results file later
#        for key in cmd:
#           cmd[key] = cmd[key].rstrip('-quit')
