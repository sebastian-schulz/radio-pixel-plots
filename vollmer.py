import math as m
import numpy as np
import matplotlib.pyplot as plt
from adaptive_convolution import AdaptiveConvolution
import configparser

from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# for Palatino and other serif fonts use:
rc('font', **{'family': 'serif', 'serif': ['T1']})
# T1 is the LaTex standard font, also used by scr familiy documents as default
rc('text', usetex=True)


class Vollmer:
    def __init__(self, pp, case):
        self.pp = pp
        if case == 'high':
            self.radio_map = self.pp.pixel_high_2d
            self.radio_pixels = self.pp.pixel_high_cut  # Not used yet
        else:  #DEFAULT
            self.radio_map = self.pp.pixel_low_2d
            self.radio_pixels = self.pp.pixel_low_cut  # Not used yet
        self.sfr_map = self.pp.pixel_sfr_2d
        self.case = case
        print('Using' + self.case + 'radio data as comparison.')

    def run(self, kernel_type):
        self.l1, self.log_phi1 = self.convolve(kernel_type, n=0)
        self.l2, self.log_phi2 = self.convolve(kernel_type, n=0.225)
        self.l3, self.log_phi3 = self.convolve(kernel_type, n=0.5)
        self.__plot(kernel_type)

    def convolve(self, kernel_type, n, k=13):
        print('Running smoothing experiment with ' + kernel_type + ' kernel and n = ', str(n))
        phi = []
        l = []
        for i in np.arange(0.1, 8.1, 0.1):
            adaptive_conv = AdaptiveConvolution(self.sfr_map, k=k, exp=n, l_0=i, sigma_0=8e-3, method=kernel_type)
            adaptive_conv.convolve()
            sfr_conv = adaptive_conv.conv_map
            l.append(i)
            print('l =\t', str(i))
            phi.append(self.__calc_phi(self.radio_map, sfr_conv))
        log_phi = []
        for i in range(len(phi)):
            log_phi.append(m.log10(phi[i]))
        print('Optimal kernel size is:', l[np.argmin(log_phi)])
        tmp = self.pp.config.get('names', 'galaxy').split(' ')
        results_ini = '../' + tmp[0] + '_' + tmp[1] + '_results.ini'
        res_out = configparser.ConfigParser()
        res_out.read(results_ini)
        sect = 'vollmer' + '_' + kernel_type
        if not res_out.has_section(sect):
            res_out.add_section(sect)
        res_out.set(sect, self.case + '_' + str(n), str(l[np.argmin(log_phi)]))
        with open(results_ini, 'w+') as configfile:
            res_out.write(configfile)
        return l, log_phi

    def __plot(self, kernel_type):
        fig, ax = plt.subplots()  # Create a figure containing a single set of axes.
        if kernel_type == 'round_gauss':
            l1 = r'Gauss kernel with n = '
        elif kernel_type == 'round_exp':
            l1 = r'Exponential kernel with n = '
        l1 = l1 + '0'
        l2 = l1 + '.225'
        l3 = l1 + '.5'
        # Plotting the resulting data
        ax.plot(self.l1, self.log_phi1, label=l1)  # Plot some data on the axes.
        ax.plot(self.l2, self.log_phi2, label=l2)  # Plot some data on the axes.
        ax.plot(self.l3, self.log_phi3, label=l3)  # Plot some data on the axes.
        ax.grid(True)
        ax.set_xlabel(r'Smoothing length $l$ in kpc')
        ax.set_ylabel(r'Goodness of fit parameter $\log(\phi)$')
        ax.legend()
        tmp_str = self.pp.config.get('values', 'freq_' + self.case)
        ax.set_title('Smoothing for ' + self.pp.config.get('names', 'galaxy') + ' w.r.t ' + tmp_str + ' MHz radio data')
        plt.savefig(kernel_type + '_' + self.case + '.png')  # .pdf
        plt.clf()  # clean for further plotting
        # Saving the data to file
        name = kernel_type + '_' + self.case + '.dat'
        file_tmp = open(name, 'w')
        file_tmp.write('# smoothing length in kpc\t goodness of fit parameter log(phi) n=0\t n=0.25\t n=0.5\n')
        for i in range(len(self.l1)):
            file_tmp.write(str(self.l1[i]) + '\t' + str(self.log_phi1[i]) + '\t' + str(self.log_phi2[i]) + '\t' + str(self.log_phi3[i]) + '\n')


    def __calc_phi(self, radio_map, sfr_map_conv):
        q = 0
        tmp = 0
        phi = 0
        for i in range(len(sfr_map_conv)):
            for j in range(len(sfr_map_conv[0])):
                q += m.pow(sfr_map_conv[i][j], 2)
                tmp += radio_map[i][j] * sfr_map_conv[i][j]
        q = q / tmp
        tmp = 0
        for i in range(len(radio_map)):
            for j in range(len(radio_map[0])):
                phi += m.pow(radio_map[i][j] - sfr_map_conv[i][j] / q, 2)
                tmp += m.pow(radio_map[i][j], 2)
        phi = phi/tmp
        return phi
