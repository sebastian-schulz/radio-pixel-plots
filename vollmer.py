from pixelplots import Pixelplots
import math as m
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['T1']})# T1 is the LaTex standard font, also used by scr familiy documents as default
rc('text', usetex=True)

from adaptive_convolution import AdaptiveConvolution

class Vollmer:
    def __init__(self, pp, case='low'):
        self.pp = pp
        if case == 'high':
            self.radio_map = pp.pixel_high_2d
            self.radio_pixels = pp.pixel_high_cut
        else:
            self.radio_map = pp.pixel_low_2d
            self.radio_pixels = pp.pixel_low_cut

        self.sfr_map = pp.pixel_sfr_2d

    def run(self, kernel_type, n=0):
        phi = []
        l = []
        for i in np.arange(0.3, 3.5, 0.1):
            adaptive_conv = AdaptiveConvolution(self.sfr_map, k=21, exp=n, l_0=i, sigma_0=8e-3, method=kernel_type)
            adaptive_conv.convolve()
            sfr_conv = adaptive_conv.conv_map
            l.append(i)
            phi.append(self.__calc_phi(self.radio_map, sfr_conv))
        #for i in range(len(l)):
        #    print(l[i], '\t', m.log10(phi[i][0]))
        log_phi = []
        for i in range(len(phi)):
            log_phi.append(m.log10(phi[i]))
        print('Optimal kernel size is:', l[np.argmin(log_phi)])
        fig, ax = plt.subplots()  # Create a figure containing a single axes.
        ax.plot(l, log_phi)  # Plot some data on the axes.
        plt.savefig(kernel_type+'_'+str(n)+'.png')
        plt.clf()  # clean for further plotting

    def __calc_phi(self, radio_map, sfr_map_conv):
        Q = 0
        tmp = 0
        phi = 0
        for i in range(len(sfr_map_conv)):
            for j in range(len(sfr_map_conv[0])):
                Q += m.pow(sfr_map_conv[i][j],2)
                tmp += radio_map[i][j] * sfr_map_conv[i][j]
        Q = Q / tmp
        tmp = 0
        for i in range(len(radio_map)):
            for j in range(len(radio_map[0])):
                phi += m.pow(radio_map[i][j] - sfr_map_conv[i][j] / Q, 2)
                tmp += m.pow(radio_map[i][j], 2)
        phi = phi/tmp
        return phi


