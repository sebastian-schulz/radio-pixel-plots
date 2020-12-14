# Pythons math package with all the constants and functions
import math as m
# Astropy library that builds 2d gauss kernels that are rotated by some agnle theta wrt the coordinate axis
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve_fft as ap_convole
from astropy.io import fits
# Library includes a zero-finder method (to find the best value for the kernel), also includes methods fit data/curves
from scipy import optimize
from calc_functions import fit_odr, calc_rms_error, convert_kpc2px, condon, convert_resolution_adv

# convert_px2kpc

# rel. calibration error for both radio and sfr maps, we assume 5%
CALIB_ERR = 0.05


class Convolution:
    def __init__(self, pp):
        self.pp = pp
        self.incl_corr_on = pp.incl_corr_on
        if self.pp.incl_corr_on:
            self.phi = m.radians(pp.config.getfloat('values', 'phi'))
            self.incl = m.radians(pp.config.getfloat('values', 'incl'))
        self.sigma_conv = pp.config.getfloat('values', 'sigma_conv')
        self.px_per_as = pp.config.getfloat('values', 'pixel_per_arcsec')
        self.distance = pp.config.getfloat('values', 'distance')

    def run(self):
        self.__fit_sigma()  # This does the actual fitting, done with optimize.fsolve
        # Create full resolution convolved sfr maps for high/low and safe them as fits
        self.data_sfr_conv_low = self.__convolve_gauss(self.sigma_low)
        self.data_sfr_conv_high = self.__convolve_gauss(self.sigma_high)
        hdul = fits.open(self.pp.config.get('names', 'sfr'))
        hdul[0].data = self.__convolve_gauss(self.sigma_low)
        tmp_str = self.pp.config.get('names', 'sfr').rstrip('.fits') + '_conv_high.fits'
        hdul.writeto(tmp_str, overwrite=True)
        hdul[0].data = self.__convolve_gauss(self.sigma_high)
        tmp_str = self.pp.config.get('names', 'sfr').rstrip('.fits') + '_conv_low.fits'
        hdul.writeto(tmp_str, overwrite=True)


    def __fit_sigma(self):
        print('Finding optimal gaussian kernel. This may take a moment...')
        self.optimal_l_low = optimize.fsolve(self.__fct_gauss_fit,  # fitting function
                                             self.sigma_conv,  # starting value
                                             args=(self.pp, False),  # additional arguments
                                             maxfev=15)  # max. no. of iterations
        self.sigma_low = self.optimal_l_low[0]
        # Diffusion length is the FWHM/2 of the final image
        self.optimal_l_low[0] = m.sqrt(m.pow(2.3548 * self.optimal_l_low[0], 2)) / 2.  # - m.pow(1.2, 2)
        print('Final value for Diffusion length:\t', '%0.3f' % self.optimal_l_low[0], 'kpc')

        self.optimal_l_high = optimize.fsolve(self.__fct_gauss_fit,
                                              self.sigma_conv,
                                              args=(self.pp, True),
                                              maxfev=15)
        self.sigma_high = self.optimal_l_high[0]
        self.optimal_l_high[0] = m.sqrt(m.pow(2.3548 * self.optimal_l_high[0], 2)) / 2.  # - m.pow(1.2, 2)
        print('Final value for Diffusion length:\t', '%0.3f' % self.optimal_l_high[0], 'kpc')

    def __fct_gauss_fit(self, sigma, pp, is_high):
        # Gaussian Kernel function for optimization purposes, because optimize searches for zeros by default
        x = self.__fct_gauss(sigma, pp, is_high)
        return x - 1.

    def __fct_gauss(self, sigma, pp, is_high):
        # Function that calculates the fit parameter a (slope) as a function of gaussian kernel width sigma
        # Convolve with gaussian:
        self.data_sfr_conv = self.__convolve_gauss(sigma)
        # Compute new spectral indices (may have changed due to different cutting)
        alpha_conv = []
        for i in range(len(pp.pixel_low)):
            alpha_conv.append(m.log10(m.fabs(pp.pixel_low[i] / pp.pixel_high[i])) /
                              m.log10(pp.config.getfloat('values', 'freq_low') /
                                      pp.config.getfloat('values', 'freq_high')))

        pixel_sfr_2d_conv = convert_resolution_adv(self.data_sfr_conv,
                                                   pp.config.getfloat('values', 'center_x'),
                                                   pp.config.getfloat('values', 'center_y'),
                                                   pp.config.getfloat('values', 'n_boxes'),
                                                   pp.px_per_box)

        pixel_sfr_conv = [item for sublist in pixel_sfr_2d_conv for item in sublist]

        pixel_sfr_conv_cut = [[], []]
        pixel_high_conv_cut = [[], [], [], []]
        pixel_low_conv_cut = [[], [], [], []]
        alpha_conv_cut = []
        # applying the same 3 sigma cut as for the non-convolved data
        for i in range(len(pixel_sfr_conv)):
            if pp.pixel_low[i] > 3. * pp.sigma['low']:
                if pp.pixel_high[i] > 3. * pp.sigma['high']:
                    if pixel_sfr_conv[i] > 3. * pp.sigma['sfr']:
                        if alpha_conv[i] <= pp.config.getfloat('boundaries', 'low'):
                            alpha_conv_cut.append(alpha_conv[i])
                            pixel_sfr_conv_cut[0].append(pixel_sfr_conv[i])
                            pixel_sfr_conv_cut[1].append(calc_rms_error(pixel_sfr_conv[i], pp.sigma['sfr'], CALIB_ERR))
                            pixel_low_conv_cut[0].append(pp.pixel_low[i])
                            pixel_low_conv_cut[1].append(calc_rms_error(pp.pixel_low[i], pp.sigma['low'], CALIB_ERR))
                            pixel_high_conv_cut[0].append(pp.pixel_high[i])
                            pixel_high_conv_cut[1].append(calc_rms_error(pp.pixel_high[i], pp.sigma['high'], CALIB_ERR))
        print('Number of points in sample after cutting:', len(alpha_conv_cut))
        self.no_cut_points = len(alpha_conv_cut)
        # Now apply the condon relation to the radio map set via parameter 'opt' and return values
        pixel_high_conv_cut[2] = condon(pixel_high_conv_cut[0], pp.config.getfloat('values', 'FWHM'),
                                             pp.config.getfloat('values', 'freq_high'))
        pixel_high_conv_cut[3] = condon(pixel_high_conv_cut[1], pp.config.getfloat('values', 'FWHM'),
                                             pp.config.getfloat('values', 'freq_high'))
        pixel_low_conv_cut[2] = condon(pixel_low_conv_cut[0], pp.config.getfloat('values', 'FWHM'),
                                            pp.config.getfloat('values', 'freq_low'))
        pixel_low_conv_cut[3] = condon(pixel_low_conv_cut[1], pp.config.getfloat('values', 'FWHM'),
                                            pp.config.getfloat('values', 'freq_low'))

        if is_high:
            self.fit_high_conv_a, self.fit_high_conv_b, _, _ = fit_odr(pixel_sfr_conv_cut[0], pixel_high_conv_cut[2], val_x_err=pixel_sfr_conv_cut[1],
                                        val_y_err=pixel_high_conv_cut[3])
            self.pixel_sfr_conv_cut_high = pixel_sfr_conv_cut
            self.pixel_high_conv_cut_high = pixel_high_conv_cut
            self.pixel_low_conv_cut_high = pixel_low_conv_cut
            self.alpha_conv_cut_high = alpha_conv_cut
            return self.fit_high_conv_a[0]

        elif not is_high:
            self.fit_low_conv_a, self.fit_low_conv_b, _, _ = fit_odr(pixel_sfr_conv_cut[0], pixel_low_conv_cut[2], val_x_err=pixel_sfr_conv_cut[1],
                                      val_y_err=pixel_low_conv_cut[3])
            self.pixel_sfr_conv_cut_low = pixel_sfr_conv_cut
            self.pixel_high_conv_cut_low = pixel_high_conv_cut
            self.pixel_low_conv_cut_low = pixel_low_conv_cut
            self.alpha_conv_cut_low = alpha_conv_cut
            return self.fit_low_conv_a[0]

    def __convolve_gauss(self, sigma_in):
        if self.pp.incl_corr_on:
            incl = self.incl
            phi = self.phi
        elif not self.pp.incl_corr_on:
            incl = 0.
            phi = 0.
        else:
            print('Inclination correction not properly set in function convolve_gauss!')
            quit()
        sigma_y = convert_kpc2px(sigma_in, self.distance, self.px_per_as)
        sigma_x = sigma_y * m.cos(incl)
        kernel2d = Gaussian2DKernel(sigma_x, y_stddev=sigma_y, theta=phi)
        res = ap_convole(self.pp.data_sfr, kernel2d, )
        return res
