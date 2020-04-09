import configparser
import os
import math as m
# Library to understand the fits format
from astropy.io import fits
from calc_functions import convert_resolution_adv, calculate_rms, calc_rms_error, condon, fit_odr

# rel. calibration error for both radio and sfr maps, we assume 5%
CALIB_ERR = 0.05


class Pixelplots:
    def __init__(self, f_name, incl_corr_on):
        self.f_name = f_name
        self.incl_corr_on = incl_corr_on
        self.config = configparser.ConfigParser()
        self.config.read(self.f_name)
        print(self.config.sections())
        for key in self.config['names']:
            print(key, ' = ', self.config['names'][key])
        for key in self.config['values']:
            print(key, ' = ', self.config['values'][key])

        # change working directory to data-path (specified in the ini file)
        os.chdir(self.config['names']['fullpath'])
        print('Current working directory is:', os.getcwd())
        print('Inclination correction is set to:', incl_corr_on)

    def run(self):
        self.__load_fits()
        self.__convert_to_1200kpc()
        self.__flatten()
        self.__calc_errors()
        self.__convert_condon()
        self.__make_fit()

    def __load_fits(self):
        # Reading in the image file and checking size as well as shape
        self.data_high = fits.getdata(self.config.get('names', 'high'))
        self.data_low = fits.getdata(self.config.get('names', 'low'))
        self.data_sfr = fits.getdata(self.config.get('names', 'sfr'))
        # if it has the shape (1,1, X, Y), SLICE it to a 2D array (x,y)
        self.data_high = self.data_high.squeeze()
        self.data_low = self.data_low.squeeze()
        self.data_sfr = self.data_sfr.squeeze()

    def __convert_to_1200kpc(self):
        x_center = self.config.getint('values', 'center_x')
        y_center = self.config.getint('values', 'center_y')
        n_boxes = self.config.getint('values', 'n_boxes')
        # number of pixels along one direction of a box to achieve 1200pc length (float)
        self.px_per_box = 1000. * self.config.getfloat('values', 'kpc') / \
                          (self.config.getfloat('values', 'distance') * m.tan(m.radians(1. / 3600.)) /
                           self.config.getfloat('values', 'pixel_per_arcsec'))

        self.pixel_high_2d = convert_resolution_adv(self.data_high, x_center, y_center, n_boxes, self.px_per_box)
        self.pixel_low_2d = convert_resolution_adv(self.data_low, x_center, y_center, n_boxes, self.px_per_box)
        self.pixel_sfr_2d = convert_resolution_adv(self.data_sfr, x_center, y_center, n_boxes, self.px_per_box)

    def __flatten(self):
        # Flatten 2D arrays in to a long list; this is easier to plot
        self.pixel_high = [item for sublist in self.pixel_high_2d for item in sublist]
        self.pixel_low = [item for sublist in self.pixel_low_2d for item in sublist]
        self.pixel_sfr = [item for sublist in self.pixel_sfr_2d for item in sublist]

    def __calc_errors(self):
        # Apply the cuts to the data based on the signal
        # create three empty lists to store the 3 sigma cut data
        #  format is [ [values], [errors] ]
        self.pixel_low_cut = [[], [], [], []]
        self.pixel_high_cut = [[], [], [], []]
        self.pixel_sfr_cut = [[], []]
        self.alpha_cut = []

        self.pixel_low_fit = [[], [], [], []]
        self.pixel_high_fit = [[], [], [], []]
        self.pixel_sfr_fit = [[], []]
        self.alpha_fit = []

        # calculate rms from the boxes set in the config file
        sigma_high = calculate_rms(self.data_high,
                                   self.config.getint('high_cutoff_box', 'center_x'),
                                   self.config.getint('high_cutoff_box', 'center_y'),
                                   self.config.getint('high_cutoff_box', 'size_x'),
                                   self.config.getint('high_cutoff_box', 'size_y'))

        sigma_low = calculate_rms(self.data_low,
                                  self.config.getint('low_cutoff_box', 'center_x'),
                                  self.config.getint('low_cutoff_box', 'center_y'),
                                  self.config.getint('low_cutoff_box', 'size_x'),
                                  self.config.getint('low_cutoff_box', 'size_y'))

        sigma_sfr = calculate_rms(self.data_sfr,
                                  self.config.getint('sfr_cutoff_box', 'center_x'),
                                  self.config.getint('sfr_cutoff_box', 'center_y'),
                                  self.config.getint('sfr_cutoff_box', 'size_x'),
                                  self.config.getint('sfr_cutoff_box', 'size_y'))

        # Define alpha array to store the spectral index (assuming a simple power law)
        alpha_tmp = []
        for i in range(len(self.pixel_high)):
            alpha_tmp.append(m.log10(m.fabs(self.pixel_low[i] / self.pixel_high[i])) / m.log10(
                self.config.getfloat('values', 'freq_low') / self.config.getfloat('values', 'freq_high')))
        # 3 sigma cutoff for all datasets
        # for each dataset, there is a corresponding set containing the errors
        # calc_rms_err adds errors with calibration error (5%) of pixel value
        # to the sigma error: m.sqrt( sigma**2 + (value*0.05)**2 )
        for i in range(len(self.pixel_low)):
            if self.pixel_low[i] > 3. * sigma_low:
                if self.pixel_high[i] > 3. * sigma_high:
                    if self.pixel_sfr[i] > 3. * sigma_sfr:
                        self.alpha_cut.append(alpha_tmp[i])
                        self.pixel_low_cut[0].append(self.pixel_low[i])
                        self.pixel_low_cut[1].append(calc_rms_error(self.pixel_low[i], sigma_low, CALIB_ERR))
                        self.pixel_high_cut[0].append(self.pixel_high[i])
                        self.pixel_high_cut[1].append(calc_rms_error(self.pixel_high[i], sigma_high, CALIB_ERR))
                        self.pixel_sfr_cut[0].append(self.pixel_sfr[i])
                        self.pixel_sfr_cut[1].append(calc_rms_error(self.pixel_sfr[i], sigma_sfr, CALIB_ERR))
                        if alpha_tmp[i] < self.config.getfloat('boundaries', 'low'):
                            self.alpha_fit.append(alpha_tmp[i])
                            self.pixel_low_fit[0].append(self.pixel_low[i])
                            self.pixel_low_fit[1].append(calc_rms_error(self.pixel_low[i], sigma_low, CALIB_ERR))
                            self.pixel_high_fit[0].append(self.pixel_high[i])
                            self.pixel_high_fit[1].append(calc_rms_error(self.pixel_high[i], sigma_high, CALIB_ERR))
                            self.pixel_sfr_fit[0].append(self.pixel_sfr[i])
                            self.pixel_sfr_fit[1].append(calc_rms_error(self.pixel_sfr[i], sigma_sfr, CALIB_ERR))

        self.sigma = {'low': sigma_low, 'high': sigma_high, 'sfr': sigma_sfr}
        print('RMS for the 3 different maps, used as sigma for 3 sigma cutoff:\n', self.sigma)

        # Some diagnostic output of how many points were cut
        print('Number of points in sample:', len(self.pixel_low))
        print('Number of points after 3 sigma cutoff:', len(self.pixel_low_cut[0]))

        mean_old = 0
        tmp = [item for sublist in self.data_low for item in sublist]
        for i in range(len(tmp)):
            mean_old += m.fabs(tmp[i])
        mean = 0
        for i in range(len(self.pixel_low)):
            mean += self.pixel_low[i]
        mean *= self.px_per_box ** 2
        print('Total total sum of lower freq. image:\t', '%0.3f' % mean_old)
        print('Total sum in the cut lower freq. map is:', '%0.3f' % mean)

    def __convert_condon(self):
        # Now convert lofar and wrst to SFR sufrace density using condon-relation
        self.pixel_low_cut[2] = condon(self.pixel_low_cut[0],
                                       self.config.getfloat('values', 'FWHM'),
                                       self.config.getfloat('values', 'freq_low'))

        self.pixel_high_cut[2] = condon(self.pixel_high_cut[0],
                                        self.config.getfloat('values', 'FWHM'),
                                        self.config.getfloat('values', 'freq_high'))

        self.pixel_low_fit[2] = condon(self.pixel_low_fit[0],
                                       self.config.getfloat('values', 'FWHM'),
                                       self.config.getfloat('values', 'freq_low'))

        self.pixel_high_fit[2] = condon(self.pixel_high_fit[0],
                                        self.config.getfloat('values', 'FWHM'),
                                        self.config.getfloat('values', 'freq_high'))

        # The same conversion is performed on the errors
        self.pixel_low_cut[3] = condon(self.pixel_low_cut[1],
                                       self.config.getfloat('values', 'FWHM'),
                                       self.config.getfloat('values', 'freq_low'))

        self.pixel_high_cut[3] = condon(self.pixel_high_cut[1],
                                        self.config.getfloat('values', 'FWHM'),
                                        self.config.getfloat('values', 'freq_high'))

        self.pixel_low_fit[3] = condon(self.pixel_low_fit[1],
                                       self.config.getfloat('values', 'FWHM'),
                                       self.config.getfloat('values', 'freq_low'))

        self.pixel_high_fit[3] = condon(self.pixel_high_fit[1],
                                        self.config.getfloat('values', 'FWHM'),
                                        self.config.getfloat('values', 'freq_high'))

    def __make_fit(self):
        self.fit_low_a, self.fit_low_b, self.fit_low_chi = fit_odr(self.pixel_sfr_fit[0],
                                                                   self.pixel_low_fit[2],
                                                                   self.pixel_sfr_fit[1],
                                                                   self.pixel_low_fit[3])

        self.fit_high_a, self.fit_high_b, self.fit_high_chi = fit_odr(self.pixel_sfr_fit[0],
                                                                      self.pixel_high_fit[2],
                                                                      self.pixel_sfr_fit[1],
                                                                      self.pixel_high_fit[3])
