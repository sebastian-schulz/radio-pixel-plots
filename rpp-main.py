# System tools (file in and output ...)

import sys
import os

# Supressing runtime warning (mainly for maximum iterations of optimize)
import warnings

from pixelplots import Pixelplots
from convolution import Convolution
from plotting import Plotting

warnings.filterwarnings("ignore")

# rel. calibration error for both radio and sfr maps, we assume 5%
CALIB_ERR = 0.05

#######################
# PROGRAM STARTS HERE #
#######################


# Read in config file from command line argument, exit if it cant be found
if len(sys.argv) == 3:
    conf_file_name = str(sys.argv[1])
    if sys.argv[2] == '0':
        incl_corr_on = False
    elif sys.argv[2] == '1':
        incl_corr_on = True
    else:
        print('Unrecognized command line option ', sys.argv[2])
        sys.exit(-1)
else:
    print('Needs exactly two arguments: path to ini file and boolean for inclination correction. Quitting...')
    sys.exit(-1)

# Read in all the settings and data from the three fits files
create_pp = Pixelplots(conf_file_name, incl_corr_on)

create_pp.run()

create_conv = Convolution(create_pp)

create_conv.run()

create_plots = Plotting(create_pp, create_conv)

create_plots.png_all()

create_plots.print_all()

create_plots.plot_all()

tmp = create_pp.config.get('names','galaxy').split(' ')
tmp2 = tmp[0]+'_'+tmp[1]+'_plots_combined.pdf'

os.system('rm *plots_combined.pdf')
os.system('pdfunite *pdf '+ tmp2 )
print("Finished!")
