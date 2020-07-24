import sys
import os
from vollmer import Vollmer
from pixelplots import Pixelplots


# Read in config file from command line argument, exit if it cant be found
if len(sys.argv) == 3:
    conf_file_name = str(sys.argv[1])
    if sys.argv[2] == 'l':
        case = 'high'
    elif sys.argv[2] == 'h':
        case = 'low'
    else:
        print('Unrecognized command line option ', sys.argv[2])
        sys.exit(-1)
else:
    print('Needs exactly two arguments: path to ini file and boolean for inclination correction. Quitting...')
    sys.exit(-1)
# Read in all the settings and data from the three fits files
create_pp = Pixelplots(conf_file_name, False)
create_pp.run()

if not os.path.isdir('./vollmer'):
    os.mkdir('./vollmer')
os.chdir('./vollmer')

compare_vollmer = Vollmer(create_pp, case)
compare_vollmer.run(kernel_type='round_exp', n=0)  # round_exp, round_gauss, elliptical not yet implemented
compare_vollmer.run(kernel_type='round_exp', n=0.225)
compare_vollmer.run(kernel_type='round_exp', n=0.5)
compare_vollmer.run(kernel_type='round_gauss', n=0)
compare_vollmer.run(kernel_type='round_gauss', n=0.225)
compare_vollmer.run(kernel_type='round_gauss', n=0.5)
