import sys
import os
from vollmer import Vollmer
from pixelplots import Pixelplots
import time
start_time = time.time()



# Read in config file from command line argument, exit if it cant be found
if len(sys.argv) == 2:
    conf_file_name = str(sys.argv[1])
#    if sys.argv[2] == 'h':
#        case = 'high'
#    elif sys.argv[2] == 'l':
#        case = 'low'
else:
    print('Needs exactly one argument: path to ini file. Quitting...')
    sys.exit(-1)
# Read in all the settings and data from the three fits files
create_pp = Pixelplots(conf_file_name, False)
create_pp.run()

if not os.path.isdir('./vollmer'):
    os.mkdir('./vollmer')
os.chdir('./vollmer')

#for i in [3,5,7,9,11,13,15,17,19,21]:
# default kernel size is 13, set in vollmer.py
compare_low = Vollmer(create_pp, 'low')
compare_low.run(kernel_type='round_exp')  # round_exp, round_gauss, elliptical not yet implemented
compare_low.run(kernel_type='round_gauss')

compare_high = Vollmer(create_pp, 'high')
compare_high.run(kernel_type='round_exp')
compare_high.run(kernel_type='round_gauss')

print("--- %s seconds ---" % (time.time() - start_time))
