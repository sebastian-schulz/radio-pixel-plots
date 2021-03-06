# For reading config files, based on key = value structure,
import configparser
# Ability to call external programs, shell commands etc
import os
# System tools (file in and output ...)
import sys

if len(sys.argv) != 2:
    print('Need path to parent directory of results!')
    sys.exit(-1)

parent = str(sys.argv[1])
# Command to search for all files ending in 'results.ini' below working directory './'
# Result is stored as a list of strings
# works recursively!
datafiles = [os.path.join(root, name)
             for root, dirs, files in os.walk(parent)
             for name in files
             if name.endswith("results.ini")]
datafiles.sort()

# Now print all results to file as tab separated values,
# starting with column descriptions
for i in range(len(datafiles)):
    print(datafiles[i])

outfile = open(parent + '/combined_results.dat', 'w')
outfile.write('Galaxy name\t' +
              'low\t' +
              'a \t' +
              'a_err\t' +
              'b \t' +
              'b_err\t' +
              'chi_sqr\t' +
              'l low\t' +
              'a err\t' +
              'high\t' +
              'a\t' +
              'a_err\t' +
              'b\t' +
              'b_err\t' +
              'chi_sqr\t' +
              'l high\t' +
              'a err\n')

# Now write the actual data
for item in datafiles:
    config = configparser.ConfigParser()
    config.read(item)
    sect = config.sections()
    outfile.write(config['name']['fname'] + '\t' +
                  '%.2f' % float(config['frequencies']['freq_low']) + '\t' +
                  '%.2f' % float(config['Low_freq_fit']['a']) + '\t' +
                  '%.2f' % float(config['Low_freq_fit']['a_err']) + '\t' +
                  '%.2f' % float(config['Low_freq_fit']['b']) + '\t' +
                  '%.2f' % float(config['Low_freq_fit']['b_err']) + '\t' +
                  '%.2f' % float(config['Low_freq_fit']['chi_sqr']) + '\t' +
                  '%.2f' % float(config['Conv_results']['sigma_l']) + '\t' +
                  '%.2f' % float(config['Conv_results']['a_err_l']) + '\t' +
                  '%.2f' % float(config['frequencies']['freq_high']) + '\t' +
                  '%.2f' % float(config['High_freq_fit']['a']) + '\t' +
                  '%.2f' % float(config['High_freq_fit']['a_err']) + '\t' +
                  '%.2f' % float(config['High_freq_fit']['b']) + '\t' +
                  '%.2f' % float(config['High_freq_fit']['b_err']) + '\t' +
                  '%.2f' % float(config['High_freq_fit']['chi_sqr']) + '\t' +
                  '%.2f' % float(config['Conv_results']['sigma_h']) + '\t' +
                  '%.2f' % float(config['Conv_results']['a_err_h']) + '\n')
