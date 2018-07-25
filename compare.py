#For reading config files, based on key = value structure,
import configparser

#Ability to call external programs, shell commands etc
import os

#System tools (file in and output ...)
import sys

#Regualar expressions in python - yay!
import re

if(len(sys.argv) != 2 ):
	print 'Need path to parent directory of results!'
	sys.exit(-1)

parent = str(sys.argv[1])
#Command to search for all files ending in 'results.ini' below working directory './'
#works recursively!
datafiles = [os.path.join(root, name)
			for root, dirs, files in os.walk(parent)
			for name in files
			if name.endswith("results.ini")]

#Result is stored as a list of strings
for i in range(len(datafiles)):
	print datafiles[i]
outfile = open(parent+'/comparison.dat', 'w')
outfile.write(	'Galaxy name\t'+
		'low\t'+
		'a \t'+
		'a_err\t'+
		'b \t'+
		'b_err\t'+
		'chi_sqr\t'+
		'high\t'+
		'a\t'+
		'a_err\t'+
		'b\t'+
		'b_err\t'+
		'chi_sqr\t'+
		'l low\t'+
		'l high\n')

#now use configparsers to read in all the files and print all their key : value pairs
for item in datafiles:
	config = configparser.ConfigParser()
	config.read(item)
	sect = config.sections()
	outfile.write(	config['name']['fname']+'\t'+'\t'+
			'%.2f' % float(config['Low_freq_fit']['a'])+'\t'+
			'%.2f' % float(config['Low_freq_fit']['a_err'])+'\t'+
			'%.2f' % float(config['Low_freq_fit']['b'])+'\t'+
			'%.2f' % float(config['Low_freq_fit']['b_err'])+'\t'+
			'%.2f' % float(config['Low_freq_fit']['chi_sqr'])+'\t'+
			'\t'+
			'%.2f' % float(config['High_freq_fit']['a'])+'\t'+
			'%.2f' % float(config['High_freq_fit']['a_err'])+'\t'+
			'%.2f' % float(config['High_freq_fit']['b'])+'\t'+
			'%.2f' % float(config['High_freq_fit']['b_err'])+'\t'+
			'%.2f' % float(config['High_freq_fit']['chi_sqr'])+'\t'+
			'%.2f' % float(config['Conv_results']['sigma_l'])+'\t'+
			'%.2f' % float(config['Conv_results']['sigma_h'])+'\n')

	#print config.sections()
	#for item in sect:
		#Two different ways to format number output
		#for key in config[item]: print key, ' = ','%.2f' % float(config[item][key])
		#for key in config[item]: print key, ' = ',round(float(config[item][key]),2)
### Now we have all the output data als configparsers lists of strings.
#to format use 



































