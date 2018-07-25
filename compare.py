#For reading config files, based on key = value structure,
import configparser

#Ability to call external programs, shell commands etc
import os

#Library to read and write csv files
import csv


#Command to search for all files ending in 'results.ini' below working directory './'
#works recursively!
datafiles = [os.path.join(root, name)
			for root, dirs, files in os.walk('./')
			for name in files
			if name.endswith("results.ini")]

#Result is stored as a list of strings
print datafiles

#now use configparsers to read in all the files and print all their key : value pairs
for item in datafiles:
	config = configparser.ConfigParser()
	config.read(item)
	print config.sections()
	sect = config.sections()
	sect.remove('name')
	for item in sect:
		#Two different ways to format number output
		for key in config[item]: print key, ' = ','%.2f' % float(config[item][key])
		#for key in config[item]: print key, ' = ',round(float(config[item][key]),2)
### Now we have all the output data als configparsers lists of strings.
#to format use 
