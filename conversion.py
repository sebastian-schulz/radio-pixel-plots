#Mighty numerical library of ptyhon
import numpy as np

#Pythons math package with all the constants and functions
import math as m

#number of pixels along one direction of a box to achieve 1200pc length (float)
def conv_px_per_box( cfg ):
	tmp = 1000. * cfg.getfloat('kpc') / ( cfg.getfloat('distance') * m.tan( m.radians(1./3600.) ) / cfg.getfloat('pixel_per_arcsec') )
	return float( tmp )

#Calculate RMS inside the box defined using DS9
def calculate_rms( data, cfg ):
	rms = 0.
	#calculate borders of the rms box in pixel, makes use of int cutting
	#therefor the box may be 1 pixel length off in any direction
	#but that only has a small impact
	y_min = cfg.getint('center_x') - cfg.getint('size_x') / 2
	y_max = cfg.getint('center_x') + cfg.getint('size_x') / 2
	x_min = cfg.getint('center_y') - cfg.getint('size_y') / 2
	x_max = cfg.getint('center_y') + cfg.getint('size_y') / 2
	n = cfg.getint('size_x') * cfg.getint('size_y')
<<<<<<< HEAD
	for i in range(int(x_min), int(x_max) ):
		for j in range(int(y_min),int(y_max) ):
=======
	for i in range( int(x_min), int(x_max) ):
		for j in range( int(y_min), int(y_max) ):
>>>>>>> 39e6c05434ac1f792f6c59921ecae3e02b27653c
			rms += m.pow( data[i][j] , 2)
	rms = m.sqrt( rms / n )
	#print('ROOT MEAN SQUARE OF BOX: ', rms)
	return rms

###Converts 2d array (image) to pixels of 1.2kpc length and cuts it to the given box size; basically produces a smaller and lower-resolution version
def convert1200( data, cfg ):
	px_per_box = conv_px_per_box( cfg )
	'''calculate starting positions (bottom left corner) in pixels'''
	s_x = cfg.getint('center_x') - cfg.getint('n_boxes') / 2 * int( px_per_box )
	s_y = cfg.getint('center_y') - cfg.getint('n_boxes') / 2 * int( px_per_box )
	#Creating empty array to store newly created image
	pixels = np.zeros(shape = ( cfg.getint('n_boxes'), cfg.getint('n_boxes') ) )

	#Loop over the 2d array. This makes use of cutting decimals by int-cutting
	#saves some computation time
	for x in range( cfg.getint('n_boxes') * int( px_per_box ) ):
		for y in range( cfg.getint('n_boxes') * int( px_per_box ) ):
			pixels[ x / int( px_per_box ), y / int( px_per_box ) ] += data[ x + s_x ][ y + s_y ]
	#Return average values (total signal divided by number of old pixels per new pixel)
	return pixels / ( px_per_box )**2

#PASSED BUGTEST (July 2018), gives the same results as DS9 with the same boxes
def convert_resolution_adv( data, cfg ):
	''' Idea: allow for pixels_per_box to be float!
assumption 1 : middle pixel for p and p' have the same position (center)
calculate everything in p distances: (n_boxes should be odd!)
bottom left x = center_x - px_per_box * n_boxes/2
same for y
loop over all p'
	loop over p inside p'
		go to bottom left corner,
		find pixel p fully inside,
		same for top right corner,
		-> calculate sum of all pixels fully inside,
        
		now borders:
		4 corners
		2 edges horizontal
		2 edges vertical
'''
	px_per_box = conv_px_per_box( cfg )
	#position of center of bottom left new pixel in old pixels
	#for some reason x and y axis are in the wrong order in the array
	img_bot_left_y = cfg.getint('center_x') - 1 - (cfg.getint('n_boxes') / 2. - 0.5) * px_per_box
	img_bot_left_x = cfg.getint('center_y') - 1 - (cfg.getint('n_boxes') / 2. - 0.5) * px_per_box
	#new pixel array
	N=0
	pixels = np.zeros(shape = ( cfg.getint('n_boxes'), cfg.getint('n_boxes') ) )
	for x in range( cfg.getint('n_boxes') ):
		for y in range( cfg.getint('n_boxes') ):
			#position of current pixel in old pixels
			img_current_x = img_bot_left_x + x * px_per_box
			img_current_y = img_bot_left_y + y * px_per_box
			#positions of the corners of the new pixels in old pixels
			pix_x_min = img_current_x - px_per_box / 2.
			pix_x_max = img_current_x + px_per_box / 2.
			pix_y_min = img_current_y - px_per_box / 2.
			pix_y_max = img_current_y + px_per_box / 2.
			#loop over the old pixels fully inside the current new pixel
			for p_x in range( ceil( pix_x_min ), floor( pix_x_max ), 1):
				for p_y in range( ceil( pix_y_min ), floor( pix_y_max ), 1):
					pixels[ x, y ] += data[ p_x, p_y ]
					N += 1
			#now loop over the edges, without the corners
			factor_max_x = m.fabs( floor(pix_x_max) - pix_x_max )
			factor_min_x = m.fabs( ceil(pix_x_min) - pix_x_min )
			factor_max_y = m.fabs( floor(pix_y_max) - pix_y_max )
			factor_min_y = m.fabs( ceil(pix_y_min) - pix_y_min )
			#if (factor_max_x > 1 or factor_max_y > 1 or factor_min_x > 1 or factor_min_y > 1):
			#	print( 'ERROR!')
			#print('x and y factors \t', factor_max_x, '\t', factor_max_y)
			for p_y in range( ceil( pix_y_min ), floor( pix_y_max ), 1):
				N += factor_min_x + factor_max_x
				pixels[ x, y ] += data[floor( pix_x_min ) ,p_y ] * factor_min_x
				pixels[ x, y ] += data[ceil( pix_x_max ) ,p_y ] * factor_max_x
			for p_x in range( ceil( pix_x_min ), floor( pix_x_max ), 1):
				N+= factor_min_y + factor_max_y
				pixels[ x, y ] += data[ p_x, floor( pix_y_min ) ] * factor_min_y
				pixels[ x, y ] += data[ p_x, ceil( pix_y_max ) ] * factor_max_y

			#finally the four corner pixels:
			pixels[ x, y ] += data[ floor(pix_x_min), floor(pix_y_min) ] *factor_min_x * factor_min_y
			pixels[ x, y ] += data[ floor(pix_x_min), ceil(pix_y_max) ] * factor_min_x *factor_max_y
			pixels[ x, y ] += data[ ceil(pix_x_max), ceil(pix_y_max) ] * factor_max_x * factor_max_y
			pixels[ x, y ] += data[ ceil(pix_x_max), floor(pix_y_min) ] * factor_max_x * factor_min_y
			N += factor_min_x * factor_min_y + factor_min_x *factor_max_y +factor_max_x * factor_max_y + factor_max_x * factor_min_y
	#if(PRINTALL==True ):
	#print( 'Total number of original pixels used in conv:\t', N)
	return pixels / (px_per_box)**2

#Custom defined ceiling and floor functions with int as return type
def ceil(x):
	return int(m.ceil(x))

def floor(x):
	return int(m.floor(x))


