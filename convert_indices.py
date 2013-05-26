import numpy as np 

def convert_indices(direction,x,y,i,j):
	'''
	Converts indices between Python and Fortran indexing, assuming that
	Python indexing begins at 0 and Fortran (for x and y) begins at 1.
	In Tracmass, the vertical indexing does begin at zero so this script
	does nothing to vertical indexing.

	Usage:
		For before a call to tracmass:
			xstart,ystart,ia,ja = convert_indices('py2f',xstart,ystart,ia,ja)
		For after a call to tracmass:
			xend,yend,iend,jend = convert_indices('f2py',xend,yend,iend,jend)
	'''

	if direction == 'py2f':
		x = x+1
		y = y+1
		i = i+1
		j = j+1
	elif direction == 'f2py':
		x = x-1
		y = y-1
		i = i-1
		j = j-1
	return x,y,i,j
