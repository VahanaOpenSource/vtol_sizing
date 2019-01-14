import numpy
def footprint_calc(spans, l_fus):


	y1 		= spans[0]*0.5
	try:
		y2 		= spans[1]*0.5
	except:
		y2 		= 0.0

	l 		= (y2*y2 + l_fus*l_fus-y1*y1)/(2.0*l_fus)

#circle radius enclosing both surfaces
	R 		= numpy.sqrt(l*l+y1*y1)

	footprint 	= 2*R
	return footprint