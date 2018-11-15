#======================================================================
# create range of blade twist, chord distribution (bilinear w/ limits)
#======================================================================
import numpy
import copy

def blade_properties():

#======================================================================
# Chord distribution
#======================================================================

#ACTUAL
#   cmin        = 0.25               # min chord/avg. chord
#   cmax        = 3.0                # max chord/avg. chord

#TESTING
   cmin        = 1.0               # min chord/avg. chord
   cmax        = 1.0                # max chord/avg. chord
   dc          = 0.25               # interval size
   nc          = numpy.floor((cmax-cmin)/dc) + 1

   ctip_array  = numpy.linspace(cmin, cmax, nc)
   croot_array = copy.copy(ctip_array)

#ACTUAL
#   ymin        = 0.3
#   ymax        = 0.8
#TESTING
   ymin        = 0.8
   ymax        = 0.8
   dy          = 0.1

   ny          = numpy.floor((ymax-ymin)/dy) + 1
   y_array     = numpy.linspace(ymin,ymax,ny)

#======================================================================
# twist array
#======================================================================

#ACTUAL
#   thmin       = -20.0              # min twist, deg, +ve nose up
#   thmax       =   2.5              # max twist, deg, +ve nose up
#TESTING
   thmin       = -10.0
   thmax       = -10.0

   dth         =   2.5              # twist angle interval size, deg

   nth         = numpy.floor((thmax - thmin)/dth) + 1
   thx_array   = numpy.linspace(thmin,thmax,nth)     
   thtip_array = numpy.linspace(thmin,thmax,nth)
   
   x_array     = copy.copy(y_array)
   
   prop_dict   = {'thx_all': thx_array, 'thtip_all': thtip_array,    
                    'x_all': x_array  ,  'ctip_all':  ctip_array,
                    'y_all': y_array  , 'croot_all': croot_array}

   return prop_dict

#======================================================================
# Compute chord at bilinear taper junction from avg. chord, root chord and tip chord 
#======================================================================

def get_cy(cbar, ctip, croot, y):

   cy          = 2.0*cbar - croot * y - ctip * (1.0 - y)
   if cy <= 0:
      print('BANG: chord cannot be negative!')
      print('mean chord: ',cbar)
      print('tip  chord: ',ctip)
      print('root chord: ',croot)
      print('taper junc: ',y)
      print('chord ther: ',cy)
      quit()
   return cy