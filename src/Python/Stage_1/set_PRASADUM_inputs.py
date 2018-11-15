#=========================================================================   
#
# Python function to setup blade properties input file for PRASADUM
# Note; the structural properties are dummy inputs since this file is 
# only meant to run a pure aerodynamics problem! Some dummy numbers are 
# inserted for the various structural properties. PRASADUM should only be
# run with the rigid rotor flag enabled for this case.
#
#=========================================================================   
#  Inputs:
#=========================================================================   
#
#  x     : where bilinear twist starts (nondiml, 0 to 1)
#  twx   : twist at nondiml location "x", deg. Positive for nose-up
#  tw1   : twist at tip of rotor blade,   deg. Positive for nose-up
#  y     : where bilinear taper starts (nondiml, 0 to 1)
#  crd0  : chord at root cut-out, ft.
#  crdy  : chord at nondiml location "y", feet
#  crd1  : chord at tip of rotor blade,  feet. 
#  rcout : root cut-out (nondiml, 0 to 1)
#  path  : Where to put the file "BladeData"
#
#=========================================================================   

import os 
import numpy
import sys
def set_csd_blade(x, twx, tw1, y, crd0, crdy, crd1, rcout, path):

#=========================================================================   
# Set file name and open file through context manager
#=========================================================================   

#   print crd0,crdy,crd1
#   quit('what is crdy?')

   filename       = os.path.join(path,'BladeData')
   with open(filename,'w') as f:

#=========================================================================   
# Write header information
#=========================================================================   

      f.write('1 \n')               # number of blades in database = 1

#=========================================================================   
# Write five dummy lines
#=========================================================================   

      write_filler(f, 5)

#=========================================================================   
# Blade numbering, airfoil id
#=========================================================================   

      f.write(" 00001               - Blade numerical identifier \n")
      f.write(" Prototype blade design                           \n")
      f.write(" 1  \n")                 # number of airfoils
      f.write(" 1  \n")                 # airfoil id: SC1095
      f.write(" 1.01e0  \n")            # end radial pos of 1st airfoil
      f.write(" 0.0e0   \n")            # Geometric pitch offset (deg)

#=========================================================================   
# Write structural properties
#=========================================================================   
      
      write_unif_prop(1.852e5,f)      # Flap bending stiffness
      write_unif_prop(1.852e5,f)      #  Lag bending stiffness
      write_unif_prop(1.852e5,f)      # Tors stiffness
      write_unif_prop(0.012e0,f)      # Mass / unit span
      write_unif_prop(0.007e0,f)      # rad. gyr. 1
      write_unif_prop(0.001e0,f)      # rad. gyr. 2

#=========================================================================   
# Prepare geometric twist information
#=========================================================================   

      twst        = numpy.zeros(3); xtwst       = numpy.zeros(3)

      xtwst[0]    = 0.0        ;    twst[0]     = 0.0       # root
      xtwst[1]    = x          ;    twst[1]     = twx       # twist change location
      xtwst[2]    = 1.0        ;    twst[2]     = tw1       # tip

#=========================================================================   
# Write twist information
#=========================================================================   

      write_aero_prop(xtwst, twst, f)

#=========================================================================   
# Inject structural properties (dummy values)
#=========================================================================   
   
      write_unif_prop(0.0, f)          # CG-EA offset
      write_unif_prop(0.0, f)          # PA-EA offset

#=========================================================================   
# Prepare chord information
#=========================================================================   

      crd         = numpy.zeros(5); xcrd        = numpy.zeros(5)

      xcrd[0]     = 0.0        ;    crd[0]      = 0.0
      xcrd[1]     = rcout-0.01 ;    crd[1]      = 0.0       # before cutout
      xcrd[2]     = rcout+0.01 ;    crd[2]      = crd0      # after  cutout
      xcrd[3]     = y          ;    crd[3]      = crdy      
      xcrd[4]     = 1.0        ;    crd[4]      = crd1

#      print crdy, crd[3]
#      print xcrd, crd
#      quit('what is crdy')
#=========================================================================   
# Write chord information
#=========================================================================   

      write_aero_prop(xcrd, crd, f)

#      print xcrd,crd
      # print crdy
#      print 'did aero chord write properly?'
#      quit()
#=========================================================================   
# Axial elongation stiffness
#=========================================================================   

      write_unif_prop(4.2672e7, f)

#=========================================================================   
# End of operations
#=========================================================================   

   return None

#=========================================================================   
#
# Python function to setup rotor properties input file for PRASADUM
#
#=========================================================================   
#  Inputs:
#=========================================================================   
#
#  R     : rotor radius, ft
#  Omega : rotor speed , rad/s
#  rcout : root cut-out (nondiml, 0 to 1)
#  path  : Where the file "RotorData" is located
#
#=========================================================================   

import os 
here        = os.getcwd()
path_loc    = os.path.normpath(os.path.join(here,'../../../py/packages/batchrun/'))
sys.path.insert(0,path_loc)

from input_swap  import input_swap

def set_csd_rotor(Nb, R, Omega, rcout, path):

   tarfile     = os.path.join(path,'RotorData')
   srcfile     = os.path.join(path,'RotorData_template')

   key_str     = ['Rotor speed','Rotor Radius','SPAR LENGTH', 'Number of blades']
   vals        = [Omega, R, R*rcout,Nb ]
   input_swap( srcfile, tarfile, key_str, vals)

   return None

#=========================================================================   
# Function to write "nlines" filler/dummy lines with '=' sign to an open
# file with handle "f"
#=========================================================================   

def write_filler(f,nlines):

   for i in xrange(0,nlines): 
      f.write('=\n')
   return None

#=========================================================================   
# Function to write dummy structural properties to blade properties file
# that has been opened with a file handle "file_handle"
#=========================================================================   

def write_unif_prop(value, f):
   
   write_filler(f,5)             # write the filler before information

   f.write(" 002 \n")             # number of data points
   f.write(" %6.4f   %13.6e    %d \n" % (0.0,value,1))
   f.write(" %6.4f   %13.6e    %d \n" % (1.0,value,1))

   return None

#=========================================================================   
# Function to write blade chord/twist information at "xval" span stations
# with corresponding values "yval" to BladeData file for PRASADUM, open 
# with file handle "f"
#=========================================================================   

def write_aero_prop(xval, yval, f):
   
   write_filler(f,5)             # write the filler before information

   f.write(" %3d \n" % len(xval))   # number of data entries
   for x,y in zip(xval,yval):
      f.write(" %6.4f   %13.6e    %d \n" % (x,y,1))
#      print x,y

#   print 'OK FOR AERO PROPERTIES??'
   return None

#=========================================================================   
# Set linearly twisted and tapered blade
#=========================================================================   

def set_linear_blade(Nb, twist, crd0, crd1, R, omega, path):

#=========================================================================   
# 15% root cut-out
#=========================================================================   

   rcout       = 0.15
   m2ft        = 1.0/0.3048
   Rft         = R*m2ft
   crdroot     = crd0*m2ft
   crdtip      = crd1*m2ft

#=========================================================================   
# Define twist values
#=========================================================================   

   tw1         =-twist                 # in deg, nose up +ve
   x           = 0.2                   # bilinear twist junction
   twx         =-twist*x               # in deg, nose up +ve

#=========================================================================   
# Define chord values
#=========================================================================   

   taper       = crdtip - crdroot
   y           = rcout + 0.2
   crdy        = crdroot + taper*y

#=========================================================================   
# Set blade properties
#=========================================================================   

   set_csd_blade(x, twx, tw1, y, crdroot, crdy, crdtip, rcout, path)

#=========================================================================   
# Set rotor properties
#=========================================================================   

   set_csd_rotor(Nb, Rft, omega, rcout, path)

#=========================================================================   
# End of operations  
#=========================================================================   
   
   return None
