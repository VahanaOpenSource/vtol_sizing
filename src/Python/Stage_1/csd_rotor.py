#=========================================================================   
# Interface for rotor comprehensive analysis : new version
#=========================================================================   
import os, sys

home                 = os.path.expanduser('~') 
batchrun_path        = home + '/Dropbox/py/packages/batchrun/'
sys.path.insert(0,batchrun_path)

from run_sim             import copy_files
from set_PRASADUM_inputs import set_csd_blade, set_csd_rotor
import os, numpy 
def set_CSD_inputs(rotor, blade, path, templ):

#====================================================================
# copy input files from template input directory to actual code 
# input directory
#====================================================================

   copy_files(templ, path)

#=========================================================================   
# Rotor properties
#=========================================================================   

   Nb          = rotor['Nb']                 # number of blades
   R           = rotor['Radius']             # rotor radius, m
   omega       = rotor['omega']              # rotor speed in cruise, rad/s
   omega_h     = rotor['omega_hover']        # rotor speed in  hover, rad/s 

#=========================================================================   
# blade properties
#=========================================================================   

   nub         = blade['flap_freq']          # hover flap freq, /rev
   mb          = blade['mass']               # blade mass, kg
   emax        = 0.05                        # max flap hinge offset

#=========================================================================   
# Compute hinge offset/radius for given hover flap frequency
#=========================================================================   

   if nub == 1.0:
      e           = 0.0
      Kbeta       = 0.0
      m           = mb/R
   elif nub < 1.0:
      print('ERROR: flap natural frequency cannot be less than 1')
      sys.exit(1)
   else:
      e           = 1/(1 + 1.5/(nub*nub-1))        # hinge offset/Radius
      if e <= emax:
         e        = e*R                            # hinge offset, meters
         Kbeta    = 0.0
         m        = mb/(R-e)                       # mass per span, kg/m
      else:
         e        = emax*R                         # hinge offset, meters
         l_bl     = R - e                          # length of blade outboard of hinge, meters
         m        = mb/(R-e)                       # in kg/m
         Ib       = m*l_bl**3/3                    # blade inertia, kg.m^2
         nub0sq   = (nub*nub - 1 - 1.5*e/l_bl)     # non-rotating flap freq. squared
         Kbeta    = Ib*omega_h**2*nub0sq           # root spring stiffness, Nm/rad

#=========================================================================   
# analytical prediction of flap natural freq.
#=========================================================================   

#   anb            = 1+3*e/2/(R-e) + Kbeta/(Ib*omega_h**2)
#   print e/R,1+3*e/2/(R-e),numpy.sqrt(anb)
#   quit()

#=========================================================================   
# Convert hinge offset and flap spring const to FPS units 
#=========================================================================   

   e              = e/0.3048                       # in feet 
   Kbeta          = Kbeta/9.8*2.2/0.3048           # in ft-lb/rad

#=========================================================================   
# At least 15% root cut-out; hinge offset is capped at 5% with spring
#=========================================================================   

   try:
      rcout       = blade['rcout']
   except:
      rcout       = 0.15

#=========================================================================   
# convert length dimensions to FPS units
#=========================================================================   

   m2ft           = 1.0/0.3048
   Rft            = R*m2ft

#=========================================================================   
# Blade geometry for aerodynamics
#=========================================================================   

#=========================================================================   
# Define chord values
# First try linear taper; if entry not found, check for bilinear
#=========================================================================   

   try:

#=========================================================================   
# NOTE: HYDRA Stages 1 and 2 (sizing) operate in SI units for chord, 
# so have to convert to FPS (it also uses linear taper)
#=========================================================================   

      crdroot     = blade['root_chord']*m2ft 
      crdtip      = blade[ 'tip_chord']*m2ft
      taper       = crdtip - crdroot
      y           = rcout + 0.2
      crdy        = crdroot + taper*y

   except:

#=========================================================================   
# NOTE: HYDRA Stage 3 operates in FPS units, no need to convert chord
# and it also uses bilinear twist definition
#=========================================================================   

      crdroot     = blade['crd0']
      crdtip      = blade['crd1']
      y           = blade['y']
      crdy        = blade['crdy']

#=========================================================================   
# Define twist values: first try linear twist, then bilinear defs
#=========================================================================   

   try:

#=========================================================================   
# NOTE: HYDRA Stages 1 and 2 use linear twist definition
# twist is measured positive nose-down for HYDRA, nose-up for 
# physics engine PRASADUM; therefore there is a sign flip. Units are deg.
#=========================================================================   

      twist       = blade['twist']
      tw1         =-twist                 # in deg, nose up +ve
      x           = 0.2                   # bilinear twist junction
      twx         =-twist*x               # in deg, nose up +ve

#=========================================================================   
# NOTE: HYDRA Stage 3 uses bilinear twist definition identical to 
# PRASADUM; no need to reverse sign of twist. Units are degrees
#=========================================================================   

   except:
      x           = blade['x']
      twx         = blade['twx']
      tw1         = blade['tw1']

#=========================================================================   
# Set blade properties
#=========================================================================   

#=========================================================================   
# old version: use for rerunning HYDRA 1.0 cases
#=========================================================================   

#   set_csd_blade(x, twx, tw1, y, crdroot, crdy, crdtip, rcout, path)

#=========================================================================   
# Set rotor properties
#=========================================================================   

#=========================================================================   
# old version: use for rerunning HYDRA 1.0 cases
#=========================================================================   

#   set_csd_rotor(Nb, Rft, omega, rcout, path)

#=========================================================================   
#new version: HYDRA 2.0
#=========================================================================   

   rotor_prop     = {'Nb': Nb, 'Radius': Rft, 'Omega': omega,              \
                     'cut_out': rcout*Rft, 'e': e, 'Kbeta': Kbeta,         \
                     'blade_wt': mb*2.2} # wt in lb for each blade
   csd_rotor(rotor_prop, path, templ)

#=========================================================================   
# Determine blade property dictionary
#=========================================================================   

#   print 'blade mass is ',m*2.2*4*(R-e),'lb'

   m              = m*0.3048/14.5939            # in slug/ft 

   blade_prop     = {'x': x, 'twx': twx, 'tw1': tw1, 'y': y, 'crd0': crdroot, \
                     'crd1': crdtip, 'crdy': crdy, 'rcout': rcout, 'mass': m}

   blade_prop['Radius']   = rotor_prop['Radius']
   blade_prop['Omega']    = rotor_prop['Omega']

   csd_blade(blade_prop, path)

#=========================================================================   
# End of operations  
#=========================================================================   
   
   return None

#=========================================================================   
#
# Python function to setup rotor properties input file for PRASADUM
#
#=========================================================================   
#  Input: path and dictionary containing following possible entries
#=========================================================================   
#
#  R     : rotor radius, ft
#  Omega : rotor speed , rad/s
#  rcout : root cut-out, ft
#  e     : hinge offset, ft
#  Kbeta : root flap spring constant, ft-lb/rad
#  path  : Where the file "RotorData" is located
#
#=========================================================================   

from input_swap  import input_swap

def csd_rotor(rotor, path, templ):

#   print 'PATH IS ',path
   tarfile     = os.path.join( path,'RotorData')
   srcfile     = os.path.join(templ,'RotorData_template')

   key_str     = []     
   vals        = []

#=========================================================================   
# Rotor radius
#=========================================================================   

   try:
      R        = rotor['Radius']  ; key_str.append('Rotor Radius'); vals.append(R)
   except:
      print('WARNING: ROTOR RADIUS NOT PART OF DICTIONARY: REVERTING TO DEFAULTS')
      pass 

#=========================================================================   
# Rotor speed
#=========================================================================   

   try:
      Omega    = rotor['Omega']   ; key_str.append('Rotor speed') ; vals.append(Omega)
   except:
      print('WARNING: ROTOR RPM NOT PART OF DICTIONARY: REVERTING TO DEFAULTS')
      pass 

#=========================================================================   
# Spar length
#=========================================================================   

   try:
      e        = rotor['e']      ; key_str.append('HINGE OFFSET'); vals.append(e)
      e2       = rotor['cut_out']; key_str.append('SPAR LENGTH') ; vals.append(e2-e)
     # print e2, e
   except:
      print('WARNING: ROTOR ROOT CUT OUT AND HINGE OFFSET NOT PART OF DICT: SWITCH TO DEFAULT')
      pass

#=========================================================================   
# Number of blades
#=========================================================================   

   try:
      Nb       = rotor['Nb']     ; key_str.append('Number of blades'); vals.append(Nb)
   except:
      print('WARNING: # BLADES NOT PART OF DICT: SWITCH TO DEFAULT')
      pass 

#=========================================================================   
# Blade root flap spring constant
#=========================================================================   

   try:
      Kbeta    = rotor['Kbeta']  ; key_str.append('Blade Flap Spring'); vals.append(Kbeta)
   except:
      print('WARNING: # FLAP SPRING CONSTANT NOT PART OF DICT: SWITCH TO DEFAULT')
      pass                      

#=========================================================================   
# blade weight in lbs
#=========================================================================   

   try:
      Wb       = rotor['blade_wt']; key_str.append('BladeWt'); vals.append(Wb)
   except:
      print('Warning: blade weight not found in dictionary; using default')
      pass 
      
#=========================================================================   
# call input swap function (AS)
#=========================================================================   

   input_swap( srcfile, tarfile, key_str, vals)

#=========================================================================   
# End of operations
#=========================================================================   

   return None

#=========================================================================   
#
# Python function to setup blade properties input file for PRASADUM
#
#=========================================================================   
#  Input: path and dictionary containing following possible entries
#=========================================================================   
#
#  x     : bilinear twist junction, nondiml span
#  twx   : twist at junction, deg
#  tw1   : twist at tip, deg
#  y     : bilinear taper junction, nondiml span
#  crd0  : chord at root cut-out, ft
#  crd1  : chord at tip, ft
#  rcout : root cut-out, nondiml span
#  mass  : mass per unit span of blade, slug/ft
#  path  : Where the file "BladeData" needs to go
#
#=========================================================================   

from set_PRASADUM_inputs import write_filler, write_unif_prop, write_aero_prop
import os 
def csd_blade(blade, path):

#=========================================================================   
# blade properties
#=========================================================================   

   x              = blade['x']
   twx            = blade['twx']
   tw1            = blade['tw1']
   y              = blade['y']
   crd0           = blade['crd0']
   crdy           = blade['crdy']
   crd1           = blade['crd1']
   rcout          = blade['rcout']
   m              = blade['mass']            #blade mass, kgs
#   m              = m*0.0685218/
   print('blade mass/span is ',m,'slug/ft')
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

      Omega    = blade['Omega']
      R        = blade['Radius']
      EIflap   = 10*m*Omega*Omega*R*R*R*R       # use normalized properties to ensure stiff blade      
      write_unif_prop(EIflap,    f)       # Flap bending stiffness
      write_unif_prop(EIflap*10, f)       #  Lag bending stiffness
      write_unif_prop(EIflap*5,  f)       # Tors stiffness
      write_unif_prop(      m, f)       # Mass / unit span
      write_unif_prop(0.007e0, f)       # rad. gyr. 1
      write_unif_prop(0.001e0, f)       # rad. gyr. 2

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
