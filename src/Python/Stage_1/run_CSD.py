#=========================================================================   
# Python code to run PRASADUM with appropriate inputs
#
#=========================================================================   
# Inputs:
#=========================================================================   
#
#        (a) Thrust     : target thrust       , Newtons
#        (b) Roll_moment: Hub rolling moment  , Newton-meters
#        (c) Speed      : Flight speed        , m/s
#        (d) shaft_tilt : Rotor fwd shaft tilt, radians
#
#=========================================================================   
# Outputs
#=========================================================================   
#
#        (a) Power      (watts)
#        (b) Hub drag   (Newtons)
#
#=========================================================================   

import os, sys, numpy

#=======================================================================
# Import path where AS routines exist, then import modules from there
#=======================================================================

home                 = os.path.expanduser('~') 
batchrun_path        = home + '/Dropbox/py/packages/batchrun/'
sys.path.insert(0,batchrun_path)
from   run_sim     import  copy_files, run_code
from input_swap    import input_swap

#=======================================================================
# Run rotor dynamics program
#=======================================================================

def run_CSD(Thrust,Roll_moment,Speed,shaft_tilt,density_alt,csd_dirs):

#=========================================================================   
# define CSD locations and copy template folder
#=========================================================================   

   csd                  = csd_dirs['prasadum']

   csd_input_path       = csd + 'Inputs/'
   csd_template_path    = csd + 'Inputs_samples/Rotor_only/HYDRA_v2/'
   csd_exec_path        = csd + 'exec/'
   csd_output_path      = csd + 'Outputs/'

#=========================================================================   
# setup and swap input for flight condition file
#=========================================================================   

   fl_file              = 'Flight_condition.input'
   source               = os.path.join(csd_template_path, fl_file)
   target               = os.path.join(csd_input_path,    fl_file)

#=======================================================================
# debug prints
#=======================================================================

#=======================================================================
# define key strings to look for in Flight condition input file
# find speed in knots, and density alt in feet 
#=======================================================================

   input_strings        = ['Forward Velocity','Altitude']
   spd                  = "%12.5f" % (Speed * 18.0 /5.0 /  1.853)
   string               = spd + '                       FORWARD VELOCITY (KT)'
   string2              = "%12.5f" % (density_alt/0.3048)         # in ft 
   string2              = string2 + '                   DENSITY ALTITUDE (FT)'
   values               = [string,string2]

   input_swap(source, target, input_strings, values)

#=======================================================================
# Set wind-tunnel targets
#=======================================================================

   wt_input_file        = 'Wind_tunnel.input'
   source               = os.path.join(csd_template_path, wt_input_file)
   target               = os.path.join(csd_input_path   , wt_input_file)

#   print source, target
#=======================================================================
# define key strings to look for
#=======================================================================

   input_strings        = [];       values = []

#=======================================================================
# Shaft tilt
#=======================================================================

   input_strings.append('Shaft tilt angle')
   alpha                = "%12.5e" % (-shaft_tilt * 180.0/numpy.pi)     # deg, back
   string0              = alpha + '                      Shaft tilt angle, +ve aft,  [deg]'
   values.append(string0)

#   print 'SENDING VALUE TO FILE',alpha
#=======================================================================
# Target thrust
#=======================================================================

   input_strings.append('Target thrust')
   thr                  = "%12.5e" % (Thrust / 9.8* 2.2)                  # lbs
   string1              = thr   + '                     Target thrust (along shaft)  [lbs]'
   values.append(string1)

#=======================================================================
# Rolling moment
#=======================================================================

   input_strings.append('Target hub rolling moment')
   roll                 = "%12.5e" % (Roll_moment / 9.8 * 2.2 / 0.3048)   # ft-lbs
   string2              = roll  + '                     Target hub rolling moment    [ft-lbs] ROLL LEFT +ve'
   values.append(string2)

#=========================================================================   
# Call input swap function 
#=========================================================================   

   input_swap(source, target, input_strings, values)

#=========================================================================   
# Run the CSD code
#=========================================================================   

   run_code(csd_exec_path, 'prasadum')

#=========================================================================   
# open Trim_02.dat and read power value in Hp 
#=========================================================================   
   
   opfile               = os.path.join(csd_output_path,'Trim_02.dat')
   with open(opfile,'r') as f:
      nlines            = 0
      for line in f:
         temp           = line.rstrip('\n').split()
         nlines         = nlines + 1
         Power          = float(temp[5])*746.0   # Shaft power in watts
         Fx             = float(temp[6])/2.2*9.8 # in Newtons, hub drag
         info           = float(temp[14])        # flag to indicate trim (1) or not (anything else)

#=========================================================================   
# Get force in X direction: Hub drag (Fx) + Thrust * sin(alpha)
#=========================================================================   

         Fx             =     Fx * numpy.cos(shaft_tilt) -                 \
                          Thrust * numpy.sin(shaft_tilt)

#=========================================================================   
# open Trim_02.dat and read power value in Hp 
#=========================================================================   
   
   opfile               = os.path.join(csd_output_path,'Trim_01.dat')
   trim_vars            =  numpy.loadtxt(opfile)
#   print len(trim_vars)
   if len(trim_vars) > 6:
      b0                   = trim_vars[6]
      b1c                  = trim_vars[7]
      b1s                  = trim_vars[8]

      if(abs(b0) > 0.2 or abs(b1c) > 0.2 or abs(b1s) > 0.2 ):
         info              = 0
   #      quit('TOO MUCH FLAPPING')
         Power             = Power + 746.0*1000.0     # add 1000 Hp penalty
#=========================================================================   
# Error trap
#=========================================================================   

   if nlines > 1:
      print('CRITICAL ERROR: CSD CODE RUN AT MULTIPLE FLGHT CONDITIONS')
      print('INSTEAD OF SINGLE CONDITION: PROGRAM TERMINATING')

   if info != 1:
      print('WARNING! CSD CODE MAY NOT HAVE TRIMMED!')
      print(info)
         

#=========================================================================   
# End of operations
#=========================================================================   

#   print Power/746.0, Fx;quit()
   return info, Power, Fx
