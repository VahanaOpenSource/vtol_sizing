#=========================================================================   
# Write file header information
#=========================================================================   

def write_header(f):
   string1              = '{:^15} {:^15} {:^15}'.format('Bilin.Twst.Jn','Jn.Twist(deg)','Tip.Twst(deg)')
   string2              = '{:^15} {:^15} {:^15}'.format('Bilin.Tapr.Jn','Rt.Chord(ft)' ,'Tip.Chord(ft)')

   string3              = '{:^15} {:^15} {:^15}'.format('Power(hp)','Fx hub (lbs)','Total Pwr (hp)')
   string               = string1 + string2 + string3 + '\n'
   f.write(string)

   return None

#=========================================================================   
# Write file header information
#=========================================================================   

def write_to_file(com, CSD_outputs, f):

#=========================================================================   
# Twist properties
#=========================================================================   

   Twist_Jn_Location    = com[0]
   Jn_Blade_Twist       = com[1]
   Tip_Twist            = com[2]

#=========================================================================   
# Taper properties
#=========================================================================   

   Taper_Jn_Location    = com[3]
   Root_Chord           = com[4]
   Tip_Chord            = com[5]

#=========================================================================   
# Flight speed (ft/s)
#=========================================================================   

   Vinf                 = CSD_outputs['Vinf']                       # in m/s

#=========================================================================   
# Prepare strings to write to file
#=========================================================================   

   string1              = '{:^15.2f} {:^15.4f} {:^15.2f}'.format(Twist_Jn_Location,  Jn_Blade_Twist,   Tip_Twist)
   string2              = '{:^15.2f} {:^15.2f} {:^15.2f}'.format(Taper_Jn_Location,      Root_Chord,   Tip_Chord)

   Power                = CSD_outputs['Power']                    # in watts
   Fxhub                = CSD_outputs['Fx']                       # in Newtons

#   print Vinf
#   print Power/746, Fxhub*Vinf/0.85/746

#=========================================================================   
# Compute   total power = Turn rotor + Power to overcome blade drag 
#=========================================================================   

   TotalPower           = (Power + Fxhub * Vinf / 0.85) / 746     # total power, Hp

#=========================================================================   
# Prepare output strings to write to file
#=========================================================================   

   string4              = '{:^15.2f} {:^15.2f} {:^15.2f}'.format(Power/746, Fxhub/9.8*2.2, TotalPower)

#=========================================================================   
# Write data to file
#=========================================================================   

   line                 = string1 + string2 + string4 + '\n'
   f.write(line)
   
   return None

#=========================================================================   
# Break line of data into two dictionaries (input and output )
#=========================================================================   

def create_subset(array):

   subset                  = {};    data = {}

#=========================================================================   
# "Input"
#=========================================================================   

   subset['Twist_Jn_Location']   = array[0]
   subset['Jn_Blade_Twist']      = array[1]
   subset['Tip_Twist']           = array[2]


   subset['Taper_Jn_Location']   = array[3] 
   subset['Root_Chord']          = array[4]
   subset['Tip_Chord']           = array[5]

#=========================================================================   
# "Output"
#=========================================================================   

   quit('HLLO: PERFORMANCE SWEEP')
#   data['Wt']              = [array[7],0]
#   data['Power']           = [array[8],0]
#   data['Radius']          = [array[9],0]
#   data['lift_off']        = [array[10],0]
#   data['span']            = [array[11],0]
#   data['RTF']             = [array[12],0]
#   data['Fuel']            = [array[13],0]
#   data['Empty']           = [array[14],0]
#   data['Wlift']           = [array[15],0]

   return subset, data

#=========================================================================   
# Call CSD solver
#=========================================================================   

from set_PRASADUM_inputs import set_csd_rotor #, set_csd_blade
from blade_properties    import get_cy
from run_CSD             import run_CSD
from csd_rotor           import set_CSD_inputs 
import os
def CSD_caller(com, CSD_params):

#=========================================================================   
# Extract quantities for rotor input file
#=========================================================================   

   Nb                      = CSD_params['Nb']
   R                       = CSD_params['Radius'] / 0.3048
   rcout                   = 0.15
   Omega                   = CSD_params['Omega_c']

#=========================================================================   
# Create rotor data input file
#=========================================================================   

   home                    = os.path.expanduser('~') 
   path                    = home + '/Dropbox/PrasadUM/Inputs/'
#   set_csd_rotor(Nb, R, Omega, rcout, path)

#=========================================================================   
# Extract quantities for blade input file
#=========================================================================   

   x                       = com[0]
   twx                     = com[1]
   tw1                     = com[2]

   y                       = com[3]
   crd0                    = com[4]
   crd1                    = com[5]

   cbar                    = CSD_params['cbar']

   crdy                    = get_cy(cbar, crd1, crd0, y)

#=========================================================================   
# Extract quantities for CSD input files: SI units!
#=========================================================================   

   T                       = CSD_params['Thrust']
   L                       = CSD_params['L']
   M                       = CSD_params['M']
   Vinf                    = CSD_params['Vinf']
   alpha_s                 = CSD_params['alpha']
   density_alt             = CSD_params['density_alt']
   hover_omega             = CSD_params['Omega_h']

#====================================================================
# Set rotor properties dictionary
#====================================================================

   rotor_properties        = {'Radius': R,      'Nb': Nb,         \
                               'omega': Omega,                    \
                               'omega_hover': hover_omega}

#====================================================================
# Blade properties dictionary
#====================================================================

#====================================================================
#old version: blade properties
#====================================================================
#   set_csd_blade(x, twx, tw1, y, crd0, crdy, crd1, rcout, path)

#====================================================================
# New version: extract from dictionary
#====================================================================

   twist                   = 10                      # nose down, deg
   nu_beta                 = CSD_params['nu_beta']
   blade_mass              = CSD_params['mbl']

#====================================================================
# Create dictionary
#====================================================================

# removed keys 'twist', 'root_chord', 'tip_chord'
#added    keys 'x', 'twx', 'tw1', 'y', 'crd0', 'crdy', crd1', 'rcout' 
#
   blade_properties        = {         'x' :    x,       'y' : y,       \
                                     'twx' :  twx,     'tw1' : tw1,     \
                                    'crd0' : crd0,    'crd1' : crd1,    \
                                    'crdy' : crdy,   'rcout' : rcout,   \
                               'flap_freq' : nu_beta,                   \
                               'mass'      : blade_mass     }

#====================================================================   
# set CSD inputs
#====================================================================   

   set_CSD_inputs(rotor_properties, blade_properties, path)

#====================================================================   
# Run CSD code
#====================================================================   

   info, Power, Fx         = run_CSD(T, L, Vinf, alpha_s, density_alt)

   if info != 1:
      Power                = 1.e6
      Fx                   = 1.e6

   CSD_outputs             = {'Power' : Power, 'Fx' : Fx, 'Vinf' : Vinf}

#=========================================================================   
# End of operations
#=========================================================================   
   
   return CSD_outputs