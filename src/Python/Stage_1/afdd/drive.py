#==================================================================
#
# DRIVE SYSTEM
#
# Gear boxes, rotor shafts, drive shafts, rotor brake
# 
#==================================================================

from numpy import pi
from conversions import *

#==================================================================
# REQUIRED INPUTS 
#==================================================================

omega_eng   = 6000  # (rpm)
#p_DSlimit   = 1000  # drive system rated power (hp)
f_rs        = 0.1    # fraction of wt for rotor shaft in (gearbox+rotor shaft) weight
#num_ds      = 1     # number of intermediate drive shafts
#f_q         = 0.03  # second rotor torque (main or tail)
imodel      = 'afdd83'
#imodel      = 'afdd00'

#==================================================================
# FIXED VALUES
#==================================================================

def drivesys_weight(vehicle_parameters):
    

#==================================================================
# unpack and process inputs
#==================================================================
   
   nrotor       = vehicle_parameters['nrotor']
   wblade       = vehicle_parameters['wblade']
   len_ds       = vehicle_parameters['len_ds']
   q_DSlimit    = vehicle_parameters['q_DSlimit']   # torque limit for Xmsn
   P_DSlimit    = vehicle_parameters['P_DSlimit']
   N_at         = vehicle_parameters['nprop']       # propeller count 
   omega        = vehicle_parameters['omega']
   fac          = vehicle_parameters['fac']

   omega        = omega*30.0/pi                     # in RPM

#==================================================================
# limit torque 
#==================================================================

#   q_DSlimit    = q_DSlimit/1.2                     # some factor based on NDARC
#==================================================================
# number of intermediate driveshafts from engine to TR: 4 for SMR
# configuration with TR driveshaft
#==================================================================

   N_ds         = 4
   N_gb         = N_ds + 1      # number of gearboxes
   if nrotor == 1:
     f_p         = 15.0
     f_q         = 3.0
   elif nrotor == 2:
     f_p         = 60.0
     f_q         = f_p
   else:
     f_p         = 10.0 + 100.0/(nrotor)
     f_q         = f_p

#==================================================================
# AFDD00 model for gear box and rotor shaft (gbrs)
#==================================================================

   if imodel=='afdd00':
     w_gbrs = ( 95.7634 * nrotor**0.38553 * 
                P_DSlimit**0.7814 * omega_eng**0.09889/omega**0.8069)

#==================================================================
# AFDD83 model 
#==================================================================

   elif imodel == 'afdd83':
     w_gbrs  = 57.72*(P_DSlimit/1.2)**0.8195*f_q**0.068*N_gb**0.0663*     \
               (omega_eng*0.001)**0.0369/(omega**0.6379)

   else:
     quit('THERE ARE MANY BUTTONS ON THE iPAD...')


#==================================================================
# the following lines are only a breakdown into gearbox (gb) and 
# rotor shaft (rs) based on a fraction "f_rs"; i.e. 
# rotor shaft weight = f_rs * weight of (rotor shaft + gearbox)
# its more for checking parts than anything else
#==================================================================

   wght_gb = (1.0-f_rs)*w_gbrs
   wght_rs = f_rs*w_gbrs

#==================================================================
# weight of driveshaft is based on maximum torque, length and #shafts
# there are some scaling factors based on other parameters
# for SMR only, tail rotor is there
#==================================================================

   if nrotor == 1:
     wght_ds   = (1.166*q_DSlimit**0.3828*len_ds**1.0455  *
                  N_ds**0.3909*(0.01*f_p)**0.2693)
   else:
     wght_ds   = 0.0 
     
#==================================================================
# added drive shaft weight for propeller, set fp = 100% i.e. DS can use 
# 100% of rated power and split equally between all propellers
#==================================================================

   if N_at > 0:
     f_p           = 100.0/N_at 
     N_ds          = N_at 
     wght_ds2 = (1.166*q_DSlimit**0.3828*len_ds**1.0455 * 
                  N_ds**0.3909*(0.01*f_p)**0.2693)

     wght_ds       = wght_ds + wght_ds2

#==================================================================
# rotor brake weight is approximated as 4% weight of blades 
# which is not YUGE; the approximation is obtained by setting Vtip
# to 700 ft/s in the rotor brake eqn. on page 158 of NDARC manual
#==================================================================
	
   wght_rb = 0.0426 * wblade 

#==================================================================
# Apply tech factors and calculate total
#==================================================================

   wght_gb = wght_gb*fac
   wght_rs = wght_rs*fac
   wght_ds = wght_ds*fac
   wght_rb = wght_rb*fac
   
   total   = wght_gb + wght_rs + wght_ds + wght_rb

#==================================================================
# end of operations: return dictionary with breakdown and totals
#==================================================================

   drive_system   = {     'gearbox': wght_gb*lb2kg, 'rotor_shafts': wght_rs*lb2kg,  \
                     'rotor_brakes': wght_rb*lb2kg,   'tail_shaft': wght_ds*lb2kg   }
   return drive_system

#====================
# derived quantities
#====================

#   omega_rotor  = V_tip/R*(60.0/(2.0*pi))
#   q_DSlimit    = p_DSlimit/omega_rotor # (hp/rpm)    
#   len_ds       = 1.2*R

#====================
# unpack input dict
#====================

#   nblade       = vehicle_parameters['nblade']
#   V_tip        = vehicle_parameters['vtip']
#   N_ds         = vehicle_parameters['N_ds']      # number of driveshafts
#   chord        = vehicle_parameters['chord']

#==================================================================
# Blade weights: used to determine rotor system brake weights
#==================================================================

   # if imodel=='afdd00':
   #    wght_blade = ( 0.0024419      * 1.0               *
   #                   nrotor         * nblade**0.53479   *
   #                   R**1.74231     * chord**0.77291    *
   #                   V_tip**0.87562 * nu_blade**2.51048 )
   # elif imodel=='afdd82':
   #    wght_blade = ( 0.02606*nrotor  * nblade**0.6592   *
   #                   R**1.3371       * chord**0.9959    *
   #                   V_tip**0.6682   * nu_blade**2.5279 )