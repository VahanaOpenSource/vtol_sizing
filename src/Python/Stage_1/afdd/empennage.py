#====================================================================
# EMPENNAGE GROUP
#
# Horizontal and vertical tail
#
# ALL UNITS IN FPS
#====================================================================
from conversions import *
f_tr             = 1.6311 #if TR on vertical tail, otherwise 1. Mostly 1.63
A_ht             = 4.0  # horizontal tail aspect ratio
A_vt             = 4.0  # vertical tail aspect ratio

def empennage_weight(vehicle_parameters, wing_wt):

    aircraftID = vehicle_parameters['aircraftID']

    wingarea    = vehicle_parameters['wingarea']      # for tilt-rotor: ft^2
    vdive       = vehicle_parameters['vdive']         # for tilt-rotor: knots
    P_DSlimit   = vehicle_parameters['pwr_installed'] # for tail rotor: hp
    N_at        = vehicle_parameters['nprop']         # for compound 
    T_at        = vehicle_parameters['T_prop']        # thruster force, lb 
    R_at        = 0.5                              # assume % rotor radius
    A_at        = pi*R_at*R_at
    fac         = vehicle_parameters['tech_factors'].empennage 
    wing        = vehicle_parameters['wing']

#====================================================================
# Vtail has same area as htail, 4.9 kg/sq.m
#====================================================================

    min_area    = 1e9 
    for ig,group in wing.groups.items():
       min_area = min(min_area,group.area)
#       quit()

#====================================================================
# Empennage weights
#====================================================================

    wght_ht   = 0.0#s_ht*(0.00395*s_ht**0.2*vdive-0.4885)
    wght_vt   = 4.9*min_area
    wght_tr = 0.0
    wght_at = 0.0

#====================================================================
# compound or propeller
#====================================================================

    if N_at > 0:
      wght_at = 0.0809484 * N_at * T_at**1.04771*(T_at/A_at)**(-0.07821)
    else:
      wght_at = 0.0
        
#====================================================================
# Apply tech factor and calculate total
#====================================================================

    wght_ht       = wght_ht*fac
    wght_vt       = wght_vt*fac
    wght_tr       = wght_tr*fac 
    wght_at       = wght_at*fac 
    
    total         = wght_ht + wght_vt + wght_tr + wght_at 

    empennage     = {     'htail': wght_ht*lb2kg,      'vtail': wght_vt*lb2kg, 
                     'tail_rotor': wght_tr*lb2kg, 'aux_thrust': wght_at*lb2kg,
                          'total':   total*lb2kg}
    return empennage 