#==================================================================
# ANTI-ICING GROUP
#
# Deicing for rotors, wings, engines
# 
#==================================================================

from numpy import pi
from conversions import *
#==================================================================
# fixed parameters
#==================================================================

k_elect = 0.5 # assume 5# of the weight of the system
k_rotor = 0.5
k_wing  = 0.5
k_air   = 0.06          # de-icing for air induction

def icing_weight(vehicle_parameters):

#==================================================================
# unpack inputs
#==================================================================

    nrotor       = vehicle_parameters['nrotor']
    R            = vehicle_parameters['radius']     # MR radius, ft
    chord        = vehicle_parameters['chord']      # rotor chord
    wght_eng     = vehicle_parameters['engine_wt']  # weight of all engines
    fac          = vehicle_parameters['tech_factors'].anti_icing
    wing         = vehicle_parameters['wing']

#==================================================================
# derived quantities
#==================================================================

    Ablades      = nrotor  *R * chord      # blade plan-form area, sq.ft 

#==================================================================
# find weight of blade de-icing wires: some const * blade area
#==================================================================

    wght_DIelect = k_elect * Ablades

#==================================================================
# heating element weights scaled 
#==================================================================

    w0           = 0.0
    for i in range(wing.ngroups):
        w0       = w0 + k_wing*wing.groups[i].span*wing.groups[i].nwings 

    wght_DIsys   = ( k_rotor *  Ablades  + 
                     w0                  + #Sum of k_wing *wing_span 
                     k_air   * wght_eng )

#==================================================================
# apply tech factor, calculate total weight and return dictionary
#==================================================================

    wght_DIelect = wght_DIelect*fac 
    wght_DIsys   = wght_DIsys*fac
    
    total_wt     = wght_DIelect + wght_DIsys 

    deicing      = {'deicing_blades': wght_DIelect*lb2kg, 
                     'deicing_bl_eq': wght_DIsys  *lb2kg, 
                     'total': total_wt    *lb2kg}

    return deicing