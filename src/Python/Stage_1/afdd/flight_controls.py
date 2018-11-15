#==================================================================
# FLIGHT CONTROL SYSTEM and HYDRAULICS
#
# Cockpit controls, automatic flight control system, system cont.
# Accounting for rotor controls and conversion controls
#
# Hydraulics for fixed wing and rotary wing flight controls,
# converstion flight controls and equipment
#==================================================================

from conversions import *

#==================================================================
# flight controls & hydraulics
#==================================================================

f_RWnb  = 1.25 # fraction rotary wing non-boosted weight
f_RWhyd = 0.4  # fraction rotary wing hydraulic weight
f_RWred = 3.0  # redundancy factor
f_CVnb  = 0.1  # fraction conversion non-boosted weight (to mass_CVmb)
f_CVmb  = 0.02 # fraction conversion boosted mech (to mass_to)
f_mbsv  = 1.3029  # ballistic survivability
f_bsv   = 1.117  # ballistic survivability
f_CVhyd = 0.1  # fraction conversion hyd. weight (to boost mech)

#==================================================================
# begin routine
#==================================================================

def flightctrl_weight(vehicle_parameters):
    
    nrotor = vehicle_parameters['nrotor'] 
    nblade = vehicle_parameters['nblade'] 
    chord  = vehicle_parameters['chord']  
    V_tip  = vehicle_parameters['vtip']   
    WMTO   = vehicle_parameters['gtow']   
    aid    = vehicle_parameters['aircraftID']
    nprop  = vehicle_parameters['nprop']
    nwing  = vehicle_parameters['nwing']
    R      = vehicle_parameters['radius']      # ft 
    fac    = vehicle_parameters['tech_factors'].flight_control 

#==================================================================
# fixed-wing controls + boost mechanisms
#==================================================================
    
    if nwing > 0:                              # including all fixed wing surfaces
        w       = 0.91*(WMTO**0.6)
    else:                                      # only Htail actuator
        s_ht    = 45.0 *(R / 26.83)**2         # horizontal tail area
        w       = 0.01735*(WMTO**0.64345)*(s_ht**0.40952)

#==================================================================
#rotary-wing flight control mechanisms
#==================================================================

    w_fc = ( 0.2873 * f_mbsv * (nrotor * nblade)**0.6257 *
             chord**1.3286 * (0.01*V_tip)**2.1129 *
             f_RWred**0.8942 )

#boosted flight control weight for rotary-wing
    wght_RWb  = ( 0.02324 * f_bsv * (nrotor*nblade)**1.0042 *
                  nrotor**0.1155 * chord**2.2296 * 
                  (0.01*V_tip)**3.1877 )

#==================================================================
#add aux thrusters also to determine controls wt
#==================================================================
    
    if nprop > 0:
        V_tip    = 600.0             # prop tip speed, ft/s 
        chord    = chord*0.2         # assume 20% of MR chord
        nrotor   = nprop             #  
        nblade   = 4                 # 4 blades/prop

        w_fc     =  w_fc + 0.2873 * f_mbsv * (nrotor * nblade)**0.6257 *   \
                  chord**1.3286 * (0.01*V_tip)**2.1129 *                   \
                  f_RWred**0.8942 

        wght_RWb = wght_RWb  + 0.02324 * f_bsv * (nrotor*nblade)**1.0042 * \
                   nrotor**0.1155 * chord**2.2296 *                        \
                   (0.01*V_tip)**3.1877 

#==================================================================
#boost mechanism weight for rotary wing
#==================================================================

    wght_RWmb = (1-f_RWhyd)*w_fc

#==================================================================
#rotary-wing non boosted flight control weight
# non-boosted control weight is a parameter value x flight control 
# weight  
#==================================================================

    wght_RWnb = f_RWnb*(1-f_RWhyd)*w_fc
   
    if aid == 2:

#==================================================================
#conversion control boost mechanism weights: fraction of GTOW 
#==================================================================

        wght_CVmb = f_CVmb*WMTO

#==================================================================
#non-boosted conversion control weight
#==================================================================

        wght_CVnb = f_CVnb*wght_CVmb
    else:
        wght_CVmb = 0.0
        wght_CVnb = 0.0

#==================================================================
#weight of avionics, rotary-wing non boosted controls, boosted controls and boost mechanisms
#==================================================================

    wght_avionics = 0.0

    all_fltcon    =  wght_RWmb + wght_RWb + wght_RWnb + wght_CVmb + wght_CVnb + wght_avionics

#add fixed wing actuator weights 

    if aid == 5:
      fx_wing_ctrl  = 0.0 
    else:
      fx_wing_ctrl  = w

    rt_wing_ctrl  = all_fltcon

#==========================================================================
# hydraulics: for rotary-wing and conversion mechanisms
# based on fraction of flight control weight and fraction of conversion 
# boost mechanism weight
#==========================================================================

    wght_RWhyd      = f_RWhyd*w_fc
    wght_CVhyd      = f_CVhyd*wght_CVmb
    
    wght_hydraulics = wght_RWhyd + wght_CVhyd

#==========================================================================
# Apply tech factor, calculate total weight
#==========================================================================

    wght_hydraulics = wght_hydraulics*fac 
    rt_wing_ctrl    = rt_wing_ctrl   *fac 
    fx_wing_ctrl    = fx_wing_ctrl   *fac 
    
    total           = wght_hydraulics + rt_wing_ctrl + fx_wing_ctrl

#==========================================================================
# use hydraulics and flight control system models for specific aircraft 
# 1,2,3,7
#==========================================================================

    flt_control = {'hydraulics' : wght_hydraulics* lb2kg, 
                   'rw_flt_ctrl': all_fltcon     * lb2kg,
                   'fw_flt_ctrl': fx_wing_ctrl   * lb2kg,
                   'total'      : total          * lb2kg }

#==========================================================================
# return dictionary with breakdown
#==========================================================================

    return flt_control