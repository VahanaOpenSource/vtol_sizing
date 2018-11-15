#=========================================================================   
# contains empty weight module for unconventional designs. the main
# AFDD empty weight model is still in the fortran. 
#=========================================================================   

from numpy import pi, exp
import os,sys

#=========================================================================   
#
# ELECTRIC MOTORS
#
#=========================================================================   


cable_density  = 1.3    # lb per foot
cooling        = 2      # lb per motor
eta_motor      = 0.97
eta_gearbox    = 0.98

def weight(aircraft, bfile):
    
    nrotors     = aircraft.nrotors

    tot_power   = aircraft.engine.ratedP * aircraft.engine.efficiency / (eta_motor*eta_gearbox)

    power       = tot_power / nrotors
    wght_motors = nrotors * 2.278*power**0.6563
    # wght_motors = nrotors * 0.5269*power**0.8983
    # wght_motors = nrotors * 2.2 * power**0.661
    # wght_motors = nrotors * 0.40199028*power**1.15382503
    # wght_motors = nrotors * 0.7907*power**0.8246
    # wght_motors = nrotors * 0.07529*power**1.226



    if(nrotors == 2):
        l_fuselage = 15.0
        wght_cables = 1.06*2*l_fuselage+nrotors*2*0.439*1.2*aircraft.rotors[0].R*1.2
    else:
        wght_cables = nrotors * 2 * aircraft.rotors[0].R * 1.414 * cable_density


    # motor cooling
    # wght_cool_m = nrotors * ( 19.1 * exp(0.007069 * power) -
    #                           7.735* exp(-0.01515 * power) )
    wght_cool_m = nrotors * (8.023    * exp(0.002331 * power) + 
                             0.005306 * exp(0.03762  * power))

    # generator cooling
    # wght_cool_g = ( 19.1 * exp(0.007069 * aircraft.engine.ratedP) -
    #                 7.735* exp(-0.01515 * aircraft.engine.ratedP) )

    wght_cool_g = ( 8.023    * exp(0.002331 * aircraft.engine.ratedP) + 
                    0.005306 * exp(0.03762  * aircraft.engine.ratedP) )

    # wght_batt   = tot_power * 0.10934 * 2.738
    wght_batt   = tot_power * 0.10934 * 3.0

    V_tip       = aircraft.rotors[0].V_tip0
    R           = aircraft.rotors[0].R
    omega_motor = 3000.0
    omega_rotor = V_tip/R*(60.0/(2.0*pi))

    wght_gb     = ( 95.7634 * 
                    power**0.78137 *
                    omega_motor**0.09889/omega_rotor**0.80686 )

    wght_dist   = 0.0
    if(aircraft.engine.__class__.__name__ == 'Hybrid_Sync'):
        wght_dist = 40.0

    if(not bfile is None):
        bfile.write("---- HYBRID SYSTEM ----\n")
        bfile.write("wght_motors %9.2f\n"%(wght_motors))
        bfile.write("wght_cables %9.2f\n"%(wght_cables)) 
        bfile.write("wght_cool_m %9.2f\n"%(wght_cool_m)) 
        bfile.write("wght_cool_g %9.2f\n"%(wght_cool_g)) 
        bfile.write("wght_batt   %9.2f\n"%(wght_batt)) 
        bfile.write("wght_gb     %9.2f\n"%(wght_gb))
        bfile.write("wght_dist   %9.2f\n"%(wght_dist))

    return ( wght_motors + wght_cables + 
             wght_cool_m + wght_cool_g + 
             wght_batt   + wght_gb + wght_dist)



#=========================================================================   
# END OF FILE
#=========================================================================   
