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

from skin_mass import skin_mass
def empennage_weight(wing, tech_factors):

    R_at        = 0.5                              # assume % rotor radius
    A_at        = pi*R_at*R_at
    fac         = tech_factors.empennage 

#====================================================================
# Vtail has same area as htail, 4.9 kg/sq.m
#====================================================================

    max_chord    = 0.0
    for ig,group in wing.groups.items():
       max_chord = max(max_chord,group.chord)

#====================================================================
# Vtail properties
#====================================================================

    chord        = max_chord*0.75
    rho_skin     = 1650.0   
    span         = 2.5*chord
    wght_vt      = skin_mass(chord, rho_skin, span)*2.2

#====================================================================
# Empennage weights
#====================================================================

    wght_ht      = 0.0
    wght_tr      = 0.0
    wght_at      = 0.0
        
#====================================================================
# Apply tech factor and calculate total
#====================================================================

    wght_ht       = wght_ht*fac
    wght_vt       = wght_vt*fac
    wght_tr       = wght_tr*fac 
    
    total         = wght_ht + wght_vt + wght_tr

    empennage     = {     'htail': wght_ht*lb2kg,      'vtail': wght_vt*lb2kg, 
                     'tail_rotor': wght_tr*lb2kg,      'total':   total*lb2kg}
    return empennage 
