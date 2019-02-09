#====================================================================
# EMPENNAGE GROUP
#
# Horizontal and vertical tail
#
# ALL UNITS IN FPS
#====================================================================
from conversions import *
f_tr             = 1.0#1.6311 #if TR on vertical tail, otherwise 1. Mostly 1.63
A_ht             = 4.0  # horizontal tail aspect ratio
A_vt             = 4.0  # vertical tail aspect ratio

def helicopter_empennage(R, P_DSlimit, Omega, WMTO, tech_factor):

   """ 
   function to calculate weight of htail, vtail and tail rotor 
   unit for a conventional helicopter

   Inputs:
   1.   R         : main rotor radius, feet 
   2. P_DSlimit   : power limit for driveshaft, hp 
   3. Omega       : main rotor speed, rad/s
   4. WMTO        : max. take-off weight, lb
   5. tech_factor : a multiplier that scales the predicted empennage weights up or down

   Output:
   1. empennage : dictionary containing weights of  htail, vtail, tail rotor 
                  and the sum of all 3 components, in kg
   """

   R_tr        = R*0.2                               # 20% of main rotor radius

# helicopter (scale linearly with rotor radius based on UH-60A design)
# SMR or COAXIAL

   s_ht = 45.0 * (R / 26.83)**2.0           # horizontal tail area
   s_vt = 32.0 * (R / 26.83)**2.0           #   vertical tail area

   wght_ht = 0.7176        * s_ht**1.1881 * A_ht**0.3173
   wght_vt = 1.0460 * f_tr * s_vt**0.9441 * A_vt**0.5332

# add horizontal stabilizer-actuator weight
   wght_fc = 0.01735*(WMTO**0.64345)*(s_ht**0.40952)

# tail rotor weight
# for some reason, NDARC *divides* p_DSlimit/1.2 to use in eqn below
   wght_tr = 1.3778 * R_tr**0.0897 * (P_DSlimit/Omega/1.2)**0.8951

# apply tech factor 
   wght_ht       = wght_ht*tech_factor
   wght_vt       = wght_vt*tech_factor
   wght_tr       = wght_tr*tech_factor   
   wght_fc       = wght_fc*tech_factor

# get total, return dictionary of empennage weight breakdown and total
   total         = wght_ht + wght_vt + wght_tr + wght_fc

   empennage     = {     'htail': wght_ht*lb2kg,      'vtail': wght_vt*lb2kg, 
                    'tail_rotor': wght_tr*lb2kg,      'fcont': wght_fc*lb2kg,
                         'total':   total*lb2kg }
   return empennage 