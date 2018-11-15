from numpy import log10
from conversions import *
#====================================================================
# FUSELAGE GROUP
#
# airframe, pressurization, crashworthiness
# 
# ALL UNITS IN IMPERIAL
#====================================================================

f_lgloc      = 1.0#1.16    # 1.1627 landing gear on fus, 1.0 otherwise
f_lgret      = 1.0#1.1437  # retractable landing gear (goes into the fuselage)
f_ramp       = 1.0     # no cargo ramp, 1.27 otherwise
f_tfold      = 0.0     # tail of fold weight 0.05 assumed
f_wfold      = 0.0     # wing and rotor fold weight
f_mar        = 0.0     # marinization weight
f_press      = 0.0     # pressurization (fraction basic body weight)
f_cw         = 0.00    # crashworthiness weight (fraction fuselage weight)
nz           = 4.0     # design ultimate load factor

a            = -2.3979 #  
b            =  1.0    # empirical factors for 
c            = -0.0866 # wetted area (source?)
d            =  0.8099 # 
imodel       = 'afdd84'       # because why not
#imodel       = 'afdd82'
def fuselage_weight(vehicle_parameters):

#====================================================================
# unpack inputs
#====================================================================
   
   aircraftID = vehicle_parameters['aircraftID']      
   gtow       = vehicle_parameters['gtow']        # lbs
   nz         = 3.5                               # loa factor
   fac        = vehicle_parameters['tech_factors'].fuselage
   wing       = vehicle_parameters['wing']
   l_fus      = vehicle_parameters['l_fus']       # ft 
   
#====================================================================
# helicopter, compound or quad rotor /derivatives
#====================================================================

   S_body = 10**(c+d*log10(gtow))
   
#====================================================================
# AFDD84 Universal model
#====================================================================
   
   wght_basic  = 25.41*f_lgloc*f_lgret*f_ramp*(gtow*0.001)**0.4879 * \
                (gtow*nz*0.001)**0.2075 * (S_body **0.1676)  *      \
                (l_fus**0.1512) 

#====================================================================
# find additional weights for folding and other mechanisms/additions
#====================================================================

   wght_tfold = f_tfold * wght_basic                  # tail folding weight
   wght_wfold = f_wfold * (wght_basic + wght_tfold)   # wing folding weight
   wght_mar   = f_mar   * wght_basic                  # marinization weight
   wght_press = f_press * wght_basic                  # pressurization weight

   wght_cw    = f_cw * ( wght_basic + wght_tfold + wght_wfold + wght_mar + wght_press )

#====================================================================
# Apply tech factors 
#====================================================================

   wght_basic = wght_basic*fac
   wght_tfold = wght_tfold*fac
   wght_mar   = wght_mar  *fac
   wght_wfold = wght_wfold*fac
   wght_press = wght_press*fac
   wght_cw    = wght_cw   *fac

   total      =  wght_basic + wght_tfold + wght_wfold + wght_mar + wght_press + wght_cw    

#====================================================================
# store kg mass in a dictionary and return the dictionary
#====================================================================

   fuselage   = {       'basic': wght_basic*lb2kg, 'tail_folding': wght_tfold*lb2kg, 
                 'marinization': wght_mar  *lb2kg, 'wing_folding': wght_wfold*lb2kg,
               'pressurization': wght_press*lb2kg,  'crashworth' : wght_cw*lb2kg,
                      'total'  : total     *lb2kg}
   return fuselage