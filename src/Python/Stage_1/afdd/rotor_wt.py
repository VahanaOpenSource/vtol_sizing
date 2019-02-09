#====================================================================
# ROTOR GROUP
#
# Blades, hub and hinge, spinner, blade fold structure
# 
# ALL UNITS IN IMPERIAL (FPS)
#====================================================================

from conversions import *
f_fold          = 0.0

#====================================================================
# begin routine
#====================================================================

def rotor_weight(vehicle_parameters):

#====================================================================
# unpack inputs
#====================================================================

   aircraftID = vehicle_parameters['aircraftID']
   nblade     = vehicle_parameters['nblade']   
   nrotor     = vehicle_parameters['nrotor']   
   R          = vehicle_parameters['radius']   
   chord      = vehicle_parameters['chord']    
   V_tip      = vehicle_parameters['vtip']    
   nu_blade   = vehicle_parameters['nu_blade'] 
   factor     = vehicle_parameters['tech_factor']

#====================================================================
# extra parameters for tilt-rotor
#====================================================================

   imodel     = 'afdd82'
   if(aircraftID == 1 or aircraftID == 3 or aircraftID == 5 or aircraftID==6 or aircraftID == 7):
     f_tilt   = 1.0
     dia_spin = 0.0
   elif(aircraftID == 2):       # tilt-rotor: has f_tilt factor in blade weight eqn
     f_tilt   = 1.1794
     dia_spin = 0.2 * R 
     imodel   = 'afdd00'
   else:
     print ('rotor.py: Unknown aircraft configuration')
     sys.exit(1)

#====================================================================
# blade weight models from statistical equations
#====================================================================

   nu_hub     = nu_blade
   
#====================================================================
# AFDD00 model: blade weight and hub weight
#====================================================================

   if imodel == 'afdd00':
     wght_blade = ( 0.0024419      * f_tilt            *
                    nrotor         * nblade**0.53479   *
                    R**1.74231     * chord**0.77291    *
                    V_tip**0.87562 * nu_blade**2.51048 )

# hub weight: for coaxial, use curve fit given in Johnson/Moodie/Yeo paper
# "Design and Performance of Lift-Offset Rotorcraft for Short-Haul Missions"
# includes the shaft weight too.. still heavier than XH-59A
     if aircraftID == 3:
       wght_hub   = ( 0.0061182 * nrotor * nblade**0.20373 *
                      R**0.60406      * V_tip**0.52803  *
                      nu_hub**1.00218 * (wght_blade/nrotor)**0.87127 )
     else:
       wght_hub   = ( 0.1837 * nrotor * nblade**0.16383 *
                      R**0.19937      * V_tip**0.06171  *
                      nu_hub**0.46203 * (wght_blade/nrotor)**1.02958 )

#====================================================================
# AFDD 82 model: blade weight and hub weight
#====================================================================

   if imodel == 'afdd82':
     wght_blade = ( 0.02606*nrotor  * nblade**0.6592   *
                    R**1.3371       * chord**0.9959    *
                    V_tip**0.6682   * nu_blade**2.5279 )

     wght_hub   = ( 0.003722*nrotor * nblade**0.2807 *
                    R**1.5377       * V_tip**0.4290  *
                    nu_hub**2.1414  * (wght_blade/nrotor)**0.5505 )
	  
#====================================================================
# spinner and folding system weights in lb
#====================================================================

   wght_spin    = 7.386*nrotor*dia_spin**2
   wght_fold    = f_fold*wght_blade

#====================================================================
# Apply tech factor for rotor group
#====================================================================

   wght_blade   = wght_blade*factor 
   wght_hub     = wght_hub  *factor 
   wght_spin    = wght_spin *factor
   wght_fold    = wght_fold *factor 
   
#====================================================================
# mass of components in kg
#====================================================================

#   total = wght_blade + wght_hub + wght_spin + wght_fold 
#   total roll-up is performed outside in empty weights

   weight = {'blades'  : wght_blade *lb2kg,
             'hub'     : wght_hub   *lb2kg,
             'spinner' : wght_spin  *lb2kg,
             'folding' : wght_fold  *lb2kg}

   mass_kg_per_unit   = (wght_blade + wght_hub + wght_spin + wght_fold)/nrotor*lb2kg

#====================================================================
# return the dictionary
#====================================================================

   return weight, mass_kg_per_unit