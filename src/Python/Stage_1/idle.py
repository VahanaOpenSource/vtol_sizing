#====================================================================
# Hover model of the aircraft
# Performs basic sizing values and obtains the power required
#====================================================================
import math
import os

#====================================================================
# Function body (segID comes in with BASE 1)
#====================================================================

class _idlemodel:
  def idlemodel(self, segID, itercount):

   icontinue      = 1

#====================================================================
# assign pointer shortcuts
#====================================================================

   mission        = self.mission
   rotor          = self.rotor
   pi             = self.constants.pi

   iseg           = segID - 1       # python counting

   segment        = mission.segment[iseg]

#====================================================================
# if first iteration (assume idle power is a fraction of total power)
# else assume that the idle power is equal to the profile power
#====================================================================
   
   if(itercount==1):
      segment.p_rotor      = 0.1*self.p_ins              # watts
      segment.p_req        = segment.p_rotor * 1.e-3     # kw
      segment.torque       = 0.01e0
      segment.thrust       = 0.01e0
      
   else:
      cd0                  = self.rotor.cd0
      cpProfile            = rotor.solidity * cd0 / 8.0

      Cptotal              = cpProfile *rotor.num

      power                = Cptotal* segment.rho*rotor.area*(0.85*rotor.tipspeed)**3

      mission.segment[iseg].p_rotor      = power                 # in watts
      mission.segment[iseg].p_req        = power*1.e-3           # in kW

      mission.segment[iseg].torque       = 0.e0
      mission.segment[iseg].thrust       = 0.e0             # N-m

   return icontinue
