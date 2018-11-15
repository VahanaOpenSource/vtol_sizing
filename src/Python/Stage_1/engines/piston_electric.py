#====================================================================
# Turboshaft (Was initially in HYDRA. Has to be updated to RPTEM)
#====================================================================

from piston_engine  import piston_engine, getSFC
from transmissions  import electric_transmission 

class piston_electric:

   def __init__(self, emp_data):

      Engines           = emp_data.Engines
      Transmission      = emp_data.Transmission

      self.eta_xmsn     = Transmission.eta 
      self.fracInstall  = Engines.frac_install          # installation losses
      self.fracPowerMCP = Engines.fracPowerMCP          # we want continuous power draw to be 85% peak  

#====================================================================
# get the weight and fuel mass required
#
# Assume losses for installation and account for pressure change
# Obtain shaft horse power required based on temp and pres corrections
#
#====================================================================

   def getFuelWeight(self,powerReq,time,flightMode,theta,delta,Pmax):

#====================================================================
# idle or hover
#====================================================================

      P_engine    = powerReq/self.eta_xmsn

#====================================================================
# altitude and temperature correction
#====================================================================

      powerTOP      = P_engine#/((1.0-KT*(theta-1.0))*(1.0+KD*(delta-1.0)) ) 
      powerTOP     *= self.fracInstall/self.fracPowerMCP # installation losses

#      print powerReq/0.746,P_engine/0.746,powerTOP/0.746
#====================================================================
# get SFC of engine, and use it to calculate fuel mass (kg) and flow rate (kg/hr)
#====================================================================

      KT            = 1.0
      KD            = 1.0

#====================================================================
# installed power is at least 4 Hp (3 kW)
#====================================================================

      if powerTOP < 4*0.746:
         powerTOP   = 4*0.746           
  
      sfc           = getSFC(theta, delta, P_engine, Pmax, KT, KD)
      massFuel      = sfc * P_engine * time / 60.0
      fuelFlowRate  = P_engine * sfc

#====================================================================
# end of operations
#====================================================================
      
      return powerTOP, sfc, massFuel, fuelFlowRate

#====================================================================
#====================================================================
#
# obtain the weight of the powerplant and associated components
#
#====================================================================
#====================================================================

   def getWeight(self, inputs): 

#====================================================================
# unpack inputs
#====================================================================

      P_ins        = inputs['powerReq']
      neng         = inputs['nengine']
      tech_factor  = inputs['tfac_eng']
#====================================================================
# get engine weight based on curve fits
#====================================================================

      mass         = piston_engine(P_ins/neng) * neng

#====================================================================
# update dictionary
#====================================================================

      engine       = {'total': mass*tech_factor}
      return engine         # dictionary, units are SI [kg]
