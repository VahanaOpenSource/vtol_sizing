#====================================================================
# Turboshaft (Was initially in HYDRA. Has to be updated to RPTEM)
#====================================================================

from turbo_engine   import turbo_engine, getSFC
from transmissions  import mechanical_transmission 
class turboshaft:

   def __init__(self, emp_data):

      Engines           = emp_data.Engines
      Transmission      = emp_data.Transmission

      self.eta_xmsn     = Transmission.eta 
      
      self.KD_mrp       = Engines.KD_mrp
      self.KD_irp       = Engines.KD_irp
      self.KD_mcp       = Engines.KD_mcp 
#
      self.KT_mrp       = Engines.KT_mrp
      self.KT_irp       = Engines.KT_irp
      self.KT_mcp       = Engines.KT_mcp
#
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
# engine power output
#====================================================================

      P_engine    = powerReq/self.eta_xmsn

#====================================================================
# idle or hover condition
#====================================================================

      if flightMode == 0 or flightMode == 1:
         KT = self.KT_irp
         KD = self.KD_irp

#====================================================================
# cruise or climb
#====================================================================

      elif flightMode == 2 or flightMode == 3:
         KT = self.KT_mcp
         KD = self.KD_mcp
      else:
         
          'Unknown... engines.py'
          sys.exit(1)

#====================================================================
# altitude and temperature correction
#====================================================================

      powerTOP      = P_engine/((1.0-KT*(theta-1.0))*(1.0+KD*(delta-1.0)) ) 
      powerTOP     *= self.fracInstall/self.fracPowerMCP # installation losses

#====================================================================
# get SFC of engine, and use it to calculate fuel mass (kg)
# and flow rate (kg/hr); SFC is in kg/kw-hr, P_engine is in kW, 
# time is in min, t/60 = hrs
#====================================================================

      sfc           = getSFC(P_engine, Pmax, delta, theta)
      massFuel      = sfc * P_engine * time / 60.0
      fuelFlowRate  = P_engine * sfc

#      print 'take-off power = ',powerTOP
#      print 'required power = ',powerReq

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

      neng           = inputs['nengine']
      P_ins          = inputs['powerReq']/neng          # in kw
      tech_factor    = inputs['tfac_eng']
      
#====================================================================
# get engine weight based on curve fits
#====================================================================

      mass           = turbo_engine(P_ins) *neng 

      total          = {'total': mass*tech_factor}
      return total         # dictionary, units are SI [kg]
