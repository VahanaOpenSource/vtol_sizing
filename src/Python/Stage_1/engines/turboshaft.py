#====================================================================
# Turboshaft (Was initially in HYDRA. Has to be updated to RPTEM)
#====================================================================

from turbo_engine   import turbo_engine, getSFC
from transmissions  import mechanical_transmission 
class turboshaft:

   def __init__(self, Engines):

      self.KD_mrp       = Engines['KD_mrp']
      self.KD_irp       = Engines['KD_irp']
      self.KD_mcp       = Engines['KD_mcp']
#
      self.KT_mrp       = Engines['KT_mrp']
      self.KT_irp       = Engines['KT_irp']
      self.KT_mcp       = Engines['KT_mcp']
#
      self.fracInstall  = Engines['frac_install']          # installation losses
      self.fracPowerMCP = Engines['fracPowerMCP']          # we want continuous power draw to be 85% peak  

#====================================================================
# get the weight and fuel mass required
#
# Assume losses for installation and account for pressure change
# Obtain shaft horse power required based on temp and pres corrections
#
#====================================================================

   def getFuelWeight(self,P_engine,time,flightMode,theta,delta,Pmax):

# idle or hover condition
      if flightMode == 'idle' or flightMode == 'hover':
         KT = self.KT_irp
         KD = self.KD_irp

# cruise or climb
      elif flightMode == 'climb' or flightMode == 'cruise':
         KT = self.KT_mcp
         KD = self.KD_mcp
      else:         
         print('you gave flightMode input as ',flightMode)
         quit('Unknown flight condition... must be 0,1,2,3: turboshaft.py')

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

      # print('inputs to getsfc',P_engine/1e3, Pmax/1e3)

      sfc           = getSFC(P_engine*0.001, Pmax*0.001, delta, theta)

#      print('after sfc calc',P_engine/1e3, Pmax/1e3, delta, theta)
#      print(sfc,P_engine/1e3);quit()
      massFuel      = sfc * P_engine*0.001 * time / 60.0
      fuelFlowRate  = P_engine * sfc*0.001
#      print 'take-off power = ',powerTOP
#      print 'required power = ',powerReq
      # print(powerTOP/1e3,sfc,massFuel,fuelFlowRate);quit()
#====================================================================
# end of operations
#====================================================================
      
      return powerTOP, sfc, massFuel, fuelFlowRate

#====================================================================
# obtain the weight of the powerplant and associated components
#====================================================================

   def getWeight(self, nengine, powerReq, tech_factor): 

      """
      this function calculates the weight of a turboshaft engine, given 
      the following input parameters
      1. nengine:       number of engines
      2. powerReq:      total installed power for this engine in kilowatts (all nengine units)
      3. tech_factor:   a multiplier to the engine weight to make it lighter or heavier
      """

      P_ins          = powerReq          # in kw, per unit
      
      mass           = turbo_engine(P_ins) *nengine

      total          = {'total': mass*tech_factor}
      return total         # dictionary, units are SI [kg]
