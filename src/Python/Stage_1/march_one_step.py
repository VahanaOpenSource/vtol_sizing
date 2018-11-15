#====================================================================
# Routine to march one step of the design
# Mission segments are defined based on the mission.data input file
#====================================================================
import copy
import fea

class _march_one_step:

#====================================================================
# to pass empty weight to log file.. need a better way to do this
# maybe update currentvalues in the fortran module
#====================================================================

   def get_empty_weight(self):
      return self.massEmptyGroup

#====================================================================
# one iteration of the sizing loop
#====================================================================

   def march_one_step(self,itercount):


#====================================================================
# assign local pointers for easy addressing
#====================================================================

      use_bemt                      = self.all_dict['sizing']['use_bemt']
      mission                       = self.mission
      wing                          = self.wing
      prop                          = self.prop
      emp                           = self.emp_data.Tech_factors
      aircraft                      = self.all_dict['aircraft']
      techfac_fus                   = emp.Weight_scaling.fuselage

      massPayload                   = mission.payload
      sizing_mode                   = mission.sizing_mode
      massTakeoff                   = self.massTakeoff
      massFuelTotal                 = 0.0
      mission.totalenergyreqd       = 0.0
      mission.maxfuelflowrate       = 0.0
      mission.max_Peng              = 0.0

#====================================================================
# Remember previous take-off mass
#====================================================================

      old_takeoff_mass     = copy.copy(self.massTakeoff)
      old_payload_mass     = copy.copy(mission.payload)
      # print(old_takeoff_mass);x=input('ok?')
#====================================================================
# March through the different segments of the missiong
# The output of each segment is the mass of fuel used, power required
# for the particular segment
#====================================================================

      nseg                 = mission.nseg
      powerInstalled       = 0.0

      icontinue            = 1

#====================================================================
# Flags to decide whether to do both sizing & performance, or only
# performance (False)
#====================================================================

      update_h    = True               # perform hover sizing + performance calculations
      update_c    = True               # perform cruise sizing + performance calculations

#====================================================================
# loop over mission segments
#====================================================================

      for i_priority in range(nseg):
         iseg              = mission.sizing_order[i_priority]
         segment           = mission.segment[iseg]
         i                 = iseg+1 
         flightmode        = segment.flightmode

#====================================================================
# Idling before take-off
#====================================================================

         if (flightmode == 'idle'):
            iflag     = self.idlemodel(i,itercount)
            icontinue = min(iflag,icontinue)

#====================================================================
# hover mode
#====================================================================

         elif (flightmode == 'hover'):
            iflag     = self.hovermodel(i, update_h, use_bemt) # hover sizing+performance
            icontinue = min(iflag,icontinue)
            update_h  = False                               # switch off sizing for subsequent hover segments

#====================================================================
# cruise mode (climb incorporated here)
#====================================================================

         elif (flightmode == 'cruise'):

            iflag     = self.cruisemodel(i, update_c, use_bemt)
            icontinue = min(iflag, icontinue)
            update_c  = False

#====================================================================
# other types of flight conditions: not yet ready!
#====================================================================

         else:
             print('flight mode is ',flightmode)
             quit('UNKNOWN FLIGHT MODE DETECTED: PROGRAM TERMINATING')

#====================================================================
# update the fuel weight (from engine models) and get the new
# segments weights (engines.py)
#====================================================================

         self.updateFuelAndSegmentWeight(i)

#====================================================================
# evaluate the total quantity of fuel used and
# total power to be installed
#====================================================================

         massFuelTotal     = massFuelTotal + segment.fuel_mass
      self.mass_fuel       = massFuelTotal

#=============================================================================
# python based afdd model ( completely trust it for full-scale design)
#=============================================================================

      empty                = self.vehicle_empty_weight(icontinue)
      self.massEmptyGroup  = empty

#=============================================================================
# if ifea==1, update fuselage group with FEA analysis
#=============================================================================

      ifea = self.all_dict['sizing']['ifea']

      if(ifea==1 and icontinue==1):

         icontinue=self.updateAirframeWeight()
#         print 'calling fea'

         # ensure FEA has converged
         if (self.current_dict['fea_optcount'][-1] < fea.currentvalues.nmaxopt):

            mass                            = self.current_dict['fea_weight_hist'][-1]
            self.massEmptyGroup['fuselage'] = mass*techfac_fus

         if icontinue != 1:
            print ('impossible airframe design')
#         print ' fea done'

#=============================================================================
# update empty mass
#=============================================================================

      massEmpty = 0;
      for key in self.massEmptyGroup:

#=============================================================================
# check if the empty weight is a dictionary, which implies
# it also carries with it the breakdown of the components, get total
#=============================================================================
#
         if isinstance(self.massEmptyGroup[key],dict):
            massEmpty += self.massEmptyGroup[key]['total']
#            print(key,self.massEmptyGroup[key]['total'])
#=============================================================================
# not a dictionary, but only a floating point
#=============================================================================

         else:
            massEmpty += self.massEmptyGroup[key]
#            print(key,self.massEmptyGroup[key])
#      print(massEmpty,powerInstalled)

#=============================================================================
# add empty weight margin: 10% of empty weight
#=============================================================================

      self.massEmptyGroup['margin']    = 0.1*massEmpty
      massEmpty                        = massEmpty + self.massEmptyGroup['margin']

#=============================================================================
# Remember empty weight roll-up as class instance
#=============================================================================

      self.massempty = massEmpty

#=============================================================================
# iteration type 1: fixed payload, variable GTOW --> converge on this one
# update takeoff mass and calculate convergence error for take off mass
#=============================================================================

      if sizing_mode == 1:
         self.massTakeoff           = self.massempty + massFuelTotal + massPayload + self.mass_battery
         error                      = abs(self.massTakeoff - old_takeoff_mass)/(old_takeoff_mass+0.01)

#=============================================================================
# iteration type 2: fixed GTOW, variable payload
# calculate payload from (take off - empty - fuel reqd)
#=============================================================================

      elif sizing_mode == 2:
         massPayload                      = self.massTakeoff - self.massempty - massFuelTotal - self.mass_battery
         error                            = abs(massPayload - old_payload_mass)/(old_payload_mass+0.01)
         mission.payload                  = massPayload     
#         print('ok?',self.massTakeoff,self.massempty,massPayload)

#=============================================================================
# for negative payload, identify Chewbacca designs
#=============================================================================

         if massPayload < 1.e-3:
            icontinue               = 0
            self.errmsg             = 'negative payload'

#=============================================================================
# if design is converged but blade design is infeasible, trigger flag to
# discontinue iterations
#=============================================================================

      if(error <= 0.001):
         icontinue = 0
         print ('Invalid design: blade stresses too high')
         self.valid  = False
#====================================================================
# run only one iteration if sizing eVTOL with no payload drop during 
# mission, when in fixed take-off weight mode; don't need more!
#====================================================================

      if(sizing_mode == 2 and mission.same_weight and self.etype == 'electric_motor'):
            error = 0.0

#====================================================================
# Convergence check 
#====================================================================
   
      converged  = self.check_converged(error, itercount, icontinue)

      return converged

# ############################################################################
# END OF FILE
# ############################################################################