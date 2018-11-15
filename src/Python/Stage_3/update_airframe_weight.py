#====================================================================
# set the beam cross sectional area and find the new airfram weight
#====================================================================

import fea
import os,sys

from create_fea_config import createAirframeConfig


class _update_airframe_weight:

   def updateAirframeWeight(self):
 
      mission       = self.mission 
      rotor         = self.rotor
      wing          = self.wing

      nz            = 3.5  # load factor
      #
      radius        = rotor.radius
      wingspan      = wing.groups[0].span

      Th            = nz * mission.segment[0].thrust      # hover  thrust 
      Tc            = nz * mission.segment[1].torque      # cruise thrust

      Qh            = nz * mission.segment[0].torque      # torque in hover 
      Qc            = nz * mission.segment[1].torque      # torque in cruise
      winglift      = nz * wing.lift
      wingdrag      = nz * wing.groups[0].drag 

#====================================================================
# set value for maximum allowable deflection
#====================================================================

      fea.currentvalues.maxallowabledeflection = 0.08   #

# ===================================================================
# transfer information to FEA data structure
# ===================================================================

      iodata         = fea.fea_types.iodata  # pointer to Fortran data structure
      iodata.new_r   = radius                # new radius 
      iodata.new_b   = wingspan              # new wing span (tip-to-tip distance)
      iodata.th      = Th                    # thrust PER ROTOR in hover 
      iodata.qh      = Qh                    # torque PER ROTOR in hover 
      iodata.qc      = Qc                    # torque PER ROTOR in cruise
      iodata.tc      = Tc                    # thrust PER ROTOR in cruise 
      iodata.l       = winglift              # lift PER WING
      iodata.d       = wingdrag              # +ve along +Z, so put -ve sign 

#====================================================================
# scale coordinates of FEA model in fortran based on new span, radius
#====================================================================

      fea.scale_coordinates()

#====================================================================
# set the forcing values and get new span if constraints are hit
#====================================================================

      fea.set_loads()

      for i in range(wing.ngroups):
            group       = wing.groups[i]
            group.span  = iodata.new_b

# ===================================================================
# update wing lift coefficient;; using modified span and same aspect
# ratio, update Cl target for next iteration 
# ===================================================================

            c           = group.span/group.aspectratio 
            area        = group.span*c                 # area of each wing 
            group.cl    = group.cl*group.area/area      

# ===================================================================
# optimize the structure (fortran routine)
# ===================================================================

      icontinue            = 1
      fea.optimizestructure(icontinue)
      
# ===================================================================
# append FEA history into self.current_dict
# ===================================================================
      
#      print self.current_dict.keys()
      if ('fea_optcount' not in self.current_dict):
         self.current_dict['fea_optcount'] = []

      if ('fea_weight_hist' not in self.current_dict):
         self.current_dict['fea_weight_hist'] = []

      optcount     = fea.currentvalues.optcount
      minoptweight = fea.currentvalues.minoptweight[0:optcount]

      self.current_dict['fea_optcount'].append(optcount)
      for j in range(optcount):
         self.current_dict['fea_weight_hist'].append(minoptweight[j])

      return icontinue

# ===================================================================
# END OF FILE
# ===================================================================
