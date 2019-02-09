#====================================================================
#
# Empty weight model
#
#====================================================================

import os,sys,numpy

from alighting       import alighting_weight    # Fraction of GTOW
from conversions     import *                   # conversion factor
#from engine          import engine_accessories  # AFDD82 
#from drive           import drivesys_weight     # AFDD83
from empennage       import empennage_weight    # AFDD82
from helicopter_emp  import helicopter_empennage
from flight_controls import flightctrl_weight   # AFDD82
from fuel_system     import fuel_system         # AFDD82
from fuselage        import fuselage_weight     # AFDD84 model
from fuselage_v2     import fuselage_v2         # GB/ZL model
from icing           import icing_weight        # refer NDARC v1_11
from rotor_wt        import rotor_weight        # AFDD82
from wing            import tiltrotor_wing_wt
from transmissions   import mechanical_transmission, electric_transmission
from emergency_sys   import emergency_sys
class _empty_weight:

   def vehicle_empty_weight(self, icontinue):
      
#====================================================================
# create links for rotor and wing
#====================================================================

        aircraft     = self.all_dict['aircraft']
        aircraftID   = aircraft['aircraftID']
        wing         = self.wing
        rotor        = self.rotor 
        prop         = self.prop 
        powerplant   = self.powerplant
        transmission = self.transmission
        fuselage     = self.fuselage
        geometry     = self.geometry 
        p_ins        = self.p_ins     # installed power, in kW
        emp_data     = self.emp_data
        tech_factors = emp_data.Tech_factors.Weight_scaling

#====================================================================
# find maximum torque and airspeed (dive speed for tiltrotor)
#====================================================================

        mission     = self.mission 
        nseg        = mission.nseg 
        Vmax        = 0.0 
        for i in range(nseg):
            V1      = mission.segment[i].cruisespeed
            Vmax = V1 if (Vmax < V1) else Vmax
  
#====================================================================
# create a dictionary with the required inputs (vehicle_parmeters)
#====================================================================

        l_fus           = geometry.fuselage_length
        Volume          = 0.0

#====================================================================
# powerplant + intake anti-icing weights
#====================================================================
      
        pp_wts          = powerplant.weight_rollup(tech_factors, mission)

#====================================================================
# fuel system weight
#====================================================================
        
        GTOW          = self.massTakeoff*kg2lb
        fuel_system   = powerplant.fuel_system(GTOW, tech_factors.fuel_system)

#=============================================================================
# Get rotor system weights
#=============================================================================

        rotor_weight  = rotor.weight_rollup(wing, fuselage, tech_factors, aircraftID)

#=============================================================================
# loop over transmission, connected rotors and find blade weights if reqd
#=============================================================================

        xmsn_weight, P, omega   = transmission.weight_rollup(tech_factors, rotor, powerplant)

#====================================================================
# place the rotors around the vehicle and estimate wire & cable weights
#====================================================================

        if not(aircraftID == 1 or rotor.nrotors == 1):
          wires, p_ins          = geometry.rotor_placement(rotor, wing, fuselage, transmission)
          powerplant.p_ins      = p_ins
          wires['power_cables'] = wires['power_cables']*self.wt_redund['wires']
          wires['total']        = wires['signal_wires'] + wires['power_cables']
        else:
          wires                 = {'total': 0.0}

#====================================================================
# Wing weight roll-up
#====================================================================

        wing_weight       = wing.weight_rollup(GTOW, tech_factors, self.wt_redund, rotor)

#====================================================================
# Airframe, empennage and landing gear weights, and volume 
#====================================================================

        Volume            = 0.0

        alighting         = alighting_weight(self.massTakeoff, tech_factors.landing_gear)

# for conventional helicopter, fuselage length = 1.3 radius
# use AFDD84 model here, as well as statistics-based interpolation for empennage
        if(aircraftID == 1):
          R_mr            = rotor.groups[0].radius
          l_fus           = R_mr*1.3
          fuselage_wt     = fuselage_weight(GTOW, tech_factors.fuselage, l_fus*m2f)
          empennage_wt    = helicopter_empennage(R_mr*m2f, P*kw2hp, omega, GTOW, tech_factors.empennage)

# for eVTOL, load paths are different, take loads based approach
        else:
          empennage_wt    = empennage_weight(wing, tech_factors)
          fuselage_wt     = fuselage_v2(GTOW, tech_factors.fuselage, wing, geometry.fuselage_length)

        Volume            = Volume + fuselage_wt['total']/1400.0 

#=============================================================================
# approximation by VT: 0.5 lb of transmission wt / unit hp of installed power
#=============================================================================

#            if (aircraftID == 5 or aircraftID == 6 or aircraftID == 7):
#               drive_system = 0.5*vparams['pwr_installed']*0.746*lb2kg   

#====================================================================
# collect empty weights into dictionary of dictionaries; all kgs
#====================================================================

        massEmptyWeight = {'wing'           : wing_weight,
                           'rotor'          : rotor_weight,
                           'empennage'      : empennage_wt,
                           'fuselage'       : fuselage_wt,
                           'alighting_gear' : alighting,
                           'fuel_system'    : fuel_system,
                           'drive_system'   : xmsn_weight,
                           'powerplant'     : pp_wts, 
                           'common_equip'   : aircraft['common_equipment']*self.wt_redund['avionics'],
                           'avionics'       : aircraft['avionics'],
                           'emergency_sys'  : emergency_sys(self.massTakeoff, tech_factors.emergency_sys),
                           'wires'          : wires}

#====================================================================
# When using our own FEA model, dont use some elements of AFDD model
# flight control system, empennage, alighting gear, engine accessories
# => common equipment
# wing load-carrying members are done with FEA-based weight model 
# need to account for filler materials
#====================================================================

        if self.all_dict['sizing']['ifea'] == 1:

            quit('not ready: FEA based sizing for airframe members')
            massEmptyWeight = {
                               'wing'           : wing_wt, 
                               'rotor'          : rotor_wt,
                               'fuselage'       : fuselage,
                               'fuel_system'    : fuel_sys, 
                               'drive_system'   : drive_system, 
                               'engine_accs'    : engine_accs,
                               'common_equip'   : aircraft['mass_common_equip'],
                               'powerplant'     : powerplant}

# return dictionary of empty weights and total volume
        return massEmptyWeight

#====================================================================
# END OF FILE
#====================================================================