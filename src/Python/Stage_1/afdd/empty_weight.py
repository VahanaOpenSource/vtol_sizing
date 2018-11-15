#====================================================================
#
# Empty weight model
#
#====================================================================

import os,sys,numpy

from alighting       import alighting_weight    # Fraction of GTOW
from conversions     import *                   # conversion factor
from engine          import engine_accessories  # AFDD82 
from drive           import drivesys_weight     # AFDD83
from empennage       import empennage_weight    # AFDD82
from flight_controls import flightctrl_weight   # AFDD82
from fuel_system     import fuel_system         # AFDD82
from fixed_wing      import fixed_wing_wt       # AFDD93 
#from fuselage        import fuselage_weight     # AFDD84
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
        aid          = aircraft['aircraftID']
        wing         = self.wing
        prop         = self.prop 
        p_ins        = self.p_ins     # installed power, in kW
        emp_data     = self.emp_data
        rotor        = self.rotor 
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
# calculate total fuel reqd for mission
#====================================================================

        fuel            = self.mass_fuel

#====================================================================
# aircraft ID
#====================================================================
  
        aircraftID      = self.all_dict['aircraft']['aircraftID']

#====================================================================
# weight for powerplant
#====================================================================

        powerplant      = self.getEngineWeight()

        l_fus           = 6.0/0.3048              # in feet

#====================================================================
# calculate costs of battery and motor
#====================================================================

        if(self.etype == 'electric_motor'):
          self.mass_battery  = self.engine.m_batt

#====================================================================
# create a dictionary with the required inputs (vehicle_parmeters)
#====================================================================
  
        vparams = {'aircraftID'    : aircraftID,
                   'flow_rate'     : self.maxfuelflowrate*kg2lb,    # lb/hr
                   'fuel'          : fuel*kg2lb,
                   'gtow'          : self.massTakeoff * kg2lb,
                   'wing'          : wing,
                   'nprop'         : prop.num,
                   'T_prop'        : prop.thrust * N2lb,
                   'nrotor'        : wing.groups[0].nrotors,
                   'nu_blade'      : rotor.groups[0].nu_beta, 
                   'pwr_installed' : p_ins * kw2hp,
                   'vdive'         : 1.3*Vmax,
                   'wingarea'      : wing.area * m2f * m2f,
                   'nwing'         : wing.nwings      ,
#                   'wing_chord'    : wing.chord * m2f     ,
                   # 'wing_AR'       : wing.aspectratio     ,
                   'wing_tip_wt'   : self.wingtip_mass*kg2lb,
                   'tech_factors'  : tech_factors,
                   'b_fus'         : emp_data.Geometry.fuselage_width*m2f,
                   'l_fus'         : l_fus, 
                   'wing'          : self.wing, 
                   'wt_redund'     : self.wt_redund}  

        Volume           = 0.0

#=============================================================================
# physics-based model for rotor weight prediction
#=============================================================================

#=============================================================================
# Loop over all rotor groups
#=============================================================================

        rotor_weight  = {}
        ngroups       = rotor.ngroups
        for i in range(ngroups):
            group     = rotor.groups[i]

#=============================================================================
# loop over wing groups that contain this rotor and find the sizing thrust, 
# sizing torque; these values are the max vals across all segments too..
#=============================================================================

            sizing_T  = 0.0
            sizing_Q  = 0.0
            for wgid in group.wing_group_ids:
                T_try    = numpy.amax(wing.groups[wgid].rotor_thrust)
                Q_try    = numpy.amax(wing.groups[wgid].rotor_torque)

                sizing_T = max(sizing_T, T_try)
                sizing_Q = max(sizing_Q, Q_try)

#=============================================================================
# Set sizing thrust and torque for the group, then perform sizing
#=============================================================================

            group.sizing_thrust   = sizing_T 
            group.sizing_torque   = sizing_Q 
            
            rotor_wt              = group.bladewt(tech_factors.rotor)

#=============================================================================
# remember weight breakdown for the group in a dictionary
#=============================================================================

            for key,value in rotor_wt.items():
               k2     = 'group'+str(i)+key 
               rotor_weight[k2] = value 

#=============================================================================
# repeat operations for de-icing weight
#=============================================================================

            anti_icing = group.icing_weight(tech_factors.anti_icing)

            for key,value in anti_icing.items():
               k2     = 'group'+str(i)+key 
               rotor_weight[k2] = value 

#====================================================================
#add up all the entries in the rotor groups
#====================================================================

        total         = 0.0
        for key,value in rotor_weight.items():
            total     = total + value 
        rotor_weight['total']   = total

#====================================================================
# Airframe, empennage and landing gear weights, and volume 
#====================================================================

        Volume            = 0.0
#        fuselage          = fuselage_weight(vparams)
        fuselage          = fuselage_v2(vparams)
        Volume            = Volume + fuselage['total']/1400.0 

        alighting         = alighting_weight(vparams)

#=============================================================================
# Flight controls: accounted for in "avionics" (common equipment weight)
#=============================================================================

        flt_control     = 0.0

#====================================================================
# engine accessories only applicable for turboshaft!
#====================================================================

        if (self.etype == 'turbo_electric' or self.etype == 'turboshaft' or self.etype == 'stuttgart'):
            engine_accs     = engine_accessories(vparams)
            Volume          = Volume + engine_accs['total']/1400.0
        else:
            engine_accs     = 0.0

#=============================================================================
# Mechanical drive systems not present for electric transmissions
# Drivetrain is counted as part of the powerplant
# NOTE: WIRE WEIGHTS ARE NOT ACCOUNTED FOR YET EXPLICITLY 
#=============================================================================

        if (self.etype=='electric_motor'):
            fuel_sys        = 0.0
            drive_system    = {'total': 0.0}              # done with battery!

#=============================================================================
# fuel handling, transmission weight for turboshaft/piston engines: AFDD models
#=============================================================================

        elif (self.etype=='turboshaft' or self.etype=='piston' or self.etype=='diesel' or   \
              self.etype=='stuttgart'):
            fuel_sys        = fuel_system(vparams)
            Volume          = Volume + fuel_sys['total']/2700.0
            drive_system    = drivesys_weight(vparams)
            Volume          = Volume + drive_system['total']/2700.0
            quit('mechanical transmissions disabled!')
#=============================================================================
# approximation by VT: 0.5 lb of transmission wt / unit hp of installed power
#=============================================================================

#            if (aircraftID == 5 or aircraftID == 6 or aircraftID == 7):
#               drive_system = 0.5*vparams['pwr_installed']*0.746*lb2kg   

#=============================================================================
# for turbo/piston with electric transmissions, use AFDD fuel handling
# there's not any mechanical drive system
#=============================================================================

        elif (self.etype=='turbo_electric' or self.etype == 'piston_electric'):
            fuel_sys        = fuel_system(vparams)            
            Volume          = Volume + fuel_sys['total']/1400.0
            drive_system    = electric_transmission(p_ins)
            drive_system    = drive_system*tech_factors.drive_system
            
            Volume          = Volume + drive_system/8960.0

#====================================================================
# Error trap
#====================================================================

        else:
            print( 'I DONT KNOW WHAT POWERPLANT IS BEING USED')
            print( 'Valid options are: ')
            print( 'electric_motor, al_air_motor')
            print( 'turboshaft, piston')
            print( 'turbo_electric, piston_electric')
            print( 'User input is ', self.etype )
            quit('CRITICAL ERROR: PROGRAM TERMINATING')

#====================================================================
# Wing weight model: use AFDD model
#====================================================================

        wing_wt, wires  = fixed_wing_wt(vparams)
        empennage       = empennage_weight(vparams, wing_wt)

#====================================================================
# collect empty weights into dictionary of dictionaries; all kgs
#====================================================================

        massEmptyWeight = {'wing'           : wing_wt,
                           'rotor'          : rotor_weight,
                           'empennage'      : empennage,
                           'fuselage'       : fuselage,
                           'alighting_gear' : alighting,
                           'fuel_system'    : fuel_sys,
                           'drive_system'   : drive_system,
                           'flight_control' : flt_control,
                           # 'anti_icing'     : anti_icing,     # moved within rotors+wings
                           'engine_accs'    : engine_accs,
                           'powerplant'     : powerplant, 
                           'common_equip'   : aircraft['common_equipment']*self.wt_redund['avionics'],
                           'avionics'       : aircraft['avionics'],
                           'emergency_sys'  : emergency_sys(vparams),
                           'wires'          : wires}
#        print(aircraft['common_equipment']);quit()
#        for k,v in sorted(massEmptyWeight.items()):
#          print(k,v)
#        quit()

#====================================================================
# Tail-sitter: no extra structure reqd. for landing gear
#====================================================================

        if (aircraftID==5):
           massEmptyWeight['alighting_gear'] = 0.0 

#====================================================================
# When using our own FEA model, dont use some elements of AFDD model
# flight control system, empennage, alighting gear, engine accessories
# => common equipment
# wing load-carrying members are done with FEA-based weight model 
# need to account for filler materials
#====================================================================

        if self.all_dict['sizing']['ifea'] == 1:

            massEmptyWeight = {
                               'wing'           : wing_wt, 
                               'rotor'          : rotor_wt,
                               'fuselage'       : fuselage,
                               'fuel_system'    : fuel_sys, 
                               'drive_system'   : drive_system, 
                               'engine_accs'    : engine_accs,
                               'common_equip'   : aircraft['mass_common_equip'],
                               'powerplant'     : powerplant}

#====================================================================
# Note: wing is only the skin + filler; "spar" is done by physics-based
# model for airframe. 
# wing weight is linearly scaled by the density of material
# obtained from 8lb QBT vehicle 
# Each section weighs 87 g (1.018 m with an AR = 4)
#
# Update by AS
# assume that wing CS has t/c =  12%, and skin thickness = 0.1% chord
# construction material is carbon fiber for skin
#====================================================================

            wing_wt       = 0.0
            for i in range(wing.ngroups):
              group       = wing.groups[i]
              chord       = group.chord
              span        = group.span
              skin_t                  =  0.001   * chord
              tmin                    =  0.0005
              if (skin_t < tmin):
                  skin_t = tmin
              skin_p                  =      2   * chord
              skin_A                  =  skin_t  * skin_p
              skin_vol                =  skin_A  * span           # cu.m
              skin_rho                = 1400.0                    # kg/cu.m
              skin_wt                 = skin_rho * skin_vol       # for one wing

              CS_t                    =     0.10 * chord
              CS_A                    = CS_t     * chord 
              CS_vol                  = CS_A     * span 
              CS_rho                  = 32.0                      # half of honeycomb, kg/cu.m
              CS_wt                   = CS_vol   * CS_rho         # for one wing 

#====================================================================
#weight of filler and skin only is called "wing"
#====================================================================

              group_wt                = group.nwings * (CS_wt + skin_wt)
              wing_wt                 = wing_wt + group_wt 

#            print 'weight breakdown so far',massEmptyWeight
#            print 'disabled wing weight thingummy',wing_wt
            massEmptyWeight['wing']   = wing_wt

#====================================================================
# return dictionary of empty weights and total volume
#====================================================================

#        Volume      = Volume * 2.5          # add 150% empty space (cooling, clearance)
  
        return massEmptyWeight

#====================================================================
# END OF FILE
#====================================================================