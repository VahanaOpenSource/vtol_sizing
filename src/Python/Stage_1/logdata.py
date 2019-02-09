#====================================================================
# save the properties of a given configuration
#====================================================================

import os,sys
import pickle
from conversions import kg2lb
#====================================================================
# can be accessed by hydraInterface
#====================================================================

class _logdata:

#====================================================================
# log file is yaml. easy to read into a post-processing script
# to evaluate the performance
#====================================================================

   def writelogdata(self,filename):

      massEmptyGroup    = self.massEmptyGroup
      powerplant        = self.powerplant 

# ===================================================================
# write class to file
# ===================================================================

      f2       = filename.rstrip('.yaml') + '.txt'
      with open(f2,'wb') as f:
         pickle.dump(self, f,protocol=2)

# ===================================================================
# myello
# ===================================================================

      std      = self.constants
      sdict    = self.all_dict['sizing']
      adict    = self.all_dict['aircraft']
      mission  = self.mission
      f=open(filename,'w')

      f.write('# ---------------------------------------------------\n')
      f.write('# Mission and Vehicle Log\n')
      f.write('# ---------------------------------------------------\n')
      f.write('\n')

# ===================================================================
# vehicle details
# ===================================================================

      f.write('vehicle:\n')
      f.write("   {:25s} {:d} \n"            .format('aircraftID:',adict['aircraftID']) )
      f.write("   {:25s} {:.2f} # [kg]\n"    .format('take_off_mass:',self.massTakeoff))
      f.write("   {:25s} {:.1f} # [kW]\n"    .format('power_installed:',self.p_ins))     # mechanical power available at rotor
      f.write("   {:25s} {:.2f} # [m] \n"    .format('D-value:',self.footprint) )
      f.write("   {:25s} {:.2f} # [Lift to Drag ratio]\n"   .format('Cruise_LbyD:',self.LbyD) )

#==============================
# max SPL and max sound distance
#==============================

#      f.write('\nAcoustics:\n')
#      f.write("   {:20s} {:5.1f} #[ft]\n"    .format('hover_altitude:',250))         # hover altitude
#      f.write("   {:20s} {:5.1f} #[ft]\n"    .format('dist_centerline:', self.maxsnd_d))
#      f.write("   {:20s} {:5.2f} #[dB]\n"    .format('maximum_SPL:',self.max_spl))   # max spl

# ===================================================================
# rotor related data
# ===================================================================

      rotor   =   self.rotor

      if rotor.nrotors > 0:
         f.write('\nRotor:\n')
         for i in range(rotor.ngroups):
            group    = rotor.groups[i]
            f.write("   set{:d}: \n"       .format(i))
            f.write("      {:20s} {:d} \n"             .format('nrotor:',group.nrotors))
            f.write("      {:20s} {:d} \n"             .format('nblade:',group.nblade) )
            f.write("      {:20s} {:.2f} # [lb/ft2]\n" .format('disk_loading:',group.diskloading/std.grav * std.kg2lb / (std.m2f**2)) )
            f.write("      {:20s} {:.2f} \n"           .format('aspect_ratio:', group.aspectratio))
            f.write("      {:20s} {:.2f} # [m]\n"      .format('radius:',group.radius))
            # print('radius is ')
            f.write("      {:20s} {:.3f} # [m]\n"      .format('chord:',group.chord ))
            f.write("      {:20s} {:.1f} # [m/s]\n"    .format('tip_speed:',group.tipspeed))
            f.write("      {:20s} {:.3f} \n"           .format('hover_FM:', group.fm))
            f.write("      {:20s} {:.3f} \n"           .format('cruise_rpm_ratio:', group.RPM_ratio))
            # f.write("      {:20s} {:.3f} \n"           .format('eta_xmsn:', self.engine.eta_xmsn))
            f.write("      {:20s} {:.5f} \n"           .format('solidity:', group.solidity))
            f.write("      {:20s} {:.3f} \n"           .format('cd0:', group.cd0))
            f.write("      {:20s} {:.3f} \n"           .format('ipf:', group.ipf))
            f.write("      {:20s} {:.3f}\n"            .format('hvr_dwld:',group.hvr_dwld))
            f.write("      {:20s} {:.3f}\n"            .format('hvr_ct_sigma:',group.ctsigma))
         # f.write("      {:20s} {:.3f} #[N-m]\n"     .format('hvr_torque:',rotor.torque))
            f.write("      {:20s} {:.3f}\n"            .format('prop_eta:',self.prop.eta) )
            f.write('\n')

# ===================================================================
# wing related data
# ===================================================================

      wing=self.wing

      if wing.ngroups > 0:
         f.write('Wings:\n')
         for i in range(wing.ngroups):
            group    = wing.groups[i]
            f.write("   set{:d}: \n"       .format(i))
            f.write("      {:20s} {:d} \n"       .format('nwing:',group.nwings) )
            f.write("      {:20s} {:.3f} # [m]\n" .format('span:',group.span))
            f.write("      {:20s} {:.3f} # [m]\n" .format('chord:',group.chord))
            f.write("      {:20s} {:.3f} # [kg, each] \n".format('structure_wt:', group.wt))
#            f.write("      {:20s} {:.3f} # [kg, each] \n".format('wires_wt:', group.wire))
            f.write("      {:20s} {:.3f} \n"       .format('aspect_ratio:', group.aspectratio))
            f.write("      {:20s} {:.3f} \n"       .format('oswald:', group.oswald))
            f.write("      {:20s} {:.3f} \n"       .format('cd0:', group.cd0))
            f.write("      {:20s} {:.3f} \n"       .format('cl:', group.cl))
            f.write("      {:20s} {:.3f} \n"       .format('lift_fraction:', group.lift_frac))
            f.write("      {:20s} {:3d} \n"        .format('rotors_per_wing:', int(group.nrotors/group.nwings)))
            f.write("      {:20s} {:3d} \n"        .format('rotor_group_id:',int(group.rotor_group_id) ))
         f.write('\n')

# ===================================================================
# propeller related data
# ===================================================================

#      if self.aircraftID != 1:
      # f.write('prop:\n')
      # f.write("   {:20s} [{:d}] \n"       .format('nprop:',adict['npropeller']) )
      # f.write("   {:20s} [{:.3f}] \n"       .format(  'eta:',adict['effpropeller']) )

#==============================
# acquisition cost breakdown
#==============================
      
      c     = self.costs 
      if(bool(c)):
         f.write('Costs:\n')
         f.write("   {:20s} [{:.6f}] # [Millions of USD]\n"    .format('Frame_acquisition:',c.acquisition/1e6))     # acquisition cost
         f.write("   {:20s} \n"    .format('acquisition_cost_breakdown:'))                              # acquisition cost breakdown
         for key in sorted(c.acq_breakdown.keys()):
               temp = c.acq_breakdown[key]
               f.write("      {:16s}: [{:15.3f}         , {:6.3f}] # [USD, % acquisition cost]\n".format(key,temp,temp/c.acquisition*100))

#==============================
# fixed cost breakdown
#==============================

         f.write("   {:20s} [{:.3f}] # [USD]\n"    .format('Fixed_operating_costs:',c.fixed_costs))     # acquisition cost
         f.write("   {:20s} \n"    .format('fixed_cost_breakdown:'))    # cost breakdown
         for key in sorted(c.fix_breakdown.keys()):
               temp = c.fix_breakdown[key]
               f.write("      {:16s}: [{:15.3f}         , {:6.3f}] # [USD, % fixed cost]\n".format(key,temp,temp/c.fixed_costs*100))

#==============================
# variable cost breakdown
#==============================

         f.write("   {:20s} [{:.3f}] # [USD/hr]\n"    .format('Variable_operating_costs:',c.variable_costs))     # variable costs, $/hr
         f.write("   {:20s} \n"    .format('variable_cost_breakdown:'))    # cost breakdown
         for key in sorted(c.var_breakdown.keys()):
               temp = c.var_breakdown[key]
               f.write("      {:16s}: [{:15.3f}         , {:6.3f}] # [USD/hr, % variable cost]\n".format(key,temp,temp/c.variable_costs*100))

#==============================
# time comparison vs. taxi
#==============================

         f.write("   {:20s} [{:15.3f}] # [minutes    ]\n"    .format('UAM_time:',c.UAM_time))           # variable costs, $/hr
         f.write("   {:20s} [{:15.3f}] # [minutes    ]\n"    .format('Taxi_time:',c.Taxi_time))         # variable costs, $/hr
         f.write("   {:20s} [{:15.3f}] # [USD        ]\n"    .format('UAM_trip_cost:',c.UAM_cost))         # variable costs, $/hr
         f.write("   {:20s} [{:15.3f}] # [USD        ]\n"    .format('Taxi_trip_cost:',c.taxi_cost))         # variable costs, $/hr
         f.write("   {:20s} [{:15.3f}] # [$/min saved]\n"    .format('Time_value:',c.dollar_per_min))   # variable costs, $/hr

#==================================================================================================
# Weight breakdown
#==================================================================================================

      unit_str    = 'kg'
      unit_mlt    = 1.0

#      unit_str    = 'lb'
#      unit_mlt    = kg2lb
      f.write('\nWeights:\n')
      f.write('   empty_weight:\n')
      f.write("      {:29s} [{:6.1f}         ,{:6.1f}] # [{:2s}, %GTOW]\n" .format('total:',self.massempty*unit_mlt, self.massempty/self.massTakeoff*100, unit_str))
      for key in sorted(massEmptyGroup):
         mtemp = massEmptyGroup[key]

# look through dictionary if it exists
         if isinstance(mtemp,dict):

            grp_tot = mtemp['total']
            if grp_tot > 1.e-3:
               f.write("      {:16s}:\n".format(key) )
               f.write("         {:25s}: [{:6.1f}         ,{:6.1f}] # [{:2s}, % GTOW]\n".format('total',grp_tot*unit_mlt,grp_tot/self.massTakeoff*100.0,unit_str))
               for key2 in sorted(mtemp):
                  mtemp2 = mtemp[key2]
                  if(mtemp2 > 1.e-3):
                     if key2 != 'total':
                        f.write("         {:25s}: [{:6.1f}] # [{:2s}]\n".format(key2,mtemp2*unit_mlt, unit_str))

         elif mtemp > 1.e-5: # check minimum mass
            f.write("      {:28s}: [{:6.1f}         ,{:6.1f}] # [{:2s}, %GTOW]\n".format(key,mtemp*unit_mlt,mtemp/self.massTakeoff*100.0,unit_str))
            
#==============================
# fuel & payload
#==============================
#      print self.battery_wt
      m_f   = self.powerplant.mass_fuel
      m_b   = self.powerplant.mass_battery
      if m_f > 0.0:
         f.write("   {:32s} [{:6.1f}         ,{:6.1f}] # [{:2s}]\n" .format('fuel:',    m_f*unit_mlt, m_f/self.massTakeoff*100,unit_str))
      if m_b > 0.0:
         f.write("   {:32s} [{:6.1f}         ,{:6.1f}] # [{:2s}]\n" .format('battery:', m_b*unit_mlt, m_b/self.massTakeoff*100,unit_str))

      m_p   = self.mission.payload
      f.write("   {:32s} [{:6.1f}         ,{:6.1f}] # [{:2s}]\n" .format('payload:', m_p*unit_mlt, m_p/self.massTakeoff*100, unit_str))

#==============================
# Volumes
#==============================

      f.write('\nVolumes:\n')
      f.write("   {:25s} {:5.2f} # [m]\n"     .format('fuselage_width:',self.geometry.fuselage_width) )
      if('pax_count' in adict):
         f.write("   {:25s} {:5d} # [people]\n"    .format('passenger_count:',adict['pax_count']))

      if(self.transmission.groups[0].type == 'electric'):
         f.write("   {:25s} {:5.3f} # [cu.m]\n"     .format('battery_volume:',powerplant.battery_vol) )

#==============================
# flat-plate area breakdown
#==============================

      f.write('\nParasitic_drag:\n')
      f_vehicle      = self.f_plate
      f.write("   {:20s} {:.3f} # [sq.m] \n"  .format('flat_plate_area:', f_vehicle))
      if adict['fdrag'] == 0:

         f.write("   {:20s} \n"  .format('flat_plate_breakdown:'))
         all_f     = adict['f_breakdown']
         for key in sorted(all_f.keys()):
            this_f = all_f[key]
            if(this_f > 1.e-3 and key != 'total'):
               f.write("      {:16s}: [{:5.1f}] # [%]\n".format(key,this_f/(f_vehicle)*100))

# ===================================================================
# mission data
# ===================================================================

      m         = self.mission
      nsegments = m.nseg
      f.write('\nMission:\n')

#==============================
# total energy
#==============================

      if(self.powerplant.groups[0].type == 'battery'):
         f.write("   {:20s} {:.2f} # [kW-hr]\n" .format('total_energy:',powerplant.ops_energy))
      f.write("   {:20s} {:d} \n".format('nsegments:',nsegments) )

      f.write("   {:20s} [".format('flight_mode:' ))
      for n in range(nsegments):
         f.write("'{:8s}'".format(m.segment[n].flightmode))
         if n < nsegments-1:
            f.write(',')
      f.write('] \n')

      f.write("   {:20s} [".format('start_alt:' ))
      for n in range(nsegments):
         f.write('{:10.0f} '.format(m.segment[n].startalt))
         if n < nsegments-1:
            f.write(',')
      f.write('] #[m]\n')

      f.write("   {:20s} [".format('end_alt:' ))
      for n in range(nsegments):
         f.write('{:10.0f} '.format(m.segment[n].endalt))
         if n < nsegments-1:
            f.write(',')
      f.write('] #[m]\n')

      f.write("   {:20s} [".format('delta_temp_ISA:' ))
      for n in range(nsegments):
         f.write('{:10.2f} '.format(m.segment[n].deltatempisa))
         if n < nsegments-1:
            f.write(',')
      f.write('] #[C]\n')

      f.write("   {:20s} [".format('rate_of_climb:' ))
      for n in range(nsegments):
         f.write('{:10.0f} '.format(m.segment[n].rateofclimb * std.m2f))
         if n < nsegments-1:
            f.write(',')
      f.write('] #[ft/min]\n')

      f.write("   {:20s} [".format('cruise_speed:' ))
      for n in range(nsegments):
         f.write('{:10.2f} '.format(m.segment[n].cruisespeed))
         if n < nsegments-1:
            f.write(',')
      f.write('] #[knots]\n')

      f.write("   {:20s} [".format('time:' ))
      for n in range(nsegments):
         f.write('{:10.2f} '.format(m.segment[n].time))
         if n < nsegments-1:
            f.write(',')
      f.write('] #[min]\n')

#==========================================================
# power required to spin rotor 
#==========================================================

      f.write("   {:20s} [".format('rotor_power_reqd:' ))
      for n in range(nsegments):
         f.write('{:10.2f} '.format(self.rotor.groups[0].p_req[n]))
         if n < nsegments-1:
            f.write(',')
      f.write('] #[kW]\n')

#==========================================================
# power draw from battery
#==========================================================

      f.write("   {:20s} [".format('battery_power_draw:' ))
      for n in range(nsegments):
         f.write('{:10.2f} '.format(m.segment[n].p_eng*0.001))
         if n < nsegments-1:
            f.write(',')
      f.write('] #[kW]\n')

#=================== 
# C-rating
#=================== 

      if self.transmission.groups[0].type == 'electric':
         f.write("   {:20s} [".format('C-rating:' ))
         for n in range(nsegments):
            f.write('{:10.2f} '.format(m.segment[n].c_rating))
            if n < nsegments-1:
               f.write(',')
         f.write('] #[1/hr]\n')

#=================== 
# cell temperature
#=================== 

         f.write("   {:20s} [".format('cell temperature:' ))
         for n in range(nsegments):
            f.write('{:10.2f} '.format(m.segment[n].cell_temp))
            if n < nsegments-1:
               f.write(',')
         f.write('] #[deg C]\n')

#=================== energy reqd
#      f.write("   {:20s} [".format('energy_segment:' ))
#      for n in range(nsegments):
#         f.write('{:6.2f} '.format(m.segment[n].energy * std.kw2hp))
#         if n < nsegments-1:
#            f.write(',')
#      f.write('] #[hp-hr]\n')

#=================== density

      f.write("   {:20s} [".format('density:' ))
      for n in range(nsegments):
         f.write('{:10.3f} '.format(m.segment[n].rho ))
         if n < nsegments-1:
            f.write(',')
      f.write('] #[kg/cu.m]\n')

#=================== segment weight

      f.write("   {:20s} [".format('segment_mass:' ))
      f.write('{:10.1f} '.format(self.massTakeoff))
      # print(self.massTakeoff)
      f.write(',')
      for n in range(nsegments-1):
         # print(m.segment[n].mass)
         f.write('{:10.1f} '.format(m.segment[n].mass ))
         if n < nsegments-2:
            f.write(',')
      f.write('] #[kg]\n')

#=================== segment weight

      f.write("   {:20s} [".format('segment_type:' ))
      for n in range(nsegments):
         f.write("'{:9s}'".format(m.segment[n].type ))
         if n < nsegments-1:
            f.write(',')
      f.write(']\n')

#=================== segment distance

      f.write("   {:20s} [".format('segment_distance:' ))
      for n in range(nsegments):
         f.write("{:11.2f}".format(m.segment[n].distance ))
         if n < nsegments-1:
            f.write(',')
      f.write(']\n')

#=================== 

#conversion factor from lb/hp-hr to kg/kw-hr
      fac        = 0.45359/0.7457

      # if(self.transmission.groups[0].type != 'electric'):
      #    f.write("   {:20s} [".format('sfc_segment:' ))
      #    for n in range(nsegments):
      #       f.write('{:6.3f} '.format(m.segment[n].sfc / fac))         # in lb/hp-hr
      #       if n < nsegments-1:
      #          f.write(',')
      #    f.write('] #[lb/hp-hr]\n')

#===================
#battery assumptions
#=================== 
      
      if('Battery' in self.all_dict['empirical']):
         Pack        = self.all_dict['empirical']['Battery']['Pack']
         f.write('\nBattery:\n')
         f.write("   {:20s} {:.2f} \n" .format('rated_capacity:',powerplant.rated_energy) )
         f.write("   {:20s} {:.3f} \n" .format('state_of_health:',Pack['SOH']) )
         f.write("   {:20s} {:.3f} \n" .format('depth_discharge:',Pack['DOD_min']) )

#=================== 
# close the file
#=================== 

      f.close()
      