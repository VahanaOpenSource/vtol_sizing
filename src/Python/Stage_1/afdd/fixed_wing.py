#====================================================================
# Ananth Sridharan, Jan 24th 2017
#
# AIRCRAFT WING GROUP: 29-1.2 SECTION IN NDARC THEORY MANUAL V1_11.PDF, pg 249
#
# AFDD 93 Model implemented by AS
# frequency-based sizing also implemented by AS
#====================================================================

f_LGloc     = 1.0                # 1.7247 if landing gear on wing, 1.0 otherwise
k_elect     = 0.5                # proportionality constants for anti-icing weight
k_rotor     = 0.5

#model       = 'statistical'
model       = 'target_freq'
fn          = 4.5                # wing spar freq, Hz
f           = 7.0                # motor mount freq, Hz
taper       = 0.8               
from conversions import *
from motor_mount_mass import motor_mount_mass as mmmass, spar_mass
import numpy
def fixed_wing_wt(vehicle_parameters, wing, rotor, motor):

   Wt          = vehicle_parameters['gtow']        # in lbs
   fac         = vehicle_parameters['tech_factors'].wing
   nr          = vehicle_parameters['nrotor']
   P           = vehicle_parameters['pwr_installed']*hp2kw
   l_fus       = vehicle_parameters['l_fus']
   b_fus       = vehicle_parameters['b_fus']*f2m
   clearance   = vehicle_parameters['clearance']      # spanwise clearance/radius from rotor to fuselage/another rotor
   red         = vehicle_parameters['wt_redund']

#====================================================================
# empirical parameters
#====================================================================

   tau_w    = 0.158                             # wing t/c (NOT TILTROTOR WING)
   Lamda    = 0.0                               # no sweep
   lamda    = 0.8                               # taper ratio
   nz       = 3.8                               # load factor 
   b_fold   = 0.0                               # fraction of span that folds (0 = not folding)

#====================================================================
# calculate weight of power wires and signal wires
#====================================================================
   
   c1          = 0.17            # kg/m of signal wire ==> based on ALpha
   c2          = 0.0057          # kg/kW/m of power wire, for (+) and (-)
   Nchannels   = 3               # number of parallel cables
   mass_sw     = c1*l_fus*Nchannels

#====================================================================
#battery to HVDB
#====================================================================

   mass_pw     = c2*(l_fus+b_fus*2)*P*red['wires'] 

#====================================================================
# expression for weight per wing, lbs
# then multiply by # wings, sum over groups
#====================================================================

   wing_wt     = 0.0
   actuator_wt = 0.0 
   area        = 0.0
   act_cost    = 0.0
   str_cost    = 0.0
   mounts_wt   = 0.0 

   # print(P);quit('ok?')
   # print(wing.ngroups)
   if wing.ngroups > 0:
      for i in range(wing.ngroups):

         group       = wing.groups[i]                 
         nwing       = group.nwings                   # number of wings
         Aw          = group.aspectratio              # aspect ratio 
         Sw          = group.area * m2f * m2f         # in sq.ft 
         fL          = group.lift_frac                # wing lift fraction for this group (0 to 1)
         W           = fL*Wt/nwing                    # thrust carried by each wing, lbs; this is both wings together!
            
#====================================================================
# first get rotor id, motor id and rotor/motor assembly masses
#====================================================================

         rgid        = group.rotor_group_id
         rg          = rotor.groups[rgid] 
         nrotors     = group.nrotors/(2*nwing)         # number of rotors along half span, each wing
         rotor_mass  = rg.mass_assembly

#====================================================================
# extract motor weight
#====================================================================

         mgid        = rg.motor_group_id 
         mg          = motor.groups[mgid]
         motor_mass  = mg.motor_mass

         P_rated     = numpy.amax(mg.p_ins)           # in kW, rated motor power

#====================================================================
# tip mass at end of rotor mount
# beam length = one rotor radius + 30% chord (where tilt axis, spar are located)
# tube radius = 15% of rotor radius
# target freq = 6.0 Hz
#====================================================================

         Mk          = motor_mass + rotor_mass 
         L_mount     = rg.radius + group.chord*0.3
         rtube       = 0.125*rg.radius*0.5

         if(model == 'statistical'):
            mount_mass  = 0.1*motor_mass 
         else:
            mount_mass  = mmmass(L_mount, Mk, rtube, f)

#====================================================================
# remember mass of all motor mounts, store in lbs
#====================================================================

         mounts_wt  = mounts_wt + 2.2*mount_mass*nrotors*(2*nwing)         

#====================================================================
# Target frequency model: we have masses of rotors, motors and 
# motor mounts; first place these rotors along the span 
# simultaneously assemble lumped mass vector
#====================================================================

         y        = []                    # spanwise location of motor mounts
         Mk       = []
         T        = []

         offset   = 0.5*b_fus + clearance*rg.radius
         Pdsum    = 0.0
         Lsum     = 0.0
         for irotor in range(int(nrotors)):
            y_rotor  = offset+rg.radius

#====================================================================
# add up rated power x distance x 2 (left+right wings) 
# for signal wire, find total wire length in wings
#====================================================================

            l_wire   = (y_rotor+L_mount)
            
            Pdsum    = Pdsum + l_wire*2*P_rated
            Lsum     = Lsum  + l_wire*2

#====================================================================
# calculate power wire and signal wire weights going out to each motor 
# remember that m_pw and m_sw are in kgs, and for one side of a wing
#====================================================================

#This effect is very small.. can be ignored.
#            m_pw     = c2*l_wire*P_rated*red['wires']
#            m_sw     = c1*l_wire*Nchannels
#            y.append(y_rotor*0.5)            # keep lumped mass @ half-span
#            Mk.append(m_pw+m_sw)             # lumped mass = signal+power wire

#====================================================================
# motor + rotor after wire, because its outboard
#====================================================================

            y.append(y_rotor)       
            Mk.append(mount_mass+ motor_mass + rotor_mass)
            T.append(numpy.amax(group.rotor_thrust[0]))
            offset   = offset + rg.radius*(2.0 + clearance)

#====================================================================
# add signal wire length of half-span x 2 (left+right) for control surface
#====================================================================
         
         Lsum        = Lsum + group.span

#====================================================================
# calculate spar mass
#====================================================================
            
         tskin       = 15e-4              # 3 layers, each 0.5mm; based on Alpha
         lbyc        = 2.1
         rho         = 1650.0             # mean density, kg/cu.m (fold in foam too)
         MbyA        = tskin*lbyc*rho

         c           = group.chord        # mean wing chord in meters
         rbar        = tau_w*c*0.5        # mean tube radius , meters   

         rroot       = 2*rbar/(1+taper)   # root spar radius
         L           = group.span*0.5     # beam length

#====================================================================
# weight to match structural frequency of wing spar
#====================================================================

         if(model == 'target_freq'):
            M_spar      = spar_mass(rroot, taper, L, MbyA, Mk, y, T, wing, fn, tau_w)
            wt          = (MbyA*group.area + M_spar*2)*2.2

#====================================================================
# wing weight, statistical method: AFDD model
#====================================================================

         else:
            wt       = 5.66411*f_LGloc* (W*0.001/cos(Lamda))**0.847 * (nz**0.3958) * (Sw**0.21754)        \
                       *sqrt(Aw)* ((1.0+taper)/tau_w)**0.09359 #* (1.0 - b_fold)**(-0.14356)

#====================================================================
# add up contributions from duplicate wings in this group
# proceed with book-keeping operations
#====================================================================

         wt          = wt * nwing
         wing_wt     = wing_wt + wt
         area        = area + Sw*nwing             # in sq.m
         group.wt    = wt*lb2kg/nwing

#====================================================================
# Add signal wire and power wire weights here; ideal place because we
# loop over wings in this function, so piggyback. 
#====================================================================

         dsignal     = c1*group.span*nwing*Nchannels
         dsignal2    = c1*Lsum*nwing*Nchannels

#====================================================================
# version 1: approximation
# version 2: exact sum
#====================================================================

         dpower      = c2*nwing*P*group.span*red['wires']
         dpower2     = c2*nwing*Pdsum*red['wires']

         mass_pw     = mass_pw + dpower2
         mass_sw     = mass_sw + dsignal2

         group.wire  = dsignal + dpower # wires in each wing, kg

#====================================================================
# for control surfaces, take total area and total weight; scaling law adapted by AS
# Note: right now redundancy is there in cost but not weight!
#====================================================================

      actuator_wt    = actuator_wt + 0.01735*(Wt**0.6435)*(area**0.40952)

#====================================================================
# weight effect of redundany actuators
#====================================================================

   actuator_wt       = actuator_wt * red['wing_flap']

#====================================================================
# tilt actuators: 7.5% of mass of wing
#====================================================================
      
   tilt_wt        = 0.075*wing_wt*red['tilt_actuator']

#====================================================================
# apply tech factor, take total weight of wing+actuators on wing
#====================================================================

   actuator_wt       = actuator_wt * fac
   wing_wt           = wing_wt     * fac
   tilt_wt           = tilt_wt     * fac 
   total             = actuator_wt + wing_wt + tilt_wt + mounts_wt

#====================================================================
#convert to kg, return dict for wing structure, actuators and wires
#====================================================================

   wt       = {'structure': wing_wt*lb2kg,'actuators':actuator_wt*lb2kg, 
               'total': total*lb2kg, 'tilters': tilt_wt*lb2kg,
               'mounts': mounts_wt*lb2kg}

   wires    = {'signal_wire': mass_sw, 'power_wire': mass_pw, \
                     'total': mass_sw + mass_pw} 

   return wt, wires