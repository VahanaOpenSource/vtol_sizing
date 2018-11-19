import sys, numpy 
from numpy import sin, cos, log10, log2, sqrt, pi
from scipy.special import jv as besselj
sys.path.insert(0,'../Stage_0/')

from conversions import *
k_elect = 0.5 # proportionality constants
k_rotor = 0.5


from blade_wt_modelv2 import blade_wt_modelv2

#====================================================================
# Motor class that remembers parameters for sizing
#====================================================================

class motors:

   def __init__(self, data, nseg):
      self.ngroups      = 0
      self.groups       = {}

      ngrp              = 0
#loop over motor groups, find design parameters
      for key in sorted(data):
         self.groups[ngrp]  = motor_group(data[key], key, nseg)
         ngrp               = ngrp + 1

      self.ngroups      = ngrp
      self.nmotors      = 0       # total number of motors
      return None

#====================================================================
# Motor group details
#====================================================================

class motor_group:

#====================================================================
# function to initialize all segments
#====================================================================

   def __init__(self, data, key, nseg):

#=======================================================================
#geometric parameters
#=======================================================================

      self.cruise_efficiency   = data['cruise_efficiency']
      self.hover_efficiency  = data['hover_efficiency']
      self.key                = key 
      self.rotor_group_ids    = []      # which rotor groups are used to size this motor?
      self.p_ins              = numpy.zeros(nseg)
      self.nmotors            = 0

#====================================================================
# Rotor class that remembers parameters for sizing
#====================================================================

class rotors:

   def __init__(self, data, nseg):
      self.ngroups      = 0
      self.nrotors      = 0         # total rotor count in aircraft
      ngrp              = 0
      self.groups       = {}
      self.ntilt        = 0 
      self.ncruise      = 0
      self.nlift        = 0
      self.Atilt        = 0.0 
      self.Acruise      = 0.0 
      self.Alift        = 0.0 

#====================================================================
#loop over rotor groups, find design parameters
#====================================================================

      for key in sorted(data):
         self.groups[ngrp]  = rotor_group(data[key], key, nseg)
         ngrp               = ngrp + 1

      self.ngroups      = ngrp

      return None
      
#====================================================================
# Python function to calculate blade drag in forward flight 
# for edgewise rotors
# Inputs
#        Vcruise (m/s)
#        VtipMax (m/s): maximum rotor section speed at advancing blade tip
#        rho     (kg/cu.m): air density
#
# Outputs
#        Blade_drag (Newtons)
#====================================================================

   def blade_drag(self, Vcruise, VtipMax, rho):

#====================================================================
# now loop over all rotor groups, and calculate rotor drag
#====================================================================

     Blade_drag           = 0.0

     for i in range(self.ngroups):
        group             = self.groups[i]
        NR                = group.nrotors


        Vtip              = group.tipspeed*group.RPM_ratio

#====================================================================
# edgewise rotor: get average profile drag coefficient of sections
#====================================================================

        if group.type == 'edgewise':
           if group.tipspeed + Vcruise > VtipMax:
              Vtip              = VtipMax - Vcruise
              group.RPM_ratio   = Vtip/group.tipspeed
           else:
              Vtip              = group.tipspeed

#====================================================================
# find average section drag coefficient at different advance ratios
#====================================================================

           cd0                  = group.cd0
           muCruise             = Vcruise/Vtip
           if muCruise > 1.0:
              cd0 = 1.5*cd0
           elif muCruise > 0.3:
              cd0 = cd0*(1.0 + (muCruise-0.3)/0.7e0)
           else:
              cd0 = cd0

#====================================================================
# find blade drag in wind direction
#====================================================================

           D_blade        = group.solidity*cd0/8.0*(3.1*muCruise)
           D_blade        = D_blade * rho * group.area * Vtip*Vtip * NR

           Blade_drag     = Blade_drag + D_blade 

#====================================================================
# tilting rotor: check helical tip mach number for limit
#====================================================================

        elif group.type == 'tilting':

          Vmax               = sqrt(Vtip*Vtip + Vcruise*Vcruise)
          if(Vmax > VtipMax):
              Vtip = sqrt(Vmax*Vmax - Vcruise*Vcruise)

          if Vtip < 0.2*group.tipspeed:
             print ('warning: SLOWED ROTOR BELOW 20% RPM')
             print ('hitting min limit: 20% HOVER TIP SPEED')
             Vtip = 0.2*group.tipspeed

#====================================================================
# unknown rotor type
#====================================================================

        else:
          quit('I dont know this rotor type')       

     return Blade_drag

#====================================================================
# Individual rotor groups are contained in the class "rotor_group"
#====================================================================

class rotor_group:

#====================================================================
# function to initialize all segments
#====================================================================

   def __init__(self, data, key, nseg):

#=======================================================================
#geometric parameters
#=======================================================================

      self.nrotors        = 0
      self.motor_group_id = -1    # which motor group is sized by this rotor
      self.wing_group_ids = []
      self.nblade         = 0
      self.set_radius     = False 
      self.span_driven    = False 
      self.set_DL         = False
      self.set_BL         = False 
      self.set_Vtip       = False
      self.set_sigma      = False
      self.type           = 'tilting'
      self.key            = key   
      self.radius         = 0.0
      self.chord          = 0.0
      self.area           = 0.0
      self.solidity       = 0.0
      self.diameter       = 0.0
      self.tipspeed       = numpy.zeros(nseg)
      self.max_Power      = 0.0
      self.thrust_share   = 'wing_based'
      
#=======================================================================
# operational parameters
#=======================================================================

      self.ctsigma        = 0.0     # in sizing segment
      self.diskloading    = 0.0     # for sizing segment
      self.nu_beta        = 1.0     # in hover, flap nat freq in /rev
      self.ipf            = 1.0
      self.cd0            = 0.01
      self.RPM_ratio      = 1.0     # cruise to hover RPM ratio

#=======================================================================
# overload factors for thrust, torque and RPM
#=======================================================================

      self.torque_scaling = 1.0
      self.T_overload     = 1.0
      self.Q_overload     = 1.0
      self.RPM_overload   = 1.0
            
#=======================================================================
#added by AS for HYDRA 2.0
#=======================================================================

      self.blade_mass   = 0.0       # in kilograms
      self.spl_hover    = 0.0       # SPL estimate
      self.ainf         = 0.0 

#===============================================================================
# size the rotor with either radius or disk loading
#===============================================================================

   def sizing(self, thrust, rho, Rmax, wing_groups, clearance, bfus):

#====================================================================
# calculate radius from disk loading
#====================================================================

     if self.set_DL:
        self.area       = thrust/self.diskloading         
        self.radius     = sqrt(self.area/pi)

#====================================================================
# radius given directly
#====================================================================

     elif self.set_radius:
        R               = self.radius

#====================================================================
# span-driven sizing
# loop over all wings this rotor appears on, find the multiplier
# for rotor radius that identifies span required
# (span - fuselage)/2 = available length along which rotors can be placed
# this value is equal to multiplier * radius, hence find radius
#====================================================================

     elif self.span_driven:
        wgids           = self.wing_group_ids
        size            = 1.0 + clearance*0.5     # clearance for rotor on each side
        radius          = 1.0e9
        for wgid in wgids:
            group       = wing_groups[wgid]
            nr          = group.nrotors/group.nwings
            multiplier  = size*float(nr) - 1.0
            Rmax        = (group.span - bfus)*0.5/multiplier
            radius      = min(radius, Rmax)

        # print('rotor radius from span driven sizing is ',radius)
        self.radius     = radius 

#====================================================================
# error message
#====================================================================

     else:
        quit('CRITICAL ERROR: EITHER SET RADIUS OR DL, or enable span-driven rotor sizing')

#====================================================================
# cap max rotor size
#====================================================================
    
     if(self.radius > Rmax):
        self.radius     = Rmax 
        print('capping rotor radius to max available')

     R                  = self.radius
#     print('rotor radius after sizing is ',R)
#set diameter
     self.diameter      = 2.0*self.radius
     A                  = pi * self.radius * self.radius 
     self.area          = A
     self.diskloading   = thrust/A

#====================================================================
# rotor blade loading or tip speed in hover 
#====================================================================

     if self.set_Vtip:
        CT            = thrust/(rho*A*self.tipspeed*self.tipspeed) 
     else:
        quit('CRITICAL ERROR: need to know tip speed for sizing')

#====================================================================
# need blade loading or solidity
#====================================================================

     if self.set_BL:
        self.solidity = CT/self.ctsigma 
     elif self.set_sigma:
        self.ctsigma  = CT/self.solidity
     else:
        quit('SET EITHER BLADE LOADING or ROTOR SOLIDITY')

#====================================================================
# Main rotor chord, SI (m)
#====================================================================

     self.aspectratio = self.nblade/(pi*self.solidity)
     self.chord       = self.radius / self.aspectratio
     self.ctsigma     = CT/self.solidity

#===============================================================================
# 15% interference penalty for induced power of coaxial
#===============================================================================

     if self.type == 'coaxial':
        self.kint     = 1.16 

     return None

#===============================================================================
# hover power calculation for the rotor
#===============================================================================

   def hover_power(self, thrust, rho, FM, use_bemt, size_switch):

#====================================================================
# Calculate thrust coefficient and profile power/induced power coeffs
#====================================================================

      CT      = thrust/(rho*self.area*self.tipspeed*self.tipspeed)

      Cpo     = self.solidity * self.cd0 / 8.0     # profile power
      Cpi     = CT*sqrt(CT*0.5)*self.kint          # ideal power + induced losses
#      print(CT,rho,self.area,self.tipspeed,self.tipspeed)

#====================================================================
# Compute power coefficient for all rotors wrt disk area of one rotor
# if using BEMT, replace aero efficiency with calibrated value
# otherwise, use momentum theory to find FM and power
#====================================================================

      if use_bemt:
         Cptotal           = Cpi/FM
         self.ipf          = (Cptotal - Cpo)/Cpi      # equiv. ipf
      else:
         if FM > 0.0:
            Cptotal        = Cpi/FM
#            print(Cpi,Cptotal);quit('ok?')
            # print(Cptotal,Cpo,Cpi)
            self.ipf       = (Cptotal - Cpo)/Cpi      # equiv. ipf
         else:                            # kappa is given
            Cptotal        = (Cpi*self.ipf + Cpo) 
            FM             = Cpi/Cptotal

#====================================================================
# calculate the FM in hover (recognize that if the mission has
# multiple hover segments, this value may change in each segment)
#====================================================================

      if size_switch:
         self.fm = FM
      Vtip       = self.tipspeed
      Phover     = Cptotal* rho*self.area*Vtip*Vtip*Vtip
      return Phover

#===============================================================================
# function to calculate RPM and thrust requirements for one rotor out conditions
#===============================================================================

   def one_rotor_out(self, wing_groups):

      nr_all            = 0
      nrmax             = 0
      for key,group in wing_groups.items():
         nr             = group.nrotors
         nr_all         = nr_all+nr           # rotors per wing 
#         nrmax          = max(nrmax,nr)       # identify min rotors in a set 
#      print(nr_all,nrmax)
      nrmax             = nr_all
      self.T_overload   = float(nrmax/(nrmax-2))+0.2    # emergency to nominal thrust ratio

#===============================================================================
# assuming this extra thrust is produced at same CT/sigma by ramping  up RPM,
# find RPM ratio to achieve this thrust overload
#===============================================================================

      self.P_overload   = 1.2*self.T_overload**1.5 

#===============================================================================
# for same FM, find torque overload
# emergency to nominal hover torque ratio
#===============================================================================

      self.RPM_overload = numpy.sqrt(self.T_overload)*1.1

      self.Q_overload   = self.P_overload/numpy.sqrt(self.T_overload)

#      print(self.P_overload,self.RPM_overload,self.Q_overload);quit()
#===============================================================================
# Compute noise from a single rotor in hover
#===============================================================================
# Based on propeller noise from "A Review of Aerodynamic Noise from Propellers,
# Rotors, and Lift Fans" JPL Tech report 32-1462
# Inputs
#  T        - thrust [N]
#  P        - power [kW], from a SINGLE rotor
#  omega    - fan rotational speed [rad/s]
#  geom     - fan geometry structure
#  air      - air property structure
#  S        - distance from observer location to noise source [m]
#  theta    - angle from thrust direction to observer location [rad]
#
# Calculates and stores:
#   spl - Sound pressure level in dB
#
#===============================================================================

# Constants

   def fanNoise(self,S,theta):

        NtoLb = 0.224809           # Newtons to lb
        WtoHp = 0.00134102         # Watts to hp
        m2f   = 1.0/0.3048         # meters to feet
        Pref  = 0.0002             # dynes/cm^2
        k     = 6.1e-27            # proportionality constant

        # print('fan noise model is switched off')

#===============================================================================
# Average blade CL
#===============================================================================

        CL      = 6 * self.ctsigma    # Ct / geom.solidity

#===============================================================================
# Thrust in pounds
#===============================================================================

        T_lb    = self.thrust[0] * NtoLb

#===============================================================================
# Parameters
#===============================================================================

        VTip    = self.tipspeed                 # tip speed, m/s
        S_ft    = S *m2f                        # Distance to source in ft
        R_ft    = self.radius * m2f             # Radius in ft
        A_ft2   = self.area* m2f*m2f            # Disc area in ft^2
        P_Hp    = self.power * WtoHp            # Power in Hp
        M_t     = VTip /self.ainf               # tip mach number
        A_b     = A_ft2 * self.solidity         # blade area, sq. ft

#        print(self.power*0.001)
#        print(self.solidity)
#        print(S_ft, R_ft, A_ft2, P_Hp, M_t, A_b)
#        quit()
#===============================================================================
# Compute rotational noise
#===============================================================================
        
        m_max   = 10                    # Maximum harmonic number
        p_m     = numpy.zeros(m_max)
        for m in range(1,m_max+1):
            p_m[m-1] = 169.3 * m * self.nblade * R_ft * M_t / (S_ft * A_ft2) *       \
                    (0.76 * P_Hp / M_t**2 - T_lb * cos(theta)) *              \
                    besselj(m * self.nblade, 0.8 * M_t * m * self.nblade * sin(theta))

#===============================================================================
# Compute total RMS sound pressure level in dynes/cm^2 = 0.1 Pa
#===============================================================================

        p = numpy.sqrt(numpy.sum(numpy.square(p_m)))

#===============================================================================
# Convert to SPL in dB's
#===============================================================================

        SPL_rotational = 20 * numpy.log10(p / Pref)

#===============================================================================
# Vortex Noise (pg. 11)
# SPL = 10 * log(k * A_b * V_0.7^6 / 10^-16) + 20 * log(CL/0.4) [dB at 300 ft]
# where
# k = constant of proportionality = 6.1e-27
# k/1e-16 = 6.1e-11
# A_b = propeller blade area [ft^2]
# V_0.7 = velocity at 0.7 * radius
#===============================================================================

#        SPL_vortex  = 10 * log10(k * A_b * (0.7 * VTip)**6 / 1e-16) + 20 * log10(CL / 0.4)
        SPL_vortex  = 10 * log10(6.1e-11 * A_b * (0.7 * VTip)**6) + 20 * log10(CL / 0.4)

# Total noise
# Adding atmospheric attenuation of 6 dB for every doubling of distance from reference
        spl = 10 * log10( 10**(SPL_rotational / 10) + 10**(SPL_vortex / 10)) - 6 * log2(S_ft/300)

        return spl

#====================================================================
# Blade weight modeL: physics-based; estimation of all major loads
#====================================================================

   def bladewt(self, tech_factor):

#====================================================================
# Loop over rotor groups to perform sizing 
#====================================================================

      nb          = self.nblade 
      Omega       = self.tipspeed/self.radius*self.RPM_overload
      Fz          = self.sizing_thrust/self.nblade*self.T_overload      
      Mz          = self.sizing_torque/self.nblade*self.Q_overload                  # torque at root 
      nz          = 1.5                           # safety factor
      nr          = self.nrotors                  # rotors in this group

      #   material    = 'Aluminum'
      #   material    = 'Titanium'
      #   material    = '090_carbon'
      material    = 'uniaxial_carbon'

#====================================================================
# Set dictionary of inputs for blade weight model
#====================================================================

      Rotor       ={ 'R': self.radius,        # Radius, m
                 'Omega': Omega,              # Rotation speed, rad/s
                 'chord': self.chord,         # chord, m 
                    'Fz': Fz,                 # Vertical lift/blade, N
                    'Mz': Mz,                 # Torque at root, N-m
                    'nz': nz,                 # load factor 
                    'nr': nr,                 # number of rotors 
                    'nb': self.nblade,        # number of blades
                   'nub': self.nu_beta,       # flap frequency
                   'mat': material,
                 'sigma': self.solidity}

#      print('hello motor',self.tipspeed,self.radius,self.RPM_overload);quit()

#      print(self.sizing_torque,self.nblade,self.Q_overload)
#      print(self.tipspeed,self.radius,self.RPM_overload)
#      quit()
#      print(Rotor);quit()
      rotor_mass  = blade_wt_modelv2(Rotor, tech_factor)

      return rotor_mass           # in kg 


#==================================================================
# ANTI-ICING WEIGHT FOR BLADES
#==================================================================
#==================================================================
# fixed parameters
#==================================================================

   def icing_weight(self, tech_factor):

#==================================================================
# unpack inputs
#==================================================================

      R            = self.radius*m2f
      chord        = self.chord*m2f

#==================================================================
# derived quantities
#==================================================================

      Ablades      = self.nrotors  *self.radius * self.chord*m2f*m2f # blade plan-form area, sq.ft 

#==================================================================
# find weight of blade de-icing wires+system: some const * blade area
#==================================================================

      wght_DIelect = k_elect * Ablades*tech_factor
      wght_DIsys   = k_rotor * Ablades*tech_factor

      deicing      = {'deicing_blades': wght_DIelect*lb2kg, 
                       'deicing_equip': wght_DIsys  *lb2kg}

      return deicing
