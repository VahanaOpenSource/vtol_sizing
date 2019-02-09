import sys, numpy 
from numpy import sin, cos, log10, log2, sqrt, pi
from scipy.special import jv as besselj
sys.path.insert(0,'../Stage_0/')

from conversions import *
k_elect = 0.125 # proportionality constants
k_rotor = 0.125


from blade_wt_modelv2 import blade_wt_modelv2

from rotor_wt    import rotor_weight
from conversions import *
from cost_class  import dict_accumulation

#====================================================================
# constants for rotor control weight estimates
#====================================================================

f_RWnb  = 1.25 # fraction rotary wing non-boosted weight
f_RWhyd = 0.4  # fraction rotary wing hydraulic weight
f_RWred = 3.0  # redundancy factor
f_mbsv  = 1.3029  # ballistic survivability
f_bsv   = 1.117  # ballistic survivability

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

#loop over rotor groups, find design parameters
      for key in sorted(data):
         self.groups[ngrp]  = rotor_group(data[key], key, nseg)
         ngrp               = ngrp + 1

      self.ngroups      = ngrp

      return None
      
#====================================================================
# Function to calculate blade drag in forward flight 
# for edgewise rotors
#====================================================================

   def blade_drag(self, Vcruise, VtipMax, rho):

     """
     Inputs
       1. Vcruise (m/s)
       2. VtipMax (m/s): maximum rotor section speed at advancing blade tip
       3. rho     (kg/cu.m): air density

     Outputs
       1. Blade_drag (Newtons)
     """

# loop over all rotor groups, and calculate rotor drag
     Blade_drag           = 0.0

     for i in range(self.ngroups):
        group             = self.groups[i]
        NR                = group.nrotors


        Vtip              = group.tipspeed*group.RPM_ratio

# edgewise rotor: get average profile drag coefficient of sections
        if group.type == 'lift':
           if group.tipspeed + Vcruise > VtipMax:
              Vtip              = VtipMax - Vcruise
              group.RPM_ratio   = Vtip/group.tipspeed

# find average section drag coefficient at different advance ratios
           cd0                  = group.cd0
           muCruise             = Vcruise/Vtip
           if muCruise > 1.0:
              cd0 = 1.5*cd0
           elif muCruise > 0.3:
              cd0 = cd0*(1.0 + (muCruise-0.3)/0.7e0)
           else:
              cd0 = cd0

           cd0    = 0.012 

# find blade drag in wind direction
           D_blade        = group.solidity*cd0/8.0*(3.1*muCruise)
           D_blade        = D_blade * rho * group.area * Vtip*Vtip * NR

           Blade_drag     = Blade_drag + D_blade 

# tilting rotor: check helical tip mach number for limit
        elif group.type == 'tilting':

          Vmax               = sqrt(Vtip*Vtip + Vcruise*Vcruise)
          if(Vmax > VtipMax):
              Vtip = sqrt(Vmax*Vmax - Vcruise*Vcruise)

          if Vtip < 0.2*group.tipspeed:
             print ('warning: SLOWED ROTOR BELOW 20% RPM')
             print ('hitting min limit: 20% HOVER TIP SPEED')
             Vtip = 0.2*group.tipspeed

# unknown rotor type
        else:
          quit('I dont know this rotor type')       

     return Blade_drag

#====================================================================
# function to accumulate weights of all rotor groups into one dictionary
#====================================================================

   def weight_rollup(self, wing, fuselage, tech_factors, aircraftID):
  
      """
      This function calculates weights of all rotor-related groups on the aircraft
      and returns a dictionary containing the breakdown and total weight
      
      Blades, hubs, actuators, anti-icing systems, flight controls, hydraulics are included.

      Inputs:
      1. wing           : class containing wing-related information
      2. fuselage       : class containing fuselage-related information
      3. tech_factors   : object containing vehicle technology factors that 
                          are weight multipliers for various empty weight groups
      4. aircraftID     : type of aircraft - 1 = elastic rotor (AFDD model), 
                          anything else = inelastic rotor (eVTOL)

      Outputs:
      1. rotor_weight   : dictionary containing breakdown of weights and total weight
                          of all rotor-related components
      """
      
      rotor_weight  = {}
      ngroups       = self.ngroups
      for i in range(ngroups):
          group     = self.groups[i]

# loop over wing groups that contain this rotor and find the sizing thrust, 
# sizing torque; these values are the max vals across all segments too..
          sizing_T  = 0.0
          sizing_Q  = 0.0
          for wgid in group.wing_group_ids:
              T_try    = numpy.amax(wing.groups[wgid].rotor_thrust)
              Q_try    = numpy.amax(wing.groups[wgid].rotor_torque)

              sizing_T = max(sizing_T, T_try)
              sizing_Q = max(sizing_Q, Q_try)
# also check for loads from fuselage-mounted rotors belonging to this group
          if(group.fuse_group_id != -1):
              T_try    = numpy.amax(fuselage.rotor_thrust)
              Q_try    = numpy.amax(fuselage.rotor_torque)

              sizing_T = max(sizing_T, T_try)
              sizing_Q = max(sizing_Q, Q_try)

# Set sizing thrust and torque for the group, then perform sizing
          group.sizing_thrust   = sizing_T 
          group.sizing_torque   = sizing_Q 

# for single MR configuration, use AFDD model            
          if(aircraftID == 1):
            rotor_wt            = group.afdd_rotor_wt(aircraftID, tech_factors.rotor)

# eVTOLs: use physics-based model
# Note: collective actuators included in the new weight model

          else:
            rotor_wt            = group.rotorwt(tech_factors.rotor)

            # print(rotor_wt);quit('ok rotor weight?')
# remember weight breakdown for the group in a dictionary
          for key,value in rotor_wt.items():
             k2     = 'group'+str(i)+key 
             rotor_weight[k2] = value 

# add de-icing systems weight
          anti_icing = group.icing_weight(tech_factors.anti_icing)

          for key,value in anti_icing.items():
             k2     = 'group'+str(i)+'anti_icing_'+key 
             rotor_weight[k2] = value 

# add flight controls for this rotor group
          if(group.nrotors == 1 or aircraftID == 1):
              flt_controls = group.conventional_rotor_control_wt()

              for key,value in flt_controls.items():
                  k2     = 'group'+str(i)+'controls_'+key 
                  rotor_weight[k2] = value 

#add up all the entries in the rotor groups
      rotor_weight['total']   = dict_accumulation(rotor_weight)

      return rotor_weight

#====================================================================
# Individual rotor groups are contained in the class "rotor_group"
#====================================================================

class rotor_group:

   def __init__(self, data, key, nseg):
      """
      function to initialize a rotor group
      Inputs:
      1.  data:     dictionary containing rotor group details 
      2.  key:      rotor group name 
      3.  nseg:     number of mission segments
      """
# rotor group parameters
      self.nrotors        = 0
      self.wing_group_ids = []
      self.fuse_group_id  = -1
      self.xmsn_group_id  = -1
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
      self.p_req          = numpy.zeros(nseg)
      self.q_req          = numpy.zeros(nseg)
      self.max_Power      = 0.0
      self.thrust_share   = 'wing_based'
      
# operational parameters
      self.ctsigma        = 0.0     # in sizing segment
      self.diskloading    = 0.0     # for sizing segment
      self.nu_beta        = 1.0     # in hover, flap nat freq in /rev
      self.ipf            = 1.0
      self.cd0            = 0.01
      self.RPM_ratio      = 1.0     # cruise to hover RPM ratio

# overload factors for thrust, torque and RPM
      self.torque_scaling = 1.0
      self.T_overload     = 1.0
      self.Q_overload     = 1.0
      self.RPM_overload   = 1.0
      self.blade_mass     = 0.0       # in kilograms
      self.spl_hover      = 0.0       # SPL estimate in hover
      self.ainf           = 0.0 

#===============================================================================
# Rotor sizing function
#===============================================================================

   def sizing(self, thrust, rho, Rmax, wing_groups, clearance, bfus):

     """
     This function sizes the rotor radius and chord with either span-driven 
     sizing, disk loading or from a given radius.

     Inputs:
     1.   thrust      : target thrust in Newtons 
     2.   rho         : hover air density, kg/cu.m 
     3.   Rmax        : maximum rotor radius, meters
     4.   wing_groups : array of wing group classes defined in wing_class.py
     5.   clearance   : fraction of rotor radius clearance between rotor planes/fuselage
     6.   bfus        : fuselage width, meters
     """

# case 1: calculate radius from disk loading
     if self.set_DL:
        self.area       = thrust/self.diskloading         
        self.radius     = sqrt(self.area/pi)

# case 2: radius given directly
     elif self.set_radius:
        R               = self.radius

# case 3: span-driven sizing
# loop over all wings this rotor appears on, find the multiplier
# for rotor radius that identifies span required
# (span - fuselage)/2 = available length along which rotors can be placed
# this value is equal to multiplier * radius, hence find radius
     elif self.span_driven:
        wgids           = self.wing_group_ids
        size            = 1.0 + clearance*0.5     # clearance for rotor on each side
        Rmin            = 1.0e9
        for wgid in wgids:
            group       = wing_groups[wgid]
            nr          = group.nrotors/group.nwings

            multiplier  = size*float(nr) - 1.0
            radius      = (group.span - bfus)*0.5/multiplier
            Rmin        = min(radius, Rmin)
            #print('group',wgid,'SPAN = ',group.span,multiplier,radius,Rmin)
        #print('rotor radius from span driven sizing is ',Rmin)
        #x1=input('?')
        self.radius     = Rmin 

# error message
     else:
        quit('CRITICAL ERROR: EITHER SET RADIUS OR DL, or enable span-driven rotor sizing')

# cap max rotor size
     if(self.radius > Rmax):
        self.radius     = Rmax 
        print('capping rotor radius to max available')

     R                  = self.radius
#     print('rotor radius after sizing is ',R)

#set diameter, area and disk loading
     self.diameter      = 2.0*self.radius
     A                  = pi * self.radius * self.radius 
     self.area          = A
     self.diskloading   = thrust/A

# rotor blade loading or tip speed in hover 
     if self.set_Vtip:
        CT            = thrust/(rho*A*self.tipspeed*self.tipspeed) 
     else:
        quit('CRITICAL ERROR: need to know tip speed for sizing')

# need blade loading or solidity
     if self.set_BL:
        self.solidity = CT/self.ctsigma 
     elif self.set_sigma:
        self.ctsigma  = CT/self.solidity
     else:
        quit('SET EITHER BLADE LOADING or ROTOR SOLIDITY')

# Main rotor chord, SI (m)
     self.aspectratio = self.nblade/(pi*self.solidity)
     self.chord       = self.radius / self.aspectratio
     self.ctsigma     = CT/self.solidity

# 15% interference penalty for induced power of coaxial
     if self.type == 'coaxial':
        self.kint     = 1.16 

     return None

#===============================================================================
# hover power calculation for the rotor
#===============================================================================

   def hover_power(self, thrust, rho, FM, use_bemt, size_switch):
      """
      this function calculates the power required to hover for a rotor
      Inputs:
      1.  thrust:       required thrust, Newtons 
      2.  rho:          ambient density, kg/cu.m 
      3.  FM:           rotor hover figure of merit 
      4.  use_bemt:     logical flag to specify if bemt is used to calculate FM 
      5.  size_switch:  logical flag to specify if this performance calculation is for a sizing segment 

      Output:
      1.  Phover:   hover shaft power in watts
      """

# Calculate thrust coefficient and profile power/induced power coeffs
      CT      = thrust/(rho*self.area*self.tipspeed*self.tipspeed)

      Cpo     = self.solidity * self.cd0 / 8.0     # profile power
      Cpi     = CT*sqrt(CT*0.5)*self.kint          # ideal power + induced losses

# Compute power coefficient for all rotors wrt disk area of one rotor
# if using BEMT, replace aero efficiency with calibrated value
# otherwise, use momentum theory to find FM and power
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

# calculate the FM in hover (recognize that if the mission has
# multiple hover segments, this value may change in each segment)
      if size_switch:
         self.fm = FM
      Vtip       = self.tipspeed
      Phover     = Cptotal* rho*self.area*Vtip*Vtip*Vtip
      return Phover

#===============================================================================
# function to calculate RPM and thrust requirements for one rotor out conditions
#===============================================================================

   def one_rotor_out(self, wing_groups):

# for electric power input, some numbers may be relevant
      if(self.input_power == 'electric'):
        nr_all            = 0
        nrmax             = 0
        for key,group in wing_groups.items():
           nr             = group.nrotors
           nr_all         = nr_all+nr           # rotors per wing 
           nrmax          = max(nrmax,nr)       # identify min rotors in a set 

        self.T_overload   = float(nrmax/(nrmax-2))+0.2    # emergency to nominal thrust ratio

# assuming this extra thrust is produced at same CT/sigma by ramping  up RPM,
# find RPM ratio to achieve this thrust overload
        self.P_overload   = 1.2*self.T_overload**1.5 

# for same FM, find torque overload
# emergency to nominal hover torque ratio
        self.RPM_overload = numpy.sqrt(self.T_overload)*1.1
        self.Q_overload   = self.P_overload/numpy.sqrt(self.T_overload)

# for mechanical input power, these numbers are not relevant

#===============================================================================
# Compute noise from a single rotor in hover
# Based on propeller noise from "A Review of Aerodynamic Noise from Propellers,
# Rotors, and Lift Fans" JPL Tech report 32-1462
#===============================================================================

   def fanNoise(self,S,theta, Thrust, power_Watts):

      """      
      Inputs
        1.S           - distance from observer location to noise source [m]
        2.theta       - angle from thrust direction to observer location [rad]
        2.Thrust      - thrust [N]
        3.power_watts - power in watts, from a SINGLE rotor

      Calculates and stores:
        spl - Sound pressure level in dB
      """

      NtoLb = 0.224809           # Newtons to lb
      WtoHp = 0.00134102         # Watts to hp
#        m2f   = 1.0/0.3048         # meters to feet
      Pref  = 0.0002             # dynes/cm^2
      k     = 6.1e-27            # proportionality constant

# Average blade CL
      CL      = 6 * self.ctsigma    # Ct / geom.solidity

# Thrust in pounds
      T_lb    = Thrust * NtoLb

# Parameters
      VTip    = self.tipspeed                 # tip speed, m/s
      S_ft    = S *m2f                        # Distance to source in ft
      R_ft    = self.radius * m2f             # Radius in ft
      A_ft2   = self.area* m2f*m2f            # Disc area in ft^2
      P_Hp    = power_Watts * WtoHp           # Power in Hp
      M_t     = VTip /self.ainf               # tip mach number
      A_b     = A_ft2 * self.solidity         # blade area, sq. ft

#        print(self.power*0.001)
#        print(self.solidity)
#        print(S_ft, R_ft, A_ft2, P_Hp, M_t, A_b)
#        quit()
# Compute rotational noise
      m_max   = 10                    # Maximum harmonic number
      p_m     = numpy.zeros(m_max)
      for m in range(1,m_max+1):
          p_m[m-1] = 169.3 * m * self.nblade * R_ft * M_t / (S_ft * A_ft2) *       \
                    (0.76 * P_Hp / M_t**2 - T_lb * cos(theta)) *              \
                      besselj(m * self.nblade, 0.8 * M_t * m * self.nblade * sin(theta))

# Compute total RMS sound pressure level in dynes/cm^2 = 0.1 Pa
      p = numpy.sqrt(numpy.sum(numpy.square(p_m)))

# Convert to SPL in dB's
      SPL_rotational = 20 * numpy.log10(p / Pref)

# Vortex Noise (pg. 11)
# SPL = 10 * log(k * A_b * V_0.7^6 / 10^-16) + 20 * log(CL/0.4) [dB at 300 ft]
# where
# k = constant of proportionality = 6.1e-27
# k/1e-16 = 6.1e-11
# A_b = propeller blade area [ft^2]
# V_0.7 = velocity at 0.7 * radius

#        SPL_vortex  = 10 * log10(k * A_b * (0.7 * VTip)**6 / 1e-16) + 20 * log10(CL / 0.4)
      SPL_vortex  = 10 * log10(6.1e-11 * A_b * (0.7 * VTip)**6) + 20 * log10(CL / 0.4)

# Total noise
# Adding atmospheric attenuation of 6 dB for every doubling of distance from reference
      spl = 10 * log10( 10**(SPL_rotational / 10) + 10**(SPL_vortex / 10)) - 6 * log2(S_ft/300)

      return spl

#====================================================================
# Blade weight model: physics-based; estimation of all major loads
#====================================================================

   def rotorwt(self, tech_factor):

      """
      This function calculates the weight of the rotor assembly (hub, blades, 
      and collective pitch actuator for eVTOLs)

      Input:
      1. tech_factor: a scaling factor that increases/decreases final predicted 
                      weight to account for materials/manufacturing improvements
      """

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

# Set dictionary of inputs for blade weight model
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

      rotor_mass, total  = blade_wt_modelv2(Rotor, tech_factor)
      self.mass_assembly = total/nr           # mass of single rotor+hub+actuator assembly, kg
      self.mass_blades   = rotor_mass['blades']
      return rotor_mass           

#==================================================================
# rotor system weight, AFDD model
#==================================================================

   def afdd_rotor_wt(self, aircraftID, tech_factor):

      """
      This function calculates the weight of the rotor system (hubs, blades, 
      spinners, folding systems) for the present rotor group
      Inputs are:
      1. aircraftID: integer flag for legacy aircraft; [1], [2] or [3]
                     [1] = Single main rotor / tail rotor helicopter 
                     [2] = Tilt-rotor/tilt-wing aircraft
                     [3] = Coaxial main rotor system
      2. tech_factor: multiplier for rotor system weights to scale results from 
         parametric models up or down
      
      Note that the inputs for the rotor weight model are in FPS

      Output is a dictionary containing the total weight of the rotor system 
      and the breakdown of this weight into blades, hubs, spinners and folding systems
      """
      inputs     = {'aircraftID'  : aircraftID,           \
                    'nblade'      : self.nblade,          \
                    'nrotor'      : self.nrotors,         \
                    'radius'      : self.radius*m2f,      \
                    'chord'       : self.chord*m2f,       \
                    'vtip'        : self.tipspeed*m2f,    \
                    'nu_blade'    : self.nu_beta,         \
                    'tech_factor' : tech_factor     }

# get breakdown for all rotors, mass of one assembly in kg
      breakdown, self.mass_assembly = rotor_weight(inputs)
      self.mass_blades   = breakdown['blades']
      return breakdown

#==================================================================
# ANTI-ICING WEIGHT FOR BLADES
#==================================================================

   def icing_weight(self, tech_factor):
      """
      This function calculates the weight of anti-icing equipment 
      for a rotor blade set. 
      Input:
      1. tech_factor: multiplier for deicing system weights to scale results from 
         parametric models up or down

      Output:
      1. deicing: dictionary with constitutent weight components for blade anti-icing
      """

      R            = self.radius*m2f
      chord        = self.chord*m2f
      Ablades      = self.nblade*self.nrotors  *self.radius * self.chord*m2f*m2f # blade plan-form area, sq.ft 

      wght_DIelect = k_elect * Ablades*tech_factor

# heating element weights 
      wght_DIsys   = k_rotor *  Ablades*tech_factor

#total weight done outside for the rotor
      # total_wt     = wght_DIelect + wght_DIsys 

      deicing      = {'blades': wght_DIelect*lb2kg, 
                       'equip': wght_DIsys  *lb2kg}

# increment assembly mass
      self.mass_assembly  = self.mass_assembly + lb2kg*(wght_DIsys + wght_DIelect)/self.nrotors
      return deicing

#==================================================================
# weight model for hydraulic control mechanisms
#==================================================================

   def conventional_rotor_control_wt(self):

      """
      this function estimates the weight of boost mechanisms, hydraulics
      and non-boosted controls for a conventional rotor with swashplates

      Output:
      1. flt_control: dictionary with weights of controls and hydraulics
      """

      nrotor = self.nrotors
      nblade = self.nblade
      chord  = self.chord*m2f
      V_tip  = self.tipspeed*m2f

#rotary-wing flight control mechanisms
      w_fc = ( 0.2873 * f_mbsv * (nrotor * nblade)**0.6257 *
               chord**1.3286 * (0.01*V_tip)**2.1129 *
               f_RWred**0.8942 )

#boosted flight control weight for rotary-wing
      wght_RWb  = ( 0.02324 * f_bsv * (nrotor*nblade)**1.0042 *
                  nrotor**0.1155 * chord**2.2296 * 
                  (0.01*V_tip)**3.1877 )

#boost mechanism weight for rotary wing
      wght_RWmb = (1-f_RWhyd)*w_fc

#rotary-wing non boosted flight control weight
# non-boosted control weight is a parameter value x flight control weight  
      wght_RWnb = f_RWnb*(1-f_RWhyd)*w_fc
   
# total control weight = boosted controls + boost mechanisms + non-boosted 
      controls  = wght_RWb + wght_RWmb + wght_RWnb 

# hydraulics: for rotary-wing and conversion mechanisms
# based on fraction of flight control weight and fraction of conversion 
# boost mechanism weight
      wght_RWhyd      = f_RWhyd*w_fc

# assemble outputs in a dictionary for return
      flt_control = {'hydraulics': wght_RWhyd*lb2kg, 
                     'mechanisms': controls  *lb2kg}

# increment rotor assembly mass
      self.mass_assembly = self.mass_assembly + lb2kg*(controls/nrotor)
      return flt_control