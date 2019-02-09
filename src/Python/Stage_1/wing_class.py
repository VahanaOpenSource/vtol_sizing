import sys, numpy 
from numpy import sin, cos, log10, log2, sqrt, pi
from scipy.special import jv as besselj
sys.path.insert(0,'../Stage_0/')
from conversions import *

from motor_mount_mass import motor_mount_mass as mmmass, spar_mass
from cost_class       import dict_accumulation

#====================================================================
# sizing constants for weight models
#====================================================================

f_LGloc     = 1.0                # 1.7247 if landing gear on wing, 1.0 otherwise
k_elect     = 0.5                # proportionality constants for anti-icing weight
k_rotor     = 0.5

#model       = 'statistical'
model       = 'target_freq'
fn          = 4.5                # wing spar freq, Hz
f           = 7.0                # motor mount freq, Hz

tau_w       = 0.158              # wing t/c ratio
Lamda       = 0.0                # no sweep
taper       = 0.8                # taper ratio = ctip/croot, linear
nz          = 3.8                # max design load factor 
b_fold      = 0.0                # fraction of span that folds (0 = not folding)

#====================================================================
# this class contains information for a fixed wing group
#====================================================================

class wing_group:

   def __init__(self, data, key, nseg):

      """
      initialization function for wing group
      """
      self.nwings         = data['nwing'][0]
      self.span           = 0.0
      self.chord          = 0.0
      self.area           = 0.0
      self.oswald         = 0.0
      self.K              = 0.0
      self.cd0            = 0.0
      self.aspectratio    = 0.0
      self.lift_frac      = 0.0                       # lift fraction carried by wing, 0 to 1, cruise
      self.stall_speed    = 80.0                      # stall speed in knots
      self.stall_speed    = self.stall_speed*0.5144   # in m/s
      self.CLmax          = 1.35

# operational parameters
      self.lift           = numpy.zeros(nseg)
      self.drag           = numpy.zeros(nseg)
      self.key            = key
      self.rotor_group_id = -1 

#remember all rotors for this wing group, including all duplicate wings in this group
      self.nrotors        = data['nrotors'][0]*self.nwings # left+right wings!

#rotor performance for units within a group
      self.rotor_thrust   = numpy.zeros(nseg)
      self.rotor_power    = numpy.zeros(nseg)
      self.motor_power    = numpy.zeros(nseg)
      self.rotor_torque   = numpy.zeros(nseg)
      self.rotor_aero_eta = 0.75*numpy.ones(nseg)

# rotor hub locations, rotor + mount assembly masses and design max thrusts: for a half-wing
      nr_halfwing         = int(self.nrotors/(self.nwings*2))
      lstore              = max(1,nr_halfwing)
      self.ycoords        = numpy.zeros(lstore) 
      self.masses         = numpy.zeros(lstore)
      self.rotor_T        = numpy.zeros(lstore)
      return None

#====================================================================
# function to size a fixed wing
#====================================================================

   def size_wing(self, lift, q, cl, df):

      """ 
      This function calculates the dimensions of a fixed wing given the 
      target lift to be supported at a given dynamic pressure and CL

      Inputs:
      1. lift  : wing lift (Newtons)
      2. q     : dynamic pressures (N/sq.m)
      3. cl    : target lift coefficient 
      4. df    : fuselage width (equivalent diameter, meters) --> for Oswald efficiency

      Outputs: None (information stored within class)
      """

      self.area     = lift/(q*cl)                        # area of each wing (left+right)
      self.span     = sqrt(self.aspectratio * self.area) # span of each wing (left+right)
      self.chord    = self.area / self.span

# Oswald efficiency calculation
# use P-Q combination from following reference to calculate oswald eff
# http://www.fzt.haw-hamburg.de/pers/Scholz/OPerA/OPerA_PUB_DLRK_12-09-10.pdf
# Equations (5), (6), (7) and (8) 

      if(self.oswald == 0.0):
         u           = 0.99                  # from Kroo [2] of above ref.
         s           = df/self.span
         s           = 1.0 - 2.0*s*s         # eqn (6)
         Q           = 1.0/(u*s)             # eqn (5)
         P           = 0.38*0.02             # eqn (7), CD0=0.02
         self.oswald = 1.0/(Q + P*pi*self.aspectratio)
      self.K         = 1.0/(pi*self.aspectratio*self.oswald)

      return None 

#====================================================================
# calculate lift for a fixed wing with CL cap
#====================================================================

   def calculate_lift(self, L_req, q, CLmax):

      """ this function calculate the lift from a fixed wing by imposing 
      a cap on max achievable CL

      Inputs:
      1. L_req    : wing_lift (N)
      2. q        : dynamic pressure (N/sq.m)
      3. CLmax    : maximum lift coefficient

      Outputs:
      1. wing_lift: wing lift in Newtons
      """

      qS        = q*self.area
      cl        = L_req/qS
      if cl > CLmax:
         cl     = CLmax

      wing_lift = qS*cl       
      return wing_lift 

#====================================================================
# function calculate wing drag coefficient * wing plan-form area
#====================================================================

   def wing_f(self, cl, gamma):
      """
      function calculate CD*S of a wing 
      Input:
      1. cl    : operating lift coefficient
      2. gamma : flight path angle (radians)

      Output:
      1. fWing : Drag/dynamic pressure = CD*S (sq.m)
      """

      CDw       = self.cd0 + self.K*cl*cl + cl*gamma  # third term is due to climb, small angle assumption
      fWing     = CDw * self.area                
      return fWing 

#====================================================================
# function to calculate weight of fixed wings
#====================================================================

   def weight_estimate(self, Wt, tech_factors, redundancy, rotor):
      """
      this function calculates the weight of a fixed wing using three different methods
      [AFDD model, frequency-based spar sizing and strength-based]

      Inputs:
      1. Wt             : take-off weight (LBS)  --> not Newtons or kgs
      2. tech_factors   : technology factor to scale output masses 
      4. redundancy     : dictionary containing details of actuator redundancies 
      5. rotor          : class containing details of the rotor

      Outputs:
      1. wt             : dictionary containing breakdown of wing structural masses (kg)
      """
# initializations
      nwing       = self.nwings
      nrotors     = self.nrotors/(2*nwing)        # number of rotors along half span, each wing

      nwing       = self.nwings                   # number of wings
      Aw          = self.aspectratio              # aspect ratio 
      Sw          = self.area * m2f * m2f         # in sq.ft 
      fL          = self.lift_frac                # wing lift fraction for this group (0 to 1)
      W           = fL*Wt/nwing                   # thrust carried by each fixed wing (lbs)
            
# size motor mounts; tip mass at end of cantilever beam
# beam length = one rotor radius + 30% chord (where tilt axis, spar are located)
# tube radius = 15% of rotor radius
# target freq = specified in mmmass function (*m*otor *m*ount *mass*)

# only one rotor type allowed on this wing
      rgid        = self.rotor_group_id
      rg          = rotor.groups[rgid] 
      Mk          = rg.mass_assembly
      L_mount     = rg.radius + self.chord*0.3
      rtube       = 0.125*rg.radius*0.5

# size the motor mounts
      mount_mass  = mmmass(L_mount, Mk, rtube, f)

# remember mass of all motor mounts, store in lbs
      mounts_wt  = kg2lb*mount_mass*nrotors*(2*nwing)         

# Target frequency model: we have masses of rotors, motors and mounts
# calculate spar mass
      tskin       = 15e-4              # 3 layers, each 0.5mm; based on Alpha
      lbyc        = 2.1
      rho         = 1650.0             # mean density, kg/cu.m (fold in foam too)
      MbyA        = tskin*lbyc*rho

      c           = self.chord         # mean wing chord in meters
      rbar        = tau_w*c*0.5        # mean tube radius , meters   

      rroot       = 2*rbar/(1+taper)   # root spar radius
      L           = self.span*0.5      # beam length

# weight to match structural frequency of wing spar
      if(model == 'target_freq'):
         y           = self.ycoords 
         Mk          = self.masses + mount_mass 
         T           = self.rotor_T
         M_spar      = spar_mass(rroot, taper, L, MbyA, Mk, y, T, fn, tau_w)
         wt          = (MbyA*self.area + M_spar*2)*2.2

# wing weight, statistical method: AFDD model
      else:
         wt       = 5.66411*f_LGloc* (W*0.001/cos(Lamda))**0.847 * (nz**0.3958) * (Sw**0.21754)        \
                    *sqrt(Aw)* ((1.0+taper)/tau_w)**0.09359 #* (1.0 - b_fold)**(-0.14356)

# add up contributions from duplicate wings in this group
      wt          = wt * nwing           # weight of all wings in group, lbs
      wing_wt     = wt
      area        = Sw*nwing             # in sq.m
      self.wt     = wt*lb2kg/nwing       # weight per wing, kg

# for control surfaces, take total area and total weight; scaling law adapted by AS
# track weight in pounds
      actuator_wt = 0.01735*(Wt**0.6435)*(area**0.40952)

# weight effect of redundany actuators
      actuator_wt = actuator_wt * redundancy['wing_flap']

# tilt actuators: 7.5% of wing weight: track in lbs
      tilt_wt     = 0.075*wing_wt*redundancy['tilt_actuator']

# apply tech factor, take total weight of wing+actuators on wing
      wing_wt     = wing_wt     * tech_factors.wing
      actuator_wt = actuator_wt * tech_factors.flight_control
      tilt_wt     = tilt_wt     * tech_factors.flight_control

#convert to kg, return dict for wing structure, actuators and wires
      wt       = {'structure': wing_wt*lb2kg,'actuators':actuator_wt*lb2kg, 
                    'tilters': tilt_wt*lb2kg,  'mounts' :  mounts_wt*lb2kg}

      return wt

#====================================================================
# class that holds data for all wing groups and methods
#====================================================================

class wings:

   def __init__(self, data={}, nseg=0):
      
      """
      initialization function for all wing groups
      """

      self.ngroups      = 0
      self.nwings       = 0
      self.area         = 0.0
      ngrp              = 0
      self.groups       = {}
#loop over wing groups, find design parameters

      if(bool(data)):
         for key in sorted(data):
#remember # of wing groups
            self.groups[ngrp]    = wing_group(data[key], key, nseg)
            self.nwings          = self.nwings + self.groups[ngrp].nwings         
            ngrp                 = ngrp + 1
         self.ngroups      = ngrp

#====================================================================
# wing weight estimation function
#====================================================================

   def weight_rollup(self, GTOW, tech_factors, redundancy, rotor):
      """
      function to calculate weight of all fixed wings

      Input:
      1. GTOW           : gross take-off weight in lbs 
      2. tech_factors   : technology factor to scale output masses
      3. redundancy     : dictionary containing details of actuator redundancies
      4. rotor          : class containing details of the rotor

      Output:
      1. weights: dictionary with breakdown and total of all wing group weights
      """

      weights     = {}

# loop over wing groups, get weight breakdown for each group
      for igroup,group in self.groups.items():
         group_weight   = group.weight_estimate(GTOW, tech_factors, redundancy, rotor)

# accumulate total and add group number tag to dictionary
         for k,v in group_weight.items():
            k2             = 'wing_group'+str(igroup)+'_'+k
            weights[k2]    = v

# store total in dictionary, return
      weights['total']     = dict_accumulation(weights) 
      return weights 

#====================================================================
# wing sizing function
#====================================================================

   def cruise_sizing(self, W, segment, size_flag, dfus):

      """
      function to calculate cruise performance
      Inputs: 
          1. W        : vehicle weight (Newtons)
          2. segment  : mission segment information (class)
          3. size_flag: indicates whether to size the wing or not
                        for this cruise condition.
          4. dfus     : fuselage width in meters (equivalent diameter)
      
      Outputs:
          1. Wings_lift : total lift from all wings (Newtons)
          2. Wings_f    : CD*S of all wings (sq.m) = drag/dynamic pressure
      """

      Wings_lift           = 0.0
      Wings_f              = 0.0
      self.max_span        = 0.0
      Vcruise              = segment.cruisespeed*kts2mps   # in m/s
      Vclimb               = segment.rateofclimb/60.0      # in m/s
      q                    = 0.5*segment.rho*Vcruise*Vcruise
      gamma                = Vclimb/Vcruise                # approx flight path angle

      for i in range(self.ngroups):
         group             = self.groups[i]
         loadingFrac       = group.lift_frac
         nwings            = group.nwings 
         wing_lift         = W*group.lift_frac/nwings

# obtain span and mean chord from wing cl
# also check for stall speed and set upper limit
         V_ratio           = Vcruise/group.stall_speed
         CLmax             = group.CLmax*(V_ratio*V_ratio)

# perform wing sizing if required
         wing_cl           = min(group.cl, CLmax)
         if size_flag:
            group.size_wing(wing_lift, q, wing_cl, dfus)

# Non-sizing cruise flight condition
# Saturate wing cl at CLmax, calculate everything else
         else:
            wing_lift      = group.calculate_lift(wing_lift, q, CLmax)

# Compute drag coefficient of wing * area = equiv. flat plate area, sqm
         fWing             = group.wing_f(wing_cl, gamma)
         Wings_f           = Wings_f + fWing*nwings 
         Wings_lift        = Wings_lift + wing_lift*nwings

         self.max_span     = max(self.max_span, group.span)

      return Wings_lift, Wings_f
