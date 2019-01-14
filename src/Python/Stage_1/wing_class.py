import sys, numpy 
from numpy import sin, cos, log10, log2, sqrt, pi
from scipy.special import jv as besselj
sys.path.insert(0,'../Stage_0/')
from dict2obj import obj
from conversions import *
#====================================================================
# fixed wing information
#====================================================================

class wing_group:

#====================================================================
# initialization
#====================================================================

   def __init__(self, data, key, nseg):
      # self              = obj(data)
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

      return None

#====================================================================
# function to size the wing given the lift "lift (N)" and
# dynamic pressure "q" (N/sq.m) as well as operating lift coefficient
# "CLmax"
# also uses "df" --> fuselage diameter, meters to calculate Oswald eff
#====================================================================

   def size_wing(self, lift, q, cl, df):
      self.area     = lift/(q*cl)                        # area of each wing (left+right)
      self.span     = sqrt(self.aspectratio * self.area) # span of each wing (left+right)
      self.chord    = self.area / self.span

#====================================================================
# use P-Q combination from following reference to calculate oswald eff
# http://www.fzt.haw-hamburg.de/pers/Scholz/OPerA/OPerA_PUB_DLRK_12-09-10.pdf
# Equations (5), (6), (7) and (8) 
#====================================================================
      
#      print('original e is ',self.oswald)
      if(self.oswald == 0.0):
         u           = 0.99                  # from Kroo [2] of above ref.
         s           = df/self.span
         s           = 1.0 - 2.0*s*s         # eqn (6)
         Q           = 1.0/(u*s)             # eqn (5)
         P           = 0.38*0.02             # eqn (7), CD0=0.02
         self.oswald = 1.0/(Q + P*pi*self.aspectratio)
      self.K         = 1.0/(pi*self.aspectratio*self.oswald)
#      print('updated e is ',self.oswald)
      return None 

#====================================================================
# calculate actual wing lift "wing_lift (N)" from the reqd. lift 
# L_req (N), dynamic pressure "q" (N/sq.m) and max lift coefficient "CLmax"
#====================================================================

   def calculate_lift(self, L_req, q, CLmax):

      qS        = q*self.area
      cl        = L_req/qS
      if cl > CLmax:
         cl     = CLmax

      wing_lift = qS*cl       
      return wing_lift 

#====================================================================
# calculate wing equivalent flat plate area for a given operating 
# lift coefficient cl and flight path (small angle) gamma
# equivalent flat-plate area defined as CD*S
#====================================================================

   def wing_f(self, cl, gamma):
      CDw       = self.cd0 + self.K*cl*cl + cl*gamma
      fWing     = CDw * self.area                
      return fWing 

#====================================================================
# bigger class for all wings in system
#====================================================================

class wings:

#====================================================================
# initialization function: called once when launching the code
#====================================================================

   def __init__(self, data, nseg):
      self.ngroups      = 0
      self.nwings       = 0
      self.area         = 0.0
      ngrp              = 0
      self.groups       = {}
#loop over wing groups, find design parameters
      for key in sorted(data):
#remember # of wing groups
         self.groups[ngrp]    = wing_group(data[key], key, nseg)
         self.nwings          = self.nwings + self.groups[ngrp].nwings         
         ngrp                 = ngrp + 1
      self.ngroups      = ngrp


#====================================================================
# function to calculate cruise performance
# inputs are 
#     vehicle weight "W" (Newtons)
#     horizontal flight speed "Vcruise" (m/s) parallel to ground
#     vertical   flight speed "Vclimb"  (m/s) normal   to ground
#     dynamic pressure "q" (N/sq.m)
#     logical "size_flag" indicates whether to size the wing or not
#     for this cruise condition.
#     fuselage width "dfus" in meters --> equivalent diameter
#
# outputs are 
#     wing CD*S "Wings_f" (sq.m)
#     total lift of all wings "Wings_lift" (Newtons) 
#====================================================================

   def cruise_sizing(self, W, segment, size_flag, dfus):

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

#====================================================================
# obtain span and mean chord from wing cl
# also check for stall speed and set upper limit
#====================================================================

         V_ratio           = Vcruise/group.stall_speed
         CLmax             = group.CLmax*(V_ratio*V_ratio)

#====================================================================
# perform wing sizing if required
#====================================================================

         wing_cl           = min(group.cl, CLmax)
         if size_flag:
            group.size_wing(wing_lift, q, wing_cl, dfus)

#====================================================================
# Non-sizing cruise flight condition
# Saturate wing cl at CLmax, calculate everything else
#====================================================================

         else:
            wing_lift      = group.calculate_lift(wing_lift, q, CLmax)

#====================================================================
# Compute drag coefficient of wing * area = equiv. flat plate area, sqm
#====================================================================

         fWing             = group.wing_f(wing_cl, gamma)
         Wings_f           = Wings_f + fWing*nwings 
         Wings_lift        = Wings_lift + wing_lift*nwings

         self.max_span     = max(self.max_span, group.span)

      return Wings_lift, Wings_f
