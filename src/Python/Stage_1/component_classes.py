import sys, numpy 
from numpy import sin, cos, log10, log2, sqrt, pi
from scipy.special import jv as besselj
sys.path.insert(0,'../Stage_0/')
from dict2obj import obj

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
#====================================================================

   def size_wing(self, lift, q, cl):
      self.area     = lift/(q*cl)                        # area of each wing (left+right)
      self.span     = sqrt(self.aspectratio * self.area) # span of each wing (left+right)
      self.chord    = self.area / self.span

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
# lift coefficient cl
# equivalent flat-plate area defined as CD*S
#====================================================================

   def wing_f(self, cl):
      CDw       = self.cd0 + self.K*cl*cl
      fWing     = CDw * self.area                
      return fWing 

#====================================================================
# bigger class for all wings in system
#====================================================================

class wings:

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
# aux thruster information
#====================================================================

class prop_sizing:

   def __init__(self, data):
      self.num          = data['npropeller']
      self.thrust       = 0.e0
      self.power        = 0.e0

#====================================================================
# aircraft properties (distributed to rotor, wing, etc later)
#====================================================================

class aircraft:
   def __init__(self, adict):

#====================================================================
# big-picture items
#====================================================================

      self.aircraftid      = int(adict['aircraftID'])

#====================================================================
# Mission
#====================================================================

      self.masscrew        =     adict['mass_crew']
      self.common_equip_wt =     adict['mass_common_equip']
#      self.fuselagewidth   =     adict['fuselage_width']
#      self.clearance       =     adict['clearance']

#====================================================================
# rotor
#====================================================================

      self.nrotor          = int(adict['nrotor'])
      self.roffset         =     adict['rotor_offset']
      self.rotorshafttilt  =     adict['rotor_shaft_tilt']

#====================================================================
# propeller
#====================================================================

      self.npropeller      =     adict['npropeller']

#====================================================================
# Engine
#====================================================================

      self.nengine         = int(adict['nengine'])

#====================================================================
# physical constants
#====================================================================

class constants:
   def __init__(self):
      self.grav    = 9.80665e0    # gravitation constant (m/s2)
      self.rho0    = 1.2256e0     # standard air density (kg/m3)
      self.f2m     = 0.3048e0     # conversion from feet to meter
      self.m2f     = 3.28084e0    # conversion from meter to feet
      self.hp2kw   = 0.7457e0     # conversion from hp to KW
      self.kw2hp   = 1.3410e0     # conversion from KW to hp
      self.kts2mps = 0.5144e0     # conversion from knots to m/s
      self.mps2kts = 1.94401e0    # conversion from m/s to knots
      self.lb2kg   = 0.45359e0    # conversion from lb to kg
      self.kg2lb   = 2.20462e0    # conversion from kg to lb
      self.mps2kph = 3.6e0        # conversion from m/s to km/hr
      self.hr2sec  = 3600e0       # conversion from hour to sec
      self.sec2hr  = 2.7777e-4    # conversion from sec to hour
      self.rpm2rps = 0.404719e0   # conversion from RPM to rad/sec
      self.in2m    = 0.0254e0     # conversion from inch to m
      self.m2in    = 3.93700e+1   # conversion from m to inch
      self.pi      = 3.14159265e0 # conversion from m to inch
      self.min2sec = 60.0e0       # conversion from mins to seconds
      self.sec2min = 1.6666666e-2 # conversion from secs to mins
      self.nm2m    = 1852.e0      # conversion from nautical miles to meters
      self.km2nm   = 0.539957e0   # conversion from kms to nautical miles
