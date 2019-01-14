import sys, numpy 
from numpy import sin, cos, log10, log2, sqrt, pi
from scipy.special import jv as besselj
sys.path.insert(0,'../Stage_0/')
from dict2obj import obj
from conversions import *

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
