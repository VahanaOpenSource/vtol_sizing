from numpy import log10, pi
from conversions import *
#====================================================================
# FUSELAGE GROUP
#
# airframe, pressurization, crashworthiness
# 
# ALL UNITS IN IMPERIAL
#====================================================================

ng      = 3.8   # Max g lift
nl      = 3.5   # Landing load factor
sf      = 1.5   # Safety factor

uni_rho     = 1660      # Uni density [kg/m^3]
uni_stress  = 450e6     # Uni ultimate stress [Pa]
bid_rho     = 1660      # BID density [kg/m^3]
bid_stress  = 275e6     # BID ultimate compressive stress [Pa]
bid_shear   = 47e6      # BID ultimate shear stress [Pa]
bid_minThk  = 0.00046   # Minimum gauge, 3 ply
bid_bearing = 400e6     # Bearing alllowable
core_minThk = 0.0064    # Core thickness
core_rho    = 52        # Core density
glue_thk    = 2.54e-4
glue_rho    = 1800
rib_thk     = 0.0015
rib_width   = 0.0254
rib_rho     = 2768
paint_thk   = 0.00015
paint_rho   = 1800
canopy_thk  = 0.003175
canopy_rho  = 1180
steel_shear = 500e6 #Bolt shear strength

arealWeight = bid_minThk*bid_rho+core_minThk*core_rho+paint_thk*paint_rho

def fuselage_v2(vehicle_parameters):

#====================================================================
# unpack inputs
#====================================================================
   
  aircraftID = vehicle_parameters['aircraftID']      
  gtow       = vehicle_parameters['gtow']        # lbs
  nz         = 3.8                               # loa factor
  fac        = vehicle_parameters['tech_factors'].fuselage
  wing       = vehicle_parameters['wing']
  l_fus      = vehicle_parameters['l_fus']       # ft 
  h_fus      = vehicle_parameters['b_fus']       # ft 
  b_fus      = h_fus                             # ft 

#====================================================================
# get max span for all wing groups
#====================================================================

  span   = 0.0
  for gname,group in wing.groups.items():
    span     = max(span, group.span)

#====================================================================
# GB/ZL fuselage model
#====================================================================
   
  length  = span*0.5
#  length  = l_fus*0.3048             # fuselage length in meters
#  width   = b_fus*0.3048             # width of fuselage in meters
#  height  = b_fus*0.3048             # height in meters
  weight  = gtow/2.2*9.81            # vehicle weight in Newtons

  width   = 1.55
  height  = 1.55

#====================================================================
# Wetted area and skin mass
#====================================================================

  Swet    = 4*pi*(((length*width/4)**1.6 + (length*height/4)**1.6 + (width*height/4)**1.6)/3)**(1/1.6)
  m_skin  = Swet*arealWeight

#====================================================================
# Bulkhead Mass
#====================================================================
  
  bulkheadMass = 4*pi*height*width/4*arealWeight

#====================================================================
# Canopy Mass
#====================================================================

  canopyMass = Swet/10*canopy_thk*canopy_rho

#====================================================================
# Keel Mass due to lift
#====================================================================

  L           = ng*weight*sf        # Lift
  M           = L/4*length*2/3      # Peak moment from tail lift
  beamWidth   = width/1.4142        # Keel width
  beamHeight  = height/1.4142       # Keel height
  A           = M*beamHeight/(4*uni_stress*(beamHeight/2)**2)
  massKeel    = A*length*uni_rho
  
#====================================================================
#Keel Mass due to torsion
#====================================================================

  M           = 0.25*L*span/2       #Wing torsion 
  A           = beamHeight*beamWidth
  t           = 0.5*M/(bid_shear*A)
  massKeel    = massKeel+2*(beamHeight+beamWidth)*t*bid_rho

#====================================================================
#Keel Mass due to landing
#====================================================================

  F           = sf*weight*nl*0.6403   # Landing force, side landing
  A           = F/steel_shear         # Requried bolt area
  d           = 2*sqrt(A/pi)          # Bolt diameter
  t           = F/(d*bid_bearing)     # Laminate thickness
  V           = pi*(20*t)**2*t/3      # Pad up volume
  massKeel    = massKeel + 4*V*bid_rho #Mass of all 4 pad ups

#====================================================================
# Apply tech factors 
#====================================================================

  fac          = fac*1.2                # add 20% for fasteners, etc

  mass_skin    = m_skin       * fac
  mass_keel    = massKeel     * fac
  mass_bulkh   = bulkheadMass * fac
  mass_canopy  = canopyMass   * fac

  total      =  mass_skin + mass_keel + mass_bulkh + mass_canopy 

#====================================================================
# store kg mass in a dictionary and return the dictionary
#====================================================================

  fuselage   = {       'skin': mass_skin, 'keel': mass_keel, 
                   'bulkhead': mass_bulkh,  'canopy' : mass_canopy,
                    'total'  : total     }
  # print('link fuselage weight model to wing lift distribution - tail lift!')
  return fuselage