import numpy
import copy 
import hydra
from spar_properties import spar_properties
from numpy import pi
#===============================================================================
# Physics based model to predict weight of rotor group
#===============================================================================

def blade_wt_modelv2(Rotor, rotor_tech_factor):

#===============================================================================
# spar properties: set based on material name
#===============================================================================

  spar        = spar_properties(Rotor['mat'].lower())

#===============================================================================
# Airfoil section properties: see test_airfoil.py
# Note: values hardcoded because they are functions of airfoil section only!
# we use 0012 for sizing
#===============================================================================

#  contour_pr  = numpy.load('naca0012.npy').item()
  contour_pr = {'Izzskin': 0.2952,   'Iyyskin': 0.00398,    'A_skin': 2.0392, 
                 'Atotal': 0.0822, 'mi1a_fill': 0.03456, 'mi1a_skin': 1.0052, 
                  'dsbyt': 2.0392, 'mi2a_fill': 0.006993}

#===============================================================================
# adjust properties iteratively: change spar thickness
#===============================================================================

  dts         = 0.0005          # minimum thickness for manufacturing
  tbyc        = 0.12            # blade thickness to chord ratio
  xRoot       = 0.1             # nondiml span length of root fitting
  Nstations   = 5               # number of spanwise stations
  betap       = 3.0             # precone angle, deg
  betap       = betap*pi/180.0  # precone angle, rad

#===============================================================================
# check if the spar material has been initialized
#===============================================================================

  if(hydra.crosssectiondata.material_init == 0):
    hydra.crosssectiondata.init_materials()

  chord       = Rotor['chord']
  R           = Rotor['R']
  Omega       = Rotor['Omega']
  nz          = Rotor['nz']
  Fz          = Rotor['Fz']

  r           = numpy.linspace(xRoot,1.0,Nstations)
  b           = 0.14*chord 
  h           = tbyc*chord
  Vtip        = Omega*R 

  c2          = chord*chord
  c3          = chord*c2
  c4          = chord*c3

#===============================================================================
# calculate masses of fixed components: leading edge protection 
#===============================================================================

  rho_LEP     = 8900.0                        # Nickel density, kg/cu.m
  l_LEP       = 0.34*chord                    
  A_LEP       = 2*l_LEP                       # area per span of leading edge strip, sq.m/m 
  t_LEP       = 0.045*0.0254*0.5              # mean thickness of coating, m, 0.045"
  m_LEP       = A_LEP*rho_LEP*t_LEP           # mass/span of Nickel plating, kg/m
  x_LEP       = 0.1*chord                     # CG of leading edge protection
  dx          = 0.25*chord - x_LEP
  I_LEP       = m_LEP*dx*dx                   # rotational inertia of LEP, kgm^2/m

#===============================================================================
# prepare constants for skin sizing
#===============================================================================

  Acs         = contour_pr['Atotal']*c2       # cross-section enclosed area, sq.m
  tmin        = 5e-4                          # 5 mm minimum
  sigma_skin  = 47e6                          # yield stress in skin 
  rho_skin    = 1660.0                        # skin material density, kg/m^3
  G_skin      = 33.0e9                        # skin shear modulus, Pa

#===============================================================================
# total torsion moment at root and reference skin thickness (shear stress based)
#===============================================================================

  Cm          = 0.02                                # PM coefficient
  c1          = 1.2256*Vtip*Vtip/6.0*c2*Cm*R
  ts_ref      = max(nz*c1/(2.0*Acs*sigma_skin), tmin)

#===============================================================================
# skin CG@ 49% chord
#===============================================================================

  l_skin      = contour_pr['A_skin']*chord
  A_skin      = l_skin*ts_ref
  m_skin      = A_skin*rho_skin         # skin mass per unit span
  x_skin      = 0.49*chord 

  Izz_skin    = contour_pr['Izzskin']*c3*ts_ref   # 
  Iyy_skin    = contour_pr['Iyyskin']*c3*ts_ref   # 
  I_skin      = (Izz_skin + Iyy_skin)*rho_skin    # skin 2nd mass M.I, kg-m^2/m

#===============================================================================
# calculate paint thickness and mass; treat it similar to the airfoil skin
#===============================================================================

  t_paint     = 0.00015               # paint
  rho_paint   = 1800
  t_ratio     = t_paint/ts_ref
  rho_ratio   = rho_paint/rho_skin

  m_paint     = m_skin*t_ratio*rho_ratio
  x_paint     = x_skin
  I_paint     = I_skin*rho_ratio*t_ratio

#===============================================================================
# calculate glue thickness and mass; treat it similar to the airfoil skin
#===============================================================================

  t_glue      = 2.54e-4               # glue thickness
  rho_glue    = 1800
  t_ratio     = t_glue/ts_ref
  rho_ratio   = rho_glue/rho_skin

  m_glue      = t_glue*rho_glue*(l_skin + 0.2*chord)
  x_glue      = x_skin 
  I_glue      = I_skin*rho_ratio*t_ratio 

#===============================================================================
# calculate filler (foam/HC) mass per span, CG at 57% chord
#===============================================================================

  rho_hc      = 52.0
  A_hc        = chord*chord*tbyc*0.4905
  m_hc        = rho_hc*A_hc             # mass per unit span, honeycomb
  x_hc        = 0.57*chord

  Ixx_hc      = contour_pr['mi2a_fill']*c4              # 2nd moment of area, entire CS
  Ixx_hc      = Ixx_hc - b*h*h*h/12.0 - h*b*b*b/12.0    # subtract out inertia from section replaced by spar
  I_hc        = Ixx_hc*rho_hc                           # rotational inertia of filler, kgm^2/m

#===============================================================================
# calculate leading edge mass per unit span to move CG to QC
#===============================================================================

  m_total     = m_LEP + m_skin + m_hc + m_paint + m_glue
  x_CG        = (m_LEP*x_LEP + m_skin*x_skin + m_hc*x_hc + m_paint*x_paint + m_glue*x_glue)/m_total

  x_LE        = 0.0
  m_LE        = m_total*(4.0*x_CG/chord - 1.0)
  I_LE        = m_LE*0.25*0.25*c2                 # 2nd mass MI about QC for LE weight

#===============================================================================
# section aggregated properties without spar
#===============================================================================

  x_CG        = (m_LE*x_LE + m_total*x_CG)/(m_LE + m_total)
  m_total     = m_total + m_LE 
  I_total     = I_hc + I_skin + I_LEP + I_LE + I_paint + I_glue

#===============================================================================
# calculate torsion natural frequency
#===============================================================================

  ds          = contour_pr['dsbyt']*chord
  J_skin      = 4*Acs*Acs*ts_ref/(ds)             # closed section polar MI, kgm^2/m
  GJ_skin     = G_skin*J_skin
  nutsq       = 1.0 + pi*pi*0.25*GJ_skin/(I_total*Vtip*Vtip)
  nu_t        = numpy.sqrt(nutsq) 
  nut_tar     = 3.3       
  if(nu_t < nut_tar):
    print('warning:torsion frequency is not matched!!',nu_t)

#===============================================================================
# calculate thickness reqd to match torsion freq. target
#===============================================================================

  # k           = I_skin/ts_ref
  # a0          = (nut_tar*nut_tar-1)*4*Vtip*Vtip/(pi*pi)
  # I_other     = I_total - I_skin
  # treq        = I_other*a0/(4*Acs*Acs/ds*G_skin - k*a0)

#===============================================================================
# loop over sections from outboard to inboard
#===============================================================================
  
  CF_out      = 0.0           # CF on outboard end
  a1          = m_total*Vtip*Vtip*0.5
  a0          = spar['rho']*(b+h)*Vtip*Vtip
  b1          = Vtip*Vtip*0.5
  M_out       = 0.0
  spar_mass   = 0.0

  for i in range(Nstations-1):

    xout      = r[Nstations-i-1]
    xin       = r[Nstations-i-2]
    dR        = (xout-xin)*R
    c0        = b1*(xout*xout-xin*xin)

#===============================================================================
# calculate centrifugal force, vertical shear and BM at this station
#===============================================================================

#===============================================================================
# CF has 3 parts: CF at outboard end, CF due to non-spar masses and 
# CF due to spar; first 2 parts are known, 3rd should be calc.
#===============================================================================

    CFo       = CF_out + a1*(xout*xout-xin*xin) 

#===============================================================================
# calculate flap bending moment from assumed quadratic variation of thrust
#===============================================================================

    FBM       = Fz*R*0.25*(xin*xin*xin*xin - 4*xin + 3)

#===============================================================================
# calculate RHS
#===============================================================================

    dRHS      = (M_out - CF_out*dR*betap - m_total*dR*betap*c0*0.5)/(h*b)
    dLHS      = c0*dR*betap*spar['rho']*(b+h)/(b*h)
    RHS       = CFo*0.5/(b+h) + FBM/(b*h) + dRHS

#===============================================================================
# calculate coefficients and solve for spar thickness
#===============================================================================

    LHS       = spar['sigma_y']/nz - a0*(xout*xout-xin*xin)*0.5/(b+h) + dLHS
    t         = max(RHS/LHS,tmin)

#===============================================================================
# calculate shear web thickness required to carry vertical load
# increase to box beam thickness if required
#===============================================================================
    
    Vz        = (1.0 - xin*xin*xin)*Fz 
    tweb      = 1.5*Vz/(h*sigma_skin/nz)

    if(t < tweb):
      t       = tweb 

#===============================================================================
# calculate section thickness, mass and CF at root end
# aggregate spar masses
#===============================================================================

    A_spar    = 2*(b*t+h*t)
    m         = A_spar*spar['rho']       # spar mass per unit span, this segment
    spar_mass = spar_mass + m*dR
    M_out     = M_out     - CF_out*dR*betap - (m+m_total)*c0*dR*betap*0.5
    CF_out    = CFo       + m*c0
    
#===============================================================================
# calculate axial stresses from bending and CF
#===============================================================================

    sigma1    = CF_out/A_spar
    Iyy       = 2*b*t*0.25*h*h

    eff_FBM   = FBM + M_out
    sigma2    = eff_FBM*0.5*h/Iyy

#===============================================================================
# prints
#===============================================================================

    # print('inboard %R = ',round(xin,2),end=' ')
    # print('axial stress ratio due to CF. = ',round(nz*sigma1/spar['sigma_y']*100,2), 
    #                         ' due to FBM = ',round(nz*sigma2/spar['sigma_y']*100,2))
    # print('spar thickness = ',round(t/chord*100,4),' % chord, or ',round(t*1000,2),' mm')
    # print('spar mass per unit span is ',round(m,2),' kg/m')
    # quit()
  m_total     = m_total*R + spar_mass

#===============================================================================
# calculate root fitting mass
#===============================================================================

  sigma_al    = 20e6                # Al 7075, fully reversing at 10^8 cycles w/ knockdowns
  rho_al      = 2800                # density, kg/cu.m
  rRoot       = 4*h                 # root fitting is circle of radius = 4*spar height

# total axial stress = due to CF + due to bending moment
# that "3" shouldn't be there, not sure why its featured
  factor      = 3.0
  RHS         = CF_out/(2.0*pi*rRoot*sigma_al) + eff_FBM/(factor*pi*rRoot*rRoot*sigma_al)
  LHS         = 1.0 - rho_al/sigma_al*Omega*Omega*(xRoot*R)*0.5
  t_Root      = RHS/LHS 
  m_Root      = 2*pi*rRoot*t_Root*(xRoot*R)*rho_al 
  m_total     = m_total + m_Root

#===============================================================================
# add a percentage for misc weights
#===============================================================================

  m_total     = m_total * 1.2

#===============================================================================
# Find blade mass in kg from wt/blade; convert to lb
#===============================================================================

  mass_blade  =  m_total*Rotor['nb']*Rotor['nr']  # blade   mass, kg 

#===============================================================================
# VPF actuator mass: scale with blade area, Omega^2
#===============================================================================

  mass_act    = 1.47*(Rotor['Omega']/228)**2*(Rotor['R']/0.75)*(Rotor['chord']/0.082)**0.5
  mass_act    = mass_act*Rotor['nr']

  mass_hub    = 4.84*Rotor['nr']*(Rotor['R']/0.75)**0.5        # scales with radius

#====================================================================
# Apply tech factor for rotor group
#====================================================================

  mass_blade  = mass_blade*rotor_tech_factor 
  mass_hub    = mass_hub  *rotor_tech_factor
  mass_act    = mass_act  *rotor_tech_factor

#===============================================================================
# return dictionary
#===============================================================================

  mass_rotor  = {'blades'  : mass_blade,
                 'hub'     : mass_hub,
                 'actuator': mass_act}

#===============================================================================
# End of operations
#===============================================================================
  
  return mass_rotor
