#=============================================================================
# estimate flat-plate area of helicopter
# fuselage and pylon drag       ==> added 20% interference and 20% roughness
# landing gear                  ==> landing gear, no interference
# main rotor hubs               ==> condensed expression
# tail rotor                    ==> condensed expression
# engine nacelles               ==> need to improve (right now fudge factor in fuselage)
# horizontal, vertical tails    ==> 8%/4% interference
# misc. items:
#   cooling losses,             ==> scaled with installed power
#   protruberances,             ==> add 15% to final answer
#   interference,               ==> add 20% to each part if not explicitly specified
#   roughness                   ==> increase Cf by 20%
#
# Inputs
#   V = flight speed in m/s
# 
#=============================================================================
from numpy import pi, sqrt, log10, exp, log, cos
from conversions import kg2lb, lb2kg, f2m,kts2mps
class _find_f:
 def find_f(self, V):

#=============================================================================
# unpack inputs
#=============================================================================

    df        = self.emp_data.Geometry.fuselage_width 
    lf        = self.wing.max_span*0.5
    f_all     = {}
    Af        = pi*df*df                     # mean cross-section area, sq ft
    fine      = lf/df                        # fineness ratio
    W         = self.massTakeoff*2.2         # weight in lbs 

#=============================================================================
# fuselage and pylon drag
#=============================================================================

    nu        = 1.48e-5                     # kinematic viscosity, m^2/s
    a         = 340.0                       # speed of sound, m/s
    M         = V/a                         # cruise mach number    
    Re        = V*lf/nu                     # Reynolds number for fuselage length
    Cf        = 0.455/(log10(Re)**2.58 * (1+0.144*M*M)**0.65)       # skin friction coefficient

#calculate form factor
    FF        = 1.05 + 0.001*fine + 1.5/fine**1.5 + 8.4/fine**3

#interference factor and extra rivets
    Fint      = 1.0
    Swet      = lf*df*4 + df*df             # wetted area, sq.m
    f_fus     = Cf*FF*Fint*Swet             # fuselage flat-plate area, sq.m

    f_all['fus'] = f_fus

#=============================================================================
# contraction length
#=============================================================================

    if V <= 82:                             # for utility, make it less streamlined
        lc    = df*1.5
    else:
        lc    = 2.0*df                      # otherwise reduce base drag

#=============================================================================
#fuselage upsweep and contraction
#=============================================================================

    df_cont   = Af*0.008*(6*(df/lc)**2.5-1)     # extra flat-plate area

    f_all['base'] = df_cont

#=============================================================================
# Loop over rotor groups
#=============================================================================

    rotor           = self.rotor 
    ngroups         = rotor.ngroups
    f_all['spin']   = 0.0
    f_all['hubs']   = 0.0
    for i in range(ngroups):
        group           = rotor.groups[i]
        radius          = group.radius
        NR              = group.nrotors
        Nb              = group.nblade

#=============================================================================
# Tilting rotors
#=============================================================================

        if group.type == 'tilting':
            rspin           = 0.15*radius
            Sf              = pi*rspin*rspin           # reference area of spinner
            f_spin          = Sf*0.15                  # drag coefficient with spinner head
            f_all['spin']   = f_all['spin'] + f_spin*NR

#=============================================================================
# main rotor hub drag: combination of center section, pitch housing, 
# blade grip and root fittings; pitch links are there but small effect overall
#=============================================================================

        else:
            fhub          = radius*radius*(6e-4 + 0.0022225*Nb)
            fhub          = 1.2*fhub                             #20% interference
            f_all['hubs'] = f_all['hubs'] + fhub*NR

#=============================================================================
# momentum drag at higher speeds (> 160 knots)
# remove for electric motor
#=============================================================================

#    if V >= 160.0:
#        f_mom           = 0.3*nengine
#        f_all['mom']    = f_mom 
#    else:
    f_mom           = 0.0

#=============================================================================
#landing gear outside
#=============================================================================

    gear_type   = 'skids'
    if gear_type == 'outwheel':
        slope     = 0.443008 if W > 1000.0 else 0.0
        y0        = 0.57549
        logf      = y0 + slope*log(W*0.001)
        fLG       = exp(logf) 
    elif gear_type == 'skids':
        slope     = 0.49759 if W > 1000.0 else 0.0
        y0        = -0.11626
        logf      = y0 + slope*log(W*0.001)
        fLG       = exp(logf) 
    elif gear_type == 'retract':        # nothing to do
        y0        = -10.0
        slope     = 0.0
        fLG       = 0.0 
    else:
        quit('unknown landing gear type')

#=============================================================================
#assume the above equations hold up to 1000 lbs. 
#below 1000 lbs, linearly scale up with weight up to asymptote
#    if W < 1000:
#        fLG       = fLG * (W/1000.0)
#=============================================================================

    fLG           = fLG*0.3048*0.3048      # convert to sq.m
    f_all['LG']   = fLG#*0.65               # effect of fairings
    # print(W,fLG)
#=============================================================================
# horizontal/vertical stabilizers
#=============================================================================

    S_ht          = 0.001
    S_vt          = 1.2                    # Vtail area, sq.m
#    print('add vtail area for drag calculation')
#    S_vt          = 25.9 * (W / 1620)       #   vertical tail area, sq.ft

#calculations below
    AR            = 2.5                     # aspect ratio
    b_ht          = sqrt(S_ht*AR)           # span of horizontal tail, ft
    b_vt          = sqrt(S_vt*AR)           # span of vertical   tail, ft 
    c_ht          = b_ht/AR
    c_vt          = b_vt/AR

    lam           = 0.8                    # ctip/croot ; taper ratio
    tbyc          = 0.12                    # thickness to chord ratio of section
    tau           = 1.0                     # thickness ratio taper

    scale         = 2*(1+0.2*tbyc)         
    xt            = 0.3
    lam           = 0.0
    coslam        = cos(lam)
    FF            = (1+0.6/xt * tbyc + 100*tbyc**4)*(1.34*M**0.18*(coslam**0.28))

    Re_ht         = V*c_ht/nu 
    Re_vt         = V*c_vt/nu 
    Cf_ht         = 0.455/( (log10(Re_ht))**2.58 * (1+0.144*M**2)**0.65)
    Cf_vt         = 0.455/( (log10(Re_vt))**2.58 * (1+0.144*M**2)**0.65)

    f_ht          = Cf_ht*FF*1.15*S_ht*scale     #Interference Factor = 1.08
    f_vt          = Cf_vt*FF*1.15*S_vt*scale     #Interference Factor = 1.04

    Cdo_ht        = f_ht/S_ht                    # somewhere around 0.01
    Cdo_vt        = f_vt/S_vt

#=============================================================================
#induced drag part: assume tail is sized so as to lift 0.2*CLopt
#=============================================================================

    flam          = 0.119 - 0.0706*lam + 0.1659*lam*lam - 0.15*lam*lam*lam + 0.0524*lam*lam*lam*lam
    e_th          = 1.0/(1.0 + flam*AR)

#switch off fuselage interference for drag; we want worst case estimate

    kef_ht        = 1# - 2*(df/b_ht)**2
    kef_vt        = 1# - 2*(df/b_vt)**2      # fuselage presence correction

    kd0           = 0.804                   # viscous correction

    if M>0.3:
        keM       = 1 - 0.00152*(M/0.3-1)**(10.82)      # mach correction
    else:
        keM       = 1.0

    e_ht          = e_th*kef_ht*kd0*keM         # Oswald efficiency factors of HT/VT
    e_vt          = e_th*kef_vt*kd0*keM         # final answers

    K_ht          = 1.0/(pi*AR*e_ht)
    K_vt          = 1.0/(pi*AR*e_vt)

#steady-state: HT, VT operate at 20% of optimum 
    CL_ht         = sqrt(Cdo_ht/K_ht)*0.5#*(2/(2+AR*e_ht))
    CL_vt         = sqrt(Cdo_vt/K_vt)*0.5#*(2/(2+AR*e_ht))

#induced drag
    CDi_ht        = CL_ht*CL_ht*K_ht  
    CDi_vt        = CL_vt*CL_vt*K_vt

#increment horizontal tail and vertical tail flat-plate areas
    f_ht          = f_ht + CDi_ht*S_ht  
    f_vt          = f_vt + CDi_vt*S_vt

    f_all['ht']   = f_ht
    f_all['vt']   = f_vt

#=============================================================================
# find cooling drag 
# rationale is that air intakes are sized by hover requirement, so larger 
# areas reqd for intakes. use simple eqn from Stepnewski & Keyes (NASA) 
#=============================================================================

#    kc            = 4.0                 # cooling system design scaling factor
#    dfe           = 2.5e-5*Pins*kc      # additional flat-plate area due to cooling losses

#not needed for electric motors; air cooled by free stream or downwash
    dfe             = 0.0
    f_all['cool']   = dfe

#====================================================================
# Add cylinder mast drag for coaxial
# Assume CD = 0.2 for mast fairing
# Frontal area from 0.2 R spacing, cyl. diameter of 0.075 R
#====================================================================

#    if (False):
#        delta_f     = 0.2*0.075*0.2*rotor.radius*rotor.radius

#=============================================================================
# add up everything, throw in another 10% for protruberances
#=============================================================================

    total           = f_all['fus'] + f_all['base'] + f_all['spin'] + f_all['hubs'] + \
                      f_all['ht']  + f_all['vt']   + f_all['LG']   + f_all['cool']

#for high speed design, remove protrusions
    f_all['prot']   = 0.1*total

    total           = total + f_all['prot'] 
    f_all['total']  = total
    return f_all, Af*lf


#====================================================================
# wrapper function to calculate parasitic drag of the airframe
# Inputs
#   Vcruise    : Cruise speed, m/s
#   update     : True/False flag to indicate whether to estimate 
#                parasitic drag or not
# Outputs:
#   fParasitic : equivalent flat-plate area, sq.m
#   Vol        : fuselage volume, cu.m #
#====================================================================

 def flat_plate_wrapper(self, segment, update):


    aircraft            = self.all_dict['aircraft']
    flatPlateFactor     = aircraft['fdrag'] 
    add_f               = segment.add_f 
    Vcruise             = segment.cruisespeed*kts2mps

#====================================================================
# Calculate from parametric equation (modified equation by AS)
#====================================================================

    if update:
        fParasitic      = flatPlateFactor*(self.massTakeoff*kg2lb*1.e-3)**0.5e0*f2m*f2m
      
#====================================================================
# estimate flat-plate area from component build-up 
#====================================================================

        if flatPlateFactor <= 1.e-3:
            f_all, Vol     = self.find_f(Vcruise)
            self.Volume    = Vol
            f_all['ext']   = add_f                  # flat-plate area due to ext. stores
            f_all['total'] = f_all['total'] + f_all['ext']
            aircraft['f_breakdown']    = f_all
            fParasitic                 = f_all['total']

#====================================================================
# update flat plate factor and flat plate area (remember in class)
#====================================================================

        self.f_plate                  = fParasitic

    return None