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
#=============================================================================
from numpy import pi, sqrt, log10, exp, log, cos

def find_f(inputs):

#=============================================================================
# unpack inputs
#=============================================================================

    W         = inputs['W']                 # weight, lbs
    R         = inputs['R']                 # rotor radius, ft  
    Vcruise   = inputs['V']                 # cruise speed, kts
    Pins      = inputs['P']                 # installed power, hp
    Nb        = inputs['Nb']                # number of blades
    f_all     = {}

#=============================================================================
# revert to defaults if some details not given 
#=============================================================================

    try:
        lf      = inputs['lf']
    except:
        lf      = 20.e0                      # fuselage length, ft

    try:
        NR      = inputs['NR']
    except:
        NR      = 1

    try:
        config  = inputs['config']
    except:
        config  = 1                     # single main rotor

    try:
        nengine = inputs['nengine']
    except:
        nengine = 1

    try:
        N_at    = inputs['Nprop']
    except:
        N_at    = 0 

#=============================================================================
# fuselage and pylon drag
#=============================================================================

    nu        = 1.72e-4                     # kinematic viscosity, ft^2/s
    a         = 1115.0                      # speed of sound, ft/s
    V         = Vcruise*1.689               # ft/s
    M         = V/a                         # cruise mach number    
    rho_pack  = 0.40                        # packing density: 20% of aluminum
    Af        = W/(7.2)/32.2/rho_pack       # equivalent cross-section area

#    Af        = 26.5/16000.0*W             # frontal area, sq.ftL scale it with uh60
    df        = 2.0*sqrt(Af/pi)*1.3         # equiv. diameter,, ft
    fine      = lf/df                       # fineness ratio
    Re        = V*lf/nu                     # Reynolds number for fuselage length
    Cf        = 0.455/(log10(Re)**2.58 * (1+0.144*M*M)**0.65)

#calculate form factor
    FF        = 1.05 + 0.001*fine + 1.5/fine**1.5 + 8.4/fine**3

#interference factor and extra rivets
    Fint      = 1.2

    Swet      = lf*df*4 + df*df

    f_fus     = Cf*FF*Fint*Swet     

    f_all['fus'] = f_fus

#=============================================================================
# contraction length
#=============================================================================

    try:
        lc    = inputs['lc']                    # length of contraction region, ft
    except:
        if V <= 160:                            # for utility, make it less streamlined
            lc    = df*0.6
        else:
            lc    = 2.0*df                      # otherwise reduce base drag

#=============================================================================
#fuselage upsweep and contraction
#=============================================================================

    df_cont   = Af*0.008*(6*(df/lc)**2.5-1)     # extra flat-plate area

    f_all['base'] = df_cont

#=============================================================================
# tilt-rotors, tilt-wings or tail-sitters
#=============================================================================

    if config == 2 or config == 5:              # need to account for pylon drag  
                                                # combine with spinner
#single engine: spinner only covers drive system
        rspin           = 0.12
        Sf              = rspin*R
        Sf              = pi*Sf*Sf*NR              # frontal area
        f_nac           = Sf*0.05                  # drag coefficient with spinner head
        f_all['nac']    = f_nac
#        print(NR,f_nac/NR/0.05/0.3048**2)
    else:

#=============================================================================
# auxiliary thrusters for smr configuration
#=============================================================================

        R_at            = 0.2*R                 # thruster radius
        rspin           = 0.2                   # spinner to rotor radius ratio
        Sf              = pi*(rspin*R_at)**2    # frontal area, sq.ft
        CD              = 0.2                   # drag coefficient of spinner
        f_nac           = Sf*CD*N_at            # spinner flat-plate area
        f_nac           = f_nac * 1.2           # add 20% for interference
        f_all['spin']   = f_nac                 # book-keep it in nacelle area

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

    try:
        gear_type   = inputs['gear_type']
    except:
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

#    print(fLG, gear_type )
#assume the above equations hold up to 1000 lbs. 
#below 1000 lbs, linearly scale up with weight up to asymptote
    if W < 1000:
        fLG       = fLG * (W/1000.0)

    fLG           = fLG*0.62
    f_all['LG']   = fLG

#=============================================================================
# horizontal/vertical stabilizers
#=============================================================================

    try:
        S_ht      = inputs['S_ht']          # h tail area, sq ft 
    except:
#        S_ht      = 0.01 * (W / 16000)      # horizontal tail area
        S_ht      = 0.0

    try:
        S_vt      = inputs['S_vt']
    except:
        S_vt      = 25.9 * (W / 1620)       #   vertical tail area, sq.ft

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
#=============================================================================
# main rotor hub drag
# combination of center section, pitch housing, blade grip and root fittings
# pitch links are there but small effect overall
#=============================================================================

    if config == 1 or config == 3:                      # edgewise configuration
        fhub          = 6e-4*R**2 + 0.0022225*Nb*R**2 
        fhub          = 1.2*fhub
        f_all['hub']  = fhub                        # 20% interference
    else:
        fhub          = 0.0

#=============================================================================
# tail rotor for smr configurations: take R = 0.2 R_MR and Nb = 4
#=============================================================================

    if config == 1:
        Nb            = 4
        fTR           = 0.005*(0.2*R)**2*(1+2*Nb)
        fTR           = 1.2*fTR
        f_all['TR']   = fTR
    else:
        fTR           = 0.0

#====================================================================
# Add cylinder mast drag for coaxial
# Assume CD = 0.2 for mast fairing
# Frontal area from 0.2 R spacing, cyl. diameter of 0.075 R
#====================================================================

    quit('OK?')
    if (aircraftID == 3):
        delta_f     = 0.2*0.075*0.2*rotor.radius*rotor.radius

#=============================================================================
# add up everything, throw in another 10% for protruberances
#=============================================================================

    total           = f_ht + f_vt + fLG + f_fus + fhub + fTR + f_nac + f_mom + \
                      df_cont + dfe 

#    print(f_ht, f_vt, fLG, f_fus, fhub, fTR, f_nac, f_mom, df_cont, dfe)

#for high speed design, remove protrusions
#    if V < 160:
    f_all['prot']   = 0.05*total
    total           = total + f_all['prot'] 

    f_all['total']  = total
    return f_all, Af*df