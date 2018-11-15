#===============================================================================
# python function to use CF and adjust spar thickness
#===============================================================================
import numpy
def evaluate_CFv2(spar, skin, Rotor, contour_pr):

#===============================================================================
# Rotor properties
#===============================================================================

    chord       = Rotor['chord']             # mean chord, m 
    R           = Rotor['R']                 # rotor radius, m 
    Omega       = Rotor['Omega']             # rotor speed, rad/s 
    Mz          = 0.0                        #  lag bending moment, root ; remove with torque offset
    My          = Rotor['Fz']*0.75*R         # flap bending moment, root

#===============================================================================
# set average coning angle
#===============================================================================

    beta        = 2*numpy.pi/180.0              # average coning angle 
    zeta        = 2*numpy.pi/180.0              # lag angle, radians 
    e_lag       = 0.05*R                        # lag hinge, meters

#===============================================================================
# derived quantities
#===============================================================================

    Vtip        = Omega * R 
    c2          = chord*chord 
    c3          = chord*c2
    c4          = chord*c3

#===============================================================================
# spar properties
#===============================================================================

    t           = spar['t']*chord               # spar thickness, m
    b           = 0.10*chord                    # spar width, m 
    E_spar      = spar['E']
    G_spar      = spar['G']
    sigma_spar  = spar['sigma_y']               # stress limit in spar
    rho_spar    = spar['rho']                   # spar material density, kg/cu.m

    h           = 0.12*chord                    # spar height, m
    A_spar      = 2*(b+h)*t                     # spar area, sq.m
    m_spar      = A_spar*rho_spar               # spar mass/span, kg/m
    EA_spar     = E_spar*A_spar                 # axial stiffness of spar

    spar['m']   = m_spar

#===============================================================================
# Spar inertia/stiffness properties
#===============================================================================

    xNA         = spar['xNA']*chord
    xc          = spar['xc']*chord              # spar center location wrt neutral axis, +ve from LE

    Izz_spar    = 0.5*h*t*b*b + t*b*b*b/6.0     # 2nd area M.I., lag bending                    
    Iyy_spar    = 0.5*b*t*h*h + t*h*h*h/6.0     # 2nd area M.I., flap bending
    Izz_spar    = Izz_spar + h*t*t*t/6.0
    Iyy_spar    = Iyy_spar + b*t*t*t/6.0

    Izz_spar    = Izz_spar + A_spar*(xNA-xc)*(xNA-xc)
    Ixx_spar    = Izz_spar + Iyy_spar           # 2nd area M.I., torsion for C.S.
    I_spar      = Ixx_spar * rho_spar           # 2nd mass M.I., torsion, kg-m^2/m

#===============================================================================
# spar torsion properties: assume closed box
#===============================================================================

    Acs         = b*h 
    intds       = 2*(b+h)                       # perimeter of box
    J_spar      = 4*Acs*Acs/(intds)*t
    GJ_spar     = G_spar*J_spar

#===============================================================================
# Find bending stiffness of spar
#===============================================================================

    EIzz_spar   = E_spar*Izz_spar
    EIyy_spar   = E_spar*Iyy_spar

#===============================================================================
# skin properties: +/- 45 deg carbon
#===============================================================================

    rho_skin    = 1600.0                          # skin density, kg/cu.m ==> carbon fiber
    G_skin      = 33.0e9                          # carbin fiber, +/- 45 deg
    E_skin      = 17.0e9                          # carbon fiber, +/- 45 deg
    sigma_skin  = 47e6                           # yield stress in skin 

#===============================================================================
# skin properties: 0/90 carbon
#===============================================================================

#    rho_skin    =  1600.0                         # skin density, kg/cu.m ==> carbon fiber
#    G_skin      =   5.0e9                         # Pa
#    E_skin      =  70.0e9                         # Pa
#    sigma_skin  = 275.0e6                         # tensile yield stress, Pa 

#===============================================================================
# Skin dimensions
#===============================================================================

    ts          = skin['t']*chord                # skin thickness, m 
    A_skin      = contour_pr['A_skin']*ts*chord  # skin area, sq.m
    m_skin      = A_skin*rho_skin                # mass / span, skin (kg/m)
    Izz_skin    = contour_pr['Izzskin']*c3*ts    # 
    Iyy_skin    = contour_pr['Iyyskin']*c3*ts    # 
    I_skin      = (Izz_skin + Iyy_skin)*rho_skin # skin 2nd mass M.I, kg-m^2/m
    EA_skin     = E_skin*A_skin
    # GJ_skin     = G_skin*(Izz_skin + Iyy_skin)

    A_cs        = contour_pr['Atotal']*c2       # assume total area is covered in honeycomb
    ds          = contour_pr['dsbyt']*chord
    GJ_skin     = 4*A_cs*A_cs/(ds/ts)*G_skin

#    print 'python version'
#    print spar['t'],chord
#    print t, ts
#===============================================================================
# Find bending stiffness of skin (N-m^2)
#===============================================================================

    EIzz_skin   = 0.0#E_skin*Izz_skin
    EIyy_skin   = 0.0#E_skin*Iyy_skin

#===============================================================================
# FInd bending stress due to flap and lag bending loads
#===============================================================================

    EIzz        = EIzz_spar# + EIzz_skin
    EIyy        = EIyy_spar# + EIyy_skin
    EA          = EA_spar  #+ EA_skin

#===============================================================================
# calculate LE wt contributions to inertia
#===============================================================================

    m_LE        = skin['m_LE']                  # mass/span, kg/m
    I_LE        = m_LE*0.125*0.125*c2           # second mass MI about QC

#===============================================================================
# Honeycomb properties
#===============================================================================

    A_hc        = A_cs - b*h - A_skin           # subtract spar, skin area
    rho_hc      = 32.0                         # filler density, kg/cu.m
    m_hc        = rho_hc * A_hc                 # mass/span, honeycomb
    Ixx_hc      = contour_pr['mi2a_fill']*c4              # 2nd moment of area, entire CS
    Ixx_hc      = Ixx_hc - b*h*h*h/12.0 - h*b*b*b/12.0    # subtract out inertia from section replaced by spar
    I_hc        = Ixx_hc*rho_hc

#===============================================================================
# add leading edge protection weight: up to 25% chord on upper, lower surface
# thickness = 0.045"
#===============================================================================

    rho_LEP     = 8900.0                        # Nickel density, kg/cu.m
    A_LEP       = 0.2*2*chord                   # area per span of leading edge strip, sq.m/m 
    t_LEP       = 0.045*0.0254*0.5              # mean thickness of coating, m 
#    t_LEP       = 0.0
    m_LEP       = A_LEP*rho_LEP*t_LEP           # mass/span of Nickel plating, kg/m
    x_LEP       = 0.125*chord                   # CG of leading edge protection

#===============================================================================
# total mass per span, kg/m without LE but use it to calculate CF loads
#===============================================================================

    m_total     = m_hc + m_spar + m_skin + m_LEP
    CF          = (m_total+m_LE)*Vtip*Vtip*0.5
#    print(ts,t) 
#    print(m_hc,m_spar,m_skin,m_LEP,m_total)
#===============================================================================
# split CF into spar and skin
#===============================================================================

    CF_spar     = CF*(m_total - m_skin)/m_total
    CF_skin     = CF*m_skin/m_total

    sxx_sp_CF   = CF_spar/A_spar
    sxx_sk_CF   = CF_skin/A_skin

#===============================================================================
# let the flap bending load be carried with "beta" deg of coning (avg)
#===============================================================================

    beta_req    = My/(CF*R*0.5)
    if(beta >= beta_req):
       beta     = beta_req*0.99
    My          = My - CF*beta*R*0.5

#===============================================================================
# reduce lag bending moment due to blade elastic lag: model as approx rigid
#===============================================================================

    zeta_req    = Mz/(e_lag*0.5*CF)
    if(zeta >= zeta_req):
        zeta    = zeta_req*0.99 

    Mz          = Mz - zeta*e_lag*0.5*CF

#===============================================================================
# find bending moment in spar and skin
#===============================================================================

    My_spar     = My*EIyy_spar/EIyy 
    My_skin     = My*EIyy_skin/EIyy 

    Mz_spar     = Mz*EIzz_spar/EIzz 
    Mz_skin     = Mz*EIzz_skin/EIzz 

#===============================================================================
# find bending stresses
#===============================================================================

    sxx_sp_flap = My_spar*0.5*h/Iyy_spar
    sxx_sk_flap = My_skin*0.05*chord/Iyy_skin

    sxx_sp_lag  = Mz_spar*0.5*b/Izz_spar
    sxx_sk_lag  = Mz_skin*0.25*chord/Izz_skin       # was 0.05; must take at LE => 0.25 chord from NA

    sxx_sp      = sxx_sp_flap + sxx_sp_lag + sxx_sp_CF + 1.0
    sxx_sk      = sxx_sk_flap + sxx_sk_lag + sxx_sk_CF + 1.0    # desingularize

#    print sxx_sp_lag/sxx_sp

#===============================================================================
# Check if spar stays inside the airfoil, i.e. max value
# h/2 = centerline location
#   t = thickness of wall
# top of spar is at h/2 + t/2
#===============================================================================

    tmax        = 2*(0.12*chord - 0.5*h)    
    tmin        = 0.0005        # min thickness for manufacturing

#===============================================================================
# find spar thickness from min stress criterion
#===============================================================================

    spar_sr     = sxx_sp/sigma_spar
    t_spar1     = t 
    if(spar_sr > 1.0):                             # if spar fails, boost its size accordingly
        t_spar1     = sxx_sp/sigma_spar*t*1.05

    if(t_spar1 <= tmin):
        t_spar1 = tmin
    if(t_spar1 >= tmax):
        t_spar1 = tmax

#===============================================================================
# find skin safety factor; manipulate spar flexural stiffness to give the 
# correct allowable bending loads
#===============================================================================

    skin_sr     = sxx_sk/sigma_skin

    t_spar2     = t 
    if(skin_sr > 1.0):
        t_spar2 = t*1.05

    if(t_spar2 <= tmin):
        t_spar2 = tmin 
    if(t_spar2 >= tmax):
        t_spar2 = tmax 

#    print t, t_spar1, t_spar2
#    print spar_sr, skin_sr

    t_spar      = max(t_spar1, t_spar2)

    spar['t']   = t_spar/chord

#===============================================================================
# Total rotational inertia about QC, kg-m^2/m; this is I_theta in ENAE 633 notes
#===============================================================================

    I_total     = I_hc + I_skin + I_spar + I_LE

#===============================================================================
# find torsion frequency with spar + skin: nonrotating, then rotating
#===============================================================================

    GJ_total    = GJ_spar + GJ_skin
    nut         = numpy.pi*0.5*numpy.sqrt(GJ_total/(I_total*Vtip*Vtip))
    nut         = numpy.sqrt(nut*nut + 1.0) 

#===============================================================================
# Update target skin thickness based on torsion frequency
#===============================================================================

    dsbyt       = contour_pr['dsbyt']
    nut_tar     = 4.3
    param       = 4.0/(numpy.pi*numpy.pi)*(nut_tar*nut_tar-1)
    GJ_tar      = param*I_total*Vtip*Vtip
    J_tar       = (GJ_tar - GJ_spar)/G_skin
    ts1         = dsbyt*chord*J_tar/(4.0*A_cs*A_cs)
    if(ts1 < tmin):          # min manufacturable limit
       ts1      = tmin       

#===============================================================================
# find total torsion moment and skin thickness reqd. to resist pitching moment
#===============================================================================

    Cm       = 0.017
    M        = 1.225*Vtip*Vtip*chord*chord*R*Cm/6.0*0.9
    M_skin   = M#*GJ_skin/GJ_total
    # print(M_skin,M)
    ts2      = M_skin/(2*sigma_skin*A_cs)            # Torsion skin thickness
    if(ts2 < tmin):
        ts2  = tmin
    skin_SF  = min(ts/ts2, sigma_skin/sxx_sk,ts/ts1)
    print(ts,ts2,ts1)
    x1 = input('eg?')
    
    # skin_SF  = sigma_skin/sxx_sk

#===============================================================================
# Also find CG of CS: get first moment of area from filler
#===============================================================================

    xcg         = contour_pr['mi1a_fill']*c3            # first moment of area for honeycomb, cu.m 

#===============================================================================
# subtract first moment of area taken up by scoop-out for spar
#===============================================================================

    xcg         = xcg  - b*h*xc
    xcg         = xcg *rho_hc                          # first mass moment for honeycomb, kg-m/m

#===============================================================================
# next: take skin and find its CG
#===============================================================================

    mi1a_skin   = contour_pr['mi1a_skin']*ts*chord*chord    # m^3: first moment of area for skin

    xcg         = xcg + rho_skin*mi1a_skin      # add moment due to skin
    xcg         = xcg + m_spar*xc               # add moment due to spar

    xcg         = xcg + m_LEP*x_LEP             # add moment due to LE protection
    xcg         = xcg/(m_total*chord)           # cg location, nondiml by chord

#===============================================================================
# find NA location
#===============================================================================

    na          = (mi1a_skin*E_skin + xc*A_spar*E_spar)/(A_skin*E_skin + A_spar*E_spar) 

#===============================================================================
# set spar location
#===============================================================================

    spar['xc']  = 0.25 + (EA_skin*0.25 - E_skin*mi1a_skin/chord)/EA_spar

    if(spar['xc'] < 0.14):
        spar['xc'] = 0.14

#===============================================================================
# set leading edge weight if required:
#===============================================================================

    if (xcg <= 0.25):
        skin['m_LE'] = 0.0
    else:

#===============================================================================
# Set leading edge weight to bring CG to QC
#===============================================================================

        m_LE          = (4*m_total*xcg - m_total)
        skin['m_LE']  = m_LE                          # mass/span, LE wts

#===============================================================================
# Recalculate cg and update total mass per span of section
#===============================================================================

    m_LE        = 0.0
    xcg         = xcg*m_total/(m_total + m_LE)
    m_total     = m_total + m_LE

#===============================================================================
# Stress ratios
#===============================================================================

    spar_SF     =  sigma_spar/sxx_sp
#    skin_SF     = 2.0

#===============================================================================
# Return sectional properties
#===============================================================================

    section     = {'m': m_total, 'Io': I_total, 'xcg': xcg, 'm_LE': m_LE, 
                   'spar_SF': spar_SF, 'skin_SF': skin_SF ,
                   'na': na}

#===============================================================================
# End of operations
#===============================================================================

    return section
