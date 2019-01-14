import numpy 
from taper_beam import taper_beam
from numpy import log, pi 
E           = 122.0e9           # Young's modulus, Pa ==> stiffness based sizing
rho         = 1650.0            # density, kg/cu.m
sigma_max   =  275.0e6          # UNIT axial limiting stress, Pa ==> for strength based sizing
tau_max     =   47.0e6          # BID ultimate shear stress [Pa] ==> for strength based sizing

nz          = 3.8 
g           = 9.81              # accln due to gravity, m/s^2
tmin        = 10e-4             # minimum gauge, m ==> 2 layers of carbon
def motor_mount_mass(L, M, r, f):

    """ 
    this function calculates the mass of a cantilever beam with a tip mass 
    by targeting a natural frequency for the structure. Input parameters are 

    1. L    = length of beam, meters
    2. M    = tip mass on beam, kg 
    3. r    = outer radius of circular cross-section, meters 
    4. f    = target natural frequency, Hz

    Based on methodology documented in the following link
    https://drive.google.com/drive/folders/1MY-0ERXqkKBujEX6o23p2v-zkhR30n9B    
    """

    wn      = f*2.0*numpy.pi        # rad/s, natural frequency

    LHS     = 1.0                   # special case solution for uniform beam
    f       = 11.0/420.0            # special case solution for uniform beam
    L4      = L*L*L*L
    A       = LHS/(3*L4)

    dRHS    = M/(9.0*L)
    coeff   = A*E*r*r/(2.0*rho*wn*wn) - f
    mass    = dRHS/coeff            # mass per unit span

    const   = 2*pi*r*rho

#min gauge check for tube thickness
    t       = mass/(const)
    t       = max(tmin,t)
    mass    = t*const
    mass    = mass*L
    return mass

#======================================================================
# function to calculate spar mass and tube thickness
# for a tapered cantilever beam
#======================================================================

def spar_mass(r1,taper,L,MbyA,Masses,xposn,Thrusts,wing,fn=6.0,tbyc=0.168):

    """
    this function calculates the mass of a cantilever beam using a structural 
    dynamics analytical solution for natural frequency with non-structural 
    mass (distributed and lumped)

    Inputs
    1. r1       = cross-section tube radius at the root, meters
    2. taper    = wing taper ratio = tip chord/root chord. [0.3 to 1.05]
    3. L        = beam length      = wing half-span, meters
    4. MbyA     = mass per unit area of non-structural masses, kg/sq.m
    5. Masses   = list of lumped masses, kgs => motor mounts, winglets, motors and rotors
    6. xposn    = list of dimensional distances for each lumped mass, meters 
    7. T        = list of rotor thrusts, Newtons
    8. wing     = wing class from HYDRA with lift, thrust, torque in all segments
    9. fn       = target natural frequency for half-wing, Hz 
    10.tbyc     = wing thickness to chord ratio

    Outputs
    1. t        = spar thickness, meters 
    2. Mspar    = spar mass, kgs

    Minimum gauge is checked at the end, and larger thickness is chosen
    """

    LHS, RHS, dRHS  = coefficients_simple(r1,taper,L,Masses,xposn,tbyc)
    # print('dRHS is ',dRHS)
    omegan          = fn*2*pi
    beta            = MbyA/tbyc
    coef            = E*pi*0.5/(omegan*omegan)*LHS-rho*pi*RHS 
    b               = beta*RHS + dRHS
    t               = b/coef
    rbar            = (1+taper)*0.5*r1
#cross-check with FEM solution: wn should match fn input value
#    fn_FEM          = taper_beam(MbyA, L, r1, t, taper, Masses, xposn)
#    print(fn_FEM,fn);quit('OK MATCH FEM?')

#calculate thickness of unidirectional carbon fiber layers
    t               = max(t, tmin)

#======================================================================
#ensure that max static stress levels in nominal flight conditions are not reached
#======================================================================

#first hover: find bending moment due to (thrust - weight)
    xnodes          = xposn
    xnodes.insert(0,0.0)
    Masses.insert(0,0.0)
    Thrusts.insert(0,0.0)
    nelem           = len(xnodes)-1        # draw elements between lumped masses
    nsample         = 5
    xpts            = numpy.zeros((nsample*nelem))
    ipt             = 0
    BM              = numpy.zeros_like(xpts)
    Shear           = numpy.zeros_like(xpts)
    for i in range(nelem):
        dx              = (xposn[i+1] - xposn[i])/float(nsample)
        for j in range(nsample):
            xpts[ipt]   = xposn[i] + dx*float(j)

#loop over lumped masses, add thrust and gravity
            for x,M,T in zip(xposn,Masses,Thrusts):
#                print('motor',x,'points',xpts[ipt])
                if x > xpts[ipt]:           # loads only due to outboard thrust
                    Fz          = T - M*g
                    BM[ipt]     = BM[ipt] + Fz*(x - xpts[ipt])
                    Shear[ipt]  = Shear[ipt] + Fz
            ipt         = ipt + 1

#calculate tube radius at sample points
    rprime          = r1*(taper-1.0)/L
    radius          = r1 + rprime*xpts

#calculate sigma = M*y/I = BM/(8 pi r^2 t) 
    sigma           = numpy.divide(BM*nz, 8*numpy.pi*t*numpy.square(radius))
    SF              = numpy.amin(sigma_max/sigma)

# assume equal number of layers for +/- 45 deg carbon fiber, calculate peak shear stress
    tshear          = t 
    CSArea          = 2*numpy.pi*radius*tshear
    tau             = 2.0*Shear/CSArea*nz
    SF2             = numpy.amin(tau_max/tau)

#======================================================================
# if minimum safety factor is less than 1.5, adjust thickness of tube
#======================================================================

    if(SF < 1.5):
        t           = t*1.5/SF 
    if(SF2 < 1.5):
        tshear      = tshear*1.5/SF2 
#    print(SF,SF2)
#======================================================================
# calculate total spar mass
#======================================================================
    
    mspar           = 2*pi*rbar*(t+tshear)*rho       # spar mass, kgs
    Mspar           = mspar*L 

    return Mspar

#======================================================================
# LHS, RHS coefficients
#======================================================================

def coefficients_simple(r1,taper,L,Masses,xposn,tbyc):

    """
    Function to calculate coefficients of KE and PE expressions for 
    tapered beam with lumped masses 

    This function is a dependency used only by spar_mass above.

    Cubic polynomials obtained from curve fits of exact expressions
    """

    r2      = taper*r1
    rprime  = (r2-r1)/L 
    ratio   = r1/r2 

#======================================================================
# alternate formulation of LHS terms
#======================================================================

    cbar    = r1/tbyc*(1+taper)             # mean chord
    AR      = 2*L/cbar                      # wing aspect ratio 
    K       = 2.0*tbyc/AR
    K3inv   = 1.0/(K*K*K)
    taper2  = taper*taper 
    taper3  = taper*taper2

#cubic approx. for LHS terms sum
    ABC2    = 1.626090567083019 - 0.482597591190902*taper + \
              2.119679981957482*taper2 - 0.595753854491954*taper3

    LHS     = ABC2*K3inv

#======================================================================
# alternate formulation for deflection coefficients
#======================================================================

    logr1   = log(r1)
    
#======================================================================
# v2.0 calculations for individual terms of RHS
#======================================================================

    scale   = L*L*K3inv/(K*K)

    RHS_C   = 0.369176557170034*taper3 + 0.151779690454288*taper2 + \
              0.269000169945435*taper  + 0.048225660574691

    RHS_v3  = RHS_C*scale

#======================================================================
# accummulate RHS
#======================================================================

    dRHS    = 0.0
    b0bar, brbar, b1bar, blbar, bmbar = bbars(taper)
    b0      = (b0bar+brbar*logr1)*K3inv
    b1      = b1bar*0.5/L*K3inv 
    bm      = bmbar*L*K*0.5*K3inv 
    bl      = blbar*K3inv
    for Mk,xk in zip(Masses,xposn):

        r       = r1 + rprime*xk
        w       = b0 + b1*xk + bm/r+bl*log(r)
        dRHS    = dRHS + Mk*0.5*w*w

    return LHS, RHS_v3, dRHS

#======================================================================
# bbar coefficients for RHS calculations
#======================================================================

def bbars(taper):

    """
    Function to calculate coefficients of deflection terms for a tapered
    cantilever beam under tip load.

    This function is a dependency used only by coefficients_simple above.
    """

    temp    = ((taper+1)/(taper-1))
    temp    = temp*temp*temp
    b1bar   = temp*(taper-1.0)*(taper-2.0)
    bmbar   = taper*temp/(taper+1)#(taper+1)**2/(taper-1)**3
    blbar   = temp                #((taper+1)/(taper-1))**3
    b0bar   = -temp*0.5*taper
    brbar   = -temp
    return b0bar, brbar, b1bar, blbar, bmbar 