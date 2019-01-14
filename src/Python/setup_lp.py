import numpy
def setup_lp(nrotor, MzSurface, W, ndof=4):
    """
    this function setups up a linear programming problem to evaluate 
    maximum overload factor to use for any rotor in case of OMI type 
    failure

    the second linear programming problem it sets up is the calculation 
    of maximum yawing moment generated with all motors working normally

    inputs are 
    1. nrotor   : number of rotors in the system
    2. MzSurface: yawing moment due to control surface deflection
    3. W        : vehicle weight in Newtons
    2. ndof     : number of independent DOF to consider = 4 (thrust, roll)
    """

#=================================================================================
# setup array sizes
#
# The unknown vector {x} has the following form 
#
# Torque rotor 1, case 1
# Torque rotor 2, case 1 
# Torque rotor 3, case 1
# ...
# ...
# Control surface yawing moment, case 1 
# Torque rotor 1, case 2 
# Torque rotor 2, case 2 
# .. 
# ..
# Control surface yawing moment, case 2 
# ..
# pattern repeats for other cases
# ..
#
#=================================================================================

    large       = 1e10                              # proxy for infinity, will be replaced with None
    ncontrols   = nrotor+1                          # number of available actuations
    ncases      = nrotor+1                          # number of cases; take nominal + failure of each rotor
#    nDV         = (nrotor+1);   nDV     = nDV*nDV+1 
#    ncols       = nDV-1                             # number of inequality constraints
                                                    # last column is max torque, added as slack variable
    ncols       = ncases*ncontrols                  # number of design variables = (# actuators x # cases + 1)
    nDV         = ncols+1

#=================================================================================
# setup lower/upper bounds for thrust overload problem
#=================================================================================

    lb          = numpy.zeros(nDV)                  # initialize storage for lower and upper bounds
    ub          = numpy.zeros_like(lb) + large      # infinite upper bound on all DV, except control surfaces
    bounds      = []                                # list of upper and lower variable bounds

    for i in range(ncases):                         # for each failure case, set appropriate limits on CS yawing moment 
        lb[nrotor+i*ncontrols] = -MzSurface
        ub[nrotor+i*ncontrols] =  MzSurface

    f       = numpy.zeros_like(lb)                  # coefficient array for DV vector that is minimized 
    f[-1]   = 1.0                                   # only max torque is minimized
    f       = numpy.transpose(f)                    # for consistency

#=================================================================================
# setup constraints for thrust overload problem
# each motor torque <= max torque => x_i <= x_end, i < end
# represented as coefficient matrices "A", "b" for [A]{x} <= b 
# b matrix is zero
#=================================================================================

    b       = numpy.zeros((ncols,1))
    A       = numpy.zeros((ncols,nDV))              # inequality constraint coefficient matrix 
    for i in range(ncols):                          
        A[i,i]  =  1.0      
        A[i,-1] = -1.0

    for i in range(ncases):                         # remove constraint for yawing moment of surfaces
        rowid           = nrotor + ncontrols*i
        A[rowid,rowid]  = 0.0

#=================================================================================
# setup equality constraints matrix
#
# [Aeq] {x} = {beq}
# 
# The equality constraints to be met are the trim conditions for 
# 4 degrees of freedom, in the following order
# Row 1: thrust; 
# Row 2: rolling moment
# Row 3: pitching moment 
# Row 4: yawing moment
#
# The individual contributions of each rotor to the vehicle loads in hover 
# are calculated using control derivatives as follows 
#
# d(Thrust)       C_T
# --------- = ------------  
# d(Torque)   C_P x radius
#
# 1. For a given rotor radius, tip speed, vehicle weight, density and  
#    rotor thrust share in hover, we can calculate C_T
# 2. From figure of merit, we get C_P from C_T
# 
# Representing thrust as "T" and torque as "Tau", this ratio above is dTdTau
#
# For the effect of this thrust on rolling/pitching moments on the body, 
# we multiply thrust by lateral offset and longitudinal offset from CG, respectively
# For the effect of individual thrusts on yawing moments, we need to take into 
# account rotor shaft mounting angles relative to the wings, because the 
# longitudinal component of rotor thrust (along body X direction) coupled with 
# the lateral offset of the rotor hub produces a yawing moment.
#
# The key assumption is that the rotor figure of merit is not changed
# significantly with thrust condition; then thrust and torque are both 
# proportional to RPM squared, and hence scale linearly with each other
#
# Represent the body loads as linear combinations of contributions from each rotor
# and all the control surfaces
#
# Fz        = Fz1 + Fz2 + Fz3 + ... + FzN 
# Mx        = Mx1 + Mx2 + Mx3 + ... + MxN 
# My        = My1 + My2 + My3 + ... + MyN 
# Mz        = Mz1 + Mz2 + Mz3 + ... + MzN + Mz_CS
#
# The control surfaces only affect yawing moment
# Because rotor loads are proportional to rotor thrust, and (under the assumption of 
# constant figure of merit) rotor torque, we can represent body loads as a matrix 
# vector product 
#
# FM    = { Fz   Mx   My  Mz}^T
# Tau   = { Tau1 Tau2 ....  TauN Mz_CS}^T
#
# {FM}  = [dFMdTau] {Tau}
#
# The "Tau" vector represents the torques of the individual rotors
#
# [dFMdTau] is a matrix of 4 rows and (Nrotor+1) columns
#
# In trim, {FM} must be equal to {-W, 0, 0, 0} ==> call this RHS as "beq"
#
# The requirement to trim the vehicle can be expressed as equality constraints
# [dFMdTau]{Tau} = {b_eq}
#
# If we're analyzing multiple cases (nominal, various rotor out), then we have 
# several sets of equality constraints. 
#
# Each case is independent, so we can assemble all the equality constraints into a 
# matrix, with diagonal blocks of [dFMdTau], and vertcal blocks of {W,0,0,0}^T
#
#=================================================================================

    neq             = ndof*ncontrols
    Aeq             = numpy.zeros((neq, nDV))
    beq             = numpy.zeros(neq)

    for i in range(ncontrols):
        beq[ndof*i] = -W

#=================================================================================
# for each failure case, switch off a rotor ==> set upper bound to zero
# for rotor torque design variable
# ..
# first     case is nominal        ==> all rotors on
# second    case is failed rotor 1 ==> row index (  ncontrols+1) ==> Python row index (ncontrols)
# third     case is failed rotor 2 ==> row index (2*ncontrols+2) ==> Python row index (2*ncontrols+1)
# fourth    case is failed rotor 3 ==> row index (3*ncontrols+3) ==> Python row index (3*ncontrols+2)
# ...
#=================================================================================

    for i in range(nrotor):
        ub[(ncontrols+1)*i+ncontrols] = 0

# assemble bounds as a list
    for i in range(nDV):
        lowerval    = lb[i]
        upperval    = ub[i]
        if(upperval == large):              # use "None" ==> Python's notation for unbounded
            upperval    = None
        bounds.append((lowerval,upperval))

#=================================================================================
# Once we solve for maximum overload factor over all failure cases, 
# we've sized the motors. Using that motor size, we then proceed to 
# estimate the maximum yawing moment that can be produced while maintaining 
# trim along the thrust direction, as well as rolling/pitching moments
#
# setup bounds for max hover authority formulation
#=================================================================================

    lb2     = numpy.zeros(ncontrols)            # do not produce negative torque
    lb2[-1] = -MzSurface                        # maximum yawing moment from control surfaces
    ub2     = numpy.ones(ncontrols)             # upper bound: will be set based on solution to first problem

    optim_details   = {'A': A, 'b': b, 'Aeq': Aeq, 'beq': beq, 'f': f, \
                       'bounds': bounds, 'lb2': lb2, 'ub2': ub2}
    return optim_details

#=================================================================================
# end of operations
#=================================================================================

def get_rotor_dirn(kk):
    """
    this function gets an array of integers to specify rotor rotation directions
    for a 4-4 tilt-wing layout
    there are 7 layouts available, best one is kk=3

    Input: 
    kk:  1 to 7, integer 

    Output:
    rotDir: array of 8 integers with entries as "+1" (CCW) or "-1" (CW)

    """
    if kk == 1:
        rotDir = [ 1,  1, -1, -1, -1, -1,  1,  1] # Option 1
    elif kk == 2:
        rotDir = [ 1, -1,  1, -1, -1,  1, -1,  1] # Option 2
    elif kk == 3:
        rotDir = [-1,  1, -1,  1,  1, -1,  1, -1] # Option 3
    elif kk == 4:
        rotDir = [-1,  1, -1,  1, -1,  1, -1,  1] # Option 2
    elif kk == 5:
        rotDir = [-1,  1,  1, -1, -1,  1,  1, -1] # Option 5
    elif kk == 6:
        rotDir = [-1,  1,  1, -1,  1, -1, -1,  1] # Option 6
    elif kk == 7:
        rotDir = [-1, -1,  1,  1,  1,  1, -1, -1] # Option 4
    else:
        print('unknown value of kk: defaulting to best known layout')
        rotDir = [-1,  1, -1,  1,  1, -1,  1, -1] # Option 3
    return rotDir 

