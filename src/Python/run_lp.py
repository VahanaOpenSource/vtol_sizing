import numpy, copy 
from scipy.optimize import linprog

def run_lp(als, rotDir, dTdTau, dMzdTau, ndof, ncontrols, \
           x_rotor, y_rotor, MzSurface, linp_dict ):

    """
    This function uses a linear programming formulation to calculate the maximum 
    thrust ratio and torque ratio for an eVTOL motor wrt nominal hover, to ensure 
    that 4dof trim can still be maintained when a motor fails in hover.

    Inputs:
    1. als       = array of rotor shaft tilt angles wrt body (-Z) axis
    2. rotDir    = array of integers with (+1) for CCW, (-1) for CW showing rotation direction
    3. dTdTau    = array of rotor thrust-to-torque ratios, 1/m
    4. ndof      = number of degrees of freedom in trim problem => must be 4
    5. ncontrols = number of independent actuators in the vehicle; each rotor is one and 
                   all control surfaces are aggregated into one
    6. x_rotor   = array of rotor longitudinal offsets from CG, meters
    7. y_rotor   = array of rotor lateral      offsets from CG, meters
    8. MzSurface = maximum yawing moment that can be generated from all surfaces, N-m
    9. linp_dict = dictionary of coefficient matrices/vectors used to setup linear programming problem
                   This dictionary must contain the following matrices and vectors 
                (a)    Aeq  = matrix of size "(ncases * ndof) x (ncontrols x ncases + 1)"  
                (b)    beq  = vector of size "(ncases * ndof) x 1"
                (c)      A  = matrix of size "(ncontrols x ncases) x (ncontrols x ncases + 1)"
                (d)      b  = vector of size "(ncontrols x ncases) x 1"
                (e) bounds  = list of size-2 tuples, number of entries = ncontrols x ncases + 1
                (f)      f  = array of coefficient of weights to minimize "x", size "ncontrols x ncases + 1"
                (g)    lb2  = array of lower bounds for maximum yawing moment problem, size "ncontrols x 1"
                (h)    ub2  = array of upper bounds for maximum yawing moment problem, size "ncontrols x 1"

    Outputs:
    1. res.x    = solution vector of linear programming problem for minimizing max. motor torque
    2. res.fun  = maximum hover torque needed for any one motor inoperative, Nm
    3. meantau  = mean hover torque, N-m 
    4. res2.fun = maximum yawing moment that can be produced in hover, Nm, no failures

    return res.x, res.fun, meantau, res2.fun
    """
#=================================================================================
# unpack matrices from linear programming dictionary
#=================================================================================

    Aeq         = linp_dict['Aeq']
    beq         = linp_dict['beq']
    A           = linp_dict['A']
    b           = linp_dict['b']
    bounds      = linp_dict['bounds']
    f           = linp_dict['f']
    lb2         = linp_dict['lb2']
    ub2         = linp_dict['ub2']

    nrotor      = len(rotDir)
    dFMdTau     = numpy.zeros((ndof,ncontrols))

#=================================================================================
# loop over rotors, calculate elements of actuator Jacobian 
# rows are: thrust (+ve down), roll moment, pitching moment and yawing moment
# body axes convention followed
#
# first "nrotor" columns correspond to individual rotors
# last column corresponds to effect of all control surfaces
#=================================================================================

    for i in range(nrotor):

#=================================================================================
# sensitivity of vertical thrust to individual rotor torques
#=================================================================================

        dFMdTau[0,i] = -dTdTau[i]*numpy.cos(als[i])                   # vertical thrust

#=================================================================================
# rolling moment sensitivity to rotor torque
#=================================================================================

        dFMdTau[1,i] =  dTdTau[i]*numpy.cos(als[i])*y_rotor[i]        # rolling  moment

#=================================================================================
# pitching moment sensitivity to rotor torque
#=================================================================================

        dFMdTau[2,i] =  dTdTau[i]*numpy.cos(als[i])*x_rotor[i]        # pitching moment

#=================================================================================
# yawing moment sensitivity to rotor torque
#=================================================================================

        dFMdTau[3,i] =  dTdTau[i]*numpy.sin(als[i])*y_rotor[i]        # yawing   moment
    dFMdTau[3,8]     = 1.0

#=================================================================================
# roll-yaw coupling due to shaft tilt wrt body z-axis
#=================================================================================

    for i in range(nrotor):
        dFMdTau[1,i] -= dMzdTau * rotDir[i] * numpy.sin(als[i])
        dFMdTau[3,i] += dMzdTau * rotDir[i] * numpy.cos(als[i])

#=================================================================================
# initialize equality constraints coefficient matrix
#=================================================================================

    for i in range(ndof):
        for j in range(ncontrols):
            Aeq[i,j]  = dFMdTau[i,j]

#=================================================================================
# Note:ignoring following nonlinear effect for the time being
# How does the yawing moment scale with the torques and control deflections as the thrust changes?
# Non-linear. There is a thrust * deflection term
#=================================================================================
            
#=================================================================================
# calculate nominal hover torques using pseudo inverse
#=================================================================================

    hoverTau    = numpy.linalg.pinv(dFMdTau[:,0:nrotor])
    hoverTau    = numpy.dot(hoverTau, beq[0:ndof])

#=================================================================================
# Run linprog of all the rotor out options
#=================================================================================

#=================================================================================
# setup equality constraint coefficient matrix to enforce trim
# done in the loop because diagonal blocks of the matrix depend on tilt angles
#=================================================================================

    for k in range(nrotor):
        tmp      = copy.copy(dFMdTau)
        tmp[:,k] = 0.0

        for i in range(ndof):
            for j in range(ncontrols):
                ioffset     = i + (k+1)*(ndof)
                joffset     = j + (k+1)*(ncontrols)
                Aeq[ioffset,joffset] = tmp[i,j]

#=================================================================================
# Solve for maximum motor size required
#=================================================================================

    res     = linprog(f, A_ub=A,b_ub=b,A_eq=Aeq, b_eq=beq, bounds=bounds,method='interior-point')
    
    if not(res.success):
        print('linear programming problem DID NOT CONVERGE!!')
        print('warning: COULD NOT FIND MAX TORQUE LIMIT!')
    meantau                 = numpy.mean(hoverTau)

#=================================================================================
# calculate max yawing moment in hover
#=================================================================================

#=================================================================================
# setup bounds for linear programming problem
#=================================================================================

    ub2     = ub2*0.0 + res.fun
    ub2[-1] = MzSurface            

    bounds2 = [] 
    for l,u in zip(lb2,ub2):
        bounds2.append((l,u))

#=================================================================================
# coefficients of dFMdTau -> used to calculate yawing moment
#=================================================================================

    f2      = dFMdTau[-1,:]
    Aeq2    = dFMdTau[0:3,:]
    beq2    = beq[0:3]
    res2    = linprog(f2,A_ub=None,b_ub=None,A_eq=Aeq2,b_eq=beq2,bounds=bounds2,method='interior-point')
#    if not(res2.success):
#        print(res2)
#        print('linear programming problem for max yawing moment DID NOT CONVERGE!!')
#        print(numpy.mean(als[0:4])*180.0/numpy.pi)
    return res.x, res.fun, meantau, res2.fun