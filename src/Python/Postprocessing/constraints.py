#====================================================================
# Calculating values used to build constraint functions
# Disk loading = DL, N/sq.m
# Tip speed    = Vtip, m/s
# sigma        = solidity 
# rho          = air density, kg/cu.m
# Lift         = Lift from each wing in a group, N 
# Vinf         = cruise speed, m/s
# CL           = wing lift coefficient 
# AR           = aspect ratio (span/chord) of a wing 
# b            = wing span, m 
# b_f          = fuselage width, m 
# R            = rotor radius, m 
# nR           = number of rotors per wing (tip to tip)
# delta        = clearance between rotor and fuselage/another rotor
#                as fraction of rotor radius (0.05 first cut)
#====================================================================
import sys 
import numpy 

sys.path.append('../Stage_1')
sys.path.append('../Stage_1/afdd')
sys.path.append('../Stage_1/engines')
sys.path.append('../Stage_3')

from hydraInterface    import hydraInterface
from conversions import *

#====================================================================
# Wing span calculation 
# wing span <= bmax is an ok design (output <= 0)
#====================================================================

def footprint_constraint(x, *args):

    design          = args[0]
    summary         = run_sizing(x,design)
    bmax            = design.constraints.max_rotor_radius            # limit, meters
    footprint       = summary[7]

    return (bmax - footprint)

#====================================================================
# Payload check: should get M_payload >= M_target payload 
# so (M_payload - M_target >= 0)
#====================================================================

def payload_constraint(x, *args):

    design     = args[0]
    summary    = run_sizing(x,design)
    mpay       = summary[5]
    mtar       = design.mission.payload_tar
    return (mpay - mtar)

#====================================================================
# constraint assembler function: dict of constraint functions
# x = design variables
#====================================================================
# x[0] = take off mass, kgs
# x[1] = ..
# x[2] = ..
# x[3] = ..
# x[4] = ..
# x[5] = ..
#====================================================================

def constraint_assembler(design):

    bnds    = []
    mTO     = design.massTakeoff

    bnds.append((0.5*mTO, 1.5*mTO))               # first DV is take-off mass

    for k in design.var_keys:
        if(k.startswith('rotor')):
            k1      = k.split('_')[-1]            # type of variable
            if(k1 == 'Vtip'):
                bnds.append((100,170))
            elif(k1 == 'ctsigma'):
                bnds.append((0.05,0.139))
            else:
                print('cannot understand rotor design variable type',k1)
                quit('CRITICAL ERROR: constraints.py, constraint assembler')
        elif(k.startswith('wing')):
            k1      = k.split('_')[-1]            # type of variable
            if(k1 == 'aspectratio'):
                bnds.append((3.0,10.0))
            elif(k1 == 'cl'):
                bnds.append((0.2,0.8))
            elif(k1 == 'liftfraction'):
#                print('preserving the same lift fraction: not adding it as design variable')
                if(design.wing.ngroups > 1):
                    lf_max  = 0.95
                    lf_min  = 0.05
                    bnds.append((lf_min, lf_max))
            else:
                print('cannot understand rotor design variable type',k1)
                quit('CRITICAL ERROR: constraints.py, constraint assembler')

#====================================================================
# constraints: payload must be what I say it is
#====================================================================

    con0    = {'type': 'ineq', 'fun': payload_constraint,   'args': (design,) }
    con1    = {'type': 'ineq', 'fun': footprint_constraint, 'args': (design,) }

#====================================================================
# solidity constraint for all rotor groups: currently deactivated
#====================================================================

#    con1    = {'type': 'ineq', 'fun': solidity_constraint,  'args': (design,) }

#====================================================================
# footprint constraint: deactivated
#====================================================================

#    con2    = {'type': 'ineq', 'fun': wing_span_constraint, 'args': (design,) }    

#====================================================================
# more than one wing group
#====================================================================

#    if design.wing.ngroups > 1:
#        con4    = {'type': 'ineq', 'fun': wing2_span_constraint,'args': (design,) }    
#        cons.append(con4)
#        cons.append(con5) 

#        bnds.append(b4)         # aspect ratio bounds 
#        bnds.append(b5)         # lift coefficient bounds 
#        bnds.append(b8)         # lift fraction    bounds

    cons    = ([con0,con1]) 

    return bnds, cons   

#====================================================================
# cost function: literally cost calculation
#====================================================================

def cost_function(x, *args):

    # print('in cost function');quit('OK?')

    design     = args[0]
    summary    = run_sizing(x, design)
    cost       = summary[6]
    # print_dvars(x,design)

#====================================================================
# error trap
#====================================================================

    if(numpy.isnan(cost)):
        print_dvars(x,design)
        quit('BANG: problem in cost function')
#====================================================================
#augment cost function with "death penalty" by doubling the cost and 
# adding $2k/hr
#====================================================================

    mpay       = summary[5]
    mtar       = design.mission.payload_tar
    bmax       = design.constraints.max_rotor_radius            # limit, meters
    footprint  = summary[7]
    if(mpay < mtar or footprint > bmax):
        cost   = 20000
#        print('TO THE DEATH!',cost,mpay,mtar,x)
#        quit()
    return cost

#====================================================================
# function to take design variables, run sizing + return summary
#====================================================================

def run_sizing(x, design):

    com        = {}
    wing       = design.wing 
    rotor      = design.rotor

#====================================================================
# design variables for rotor
#====================================================================

    for ikey,key in enumerate(design.var_keys):
        com[key]    = x[ikey+1]

#========================================
# run sizing
#========================================

    design.com.update(com)
    design.modify_case(design.com)

#========================================
# run in fixed take-off weight mode
#========================================

    mission                 = design.mission                # pointer to class
    mission.sizing_mode     = 2                             # fixed take-off weight sizing mode
    mission.gtow_target     = x[0]                          # set target take-off mass
    mTO                     = mission.mass_takeoff_guess()  # 
    design.massTakeoff      = mTO                           # set take-off mass initial condition

#==================================================
# calculate payload, cost & big-picture parameters
#==================================================

    design.kick_off()
    summary    = design.get_essentials()
    return summary

#======================================================================
# function to print design variables
#======================================================================

def print_dvars(x, design):
    bnds    = design.bounds
    nDV     = design.nDV
    print('\n=================================================')
    print('DESIGN VARIABLES, final solution and bounds:       ')
    print('=================================================\n')
    print('{:30s} {:12.3f} {:12.3f} to {:12.3f}'.format('x[0]: take off mass (kg)',x[0],bnds[0][0],bnds[0][1]))
    for i in range(1,nDV):
        print('{:30s} {:12.3f} {:12.3f} to {:12.3f}'.format('x['+str(i)+']: '+design.var_keys[i-1],x[i],bnds[i][0],bnds[i][1]))


    return None 