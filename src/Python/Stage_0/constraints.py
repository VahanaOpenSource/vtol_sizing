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
from copy import copy

sys.path.append('../Stage_1')
sys.path.append('../Stage_1/afdd')
sys.path.append('../Stage_1/engines')
sys.path.append('../Stage_3')

from hydraInterface    import hydraInterface
from conversions import *

#====================================================================
# wing span calculator for first wing group
#====================================================================

def wing_span_calc(x, design):

    nwing           = design.wing.groups[0].nwings                  # wings in group
    rho             = design.mission.segment[0].rho                 # density, kg/cu.m
    Vinf            = design.mission.segment[1].cruisespeed*kts2mps # V_inf  , m/s
    q               = 0.5*rho*Vinf*Vinf
    AR              = x[4]
    CL              = x[5]

#====================================================================
# if there is more than one group, lift fraction is specified for the 
# second group as a design variable; use it to calculate lift frac 
# for the first group; assume there are only two groups
#====================================================================

    if(design.wing.ngroups > 1):
        lf          = 1.0 - x[8]
    else:
        lf          = 1.0

    Lift            = x[3]*grav/nwing*lf        # lift, Newtons
    S               = Lift/(q*CL)               # wing area, sq.m
    wing_span       = sqrt(S*AR)                # wing span, sq.m

    # print('===> group 0;',wing_span)
    # print(Lift, x[3],x[8],lf)

    return wing_span    

#====================================================================
# wing span calculator
#====================================================================

def wing2_span_calc(x, design):

    nwing           = design.wing.groups[1].nwings
    rho             = design.mission.segment[0].rho                 # density, kg/cu.m
    Vinf            = design.mission.segment[1].cruisespeed*kts2mps # V_inf  , m/s
    q               = 0.5*rho*Vinf*Vinf

    lf              = x[8]                    # lift fraction specified for 2nd wing group
    AR              = x[6]                    # aspect ratio,     DV
    CL              = x[7]                    # lift coefficient, DV
    Lift            = x[3]*grav/nwing*lf      # lift, Newtons

    S               = Lift/(q*CL)             # wing area, sq.m                                    
    wing_span       = sqrt(S*AR)              # span, sq.m

    return wing_span

#====================================================================
# Wing span calculation 
# wing span <= bmax is an ok design (output <= 0)
#====================================================================

def wing_span_constraint(x, *args):

    design          = args[0]
    # print(dir(design));x2 = input('yeahhh?')
    bmax            = design.constraints.max_rotor_radius            # limit, meters
    wing_span       = wing_span_calc(x,design)
    # print(wing_span-bmax)
    # print('wing span constraint here: ',wing_span, bmax)
    return (bmax - wing_span)

#====================================================================
# Wing span calculation 
# wing span <= bmax is an ok design (output <= 0)
#====================================================================

def wing2_span_constraint(x, *args):

    design          = args[0]
    # print(dir(design));x2 = input('yeahhh?')
    bmax            = design.constraints.max_rotor_radius            # limit, meters
    wing_span       = wing2_span_calc(x,design)

    return (bmax - wing_span)

#====================================================================
# clearance parameter for rotors, wings and fuselage to make sure 
# everything fits 
# clearance >= 0 is ok design, so (-clearance) <= 0
#====================================================================

def clearance2_constraint(x, *args):

    design          = args[0]
    # print(design);x1=input('ehh?')
    summary         = run_sizing(x, design)
    R               = summary['Radius']     
    bmax            = design.constraints.max_rotor_radius
    b               = wing2_span_calc(x, design)
    geom            = design.emp_data.Geometry
    bf              = geom.fuselage_width
    delta           = geom.clearance 
    nR              = design.wing.groups[1].nrotors

    clearance       = 0.5*(b - bf) - R*( (1+delta)*nR - 1)
    # print(0.5*(b-bf), R*( (1+delta)*nR - 1),-clearance)
    return clearance

#====================================================================
# clearance parameter for rotors, wings and fuselage to make sure 
# everything fits 
# clearance >= 0 is ok design, so (-clearance) <= 0
#====================================================================

def clearance_constraint(x, *args):

    design          = args[0]
    bmax            = design.constraints.max_rotor_radius
    summary         = run_sizing(x, design)
    R               = summary['Radius']     
    b               = wing_span_calc(x, design)
    geom            = design.emp_data.Geometry
    bf              = geom.fuselage_width
    size            = 1 + geom.clearance
    nR              = design.wing.groups[0].nrotors

    clearance       = 0.5*(b - bf) - R*( size*nR - 1)
    # print('=====> wing clearance constraint here: ',b, bf, clearance)
    # print('radius = ',R)
    # print(0.5*(b-bf), R*( (1+delta)*nR - 1),-clearance)
    return clearance

#====================================================================
# adapter functions to take design variable vector "x" as an input 
# and call the constraint function for blade loading CT/sigma
#====================================================================

def CTsig_constraint(x, *args):

    design     = args[0]
    DL          = x[0]*lb2N/(f2m*f2m)               # in N/sq.m
    rho         = design.mission.segment[0].rho     # density, kg/cu.m
    Vtip        = x[1]                              # tip speed, m/s
    sigma       = x[2]                              # solidity
    maxval      = design.constraints.max_ct_sigma   # maximum CT/sigma 
    
    CTsigma     = DL/(rho*Vtip*Vtip*sigma)

    # print('in CT/sig function',CTsigma);x1=input('OK?')
    return (maxval - CTsigma)

#====================================================================
# Payload check: should get M_payload >= M_target payload 
# so (M_tar - M_pay <= 0)
#====================================================================

def payload_constraint(x, *args):

    design     = args[0]
    summary    = run_sizing(x,design)
    mpay       = summary['Payload']
    mtar       = design.mission.payload_tar
    return (mpay - mtar)

#====================================================================
# constraint assembler function: dict of constraint functions
# x = design variables, contains following in an array
#     for identical tandem wings (Alpha-like) 
#====================================================================
# x[0] = Disk loading, lb/sq.ft
# x[1] = Tip speed, m/s
# x[2] = rotor solidity
# x[3] = take off mass, kgs
# x[4] = wing aspect ratio
# x[5] = wing CL
#------
# for one wing group, that's it; for multiple wing groups, keep reading
#-----
# x[6] = wing aspect ratio,  second group
# x[7] = wing CL,            second group
# x[8] = wing lift fraction, second group
#
# if there is another group in this design, then the lift fraction 
# for the second group is an input; lift for the first group in 
# cruise condition is automatically determined as 1-x[8], assuming
# that all lift is carried by fixed wings in cruise
#====================================================================

def constraint_assembler(design):

#====================================================================
# design variable min/max limits
#====================================================================

    bnds    = []

#====================================================================
# the following bounds are the same regardless of # wings
#====================================================================

    mTO     = design.massTakeoff

    b0      = (2.0        ,   18.0 )                # disk loading bounds
    b1      = (50         ,  172.0 )                # tip speed bounds
    b2      = ( 0.01      ,    0.17)                # solidity  bounds
    b3      = ( 0.5*mTO   ,1.5*mTO )                # take-off mass: 10% starting to 4 tons

    bnds.append(b0)
    bnds.append(b1)
    bnds.append(b2)
    bnds.append(b3)

#====================================================================
# wing design variables: varies with # wings
#====================================================================

    b4      = ( 2  , 14.0 )                # wing aspect ratio
    b5      = ( 0.2 , 0.7 )               # wing CL   bounds

    bnds.append(b4)
    bnds.append(b5)

#====================================================================
# if there is another wing, add its AR, CL and lift fraction
#====================================================================

    b8      =   (0.1, 0.5)                 # lift fraction can be between 0 and 1

#====================================================================
# constraints
#====================================================================

    con0    = {'type': 'ineq', 'fun': payload_constraint,   'args': (design,) }
    con1    = {'type': 'ineq', 'fun': CTsig_constraint,     'args': (design,) }
    con2    = {'type': 'ineq', 'fun': wing_span_constraint, 'args': (design,) }    
    con3    = {'type': 'ineq', 'fun': clearance_constraint, 'args': (design,) }

    cons    = ([con0,con1,con2,con3]) 

#====================================================================
# more than one wing group
#====================================================================

    if design.wing.ngroups > 1:
        con4    = {'type': 'ineq', 'fun': wing2_span_constraint,'args': (design,) }    
        con5    = {'type': 'ineq', 'fun': clearance2_constraint,'args': (design,) }
        cons.append(con4)
        cons.append(con5) 

        bnds.append(b4)         # aspect ratio bounds 
        bnds.append(b5)         # lift coefficient bounds 
        bnds.append(b8)         # lift fraction    bounds

    return bnds, cons   

#====================================================================
# cost function: literally cost calculation
#====================================================================

def cost_function(x, *args):

    # print('in cost function');quit('OK?')
    design     = args[0]
    summary    = run_sizing(x, design)

    return summary['Cost'] #+ summary['Weight'])

#====================================================================
# function to take design variables and return summary dictionary
#====================================================================

def run_sizing(x, design):

    # print(design.com)
    com        = {'disk_loading': x[0], 'vtip': x[1],  'solidity': x[2], 
                  'wing_group0_aspectratio': x[4],'wing_group0_cl': x[5]}

#========================================
# assign wing group #2 DVs if applicable
#========================================

    if(design.wing.ngroups > 1):
        com['wing_group1_aspectratio']  = x[6]
        com['wing_group1_cl']           = x[7]
        com['wing_group1_liftfraction'] = x[8]

#========================================
# run sizing
#========================================

    mission    = design.mission

    design.com.update(com)
    design.modify_case(design.com)

#========================================
# run in fixed take-off weight mode
#========================================

    mission.sizing_mode     = 2
    mission.gtow_target     = x[3]
    mTO                     = mission.mass_takeoff_guess()
    design.massTakeoff      = mTO

    # print(x[3],design.massTakeoff,mTO); x1 = input('OK?')
#========================================
# calculate payload
#========================================

    design.kick_off()
    summary    = design.get_essentials()
    # print(x)
    # print('hello',a)
    # x1 = input('yes?')
    return summary