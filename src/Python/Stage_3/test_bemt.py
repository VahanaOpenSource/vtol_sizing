import bemt, numpy

def run_bemt_once(print_flag):
#======================================================================
# populate inputs to bemt
#======================================================================

    interface               = bemt.bemt_interface
    inputs                  = interface.inputs
    inputs.nb               = 2                      # number of blades
    inputs.radius           = 0.3e0                  # radius, m
    inputs.solidity         = 0.12e0                 # blade solidity
    inputs.omegah           = 3000.0*numpy.pi/30.0   # hover rotor speed, rad/s
    inputs.use_tables       = 1                      # 1 for True, 0 for False
    inputs.cl_alpha         = 5.73e0                 # 1 for True, 0 for False
    inputs.cdo              = 0.028e0                # 1 for True, 0 for False
    inputs.sweep_rpm        = 0                      # 1 = test all cruise RPM
        
#======================================================================
# set flight condition details
#======================================================================

    inputs.n_flight         = 2
    flt                     = interface.flt

#======================================================================
# 5 min hover segment, MSL
#======================================================================

    flt[0].time             =  300.0e0          # segment time, seconds
    flt[0].vc               =    0.0e0          # climb velocity, m/s
    flt[0].thrust           =   58.0e0          # thrust, N
    flt[0].rho              =    1.225e0        # density, kg/cu.m
    flt[0].altitude         =    0.00e0         # altitude in meters

#======================================================================
# 30 min cruise segment, MSL
#======================================================================

    flt[1].time             = 1800.0e0          # 30 min segment
    flt[1].vc               =   18.0e0          # climb velocity, m/s
    flt[1].thrust           =   11.6e0          # thrust, N
    flt[1].rho              =    1.225e0        # density, kg/cu.m
    flt[1].altitude         =    0.00e0         # altitude in meters

#======================================================================
# call bemt analysis: calculate performance
#======================================================================

    bemt.setup_bemt()

#======================================================================
# screen output for results
#======================================================================

    if print_flag:
        for i in xrange(inputs.n_flight):
            print ('================================')
            print ('Flight condition # %3d' % (i+1))
            print ('================================')
            print ('Root collective  is %12.3f degrees' % flt[i].coll)
            print ('Rotor efficiency is %12.3f' % (flt[i].efficiency))
            print ('Power required   is %12.3f watts' % flt[i].power)
            print (' ')

        print ('optimal design is ')
        print (bemt.bemt_interface.best_design)

    print ('\nInitialized bemt for sizing calculations')

    return None 

#======================================================================
# running bemt once initializes the airfoil table and runs a sample case
#======================================================================
