#======================================================================
# This function calls the BEMT module in HYDRA and calculates power 
# required for a given RPM
#======================================================================

import bemt, numpy 

class _bemt_model():

  def bemt_model(self):

#======================================================================
# create hydra pointers
#======================================================================

    aircraft                = self.all_dict['aircraft']
    rotor                   = self.rotor
    segment                 = self.mission.segment 
    std                     = self.constants
    kts2mps                 = std.kts2mps

#======================================================================
# populate inputs to bemt
#======================================================================    

    interface               = bemt.bemt_interface
    inputs                  = interface.inputs
    inputs.mode_operation   = self.bemt_mode

    inputs.nb               = rotor.nblade                  # number of blades 
    inputs.radius           = rotor.radius                  # radius, m
    inputs.solidity         = rotor.solidity                # blade solidity
    inputs.omegah           = rotor.tipspeed/rotor.radius   # hover rotor speed, rad/s
    inputs.sweep_rpm        = 0                             # only ONE rpm
    inputs.rpm_ratio        = aircraft['cruise_rpm_ratio']
    # print 'MYELLO ',inputs.rpm_ratio 
    inputs.use_tables       = 1                      # 1 for True, 0 for False
    inputs.cdo              = 0.028e0                # 1 for True, 0 for False
        
#======================================================================
# set flight condition details
#======================================================================

    N                       = self.mission.nseg
    inputs.n_flight         = N
    flt                     = interface.flt

#======================================================================
# 5 min hover segment, MSL
#======================================================================

    for i in range(N):
        flt[i].time             =  segment[i].time       # segment time, seconds
        flt[i].vc               =  segment[i].cruisespeed*kts2mps                  # climb velocity, m/s
        flt[i].thrust           =  segment[i].thrust          # thrust, N
        flt[i].rho              =  segment[i].rho             # density, kg/cu.m
        flt[i].altitude         =  segment[i].startalt      # altitude in meters

#        print flt[i]

#        print ('---------------------------')
#        print ('Flight condition # %3d' % (i+1))
#        print ('---------------------------')
#        print ('Thrust required  is %12.3f lbs' % flt[i].thrust)
#        print ('Rotor speed      is  %9.1f   percent of hover speed' % (flt[i].omega/inputs.omegah*100))
#        print ('time in segment is %12.3f min' % flt[i].time)
#======================================================================
# call bemt analysis: calculate performance
#======================================================================

    print ('running bemt')
    bemt.setup_bemt()
    print ('ran bemt')
#======================================================================
# store segment aero efficiency for HYDRA resizing
#======================================================================

#======================================================================
# screen output for results
#======================================================================

    self.valid       = (bemt.bemt_interface.best_design.valid == 1)

    if self.valid and self.bemt_mode != 2:
        print ('\n')
        for i in range(N):
            # print 'changing efficiency to ',flt[i].efficiency
            eff                     = flt[i].efficiency
            if eff > 0.1:
                segment[i].aero_eta = flt[i].efficiency 
            else:
                self.valid               = False 
                print ('less than 10% rotor efficiency: FUHGEDDABOUDIT')
                
            if self.bemt_mode == 0:
                print ('---------------------------')
                print ('Flight condition # %3d' % (i+1))
                print ('---------------------------')
                print ('Thrust required  is %12.3f lbs' % flt[i].thrust)
                print ('Root collective  is %12.3f degrees' % flt[i].coll)
                print ('Rotor efficiency is %12.3f' % (flt[i].efficiency))
                print ('Rotor speed      is  %9.1f   percent of hover speed' % (flt[i].omega/inputs.omegah*100))
                print ('Power required   is %12.3f hp' % (rotor.num*flt[i].power/746))
                print ('time in segment is %12.3f min' % flt[i].time)
                print (' ')

            print ('optimal design is ')
        print (bemt.bemt_interface.best_design)
    else:
        if(self.bemt_mode !=2):
            print ('BEMT INDICATED FAILURE TO TRIM ROTOR OR MATCH EFFICIENCY TARGETS')
    #     for i in xrange(N):
    #         segment[i].aero_eta     = 0.05 
    #     valid     = True 

#======================================================================
# end of operations
#======================================================================

    return None         # 1 for true, 0 for false