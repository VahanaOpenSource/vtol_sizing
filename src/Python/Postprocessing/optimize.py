#====================================================================
# Python program to call SCIPY's optimizer
#====================================================================

#====================================================================
# define locations and import modules
#====================================================================

import os, numpy, math, pickle, matplotlib, sys, time

sys.path.append('../Stage_0')
sys.path.append('../Stage_1')
sys.path.append('../Stage_1/afdd')
sys.path.append('../Stage_1/engines')
sys.path.append('../Stage_3')

from hydraInterface     import hydraInterface
from constraints        import *
from scipy.optimize     import minimize
from assign_rank        import assign_rank
from extract_plot_data  import write_to_file 

#====================================================================
# MPI related variables
#====================================================================

from mpi4py import MPI
comm     = MPI.COMM_WORLD               # communicator
rank     = comm.Get_rank()              # process rank
nprocs   = comm.Get_size()              # number of processes

#====================================================================
# function to read data from best_designs.dat file to identify 
# valid best designs and run optimization
#====================================================================

def optimize(path):

    t1          = time.time()
    fname       = path + '/best_design.dat' 

    with open(fname,'r') as f:
        line    = next(f)

#====================================================================
# determine file names
# take header and duplicate in optim_summary.dat
#====================================================================
    
    if(rank == 0):
        fname2  = path + '/optim_summary.dat' 
        with open(fname2,'w') as f:
            f.write(line)
    else:
        fname2  = path + '/optim_p' + str(int(rank)) + '.dat' 
        with open(fname2,'w') as f:
            print('flushed file',fname2)

#====================================================================
# get list of designs and assign process rank for each
#====================================================================

    data        = numpy.loadtxt(fname, skiprows=1)

#====================================================================
# special processing for single dimensional arrays (only 1 valid design)
#====================================================================

    sh          = numpy.shape(data)
    if(len(sh) == 1):
        designs     = data.reshape(sh[0],1)
        ndesigns    = 1
    else:
        designs     = list(data[:,0])
        ndesigns    = len(designs)

    design_id   = str(int(designs[0]))
    pids        = assign_rank(ndesigns, nprocs)

    if(rank == 0):
        print('found ',ndesigns, 'valid designs to optimize around')
    t1          = time.time()

#====================================================================
# loop over designs handled by this process, run scipy optimization
#====================================================================

    for i in range(ndesigns):
        if(rank == pids[i]):
            optimize_single(path, int(designs[i]), fname2)

    t2          = time.time()

#====================================================================
# Collocate summary_p*.dat lines into one file
#====================================================================

    comm.Barrier()

    if(rank == 0):
        print('time taken is {:12.2f} seconds'.format(t2-t1))

    if(rank == 0 and nprocs > 1):
        with open(fname2,'a') as f:
            for i in range(1,nprocs):
                fname   = path + '/optim_p' + str(i) + '.dat'
                with open(fname,'r') as origin:
                    for line in origin:
                        f.write(line)
                os.remove(fname)
        t3        = time.time()

        print ('File reconstruction time taken is ',t3-t2,'seconds')

#====================================================================
# filter out very similar designs: defined as metrics being within 1%
# To do!
#====================================================================

#====================================================================
# End of operations
#====================================================================

    return None

#====================================================================
# get design ID, determine file name to read data from
#====================================================================

def optimize_single(path, design_id, fname_summary):

    logdir      = path    + '/output/logs/'
    fname       = logdir  + '/log' + str(design_id) + '.txt'
    with open(fname,'rb') as f:
       design = pickle.load(f)
    
    design.kick_off()
    baseline          = design.get_essentials()
    nvars             = len(baseline)
    design.bemt_mode  = 1          # fix design, calculate performance
    dvars             = list(baseline.keys())
    com               = design.com 

#====================================================================
# setup initial vector of design variables
#====================================================================

    if(design.wing.ngroups > 1):
        nDV           = 9
    else:
        nDV           = 6

    x                 = numpy.zeros(nDV)
    x[0]              = com['disk_loading']
    x[1]              = com['vtip']
    x[2]              = com['solidity']
    x[3]              = design.massTakeoff
    x[4]              = com['wing_group0_aspectratio']
    x[5]              = com['wing_group0_cl']

    if(design.wing.ngroups > 1):
        x[6]          = com['wing_group1_aspectratio']
        x[7]          = com['wing_group1_cl']
        x[8]          = com['wing_group1_liftfraction']

    a                 = run_sizing(x, design)

#====================================================================
# setup and run optimization problem
#====================================================================

    bnds, cons 		  = constraint_assembler(design)
    try:
        options           = {'eps': 1e-7}
        solution 	      = minimize(cost_function,x,args=(design,), 	\
                            method='SLSQP',bounds=bnds,constraints=cons,options=options)

#        print(solution)
        xsoln             = solution.x;

        a                 = run_sizing(xsoln, design)
#        for con in cons:
#            print('constraint value = ',con['fun'](xsoln,design))
#        print(a)
        if(a['valid']):
            print('Design ' + str(design_id) + ' Cost = ',round(a['Cost'],2))
            fname             = logdir  + '/optim_log' + str(design_id) + '.yaml'
            design.writelogdata(fname)
        else:
            print('Design ' + str(design_id) + ' Optimization failed to identify valid solution!' + '\n' + design.errmsg)
    except:
        print('encountered error')
        pass

    with open(fname_summary,'a') as f:
        write_to_file(design.com, design.summary, f)

    return None