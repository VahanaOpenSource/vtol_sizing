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
from scipy.optimize     import differential_evolution as diffl_evo
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
# parse design vars and cost data to find unique designs
#====================================================================

def parse_designs(fnames, design):

    data    = numpy.loadtxt(fnames['dvars'])

#====================================================================
# first column is design id, next columns are DVs, final column is cost
#====================================================================

    ids     = data[:,0].astype(int)
    XandC   = data[:,1:]

    nr,nc   = numpy.shape(XandC)
    nX      = nc-1                  #number of dvars = # cols - 2 

#====================================================================
# loop over designs: for each row, check all other rows
#====================================================================
    
    uflag   = numpy.ones((nr,1),dtype=bool)        # indicate if unique

    for i in range(nr-1):
        if(uflag[i]):
            print('checking against unique row ',i)
            xi      = XandC[i,:]
            for j in range(i+1,nr):
                if(uflag[j]):
                    xj  = XandC[j,:]
                    dx  = numpy.amax(numpy.divide(xj-xi,xi))
                    if(dx < 0.02):
                        uflag[j]    = False 
                        print('row # ',j,'is not unique')
    uids    = [i for i,flag in enumerate(uflag) if flag]

#====================================================================
# write identified unique designs to file 
#====================================================================

    with open(fnames['summary'],'r') as f:
        line    = next(f)

    with open(fnames['summary'],'w') as f:
        f.write(line)

#====================================================================
# run sizing one last time
#====================================================================

    with open(fnames['dvars'],'w') as dvars, open(fnames['summary'],'a') as summary:
        for uid in uids:

            x                   = XandC[uid,:nX]
            design_id           = str(ids[uid])
            design.com['icom']  = ids[uid]

#====================================================================
# using the row ids, run the designs one last time, store details
#====================================================================

            a               = run_sizing(x, design)

            if(a[0] == 1):
                print('{:8s} {:6s} {:12.2f} USD/hr'.format('Design # ',design_id, round(a[6],2)))
                fname             = 'output/logs/optim_log' + design_id + '.yaml'
                design.writelogdata(fname)

#write design variables and summary
                write_dvars(dvars, design_id, x, XandC[uid,-1])
                write_to_file(design.com, design.summary, summary)

    return uids         # these are unique row identifiers!!

#====================================================================
# function to read data from best_designs.dat file to identify 
# valid best designs and run optimization
#====================================================================

def optimize(path,method=None):

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
        fname3  = path + '/optim_dvars_summary.dat' 
        with open(fname2,'w') as f:
            f.write(line)
        with open(fname3,'w') as f:
            print('flushed file ',fname3)
    else:
        fname2  = path + '/optim_p' + str(int(rank)) + '.dat' 
        fname3  = path + '/optim_dvars_p' + str(int(rank)) + '.dat' 
        with open(fname2,'w') as f:
            print('flushed file ',fname2)
        with open(fname3,'w') as f:
            print('flushed file ',fname3)

    fnames      = {'dvars': fname3, 'summary': fname2}

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

    nmax        = int(max(designs))+1
    design_id   = str(int(designs[0]))
    pids        = assign_rank(ndesigns, nprocs)

    if(rank == 0):
        print('found ',ndesigns, 'valid designs to optimize around')
    t1          = time.time()

#====================================================================
# set number of designs to investigate: 
# all for gradient-based 
# one for differential evolution
#====================================================================

    if(rank == 0):
        design  = optimize_single(path, design_id, fnames, 'diffl_evo', str(nmax))

#====================================================================
# loop over designs handled by this process, run scipy optimization
#====================================================================

    for i in range(ndesigns):
        if(rank == pids[i]):
            design = optimize_single(path, str(int(designs[i])), fnames, 'gradient_based', str(int(designs[i])))

    t2          = time.time()

#====================================================================
# Collocate summary_p*.dat lines into one file
#====================================================================

    comm.Barrier()

    if(rank == 0):
        print('time taken is {:12.2f} seconds'.format(t2-t1))

    if(rank == 0 and nprocs > 1):

        t2          = time.time()

        with open(fname2,'a') as f:
            for i in range(1,nprocs):
                fname   = path + '/optim_p' + str(i) + '.dat'
                with open(fname,'r') as origin:
                    for line in origin:
                        f.write(line)
                os.remove(fname)

        with open(fname3,'a') as f:
            for i in range(1,nprocs):
                fname   = path + '/optim_dvars_p' + str(i) + '.dat'
                with open(fname,'r') as origin:
                    for line in origin:
                        f.write(line)
                os.remove(fname)

        t3        = time.time()

        print ('File reconstruction time is ',t3-t2,'seconds')

#====================================================================
# filter out very similar designs: defined when max DV and cost 
# differences are within 0.5% 
#====================================================================

    if(rank == 0):
        print('parsing designs for unique values')
        uids = parse_designs(fnames, design)

    return None

#====================================================================
# load design from file
#====================================================================

def load_design(path, design_id):

    logdir      = path    + '/output/logs/'
    fname       = logdir  + '/log' + design_id + '.txt'
    with open(fname,'rb') as f:
       design = pickle.load(f)
    
    return design

#====================================================================
# get design ID, determine file name to read data from
#====================================================================

def optimize_single(path, design_id, fnames, method, out_id):

    design      = load_design(path, design_id)
    
    design.kick_off()
    baseline          = design.get_essentials()
    nvars             = len(baseline)
    design.bemt_mode  = 1          # fix design, calculate performance
    com               = design.com 

#====================================================================
# iterate over "com" dictionary ==> this has the design variables
#====================================================================

    collected_keys    = []
    for k,v in sorted(com.items()):
        # k             = k1.lower()
        # print(k)
        if(k.startswith('rotor')):
            if(not(k.endswith('flap_freq') or k.endswith('Nb') or k.endswith('cruise_rpm_ratio') )):
                collected_keys.append(k)
#                if(rank == 0):
#                    print('found a rotor design variable',k)
        elif(k.startswith('wing')):
            if(not(k.endswith('nrotors') or k.endswith('nwing') or k.endswith('cruise_rpm_ratio') )):
                if(design.wing.ngroups == 1 and k.endswith('liftfraction')):
                    if(rank == 0):
                        print('NOT ADDING WING LIFT FRACTION FOR SINGLE WING GROUP CONFIGURATIONS')
                else:
#                    if(k.endswith('liftfraction')):
#                        if(rank == 0):
#                            print('multiple wings: not using lift fraction as DV')
#                    else:
                         collected_keys.append(k)
#                        if(rank == 0):
#                           print('found a wing  design variable',k) 
 
    design.var_keys     = collected_keys
    nDV                 = len(design.var_keys)+1
    bnds, cons          = constraint_assembler(design)
    design.bounds       = bnds 
    design.nDV          = nDV

#====================================================================
# setup initial vector of design variables
#====================================================================

    x                 = numpy.zeros(nDV)
    x[0]              = design.massTakeoff

    for i in range(1,nDV):
        k             = design.var_keys[i-1]
        x[i]          = com[k]

#====================================================================
# perform sizing for baseline design
#====================================================================

    a                 = run_sizing(x, design)

#====================================================================
# gradient based optimizer
#====================================================================

    if method == 'gradient_based':
        options           = {'eps': 1e-7}
        solution 	      = minimize(cost_function,x,args=(design,), 	\
                        method='SLSQP',bounds=bnds,constraints=cons,options=options)

#====================================================================
# differential evolution
#====================================================================

    else:
        options           = {'eps': 1e-7}
        print('CHOSEN DIFFERENTIAL EVOLUTION BASED OPTIMIZER')
        solution          = diffl_evo(cost_function,bnds,args=(design,),    \
                            strategy='best1bin',popsize=20,seed=2,disp=True,\
                            maxiter=200)

#====================================================================
# postprocess results
#====================================================================

    xsoln             = solution.x;
    a                 = run_sizing(xsoln, design)

#====================================================================
# for valid designs, write log data
#====================================================================

    if(a[0] == 1):

        x              = xsoln

#====================================================================
# print design variables, starting design ID and cost to file
#====================================================================
            
        with open(fnames['dvars'],'a') as f:
            write_dvars(f, out_id, x,a[6])

# print summary message to screen for process rank 0
        if(rank == 0):
            print_dvars(x,design)

#====================================================================
# print failure message
#====================================================================

    else:
        print('Design ' + out_id + ' Optimization failed to identify valid solution!' + '\n' + design.errmsg)

    return design

#====================================================================
# function to write design variables to file
#====================================================================

def write_dvars(f, design_id, x, Cost):
    f.write("   {:8s} "    .format(design_id))

    nDV         = len(x)
    for i in range(nDV):
        f.write("   {:15.6f} "    .format(x[i]))

    f.write("   {:15.6f} \n"    .format(Cost))

    return None