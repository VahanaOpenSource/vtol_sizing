#====================================================================
# Python program to allocate process rank array for "N" length 
# vector 
#====================================================================
import numpy
def assign_rank(Ncases, Nprocs):


    ncases      = int(Ncases/Nprocs)                        # cases per process
    pid_array   = numpy.ones(Ncases).astype(int)*(Nprocs-1)  # assign default rank to last process

#====================================================================
# Loop over case numbers
#====================================================================

    for icase in range(Ncases):

        found         = False 

#====================================================================
# loop over process rank 
#====================================================================
        
        for prank in range(Nprocs):

#====================================================================
# smallest and largest case ids covered by this process
#====================================================================

            id_start  = ncases*prank                       # starting ID
            id_end    = id_start + ncases - 1              # ending   ID 
            if (icase >= id_start and icase <= id_end):      #
                pid_array[icase] = prank 
                found            = True
            
            if found:
                break 

#====================================================================
# End of operations
#====================================================================

    return pid_array