#====================================================================
# Routine to iterate multiple designs 
# Accommodates changes in AR, nblade, tip speed etc.
#====================================================================

#====================================================================
# Access fortran routines
#====================================================================

import numpy as np

from   mpi4py import MPI
comm     = MPI.COMM_WORLD
rank     = comm.Get_rank()
nprocs   = comm.Get_size()

#====================================================================
# Begin function
#====================================================================

class _iterate_design:

   def iterate_design(self,blade, run_CSD):
            
#====================================================================
# iterate until convergence (or max steps)
#====================================================================
   
      mission           = self.mission
      sizing_mode       = mission.sizing_mode
      Mmax              = self.constraints.max_gtow

#====================================================================
# for sizing mode = 2 (fixed weight, variable payload) nothing special
# use fixed point iterations and continuously check convergence
#====================================================================

      if (sizing_mode == 2):
         self.converge_mode_2(blade, run_CSD)

#====================================================================
# fixed payload mode: change things slightly
# run in fixed take off weight mode, calculate payload and adjust 
# takeoff weight to get the right payload
#====================================================================

      else: 

#====================================================================
# set initial weight (guideline) and get corresponding payload
# 
#====================================================================
         
         try:
            M0                =  self.mission.guess_weight
         except:
            M0                = 6.0*mission.payload_tar     

         dMdP                 = 2.0         # inverse of payload fraction for iterations
         err                  = 100 
         P0                   = self.get_payload(M0, blade, run_CSD)       # get payload
         Ptar                 = mission.payload_tar
         err                  = abs(Ptar - P0)
         iterc                = 1
         failed               = False 
         while(err > 0.05):
            deltaP            = Ptar - P0 
            M1                = M0 + dMdP*deltaP
#            print(Ptar,P0,deltaP)
#            print('adjusting weight','slope dMdP',dMdP,'new M = ',M0,'increment',dMdP*deltaP)
            P1                = self.get_payload(M1, blade, run_CSD)

            dMdP              = (M1-M0)/(P1-P0) 
            dMdP              = dMdP if dMdP > 0 else 2.5

            err               = abs(Ptar - P1)
            M0                = M1 
            P0                = P1
            iterc             = iterc + 1

#====================================================================
# Exceeded # iterations or weight ballooning triggered: halt iterations
#====================================================================

            if(iterc > 100 or M0 > 1.05*Mmax):
               failed         = True 
               err            = 0.0 
               self.err_msg   = 'could not converge weight with new method'

#====================================================================
# if convergence failed, try the old way
#====================================================================

         if(failed):
            print('faster convergence failed; trying original method')
            mission.payload   = mission.payload_tar 
            self.converge_mode_2(blade,run_CSD)
#         print(Ptar,P0, 'iterations = ',iterc);quit()

#====================================================================
# End of operations
#====================================================================

      return None 

#====================================================================
# This function obtains a converged payload for a given take-off 
# weight
#====================================================================

   def converge_mode_2(self, blade, run_CSD):

      itercount         = 0 
      converged         = False 
      while not converged:
         itercount   = itercount + 1          
         converged   = self.march_one_step(itercount)

#      print('performed # iterations = ',itercount)
      return None

#====================================================================
# calculate payload for a given take-off weight
#====================================================================

   def get_payload(self, M, blade, run_CSD):

#====================================================================
# create shortcuts, remmeber machine state
#====================================================================

      mission                 = self.mission
      smode                   = mission.sizing_mode

#====================================================================
# set fixed take off weight mode and run sizing
#====================================================================

      mission.sizing_mode     =  2 
      mission.gtow_target     = M
      self.massTakeoff        = M# mission.mass_takeoff_guess(M)
      self.converge_mode_2(blade, run_CSD)

#====================================================================
# restore machine state, return payload
#====================================================================

      mission.sizing_mode     = smode 
      return mission.payload 

