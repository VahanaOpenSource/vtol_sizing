!
! main program
!
! bharath govindarajan
!
program fea

   implicit none

   integer, parameter :: ioptimize = 0

   ! read in inputs
   call readInputs

   ! preprocess
   call preprocess

   ! iteratively solve for optimal structure (if required)
   call optimizeStructure(ioptimize)
   
   ! write outputs
   call writeOutputs

end program fea
