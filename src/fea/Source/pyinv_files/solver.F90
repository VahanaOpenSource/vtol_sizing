!================================================================================
!
! solver for the sparse matrix
!
!================================================================================

subroutine solver(failure)

   use layout
   use currentvalues


   implicit none

   logical, intent(inout)  :: failure
!   integer, intent(in)  :: n 
   integer        :: j!,k,nodeid,rowid, iter
   integer        :: info, loadingid, ipiv(ndof_red)
!   real(kind=8)   :: disp, work(n*nloading)
!   real(kind=4)   :: swork(n*(n+nloading))

!================================================================================
! begin executable code 
!================================================================================

!================================================================================
! simultaneous leq solver ==> DA BEST
!================================================================================

   call DGESV( ndof_red, nloading, K_red, ndof_red, ipiv, rhs_red, ndof_red, info)
   if(info .ne. 0) then
      write(*,*) 'solver.F90: matrix inversion failed. info=0'
      failure = .true.
      return
   endif

!================================================================================
! fill the full xdof array
!================================================================================

   do loadingid = 1,nloading
      do j = 1,ndisp
         xdof(loadingid,extdisp_dofid(j)) = extdisp(j)
      enddo
   !
      do j = 1,ndof_red
         xdof(loadingid,row_red(j)) = rhs_red(j,loadingid)
      enddo
   end do 

end subroutine solver