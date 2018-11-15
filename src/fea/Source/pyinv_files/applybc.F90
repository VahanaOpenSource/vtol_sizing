!================================================================================
! apply essential bc to stiffness matrix
! also set up vectors for right hand side (loading) for all cases
!================================================================================
subroutine applybc()
   
   use layout
   use currentvalues
   
   implicit none
   
   integer        :: j,jj,i

!================================================================================
! begin executable code
!================================================================================

!================================================================================
! create initial force vector (fx,fy,fz,mx,my,mz for each node)
! done by set_loads routine in pyint_files/one_time_transfer.F90 
!================================================================================

!   rhs_vec = 0.d0
!   do i = 1, nloading
!      do j = 1,nforce(i)
!         rhs_vec(extforce_dofid(i,j),i) = extforce(i,j)
!      end do
!   end do  

!================================================================================
! find reduced rhs_vec
!transferred to set_loads routine in pyint_files/one_time_transfer.F90 
!================================================================================

   do i = 1, nloading
      do jj = 1,ndof_red                            ! loop over free DOF
         j                 = row_red(jj)
         rhs_red(jj,i)     = rhs_vec(j,i)
      enddo
   end do

!================================================================================
! end of operations
!================================================================================

end subroutine applybc