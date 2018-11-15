!================================================================================
! preprocessor to setup Finite-element analysis
! create a map to represent all the elements of a particular member
!================================================================================

subroutine preprocess

   use layout
   use material

   implicit none

   integer ::  j, k, kk, memID

!================================================================================
! begin executable code
!================================================================================

   nmember = 0
   do j = 1,nelem
      nmember = max(nmember,memberID(j))
   enddo
!
   allocate(nelem2member(nmember))
!
   nelem2member = 0
   do j = 1,nelem
      memID = memberID(j)
      nelem2member(memID) = nelem2member(memID) + 1
   enddo
   !
   maxnelem2member = 0
   do j = 1,nmember
      maxnelem2member = max(maxnelem2member, nelem2member(j))
   enddo
   !
   allocate(elem2member(maxnelem2member,nmember))
   elem2member  = -1
   !
   nelem2member = 0
   do j = 1,nelem
      memID = memberID(j)
      
      nelem2member(memID)                    = nelem2member(memID) + 1
      elem2member(nelem2member(memID),memID) = j

   enddo
   
end subroutine preprocess

!================================================================================
! find rows to eliminate and map from reduced rows to original rows
!================================================================================

subroutine assign_indices
   use layout
   use currentvalues
   implicit none 

!================================================================================
! local variables
!================================================================================

   integer     :: kk, j, jj

!================================================================================
! begin executable code
!================================================================================

   kk       = 1                  ! running counter
   bc_flag  = 0                  ! 1 = boundary condition, 0 = free node
   do j = 1,ndof                 ! all rows in original matrix

      bc_flag(j)           = kk  ! maps original matrix location -> reduced matrix
      do jj = 1,ndisp            ! all BC
         if (extdisp_dofid(jj) == j) then    ! if current DOF shd be eliminated
            bc_flag(j)     = 0
            exit                 ! exit this loop if boundary condition identified
         endif
      enddo
!
      if(bc_flag(j) .ne. 0) then      ! if current dof shd be preserved=
        row_red(kk) = j             ! map index in reduced matrix to original matrix 
         kk         = kk + 1        ! increment red matrix index
      end if
   enddo
   
!================================================================================
! end of operations
!================================================================================

end subroutine assign_indices