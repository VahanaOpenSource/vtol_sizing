!==================================================================================
! routine to compute boundary condition indices and corresponding DOF in matrices
!==================================================================================

    subroutine index_bc(disp_id)

    use layout

    implicit none 

!==================================================================================
! inputs
!==================================================================================

    integer     , intent(in)    :: disp_id(ndisp,3)

!==================================================================================
! local variables
!==================================================================================

    integer                     :: i, nodeid, dofid, nbc 

!==================================================================================
! begin executable code
!==================================================================================

    do i = 1, ndisp
      nodeid           = disp_id(i,1)       ! id # of the node 
      dofid            = disp_id(i,2)       ! which DOF to constrain
      extdisp_dofid(i) = ndof_beam * (nodeid-1) + dofid     ! index in the matrix 

    end do 

!==================================================================================
! also zero out external forces
!==================================================================================

    extforce           = 0.d0

    end subroutine    

!==================================================================================
! end of operations
!==================================================================================


!==================================================================================
! routine to compute loading indices for FEA
!==================================================================================

    subroutine index_loading(il, nforces, force_id)

    use layout

    implicit none 

!==================================================================================
! inputs
!==================================================================================

    integer     , intent(in)    :: il                   ! which flight condition
    integer     , intent(in)    :: nforces              ! how many forces are applied
    real(kind=8), intent(in)    :: force_id(nforces,3)  ! node ids etc for applying loads

!==================================================================================
! local variables
!==================================================================================

    integer                     :: i, nodeid, dofid, nbc 
!    real(kind=8)                :: dload 

!==================================================================================
! begin executable code
!==================================================================================

    extforce(il,:)          = 0.d0
    do i = 1, nforces
      nodeid                = int(force_id(i,1))   ! id # of the node 
      dofid                 = int(force_id(i,2))  ! which DOF is loaded

      extforce(il,i)        = extforce(il,i) + force_id(i,3)
      extforce_dofid(il,i)  = ndof_beam * (nodeid-1) + dofid   ! index in the matrix 

!      write(*,*) il,i,extforce(il,i)
    end do 

    nforce(il)              = nforces
    end subroutine    

!==================================================================================
! end of operations
!==================================================================================