#define type(X) Type(X), target

!==================================================================================
! data structure containing inputs for FEA
!==================================================================================

module fea_types 

   implicit none
   integer, parameter :: max_nodes      = 100
   integer, parameter :: mx_wingnode    = 10
   integer, parameter :: mx_rotors      = 20
   integer, parameter :: mx_wings       = 4

   type :: interface_def

!==================================================================================
! run-time options 
!==================================================================================

      logical         :: get_modes      = .false.

!==================================================================================
! rotor-related information
!==================================================================================

      integer         :: id_thrust, id_torque
      integer         :: nrotors, id_rotor(mx_rotors)  ! # rotors and rotor nodes
      real(kind=8)    :: dirn(mx_rotors)               ! direction of rotation (+1 = ccw, -1 = cw)
      real(kind=8)    :: Th, Qh, Tc, Qc                ! hover thrust and torque, cruise thrust and torque 

!==================================================================================
! wing-related information
!==================================================================================

      integer         :: nwings                        ! number of wings
      integer         :: nwing_elem(mx_wings)          ! number of elements in each wing
      integer         :: iw_elem(mx_wingnode,mx_wings) ! element IDs in all wings
      real(kind=8)    :: L, D                          ! wing lift and drag, newtons

!==================================================================================
! to scale up model from prototype python code
!==================================================================================

      integer         :: nnodes                        ! total # nodes in model
      integer         :: nwing_nodes    = 0            ! # nodes b/w wing tip and mount
      integer         :: norig_nodes    = 0            ! # wing tip locations
      integer         :: iw_nodes(3,mx_wingnode)       ! intermediate wing nodes
      integer         :: orig_nodes(mx_wingnode)       ! wing tip nodes 
      real(kind=8)    :: wtw_nodes(2,mx_wingnode)      ! interpolation weights 

!==================================================================================
! scaling parameters are here
!==================================================================================

      real(kind=8)    :: ref_R, new_R                  ! reference radius and new rotor radius 
      real(kind=8)    :: ref_b, new_b                  ! reference radius and new rotor radius 
      real(kind=8)    :: norm_xyz(3,max_nodes)         ! coordinates of points (normalized)

   end type

!==================================================================================
! actual data is here
!==================================================================================

    type(interface_def) :: iodata

end module