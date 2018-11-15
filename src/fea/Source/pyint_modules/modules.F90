!==================================================================================
! list of modules
!==================================================================================

module layout

   implicit none

   integer, parameter        :: ndof_beam = 6
   integer, parameter        :: nmodes = 8
   integer                   :: nnodes,nelem,nmoment,ndisp,ndof,nloading,maxnforce
   integer, allocatable      :: extforce_dofid(:,:),extdisp_dofid(:)
   real(kind=8), allocatable :: nodepos(:,:)
   integer, allocatable      :: conn(:,:), memberID(:),nmember
   integer, allocatable      :: nforce(:)
   real(kind=8), allocatable :: extforce(:,:),extdisp(:)
   integer, allocatable      :: maxnelem2member,nelem2member(:),elem2member(:,:)
contains

!==================================================================================
! initializes FEA layout
!==================================================================================

   subroutine initializeLayout()

   implicit none
      !
      allocate(nodepos(3,nnodes))
      !
      ndof = ndof_beam * nnodes
      !
      allocate(conn(2,nelem),memberID(nelem))
      !
      allocate(nforce(nloading))
      allocate(extforce(nloading,maxnforce))
      allocate(extforce_dofid(nloading,maxnforce))
      !
      allocate(extdisp_dofid(ndisp),extdisp(ndisp))

   end subroutine initializeLayout

end module layout

!==================================================================================
! material properties
!==================================================================================

module material

   implicit none

   integer, parameter        :: nmatprop=12
   real(kind=8), allocatable :: matprop(:,:)
   character(100)            :: cstype

contains

   subroutine initializeMaterial(nelem)

      implicit none
      integer, intent(in) :: nelem

      allocate(matprop(nmatprop,nelem))
      matprop = 0.d0

   end subroutine initializeMaterial

end module material

!==================================================================================
! values stored at every iteration
!==================================================================================
!
module currentvalues

   implicit none

   integer                   :: ndof_red,half_ndof_beam
!   integer, allocatable      :: nodemap_red(:)
   integer, allocatable      :: row_red(:), bc_flag(:)  
   real(kind=8), allocatable :: K_global(:,:), K_red(:,:)  ! stiffness matrix
   real(kind=8), allocatable :: M_global(:,:), M_red(:,:)  ! mass matrix
   real(kind=8), allocatable :: rhs_vec(:,:),  rhs_red(:,:)! rhs vector
   real(kind=8), allocatable :: xdof(:,:)!,   xdof_red(:,:) ! DOF vector
   !
   ! stresses and strains
   !
   real(kind=8),parameter    :: fBuckling=5.d0
   real(kind=8)              :: maxBuckling
   real(kind=8), allocatable :: strain(:)                  ! strain in the elements
   real(kind=8), allocatable :: stress(:), stressFOS(:)    ! stress in the elements
   real(kind=8), allocatable :: bucklingFOS(:)             ! buckling FOS
   real(kind=8), allocatable :: weightElem(:), totalWeight ! weight 
   real(kind=8), allocatable :: memStress(:), memFOS(:), memDeflection(:)
   real(kind=8)              :: minFOS, maxDeflection
   real(kind=8), allocatable :: deflection(:)

!==================================================================================
! optimizer values
!==================================================================================

   integer, parameter         :: nmaxopt = 100
   integer                    :: optcount
   real(kind=8), allocatable  :: minOptWeight(:), minOptStressFOS(:),maxOptDeflection(:)
   real(kind=8)               :: maxAllowableDeflection
   integer,allocatable        :: idefFlag(:),istressFlag(:)

!==================================================================================
! natural response
!==================================================================================

   real(kind=8),allocatable   :: modes(:,:),freqs(:)
   real(kind=8),allocatable   :: modeShapes(:,:),modeFreqs(:)

contains

!==================================================================================
! to initialize calculations
!==================================================================================

   subroutine initializecvals
      use layout, only : nmember,ndof,nnodes,nmodes,nelem,nloading,ndisp,ndof_beam
      implicit none

      ndof_red       = ndof-ndisp
      half_ndof_beam = 0.5*ndof_beam

!==================================================================================
! from optimizeStructure
!==================================================================================

      allocate(minOptWeight(nmaxopt),minOptStressFOS(nmaxopt))
      allocate(maxOptDeflection(nmaxopt))
      allocate(idefFlag(nmember))
      allocate(istressFlag(nmember))

!==================================================================================
! apply bc
!==================================================================================

      allocate(rhs_vec(ndof,nloading))
      allocate(K_red(ndof_red,ndof_red),rhs_red(ndof_red,nloading))
      allocate(M_red(ndof_red,ndof_red))

!==================================================================================
! assemble
!==================================================================================

      allocate(K_global(ndof,ndof),M_global(ndof,ndof))

!==================================================================================
! calculate loads
!==================================================================================

      allocate(stress(nelem),strain(nelem),bucklingFOS(nelem))
      allocate(stressFOS(nelem),weightElem(nelem))
      allocate(deflection(nnodes))

!==================================================================================
! reset loads
!==================================================================================

      allocate(xdof(nloading,ndof))
      allocate(memStress(nmember))
      allocate(memFOS(nmember))
      allocate(memDeflection(nmember))

!==================================================================================
! extractnaturalresponse
!==================================================================================

      allocate(modeShapes(nmodes,half_ndof_beam*nnodes),modeFreqs(nmodes))

!==================================================================================
! indices to map from reduced to full stiffness matrix
!==================================================================================

      allocate(row_red(ndof_red))
      allocate(bc_flag(ndof))

!==================================================================================
! set normalized deflection constraint
!==================================================================================

      maxAllowableDeflection = 0.1d0            ! default value 

   end subroutine initializecvals

end module currentvalues
