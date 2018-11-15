!====================================================================
! This module contains the definitions of all the data structures
! that are private to bemt analysis, i.e. don't need python wrapping 
! note: give all parameters in consistent units!
!====================================================================

!====================================================================
! Definition of data structures  
!====================================================================

module bemt_data_structures

    use precision

!====================================================================
! rotor geometry details (use in calculations)
!====================================================================

    type :: rotor_def
        logical                         :: use_tables = .true.
        integer                         :: Nb                ! get from inputs
        integer                         :: Nomega = mxOmega
        real(kind=rdp)                  :: radius, solidity  ! get from inputs 
        real(kind=rdp)                  :: omega, th0, power
        real(kind=rdp)                  :: rcout = 0.1d0     ! root cut out/radius, nondiml
        real(kind=rdp)                  :: CL_alpha
        real(kind=rdp)                  :: Cdo
        real(kind=rdp)                  :: th_75
        real(kind=rdp), dimension(Nseg) :: r, chord, twist   ! calculate from design
        real(kind=rdp)                  :: eta
    end type

!====================================================================
! data structure containing design variables (use to generate cases)
!====================================================================

    type :: dvars_def
        integer         :: nOmega
        real(kind=rdp)  :: taper(ntaper), thx(nthx) 
        real(kind=rdp)  :: omega(mxOmega), thtip(nthtip), x(nx)
    end type

end module
