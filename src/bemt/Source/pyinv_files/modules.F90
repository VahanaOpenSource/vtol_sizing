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
! data structure containing design variables (use to generate cases)
!====================================================================

    type :: dvars_def
        integer         :: nOmega
        real(kind=rdp)  :: taper(ntaper), thx(nthx) 
        real(kind=rdp)  :: omega(mxOmega), thtip(nthtip), x(nx)
    end type

end module
