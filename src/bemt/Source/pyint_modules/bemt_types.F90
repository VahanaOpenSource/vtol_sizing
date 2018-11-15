!====================================================================
! This module contains the definitions of all the data structures
! that need access from Python
!====================================================================

#define type(x) Type(x), target
module bemt_types

    use precision

!====================================================================
! flight condition
!====================================================================

    type :: flt_def
        real(kind=rdp)      :: thrust, Vc, time, rho, altitude          ! inputs
        real(kind=rdp)      :: power, omega, coll, efficiency           ! outputs

!for RPM sweep with fixed collective
        integer             :: nvalid_RPMs(2)
        real(kind=rdp)      ::  RPM_sweep(NRPM,2)
        real(kind=rdp)      :: T_RPMsweep(NRPM,2), P_RPMsweep(NRPM,2)   ! thrust, power for RPM sweep
    end type

!====================================================================
! one design case, performance at all flight conditions given 
!====================================================================

    type :: case_def
        logical             :: valid = .false.
        real(kind=rdp)      :: taper, thx, thtip, x, omega
        real(kind=rdp)      :: energy = zero, FMmin = one, etamin = one
    end type

!====================================================================
! numbers that flow from sizing to bemt and vice versa: Nb
! & flight conditions, radius, thrust/speed/time in each condition
!====================================================================

    type :: inter_def
        logical             :: use_tables = .false.     ! airfoil tables
        logical             :: sweep_rpm  = .true.      ! switch: rpm sweep (cruise)
        integer             :: mode_operation=0         ! default mode
        real(kind=rdp)      :: CL_alpha = 5.73d0        ! linear airfoil properties
        real(kind=rdp)      :: Cdo = 0.01d0
        integer             :: Nb, N_flight
        real(kind=rdp)      :: Radius, solidity, omegaH
        real(kind=rdp)      :: rpm_ratio = 0.5d0        ! cruise to hover RPM ratio
    end type 

end module

!====================================================================
!actual data is here
!====================================================================

module bemt_interface
    use bemt_types
    type(inter_def)     :: Inputs
    type(case_def)      :: best_design
    type(flt_def)       :: flt(nfc)
end module
