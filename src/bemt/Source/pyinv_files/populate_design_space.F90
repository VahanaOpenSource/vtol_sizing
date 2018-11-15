!===========================================================================
!> fortran routine to populate the design space for BEMT parametric study
!===========================================================================

subroutine populate_design_space(inputs, Rotor, dvars, nd)

    use bemt_data_structures
    use bemt_types
    implicit none

!===========================================================================
! input/output
!===========================================================================

    type(inter_def)     , intent(in)      :: inputs
    type(Rotor_def)     , intent(inout)   :: Rotor   
    type(dvars_def)     , intent(inout)   :: dvars 

!===========================================================================
! output: total # of designs
!===========================================================================

    integer             , intent(out)     :: nd

!===========================================================================
! local variables
!===========================================================================

    integer             :: i, N
    real(kind=rdp)      :: dr, Tmax
    real(kind=rdp)      :: max_omega, min_omega, domega, vmin, vmax, dv

!===========================================================================
! begin executable code
!===========================================================================

    N                 = inputs % N_flight

!===========================================================================
! transfer geometry-related inputs to BEMT data structures
!===========================================================================

    Rotor % radius     = inputs % radius 
    Rotor % Nb         = inputs % Nb
    Rotor % solidity   = inputs % solidity

!===========================================================================
! generate spanwise stations
!===========================================================================

    dr          = (one - Rotor % rcout)/real(Nseg)
    do i = 1, Nseg
       Rotor % r(i)    = Rotor % rcout + (real(i) - half)*dr 
    end do 

!===========================================================================
! populate design variables
!===========================================================================

!===========================================================================
! taper: from 1 to 3
!===========================================================================

    vmin    = one 
    vmax    = three 
    dv      = (vmax-vmin)/real(ntaper-1)
    do i = 1,ntaper
       dvars % taper(i) = vmin + real(i-1)*dv
    end do 

!===========================================================================
! twist at bilinear junction: 5 to 20 deg (defined +ve for nose down)
!===========================================================================

    vmin    = 8.0d0
    vmax    = 20.0d0 
    dv      = (vmax-vmin)/real(nthx-1)
    do i = 1,nthx
       dvars % thx(i)   = vmin + real(i-1)*dv
    end do 

!===========================================================================
! twist at blade tip: 10 to 35 deg (defined positive nose down)
!===========================================================================

    vmin    = 10.0d0 
    vmax    = 30.0d0
    dv      = (vmax-vmin)/real(nthtip-1)
    do i = 1,nthtip
       dvars % thtip(i) = vmin + real(i-1)*dv
    end do 

!===========================================================================
! bilinear twist junction location: 0.4 to 0.7 of span
!===========================================================================

    vmin    = 0.4d0 
    vmax    = 0.7d0 
    dv      = (vmax - vmin)/real(nx-1)
    do i = 1, nx
       dvars % x(i)     = vmin + real(i-1)*dv
    end do 

!===========================================================================
! If cruise RPM is unknown, then sweep from 100% -> 50% of hover RPM 
!===========================================================================

    if(inputs % sweep_RPM) then
        max_omega   = 1.0d0 * inputs % omegaH          
        min_omega   = 0.5d0 * max_omega
        domega      = (max_omega - min_omega)/real(mxomega-1)

        do i = 1, mxomega
           dvars % omega(i) = min_omega + real(i-1)*domega
        end do 
        dvars % nOmega      = mxOmega

!===========================================================================
! cruise RPM known: run only one RPM for geometry sweep
!===========================================================================

    else 
        dvars % omega(1)    = inputs % omegaH * inputs % rpm_ratio 
        dvars % nOmega      = 1 
    end if 

!===========================================================================
! find total # of geometries to evaluate per flight condition
!===========================================================================

    nd          = ntaper * nthx * nthtip * nx

!===========================================================================
! end of operations
!===========================================================================

return
end subroutine