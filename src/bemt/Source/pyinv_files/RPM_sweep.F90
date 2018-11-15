!===========================================================================
! fortran routine to perform RPM sweep and obtain power, thrust prediction
! for a range of 
!===========================================================================

  subroutine RPM_sweep(flt, Rotor, k)

    use precision
    use bemt_data_structures
    use bemt_types
    implicit none

!===========================================================================
! inputs
!===========================================================================

    integer         , intent(in)        :: k

!===========================================================================
! Input/Output
!===========================================================================

    type(flt_def)   , intent(inout)     :: flt
    type(Rotor_def) , intent(inout)     :: Rotor   

!===========================================================================
! local variables
!===========================================================================

    logical           :: converged 
    integer           :: i, nvalid, j
    real(kind=rdp)    :: omega, omin, omax, domega, omega_ref, T,P

!===========================================================================
! begin executable code 
!===========================================================================

    nvalid            = 0

!===========================================================================
! calling fixed-point iteration for BEMT, Ananth ishtyle
!===========================================================================

    omega_ref         = Rotor % Omega
    omin              = 0.1d0*omega_ref
    omax              = 2.3d0*omega_ref
    domega            = (omax-omin)/real(nRPM-1)
    j                 = 0

    do i = 1,nRPM
      Omega  = omin + domega*real(i-1)

!===========================================================================
! converge for inflow OGE, find thrust and power at this pitch angle
!===========================================================================

      Rotor % Omega       = Omega
      call converge_inflow(flt, Rotor, T, P, converged)

      if(converged) then 
         j                      = j + 1
         flt % T_RPMsweep(j,k)  = T 
         flt % P_RPMsweep(j,k)  = P
         flt % RPM_sweep(j,k)   = Omega 
         flt % nvalid_RPMs(k)   = j
      end if 

!===========================================================================
! check if found low and high limits
!===========================================================================

    end do 

!restore reference rotor RPM
    Rotor % Omega     = omega_ref 

!===========================================================================
! end of operations
!===========================================================================

return
end subroutine RPM_sweep