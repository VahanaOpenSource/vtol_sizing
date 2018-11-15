!===========================================================================
! fortran routine to perform collective sweep and obtain power prediction
! for a given thrust & RPM using BEMT
!===========================================================================

subroutine trim_thrust(flt, Rotor, found_soln)

    use precision
    use bemt_data_structures
    use bemt_types
    implicit none

!===========================================================================
! inputs
!===========================================================================

    type(flt_def)   , intent(in)        :: flt

!===========================================================================
! Input/Output
!===========================================================================

    logical         , intent(out)       :: found_soln
    type(Rotor_def) , intent(inout)     :: Rotor   

!===========================================================================
! local variables
!===========================================================================

    logical           :: done, converged
    integer           :: i, iter, strikes
    real(kind=rdp)    :: flt_T, dT, dth, T, P, A, Vtip, R, sigma, Pideal
    real(kind=rdp)    :: dCTdth, CT, dCT, Omega, old_T, new_T, old_th0, new_th0

!===========================================================================
! begin executable code 
!===========================================================================

!===========================================================================
! try simple solution first: set initial guess for collective 
!===========================================================================

    Rotor % Power     = 1.d16                       ! some value
    Rotor % th0       = -90.0d0
    R                 = Rotor % radius 
    Omega             = Rotor % omega 
    vtip              = Omega * R 
    A                 = pi * R * R

!===========================================================================
! from CT target, generate guess for collective (in radians)
!===========================================================================

    flt_T             = flt % Thrust 
    CT                = flt_T/(flt % rho * A * Vtip * Vtip )
    sigma             = Rotor % solidity

    Rotor % th0       = 6.0d0*CT/(sigma*Rotor % CL_alpha) +                 &
                        3.0d0*sqrt(CT*half)*half - Rotor % th_75*d2r
    Rotor % th0       = Rotor % th0 * r2d  
    dCTdth            = one/(6.0d0/(sigma*Rotor % CL_alpha) + 0.25d0/sqrt(two*CT))

!===========================================================================
! calling fixed-point iteration for BEMT, Ananth ishtyle
!===========================================================================

    found_soln        = .false. ; iter = 0; strikes = 0
    done              = .false.
    new_T             = zero; new_th0  = zero 
    do while (.not. done)

      iter            = iter + 1            

!===========================================================================
! remember previous solution
!===========================================================================

      if(iter .ge. 2) then 
          old_T       = new_T
          old_th0     = new_th0 
      end if 

!===========================================================================
! converge for inflow OGE, find thrust and power at this pitch angle
!===========================================================================

      call converge_inflow(flt, Rotor, T, P, converged)

!===========================================================================
! remember solution for comparison later
!===========================================================================

      new_T           = T
      new_th0         = Rotor % th0 

!===========================================================================
! find collective update
!===========================================================================

      dT              = flt_T - T 
      dCT             = dT/(flt % rho* A * Vtip * Vtip)
      dth             = dCT/dCTdth

      dth             = min(dth, 5.0d0)
      dth             = max(dth,-5.0d0)

      Rotor % th0     = Rotor % th0 + 0.8d0*dth*r2d 

!===========================================================================
! check for convergence
!===========================================================================

      if(dT .le. 0.005d0*flt_T) then 
        found_soln    = .true.
        done          = .true.
        Rotor % power = P 
        found_soln    = .true.
      end if 

!===========================================================================
! if thrust starts decreasing even with increase in collective, go to 
! collective sweep 
!===========================================================================

      if(iter .ge. 2) then
         if((new_T - old_T)*(new_th0 - old_th0) .lt. zero) then 
            strikes       = strikes + 1
!            write(*,*) new_T, old_T, new_th0, old_th0
!            write(*,*) sign(one, new_T-old_T), sign(one, new_th0 - old_th0)
         else 
            strikes       = 0
         end if 
      end if 

!===========================================================================
! if this trimming method doesnt work, try collective sweep (brute force)
!===========================================================================

      if (iter .ge. 50 .or. Rotor % th0 .ge. 85.0d0 .or. strikes .ge. 10) then 
!        write(*,*) 'triggered error: reverting to collective sweep'
        done        = .true.
!        write(*,*) 'used ',strikes,' # of strikes'
      end if 

    end do 

!===========================================================================
! if not converged in this trimming mode mode, use collective sweep 
!===========================================================================

    if (.not. found_soln) then
        call collective_sweep(flt, Rotor, found_soln)
    end if 

!===========================================================================
! calculate ideal power
!===========================================================================

    if(flt % Vc .ge. 1.d-3) then 
       Pideal       = flt_T * flt % Vc
    else 
       Pideal       = (flt_T ** 1.5d0) / sqrt(two*flt % rho * A)
    end if 
    Rotor % eta     = Pideal/Rotor % power

!===========================================================================
! end of operations
!===========================================================================

return
end subroutine