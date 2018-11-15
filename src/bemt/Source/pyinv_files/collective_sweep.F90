!===========================================================================
! fortran routine to perform collective sweep and obtain power prediction
! for a given thrust & RPM using BEMT
!===========================================================================

subroutine collective_sweep(flt, Rotor, found_soln)

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

    logical           :: found_LB, found_UB
    integer           :: Nmax, i 
    integer, parameter:: ncoll = 46
    real(kind=rdp)    :: th0(ncoll), Thrust(ncoll), Power(ncoll)
    real(kind=rdp)    :: flt_T, dP, dT, dth, dT1, dT2
    logical           :: converged 

!===========================================================================
! begin executable code 
!===========================================================================

    Power             = 1.d16                       ! some value
    Rotor % th0       = -90.0d0
    found_LB          = .false.
    found_UB          = .false.
    Nmax              = 1

!===========================================================================
! remember target thrust
!===========================================================================

    flt_T     = flt % Thrust 

!===========================================================================
! calling fixed-point iteration for BEMT, Ananth ishtyle
!===========================================================================

    dth       = 90.0d0/real(ncoll-1)
    do i = 1,ncoll
      th0(i)  = real(i-1)*dth
    end do 

    do i = 1,ncoll
      Rotor % th0     = th0(i)                  ! set collective pitch in degrees
!===========================================================================
! converge for inflow OGE, find thrust and power at this pitch angle
!===========================================================================

      call converge_inflow(flt, Rotor, Thrust(i), Power(i), converged)

!===========================================================================
! check if found low and high limits
!===========================================================================

      if (Thrust(i) .ge. flt_T .and. converged) then 
        found_UB = .true.
      end if 

      if (Thrust(i) .lt. flt_T .and. converged) then
        found_LB  = .true.
      end if

      if (found_LB .and. found_UB) then
        Nmax = i
        found_soln = .true.
        !write(*,*) 'found limits; interpolating',i
        exit
      end if

      !write(*,*) i,th0(i),thrust(i)
    end do 

!===========================================================================
! compute collective and power required by interpolation
!===========================================================================

  if (found_lB .and. found_UB .and. Nmax .ge. 2) then

    do i = 1, Nmax-1
       dT1    = Thrust(i)   - flt_T
       dT2    = Thrust(i+1) - flt_T

!       write(*,*) i,Power(i),Power(i+1)!Power(i)
       if (dT1 * dT2 .le. zero) then
          dT          = Thrust(i+1) - Thrust(i)             ! dx
          dP          = Power(i+1)  - Power(i)
          dth         = th0(i+1)    - th0(i) 

          Rotor % Power  = Power(i) +  dP/dT*(flt_T - Thrust(i)) 
          Rotor % th0    =   th0(i) + dth/dT*(flt_T - Thrust(i)) 

!          write(*,*) 'analyzed data and interpolated',i,th0(i)
!          write(*,*) dP, dT, dth

          exit
       end if
     end do  
  end if

!  write(*,*) Rotor % Power, rotor % th0
!===========================================================================
! end of operations
!===========================================================================

return
end subroutine