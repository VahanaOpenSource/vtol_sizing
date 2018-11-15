!===========================================================================
! fortran routine to convert solidity, bilinear twist and linear taper
! to chord, twist distribution along span 
!===========================================================================

subroutine interpret_geometry(design, Rotor)

    use bemt_data_structures
    use bemt_interface
    implicit none

!===========================================================================
! input/output
!===========================================================================

    type(case_def)   , intent(in)    :: design
    type(Rotor_def)  , intent(inout) :: Rotor 

!===========================================================================
! local variables
!===========================================================================

    integer             :: i 
    real(kind=rdp)      :: radius, cbar, ctip, croot, th_offset

!===========================================================================
! begin executable code
! note: we calculate cbar from solidity assuming no root cutout
!===========================================================================

    radius      = Rotor % radius 
    cbar        = Rotor % solidity * pi * radius / real(rotor % Nb)

!===========================================================================
! for chord distribution: array
!===========================================================================

    ctip  = two*cbar/(one + design % taper)
    croot = design % taper*ctip

    do i = 1,Nseg
        Rotor % chord(i)  = croot + (ctip - croot)*Rotor % r(i)
    end do 

!===========================================================================
! for twist distribution: array
!===========================================================================

    do i = 1, Nseg 
        if (Rotor % r(i) .le. design % x) then      ! inboard of bilinear junction 

            Rotor % twist(i) = - Rotor % r(i) * design % thx/design % x         
        else
            Rotor % twist(i) = -design % thx -                                  &
                                (design % thtip - design % thx)*(Rotor % r(i) - design % x)/(one - design % x)  
        end if 
    end do 

!===========================================================================
! find twist at 75% span 
!===========================================================================

    if(0.75d0 .le. design % x) then 
        Rotor % th_75   = - 0.75d0 * design % thx / design % x
    else 
        Rotor % th_75   = -design % thx -                               &
                            (design % thtip - design % thx)*(0.75d0 - design % x)/(one - design % x)  
    end if 

!===========================================================================
! offset twist distribution so that 75% span has zero geometric twist
!===========================================================================

    th_offset           = Rotor % th_75 
    Rotor % twist       = Rotor % twist - th_offset
    Rotor % th_75       = Rotor % th_75 - th_offset 

!===========================================================================
! debug prints
!===========================================================================

    ! do i = 1, Nseg
    !     write(*,*) i, Rotor % r(i), Rotor % chord(i), Rotor % twist(i)
    ! end do 
    ! stop 'ok?'
!===========================================================================
! end of operations
!===========================================================================

return
end subroutine