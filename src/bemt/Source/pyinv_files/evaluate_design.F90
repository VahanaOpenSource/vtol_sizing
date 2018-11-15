!===========================================================================
! fortran routine to perform collective sweep and obtain power prediction
! for a given thrust & RPM using BEMT
!===========================================================================

subroutine evaluate_design(inputs, flt, design, rotor, dvars, save_power)

    use bemt_data_structures
    use bemt_types
    implicit none

!===========================================================================
! input/output
!===========================================================================

    logical         , intent(in)        :: save_power


    type(dvars_def) , intent(in)        :: dvars
    type(inter_def) , intent(in)        :: inputs 
    type(flt_def)   , intent(inout)     :: flt(nfc)
    type(Rotor_def) , intent(inout)     :: Rotor 
    type(case_def)  , intent(inout)     :: design 

!===========================================================================
! local variables
!===========================================================================

    logical         :: valid_soln
    integer         :: i, k, nvalid, idmin 
    real(kind=rdp)  :: omega(mxOmega), Power(mxOmega), Pmin, th0(mxOmega)

!===========================================================================
! build rotor chord and twist distribution from design variables 
!===========================================================================

    call interpret_geometry(design, Rotor)

    design % valid      = .true.
    design % energy     = zero

!===========================================================================
! if outboard section is more nose-up than inboard section, exit 
!===========================================================================
    
    if (design % thx .gt. design % thtip) then 
        design % valid  = .false.
        return 
    end if 

!===========================================================================
! loop over flight conditions
!===========================================================================

    do i = 1, inputs % N_flight

!===========================================================================
! if its hover, use hover rpm (input) and obtain power for target thrust
! hover is when climb velocity is very small in magnitude
!===========================================================================

        if (abs(flt(i) % Vc) .le. 1.d-3) then 

            Rotor % omega       = inputs % omegaH
            call trim_thrust(flt(i), Rotor, valid_soln)
!            write(*,*) 'valid soln?',valid_soln
!===========================================================================
! for invalid solutions, don't run other flight conditions 
! if the design can't hover, don't consider it!
!===========================================================================

            if (.not. valid_soln) then  
                design % valid  = .false.
!                write(*,*) 'HAAAA: PAM?' 
                exit 
            else 
                design % energy = design % energy + Rotor % power * flt(i) % time

                if (save_power) then 
                    flt(i) % power      = Rotor % Power 
                    flt(i) % omega      = Rotor % omega 
                    flt(i) % coll       = Rotor % th0
                    flt(i) % efficiency = Rotor % eta
                end if 
!remember minimum FM
                if (Rotor % eta .le. design % FMmin) then 
                    design % FMmin  = Rotor % eta
                end if  

            end if 

!===========================================================================
! if its climb/cruise in axial flight, also check various RPMs
!===========================================================================
    
        else
            nvalid = 0 
            do k = 1, dvars % nOmega
                Rotor % omega   = dvars % omega(k)
                call trim_thrust(flt(i), Rotor, valid_soln)
                if(valid_soln) then
                    nvalid           = nvalid + 1 
                    Power(nvalid)    = Rotor % power
                    omega(nvalid)    = Rotor % omega
                    th0(nvalid)      = Rotor % th0      ! collective, deg 

!===========================================================================
! if there is a case working at lower RPM and it has lesser power, use it
!===========================================================================

                    if(nvalid .ge. 2) then 
                        if (Power(nvalid) .lt. Power(nvalid-1)) then 
                            exit
                        end if 
                    end if 
                end if
            end do 

            if (nvalid == 0) then
                design % valid = .false. 
                exit 
            else 
                Pmin = Power(1); idmin   = 1; 
                do k = 2, nvalid
                    if (Pmin .ge. Power(k)) then 
                        Pmin             = Power(k); idmin  = k
                    end if 
                end do                     

!===========================================================================
! rerun optimal RPM case if save mode is enabled
!===========================================================================

                if (save_power) then 
                    Rotor % omega       = omega(idmin)
                    call trim_thrust(flt(i), Rotor, valid_soln)
                    flt(i) % power      = Rotor % power 
                    flt(i) % omega      = Rotor % omega
                    flt(i) % coll       = Rotor % th0 
                    flt(i) % efficiency = Rotor % eta
                    design % omega      = Rotor % omega 
                end if 
                design % energy     = design % energy + Pmin * flt(i) % time

!remember minimum FM
                if ( Rotor % eta .le. design % etamin) then 
                    design % etamin  = Rotor % eta
                end if  

            end if              ! end if valid design
        end if                  ! end if in climb

    end do                      ! loop over flight conditions

!===========================================================================
! add filter for minimum FM, eta_p
!===========================================================================
    
    if (design % FMmin .le. 0.75d0 .or. design % etamin .le. 0.8d0) then 
        design % valid = .false. 
    end if 

!    write(*,*) 'design energy is ',design % energy, design % valid
return
end subroutine

!===========================================================================
! fortran routine to find the best design in a set based on energy usage
!===========================================================================

subroutine find_top_design(designs, nd, id)

use bemt_types 
!===========================================================================
! input/output
!===========================================================================

    integer         , intent(in)    :: nd 
    type(case_def)  , intent(inout) :: designs(nd) 

    integer         , intent(out)   :: id 

!===========================================================================
! local variables
!===========================================================================

    logical             :: found_base 
    integer             :: i
    real(kind=rdp)      :: e_min

!===========================================================================
! begin executable code 
!===========================================================================

    found_base = .false.
    id         = -1
    e_min      = 1.d15
    
!===========================================================================
! loop over all designs
!===========================================================================

    do i = 1, nd 

!===========================================================================
! if found a valid design and current design is better, set it as best 
!===========================================================================

        if (found_base) then
            if(designs(i) % valid) then 
                if(designs(i) % energy .le. e_min) then 
                    id      = i 
                    e_min   = designs(i) % energy 
                end if     
            end if 

!===========================================================================
! if we havent established a best, set it to first valid design in list
!===========================================================================

        else 
            if(designs(i) % valid) then 
                found_base  = .true. 
                id          = i
                e_min       = designs(i) % energy 
            end if 
        end if 
    end do 

!===========================================================================
! end of operations
!===========================================================================

return 
end subroutine