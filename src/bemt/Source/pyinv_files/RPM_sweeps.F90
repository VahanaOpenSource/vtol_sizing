!===========================================================================
! fortran routine to perform RPM sweep and store power, thrust arrays 
! for several flight conditions
!===========================================================================

subroutine rpm_sweeps(inputs, flt, design, rotor, dvars, save_power)

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

!===========================================================================
! build rotor chord and twist distribution from design variables 
!===========================================================================

    call interpret_geometry(design, Rotor)

!===========================================================================
! set reference rotor RPM and collective; first will change, second is fixed
!===========================================================================

    Rotor % omega   = inputs % omegaH

!===========================================================================
! loop over flight conditions, set root collective based on 
! hover setting flt(1) % coll or variable collective, flt(i) % coll 
!===========================================================================

    do i = 1, inputs % N_flight

        Rotor % th0     = flt(1) % coll 
        call RPM_sweep(flt(i), Rotor, 1)

        Rotor % th0     = flt(i) % coll
        call RPM_sweep(flt(i), Rotor, 2)
    end do 

!===========================================================================
! end of operations
!===========================================================================

return
end subroutine rpm_sweeps
