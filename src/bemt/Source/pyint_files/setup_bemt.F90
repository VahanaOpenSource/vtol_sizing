!===========================================================================
! fortran routine to perform collective sweep and obtain power prediction
! for a given thrust & RPM using BEMT
!===========================================================================

subroutine setup_bemt()

    use bemt_data_structures
    use bemt_interface
    implicit none

!===========================================================================
! local variables
!===========================================================================

    integer                 :: j, nd, omp_get_thread_num, omp_get_max_threads
    type(Rotor_def)         :: Rotor 
    type(dvars_def)         :: design_vars

!all designs 
    type(case_def)       , allocatable  :: all_cases(:) 
    type(case_def)                      :: design 

!===========================================================================
! begin executable code
!===========================================================================

    Rotor % use_tables      = Inputs % use_tables 
    Rotor % Cdo             = Inputs % Cdo 
    Rotor % CL_alpha        = Inputs % CL_alpha 

    if(Rotor % use_tables) call setup_tables()

!===========================================================================
! generate blade twist and taper arrays
!===========================================================================

    call populate_design_space(Inputs, Rotor, design_vars, nd)

!===========================================================================
! check mode of operation: default mode is parameter sweep of blade geometry
! i.e., try all twist/taper combinations and figure out which has least fuel
! however.. another mode of operation exists; we specify the design coming 
! in, and just calculate the performance
!===========================================================================

    if (Inputs % mode_operation == 1) then 
        write(6,*) 
        write(6,*) 'running BEMT analysis in fixed RPM, variable collective mode'
        write(6,*) 'evaluating performance of best design'
        call evaluate_design(inputs, flt, best_design, Rotor, design_vars, .true.)
        return 
    elseif(Inputs % mode_operation == 2) then 
        write(6,*) 
        write(6,*) 'running BEMT analysis in fixed collective, variable-RPM mode'
        write(6,*) 'computing thrust-power curve with fixed-pitch propellers'        
        call RPM_sweeps(inputs, flt, best_design, Rotor, design_vars, .true.)
    else
        write(6,*) 'running BEMT analysis with blade geometry parameter sweep'
    end if 

    if (allocated(all_cases)) deallocate(all_cases); allocate(all_cases(nd))

!===========================================================================
! loop over all geometries and evaluate the design 
! the metric used is "energy" need to power the design through all flight 
! conditions given 
!===========================================================================

!===========================================================================
! generate factorial combinations of designs
!===========================================================================

    call factorial_combinations(all_cases, nd, design_vars)

!===========================================================================
! evaluate all designs and return best design if relevant
!===========================================================================

    call evalute_designs(inputs, flt, all_cases, Rotor, design_vars, nd, best_design)

!===========================================================================
! deallocate memory 
!===========================================================================

    if (allocated(all_cases)) deallocate(all_cases)

!===========================================================================
! end of operations
!===========================================================================

return
end subroutine
