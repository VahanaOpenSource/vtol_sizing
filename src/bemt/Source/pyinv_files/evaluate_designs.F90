!===========================================================================
!> This subroutine calculates performance for all blade designs
!===========================================================================

    subroutine evalute_designs(inputs, flt, all_cases, Rotor, design_vars, nd, best_design)

    use bemt_data_structures
    use bemt_types
    implicit none

!===========================================================================
! input/output
!===========================================================================

    integer         , intent(in)        :: nd 
    type(dvars_def) , intent(in)        :: design_vars
    type(inter_def) , intent(in)        :: inputs 
    type(flt_def)   , intent(inout)     :: flt(nfc)
    type(Rotor_def) , intent(inout)     :: Rotor 
    type(case_def)  , intent(inout)     :: best_design, all_cases(nd) 

!===========================================================================
! local variables
!===========================================================================

    integer             :: j, omp_get_thread_num, omp_get_max_threads


!===========================================================================
! performance calculations
!===========================================================================

#ifdef openmp 
    call omp_set_num_threads(4)
    write(*,*) 'engaging openmp parallelization for bemt cases'
!$OMP PARALLEL DO DEFAULT(SHARED) FIRSTPRIVATE(Rotor) PRIVATE(j)
#endif

    do j = 1,nd
!        write(*,*) 'hello from thread # ',omp_get_thread_num()
        call evaluate_design(inputs, flt, all_cases(j), Rotor, design_vars, .false.)
    end do 

#ifdef openmp 
!$OMP END PARALLEL DO
#endif

!===========================================================================
! find the design that uses minimum energy
!===========================================================================

    call find_top_design(all_cases(1:nd), nd, j)

!===========================================================================
! if there is NO valid design, trigger "valid" flag
!===========================================================================

    if (j == -1) then 
        best_design % valid = .false.
    else

!===========================================================================
! evaluate top design, and save the power 
!===========================================================================

        best_design     = all_cases(j)
        if (best_design % valid) then
            call evaluate_design(inputs, flt, best_design, Rotor, design_vars, .true.)
        end if 
    end if 

!    write(*,*) 'is design valid?',best_design % valid

!===========================================================================
! end of operations
!===========================================================================

    return 
    end subroutine evalute_designs