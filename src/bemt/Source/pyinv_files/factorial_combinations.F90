!===========================================================================
!> fortran routine to generate factorial combinations of design variables
!===========================================================================

subroutine factorial_combinations(cases, nd, dvars)

    use bemt_types
    use bemt_data_structures
    implicit none

!===========================================================================
! input/output
!===========================================================================

    integer         , intent(in)      :: nd
    type(dvars_def) , intent(in)      :: dvars 
    type(case_def)  , intent(inout)   :: cases(nd) 

!===========================================================================
! local variables
!===========================================================================

    integer     :: i, j, k, l, m, ic 

!===========================================================================
! begin executable code
!===========================================================================

!===========================================================================
! populate design variables
!===========================================================================

    ic   = 0
    do                  i = 1, ntaper
        do              j = 1, nthx 
            do          k = 1, nthtip
                do      l = 1, nx 
                        
                    ic = ic + 1
                    cases(ic) % taper = dvars % taper(i) 
                    cases(ic) % thx   = dvars % thx(j)   
                    cases(ic) % thtip = dvars % thtip(k) 
                    cases(ic) % x     = dvars % x(l)

                end do 
            end do 
        end do 
    end do 

!===========================================================================
! end of operations
!===========================================================================

return
end subroutine