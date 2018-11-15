!================================================================================
!> This subroutine is a specialized linear equation solver for symmetric matrices
!>
!> The "n" equations we are solving are of the form 
!>             
!>                A x = b
!> 
!> and multiple right hand sides (m). The size of the coefficient matrix is 
!>
!> C(m+n,n)    , where 
!> n           is the number of linear equations and
!> m           is the number of RHS vectors 
!> x(n,m)      is the solution vector ("m columns of size n")
!> 
!> On input, the matrix is partitioned in the following form 
!>           _         _
!>          |           |
!>          |    A      |
!> C     =  |           |
!>          |    b^T    |
!>          |_         _|
!>
!> i.e. the first "n" rows of the input matrix contain the coefficient matrix 
!> and each subsequent row of "C" corresponds to a column of the right hand side 
!> b_i, expressed in row vector form
!>
!> on output, X contains the solution in multiple columns
!> on output, C contains the row-reduced matrix
!================================================================================

   subroutine leqsolver_spl(C, n, m, x)

   implicit none

!================================================================================
! input/output
!================================================================================

   integer     , intent(in)      :: m, n 
   real(kind=8), intent(inout)   :: C(n+m,n)
   real(kind=8), intent(out)     :: x(n,m)

!================================================================================
! local variables
!================================================================================

   integer        :: i, j, k 
   real(kind=8)   :: f1, diag(n), factor

!================================================================================
! begin executable code
!================================================================================

#ifdef openmp
!   call omp_set_num_threads(5)
#endif 
   do k = 1,n              ! all reduction rows 

!================================================================================
! remember inverse of diagonal entry
!================================================================================

      f1             = 1.d00/C(k,k)

!================================================================================
! loop over all other equations and reduce entries
!================================================================================

      do i = 1,n
         if(i .ne. k) then
            factor = C(k,i)*f1
!!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j)
            do j = k,n+m
               C(j,i)  = C(j,i) - factor*C(j,k)
            end do 
!!$OMP END PARALLEL DO 
         end if 
      end do 
      diag(k)  = f1
   end do 

   do i = 1, m
      do j = 1,n
         x(j,i) = C(n+i,j)*diag(j)
      end do 
   end do

!================================================================================
! end of operations
!================================================================================
   
   return
   end subroutine