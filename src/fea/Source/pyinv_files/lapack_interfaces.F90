!======================================================================
!These subroutines provide interfaces to the netlib canned procedures
!that perform standard linear-algebra operations on matrices 
!======================================================================

! 1 - matrix inversion                                      DGECO/DGEDI
! 2 - simultaneous linear equation solver                   DGESV

      subroutine INVERSE(A,LDA,N,RCOND,DETERM)
      implicit none

!======================================================================
! Inputs
!======================================================================

      integer        , intent(in)   :: N, LDA

!======================================================================
! Input/Output
!======================================================================
      
      real(kind=8) , intent(inout):: A(LDA, N)

!======================================================================
! Outputs
!======================================================================
      
      real(kind=8) , intent(inout):: rcond, determ
      
!======================================================================
! Local variables
!======================================================================

      INTEGER                       :: JOB, IPVT(N)
      real(kind=8)                :: DET(2)
      real(kind=8)                :: WORK2(N) ,WORK(N)     

!======================================================================
! Begin executable code
!======================================================================
      
!set job to 11 to get both determinant and inverse
      JOB=11      

!double precision version
!#ifdef usedouble

      CALL DGECO(A,LDA,N,IPVT,RCOND,WORK2)
      CALL DGEDI(A,LDA,N,IPVT,DET,WORK,JOB)

!single precision version
!      call sgeco(A,LDA,N,IPVT,RCOND,WORK2)
!      CALL SGEDI(A,LDA,N,IPVT,DET,WORK,JOB)
!#endif

      DETERM=DET(1) * 10.0**DET(2)
      
      return
      end subroutine
      
