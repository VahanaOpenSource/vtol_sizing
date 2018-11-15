! ###################################################################
!
! code to extract the natural frequencies and modes
!
! ###################################################################

subroutine extractNaturalResponse

   use layout 
   use currentValues


   implicit none
   
   integer :: j,k,kk
   integer :: key(ndof), ierr!, half_ndof_beam

   real(kind=8) :: alfr(ndof),alfi(ndof),beta(ndof)
   real(kind=8) :: srtf(ndof),evec(ndof,ndof)

   if(allocated(modes)) deallocate(modes)
   if(allocated(freqs)) deallocate(freqs)
   allocate(modes(ndof,ndof),freqs(ndof))

! ===================================================================
! compute eigenvalues/eigenvectors
! ===================================================================

   call rgg(ndof,ndof,K_global,M_global,alfr,alfi,beta,1,evec,ierr)
 
! ===================================================================
! compute frequencies 
! ===================================================================

   do j = 1,ndof
      freqs(j) = sqrt(alfr(j)/beta(j))
      key(j)   = j
   enddo

! ===================================================================
! sort by order of increasing frequencies
! ===================================================================

   call DLASRT2('i',ndof,freqs,key,ierr)
   
! ===================================================================
! rerrange eigenvector order (incresing freq)
! ===================================================================

   do j = 1,ndof
      modes(:,j) = evec(:,key(j))
   enddo

! ===================================================================
! compute the mode shapes
! ===================================================================

   !half_ndof_beam = 0.5*ndof_beam
   !if(allocated(modeShapes)) deallocate(modeShapes)
   !if(allocated(modeFreqs))  deallocate(modeFreqs)
   !allocate(modeShapes(nmodes,half_ndof_beam*nnodes),modeFreqs(nmodes))

   kk = 0
   do j = 1,ndof
      
      if( freqs(j) < 1.e-4) cycle;

      ! mode counter
      kk=kk+1

      ! mode frequencies
      modeFreqs(kk) = freqs(j)

      ! mode shapes 2nd dimensions (n1-x,n1-y,n1-z | n2-x,n2-y,n2-z | ... )
      do k = 1,nnodes
         modeShapes(kk, (k-1)*half_ndof_beam + 1) = modes((k-1)*ndof_beam + 1, j)
         modeShapes(kk, (k-1)*half_ndof_beam + 2) = modes((k-1)*ndof_beam + 2, j)
         modeShapes(kk, (k-1)*half_ndof_beam + 3) = modes((k-1)*ndof_beam + 3, j)
      enddo

      if (kk == nmodes) exit;
   enddo

end subroutine extractNaturalResponse

