!===========================================================================
! assemble global stiffness matrix 
!===========================================================================

subroutine assemble

   use layout
   use currentvalues
   use fea_types
   implicit none
   
   integer      :: j,k,jr,kr,g1,g2, ie
   integer      :: node1,node2,off1,off2
   integer      :: nn, nred, omp_get_thread_num
   real(kind=8), dimension(2*ndof_beam,2*ndof_beam,nelem) :: K_all, M_all

!====================================================================================
! begin executable code
!====================================================================================

!====================================================================================
! find element M,K matrices
!====================================================================================

!   do j = 1, nelem 
!      call beam3d(j,K_all(:,:,j),M_all(:,:,j))
!   end do       

!====================================================================================
! old style: assemble global then apply BC
!====================================================================================

   K_red       = 0.d0 
   nred        = ndof_red
   nn          = 2*ndof_beam              ! size of element matrices

!====================================================================================
! version for only reduced stiffness
!====================================================================================

#ifdef openmp 
   call omp_set_num_threads(4)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ie,node1,node2,off1,off2,k,g2,kr,j,g1,jr) &
!$OMP REDUCTION(+:K_red) 
#endif
   do ie = 1,nelem                  ! loop over elements
!      write(*,*) 'fea: hello from thread # ',omp_get_thread_num()
      node1 = conn(1,ie)            ! node 1 index 
      node2 = conn(2,ie)            ! node 2 index
      !
      off1  = 6*(node1-1)           ! offset index for node 1
      off2  = 6*(node2-1)           ! offset index for node 2

!====================================================================================
! find beam matrices
!====================================================================================

      call beam3d(ie,K_all(:,:,ie),M_all(:,:,ie))

!====================================================================================
! assemble global matrix from local matrix
!====================================================================================

      do k = 1,nn

!====================================================================================
! perform matrix reduction based on BC here
!====================================================================================

         g2 = off1 + k
         if (k > ndof_beam) g2 = off2 + (k-ndof_beam)

         kr       = bc_flag(g2)           ! location in reduced matrix 
         if(kr .ne. 0) then

            do j = 1,nn

               g1 = off1 + j
               if (j > ndof_beam) g1 = off2 + (j-ndof_beam)
               jr  =  bc_flag(g1)      ! storage index in reduced matrix

               if(jr .ne. 0) then
                  K_red(jr,kr)    = K_red(jr,kr) + K_all(j,k,ie)
               end if 
            end do
         end if 
      end do
   end do
#ifdef openmp 
!$OMP END PARALLEL DO 
#endif
!====================================================================================
! version for global M, K (mode extraction)
!====================================================================================

   if(iodata % get_modes) then 
      K_global    = 0.d0 
      M_global    = 0.d0
!      M_red       = 0.d0

      do ie = 1,nelem                  ! loop over elements
         node1 = conn(1,ie)            ! node 1 index 
         node2 = conn(2,ie)            ! node 2 index
         !
         off1  = 6*(node1-1)           ! offset index for node 1
         off2  = 6*(node2-1)           ! offset index for node 2

!====================================================================================
! assemble global matrix from local matrix
!====================================================================================

         do k = 1,nn

!====================================================================================
! perform matrix reduction based on BC here
!====================================================================================

            g2 = off1 + k
            if (k > ndof_beam) g2 = off2 + (k-ndof_beam)

            do j = 1,nn

               g1 = off1 + j
               if (j > ndof_beam) g1 = off2 + (j-ndof_beam)

               K_global(g1,g2) = K_global(g1,g2) + K_all(j,k,ie)
               M_global(g1,g2) = M_global(g1,g2) + M_all(j,k,ie)

            end do
         end do         ! 
      end do
   end if 
end subroutine assemble
