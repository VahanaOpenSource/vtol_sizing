!
! routine to optimize structure
!
subroutine optimizeStructure(ioptimize)

   use layout, only : nmodes,nloading,nmember
   use currentvalues
   use fea_types
   implicit none
   !
   integer, intent(inout) :: ioptimize
   !
   integer        :: n
   integer        :: icount, countmax, iconstraint,iconverge
   real(kind=8)   :: error
   logical        :: failure 

! ===================================================================
! allocate
! ===================================================================

   !if(allocated(minOptWeight))     deallocate(minOptWeight)
   !if(allocated(minOptStressFOS))  deallocate(minOptStressFOS)
   !if(allocated(maxOptDeflection)) deallocate(maxOptDeflection)
   !if(allocated(idefFlag))         deallocate(idefFlag)
   !if(allocated(istressFlag))      deallocate(istressFlag)
   !allocate(minOptWeight(nmaxopt),minOptStressFOS(nmaxopt),maxOptDeflection(nmaxopt))
   !allocate(idefFlag(nmember))
   !allocate(istressFlag(nmember))

! ===================================================================
! if no optimization is required
! assemble the matrices, apply BC
! solver the system and calculate loads
! ===================================================================

   icount      = 1
   failure     = .false.
   call assemble()

!====================================================================
! for modal analysis, we have global M, K assembled, exit
!====================================================================

   if(iodata % get_modes) return

!====================================================================
! if still here, we're doing static analysis/optimization
!====================================================================

   call resetLoads(icount)
   call applybc()
   call solver(failure)
   if(failure) then
      ioptimize = -1
      return 
   end if  

! ===================================================================
! loop over number of loading conditions
! ===================================================================

   do n = 1,nloading
      call calculateLoads(n,icount)
   enddo

   iconverge = 0
   icount    = icount + 1
   countmax  = nmaxopt

! ===================================================================
! if optimization of the structure is required
! ===================================================================

   if(ioptimize == 1) then

      do while (icount < countmax+1 .and. iconverge==0) 

! ===================================================================
! update the cross section characteristics 
! ===================================================================

         call updateStructure(iconstraint)

! ===================================================================
! convergence exit condition
! this implies that the deflection is above the specified limit and
! the bending FOS is within the specified band. However, there is the
! case of 0 loads, where there is a minimum size below which the 
! members will not get sized. Because the FOS is not within the band
! the code will never converge.
! ===================================================================

         if(iconstraint==1) then
            iconverge=1
            exit
         endif

! ===================================================================
! reassemble stiffness and mass matrix, solve for dof's
! ===================================================================

         call assemble

! ===================================================================
! loop through the loading cases
! ===================================================================

         call resetLoads(icount)
         call applybc()
         call solver(failure)

         if(failure) then
            ioptimize = -1
            return 
         end if  

         do n = 1,nloading

            call calculateLoads(n,icount)

         enddo

! ===================================================================
! this condition effectively checks if all the stress and
! deflection conditions have been met and the weight
! between iterations do not change by 1%
!
! failsafe implemented to exit at countmax
! ===================================================================

         if(icount > 1)then
            if(sum(idefFlag)==nmember .and. sum(istressFlag)==nmember)then
               error=abs(minOptWeight(icount) - minOptWeight(icount-1))/minOptWeight(icount-1)
               if(error < 0.05d0 .or. icount == countmax) then
                  iconverge=1
                  exit
               endif
            endif
         endif
         icount = icount + 1
      enddo

      optcount = icount-1
   endif

end subroutine optimizeStructure
