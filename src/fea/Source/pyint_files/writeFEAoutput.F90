! ###################################################################
!
! outputs
!
! ###################################################################

subroutine writeFEAoutputs

   use currentvalues
   use layout
   use material

   implicit none

   integer n,j,jst,jend,offset,k
   !integer half_ndof_beam
   !
   ! dof.dat
   !
   open(unit=1,file='./output/dof.dat')
   write(1,*) '# number of loading conditions'
   write(1,*) nloading

   do n = 1,nloading
      write(1,*) '# NodeID DOF (u, v, w, thx, thy, thz)'
      do j = 1,nnodes
         jst  = ndof_beam*(j-1)+1
         jend = ndof_beam*(j-1)+6
         write(1,'(I4,6E16.6)') j,xdof(n,jst:jend)

      enddo
   enddo
   close(1)
   !
   ! stress.dat
   !
   open(unit=1,file='./output/stress.dat')
   write(1,*) '# ElemID Stress-strain (bucklingFOS, VM-Stress (MPa), yieldFOS)'
   do j = 1,nelem
      write(1,'(I4,3E16.6)') j,bucklingFOS(j),stress(j)*1.e-6,stressFOS(j)
   enddo
   close(1)
   !
   ! optimize.dat
   !
   open(unit=1,file='./output/optimize.dat')
   write(1,*) '# Iteration count, min FOS, min weight (kg), max deflection (m)'
   write(1,*) optcount
   do j = 1,optcount
      write(1,'(I4,3E16.6)') j,minOptStressFOS(j), minOptWeight(j), maxOptDeflection(j)
   enddo
   close(1)
   !
   ! modeshapes.dat
   !
   half_ndof_beam = 0.5*ndof_beam
   open(unit=1,file='./output/modeshapes.dat')
   write(1,*) '# something .. '
   write(1,*) nmodes
   do j = 1,nnodes*half_ndof_beam

      if (j==1) then
         do k = 1,nmodes
            write(1,'(E16.6,$)') modeFreqs(k)
         enddo
         write(1,*)
      endif

      do k = 1,nmodes
         write(1,'(E16.6,$)') modeShapes(k,j)
      enddo
      write(1,*)
   enddo
   close(1)
   !
   ! dimension_member.dat
   !
   open(unit=1,file='./output/dimension_member.dat')
   write(1,*) '# nmember, outer dimension [m]'
   write(1,*) nmember
   do j = 1,nmember
      write(1,'(E16.6)') matprop(10,elem2member(1,j))
   enddo
   close(1)
   


end subroutine writeFEAoutputs
