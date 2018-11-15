! ###################################################################
!
! input-output routines
!
! ###################################################################

subroutine readInputs

   use layout
   use material

   implicit none

   integer      :: j,k,n
   integer      :: nodeid,dofid,maxnforce
   real(kind=8) :: dummy

   !
   ! read input file
   !
   open(unit=1,file="../inputs/config.inp")
 
   read(1,*)
   read(1,*) nnodes,nelem
   !
   if(allocated(nodepos)) deallocate(nodepos)
   allocate(nodepos(3,nnodes))
   !
   ndof = ndof_beam * nnodes
   !
   ! node position
   read(1,*)
   do j = 1,nnodes
      read(1,*) nodepos(1,j),nodepos(2,j),nodepos(3,j)
   enddo
   !
   if(allocated(conn))     deallocate(conn)
   if(allocated(memberID)) deallocate(memberID)
   allocate(conn(2,nelem),memberID(nelem))
   !
   ! element connectivity
   !
   read(1,*)
   do j = 1,nelem
      read(1,*) conn(1,j),conn(2,j),memberID(j)
   enddo
   !
   ! read in number of loading conditions and maxnforce
   !
   read(1,*) 
   read(1,*) nloading,maxnforce
   if(allocated(extforce_dofid)) deallocate(extforce_dofid)
   if(allocated(extforce))       deallocate(extforce)
   if(allocated(nforce))         deallocate(nforce) 
   allocate(nforce(nloading))
   allocate(extforce(nloading,maxnforce),extforce_dofid(nloading,maxnforce))
   !
   ! loop through the number of loadings
   !
   do n = 1,nloading
      !
      ! force/moment boundary condition
      !
      read(1,*)
      read (1,*) nforce(n)
      !
      extforce(n,:)=0.d0
      !
      do j = 1,nforce(n)
         read(1,*) nodeid, dofid, dummy

         extforce(n,j)       = extforce(n,j) + dummy

         extforce_dofid(n,j) = ndof_beam * (nodeid-1) + dofid
      enddo
   enddo
   !
   ! read displacement boundary condition
   !
   read(1,*)
   read(1,*) ndisp
   if(allocated(extdisp_dofid)) deallocate(extdisp_dofid)
   if(allocated(extdisp)) deallocate(extdisp)
   allocate(extdisp(ndisp),extdisp_dofid(ndisp))
   !
   do j = 1,ndisp
      read(1,*) nodeid,dofid,extdisp(j)
      extdisp_dofid(j) = ndof_beam * (nodeid-1) + dofid
   enddo
   !
   ! read element properties
   !
   if(allocated(matprop)) deallocate(matprop)
   allocate(matprop(nmatprop,nelem))
   !
   read(1,*)
   do j = 1,nmatprop
      read(1,*) dummy
      matprop(j,:) = dummy
   enddo
   !
   ! read cross section type
   !
   read(1,*)
   read(1,*) cstype
   !
   close(1)

end subroutine readInputs

! ###################################################################
!
! outputs
!
! ###################################################################

subroutine writeOutputs

   use currentvalues
   use layout

   implicit none

   integer n,j,jst,jend,offset,k
   integer half_ndof_beam
   !
   ! dof.dat
   !
   open(unit=1,file='../outputs/dof.dat')
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
   open(unit=1,file='../outputs/stress.dat')
   write(1,*) '# ElemID Stress-strain (bucklingFOS, VM-Stress (MPa), yieldFOS)'
   do j = 1,nelem
      write(1,'(I4,3E16.6)') j,bucklingFOS(j),stress(j)*1.e-6,stressFOS(j)
   enddo
   close(1)
   !
   ! optimize.dat
   !
   open(unit=1,file='../outputs/optimize.dat')
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
   open(unit=1,file='../outputs/modeshapes.dat')
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
   


end subroutine writeOutputs
