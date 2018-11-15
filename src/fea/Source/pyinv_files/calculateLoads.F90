! ###################################################################
!
! calculate loads
!
! ###################################################################

subroutine calculateLoads(loadingid,icount)
   
   use currentvalues
   use layout
   use material

   implicit none
   !
   integer, intent(in) :: icount,loadingid
   !
   integer :: i,k,off,kk,ii
   integer :: j1,j2,j3
   integer :: matid,node1,node2
   integer :: ist,iend,elemID

   real(kind=8) :: l1,l2,l3,m1,m2,m3,n1,n2,n3
   real(kind=8) :: L,invl,invl2,invl3,sum
   real(kind=8) :: T(12,12),ldof(12),ldoftmp(12)
   real(kind=8) :: area,E,G,J,phi_z,phi_y,Iy,Iz,ky,kz,rho,nu
   real(kind=8) :: V1,V2,V3,M11,M2l,M2r,M3l,M3r,bucklingLoad,maxy,maxz
   real(kind=8) :: sig(6,6),yy,zz,sigv,maxsigv,sigYield
   !
   real(kind=8) :: newRadius, d2, d1
   real(kind=8) :: vectors(3,3)
   !
   real(kind=8), parameter :: pi = 4.d0*atan(1.d0)

! ===================================================================   
! allocatable
! ===================================================================   

   !if(allocated(stress))      deallocate(stress)
   !if(allocated(stressFOS))   deallocate(stressFOS)
   !if(allocated(strain))      deallocate(strain)
   !if(allocated(bucklingFOS)) deallocate(bucklingFOS)
   !if(allocated(weightElem))  deallocate(weightElem)
   !if(allocated(deflection))  deallocate(deflection)
   !allocate(stress(nelem),strain(nelem),bucklingFOS(nelem),stressFOS(nelem),weightElem(nelem))
   !allocate(deflection(nnodes))

! ===================================================================   
! loop through the elements
! ===================================================================   

   totalWeight  = 0.d0

   do i = 1,nelem
      E        = matprop( 1,i)
      nu       = matprop( 2,i)
      rho      = matprop( 3,i)
      area     = matprop( 4,i)
      Iy       = matprop( 5,i)
      Iz       = matprop( 6,i)
      J        = matprop( 7,i)
      ky       = matprop( 8,i)
      kz       = matprop( 9,i)
      maxy     = matprop(10,i)
      maxz     = matprop(11,i)
      sigYield = matprop(12,i)
      !
      G        = E/(2.d0*(1.d0+nu))
      !
      node1    = conn(1,i)
      node2    = conn(2,i)

! ===================================================================   
!find distances of nodes from origin
! ===================================================================   

      d1       = nodepos(1,node1)*nodepos(1,node1) +           &
                 nodepos(2,node1)*nodepos(2,node1) +           &
                 nodepos(3,node1)*nodepos(3,node1)
      d1       = sqrt(d1)

      d1       = d1 + 1.d-6 

      d2       = nodepos(1,node2)*nodepos(1,node2) +           &
                 nodepos(2,node2)*nodepos(2,node2) +           &
                 nodepos(3,node2)*nodepos(3,node2)
      d2       = sqrt(d2)

      d2       = d2 + 1.d-6

! ===================================================================   
! find element length
! ===================================================================   

      L = (nodepos(1,node2)-nodepos(1,node1))*(nodepos(1,node2)-nodepos(1,node1)) &
        + (nodepos(2,node2)-nodepos(2,node1))*(nodepos(2,node2)-nodepos(2,node1)) &
        + (nodepos(3,node2)-nodepos(3,node1))*(nodepos(3,node2)-nodepos(3,node1))
      L = sqrt(L)
      !
      invl   = 1.d0/L
      invl2  = 1.d0/(L*L)
      invl3  = 1.d0/(L*L*L)

! ===================================================================   
! find local coordinates
! ===================================================================   

      call findLocalCoordinates(node1,node2,vectors)
      !
      l1 = vectors(1,1)
      m1 = vectors(1,2)
      n1 = vectors(1,3)
      !
      l2 = vectors(2,1)
      m2 = vectors(2,2)
      n2 = vectors(2,3)
      !
      l3 = vectors(3,1)
      m3 = vectors(3,2)
      n3 = vectors(3,3)
      
! ===================================================================   
! define 12x12 transformation matrix
! ===================================================================   

      T = 0.d0
      do k = 1,4
         off = 3*(k-1)
         T(off+1,off+1) = l1
         T(off+1,off+2) = m1
         T(off+1,off+3) = n1
         T(off+2,off+1) = l2
         T(off+2,off+2) = m2
         T(off+2,off+3) = n2
         T(off+3,off+1) = l3
         T(off+3,off+2) = m3
         T(off+3,off+3) = n3
      enddo
      !
      ist            = ndof_beam*(node1-1) + 1 
      iend           = ndof_beam*(node1-1) + 6
      ldoftmp(1:6)   = xdof(loadingid,ist:iend)
      !
      ist            = ndof_beam*(node2-1) + 1 
      iend           = ndof_beam*(node2-1) + 6
      ldoftmp(7:12)  = xdof(loadingid,ist:iend)
      !
      do k = 1,2*ndof_beam
         sum = 0.d0
         do kk = 1,2*ndof_beam
            sum = sum + T(k,kk) * ldoftmp(kk)
         enddo
         ldof(k) = sum
      enddo

! ===================================================================   
! calculate nodal defletions
! ===================================================================   

      deflection(node1) = sqrt(ldof(1)*ldof(1) + ldof(2)*ldof(2) + ldof(3)*ldof(3))
      deflection(node2) = sqrt(ldof(7)*ldof(7) + ldof(8)*ldof(8) + ldof(9)*ldof(9))

!=========================================================================  
! normalize deflections by distance from center of vehicle (set at origin!!)
! this way, we get a "normalized" deflection for all scales of construction 
! rather than beging tied to a dimensional displacement value
!=========================================================================

      if( abs(d1) .ge. 1.d-5) deflection(node1) = deflection(node1)/d1 
      if( abs(d2) .ge. 1.d-5) deflection(node2) = deflection(node2)/d2 

! ===================================================================   
! calculate nodal forces
! ===================================================================   

      V1  = E*area*(ldof(7)-ldof(1))/L
      V2  = E*Iy*( 12.d0*invl3*ldof(3) - 6.d0*invl2*ldof(5) - 12.d0*invl3*ldof(9) - 6.d0*invl2*ldof(11) )
      V3  = E*Iz*( 12.d0*invl3*ldof(2) - 6.d0*invl2*ldof(6) - 12.d0*invl3*ldof(8) - 6.d0*invl2*ldof(12) )
      !
      M11  = G*J*(ldof(10)-ldof(4))/L
      M2l = E*Iy*( -6.d0*invl2*ldof(3) + 4.d0*invl*ldof(5) + 6.d0*invl2*ldof(9) + 2.d0*invl*ldof(11) )
      M2r = E*Iy*(  6.d0*invl2*ldof(3) - 2.d0*invl*ldof(5) - 6.d0*invl2*ldof(9) - 4.d0*invl*ldof(11) )
      M3l = E*Iy*( -6.d0*invl2*ldof(2) - 4.d0*invl*ldof(6) + 6.d0*invl2*ldof(8) - 2.d0*invl*ldof(12) )
      M3r = E*Iy*(  6.d0*invl2*ldof(2) + 2.d0*invl*ldof(6) - 6.d0*invl2*ldof(8) + 4.d0*invl*ldof(12) )


      !print*,'elem: ',i,'forces: ',M11,M2l,M2r,M3l,M3r

! ===================================================================   
! find buckling load
! ===================================================================   

      bucklingLoad = pi*pi*E*Iy/(fBuckling*fBuckling)*invl2
      !
      !if (V1 > bucklingLoad) then
      !   print*,'calculateLoads.F90: beam element failed because of buckling'
      !   maxBuckling = -1.d0
      !else
!modified by AS: desingularize
         maxBuckling = max(maxBuckling,bucklingLoad/(abs(V1)+1.d-2))
      !endif
      !
      bucklingFOS(i) = maxBuckling

!      if(V1 .le. -1.d-2) write(*,*) 'hi',V1,bucklingLoad,bucklingFOS(i)
! ===================================================================   
! von-Mises stress (check on 4 points on the cross section)
! ===================================================================   

      sig=0.d0
      maxsigv=0.d0
      do ii = 1,2 ! for nodes at either end
         do k = 1,4  ! 4 points at each cross-section
            select case (k)
               case (1)
                  yy =  maxy; zz =  0.d0
               case (2)
                  yy = 0.d0;  zz =  maxz
               case (3)
                  yy = -maxy; zz =  0.d0
               case (4)
                  yy = 0.d0;  zz = -maxz
            end select
            !
            if(ii==1) sig(1,1) = V1/area - M3l*yy/Iz - M2l*zz/Iy
            if(ii==2) sig(1,1) = V1/area - M3r*yy/Iz - M2r*zz/Iy


            ! uni-axial stress assumption
            sig(2,2) = 0.d0
            sig(3,3) = 0.d0
            sig(2,3) = 0.d0

            sigv = (sig(1,1)-sig(2,2))*(sig(1,1)-sig(2,2)) &
                 + (sig(2,2)-sig(3,3))*(sig(2,2)-sig(3,3)) &
                 + (sig(3,3)-sig(1,1))*(sig(3,3)-sig(1,1)) &
                 + 6.d0*(sig(2,3)*sig(2,3) + sig(1,2)*sig(1,2) + sig(1,3)*sig(1,3))

            sigv = sqrt(sigv)


            maxsigv = max(maxsigv,sigv)

         enddo
      enddo


      ! store only the von-mises stress
      stress(i) = maxsigv                    ! magnitude of stress
        
      ! stress FOS of safety in each element
      stressFOS(i) = sigYield/maxsigv        ! always positive

!=========================================================================
!modification by AS
! for compressive loads (V1 < 0), consider buckling also
! choose lower of  buckling safety factor or Von.Mises stress based SF
!=========================================================================

      if(V1 .lt. -1.d-2) then        !compressive load of 0.01N  
         if( stressFOS(i) .ge. bucklingFOS(i) ) then 
            stressFOS(i) = bucklingFOS(i)
         end if 
      end if 

!=========================================================================
! weight of each element
!=========================================================================

      weightElem(i) = rho * area * L
      totalWeight   = totalWeight + weightElem(i)

   enddo

! ===================================================================   
! store the max stresses within each member
! ===================================================================   
   !
   do i = 1,nmember
      do k = 1,nelem2member(i)
         elemID       = elem2member(k,i)
         memStress(i) = max(memStress(i),stress(elemID))
         memFOS(i)    = min(memFOS(i),stressFOS(elemID))
         !
         ! check both nodes of a given element
         !
         memDeflection(i) = max(memDeflection(i), deflection(conn(1,elemID)))
         memDeflection(i) = max(memDeflection(i), deflection(conn(2,elemID)))

      enddo
      
      minFOS        = min(memFOS(i),minFOS)
      maxDeflection = max(memDeflection(i),maxDeflection)


   enddo

! ===================================================================   
! save values to opt array
! ===================================================================   

   minOptWeight(icount) = totalWeight
   do i = 1,nmember
      minOptStressFOS(icount)  = min(minOptStressFOS(icount),memFOS(i))
      maxOptDeflection(icount) = max(maxOptDeflection(icount),memDeflection(i))
   enddo

   if(.false.)then
      write(*,'(A)') ' --------------------------------------------------------------------'
      do i = 1,nmember
         !write(*,'(A,I3,A,E16.6,A,E16.6,A,E16.6)') '  Member ' , i , ' || Minimum stress FOS: ',memFOS(i), &
         write(*,'(A,I3,A,E16.6,A,E16.6,A,E16.6)') '  Member ' , i , ' || Minimum stress FOS: ',memFOS(i), &
            ' || Deflection: ',memDeflection(i),' || Buckling FOS: ',bucklingFOS(i)
      enddo
      write(*,'(A,F12.3,A)') ' Total weight:       ',totalWeight,' kg'
      write(*,'(A)') ' --------------------------------------------------------------------'
      write(*,*) ''
   endif
   



end subroutine calculateLoads

! ###################################################################
!
! reset loads
!
! ###################################################################

subroutine resetLoads(icount)

   use layout
   use currentvalues

   implicit none

   integer, intent(in) :: icount
   !
   ! moved from solver
   !
   !if(allocated(xdof)) deallocate(xdof)
   !allocate(xdof(nloading,ndof))
   !!
   !! moved from calculate loads
   !!
   !if(allocated(memStress))     deallocate(memStress)
   !if(allocated(memFOS))        deallocate(memFOS)
   !if(allocated(memDeflection)) deallocate(memDeflection)

   !allocate(memStress(nmember),memFOS(nmember))
   !allocate(memDeflection(nmember))
   !
   ! calculate the maximum stress within each member
   ! and the minimum FOS within each member
   !
   memStress = 0.d0
   memFOS    = 10000.d0
   memDeflection = 0.d0
   !
   minFOS        = 100000.d0
   maxDeflection = 0.d0
   !
   minOptStressFOS(icount) = 1000.d0
   maxOptDeflection(icount) = 0.d0
   !
   maxBuckling  = 0.d0



end subroutine resetLoads


