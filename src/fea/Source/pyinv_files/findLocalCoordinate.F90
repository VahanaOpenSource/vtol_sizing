! ###################################################################
!
! findLocalCoordinates
!
! ###################################################################

subroutine findLocalCoordinates(node1,node2,vectors)

   use layout

   implicit none

   integer, intent(in)  :: node1,node2
   real(kind=8)         :: vectors(3,3)
   !
   integer              :: j,k,jj,kk
   integer              :: nzeros,zeroid(3),izero,ii(2)
   real(kind=8)         :: L,invl
   real(kind=8)         :: a1,b1,c1,p1(3),p2(3),v1(3),v2(3),v3(3),d, dp(3)
   !
   p1(1) = nodepos(1,node1)
   p1(2) = nodepos(2,node1)
   p1(3) = nodepos(3,node1)
   !
   p2(1) = nodepos(1,node2)
   p2(2) = nodepos(2,node2)
   p2(3) = nodepos(3,node2)
   
   dp    = p2 - p1
   !
!   L = (p2(1)-p1(1))*(p2(1)-p1(1)) &
!     + (p2(2)-p1(2))*(p2(2)-p1(2)) &
!     + (p2(3)-p1(3))*(p2(3)-p1(3))
   L  = dp(1)*dp(1) + dp(2)*dp(2) + dp(3)*dp(3)
   L = sqrt(L)
   !
   invl   = 1.d0/L
   !
   ! v1 vector
   !
!   v1(1) = (p2(1)-p1(1))*invl
!   v1(2) = (p2(2)-p1(2))*invl
!   v1(3) = (p2(3)-p1(3))*invl

   v1    = dp*invl
   !
   ! value of d in equation of plane (ax+by+cz=d)
   !
   d = v1(1)*p1(1) + v1(2)*p1(2) + v1(3)*p1(3)
   !
   zeroid = 0
   if(abs(v1(1)) < 1.d-12) zeroid(1)=1
   if(abs(v1(2)) < 1.d-12) zeroid(2)=1
   if(abs(v1(3)) < 1.d-12) zeroid(3)=1
   !
   nzeros = zeroid(1) + zeroid(2) + zeroid(3)
   !
   !print*,'nzeros',nzeros
   !print*,'zeroid',zeroid
   !
   ! case where plane is along 2 axes (2 zeros, e.g. z = 0)
   !
   v2=-1.d0
   if (nzeros == 2) then
   
      kk = 1
      do j = 1,3
         if (zeroid(j) == 0) then
            izero = j
         else
            ii(kk) = j ! ii contains the indices of the two non-zero values
            kk = kk + 1
         endif
      enddo
      !
      ! find vector v2
      !
      v2(izero) = d/v1(izero)
      v3(izero) = d/v1(izero)
      !
      !v2(ii(1)) = 1.d0
      !v2(ii(2)) = p1(ii(2)) + sqrt(1.d0 - (v2(ii(1)) - p1(ii(1)))**2)
      !
      v2(ii(1)) = p1(ii(1))
      v2(ii(2)) = 1 + p1(ii(2))
      
      !
      v2(ii(1)) = v2(ii(1)) - p1(ii(1))
      v2(ii(2)) = v2(ii(2)) - p1(ii(2))
      v2(izero) = v2(izero) - p1(izero)
      !
      L    = sqrt(v2(1)*v2(1) + v2(2)*v2(2) + v2(3)*v2(3))
      invl = 1.d0/L
      !
      v2(1) = v2(1)*invl
      v2(2) = v2(2)*invl
      v2(3) = v2(3)*invl
      !
      ! find vector v3
      !
      v3(1) =   v1(2)*v2(3) - v1(3)*v2(2)
      v3(2) = -(v1(1)*v2(3) - v1(3)*v2(1))
      v3(3) =   v1(1)*v2(2) - v1(2)*v2(1)
      L    = sqrt(v3(1)*v3(1) + v3(2)*v3(2) + v3(3)*v3(3))
      invl = 1.d0/L
      !
      v3(1) = v3(1)*invl
      v3(2) = v3(2)*invl
      v3(3) = v3(3)*invl
      !
      !print*,'p1: ',p1
      !print*,'v1: ',v1
      !print*,'v2: ',v2
      !print*,'v3: ',v3
      !print*,'-----------------------------'

   else if (nzeros==1) then

      if (zeroid(3) == 1) then
         v1(1) =  (nodepos(1,node2)-nodepos(1,node1))*invl
         v1(2) =  (nodepos(2,node2)-nodepos(2,node1))*invl
         v1(3) = 0.d0

         v2(1) = -(nodepos(2,node2)-nodepos(2,node1))*invl
         v2(2) =  (nodepos(1,node2)-nodepos(1,node1))*invl
         v2(3) = 0.d0

         v3(1) = 0.d0
         v3(2) = 0.d0
         v3(3) = 1.d0
      else if (zeroid(1) == 1) then
         v1(1)  = 1.d0
         v1(2)  = 0.d0
         v1(3)  = 0.d0
               
         v2(1)  = 0.d0
         v2(2)  =  (nodepos(1,node2)-nodepos(1,node1))*invl
         v2(3)  =  (nodepos(2,node2)-nodepos(2,node1))*invl
               
         v3(1)  = 0.d0
         v3(2)  = -(nodepos(2,node2)-nodepos(2,node1))*invl
         v3(3)  =  (nodepos(1,node2)-nodepos(1,node1))*invl

      else if (zeroid(2) == 1) then
         v1(1) =  (nodepos(1,node2)-nodepos(1,node1))*invl
         v1(2) = 0.d0
         v1(3) =  (nodepos(2,node2)-nodepos(2,node1))*invl
              
         v2(1) = 0.d0
         v2(2) = 1.d0
         v2(3) = 0.d0
              
         v3(1) = -(nodepos(2,node2)-nodepos(2,node1))*invl
         v3(2) = 0.d0
         v3(3) =  (nodepos(1,node2)-nodepos(1,node1))*invl
      else
         stop 'beam3d.F90: Should not be here...'
      endif



   else
      stop 'beam3d.F90: not yet implemented'
   endif
 
   vectors(1,1) = v1(1)
   vectors(1,2) = v1(2)
   vectors(1,3) = v1(3)
   !
   vectors(2,1) = v2(1)
   vectors(2,2) = v2(2)
   vectors(2,3) = v2(3)
   !
   vectors(3,1) = v3(1)
   vectors(3,2) = v3(2)
   vectors(3,3) = v3(3)


   !vectors(:,1) = (/ v1(1), v1(2), v1(3) /)
   !vectors(:,2) = (/ v2(1), v2(2), v2(3) /)
   !vectors(:,3) = (/ v3(1), v3(2), v3(3) /)

   !print*,'vectors 1 --- ',vectors(:,1)
   !print*,'vectors 2 --- ',vectors(:,2)
   !print*,'vectors 3 --- ',vectors(:,3)
   !print*,'===================================='


   ! !
   ! px = nodepos(1,node1)
   ! py = nodepos(2,node1)
   ! pz = nodepos(3,node1)
   ! !
   ! ! x' vector
   ! !
   ! a1 = (nodepos(1,node2)-nodepos(1,node1))*invl
   ! b1 = (nodepos(2,node2)-nodepos(2,node1))*invl
   ! c1 = (nodepos(3,node2)-nodepos(3,node1))*invl
   ! !
   ! ! y' vector
   ! !
   ! if(abs(a1) < 1.d-12) a2 = 0.d0
   ! if(abs(b1) < 1.d-12) b2 = 0.d0
   ! if(abs(c1) < 1.d-12) c2 = 0.d0

   ! if (abs(a1) > 1.d-12) then
   !    b2 = 1.d0; c2 = 1.d0
   !    a2 = 1.d0/a1 * (a1*px +b1*py+c1*pz - b1*b2 - c1*c2)
   ! else if (abs(b1) > 1.d-12) then
   !    a2 = 1.d0; c2 = 1.d0
   !    b2 = 1.d0/b1 * (a1*px +b1*py+c1*pz - a1*a2 - c1*c2)
   ! else if (abs(c1) > 1.d-12) then
   !    a2 = 1.d0; b2 = 1.d0
   !    c2 = 1.d0/c1 * (a1*px +b1*py+c1*pz - a1*a2 - b1*b2)
   ! else
   !     stop 'Error: beam3d.F90. Check y vector'
   ! endif
   ! !
   ! tempsum = sqrt(a1*a2+b2*b2+c2*c2)
   ! tempsum = 1.d0/tempsum
   ! !
   ! a2 = a2*tempsum
   ! b2 = a2*tempsum
   ! c2 = a2*tempsum
   ! !
   ! print*,'------- >',px,py,pz
   ! print*,'------- >',a1,b1,c1
   ! print*,'------- >',a2,b2,c2
   ! print*,'------------------------------------------'
   ! !
   ! ! z' vector
   ! !



















end subroutine findLocalCoordinates









