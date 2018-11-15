!
! find global stiffness and mass matrices for a 3-dof 3d beam element 
!
subroutine beam3d(elemid,K_loc,M_loc)

   use layout
   use material
   use fea_types
   implicit none


   integer,       intent(in) :: elemid
   real(kind=8), intent(out) :: K_loc(12,12), M_loc(12,12)


   integer :: i,k,off
   integer :: j1,j2,j3
   integer :: matid,node1,node2
   real(kind=8) :: l1,l2,l3,m1,m2,m3,n1,n2,n3
   real(kind=8) :: L,invl
   real(kind=8) :: T(12,12),mtx1(12,12)
   real(kind=8) :: X,Y1,Y2,Y3,Y4,Z1,Z2,Z3,Z4,S
   real(kind=8) :: sum,tempsum

   real(kind=8) :: area,E,G,J,phi_z,phi_y,Iy,Iz,ky,kz,rho,nu

   real(kind=8) :: vectors(3,3),temp

   ! get material properties for this particular element
   E    = matprop(1,elemid)
   nu   = matprop(2,elemid)
   rho  = matprop(3,elemid)
   area = matprop(4,elemid)
   Iy   = matprop(5,elemid)
   Iz   = matprop(6,elemid)
   J    = matprop(7,elemid)
   ky   = matprop(8,elemid)
   kz   = matprop(9,elemid)
   !
   G = E/(2.d0*(1.d0+nu))
   !
   node1 = conn(1,elemid)
   node2 = conn(2,elemid)
   !
   L = (nodepos(1,node2)-nodepos(1,node1))*(nodepos(1,node2)-nodepos(1,node1)) &
     + (nodepos(2,node2)-nodepos(2,node1))*(nodepos(2,node2)-nodepos(2,node1)) &
     + (nodepos(3,node2)-nodepos(3,node1))*(nodepos(3,node2)-nodepos(3,node1))
   L = sqrt(L)
   !
   invl   = 1.d0/L
   !
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
   
   !l1 =  (nodepos(1,node2)-nodepos(1,node1))*invl
   !l2 = -(nodepos(2,node2)-nodepos(2,node1))*invl
   !l3 = 0.d0

   !m1 =  (nodepos(2,node2)-nodepos(2,node1))*invl
   !m2 =  (nodepos(1,node2)-nodepos(1,node1))*invl
   !m3 = 0.d0

   !n1 = 0.d0
   !n2 = 0.d0
   !n3 = 1.d0


   ! 
   ! define 12x12 transformation matrix
   ! 
   T = 0.d0
   do i = 1,4
      off = 3*(i-1)
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
! define stiffness matrix terms
!
!   phi_y = 12.d0*E*Iz*ky/(area*G*L*L)
!   phi_z = 12.d0*E*Iz*kz/(area*G*L*L)

   temp     = 12.d0*E/(area*G*L*L)
   phi_y = Iz*ky*temp!/(area*G*L*L)
   phi_z = Iz*kz*temp!/(area*G*L*L)

   !
   X     = E*area*invL                 ! EA/L

!   Y1    = 12.d0*E*Iz/((1.d0+phi_y)*L*L*L)
!   Y2    =  6.d0*E*Iz/((1.d0+phi_y)*L*L)
!   Z1    = 12.d0*E*Iy/((1.d0+phi_z)*L*L*L)
!   Z2    =  6.d0*E*Iy/((1.d0+phi_z)*L*L)

   temp  = E*Iz/((1.d0+phi_y)*L*L)
   Y2    = 6.d0*temp
   Y1    = 2.0d0*Y2*invL 

   temp  = temp*L
   Y3    = (4.d0+phi_y)*temp      !temp = E*Iz/((1.d0+phi_y)*L)
   Y4    = (2.d0-phi_y)*temp      !E*Iz/((1.d0+phi_y)*L)


   temp  = E*Iy/((1.d0+phi_z)*L*L)
   Z2    = 6.d0*temp              !temp = E*Iy/((1.d0+phi_z)*L*L)
   Z1    = 2.d0*Z2*invL
   
   temp  = temp*L
   Z3    = (4.d0+phi_z)*temp      !temp = E*Iy/((1.d0+phi_z)*L)
   Z4    = (2.d0-phi_z)*temp      !E*Iy/((1.d0+phi_z)*L)
   S     = G*J*invL               !Linv = 1/L
   !
   ! define local stiffness matrix
   !
   K_loc( 1,:) = (/    X, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,   -X, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
   K_loc( 2,:) = (/ 0.d0,   Y1, 0.d0, 0.d0, 0.d0,   Y2, 0.d0,  -Y1, 0.d0, 0.d0, 0.d0,   Y2 /)
   K_loc( 3,:) = (/ 0.d0, 0.d0,   Z1, 0.d0,  -Z2, 0.d0, 0.d0, 0.d0,  -Z1, 0.d0,  -Z2, 0.d0 /)
   K_loc( 4,:) = (/ 0.d0, 0.d0, 0.d0,    S, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,   -S, 0.d0, 0.d0 /)
   K_loc( 5,:) = (/ 0.d0, 0.d0,  -Z2, 0.d0,   Z3, 0.d0, 0.d0, 0.d0,   Z2, 0.d0,   Z4, 0.d0 /)
   K_loc( 6,:) = (/ 0.d0,   Y2, 0.d0, 0.d0, 0.d0,   Y3, 0.d0,  -Y2, 0.d0, 0.d0, 0.d0,   Y4 /)
   K_loc( 7,:) = (/   -X, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,    X, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
   K_loc( 8,:) = (/ 0.d0,  -Y1, 0.d0, 0.d0, 0.d0,  -Y2, 0.d0,   Y1, 0.d0, 0.d0, 0.d0,  -Y2 /)
   K_loc( 9,:) = (/ 0.d0, 0.d0,  -Z1, 0.d0,   Z2, 0.d0, 0.d0, 0.d0,   Z1, 0.d0,   Z2, 0.d0 /)
   K_loc(10,:) = (/ 0.d0, 0.d0, 0.d0,   -S, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,    S, 0.d0, 0.d0 /)
   K_loc(11,:) = (/ 0.d0, 0.d0,  -Z2, 0.d0,   Z4, 0.d0, 0.d0, 0.d0,   Z2, 0.d0,   Z3, 0.d0 /)
   K_loc(12,:) = (/ 0.d0,   Y2, 0.d0, 0.d0, 0.d0,   Y4, 0.d0,  -Y2, 0.d0, 0.d0, 0.d0,   Y3 /)

!
! calculate the global stiffness matrix
! 
! mtx1 = K_loc*T
   mtx1=0.d0
   do j1 = 1,12
      do j2 = 1,12
         sum=0.d0
         do j3 = 1,12
            sum = sum + K_loc(j3,j1)*T(j3,j2)      ! use symmetry to switch j3,j1
         enddo
         mtx1(j1,j2) = sum
      enddo
   enddo
!
! K_loc = T' * mtx1
!
   do j1 = 1,12
      do j2 = 1,12
         sum=0.d0
         do j3 = 1,12
            sum = sum + T(j3,j1)*mtx1(j3,j2)
         enddo
         K_loc(j1,j2) = sum
      enddo
   enddo

!=========================================================================================
! perform mass matrix calculations only if modes flag is on
!=========================================================================================

   if(iodata % get_modes) then

   !
   ! define local mass matrix
   !
      temp        = J/area*70.d0 
      Y1          = 4.d0*L*L 
      Y2          = 22.d0*L
      Y3          = 13.d0*L
      Y4          =  3.d0*L*L
      M_loc( 1,:) = (/ 140.d0,     0.d0,     0.d0,          0.d0,      0.d0,      0.d0,  70.d0,&
                         0.d0,     0.d0,          0.d0,      0.d0,      0.d0  /)
      M_loc( 2,:) = (/   0.d0,   156.d0,     0.d0,          0.d0,      0.d0,   Y2,   0.d0,&
                        54.d0,     0.d0,          0.d0,      0.d0,  -Y3  /)
      M_loc( 3,:) = (/   0.d0,     0.d0,   156.d0,          0.d0,  -Y2,      0.d0,   0.d0,&
                         0.d0,    54.d0,          0.d0,   Y3,      0.d0  /)
      M_loc( 4,:) = (/   0.d0,     0.d0,     0.d0, 2.d0*temp,        0.d0,      0.d0,   0.d0,&
                         0.d0,     0.d0,  temp,        0.d0,      0.d0  /)
      M_loc( 5,:) = (/   0.d0,     0.d0, -Y2,          0.d0,  Y1,      0.d0,   0.d0,&
                         0.d0, -Y3,          0.d0, -Y4,      0.d0  /)
      M_loc( 6,:) = (/   0.d0,  Y2,     0.d0,          0.d0,      0.d0,  Y1,   0.d0,&
                      Y3,     0.d0,          0.d0,      0.d0, -Y4  /)
      M_loc( 7,:) = (/  70.d0,     0.d0,     0.d0,          0.d0,      0.d0,      0.d0, 140.d0,& 
                         0.d0,     0.d0,          0.d0,      0.d0,      0.d0  /)
      M_loc( 8,:) = (/   0.d0,    54.d0,     0.d0,          0.d0,      0.d0,   Y3,   0.d0,&
                       156.d0,     0.d0,          0.d0,      0.d0,  -Y2  /)
      M_loc( 9,:) = (/   0.d0,     0.d0,    54.d0,          0.d0,  -Y3,      0.d0,   0.d0,&
                         0.d0,   156.d0,          0.d0,   Y2,      0.d0  /)
      M_loc(10,:) = (/   0.d0,     0.d0,     0.d0,  temp,        0.d0,      0.d0,   0.d0,& 
                         0.d0,     0.d0, 2.d0*temp,        0.d0,      0.d0  /)
      M_loc(11,:) = (/   0.d0,     0.d0,  Y3,          0.d0, -Y4,      0.d0,   0.d0,&
                         0.d0,  Y2,          0.d0,  Y1,      0.d0  /)
      M_loc(12,:) = (/   0.d0, -Y3,     0.d0,          0.d0,      0.d0, -Y4,   0.d0,& 
                     -Y2,     0.d0,          0.d0,      0.d0,  Y1  /)
   ! 
      do i=1,12
         do k=1,12
            M_loc(i,k) = M_loc(i,k)*rho*area*L/420.d0
         enddo
      enddo
!
! calculate the global mass matrix
!
! mtx1 = M_loc*T
      mtx1=0.d0
      do j1 = 1,12
         do j2 = 1,12
            sum=0.d0
            do j3 = 1,12
               sum = sum + M_loc(j1,j3)*T(j3,j2)
            enddo
            mtx1(j1,j2) = sum
         enddo
      enddo
!
! M_loc = T' * mtx1
!
      do j1 = 1,12
         do j2 = 1,12
            sum=0.d0
            do j3 = 1,12
               sum = sum + T(j3,j1)*mtx1(j3,j2)
            enddo
            M_loc(j1,j2) = sum
         enddo
      enddo

   end if 

!=========================================================================================
! end of operations
!=========================================================================================

end subroutine beam3d
