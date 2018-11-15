!==================================================================================
! routine to transfer normalized coordinates to fortran
!==================================================================================

!==================================================================================
! routine to transfer normalized coordinates to fortran storage
!==================================================================================

subroutine transfer_norm_xyz(nnodes, rxyz)

   use fea_types
   implicit none 

!==================================================================================
! inputs
!==================================================================================

   integer        , intent(in)   :: nnodes 
   real(kind=8)   , intent(in)   :: rxyz(3,nnodes)

!==================================================================================
! local variables
!==================================================================================

   integer        :: i

!==================================================================================
! begin executable code
!==================================================================================

   iodata % nnodes = nnodes               ! remember # of nodes

   do i = 1,nnodes 
      iodata % norm_xyz(1,i)  = rxyz(1,i)
      iodata % norm_xyz(2,i)  = rxyz(2,i)
      iodata % norm_xyz(3,i)  = rxyz(3,i)
   end do 

!==================================================================================
! end of operations
!==================================================================================

return 
end subroutine

!==================================================================================
! routine to scale normalized coordinates 
!==================================================================================

subroutine scale_coordinates()

   use fea_types
   use layout
   implicit none 

!==================================================================================
! inputs
!==================================================================================

!==================================================================================
! local variables
!==================================================================================

   integer        :: i, n, j, il, ir
   real(kind=8)   :: scale_factor, wtl, wtr, span, rmnt 

!==================================================================================
! begin executable code
!==================================================================================

!==============================================================================
!get normalized half span
!==============================================================================

   span  = iodata % new_b * 0.5d0 * 0.1905d0/iodata % new_r

!==============================================================================
! scale the rotor mount point location with the wing span
! performed to prevent excessive bending moments at the root
!==============================================================================

   rmnt = 0.2329 

!==============================================================================
! disabled: changing the wing span to be at least rotor diameter + clearance
! reasoning is that the designs are eliminated later!
!==============================================================================

   if (span .lt. rmnt) then                  ! nondiml span > 0.23 R
      span            = 1.1d0 * rmnt 
      iodata % new_b  = iodata % new_r*span / (0.5d0*0.1905d0)
   end if 

!==================================================================================
! set location of all nodes based on nondimensional model (created in python)
!==================================================================================

   scale_factor      = iodata % new_R/iodata % ref_R
   n                 = iodata % nnodes

   do i = 1, n
      nodepos(1,i) = iodata % norm_xyz(1,i)*scale_factor
      nodepos(2,i) = iodata % norm_xyz(2,i)*scale_factor
      nodepos(3,i) = iodata % norm_xyz(3,i)*scale_factor
   end do 

!==================================================================================
! for "original nodes" (i.e. wing tips), set scale factor to span
!==================================================================================

   scale_factor      = iodata % new_b/iodata % ref_b
   n                 = iodata % norig_nodes 
   do i = 1, n 
      j              = iodata % orig_nodes(i)
      nodepos(1,j)   = iodata % norm_xyz(1,j)*scale_factor
!      write(*,*) 'modifying node #',j 
!      write(*,*) iodata % norm_xyz(1,j), scale_factor 
   end do 

!==================================================================================
! for "wing nodes" (i.e. between wing and rotor mounts for example) obtain by 
! interpolation
!==================================================================================

   scale_factor      = iodata % new_b/iodata % ref_b
   n                 = iodata % nwing_nodes 
   do i = 1, n 
      j              = iodata % iw_nodes(1,i)      ! actual node ID to interpolate
      il             = iodata % iw_nodes(2,i)      ! left reference node id 
      ir             = iodata % iw_nodes(3,i)      ! right reference node id 

      wtl            = iodata % wtw_nodes(1,i)     ! left interp weight 
      wtr            = iodata % wtw_nodes(2,i)     ! right interp weight 

!      write(*,*) 'we want for node # ',j
!      write(*,*) 'between nodes ',il,ir
      nodepos(:,j)   = nodepos(:,il)*wtl + nodepos(:,ir)*wtr
   end do 

!==================================================================================
! end of operations
!==================================================================================

return 
end subroutine

!==================================================================================
! routine to set external loads for FEA
!==================================================================================

subroutine set_loads()

   use fea_types
   use currentvalues
   use layout
   implicit none 

!==================================================================================
! inputs
!==================================================================================

!==================================================================================
! local variables
!==================================================================================

   integer        :: i, j, iflt, iw, jl, jr, n, idof
   real(kind=8)   :: l_elem, dx(3), p, q, cosx

!==================================================================================
! begin executable code
!==================================================================================

   iflt           = 1
   rhs_vec        = 0.0d0

!==================================================================================
! Note: for quad bi plane, id_thrust = 3, id_torque = 6, id_lift = 2
!==================================================================================
   
   if(nloading .lt. 1) stop 'error: must have at least one flight condition'

!==================================================================================
! set hover loads
!==================================================================================

   do i = 1,iodata % nrotors
      j                    = iodata % id_rotor(i)

!==================================================================================
! thrust loading in hover
!==================================================================================

      jl                = ndof_beam*(j-1) + iodata % id_thrust    ! location in forcing vector
      rhs_vec(jl,iflt)  = rhs_vec(jl,iflt) + iodata % Th

!==================================================================================
! torque loading in hover
!==================================================================================

      jl                = ndof_beam*(j-1) + iodata % id_torque
      rhs_vec(j,iflt)   = rhs_vec(j,iflt) + iodata % Qh * iodata % dirn(i)
   end do 

   if(nloading .le. 1) return

!==================================================================================
! set cruise loads
!==================================================================================

   iflt     = iflt + 1
   
   do i = 1,iodata % nrotors
      j     = iodata % id_rotor(i)

!==================================================================================
! thrust loading
!==================================================================================

      jl                = ndof_beam*(j-1) + iodata % id_thrust    ! location in forcing vector
      rhs_vec(jl,iflt)  = rhs_vec(jl,iflt) + iodata % Tc

!==================================================================================
! torque loading 
!==================================================================================

      jl                = ndof_beam*(j-1) + iodata % id_torque    ! location in forcing vector
      rhs_vec(jl,iflt)  = rhs_vec(jl,iflt) + iodata % Qc * iodata % dirn(i)

   end do 

!==================================================================================
!   write(*,*) '================'
!   write(*,*) 'ROTORS DONE'
!   write(*,*) '================'

!   write(*,*) '================'
!   write(*,*) 'WINGS STARTING'
!   write(*,*) '================'
!==================================================================================

   p       = iodata % L/iodata % new_b       ! wing load/span
   q       = iodata % D/iodata % new_b       ! wing drag/span (acts along -Z axis)
   do iw = 1, iodata % nwings
      n    = iodata % nwing_elem(iw)
      do i = 1, n                            ! loop over wing elements        

!==================================================================================
! find element ID and left/right node IDs for the element
!==================================================================================

         j        = iodata % iw_elem(i,iw)   ! element index

         jl       = conn(1,j)                ! left node ID
         jr       = conn(2,j)                ! right node ID 

!==================================================================================
! find element length and cosine projection along dimension 1
! the logic is that if the element left -> right direction is along +x, we have 
! -ve flap moment at the right end and +ve moment at left end due to lift   
!==================================================================================

         dx       = nodepos(:,jr) - nodepos(:,jl)                    ! vector connecting elements
         l_elem   = sqrt(dx(1)*dx(1) + dx(2)*dx(2) + dx(3)*dx(3))    ! element length
         cosx     = dx(1)/l_elem                                     ! cosine projection
            
!==================================================================================
! apply concentrated lift load on left end (set jl = offset, index = offset+2 along Y)
!==================================================================================

         jl                = ndof_beam*(jl-1) + 2
         rhs_vec(jl,iflt)  = rhs_vec(jl,iflt) + p*l_elem*0.5d0

!==================================================================================
! apply concentrated drag load on left end (offset + 3 = force along Z)
!==================================================================================

         jl                = jl + 1
         rhs_vec(jl,iflt)  = rhs_vec(jl,iflt) - q*l_elem*0.5d0

!==================================================================================
! apply concentrated moment on left end (offset + 5 = lead-lag bending moment) 
!==================================================================================

         jl                = jl + 2
         rhs_vec(jl,iflt)  = rhs_vec(jl,iflt) - q*l_elem*l_elem/12.0d0*cosx

!==================================================================================
! apply concentrated moment on left end (offset + 6 = flap bending moment) 
!==================================================================================

         jl                = jl + 1
         rhs_vec(jl,iflt)  = rhs_vec(jl,iflt) + p*l_elem*l_elem/12.0d0*cosx

!==================================================================================
! apply concentrated load on right end
!==================================================================================

         jr                = ndof_beam*(jr-1) + 2
         rhs_vec(jr,iflt)  = rhs_vec(jr,iflt) + p*l_elem*0.5d0

!==================================================================================
! apply concentrated drag load on left end (offset + 3 = force along Z)
!==================================================================================

         jr                = jr + 1
         rhs_vec(jr,iflt)  = rhs_vec(jr,iflt) - q*l_elem*0.5d0

!==================================================================================
! apply concentrated moment on left end (offset + 5 = lead-lag bending moment) 
!==================================================================================

         jr                = jr + 2
         rhs_vec(jr,iflt)  = rhs_vec(jr,iflt) + q*l_elem*l_elem/12.0d0*cosx

!==================================================================================
! apply concentrated moment on right end   ( index 2 = lift, 6 = flap bending moment) 
!==================================================================================

         jr                = jr + 1
         rhs_vec(jr,iflt)  = rhs_vec(jr,iflt) - p*l_elem*l_elem/12.0d0*cosx

      end do  
   end do 

!==================================================================================
! print out RHS vector (non-zero entries)
!==================================================================================
   
   ! do iflt = 1, nloading
   !    write(*,*) 'FLIGHT CONDITION # ',iflt
   !    do i = 1, nnodes
   !       do j = 1,ndof_beam
   !          jl    = (i-1)*ndof_beam + j         ! 
   !          if(abs(rhs_vec(jl,iflt)) .ge. 1.d-5) then 
   !             write(*,*) 'Node',i,'DOF',j,'Loading',rhs_vec(jl,iflt)
   !          end if
   !       end do  
   !    end do 
   ! end do 
   ! stop 'ok?'
!==================================================================================
! end of operations
!==================================================================================

return 
end subroutine