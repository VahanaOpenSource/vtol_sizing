! !==============================================================================
! ! this subroutine creates the airframe configuration based on quad bi plane 
! ! parameters
! !==============================================================================

!   subroutine create_airframe(params)

!   use fea_types
!   implicit none 

! !==============================================================================
! ! inputs
! !==============================================================================
   
!   type(airframe_params_def), intent(inout) :: params
  
! !==============================================================================
! ! local variables
! !==============================================================================
  
!   real(kind=8)  :: radius, span, rmnt, wingspan, thrust, hvr_torque, ff_torque, &
!                    drag, winglift, bbyd, bbyr, spanbyr, c1

!   real(kind=8)  :: rxyz0(21,3)        ! location of points (nondiml)

!   logical       :: ihover, iforward
!   integer       :: base, nbeamelem
!   integer       :: conn0(26,3)

! !for applying loads to nodes
!   integer             :: maxnforce
 
!   integer             :: force_id(2,nmax_loads,nmax_flight)    ! nodes, dof, flight conditions
!   integer             :: nloads(nmax_flight)                   ! loads per flight condition
!   real(kind=8)        :: ploads(nmax_loads,nmax_flight)        ! point loads at nodes 

!   integer             :: idf, nrotors, rotor_nodes(nrmax)
!   real(kind=8)        :: dirn(nrmax)

! !to apply boundary conditions 
!   integer             :: dispid(6,2)     ! max 6 BCs, i.e. one node 

! !find total # nodes 
!   integer             :: nnodes, nnodes0, nelem, nelem0

! !==============================================================================
! !begin executable code: interpret input parameters
! !==============================================================================

!   radius     = params % R
!   wingspan   = params % b
!   thrust     = params % T
!   hvr_torque = params % Qh
!   ff_torque  = params % Qc
!   drag       = params % D
!   winglift   = params % L 

! !==============================================================================
! ! interpret inputs 
! !==============================================================================

!   BASE       = 1
!   nbeamelem  = 3
!   ihover     = .true.
!   iforward   = .true.

! !==============================================================================
! ! scale the wing size(b/2) based on the wingspan
! !==============================================================================

!   bbyd       = 1.1d0
!   bbyr       = 0.1905d0
!   spanbyr    = 0.5d0*bbyr
!   span       = wingspan/radius * spanbyr      ! normalized span

! !==============================================================================
! ! scale the rotor mount point location with the wing span
! ! performed to prevent excessive bending moments at the root
! !==============================================================================

!   rmnt       = 0.2329d0  

! !==============================================================================
! !wingspan_new = wingspan
! !==============================================================================

!   if (span .lt. bbyd*rmnt) then        ! don't use a smaller span/radius than AIREZ
      
! !==============================================================================
! ! re-size span
! !==============================================================================

!     span         = bbyd * rmnt            ! normalized span/radius 
!     wingspan     = radius*span / spanbyr  ! dimensionalize using radius
!     params % b   = wingspan
!   end if 

! !==============================================================================
! ! location of nodes (normalized by radius) - quad bi 
! !==============================================================================

!   c1          = 0.09d0
!   c2          = 0.21d0
!   c3          = 0.14d0

!   rxyz0( 1,1) = -span;  rxyz0( 1,2) = rmnt;  rxyz0( 1,3) = 0.0d0
!   rxyz0( 2,1) = -rmnt;  rxyz0( 2,2) = rmnt;  rxyz0( 2,3) = 0.0d0
!   rxyz0( 3,1) =  rmnt;  rxyz0( 3,2) = rmnt;  rxyz0( 3,3) = 0.0d0
!   rxyz0( 4,1) =  span;  rxyz0( 4,2) = rmnt;  rxyz0( 4,3) = 0.0d0
!   rxyz0( 5,1) =  0.0d0; rxyz0( 5,2) = 0.d0;  rxyz0( 5,3) = 0.0d0
!   rxyz0( 6,1) = -span;  rxyz0( 6,2) =-rmnt;  rxyz0( 6,3) = 0.0d0
!   rxyz0( 7,1) = -rmnt;  rxyz0( 7,2) =-rmnt;  rxyz0( 7,3) = 0.0d0
!   rxyz0( 8,1) =  rmnt;  rxyz0( 8,2) =-rmnt;  rxyz0( 8,3) = 0.0d0
!   rxyz0( 9,1) =  span;  rxyz0( 9,2) =-rmnt;  rxyz0( 9,3) = 0.0d0
!   rxyz0(10,1) =   -c1;  rxyz0(10,2) =   c1;  rxyz0(10,3) = 0.0d0
!   rxyz0(11,1) =    c1;  rxyz0(11,2) =   c1;  rxyz0(11,3) = 0.0d0
!   rxyz0(12,1) =   -c1;  rxyz0(12,2) =  -c1;  rxyz0(12,3) = 0.0d0
!   rxyz0(13,1) =    c1;  rxyz0(13,2) =  -c1;  rxyz0(13,3) = 0.0d0
!   rxyz0(14,1) =   -c1;  rxyz0(14,2) =   c1;  rxyz0(14,3) = -c2
!   rxyz0(15,1) =    c1;  rxyz0(15,2) =   c1;  rxyz0(15,3) = -c2
!   rxyz0(16,1) =   -c1;  rxyz0(16,2) =  -c1;  rxyz0(16,3) = -c2
!   rxyz0(17,1) =    c1;  rxyz0(17,2) =  -c1;  rxyz0(17,3) = -c2
!   rxyz0(18,1) = -rmnt;  rxyz0(18,2) = rmnt;  rxyz0(18,3) =  c3
!   rxyz0(19,1) =  rmnt;  rxyz0(19,2) = rmnt;  rxyz0(19,3) =  c3
!   rxyz0(20,1) = -rmnt;  rxyz0(20,2) =-rmnt;  rxyz0(20,3) =  c3
!   rxyz0(21,1) =  rmnt;  rxyz0(21,2) =-rmnt;  rxyz0(21,3) =  c3

! !==============================================================================
! ! scale up based on AirEZ dimensions
! !==============================================================================

!   rxyz0      = rxyz0*radius/bbyr

! !==============================================================================
! ! define element connectivity
! !==============================================================================

!   conn0( 1,1)  =  1; conn0( 1,2) =  2; conn0( 1,3) = 2
!   conn0( 2,1)  =  2; conn0( 2,2) =  3; conn0( 2,3) = 0
!   conn0( 3,1)  =  3; conn0( 3,2) =  4; conn0( 3,3) = 2
!   conn0( 4,1)  =  6; conn0( 4,2) =  7; conn0( 4,3) = 2
!   conn0( 5,1)  =  7; conn0( 5,2) =  8; conn0( 5,3) = 0
!   conn0( 6,1)  =  8; conn0( 6,2) =  9; conn0( 6,3) = 2
!   conn0( 7,1)  =  2; conn0( 7,2) = 10; conn0( 7,3) = 2
!   conn0( 8,1)  = 10; conn0( 8,2) =  5; conn0( 8,3) = 2
!   conn0( 9,1)  =  5; conn0( 9,2) = 13; conn0( 9,3) = 2
!   conn0(10,1)  = 13; conn0(10,2) =  8; conn0(10,3) = 2
!   conn0(11,1)  =  3; conn0(11,2) = 11; conn0(11,3) = 2
!   conn0(12,1)  = 11; conn0(12,2) =  5; conn0(12,3) = 2
!   conn0(13,1)  =  5; conn0(13,2) = 12; conn0(13,3) = 2
!   conn0(14,1)  = 12; conn0(14,2) =  7; conn0(14,3) = 2
!   conn0(15,1)  = 10; conn0(15,2) = 14; conn0(15,3) = 0
!   conn0(16,1)  = 11; conn0(16,2) = 15; conn0(16,3) = 0
!   conn0(17,1)  = 12; conn0(17,2) = 16; conn0(17,3) = 0
!   conn0(18,1)  = 13; conn0(18,2) = 17; conn0(18,3) = 0
!   conn0(19,1)  = 18; conn0(19,2) =  2; conn0(19,3) = 2
!   conn0(20,1)  = 19; conn0(20,2) =  3; conn0(20,3) = 2
!   conn0(21,1)  = 20; conn0(21,2) =  7; conn0(21,3) = 2
!   conn0(22,1)  = 21; conn0(22,2) =  8; conn0(22,3) = 2
!   conn0(23,1)  = 14; conn0(23,2) = 15; conn0(23,3) = 2
!   conn0(24,1)  = 15; conn0(24,2) = 17; conn0(24,3) = 2
!   conn0(25,1)  = 17; conn0(25,2) = 16; conn0(25,3) = 2
!   conn0(26,1)  = 16; conn0(26,2) = 14; conn0(26,3) = 2

!   nnodes       = 
! !==============================================================================
! ! types of different loadings to design for (hover, forward flight)
! !==============================================================================

!   if(ihover .and. (.not. iforward)) then                  
!     nloadings = 1
!   elseif( (.not. ihover) .and. iforward) then 
!     nloadings = 1
!   elseif (ihover .and. iforward) then 
!     nloadings = 2
!   else
!     write(*,*) 'Incorrect loading condition...',ihover,iforward
!     stop 'critical error : must have at least one loading condition'
!   end if 

! !==============================================================================
! ! determine where the rotors are
! !==============================================================================

!   nrotors     = 4               ! make it input later!
!   if (nrotors .gt. nrmax) stop 'airframe_config.F90 in FEA: increase rotor storage'

!   rotor_nodes(1)  = 18    ; dirn(1)   = 1.d0      ! ccw rotor 
!   rotor_nodes(2)  = 19    ; dirn(2)   =-1.d0      ! cw  rotor 
!   rotor_nodes(3)  = 20    ; dirn(3)   =-1.d0      ! cw  rotor 
!   rotor_nodes(4)  = 21    ; dirn(4)   = 1.d0      ! ccw rotor 

! !==============================================================================
! ! store loads in all flight conditions
! !==============================================================================

!   iflt          = 0               ! flight condition counter

! !==============================================================================
! ! force/moment loading information is stored in an array with 
! ! 2 columns per entry (nodeid, dofid, flight condition id)
! ! forcing stored separately (loads, flight condition id)
! !==============================================================================

! !==============================================================================
! ! apply rotor hub loads: thrust and torque
! !==============================================================================

!   if(ihover) then 
!     iflt        = iflt + 1
!     idf         = 0               ! loading counter

!     do i = 1, nrotors
!       idf                   = idf + 1
!       force_id(1,idf,iflt)  = rotor_nodes(i)
!       force_id(2,idf,iflt)  = 3
!       ploads(idf,iflt)      = thrust

!       idf                   = idf + 1
!       force_id(1,idf,iflt)  = rotor_nodes(i)
!       force_id(2,idf,iflt)  = 6
!       ploads(idf,iflt)      = hvr_torque*dirn(i)
!     end do 
!   end if 

! !==============================================================================
! ! remember loads in cruise
! !==============================================================================

!   if(iforward) then 
!     iflt        = iflt + 1
!     idf         = 0               ! loading counter

!     do i = 1, nrotors
!       idf                   = idf + 1
!       force_id(1,idf,iflt)  = rotor_nodes(i)
!       force_id(2,idf,iflt)  = 3
!       ploads(idf,iflt)      = drag

!       idf                   = idf + 1
!       force_id(1,idf,iflt)  = rotor_nodes(i)
!       force_id(2,idf,iflt)  = 6
!       ploads(idf,iflt)      = ff_torque*dirn(i)
!     end do 
!   end if 

!   nloading    = iflt        ! keep track of how many flight conditions there are

! !==============================================================================
! !
! !==============================================================================
      
! !flag to compute nodal loads 
!   integer         :: nloads(nmax_flight)
!   integer         :: id_mem(nmax_dist_mem,nmax_dist,nmax_flight)
!   real(kind=rdp)  :: load_val(nmax_dist,nmax_flight)

!   dist_flag(1)    = .false.       ! no distributed loads in hover
!   dist_flag(2)    = .true.        ! yes distributed loads in cruise from wing

! !==============================================================================
! ! distributed loads (member ID, dofid, load (N/m))
! !==============================================================================
   
! ! forward flight
!   if(iforward) then 

!     iflt              = 2
!     load_val(1,iflt)  = winglift      ! first  dist load, second flight condition
!     load_val(2,iflt)  = winglift      ! second dist load, second flight condition

!     imem  = 0             ! member id 
!     iload = 1             ! dist load id 

!     imem  = imem + 1;   id_mem(imem,iload,iflt)   = 1
!     imem  = imem + 1;   id_mem(imem,iload,iflt)   = 2
!     imem  = imem + 1;   id_mem(imem,iload,iflt)   = 3

!     distloadid = np.array([
!                            [[1,2,3],2,winglift], ! simulating lift from wing
!                            [[4,5,6],2,winglift]
!                          ])
    
!     ! append distributed load
!     alldistload.append(distloadid)
!   end if 
! !====================================================================
! ! set # of boundary conditions
! ! also specify the constrained DOF (column 1) and the node id (col 2)
! !====================================================================

!   do i = 1,6
!     dispid(i,1) = 5 
!     dispid(i,2) = i 
!   end do 

! ! the rest of the code takes the above inputs and creates
! ! the input file necessary for the fea code
! ! no user intervention required
!    !
   
!    ! ===================================================================
!    ! given the above rxyz and conn
!    ! create the node and conn array assuming nbeamelem
!    ! ===================================================================
   
!    conn0 -= BASE ! convert to base 0 from base 1
   
!    nnodes0 = len(rxyz0)
!    nelem0  = len(conn0)
   
!    nnodes = nnodes0 + nelem0*(nbeamelem-1)
!    nelem  = nelem0 * nbeamelem
   
!    rxyz     = np.zeros([nnodes,3])
!    conn     = np.zeros([nelem,2])
!    memberid = np.zeros(nelem)
!    !
!    ! copy positions
!    !
!    for j in range(nnodes0):
!       rxyz[j][:] = rxyz0[j][:]
!    !
!    ! create member to node array
!    !
!    maxconn      = 10
!    nmembers     = len(conn0)
!    member2node  = np.empty([nmembers,maxconn])
!    member2node.fill(-1)
!    nmember2node = np.zeros(nmembers)
!    !
!    ! create the position and conn arrays 
!    !
!    inode = nnodes0
!    ielem = 0
   
!    for j in range(nelem0):
!       n1 = conn0[j][0]
!       n2 = conn0[j][1]
!       !
!       r1 = rxyz0[n1][:]
!       r2 = rxyz0[n2][:]
   
!       if(conn0[j][2] == -1):
!          nn = nbeamelem - 1
!       else:
!          nn = conn0[j][2]
   
!    !
!    ! 1 beam element
!    !
!       if(nn==0):
         
!          memid           = j
   
!          conn[ielem][0]  = int(n1)
!          conn[ielem][1]  = int(n2)
!          memberid[ielem] = memid
         
!          ! add conn[ielem][0] to member2node
!          if (conn[ielem][0] in member2node[memid][:]) == False:
!             member2node[memid][nmember2node[memid]] = conn[ielem][0]
!             nmember2node[memid] += 1
   
!          ! add conn[ielem][1] to member2node
!          if (conn[ielem][1] in member2node[memid][:]) == False:
!             member2node[memid][nmember2node[memid]] = conn[ielem][1]
!             nmember2node[memid] += 1
   
!          ielem += 1
   
!    ! 
!    ! more than 1 beam element
!    ! 
!       else:
!          ! add additional points
!          for k in range(nn):
!             rxyz[inode][0] = r1[0] + (r2[0]-r1[0])/(nn+1) * (k+1)
!             rxyz[inode][1] = r1[1] + (r2[1]-r1[1])/(nn+1) * (k+1)
!             rxyz[inode][2] = r1[2] + (r2[2]-r1[2])/(nn+1) * (k+1)
            
!             nn1 = inode - 1
!             nn2 = inode + 0
   
!             if(k==0):
!                nn1 = n1
   
!             conn[ielem][0]  = int(nn1)
!             conn[ielem][1]  = int(nn2)
!             memberid[ielem] = j
            
!             memid = int(memberid[ielem])
   
!             ! add conn[ielem][0] to member2node
!             if (conn[ielem][0] in member2node[memid][:]) == False:
!                member2node[memid][int(nmember2node[memid])] = conn[ielem][0]
!                nmember2node[memid] += 1
   
!             ! add conn[ielem][1] to member2node
!             if (conn[ielem][1] in member2node[memid][:]) == False:
!                member2node[memid][int(nmember2node[memid])] = conn[ielem][1]
!                nmember2node[memid] += 1
   
   
!             ielem += 1
!             inode += 1
   
!          conn[ielem][0]  = int(inode-1)
!          conn[ielem][1]  = int(n2)
!          memberid[ielem] = j
   
!          memid = int(memberid[ielem])
   
!          ! add conn[ielem][0] to member2node
!          if (conn[ielem][0] in member2node[memid][:]) == False:
!             member2node[memid][int(nmember2node[memid])] = conn[ielem][0]
!             nmember2node[memid] += 1
   
!          ! add conn[ielem][1] to member2node
!          if (conn[ielem][1] in member2node[memid][:]) == False:
!             member2node[memid][int(nmember2node[memid])] = conn[ielem][1]
!             nmember2node[memid] += 1
   
!          ielem += 1
   
!    ! 
!    ! loop over the number of loading conditins
!    !
! !   forceOut = []     ! AS converted it to a dictionary

!    maxnforce      = 0
!    nforces        = []
     
!    for n in range(nloadings):
   
! !      print 'n is ',n 
!       forceid    = allforce[n]
!       distloadid = alldistload[n]
   
!       !
!       ! logic to find the nodes to apply the distributed loads
!       ! 
!       ! loop through the number of distributed loads
!       ndistload = len(distloadid)

! !      print 'n is ',n,' at the top-ish'

!       for j in range(ndistload):
      
!          member_length = 0.0
!          nmember_load  = len(distloadid[j][0])
!          nodelist      = np.zeros(0,dtype=int)

! !         print 'n is ',n,' at the top-middle'
      
!          ! loop through the members for each distributed load
!          for k in range(nmember_load):
!             memid = int(distloadid[j][0][k]) - BASE
      
!             ! add length to member
!             n1 = conn0[memid][0]
!             n2 = conn0[memid][1]
!             member_length += np.sqrt(  (rxyz[n1,0]-rxyz[n2,0])*(rxyz[n1,0]-rxyz[n2,0])  \
!                                      + (rxyz[n1,1]-rxyz[n2,1])*(rxyz[n1,1]-rxyz[n2,1])  \
!                                      + (rxyz[n1,2]-rxyz[n2,2])*(rxyz[n1,2]-rxyz[n2,2]))
      
!             ! find a list of all nodes belonging to the members
!             for node in member2node[memid][:]:
!                node = int(node)
!                if (node >= 0 and ((node in nodelist)==False)):
!                   nodelist = np.append(nodelist,node)
      
!          ! apply the distributed load (defined as force/length)
!          ! to the nodes in nodelist. Must multiply by the length
!          ! of the member to ensure compatible dimensions
!          !print 'nodelist',nodelist
!          !for node in nodelist:
!          !   print 'node',node,'position',rxyz[node,:]
!          !
!          ! currently nodes are arranged in ascending order of 'x'
!          ! which is convienient..
!          !
!          nnodes = len(nodelist)
!          midpoints = np.zeros([nnodes-1,3])
!          for kk in range(nnodes-1):
!             n1 = nodelist[kk]
!             n2 = nodelist[kk+1]
!             midpoints[kk,:] = 0.5*(rxyz[n1,:] + rxyz[n2,:])
         
!          integ_limits = np.zeros([nnodes,2])

!          !print '\n\nmidpoints',midpoints

!          !
!          ! find limits of integration
!          !
!          for kk in range(nnodes):

!             if kk==0:
!                n1 = nodelist[kk]
!                xstart = rxyz[n1,:]
!             else:
!                xstart = midpoints[kk-1,:]
!             !
!             if kk == nnodes-1:
!                n1 = nodelist[kk]
!                xend = rxyz[n1,:]
!             else:
!                xend   = midpoints[kk,:]
            
!             integ_limits[kk,0] = xstart[0]
!             integ_limits[kk,1] = xend[0]
!          !print '\n\ninteglimits',integ_limits
!          !
!          ! find area under curve
!          !
! !         print 'n is ',n,' at the middle'
!          npts = 10
!          totarea = 0
!          area_node = np.zeros(nnodes)
!          for kk in range(nnodes):
!             xstart = integ_limits[kk,0]
!             xend   = integ_limits[kk,1]
            
!             deltax = (xend-xstart)/npts
            
!             ! rectangular rule
!             area = 0
!             for n1 in range(npts):
!                x1 = xstart + n1 * deltax
!                x2 = x1 + deltax
!                xmid = 0.5*(x1+x2)
               
!                !print 'TERM..',2*xmid/wingspan
!                area += np.sqrt(1.0 - (2*xmid/wingspan)**2) * deltax

!             area_node[kk] = area
!             totarea += area
!          !
!          ! find L0
!          !
!          L0 = 4*winglift/(wingspan*np.pi)
!          !
!          dofid            = distloadid[j][1]
!          force_per_length = distloadid[j][2]
               
!          ! elliptical distribution of force
!          kk = -1
!          temp =0
!          for node in nodelist:
!             kk+=1
!             force_in_node = L0 * area_node[kk]
!             forceid = np.vstack((forceid,(int(node)+BASE,int(dofid),force_in_node)))
!             temp += force_in_node

! !      print 'n is ',n,' at the bottom'
!       forceOut[str(n)] = forceid

!       nforce      = len(forceid)
!       nforces.append(nforce) 
!       maxnforce   = max(maxnforce,nforce)
   
!    nforces = np.asarray(nforces)

! ! ===================================================================
! ! write out to file
! ! ===================================================================
   
!    !convert back to base 1
!    conn     += BASE
!    memberid += BASE
   
!    nnodes         = inode
!    nelem          = ielem
!    !     
!    ndisp = len(dispid)
