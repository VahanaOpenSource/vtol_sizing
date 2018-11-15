! ###################################################################
!
! update cross sectional properties
! 
! ###################################################################

subroutine updateStructure(iconstraint)

   use layout
   use currentValues
   use material

   implicit none
   !
   integer, intent(out) :: iconstraint


   real(kind=8), parameter :: par_FOS  = 1.5d0
   real(kind=8), parameter :: band_FOS = 0.1d0
   !
   !real(kind=8), parameter :: par_DEF  = 10.d0!0.15d0
   real(kind=8) :: par_DEF
   !
   integer i,k,kk,elemID
   integer :: iflag(2)

   real(kind=8)   :: outerDim1, outerDim2,newOuterDim,maxy,sigYield
   real(kind=8)   :: matpropTmp(nmatprop)

   real(kind=8)   :: minPAR_FOS,maxPAR_FOS, drmax

   idefFlag=0
   istressFlag=0

! ===================================================================
! set maximum allowable deflection from rotor radius set in 
! update_airframe_weight.py in Python/Stage_3/..
! ===================================================================

   par_DEF     = maxAllowableDeflection       ! 15% of radius 

! ===================================================================
! min and max FOS
! ===================================================================

   minPAR_FOS = par_FOS - band_FOS
   maxPAR_FOS = par_FOS + band_FOS

! ===================================================================
! copying common material properties
! ===================================================================

   matpropTmp     = matprop(:,1)

! ===================================================================
! maximum deflection ratio and scaling factor for CS dimension
! ===================================================================

   drmax          = maxval(memDeflection(1:nmember))/par_DEF

   if(drmax .ge. 0.995d0 .or. drmax .le. 0.65d0) then  ! dont make it super stiff either
      drmax          = drmax**0.25d0 
   else 
      drmax          = 1.0d0
   end if 

! ===================================================================
! loop over the members
! ===================================================================

   do i = 1,nmember
   
! ===================================================================
! save maxy and sigYield for this member
! pick any element (1st for instance)
! ===================================================================

      elemID         = elem2member(1,i)
      maxy           = matprop(10,elemID)
      sigYield       = matprop(12,elemID)

! ===================================================================
! new outer radius
! ===================================================================

      outerDim1      = 0.d0
      outerDim2      = 0.d0
      iflag          = 0
     
      newOuterDim    = maxy

! ===================================================================
! check for factor of safety wrt yield stress
! NOTE: maxy * (memStress(i) * par_FOS/sigYield) ** 0.25
! should be 0.33d0 barry
! ===================================================================

      ! if(memFOS(i) < minPAR_FOS) then
      !    outerDim1   = 1.1d0 * maxy
      !    iflag(1)    = 1
      ! else if (memFOS(i) > maxPAR_FOS) then

      !    ! outer diameter cannot be smaller than 10 mm
      !    outerDim1 = max(0.9d0 * maxy, 0.01d0)
      !    !outerDim1 =0.9d0 * maxy
      !    iflag(1)  = 1
      ! endif

!====================================================================
!new code : AS
! we're solving statics problem. Factor of safety for prescribed 
! external loads is proportional to CS dimension cubed, 
! and stress is inversely proportional to cube of CS dimension for 
! hollow circle/solid square with assumed wall thickness/side ratio 
!====================================================================

      if (abs(memFOS(i) - par_FOS) .ge. band_FOS) then 
         outerDim1 = par_FOS/memFOS(i)
         outerDim1 = maxy*(outerDim1**(1.0d0/3.0d0))
         iflag(1)  = 1
      else 
         outerDim1 = maxy 
      end if  

!====================================================================
! track member deflection/allowed defln ratio   
!====================================================================

      outerdim2    = maxy*drmax 

! ===================================================================
! if neither of the above conditions are violated, we're ok 
! otherwise change the cross section
! ===================================================================
      
      newOuterDim  = max(outerDim1,outerDim2)
!      write(*,*) i,newOuterDim,memFOS(i),drmax
! ===================================================================
! update cross sectional properties
! ===================================================================

      call crossSectionProperties(newOuterDim,matpropTmp)      

! ===================================================================
! copy new cross sectional information to all
! elements of a given member
! ===================================================================

      do k = 1,nelem2member(i)
         elemID = elem2member(k,i)
         matprop(:,elemID) = matpropTmp(:)
      enddo

! ===================================================================
! set idefFlag and istressFlag for the given member
! ===================================================================

      idefFlag(i)=iflag(1)
      if(memFOS(i) > maxPAR_FOS) istressFlag(i)=1

   enddo

! ===================================================================
! check iconstraint
! ===================================================================

   iconstraint = 0

   if( (minFOS < maxPAR_FOS) .and. (minFOS > minPAR_FOS) .and. maxDeflection < par_DEF) then
      iconstraint=1
   endif

   !do i = 1,nmember
   !   if( (memFOS(i) < minFOS) .or. memFOS(i) > maxFOS .or. memDeflection(i) > par_DEF) then
   !      iconstraint = 0
   !      exit
   !   endif
   !enddo

end subroutine updateStructure

! ###################################################################
!
! get geometric parameters
!
! ###################################################################

subroutine crossSectionProperties(maxOuterDimension,matpropTmp)

   use material, only : nmatprop,cstype

   implicit none

   real(kind=8), intent(in)    :: maxOuterDimension
   real(kind=8), intent(inout) :: matpropTmp(nmatprop)
   !
   real(kind=8), parameter :: pi = 4.d0*atan(1.d0)
   !
   real(kind=8) :: r1,r2,area,Iy,Iz,J,maxy,maxz,side

! ===================================================================
! hollow circle cross section
! ===================================================================

   if (cstype=='hollowcircle') then
      
      r2   = maxOuterDimension
      !
      r1   = 0.9d0*r2
      area = pi*(r2*r2-r1*r1)
      J    = 0.5d0*pi*(r2*r2*r2*r2 - r1*r1*r1*r1)
      Iy   = 0.5d0*J
      Iz   = Iy
      maxy = r2
      maxz = r2

! ===================================================================
! solid square cross section
! ===================================================================

   else if (cstype=='solidsquare') then


      side = maxOuterDimension
      area = side*side
      Iy   = (area*area) / 12.0d0
      Iz   = Iy
      J    = 2.0d0*Iy

      maxy = side
      maxz = side

! ===================================================================
! incorrect option. stop
! ===================================================================

   else
      stop 'incorrect cross-section type'
   endif

! ===================================================================
! copy back into matpropTmp
! ===================================================================

   matpropTmp( 4) = area
   matpropTmp( 5) = Iy
   matpropTmp( 6) = Iz
   matpropTmp( 7) = J
   matpropTmp(10) = maxy
   matpropTmp(11) = maxz
 


end subroutine crossSectionProperties
