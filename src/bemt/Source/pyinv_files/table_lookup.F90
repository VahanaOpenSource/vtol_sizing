!=======================================================================
!> Given the operating angle of attack and Re, this subroutine finds the
!> lift and drag  coefficients using interpolation/table look up
!=======================================================================

subroutine table_lookup( xnd, alpha, Mach, Reynolds, Cl, Cd)
      
   use airfoil_tables
   implicit none

!=======================================================================
! inputs
!=======================================================================

   real(kind=rdp), intent(in)  :: xnd, alpha, Mach, Reynolds

!=======================================================================
! outputs
!=======================================================================

   real(kind=rdp), intent(out) :: Cl, Cd
      
!=======================================================================
! local variables
!=======================================================================

   logical                     :: foundAirfoil
   integer                     :: iAirfoil
   real(kind=rdp)              :: xend
   real(kind=rdp)              :: Mach_or_Reynolds

!=======================================================================
! Begin executable code
!=======================================================================

   iAirfoil           = 1
   foundAirfoil       = .false.

!=======================================================================
! keep searching while airfoils are not identified
! find where this airfoil ends, see if the point is inboard of that limit
! if yes, we know what airfoil it is; if not, keep searching 
!=======================================================================

   do while(.not. foundAirfoil)
      xend            = Blade % span_end(iAirfoil)
      if(xnd .le. xend) then
         foundAirfoil = .true.
      else
         iAirfoil     = iAirfoil + 1
      end if
   end do

!=======================================================================
!map airfoil numerical identifier to storage location in memory
!=======================================================================

   iAirfoil           = Blade % airfoil_id(iAirfoil)

!=======================================================================
! Choose Reynolds or Mach number for look-up tables
!=======================================================================
  
   if (Tables(iAirfoil) % lookup_method == 1) then
      Mach_or_Reynolds     = Mach
   else 
      Mach_or_Reynolds     = Reynolds
   end if 

!=======================================================================
! Find Cl, Cd, Cm for the identified airfoil
!=======================================================================

   call interpolate_tables(alpha, Mach_or_Reynolds, Tables(iAirfoil), Cl, Cd)

!=======================================================================
! End of operations
!=======================================================================
      
   return
   end subroutine table_lookup
