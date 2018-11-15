!=======================================================================
!> This subroutine computes the lift/drag coefficient of an airfoil 
!> using table look-up. 
! though input says Mach number, it can also accommodate Reynolds #
!=======================================================================

      subroutine interpolate_tables(alpha, Machnum, Airfoil, Cl, Cd) 
      use airfoil_tables
      implicit none

!=======================================================================
!                           inputs
!=======================================================================

      real(kind=rdp)    , intent(in) :: alpha, Machnum
      type(AirfoilDef)  , intent(in) :: Airfoil

!=======================================================================
!                          outputs
!=======================================================================

      real(kind=rdp)    , intent(out):: Cl, Cd 

!=======================================================================
!                       local variables
!=======================================================================

      integer        :: index_ofst, index_left, index_right, Mach_ofst, &
                         Mach_left, Mach_right, Nmax, Nalpha
      real(kind=rdp) :: delta_alfa, x1, x2, y1, y2, x3, slope,          &
                        delta_Mach, y11, y12, y21, y22, mach,           &
                        Mmax, Mmin, z11, z12, z21, z22, z1, z2, Cla, Cda

      real(kind=rdp) :: amin, amax, temp 

!=======================================================================
!                    Begin executable code
!=======================================================================

!=======================================================================
!          Apply max cutoff for mach number based on tables
!=======================================================================

      Nmax   = Airfoil % nMachCLlow 
      Nalpha = Airfoil % nalphaClhigh

      Mmax   = Airfoil % MachCLlow(Nmax) - 0.001d0    ! add buffer for num. safety
      Mmin   = Airfoil % MachCLlow(1)    + 0.001d0    ! add buffer for num. safety

!      amin   = Airfoil % alphaCLlow(1)          ! min AoA 
!      amax   = Airfoil % alphaCLhigh(Nalpha)    ! max AoA

      Mach   = min(Mmax, Machnum)
      Mach   = max(Mmin, Mach)

!=======================================================================
! Compute lift coefficient
!=======================================================================

      if (alpha.gt.Airfoil % CLHiPositiveAlpha .or.         &
          alpha.lt.Airfoil % CLHiNegativeAlpha) then
 
!=======================================================================
! high alpha table :no mach number variation + regular intervals  
!=======================================================================

         delta_alfa = Airfoil % AlphaCLhigh(2) - Airfoil % AlphaCLhigh(1)
         index_ofst = floor(-Airfoil % AlphaCLhigh(1)/delta_alfa)

         index_left  = floor(alpha/delta_alfa) + 1 + index_ofst
         index_right = index_left + 1
         x1          = Airfoil % AlphaCLhigh(index_left)
         x2          = Airfoil % AlphaCLhigh(index_right)
         y1          = Airfoil % CLHigh(index_left,1)
         y2          = Airfoil % CLHigh(index_right,1)
         x3          = alpha
         slope       = (y2-y1)/(x2-x1)
         Cl          = y1 + slope*(x3-x1)

!find drag also
         y1          = Airfoil % CDHigh(index_left,1)
         y2          = Airfoil % CDHigh(index_right,1)
         slope       = (y2-y1)/(x2-x1)
         Cd          = y1 + slope*(x3-x1)

      else

!=======================================================================
!    we are using the low alpha table for lift and pitching moment
!=======================================================================

         delta_alfa  = Airfoil % AlphaCLlow(2) -Airfoil % AlphaCLlow(1)

         delta_alfa  = 1.d0/delta_alfa       
         index_ofst  = floor(-Airfoil % AlphaCLlow(1) * delta_alfa)
         index_left  = floor(alpha*delta_alfa) + 1 + index_ofst
         index_right = index_left + 1

!=======================================================================
!found aoa indices for bracketing, now get mach indices for bracketing
!=======================================================================

         delta_Mach  = Airfoil % MachCLlow(2) -                      &
                       Airfoil % MachCLlow(1)

         Mach_ofst   = 1 - floor(Airfoil % MachCLlow(1)/delta_Mach)

         Mach_left   = floor(Mach/delta_Mach)+ Mach_ofst
         Mach_right  = Mach_left + 1

!=======================================================================
!                  get Cl(a,M) at 4 locations
!                      a1,M1    a2, M1
!                      a1,M2    a2, M2
!=======================================================================

         y11         = Airfoil % CLlow(index_left, Mach_left)
         y12         = Airfoil % CLlow(index_left, Mach_right)
         y21         = Airfoil % CLlow(index_right,Mach_left)
         y22         = Airfoil % CLlow(index_right,Mach_right)

!=======================================================================
! find drag values from table
!=======================================================================

         z11         = Airfoil % CDlow(index_left, Mach_left)
         z12         = Airfoil % CDlow(index_left, Mach_right)
         z21         = Airfoil % CDlow(index_right,Mach_left)
         z22         = Airfoil % CDlow(index_right,Mach_right)

!=======================================================================
!              Interpolate along Mach # to get
!                    a1,M      a2, M
!=======================================================================

         x3          = Mach
         x1          = Airfoil % MachCLlow(Mach_left)       ! low  Re
         x2          = Airfoil % MachCLlow(Mach_right)      ! high Re

         temp        = (x3-x1)/(x2-x1)
         y1          = y11 + (y12-y11)*temp       ! lift @ low Re
         y2          = y21 + (y22-y21)*temp      ! lift @ hi  Re

         z1          = z11 + (z12-z11)*temp      ! drag @ low Re 
         z2          = z21 + (z22-z21)*temp      ! drag @ hi  Re 

!=======================================================================
!              now interpolate along alpha
!=======================================================================

         x1          = Airfoil % AlphaCLlow(index_left)     ! low alpha 
         x3          = alpha                                ! aoa 

         delta_alfa  = delta_alfa * (x3-x1)
         Cl          = y1 + (y2-y1)*delta_alfa!*(x3-x1)
         Cd          = z1 + (z2-z1)*delta_alfa!*(x3-x1)

      end if


!=======================================================================
! End of operations
!=======================================================================
   return
   end subroutine interpolate_tables