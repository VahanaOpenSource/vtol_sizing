!=======================================================================
!>       The current 'airfoil' setup is intended to match the format of 
!>    the tables in NASA CR-166309. The airfoil properties are given as
!>     2D tables with angle of attack and Mach number as independent 
!>    variables.  Two different tables are provided for each property, 
!>    namely, for 'low' and 'high' angles of attack.  
!=======================================================================

   module airfoil_tables
      use precision

!=======================================================================
! blade data structure : detailed aerodynamic design
!=======================================================================
      
      type :: BladeDef 
         logical          :: populated = .false.
         integer          :: Nairfoils
         real(kind=rdp)   :: airfoil_id(mxafoil)
         real(kind=rdp)   :: span_end(mxafoil)
      end type

!=======================================================================
! airfoil data structure
!=======================================================================

      TYPE :: AirfoilDef

         CHARACTER (LEN=30)   :: AirfoilName

         INTEGER                             :: AirfoilNumId         !Airfoil number index
         INTEGER                             :: nAlphaCLlow,  nMachCLlow
         INTEGER                             :: nAlphaCLhigh, nMachCLhigh
         INTEGER                             :: nAlphaCDlow,  nMachCDlow
         INTEGER                             :: nAlphaCDhigh, nMachCDhigh
      
         integer                             :: lookup_method        ! 1 = Mach + AoA, 2 = Reynolds + AoA 

!=======================================================================
! Angles at which to switch between low angle and high angle tables
!=======================================================================

         real(kind=rdp)                      :: CLHiPositiveAlpha    !   for positive alphas
         real(kind=rdp)                      :: CLHiNegativeAlpha    !   for negative alphas
         real(kind=rdp)                      :: CDHiPositiveAlpha    !   for positive alphas
         real(kind=rdp)                      :: CDHiNegativeAlpha    !   for negative alphas

!=======================================================================
! Mach, AoA at which data is available
!=======================================================================

         real(kind=rdp),  DIMENSION(11)      :: MachCLlow,  MachCDlow
         real(kind=rdp),  DIMENSION(5)       :: MachCLhigh, MachCDHigh

         real(kind=rdp),  DIMENSION(65)      :: AlphaCLlow,  AlphaCDlow
         real(kind=rdp),  DIMENSION(181)     :: AlphaCLhigh, AlphaCDHigh

!=======================================================================
! Tabulated coefficients: low and high angle ranges
!=======================================================================

         real(kind=rdp),  DIMENSION(65,11)   :: CLlow,  CDlow
         real(kind=rdp),  DIMENSION(181,1)   :: CLHigh, CDhigh

      end type AirfoilDef

!=======================================================================
! This module also contains "mxafoil" airfoils
!=======================================================================

      INTEGER          :: store_index(20)    ! maps airfoil id to storage id
      type(AirfoilDef) :: Tables(6)          ! at most 6 airfoils in a blade 
      type(BladeDef)   :: Blade              ! stores details like spanwise end position

   end module airfoil_tables