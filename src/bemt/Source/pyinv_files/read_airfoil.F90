!======================================================================
!> Reads data tables for a rotor airfoil sections
!======================================================================

      subroutine read_AirfoilData(ThisAirfoil) 
      use airfoil_tables
      implicit none 

!======================================================================
!                             Outputs
!======================================================================

      type(AirfoilDef), intent(out) :: ThisAirfoil
      
!======================================================================
!                       Local Variables
!======================================================================

      character*80            :: EchoLine

      integer                 :: i, jmach, ialpha

      real(kind=rdp)          :: alphamax, alphamin 

!======================================================================
!                     Begin executable code
!======================================================================

!skip next five lines
      do i=1,5
         read (20, *, end=240, err=240) EchoLine
      end do
      
      read (20, 9040, end=260, err=260) ThisAirfoil % AirfoilNumId
      read (20, 9060, end=280, err=280) ThisAirfoil % AirfoilName

!======================================================================
!          Lift coefficient - low angle table
!======================================================================

!skip next five lines
      do i=1,5
         read (20, *, end=320, err=320) EchoLine
      end do

!======================================================================
! Read number of Mach/Re number data points, # of angle of attack values
!======================================================================
  
      read (20, 9080, end=340, err=340) ThisAirfoil % nMachCLlow, &
                                        ThisAirfoil % nAlphaCLlow

!======================================================================
! Default min/max limits for angle of attack
!======================================================================
  
      alphaMax                    = -91.d0
      alphaMin                    =  91.d0
      ThisAirfoil % lookup_method = 1

!======================================================================
! Loop over Mach/Reynolds numbers
!======================================================================

      do jMach = 1,ThisAirfoil % nMachCLlow 

!======================================================================
! Loop over angles of attack
!======================================================================

         do iAlpha = 1,ThisAirfoil % nAlphaCLlow

!            write(*,*) 'in lift',ialpha,jmach
            read (20, 9100, end=360, err=360)                           &
                ThisAirfoil % MachCLlow(jMach),                         &
                ThisAirfoil % AlphaCLlow(iAlpha),                       &
                ThisAirfoil % CLlow(iAlpha,jMach) 

            if (ThisAirfoil % AlphaCLlow(iAlpha).lt.alphaMin)       &
                   alphaMin = ThisAirfoil % AlphaCLlow(iAlpha)

            if (ThisAirfoil % AlphaCLlow(iAlpha).gt.alphaMax)       &
                   alphaMax = ThisAirfoil % AlphaCLlow(iAlpha)

!======================================================================
! If the "Mach" number is more than 2, interpret it as Reynolds #
! Once this trigger is encountered once, its used for all airfoil 
! properties, i.e. all "mach" number entries are interpreted as Re
!======================================================================

            if (ThisAirfoil % MachCLlow(jMach) .gt. two) then
              ThisAirfoil % lookup_method = 2
            end if 
         end do
      end do

      ThisAirfoil % CLHiPositiveAlpha = alphaMax
      ThisAirfoil % CLHiNegativeAlpha = alphaMin
  
!======================================================================
! Lift coefficient - high angle table
!======================================================================

!skip next five lines
      do i=1,5
         read (20, *, end=400, err=400) EchoLine
      end do
  
      read (20, 9080, end=440, err=440) ThisAirfoil % nMachCLhigh, &
                                        ThisAirfoil % nAlphaCLhigh
  
      do jMach = 1,ThisAirfoil % nMachCLhigh
         do iAlpha = 1,ThisAirfoil % nAlphaCLhigh

            read (20, 9100, end=480, err=480)               &
                ThisAirfoil % MachCLhigh(jMach),            &
                ThisAirfoil % AlphaCLhigh(iAlpha),          &
                ThisAirfoil % CLhigh(iAlpha,jMach) 

         end do
      end do

!======================================================================
!     Drag coefficient - low angle table
!======================================================================

!skip next five lines

      do i=1,5
         read (20, *, end=520, err=520) EchoLine
      end do
  
      read (20, 9080, end=560, err=560) ThisAirfoil % nMachCDlow,       &
                                        ThisAirfoil % nAlphaCDlow
  
      alphaMax = -9999.d0
      alphaMin =  9999.d0

      do jMach = 1,ThisAirfoil % nMachCDlow 
         do iAlpha = 1,ThisAirfoil % nAlphaCDlow

!            write(*,*) 'in drag'
            read (20, 9100, end=360, err=360)                           &
                ThisAirfoil % MachCDlow(jMach),                         &
                ThisAirfoil % AlphaCDlow(iAlpha),                       &
                ThisAirfoil % CDlow(iAlpha,jMach) 

                if (ThisAirfoil % AlphaCDlow(iAlpha).lt.alphaMin)       &
                       alphaMin = ThisAirfoil % AlphaCDlow(iAlpha)
                if (ThisAirfoil % AlphaCDlow(iAlpha).gt.alphaMax)       &
                       alphaMax = ThisAirfoil % AlphaCDlow(iAlpha)

         end do
      end do

      ThisAirfoil % CDHiPositiveAlpha = alphaMax
      ThisAirfoil % CDHiNegativeAlpha = alphaMin
  
!======================================================================
!             drag coefficient - high angle table
!======================================================================

!skip next five lines
      do i=1,5
         read (20, *, end=640, err=640) EchoLine
      end do
  
      read (20, 9080, end=680, err=680) ThisAirfoil % nMachCDhigh, &
                                        ThisAirfoil % nAlphaCDhigh
  
      do jMach = 1,ThisAirfoil % nMachCDhigh
         do iAlpha = 1,ThisAirfoil % nAlphaCDhigh

            read (20, 9100, end=720, err=720)         &
                ThisAirfoil % MachCDhigh(jMach),      &
                ThisAirfoil % AlphaCDhigh(iAlpha),    &
                ThisAirfoil % CDhigh(iAlpha,jMach) 

         end do
      end do

!======================================================================
!                       End of Operations
!======================================================================

#ifdef display_off
#else
!      write(*,'(8x,A,A)') 'I finished reading an airfoil called ... ',&
!                            ThisAirfoil % AirfoilName(1:20)
#endif
      return

!======================================================================
!                       Error Messages
!======================================================================

 9100 format (2x, f12.3, 2x, f8.3, 2x, d13.6)
 9080 format (2x, i4, 2x, i4)
 9060 format (1x, A30)
 9040 format (1x, i5)
  240 stop 'RAFDT error 1'
  260 stop 'RAFDT error 2'
  280 stop 'RAFDT error 3'
  320 stop 'RAFDT error 4'
  340 stop 'RAFDT error 5'
  360 stop 'RAFDT error 6'
  400 stop 'RAFDT error 7'
  440 stop 'RAFDT error 8'
  480 stop 'RAFDT error 9'
  520 stop 'RAFDT error 10'
  560 stop 'RAFDT error 11'
  600 stop 'RAFDT error 12'
  640 stop 'RAFDT error 13'
  680 stop 'RAFDT error 14'
  720 stop 'RAFDT error 15'
  760 stop 'RAFDT error 16'
  800 stop 'error reading pitching moment tables (low angles)'
  840 stop 'RAFDT error 18'
  880 stop 'RAFDT error 19'
  920 stop 'RAFDT error 20'
  960 stop 'RAFDT error 21'
      
!======================================================================

      end subroutine read_AirfoilData
