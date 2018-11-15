!=======================================================================
!> this subroutine reads airfoil tables and places it in a module for 
!> performing table look-up for lift and drag
!=======================================================================

subroutine setup_tables()
   use airfoil_tables
   implicit none

!=======================================================================
! local variables
!=======================================================================

   integer           :: i, icount, i_data, N_all, k 
   logical           :: exist , found_airfoil(mxafoil)
   type(AirfoilDef)  :: temp 

!=======================================================================
! Begin executable code
!=======================================================================

!=======================================================================
! check if already setup: return and do nothing
!=======================================================================

   if(Blade % populated) return 

!=======================================================================
! otherwise, read the input files
!=======================================================================

   open(20,file='BladeData',status='old',err=22)

!=======================================================================
! # of airfoils in this rotor design  
!=======================================================================

   read(20,*) Blade % Nairfoils

!=======================================================================
! read spanwise end position of each of these airfoils, along with airfoil
! "type" identification 
!=======================================================================

   do i = 1, Blade % Nairfoils  
      read(20,*,err=22) Blade % airfoil_id(i)
      read(20,*,err=22) Blade % span_end(i)
      found_airfoil(i)     = .false. 
   end do 
   close(20)

!=======================================================================
! Open airfoil details
!=======================================================================

   open(20,file='AirfoilData',status='old',err=23)

!=======================================================================
! the file is open; now read total # of airfoils
!=======================================================================

   read(20,*) N_all
   icount               = 0

!=======================================================================
! loop over and read each airfoil in the database 
!=======================================================================

   do i = 1, N_all         
      call read_AirfoilData(temp)

      i_data            = temp % AirfoilNumId

!=======================================================================
! see whether this airfoil matches any of the airfoils in the blade
! if yes, store the airfoil in memory, and remember the storage index 
! for the corresponding airfoil
!=======================================================================

      do k = 1, Blade % Nairfoils               
         if (i_data == Blade % airfoil_id(k)) then 
            icount                  = icount + 1 
            Tables(icount)          = temp
            Blade % airfoil_id(k)   = icount 
            found_airfoil(k)        = .true.
            write(*,*) 'remembering airfoil ',temp % AirfoilName, ' in storage location', icount 
         end if 
      end do 
   end do 

   do i = 1, Blade % Nairfoils
      if(.not. found_airfoil(i)) then 
         write(*,*) 'critical error: did not find, in database, airfoil # ',i 
         stop 'program halting'
      end if 
   end do 

!=======================================================================
! trip flag so that we dont read input files again on repeated bemt calls
!=======================================================================

   Blade % populated = .true.
!=======================================================================
! End of operations
!=======================================================================
      
   return
22 stop 'CRITICAL PROGRAM ERROR: BladeData file incomplete'
23 stop 'CRITICAL PROGRAM ERROR: HALTING'
   end subroutine setup_tables