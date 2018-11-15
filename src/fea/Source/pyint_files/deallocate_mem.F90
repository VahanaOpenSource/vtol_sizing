!==================================================================================
! routine to perform deallocation of memory
!==================================================================================

subroutine deallocate_variables()

   use currentvalues
   use layout
   use material

   implicit none 

!from preprocessor
   if(allocated(nelem2member)) deallocate(nelem2member)
   if(allocated(elem2member)) deallocate(elem2member)

!from layout module
   if(allocated(nodepos))        deallocate(nodepos)
   if(allocated(conn))           deallocate(conn)
   if(allocated(memberID))       deallocate(memberID)
   if(allocated(extforce_dofid)) deallocate(extforce_dofid)
   if(allocated(extforce))       deallocate(extforce)
   if(allocated(nforce))         deallocate(nforce)
   if(allocated(extdisp_dofid)) deallocate(extdisp_dofid)
   if(allocated(extdisp))       deallocate(extdisp)

!from material properties
   if(allocated(matprop)) deallocate(matprop)

!from initialize currentvalues
      if(allocated(minOptWeight))     deallocate(minOptWeight)
      if(allocated(minOptStressFOS))  deallocate(minOptStressFOS)
      if(allocated(maxOptDeflection)) deallocate(maxOptDeflection)
      if(allocated(idefFlag))         deallocate(idefFlag)
      if(allocated(istressFlag))      deallocate(istressFlag)

!bc related
      if(allocated(rhs_vec)) deallocate(rhs_vec)
      if(allocated(K_red))   deallocate(K_red)
      if(allocated(M_red))   deallocate(M_red)
      if(allocated(rhs_red)) deallocate(rhs_red)

!eigenvalue solution
      if(allocated(K_global)) deallocate(K_global)
      if(allocated(M_global)) deallocate(M_global)

!stress values
      if(allocated(stress))      deallocate(stress)
      if(allocated(stressFOS))   deallocate(stressFOS)
      if(allocated(strain))      deallocate(strain)
      if(allocated(bucklingFOS)) deallocate(bucklingFOS)
      if(allocated(weightElem))  deallocate(weightElem)
      if(allocated(deflection))  deallocate(deflection)

!solution for statics
      if(allocated(xdof))          deallocate(xdof)
      if(allocated(memStress))     deallocate(memStress)
      if(allocated(memFOS))        deallocate(memFOS)
      if(allocated(memDeflection)) deallocate(memDeflection)

!for mode shape extraction 
      if(allocated(modeShapes)) deallocate(modeShapes)
      if(allocated(modeFreqs))  deallocate(modeFreqs)

!boundary condition indices
      if(allocated(row_red)) deallocate(row_red)
      if(allocated(bc_flag)) deallocate(bc_flag)    

end subroutine
