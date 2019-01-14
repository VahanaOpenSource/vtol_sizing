from numpy import zeros,sqrt,linspace,unique
from scipy import linalg

#======================================================================
# python function to calculate element mass and stiffness matrices for
# vertical bending given beam mass per unit span "m" (mass/length) and 
# flexural stiffness "EI" (Force length^2) and element size "le" (length)
#======================================================================

def element_matrices(m,EI,le):

	M 	 	= zeros((4,4))
	K 	 	= zeros((4,4))

#======================================================================
# mass matrix entries for element with uniform mass along its length
#======================================================================

	M[0,0] 	= m*13.0*le/35.0
	M[0,1] 	= m*11.0*le*le/210.0
	M[0,2]  = m*9.0*le/70.0
	M[0,3] 	= m*(-13.0)*le*le/420.0

	M[1,0] 	= M[0,1]
	M[1,1] 	= m*le*le*le/105.0
	M[1,2] 	=-M[0,3]
	M[1,3]  = m*(-le*le*le)/140.0

	M[2,0]  = M[0,2]
	M[2,1]  = M[1,2]
	M[2,2]  = M[0,0]
	M[2,3]  =-M[0,1]

	M[3,0]  = M[0,3]
	M[3,1]  = M[1,3]
	M[3,2]  = M[2,3]
	M[3,3]  = M[1,1]

#======================================================================
# stiffness matrix entries
#======================================================================

	linv 	= 1.0/le 
	linsq  	= linv*linv
	lincu  	= linsq*linv
	K[0,0] 	= EI*12.0*lincu 
	K[0,1]  = EI* 6.0*linsq
	K[0,2] 	=-K[0,0]
	K[0,3] 	= K[0,1]

	K[1,0] 	= K[0,1]
	K[1,1] 	= EI*4.0*linv
	K[1,2]  =-K[1,0]
	K[1,3] 	= 0.5*K[1,1]

	K[2,0] 	= K[0,2]
	K[2,1] 	= K[1,2]
	K[2,2] 	= K[0,0]
	K[2,3] 	= K[2,1]

	K[3,0] 	= K[0,3]
	K[3,1] 	= K[1,3]
	K[3,2] 	= K[2,3]
	K[3,3] 	= K[1,1]

	return M,K

#======================================================================
# python function to assemble matrices for cantilever beam vibration analysis
#======================================================================

def assemble(A_elem):

	ne 		= len(A_elem)
	ndof 	= 2*ne+2
	nbc 	= 2
	nred 	= ndof-nbc
	Ared 	= zeros((nred,nred)) 		# for a cantilever beam, 2n dof
	for k in range(ne):
		for i in range(4):
			iglobal 	= i + 2*k - nbc
			for j in range(4):
				jglobal = j + 2*k - nbc

				if(iglobal >= 0 and jglobal >= 0):
					Ared[iglobal,jglobal] = Ared[iglobal,jglobal] + A_elem[k][i,j]

	return Ared 

#======================================================================
# python function to generate elements so that lumped masses are placed
# at the right ends of certain elements; also tag those elements so we
# know which elements we have to augment the mass matrix for
#======================================================================

def generate_elements(xposn,L,reqd_node=0.0):

	nodes  	= list(linspace(0,L,6))
	for x in xposn:
 		nodes.append(x)
	nodes.append(reqd_node)
	nodes 	= sorted(unique(nodes))
	ne 		= len(nodes)-1
	le 		= [] 
	ie 		= [] 						
	for i in range(ne): 			 	# loop over elements
		le.append(nodes[i+1]-nodes[i])
		if(nodes[i+1] in xposn):
			ie.append(i) 				# tag element with tip mass at right end
	# print(nodes)
	# print(le)
	# print(ie)
	# quit()
	return le,ie

#======================================================================
# python function to build M,K matrices and return natural frequencies
# for a nonrotating cantilever beam 
# inputs are element lengths "le", element mass/span "m" and element 
# flexural stiffnesses "EI" as 3 lists
# "elem_id" = element indices which have lumped masses "Masses" at the right 
# ends of those elements
#======================================================================

def beam_freq_FEM(m_list, EI_list, le_list, elem_id, Masses):
	Me      = []
	Ke 		= []

#======================================================================
# first build the regular matrix
#======================================================================
	# print('have to add tip mass for this element:',elem_id)
	ie 		= -1
	im 		= 0 						# counter for lumped masses
	# print(m_list,EI_list,le_list)
	for m,EI,le in zip(m_list,EI_list,le_list):
		ie  	= ie + 1
		M,K 	= element_matrices(m,EI,le)

#======================================================================
# append mass matrix with lumped mass at right end: modify 2,2 entry
#======================================================================

		if(ie in elem_id):
			M[2,2] 	= M[2,2] + Masses[im]
			im 		= im + 1
		Me.append(M)
		Ke.append(K)

#======================================================================
# get assembled and reduced matrices for eigenvalue solution
#======================================================================

	wn 		= nat_freq(Me,Ke)

	return wn

#======================================================================
# get natural frequencies from element matrices
#======================================================================

def nat_freq(Me, Ke):
	Mred  	= assemble(Me)
	Kred 	= assemble(Ke)
	a,b  	= linalg.eig(Kred, Mred)
	a 		= sqrt(abs(a))[-1:0:-1] 	 # sorted in ascending order
	return a 