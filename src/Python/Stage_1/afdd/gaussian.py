#=======================================================================
# quadrature weight/location for numerical integration along one dimension
#=======================================================================

import numpy
def gaussian(NG):

   xg             = numpy.zeros(NG)
   wg             = numpy.zeros_like(xg)

#=======================================================================
# generate initial guess values
#=======================================================================

   M = int((NG + 1)/2)		#max # of gauss points to generate (symmetry)
      
   for i in range(0,M):
      Z    = numpy.cos(numpy.pi*(i + 1.0 - 0.25)/(NG + 0.5))
      eps  = 1.0
      while abs(eps) > 1e-10:
         P1 = 1.0          # first   Legendre polynomial
         P2 = 0.0          # second  Legendre polynomial
         for j in range(0,NG):
            P3 = P2
            P2 = P1

#=======================================================================
# n P_n(z) = (2*n-1) z P_(n-1) (z) - (n-1) P_(n-2) (z)
#=======================================================================

            P1 = ((2.0*j + 1.0)*Z*P2 - j*P3)/(j+1)

#=======================================================================
#  d/dx[p_n(z)] = [z P_(n) (z) - P_(n-1) (z)] * n/( z^2-1)
#=======================================================================

         PP = NG*(Z*P1 - P2)/(Z*Z - 1.0) 		#slope of Pn(z)
         Z1 = Z		               

#=======================================================================
#Use Newton-Raphson method to find the "true" zero of Legendre polynomial
#=======================================================================

         Z           = Z1 - P1/PP
         eps         = abs(Z-Z1)
#=======================================================================
# Store converged results in array
#=======================================================================

      xg[i]          = -Z
      xg[NG - 1 - i] = Z
      wg[i]          = 2.0/((1.0 - Z*Z)*PP*PP)	#Gauss-Legendre weight value
      wg[NG - 1 - i] = wg[i]

#=======================================================================
# change limits of integration to be [0,1]
#=======================================================================

   xg             = 0.5*(xg+1.0)
   wg             = 0.5*wg

   return xg, wg