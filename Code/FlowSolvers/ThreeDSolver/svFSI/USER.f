!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!      
!--------------------------------------------------------------------
!      
!     This subroutine initializes the parameters, you need to read the
!     comments and specify them carefuly
!      
!--------------------------------------------------------------------
!     This is an example of Windkessel model (RCR) coupled to a Neumann
!     BC

      SUBROUTINE cplBC_INI(nFaces, nTimeSteps, qConv, pConv, face)
      
      IMPLICIT NONE
      
      INCLUDE "cplBC.h"
      
      INTEGER, INTENT(IN) :: nFaces
      INTEGER, INTENT(OUT) :: nTimeSteps
      REAL(KIND=8), INTENT(OUT) :: pConv, qConv
      TYPE(cplFaceType), INTENT(OUT) :: face(nFaces)

!--------------------------------------------------------------------
!     Block for your inputs
!     For instance if pressure in 3D solver is in cgs and here mmHg
!     pConv=1334=1.334D3, also if flowrate in 3D solver is mm^3/s and 
!     here is mL/s qConv=1000=1D3. In the case both solver are using 
!     same unites you can set these two conversion coefficients to 1D0
      pConv = 1.334D3
      qConv = 1D0

!     Number of time step between N and N+alpha
      nTimeSteps = 10

!     The name, bGrp and x pointer are set here
      INCLUDE "faces.f"

      END SUBROUTINE cplBC_INI

!####################################################################
!     Here you should find the f_i=dx_i/dt, based on the following 
!     parameters:
!     current x_i:         x(i)
!     Flowrates of faces:  Q(i)
!     Pressures of faces:  P(i)

      SUBROUTINE cplBC_FINDF(nFaces, nX, x, f, Q, P, offset, nXprint, 
     2   Xprint, t)
      
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nFaces, nX, nXprint
      REAL(KIND=8), INTENT(IN) :: Q(nFaces), P(nFaces), t
      REAL(KIND=8), INTENT(OUT) :: f(nX), offset(nFaces)
      REAL(KIND=8), INTENT(INOUT) :: x(nX), Xprint(nXprint)
      
!     These are the dumy variables 
      INCLUDE "parameters.f"

!     The main body of equations
      f(1) = (Q(1) - x(1)/Rd(1))/C(1)
      offset(1) = Q(1)*Rp(1)
      
      Xprint(1) = t
      Xprint(2) = Q(1)
      
      RETURN
      END SUBROUTINE cplBC_FINDF
