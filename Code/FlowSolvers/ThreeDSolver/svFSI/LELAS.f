!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!--------------------------------------------------------------------
!      
!     This is for solving linear elasticity
!      
!--------------------------------------------------------------------

      PURE SUBROUTINE LELAS3D (eNoN, w, N, Nx, al, dl, lR, lK)

      USE COMMOD

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(nsd,eNoN), 
     2   al(tDof,eNoN), dl(tDof,eNoN)
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,eNoN), lK(dof*dof,eNoN,eNoN)

      INTEGER a, b, i, j, k
      REAL(KIND=8) NxdNx, rho, elM, nu, lambda, mu, divD, T1, amd, wl, 
     2   lDm, ed(6), ud(nsd), f(nsd)
      
      rho  = eq(cEq)%dmn(cDmn)%prop(solid_density)
      elM  = eq(cEq)%dmn(cDmn)%prop(elasticity_modulus)
      nu   = eq(cEq)%dmn(cDmn)%prop(poisson_ratio)
      f(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      f(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      f(3) = eq(cEq)%dmn(cDmn)%prop(f_z) 
      i    = eq(cEq)%s
      j    = i + 1
      k    = j + 1
      
      lambda = elM*nu/(1D0 + nu)/(1D0 - 2D0*nu)
      mu     = elM/2D0/(1D0 + nu)
      lDm    = lambda/mu
      T1     = eq(cEq)%af*eq(cEq)%beta*dt*dt
      amd    = eq(cEq)%am/T1*rho
      wl     = w*T1*mu

      ed = 0D0
      ud = -f
      DO a=1, eNoN
         ud(1) = ud(1) + N(a)*al(i,a)
         ud(2) = ud(2) + N(a)*al(j,a)
         ud(3) = ud(3) + N(a)*al(k,a)
               
         ed(1) = ed(1) + Nx(1,a)*dl(i,a)
         ed(2) = ed(2) + Nx(2,a)*dl(j,a)
         ed(3) = ed(3) + Nx(3,a)*dl(k,a)
         ed(4) = ed(4) + Nx(3,a)*dl(j,a) + Nx(2,a)*dl(k,a)
         ed(5) = ed(5) + Nx(3,a)*dl(i,a) + Nx(1,a)*dl(k,a)
         ed(6) = ed(6) + Nx(2,a)*dl(i,a) + Nx(1,a)*dl(j,a)
      END DO
      divD = lambda*(ed(1) + ed(2) + ed(3))
 
      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*(rho*N(a)*ud(1) + mu*Nx(2,a)*ed(6) 
     2      + mu*Nx(3,a)*ed(5) + Nx(1,a)*(2D0*mu*ed(1) + divD))

         lR(2,a) = lR(2,a) + w*(rho*N(a)*ud(2) + mu*Nx(1,a)*ed(6) 
     2      + mu*Nx(3,a)*ed(4) + Nx(2,a)*(2D0*mu*ed(2) + divD))
         
         lR(3,a) = lR(3,a) + w*(rho*N(a)*ud(3) + mu*Nx(1,a)*ed(5) 
     2      + mu*Nx(2,a)*ed(4) + Nx(3,a)*(2D0*mu*ed(3) + divD))

         DO b=1, eNoN
            NxdNx = Nx(1,a)*Nx(1,b) + Nx(2,a)*Nx(2,b) + Nx(3,a)*Nx(3,b)

            T1 = amd*N(a)*N(b)/mu + NxdNx
            
            lK(1,a,b) = lK(1,a,b) + wl*(T1 
     2         + (1D0 + lDm)*Nx(1,a)*Nx(1,b))

            lK(2,a,b) = lK(2,a,b) + wl*(lDm*Nx(1,a)*Nx(2,b) 
     2         + Nx(2,a)*Nx(1,b))

            lK(3,a,b) = lK(3,a,b) + wl*(lDm*Nx(1,a)*Nx(3,b) 
     2         + Nx(3,a)*Nx(1,b))

            lK(dof+1,a,b) = lK(dof+1,a,b) + wl*(lDm*Nx(2,a)*Nx(1,b) 
     2         + Nx(1,a)*Nx(2,b))

            lK(dof+2,a,b) = lK(dof+2,a,b) + wl*(T1
     2         + (1D0 + lDm)*Nx(2,a)*Nx(2,b))
            
            lK(dof+3,a,b) = lK(dof+3,a,b) + wl*(lDm*Nx(2,a)*Nx(3,b) 
     2         + Nx(3,a)*Nx(2,b))

            lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + wl*(lDm*Nx(3,a)*Nx(1,b)
     2         + Nx(1,a)*Nx(3,b))

            lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + wl*(lDm*Nx(3,a)*Nx(2,b) 
     2         + Nx(2,a)*Nx(3,b))
            
            lK(2*dof+3,a,b) = lK(2*dof+3,a,b) + wl*(T1 
     2         + (1D0 + lDm)*Nx(3,a)*Nx(3,b))
         END DO
      END DO

      RETURN
      END SUBROUTINE LELAS3D

!====================================================================
!     This is for boundary elasticity equation
      PURE SUBROUTINE BLELAS (eNoN, w, N, h, nV, lR)

      USE COMMOD

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), h, nV(nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,eNoN)

      INTEGER a, i
      REAL(KIND=8) hc(nsd)

      hc = h*nV
      DO a=1,eNoN
         DO i=1, nsd
            lR(i,a) = lR(i,a) - w*N(a)*hc(i)
         END DO
      END DO

      RETURN
      END SUBROUTINE BLELAS

!####################################################################
      
      PURE SUBROUTINE LELAS2D (eNoN, w, N, Nx, al, dl, lR, lK)

      USE COMMOD

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(nsd,eNoN), 
     2   al(tDof,eNoN), dl(tDof,eNoN)
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,eNoN), lK(dof*dof,eNoN,eNoN)

      INTEGER a, b, i, j
      REAL(KIND=8) NxdNx, rho, elM, nu, lambda, mu, divD, T1, amd, wl, 
     2   lDm, ed(3), ud(nsd), f(nsd)
      
      rho  = eq(cEq)%dmn(cDmn)%prop(solid_density)
      elM  = eq(cEq)%dmn(cDmn)%prop(elasticity_modulus)
      nu   = eq(cEq)%dmn(cDmn)%prop(poisson_ratio)
      f(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      f(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      i    = eq(cEq)%s
      j    = i + 1
      
      lambda = elM*nu/(1D0 + nu)/(1D0 - 2D0*nu)
      mu     = elM/2D0/(1D0 + nu)
      lDm    = lambda/mu
      T1     = eq(cEq)%af*eq(cEq)%beta*dt*dt
      amd    = eq(cEq)%am/T1*rho
      wl     = w*T1*mu

      ed = 0D0
      ud = -f
      DO a=1, eNoN
         ud(1) = ud(1) + N(a)*al(i,a)
         ud(2) = ud(2) + N(a)*al(j,a)
         
         ed(1) = ed(1) + Nx(1,a)*dl(i,a)
         ed(2) = ed(2) + Nx(2,a)*dl(j,a)
         ed(3) = ed(3) + Nx(2,a)*dl(i,a) + Nx(1,a)*dl(j,a)
      END DO
      divD = lambda*(ed(1) + ed(2))

      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*(rho*N(a)*ud(1) + mu*Nx(2,a)*ed(3) 
     2                     + Nx(1,a)*(2D0*mu*ed(1) + divD))

         lR(2,a) = lR(2,a) + w*(rho*N(a)*ud(2) + mu*Nx(1,a)*ed(3) 
     2                     + Nx(2,a)*(2D0*mu*ed(2) + divD))
         
         DO b=1, eNoN
            NxdNx = Nx(1,a)*Nx(1,b) + Nx(2,a)*Nx(2,b)
            
            T1 = amd*N(a)*N(b)/mu + NxdNx
            
            lK(1,a,b) = lK(1,a,b) + wl*(T1 
     2         + (1D0 + lDm)*Nx(1,a)*Nx(1,b))

            lK(2,a,b) = lK(2,a,b) + wl*(lDm*Nx(1,a)*Nx(2,b) 
     2         + Nx(2,a)*Nx(1,b))

            lK(dof+1,a,b) = lK(dof+1,a,b) + wl*(lDm*Nx(2,a)*Nx(1,b) 
     2         + Nx(1,a)*Nx(2,b))

            lK(dof+2,a,b) = lK(dof+2,a,b) + wl*(T1
     2         + (1D0 + lDm)*Nx(2,a)*Nx(2,b))
         END DO
      END DO

      RETURN
      END SUBROUTINE LELAS2D

