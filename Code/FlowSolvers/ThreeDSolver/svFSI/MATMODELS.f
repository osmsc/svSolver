!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!--------------------------------------------------------------------
!      
!     This routines is for material models for structural dynamics.
!      
!--------------------------------------------------------------------

!####################################################################
!     St.Venant-Kirchhoff model
      SUBROUTINE PK2_stVenantKirchhoff (elM, nu, F, S, CC)
      
      USE COMMOD
      USE ALLFUN
      
      IMPLICIT NONE
      
      REAL(KIND=8), INTENT(IN) :: elM, nu, F(nsd,nsd)
      REAL(KIND=8), INTENT(OUT) :: S(nsd,nsd), CC(nsd,nsd,nsd,nsd)
      
      REAL(KIND=8) :: lam, mu, kap, trE, Im(nsd,nsd), C(nsd,nsd),
     2   E(nsd,nsd)
      
      lam = elM*nu/(1D0 + nu)/(1D0 - 2D0*nu)
      mu  = elM/2D0/(1D0 + nu)
      kap = lam + 2D0/3D0*mu
      Im  = MAT_ID(nsd)
      C   = MATMUL(TRANSPOSE(F), F)
      E   = 5D-1 * (C - Im)
      trE = MAT_TRACE(E, nsd)
      
      S   = 2D0*mu*E + lam*trE*Im
      CC  = lam*TEN_DYADPROD(Im, Im, nsd) + 2D0*mu*TEN_IDs(nsd)
      
      RETURN
      END SUBROUTINE PK2_stVenantKirchhoff
!--------------------------------------------------------------------
!     Modified St.Venant-Kirchhoff (Simo85 penalty) model
      SUBROUTINE PK2_stVenantKirchhoffSimo85 (elM, nu, F, S, CC)
      
      USE COMMOD
      USE ALLFUN
      
      IMPLICIT NONE
      
      REAL(KIND=8), INTENT(IN) :: elM, nu, F(nsd,nsd)
      REAL(KIND=8), INTENT(OUT) :: S(nsd,nsd), CC(nsd,nsd,nsd,nsd)
      
      REAL(KIND=8) :: lam, mu, kap, lnJ, Im(nsd,nsd), C(nsd,nsd),
     2   Ci(nsd,nsd)
      
      lam = elM*nu/(1D0 + nu)/(1D0 - 2D0*nu)
      mu  = elM/2D0/(1D0 + nu)
      kap = lam + 2D0/3D0*mu
      Im  = MAT_ID(nsd)
      lnJ = LOG(MAT_DET(F, nsd))
      C   = MATMUL(TRANSPOSE(F), F)
      Ci  = MAT_INV(C, nsd)
      
      S   = kap*lnJ*Ci + mu*(C-Im)
      CC  = kap * ( -2D0*lnJ*TEN_SYMMPROD(Ci, Ci, nsd) + 
     2   TEN_DYADPROD(Ci, Ci, nsd) ) + 2D0*mu*TEN_IDs(nsd)
      
      RETURN
      END SUBROUTINE PK2_stVenantKirchhoffSimo85
!--------------------------------------------------------------------
!     Neo-Hookean model
      SUBROUTINE PK2_neoHookean (elM, nu, F, S, CC)
      
      USE COMMOD
      USE ALLFUN
      
      IMPLICIT NONE
      
      REAL(KIND=8), INTENT(IN) :: elM, nu, F(nsd,nsd)
      REAL(KIND=8), INTENT(OUT) :: S(nsd,nsd), CC(nsd,nsd,nsd,nsd)
      
      REAL(KIND=8) :: lam, mu, kap, J, J2, J23, tC3, tC9, Im(nsd,nsd),
     2   C(nsd,nsd), Ci(nsd,nsd)
      REAL(KIND=8), DIMENSION(nsd,nsd,nsd,nsd) :: T1, T2, T3
      
      lam = elM*nu/(1D0 + nu)/(1D0 - 2D0*nu)
      mu  = elM/2D0/(1D0 + nu)
      kap = lam + 2D0/3D0*mu
      Im  = MAT_ID(nsd)
      J   = MAT_DET(F, nsd)
      J2  = J**2
      J23 = J**(-2D0/3D0)
      C   = MATMUL(TRANSPOSE(F), F)
      Ci  = MAT_INV(C, nsd)
      tC3 = MAT_TRACE(C, nsd) / 3D0
      tC9 = tC3/3D0
      
      S   = mu*J23*(Im - tC3*Ci) + 5D-1*kap*(J2-1D0)*Ci
      T1  = (2D0*mu*J23*tC9 + kap*J2) * TEN_DYADPROD(Ci, Ci, nsd)
      T2  = (2D0*mu*J23*tC3 - kap*(J2-1D0)) * TEN_SYMMPROD(Ci, Ci, nsd)
      T3  = -2D0*mu*J23* ( TEN_DYADPROD(Im, Ci, nsd) +
     2   TEN_DYADPROD(Ci, Im, nsd) ) / 3D0
      CC  = T1 + T2 + T3
      
      RETURN
      END SUBROUTINE PK2_neoHookean
!####################################################################
