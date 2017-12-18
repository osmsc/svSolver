!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!--------------------------------------------------------------------
!      
!     This routine contains multiple functions that are generally
!     desined to interface with user.
!      
!--------------------------------------------------------------------
!     This is to find if any exception occures, commenting this out
!     untill fortran2003 is standard in all compilers
      SUBROUTINE EXCEPTIONS

c      USE IEEE_EXCEPTIONS 
      USE COMMOD, ONLY: stdL

      IMPLICIT NONE

!     I am assuming nExCk is 5, if a compilation error occured, you need
!     to adjust the following arrays accordingely
c      INTEGER, PARAMETER :: nExCk = SIZE(IEEE_ALL) 
c      LOGICAL, PARAMETER :: check(nExCk) = 
c     2   (/.TRUE.,.TRUE.,.TRUE.,.TRUE.,.FALSE./)
c      CHARACTER(LEN=stdL), PARAMETER :: ieWarn(nExCk) = 
c     2   (/"Overflow", "Divide by zero", "Invalid arithmetic operation",
c     3   "Underflow", "Inexact operation"/)
c
c      LOGICAL, SAVE :: iniSet = .TRUE., sprtFlag(nExCk)
c
c      LOGICAL fg
c      INTEGER i
c      REAL(KIND=8) r
c
c      IF (iniSet) THEN
c         iniSet = .FALSE.
c         DO i=1, nExCk
c            IF (.NOT.check(i)) CYCLE
c            sprtFlag(i) = IEEE_SUPPORT_FLAG(IEEE_ALL(i), r)
c         END DO
c         
c         fg = .FALSE.
c         DO i=1, nExCk
c            IF (.NOT.check(i)) CYCLE
c            IF (sprtFlag(i)) CALL IEEE_SET_HALTING_MODE(IEEE_ALL(i), fg)
c         END DO
c      END IF
c
c      DO i=1, nExCk
c         IF (.NOT.check(i)) CYCLE
c         IF (sprtFlag(i)) THEN
c            CALL IEEE_GET_FLAG(IEEE_ALL(i), fg)
c            IF (fg) THEN
c               CALL WARNING(ieWarn(i))
c               CALL IEEE_SET_FLAG(IEEE_ALL(i), .FALSE.)
c            END IF
c         END IF
c      END DO

      RETURN
      END SUBROUTINE EXCEPTIONS

!####################################################################
!     Prepares the output of mupfes to the standard output.
      SUBROUTINE OUTRESULT(timeP, co, iEq)

      USE COMMOD
      USE ALLFUN

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: co, iEq
      REAL(KIND=8), INTENT(INOUT) :: timeP(3)

      CHARACTER(LEN=*), PARAMETER :: sepLine = "---------------------"//
     2   "-------------------------------------"

      INTEGER fid, i
      REAL(KIND=8) tmp, tmp2
      CHARACTER c1, c2
      CHARACTER(LEN=stdL) sOut

      IF (cm%slv()) RETURN
      
      fid = 1
      tmp = CPUT()
      
      IF (co .EQ. 1) THEN
         timeP(1) = tmp - timeP(1)
         timeP(2) = 0D0
         std = " "
         std = sepLine
         std = "Eq     N-i     T      dB   Ri/R0    R/Ri     lsIt  "//
     2      "dB  %t"
         IF (nEq .EQ. 1) std = sepLine
         RETURN
      END IF

      IF (nEq.GT.1 .AND. iEq.EQ.1 .AND. eq(iEq)%itr.EQ.1) THEN
         std = sepLine
      END IF

      c1 = " "
      IF (co .EQ. 3) c1 = "s"
      timeP(3) = tmp - timeP(1)
      sOut = eq(iEq)%sym//" "//STR(cTS,5)//"-"//STR(eq(iEq)%itr)//c1//
     2   " "//STR(timeP(3),6)

      IF (ISZERO(eq(iEq)%iNorm)) THEN
         tmp  = 1D0
         tmp2 = 1D0
         i    = 0
      ELSE
         tmp  = eq(iEq)%ls%RI%iNorm/eq(iEq)%iNorm
         tmp2 = eq(iEq)%ls%RI%fNorm/eq(iEq)%ls%RI%iNorm
         i    = INT(2D1*LOG10(tmp/eq(iEq)%pNorm))
      END IF

      IF (i .GT. 20) THEN
         c1 = "!"; c2 = "!"
      ELSE
         c1 = "["; c2 = "]"
      END IF
      sOut = TRIM(sOut)//"  "//c1//STR(i,4)//" "//STR(tmp,7)//" "//
     2   STR(tmp2,7)//c2

      IF (ISZERO(timeP(3),timeP(2))) timeP(3) = (1D0+eps)*timeP(2) + eps
      tmp = 1D2*eq(iEq)%ls%RI%callD/(timeP(3) - timeP(2))
      timeP(2) = timeP(3)
      IF (ABS(tmp) .GT. 1D2) tmp = 1D2

      IF (eq(iEq)%ls%RI%suc) THEN
         c1 = "["; c2 = "]"
      ELSE
         c1 = "!"; c2 = "!"
      END IF
      sOut = TRIM(sOut)//"  "//c1//STR(eq(iEq)%ls%RI%itr,4)//" "//
     2   STR(NINT(eq(iEq)%ls%RI%dB),3)//" "//STR(NINT(tmp),3)//c2

      IF (nEq .GT. 1) THEN
         std = CLR(sOut,iEq)
         IF (sepOutput .AND. iEq.NE.cEq) std = "--"
      ELSE
         std = sOut
      END IF
      
      RETURN
      END SUBROUTINE OUTRESULT


