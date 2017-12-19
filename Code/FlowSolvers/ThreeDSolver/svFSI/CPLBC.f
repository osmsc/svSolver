!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!--------------------------------------------------------------------
!
!     Here data are received from 3D domain and it setup the data
!     required for integration of ODE's inside cplBC_FINDF
!
!--------------------------------------------------------------------
!
      PROGRAM cplBC_Integ_X

      IMPLICIT NONE

      INCLUDE "cplBC.h"

      INTEGER nTimeSteps, i, j, k, n, nFaces, nX, v, fid, faIn, nXprint
      REAL(KIND=8) pConv, qConv, dt, temp, dt_3D, t, tt
      CHARACTER(LEN=128) cplBC_commu_name

      REAL(KIND=8), ALLOCATABLE :: x(:), Q(:), P(:), f(:,:), offset(:),
     2   xo(:), Xprint(:)
      TYPE(cplFaceType), ALLOCATABLE :: face(:), inFace(:)

      i = COMMAND_ARGUMENT_COUNT()
      IF (i .EQ. 0) THEN
        STOP "0D-3D communication file name is required as an argument"
      ELSEIF (i .GT. 1) THEN
         STOP "Too many arguments"
      END IF
      CALL GETARG(1,cplBC_commu_name)

      fid = 1
!     Reading data required by 0D domain. This includes comparison data.
      OPEN(fid, FILE=cplBC_commu_name, STATUS='OLD', FORM='UNFORMATTED')
      READ(fid) v
      IF (v .NE. cplBCVersion) THEN
         STOP "CPLBC-ERR: 0D code version not compatible with mupfes"
      END IF
      READ(fid) nFaces
      IF (nFaces .EQ. 0) THEN
         STOP "CPLBC-ERROR: 0D code is called with nFaces = 0"
      END IF
      READ(fid) nX
      READ(fid) dt_3D
      READ(fid) t
      ALLOCATE (face(nFaces), inFace(nFaces), P(nFaces), Q(nFaces),     &
     &   f(nX,4), x(nX), offset(nFaces), xo(nX))

      READ(fid) x
      DO faIn=1, nFaces
         READ(fid) inFace(faIn)%bGrp
         READ(fid) inFace(faIn)%Qo
         READ(fid) inFace(faIn)%Qn
         READ(fid) inFace(faIn)%Po
         READ(fid) inFace(faIn)%Pn
         READ(fid) inFace(faIn)%name
      END DO
      CLOSE(fid)

      CALL cplBC_INI(nFaces, nTimeSteps, qConv, pConv, face)

      IF (inFace(1)%eqv .EQ. 0) THEN
         DO i=1, nFaces
            DO j=1, nFaces
               IF (inFace(i)%name .EQ. face(j)%name) THEN
                  inFace(i)%eqv = j
                  inFace(i)%xPtr = face(j)%xPtr
                  IF (inFace(i)%bGrp .NE. face(j)%bGrp) THEN
                     PRINT *, "CPLBC-ERROR: face ", i,                  &
     &                  " bGrp is different in 0D and 3D domains"
                     STOP
                  END IF
                  EXIT
               END IF
            END DO
            IF (j .GT. nFaces) THEN
               PRINT *, "CPLBC-ERROR: face "//TRIM(inFace(i)%name)//    &
     &            " not found in cplBC"
               STOP
            END IF
         END DO
      END IF

      DO i=1, nFaces
         j = inFace(i)%eqv
         face(j)%eqv = i
         face(j)%Qo = inFace(i)%Qo/qConv
         face(j)%Qn = inFace(i)%Qn/qConv
         face(j)%Po = inFace(i)%Po/pConv
         face(j)%Pn = inFace(i)%Pn/pConv
      END DO

      nXprint = 2
      ALLOCATE(Xprint(nXprint))

!--------------------------------------------------------------------
!     Setting up the system of equations
      f      = 0D0
      dt     = dt_3D/REAL(nTimeSteps,8)
      offset = 0D0
      DO n=1, nTimeSteps
         xo = x
         DO i=1, 4
            temp = (REAL(n-1,8) + REAL(i-1,8)/3D0)/REAL(nTimeSteps,8)

            Q = face%Qo + (face%Qn - face%Qo)*temp
            P = face%Po + (face%Pn - face%Po)*temp
            DO j=1, nX
               IF (x(j) .NE. x(j)) THEN
                  PRINT *, "CPLBC-ERROR: NaN detected in cplBC"
                  PRINT *, "Q = ", Q
                  PRINT *, "P = ", p
                  PRINT *, "x = ", x
                  STOP
               END IF
            END DO

            SELECT CASE(i)
            CASE(1)
               tt = t
            CASE(2)
               tt = t + dt/3D0
            CASE(3)
               tt = t + dt*2D0/3D0
            CASE(4)
               tt = t + dt
            END SELECT

            CALL cplBC_FINDF(nFaces, nX, x, f(:,i), Q, P, offset,
     2         nXprint, Xprint, tt)

            SELECT CASE(i)
            CASE(1)
               xo = x
               x = xo + dt*f(:,1)/3D0
            CASE(2)
               x = xo - dt*(f(:,1)/3D0 - f(:,2))
            CASE(3)
               x = xo + dt*(f(:,1) - f(:,2) + f(:,3))
            CASE(4)
               x = xo + dt*(f(:,1) + 3D0*f(:,2) + 3D0*f(:,3)+f(:,4))/8D0
            END SELECT
         END DO
         t = t + dt
      END DO

      DO i=1, nFaces
         j = inFace(i)%xPtr
         k = inFace(i)%eqv
         IF (inFace(i)%bGrp .EQ. cplBC_Dir) THEN
            inFace(i)%y = qConv*(x(j) + offset(k))
         ELSE
            inFace(i)%y = pConv*(x(j) + offset(k))
         END IF
      END DO

!     Writing the data required by the 3D domain
      OPEN(fid, FILE=cplBC_commu_name, STATUS='OLD', FORM='UNFORMATTED')
      WRITE(fid) x
      DO faIn=1, nFaces
         WRITE(fid) inFace(faIn)%y
      END DO
      CLOSE(fid)

      OPEN(fid, FILE='AllData',STATUS='UNKNOWN',POSITION='APPEND')
      DO j=1, nXprint
         WRITE(fid,'(ES14.6E2,A)',ADVANCE='NO') Xprint(j),' '
      END DO
      WRITE(fid,'(A)')
      CLOSE(fid)

      DEALLOCATE(face, inFace, P, Q, f, x, offset, xo, Xprint)

      END PROGRAM cplBC_Integ_X
