!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!--------------------------------------------------------------------
!      
!     Creating the data structure and assembling LHS sparse matrix.
!      
!--------------------------------------------------------------------

      SUBROUTINE LHSA(nnz)

      USE COMMOD
      USE ALLFUN

      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: nnz

      LOGICAL flag
      INTEGER i, j, a, b, e, mnnzeic, rowN, colN, iM

      INTEGER, ALLOCATABLE :: uInd(:,:)

      flag = .FALSE.
      mnnzeic = 10*MAXVAL(msh%eNoN)
 003  mnnzeic = mnnzeic + MAX(5,mnnzeic/5)

      IF (ALLOCATED(uInd)) DEALLOCATE(uInd)
      ALLOCATE (uInd(mnnzeic,tnNo))
      uInd = 0
      DO iM=1, nMsh
         DO e=1, msh(iM)%nEl
            DO a=1, msh(iM)%eNoN
               rowN = msh(iM)%IEN(a,e)
               DO b=1, msh(iM)%eNoN
                  colN = msh(iM)%IEN(b,e)
                  DO i=1, mnnzeic
!     If current entry is zero, then  fill it up
                     IF (uInd(i,rowN) .EQ. 0) THEN
                        uInd(i,rowN) = colN
                        EXIT
                     END IF
!     If current entry is still smaller, keep going
                     IF (colN .GT. uInd(i,rowN)) CYCLE
!     If there is already a similar entry, no point to add it again
                     IF (colN .EQ. uInd(i,rowN)) EXIT
!     If we are this point, then then the current entry is bigger. Hence
!     we need to shift all the entries from here to the end of the list.
!     If list is full, we request a larger list, otherwise we
!     shift and add the item at the current entry position.
                     IF (uInd(mnnzeic,rowN) .NE. 0) GOTO 003
                     DO j=mnnzeic, i+1, -1
                        uInd(j,rowN) = uInd(j-1,rowN)
                     END DO
                     uInd(i,rowN) = colN
                     EXIT
                  END DO
!     In the following case the list is too short.
                  IF (i .GT. mnnzeic) THEN
                     flag = .TRUE.
                     GOTO 003
                  END IF
               END DO
            END DO
         END DO
      END DO
      IF (flag) std = "max(nnz) is increased to: "//mnnzeic

!--------------------------------------------------------------------
!     Finding number of non-zeros in colPtr vector
      nnz = 0
      DO rowN=1, tnNo
         IF (uInd(1,rowN) .EQ. 0) THEN
            err = "Node "//rowN//" is isolated"
         END IF
         DO i = 1, mnnzeic
            IF (uInd(i,rowN) .NE. 0) THEN
               nnz = nnz + 1
            END IF
         END DO  
      END DO

!--------------------------------------------------------------------
!     Now constructing compact form of rowPtr and colPtr
      ALLOCATE (colPtr(nnz), rowPtr(tnNo+1))
      j  = 1
      rowPtr(1) = 1
      DO rowN=1, tnNo
         DO i=1, mnnzeic
            IF (uInd(i,rowN) .NE. 0) THEN
               colPtr(j) = uInd(i,rowN)
               j = j + 1
            END IF
         END DO
         rowPtr(rowN+1) = j
      END DO
      DEALLOCATE (uInd)

      RETURN
      END SUBROUTINE LHSA

!####################################################################
!     This subroutine assembels the element stiffness matrix into the
!     global stiffness matrix (Val sparse matrix formatted as a vector)
      SUBROUTINE DOASSEM (d, eqN, lK, lR)

      USE COMMOD, ONLY: dof, rowPtr, colPtr, R, Val

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: d, eqN(d)
      REAL(KIND=8), INTENT(IN) :: lK(dof*dof,d,d), lR(dof,d)

      INTEGER a, b, ptr, rowN, colN, left, right

      DO a=1, d
         rowN = eqN(a)
         R(:,rowN) = R(:,rowN) + lR(:,a)
         DO b=1, d
            colN = eqN(b)
            left  = rowPtr(rowN)
            right = rowPtr(rowN+1)
            ptr   = (right + left)/2
            DO WHILE (colN .NE. colPtr(ptr))
               IF (colN .GT. colPtr(ptr)) THEN
                  left  = ptr
               ELSE
                  right = ptr
               END IF
               ptr = (right + left)/2
            END DO
            Val(:,ptr) = Val(:,ptr) + lK(:,a,b)
         END DO
      END DO

      RETURN
      END SUBROUTINE DOASSEM
