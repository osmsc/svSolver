!--------------------------------------------------------------------
!     Created by Mahdi Esmaily Moghadam
!     contact memt63@gmail.com for reporting the bugs.
!--------------------------------------------------------------------
!
!     UC Copyright Notice
!     This software is Copyright ©2012 The Regents of the University of
!     California. All Rights Reserved.
!
!     Permission to copy and modify this software and its documentation
!     for educational, research and non-profit purposes, without fee, 
!     and without a written agreement is hereby granted, provided that
!     the above copyright notice, this paragraph and the following three
!     paragraphs appear in all copies.
!
!     Permission to make commercial use of this software may be obtained
!     by contacting:
!     Technology Transfer Office
!     9500 Gilman Drive, Mail Code 0910
!     University of California
!     La Jolla, CA 92093-0910
!     (858) 534-5815
!     invent@ucsd.edu
!
!     This software program and documentation are copyrighted by The
!     Regents of the University of California. The software program and
!     documentation are supplied "as is", without any accompanying
!     services from The Regents. The Regents does not warrant that the
!     operation of the program will be uninterrupted or error-free. The
!     end-user understands that the program was developed for research
!     purposes and is advised not to rely exclusively on the program for
!     any reason.
!
!     IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY 
!     PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL 
!     DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS 
!     SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF 
!     CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
!     THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY 
!     WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
!     OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE 
!     SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE 
!     UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE 
!     MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
!
!--------------------------------------------------------------------
!     Routines prototypes are declared here. This provides the details 
!     of how to call these functions
!--------------------------------------------------------------------
      
      INTERFACE
         SUBROUTINE svLS_COMMU_CREATE(commu, commi)
            INCLUDE "svLS_STRUCT.h"
            TYPE(svLS_commuType), INTENT(INOUT) :: commu
            INTEGER, INTENT(IN) :: commi
         END SUBROUTINE svLS_COMMU_CREATE
         SUBROUTINE svLS_COMMU_FREE(commu)
            INCLUDE "svLS_STRUCT.h"
            TYPE(svLS_commuType), INTENT(INOUT) :: commu
         END SUBROUTINE svLS_COMMU_FREE
         
         SUBROUTINE svLS_LHS_CREATE(lhs, commu, gnNo, nNo, nnz,        &
     &      gNodes, rowPtr, colPtr, nFaces)
            INCLUDE "svLS_STRUCT.h"
            TYPE(svLS_lhsType), INTENT(INOUT) :: lhs
            TYPE(svLS_commuType), INTENT(IN) :: commu
            INTEGER, INTENT(IN) :: gnNo, nNo, nnz
            INTEGER, INTENT(IN) :: gNodes(nNo), rowPtr(nNo+1),          &
     &         colPtr(nnz)
            INTEGER, INTENT(IN) :: nFaces
         END SUBROUTINE svLS_LHS_CREATE
         SUBROUTINE svLS_LHS_FREE(lhs)
            INCLUDE "svLS_STRUCT.h"
            TYPE(svLS_lhsType), INTENT(INOUT) :: lhs
         END SUBROUTINE svLS_LHS_FREE

         SUBROUTINE svLS_LS_CREATE(ls, LS_type, relTol, absTol, maxItr,&
     &      dimKry, relTolIn, absTolIn, maxItrIn)
            INCLUDE "svLS_STRUCT.h"
            TYPE(svLS_lsType), INTENT(INOUT) :: ls
            INTEGER, INTENT(IN) :: LS_type
            REAL(KIND=8), INTENT(IN), OPTIONAL :: relTol, absTol,       &
     &         relTolIn(2), absTolIn(2)
            INTEGER, INTENT(IN), OPTIONAL :: maxItr, dimKry, maxItrIn(2)
         END SUBROUTINE svLS_LS_CREATE
         SUBROUTINE svLS_LS_FREE(ls)
            INCLUDE "svLS_STRUCT.h"
            TYPE(svLS_lsType), INTENT(INOUT) :: ls
         END SUBROUTINE svLS_LS_FREE
         
         SUBROUTINE svLS_BC_CREATE(lhs, faIn, nNo, dof, BC_type,       &
     &         gNodes, Val)
            INCLUDE "svLS_STRUCT.h"
            TYPE(svLS_lhsType), INTENT(INOUT) :: lhs
            INTEGER, INTENT(IN) :: faIn, nNo, dof
            INTEGER, INTENT(IN) :: BC_type
            INTEGER, INTENT(IN) :: gNodes(nNo)
            REAL(KIND=8), INTENT(IN), OPTIONAL :: Val(dof,nNo)
         END SUBROUTINE svLS_BC_CREATE
         SUBROUTINE svLS_BC_FREE(lhs, faIn)
            INCLUDE "svLS_STRUCT.h"
            TYPE(svLS_lhsType), INTENT(INOUT) :: lhs
            INTEGER, INTENT(IN) :: faIn
         END SUBROUTINE svLS_BC_FREE
         
         SUBROUTINE svLS_SOLVE (lhs, ls, dof, Ri, Val, isS, incL, res)
            INCLUDE "svLS_STRUCT.h"
            TYPE(svLS_lhsType), INTENT(INOUT) :: lhs
            TYPE(svLS_lsType), INTENT(INOUT) :: ls
            INTEGER, INTENT(IN) :: dof
            REAL(KIND=8), INTENT(INOUT) :: Ri(dof,lhs%nNo)
            REAL(KIND=8), INTENT(IN) :: Val(dof*dof,lhs%nnz)
            LOGICAL, INTENT(IN), OPTIONAL :: isS(lhs%nNo)
            INTEGER, INTENT(IN), OPTIONAL :: incL(lhs%nFaces)
            REAL(KIND=8), INTENT(IN), OPTIONAL :: res(lhs%nFaces)
         END SUBROUTINE svLS_SOLVE
      END INTERFACE

      INTERFACE svLS_INP
         FUNCTION svLS_DOTV(dof, nNo, commu, U, V)
            INCLUDE "svLS_STRUCT.h"
            INTEGER, INTENT(IN) :: dof, nNo
            TYPE(svLS_commuType), INTENT(IN) :: commu
            REAL(KIND=8), INTENT(IN) :: V(dof,nNo), U(dof,nNo)
            REAL(KIND=8) svLS_DOTV
         END FUNCTION svLS_DOTV
         FUNCTION svLS_DOTS(nNo, commu, U, V)
            INCLUDE "svLS_STRUCT.h"
            INTEGER, INTENT(IN) :: nNo
            TYPE(svLS_commuType), INTENT(IN) :: commu
            REAL(KIND=8), INTENT(IN) :: V(nNo), U(nNo)
            REAL(KIND=8) svLS_DOTS
         END FUNCTION svLS_DOTS
         FUNCTION svLS_NORMV(dof, nNo, commu, U)
            INCLUDE "svLS_STRUCT.h"
            INTEGER, INTENT(IN) :: dof, nNo
            TYPE(svLS_commuType), INTENT(IN) :: commu
            REAL(KIND=8), INTENT(IN) :: U(dof,nNo)
            REAL(KIND=8) svLS_NORMV
         END FUNCTION svLS_NORMV
         FUNCTION svLS_NORMS(nNo, commu, U)
            INCLUDE "svLS_STRUCT.h"
            INTEGER, INTENT(IN) :: nNo
            TYPE(svLS_commuType), INTENT(IN) :: commu
            REAL(KIND=8), INTENT(IN) :: U(nNo)
            REAL(KIND=8) svLS_NORMS
         END FUNCTION svLS_NORMS
      END INTERFACE svLS_INP

      INTERFACE svLS_COMMU
         SUBROUTINE svLS_COMMUV(lhs, dof, R)
            INCLUDE "svLS_STRUCT.h"
            TYPE(svLS_lhsType), INTENT(INOUT) :: lhs
            INTEGER, INTENT(IN) :: dof
            REAL(KIND=8), INTENT(INOUT) :: R(dof,lhs%nNo)
         END SUBROUTINE svLS_COMMUV
         SUBROUTINE svLS_COMMUS(lhs, R)
            INCLUDE "svLS_STRUCT.h"
            TYPE(svLS_lhsType), INTENT(INOUT) :: lhs
            REAL(KIND=8), INTENT(INOUT) :: R(lhs%nNo)
         END SUBROUTINE svLS_COMMUS
      END INTERFACE svLS_COMMU

      INTERFACE svLS_SPARMUL
         SUBROUTINE svLS_SPARMULVV(lhs, rowPtr, colPtr, dof, K, U, KU)
            INCLUDE "svLS_STRUCT.h"
            TYPE(svLS_lhsType), INTENT(INOUT) :: lhs
            INTEGER, INTENT(IN) :: rowPtr(2,lhs%nNo), colPtr(lhs%nnz)
            INTEGER, INTENT(IN) :: dof
            REAL(KIND=8), INTENT(IN) :: K(dof*dof,lhs%nnz),             &
     &         U(dof,lhs%nNo)
            REAL(KIND=8), INTENT(OUT) :: KU(dof,lhs%nNo)
         END SUBROUTINE svLS_SPARMULVV
         SUBROUTINE svLS_SPARMULVS(lhs, rowPtr, colPtr, dof, K, U, KU)
            INCLUDE "svLS_STRUCT.h"
            TYPE(svLS_lhsType), INTENT(INOUT) :: lhs
            INTEGER, INTENT(IN) :: rowPtr(2,lhs%nNo), colPtr(lhs%nnz)
            INTEGER, INTENT(IN) :: dof
            REAL(KIND=8), INTENT(IN) :: K(dof,lhs%nnz), U(dof,lhs%nNo)
            REAL(KIND=8), INTENT(OUT) :: KU(lhs%nNo)
         END SUBROUTINE svLS_SPARMULVS
         SUBROUTINE svLS_SPARMULSV(lhs, rowPtr, colPtr, dof, K, U, KU)
            INCLUDE "svLS_STRUCT.h"
            TYPE(svLS_lhsType), INTENT(INOUT) :: lhs
            INTEGER, INTENT(IN) :: rowPtr(2,lhs%nNo), colPtr(lhs%nnz)
            INTEGER, INTENT(IN) :: dof
            REAL(KIND=8), INTENT(IN) :: K(dof,lhs%nnz), U(lhs%nNo)
            REAL(KIND=8), INTENT(OUT) :: KU(dof,lhs%nNo)
         END SUBROUTINE svLS_SPARMULSV
         SUBROUTINE svLS_SPARMULSS(lhs, rowPtr, colPtr, K, U, KU)
            INCLUDE "svLS_STRUCT.h"
            TYPE(svLS_lhsType), INTENT(INOUT) :: lhs
            INTEGER, INTENT(IN) :: rowPtr(2,lhs%nNo), colPtr(lhs%nnz)
            REAL(KIND=8), INTENT(IN) :: K(lhs%nnz), U(lhs%nNo)
            REAL(KIND=8), INTENT(OUT) :: KU(lhs%nNo)
         END SUBROUTINE svLS_SPARMULSS
      END INTERFACE svLS_SPARMUL

      INTERFACE svLS_CPUT
         FUNCTION svLS_CPUT()
            REAL(KIND=8) svLS_CPUT
         END FUNCTION svLS_CPUT
      END INTERFACE svLS_CPUT

      INTERFACE svLS_HRCPUT
         FUNCTION svLS_HRCPUT()
            REAL(KIND=8) svLS_HRCPUT
         END FUNCTION svLS_HRCPUT
      END INTERFACE svLS_HRCPUT
