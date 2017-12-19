!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!      
!--------------------------------------------------------------------
!      
!     This is a header file for the data structure used for 
!     communication with mupfes
!      
!--------------------------------------------------------------------

      INTEGER, PARAMETER :: cplBC_Dir = 66112, cplBC_Neu = 66113,
     2   cplBCVersion = 8
     
      TYPE cplFaceType
         SEQUENCE
         INTEGER bGrp            ! (IN)  geBC_Dir/genBC_Neu
         INTEGER xPtr            ! (USE) pointer to x
         INTEGER :: eqv = 0      ! (USE) internal genBC use
         INTEGER reserved        ! ( - ) reserved for alignment
         REAL(KIND=8) Qo         ! (IN)  flow rate at t
         REAL(KIND=8) Qn         ! (IN)  flow rate at t+dt
         REAL(KIND=8) Po         ! (IN)  pressure at t
         REAL(KIND=8) Pn         ! (IN)  pressure at t+dt
         REAL(KIND=8) y          ! (OUT) imposed flow/pressure
         CHARACTER(LEN=128) name ! (IN)  name of the face
      END TYPE cplFaceType
     
