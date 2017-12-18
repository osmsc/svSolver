!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!--------------------------------------------------------------------
!      
!     This is to initialize or finalize mupfes variables/structures. 
!      
!--------------------------------------------------------------------
      
      SUBROUTINE INITIALIZE(timeP)
      
      USE COMMOD
      USE ALLFUN
      
      IMPLICIT NONE
      
      REAL(KIND=8), INTENT(OUT) :: timeP(3)
      
      LOGICAL :: flag
      INTEGER :: iEq, ierr, gnnz, nnz, iDmn
      CHARACTER(LEN=stdL) :: fName
      REAL(KIND=8) :: am
      TYPE(memLS_commuType) :: communicator
      
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: s
      
      INTEGER :: iM, i
      CHARACTER(LEN=stdL) :: sTmp, fTmp
      
      tDof     = 0
      dFlag    = .FALSE.
      nFacesLS = SUM(eq%nBc)
      IF (mvMsh) nFacesLS = nFacesLS + 1
      IF (ANY(eq(1)%bc%cplBCptr .NE. 0)) cplBC%coupled = .TRUE.
      
      DO iEq=1, nEq
!     This would be the default value of am for first order equations
         eq(iEq)%am = 5D-1*(3D0 - eq(iEq)%roInf) / 
     2                (1D0 + eq(iEq)%roInf)
!     If the equation is second order, am is calculated as follows
         am = (2D0 - eq(iEq)%roInf)/(1D0 + eq(iEq)%roInf)
         SELECT CASE (eq(iEq)%phys)
         CASE (phys_fluid)
            eq(iEq)%dof = nsd + 1
            eq(iEq)%sym = 'NS'
         CASE (phys_heatF)
            eq(iEq)%dof = 1
            eq(iEq)%sym = 'HF'
         CASE (phys_heatS)
            eq(iEq)%dof = 1
            eq(iEq)%sym = 'HS'
         CASE (phys_lElas)
            dFlag = .TRUE.
            eq(iEq)%dof = nsd
            eq(iEq)%am  = am
            eq(iEq)%sym = 'LE'
         CASE (phys_struct)
            dFlag = .TRUE.
            eq(iEq)%dof = nsd
            eq(iEq)%am  = am
            eq(iEq)%sym = 'ST'
         CASE (phys_FSI)
            dFlag = .TRUE.
            eq(iEq)%dof = nsd + 1
            eq(iEq)%sym = 'FS'
         CASE (phys_mesh)
            dFlag = .TRUE.
            eq(iEq)%dof = nsd
            eq(iEq)%am  = am
            eq(iEq)%sym = 'MS'
         CASE (phys_BBO)
            eq(iEq)%dof = nsd
            eq(iEq)%sym = 'BB'
         CASE DEFAULT
            err = "Equation type "//eq(iEq)%sym//" is not defined"
         END SELECT
         
         eq(iEq)%pNorm = HUGE(eq(iEq)%pNorm)
         eq(iEq)%af    = 1D0/(1D0 + eq(iEq)%roInf)
         eq(iEq)%beta  = 25D-2*(1D0 + eq(iEq)%am - eq(iEq)%af)**2D0
         eq(iEq)%gam   = 5D-1 + eq(iEq)%am - eq(iEq)%af
         eq(iEq)%s     = tDof + 1
         eq(iEq)%e     = tDof + eq(iEq)%dof
         tDof          = eq(iEq)%e
      END DO
      ierr = 0; IF (dFlag) ierr = 1
      i = 0; IF (cplBC%coupled) i = cplBC%nX
      
      stamp = (/cm%np(), nEq, nMsh, tnNo, i, tDof, ierr, version/)
      
!     Calculating the record length
      i = 2
      IF (dFlag) i = 3
      IF (rmsh%isReqd .AND. saveAve) i = 6
      i = 4*(1+SIZE(stamp)) + 8*(2+nEq+cplBC%nX+i*tDof*tnNo)
      IF (cm%seq()) THEN
         recLn = i
      ELSE
         CALL MPI_ALLREDUCE(i, recLn, 1, mpint, MPI_MAX, cm%com(), ierr)
      END IF
      
      ALLOCATE(Ao(tDof,tnNo), An(tDof,tnNo), Yo(tDof,tnNo), 
     2   Yn(tDof,tnNo), Do(tDof,tnNo), Dn(tDof,tnNo))
      
      IF (.NOT.resetSim) THEN
         IF (.NOT.ALLOCATED(rmsh%flag)) ALLOCATE(rmsh%flag(nMsh))
         rmsh%flag(:) = .FALSE.
         rmsh%fTS = rmsh%freq
         IF (rmsh%isReqd) THEN
            ALLOCATE(rmsh%A0(tDof,tnNo))
            ALLOCATE(rmsh%Y0(tDof,tnNo))
            ALLOCATE(rmsh%D0(tDof,tnNo))
            ALLOCATE(rmsh%iNorm(nEq))
            IF (saveAve) THEN
               ALLOCATE(rmsh%Aav(tDof,tnNo))
               ALLOCATE(rmsh%Yav(tDof,tnNo))
               ALLOCATE(rmsh%Dav(tDof,tnNo))
            END IF
         END IF
      END IF
      
      std = " Constructing stiffness matrix sparse structure"
      CALL LHSA(nnz)
       
      gnnz = nnz
      CALL MPI_ALLREDUCE(nnz, gnnz, 1, mpint, MPI_SUM, cm%com(), ierr)
      
      IF (resetSim) THEN
         IF (communicator%foC) CALL memLS_COMMU_FREE(communicator)
         IF (lhs%foC) CALL memLS_LHS_FREE(lhs)
      END IF ! resetSim
      
      std = " Total number of non-zeros in the LHS matrix: "//gnnz
      dbg = "Calling memLS_COMMU_CREATE"
      CALL memLS_COMMU_CREATE(communicator, cm%com())
      
      dbg = "Calling memLS_LHS_CREATE"
      CALL memLS_LHS_CREATE(lhs, communicator, gtnNo, tnNo, nnz, ltg,
     2   rowPtr, colPtr, nFacesLS)
      
      IF (.NOT.resetSim) THEN
         IF (stFileFlag) THEN
            INQUIRE (FILE=iniFilePath, EXIST=flag)
            IF (flag) THEN
               i = LEN(TRIM(iniFilePath))
               IF (iniFilePath(i-2:i) .EQ. "bin") THEN
                  CALL INITFROMBIN(iniFilePath)
               ELSE
                  CALL INITFROMVTK(iniFilePath)
               END IF
            ELSE
               fName = TRIM(stFileName)//"_last.bin"
               INQUIRE (FILE=fName, EXIST=flag)
               IF (flag) THEN
                  CALL INITFROMBIN(fName)
               ELSE
                  IF (cm%mas()) wrn = TRIM(fName)//" can not be opened"
                  CALL ZEROINIT
               END IF
            END IF
            IF (rmsh%isReqd) THEN
               rmsh%fTS = (cTS/rmsh%fTS + 1)*rmsh%freq
               rmsh%rTS = cTS
               rmsh%time = time
               rmsh%iNorm(:) = eq(:)%iNorm
               rmsh%A0(:,:) = Ao(:,:)
               rmsh%Y0(:,:) = Yo(:,:)
               rmsh%D0(:,:) = Do(:,:)
               IF (saveAve .AND. zeroAve) THEN
                  rmsh%Aav = 0D0
                  rmsh%Yav = 0D0
                  rmsh%Dav = 0D0
               END IF
            END IF
         ELSE
            CALL ZEROINIT
         END IF ! stFileFlag
         rsTS = cTS
      ELSE
         cTS  = rmsh%rTS
         time = rmsh%time
         eq(:)%iNorm = rmsh%iNorm(:)
         Ao = LOCAL(rmsh%A0)
         Yo = LOCAL(rmsh%Y0)
         Do = LOCAL(rmsh%D0)
         IF (saveAve) THEN
            rmsh%A0 = rmsh%Aav
            rmsh%Y0 = rmsh%Yav
            rmsh%D0 = rmsh%Dav
            DEALLOCATE(rmsh%Aav, rmsh%Yav, rmsh%Dav)
            ALLOCATE(rmsh%Aav(tDof,tnNo))
            ALLOCATE(rmsh%Yav(tDof,tnNo))
            ALLOCATE(rmsh%Dav(tDof,tnNo))
            rmsh%Aav = LOCAL(rmsh%A0)
            rmsh%Yav = LOCAL(rmsh%Y0)
            rmsh%Dav = LOCAL(rmsh%D0)
         END IF
         DEALLOCATE(rmsh%A0,rmsh%Y0,rmsh%D0)
         ALLOCATE(rmsh%A0(tDof,tnNo)); rmsh%A0(:,:) = Ao(:,:)
         ALLOCATE(rmsh%Y0(tDof,tnNo)); rmsh%Y0(:,:) = Yo(:,:)
         ALLOCATE(rmsh%D0(tDof,tnNo)); rmsh%D0(:,:) = Do(:,:)
      END IF ! resetSim
      
      DO iM=1, nMsh
         IF (cm%mas()) THEN
            fTmp = TRIM(appPath)//".partitioning_"//
     2         TRIM(msh(iM)%name)//".bin"
            sTmp = TRIM(appPath)//".partitioning_"//
     2         TRIM(msh(iM)%name)//"_"//STR(cTS)//".bin"
            INQUIRE(FILE=TRIM(fTmp), EXIST=flag)
            IF (flag) THEN
               sTmp = "cp  "//TRIM(fTmp)//" "//TRIM(sTmp)
               CALL SYSTEM(TRIM(sTmp))
            END IF
         END IF
      END DO
      
!     Calculating the volume of each domain
      ALLOCATE(s(1,tnNo))
      s = 1D0
      DO iEq=1, nEq
         DO iDmn=1, eq(iEq)%nDmn
            eq(iEq)%dmn(iDmn)%v = Integ(eq(iEq)%dmn(iDmn)%Id, s, 1, 1)
            IF (ISZERO(eq(iEq)%dmn(iDmn)%v)) wrn = "Volume of "//
     2         "domain "//iDmn//" of equation "//iEq//" is zero"
            IF (eq(iEq)%dmn(iDmn)%phys .EQ. phys_struct) THEN
               CALL TEN_INIT(nsd)
            END IF
         END DO
      END DO
      
!     Predicting new variables
      CALL PICP
      
!     Preparing faces and BCs
      CALL BAFINI()
      
!     Making sure the old solution satisfies BCs
      CALL SETBCDIR(Ao, Yo, Do)
      
!     Preparing TXT files
      CALL TXT(.TRUE.)
      
!     Printing the first line and initializing timeP
      CALL OUTRESULT(timeP, 1, 1)
      rmsh%flag(:) = .FALSE.
      resetSim = .FALSE.
      
      RETURN
      CONTAINS
!--------------------------------------------------------------------
!     Initializing accelaration, velocity and displacement to zero
      SUBROUTINE ZEROINIT
      
      IMPLICIT NONE
      
      std = " Initializing state variables to zero"
      
!     This cTS corresponds to old variables. As soon as incrementing it
!     by one, it will be associated to new variables.
      cTS      = 0
      time     = 0D0
      timeP(1) = 0D0
      eq%iNorm = 0D0
      
      Ao = 0D0
      Yo = 0D0
      Do = 0D0
      
      IF (rmsh%isReqd) THEN
         rmsh%A0 = 0D0
         rmsh%Y0 = 0D0
         rmsh%D0 = 0D0
         IF (saveAve) THEN
            rmsh%Aav = 0D0
            rmsh%Yav = 0D0
            rmsh%Dav = 0D0
         END IF
      END IF
      
      RETURN
      END SUBROUTINE ZEROINIT
      
!--------------------------------------------------------------------
!     Using the saved VTK files for initialization
      SUBROUTINE INITFROMVTK(fName)
      
      IMPLICIT NONE
      
      CHARACTER(LEN=stdL), INTENT(IN) :: fName
      
      INTEGER l
      
      REAL(KIND=8), ALLOCATABLE :: tmpA(:,:), tmpY(:,:), tmpD(:,:)
      
      err = "Initialization from vtk is deprecated"
      
      l = LEN(TRIM(fName))
      IF (fName(l-2:l) .NE. "vtk") err = "Format of <"//
     2   TRIM(fName)//"> is not recognized"
      
      IF (ANY(msh%eType .EQ. eType_NRB)) err = "Only isoparametric"
     2   //" meshes can be read from VTK files"
      
      IF (nMsh .GT. 1) wrn = "Reading from VTK will fail in"//
     2   " presence of projected boundaries"
      std = " Initializing from "//fName
      
      cTS      = 0
      time     = 0D0
      timeP(1) = 0D0
      eq%iNorm = 0D0
      
      IF (cm%mas()) THEN
         ALLOCATE(tmpA(tDof,gtnNo), tmpY(tDof,gtnNo), tmpD(tDof,gtnNo))
!         CALL READVTK(fName, tmpA, tmpY, tmpD)
      ELSE
         ALLOCATE(tmpA(0,0), tmpY(0,0), tmpD(0,0))
      END IF
      
      Ao = LOCAL(tmpA)
      Yo = LOCAL(tmpY)
      Do = LOCAL(tmpD)
      
      RETURN
      END SUBROUTINE INITFROMVTK
      
!--------------------------------------------------------------------
!     Using the MUPFES specific format binary file for initialization
      SUBROUTINE INITFROMBIN(fName)
      
      IMPLICIT NONE
      
      CHARACTER(LEN=stdL), INTENT(IN) :: fName
      
      INTEGER, PARAMETER :: fid = 1
      
      INTEGER tStamp(SIZE(stamp))
      
      std = " Initializing from "//fName
      OPEN(fid, FILE=fName, ACCESS='DIRECT', RECL=recLn)
      IF (dFlag) THEN
         IF (rmsh%isReqd .AND. saveAve) THEN
            READ(fid,REC=cm%tF()) tStamp, cTS, time, timeP(1), eq%iNorm, 
     2         cplBC%xo, Yo, Ao, Do, rmsh%Aav, rmsh%Yav, rmsh%Dav
         ELSE
            READ(fid,REC=cm%tF()) tStamp, cTS, time, timeP(1), eq%iNorm, 
     2         cplBC%xo, Yo, Ao, Do
         END IF
      ELSE
         READ(fid,REC=cm%tF()) tStamp, cTS, time, timeP(1), eq%iNorm, 
     2      cplBC%xo, Yo, Ao
      END IF
      CLOSE(fid)
      
!     First checking all variables on master processor, since on the
!     other processor data will be shifted due to any change on the
!     sizes
      IF (cm%mas()) THEN
         IF (tStamp(1).NE.stamp(1)) err = "Number of processors <"//
     2      tStamp(1)//"> does not match with "//
     3      TRIM(fName)//" <"//stamp(1)//">"
         IF (tStamp(2).NE.stamp(2)) err = "Number of equations <"//
     2      tStamp(2)//"> does not match with "//
     3      TRIM(fName)//" <"//stamp(2)//">"
         IF (tStamp(3).NE.stamp(3)) err = "Number of 0D unknowns"//
     2      tStamp(3)//"> does not match with "//
     3      TRIM(fName)//" <"//stamp(3)//">"
         IF (tStamp(4).NE.stamp(4)) err = "Number of total dof <"//
     2      tStamp(4)//"> does not match with "//
     3      TRIM(fName)//" <"//stamp(4)//">"
         IF (tStamp(5).NE.stamp(5)) err = "Number of variables"//
     2      tStamp(5)//"> does not match with "//
     3      TRIM(fName)//" <"//stamp(5)//">"
         IF (tStamp(6).NE.stamp(6)) err = "Number of elements"//
     2      tStamp(6)//"> does not match with "//
     3      TRIM(fName)//" <"//stamp(6)//">"
         IF (tStamp(7).NE.stamp(7)) err = "Number of nodes"//
     2      tStamp(7)//"> does not match with "//
     3      TRIM(fName)//" <"//stamp(7)//">"
         IF (tStamp(8).NE.stamp(8)) err = "Version of solver"//
     2      tStamp(8)//"> does not match with "//
     3      TRIM(fName)//" <"//stamp(8)//">"
      END IF
      
      CALL MPI_BARRIER(cm%com(), ierr)
      IF (ANY(tStamp.NE.stamp)) err = "Simulation stamp"//
     2   " does not match with "//fName
      
      RETURN
      END SUBROUTINE INITFROMBIN
      
      END SUBROUTINE INITIALIZE
!####################################################################
      SUBROUTINE FINALIZE
      
      USE COMMOD
      USE ALLFUN
      
      IMPLICIT NONE
      
      INTEGER iM, iEq
      
!     Deallocating meshes
      IF (ALLOCATED(msh)) THEN
         DO iM=1, nMsh
            CALL DESTROY(msh(iM))
         END DO
         DEALLOCATE(msh)
      END IF
      
!     Deallocating equations
      IF (ALLOCATED(eq)) THEN
         DO iEq=1, nEq
            CALL DESTROY(eq(iEq))
         END DO
         DEALLOCATE(eq)
      END IF
      
!     Deallocating sparse matrix structures
      IF(lhs%foc) CALL memLS_LHS_FREE(lhs)
      
!     Closing the output channels
      CALL std%close()
      CALL wrn%close()
      CALL err%close()
      CALL dbg%close()
      
      RETURN
      END SUBROUTINE FINALIZE
