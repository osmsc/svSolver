!     Copyright, 2013
!     Mahdi Esmaily Moghadam

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!--------------------------------------------------------------------
!
!     Main routine that contains the calls to all major routines and
!     general structure of the code.
!
!--------------------------------------------------------------------
      
      PROGRAM MAIN
      
      USE COMMOD
      USE ALLFUN
      
      IMPLICIT NONE
      
      LOGICAL l1, l2, l3, l4, l5
      INTEGER i, a, Ac, e, ierr, iEqOld, iBc, eNoN, iM, j, fid
c      INTEGER OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
      REAL(KIND=8) timeP(3)
      CHARACTER(LEN=stdL) fName, tmpS
      
      LOGICAL, ALLOCATABLE :: isS(:)
      INTEGER, ALLOCATABLE :: ptr(:), incL(:)
      REAL(KIND=8), ALLOCATABLE :: xl(:,:), Ag(:,:), al(:,:), Yg(:,:),
     2   yl(:,:), Dg(:,:), dl(:,:), dol(:,:), fNl(:,:), res(:)
      
!--------------------------------------------------------------------
      l1 = .FALSE.
      l2 = .FALSE.
      l3 = .FALSE.
      l4 = .FALSE.
      l5 = .FALSE.
      
      fid       = 27
      savedOnce = .FALSE.
      CALL MPI_INIT(i)
      CALL cm%new(MPI_COMM_WORLD)
      
!     Initiating the exception tracing
      CALL EXCEPTIONS
      
      resetSim  = .FALSE.
      rmsh%cntr = 0
      
!     Reading the user-defined parameters from foo.inp
 101  CALL READFILES
      
!     Doing the partitioning and distributing the data to the all
!     Processors
      CALL DISTRIBUTE
      
!     Initializing the solution vectors and constructing LHS matrix
!     format
      CALL INITIALIZE(timeP)
      
      dbg = 'Allocating intermediate variables'
      ALLOCATE(Ag(tDof,tnNo), Yg(tDof,tnNo), Dg(tDof,tnNo),
     3   res(nFacesLS), incL(nFacesLS), isS(tnNo))
      isS = .FALSE.
      DO Ac=1, tnNo
         IF (ISDOMAIN(1, Ac, phys_struct)) isS(Ac) = .TRUE.
      END DO
      
!--------------------------------------------------------------------
!     Outer loop for marching in time. When entring this loop, all old
!     variables are completely set and satisfy BCs.
      IF (cTS .LE. nITS) dt = dt/1D1
      DO
!     Adjusting the time step size once initialization stage is over
         IF (cTS .EQ. nITS) THEN
            dt = dt*1D1
            std = "New time step size: "//dt
         END IF
!     Incrementing time step, hence cTS will be associated with new
!     variables, i.e. An, Yn, and Dn
         cTS    = cTS + 1
         time   = time + dt
         cEq    = 1
         eq%itr = 0
         eq%ok  = .FALSE.
            
!     Compute mesh properties to check if remeshing is required
         IF (mvMsh .AND. rmsh%isReqd) THEN
            CALL CALCMESHPROPS(nMsh, msh)
            IF (resetSim) EXIT
         END IF
         
!     Predictor step
         CALL PICP
         CALL SETBCDIR(An, Yn, Dn)
         
!     Inner loop for iteration
         DO
            iEqOld = cEq
            
            IF (cplBC%coupled .AND. cEq.EQ.1) THEN
               CALL SETBCCPL
               CALL SETBCDIR(An, Yn, Dn)
            END IF
!     Initiator step
            CALL PICI(Ag, Yg, Dg)
            
            dbg = 'Allocating the RHS and LHS'
            IF (ALLOCATED(R)) THEN
               IF (SIZE(R(:,1)) .NE. dof) THEN
                  DEALLOCATE(R, Val)
                  ALLOCATE (R(dof,tnNo), Val(dof*dof,lhs%nnz))
               END IF
            ELSE
               ALLOCATE (R(dof,tnNo), Val(dof*dof,lhs%nnz))
            END IF
            
!$OMP DO SCHEDULE(GUIDED,mpBs)
            DO i=1, lhs%nnz
               Val(:,i) = 0D0
            END DO
!$OMP END DO
!$OMP DO SCHEDULE(GUIDED,mpBs)
            DO a=1, tnNo
               R(:,a) = 0D0
            END DO
!$OMP END DO
            
            dbg = "Assembling equation <"//eq(cEq)%sym//">"
            DO iM=1, nMsh
               eNoN = msh(iM)%eNoN
               i    = nsd + 2
               j    = 2*nsd + 1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, e, a, Ac, al, yl, dl,
!$OMP&   xl, fNl, dc, ptr, cDmn)
               ALLOCATE(al(tDof,eNoN), yl(tDof,eNoN), dl(tDof,eNoN),
     2            dol(nsd,eNoN), xl(nsd,eNoN), fNl(nsd,eNoN), ptr(eNoN))
!$OMP DO SCHEDULE(GUIDED,mpBs)
               DO e=1, msh(iM)%nEl
                  fNl = 0D0
                  DO a=1, eNoN
                     Ac      = msh(iM)%IEN(a,e)
                     ptr(a)  = Ac
                     xl(:,a) = x (:,Ac)
                     al(:,a) = Ag(:,Ac)
                     yl(:,a) = Yg(:,Ac)
                     dl(:,a) = Dg(:,Ac)
                     IF (mvMsh) THEN
                        dol(:,a) = Do(i:j,Ac)
                     END IF
                     IF (ALLOCATED(fN)) fNl(:,a) = fN(:,Ac)
                  END DO
!     Add contribution of current equation to the LHS/RHS
                  CALL CONSTRUCT(msh(iM), al, yl, dl, dol, xl, fNl, 
     2               ptr, e)
               END DO
               DEALLOCATE(al, yl, dl, xl, dol, fNl, ptr)
               dbg = "Mesh "//iM//" is assembled"
!$OMP END DO
            END DO
            
!     Constructing the element stiffness matrix for boundaries
            CALL SETBCNEU(Yg, Dg)
            incL = 0
            IF (eq(cEq)%phys .EQ. phys_mesh) incL(nFacesLS) = 1
            DO iBc=1, eq(cEq)%nBc
               i = eq(cEq)%bc(iBc)%lsPtr
               IF (i .NE. 0) THEN
                  res(i) = eq(cEq)%gam*dt*eq(cEq)%bc(iBc)%r
                  incL(i) = 1
               END IF
            END DO
!$OMP END PARALLEL

            dbg = "Solving equation <"//eq(cEq)%sym//">"
!     Solving the linear system of equations
! *** HACK!! NATE *** isS doesn't exist in svLS, so calling dummy routine here            
!            CALL svLS_SOLVE(lhs, eq(cEq)%ls, dof, R, Val, isS,
!     2         incL, res)
            CALL svLS_SOLVE(lhs, eq(cEq)%ls, dof, R, Val,
     2         incL, res)

         
!     Solution is obtained, now updating (Corrector)
            CALL PICC
            
!     Checking for exceptions
            CALL EXCEPTIONS
            
!     Writing out the time passed, residual, and etc.
            IF (ALL(eq%ok)) EXIT
            CALL OUTRESULT(timeP, 2, iEqOld)
         END DO
!     End of inner loop
      
!     Saving the TXT files containing average and fluxes
         CALL TXT(.FALSE.)
         
         IF (rmsh%isReqd) THEN
            l1 = MOD(cTS,rmsh%cpVar) .EQ. 0
            IF (l1) THEN
               rmsh%rTS = cTS-1
               rmsh%time = time-dt
               rmsh%iNorm(:) = eq(:)%iNorm
               rmsh%A0(:,:) = Ao(:,:)
               rmsh%Y0(:,:) = Yo(:,:)
               rmsh%D0(:,:) = Do(:,:)
            END IF
         END IF
         
         IF (cm%mas()) INQUIRE(FILE=stopTrigName, EXIST=l1)
         CALL cm%bcast(l1)
         l2 = cTS .GE. nTS
         l3 = MOD(cTS,stFileIncr) .EQ. 0
!     Saving the result to restart simulation in the future
         IF (l1 .OR. l2 .OR. l3) THEN
            fName = TRIM(stFileName)//"_last.bin"
            tmpS  = fName
            IF (.NOT.stFileRepl) THEN
               WRITE(fName,'(I3.3)') cTS
               IF (cTS .GE. 1000) fName = STR(cTS)
               fName = TRIM(stFileName)//"_"//TRIM(fName)//".bin"
            END IF
            IF (cm%mas()) THEN
               OPEN(fid, FILE=TRIM(fName))
               CLOSE(fid, STATUS='DELETE')
            END IF
!     This call is to block all processors
            CALL cm%bcast(l1)
            OPEN(fid, FILE=TRIM(fName), ACCESS='DIRECT', RECL=recLn)
            IF (dFlag) THEN
               IF (rmsh%isReqd .AND. saveAve) THEN
                  WRITE(fid, REC=cm%tF()) stamp, cTS, time,
     2               CPUT()-timeP(1), eq%iNorm, cplBC%xn, Yn, An, Dn,
     3               rmsh%Aav, rmsh%Yav, rmsh%Dav
               ELSE
                  WRITE(fid, REC=cm%tF()) stamp, cTS, time,
     2               CPUT()-timeP(1), eq%iNorm, cplBC%xn, Yn, An, Dn
               END IF
            ELSE
               WRITE(fid, REC=cm%tF()) stamp, cTS, time,
     2            CPUT()-timeP(1), eq%iNorm, cplBC%xn, Yn, An
            END IF
            CLOSE(fid)
            IF (.NOT.stFileRepl .AND. cm%mas()) THEN
               CALL SYSTEM("ln -f "//TRIM(fName)//" "//TRIM(tmpS))
            END IF
         END IF
         
         l3 = MOD(cTS,saveIncr) .EQ. 0
         l4 = saveFormat .NE. saveF_none
         l5 = cTS .GE. saveATS
!     Writing results into the disk with VTU format
         IF (l3 .AND. l4 .AND. l5) THEN
            CALL OUTRESULT(timeP, 3, iEqOld)
            CALL WRITEVTUS(An, Yn, Dn)
            IF (rmsh%isReqd .AND. saveAve) THEN
               rmsh%Aav = rmsh%Aav + An
               rmsh%Yav = rmsh%Yav + Yn
               rmsh%Dav = rmsh%Dav + Dn
            END IF
         ELSE
            CALL OUTRESULT(timeP, 2, iEqOld)
         END IF
         
!     Exiting outer loop if l1 or l2 happens
         IF (l1 .OR. l2) EXIT
         
!     Solution is stored here before replacing it at next time step
         Ao = An
         Yo = Yn
         IF (dFlag) Do = Dn
         cplBC%xo = cplBC%xn
      END DO
!     End of outer loop
      
      IF (resetSim) THEN
         CALL REMESHRESTART(timeP)
         DEALLOCATE(Ag, Yg, Dg, incL, res, isS)
         GOTO 101
      END IF
      
      IF (l2 .AND. saveAve) CALL CALCAVE
      
      CALL FINALIZE
      CALL MPI_FINALIZE(ierr)
      
      END PROGRAM MAIN
      
!--------------------------------------------------------------------
      
      SUBROUTINE STOPSIM()
      
      CALL FINALIZE
c      CALL cm%fStop()
      STOP "MPI is forced to stop by a fatal error"
      
      END SUBROUTINE STOPSIM
