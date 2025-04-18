      MODULE pointerAssignments
      
      IMPLICIT NONE

      TYPE :: allocVal
         CHARACTER(LEN=:), ALLOCATABLE :: val
      END TYPE allocVal
      
      CONTAINS
      
      SUBROUTINE assignINT(pnt,val)
      INTEGER, POINTER, INTENT(INOUT) :: pnt
      CHARACTER(LEN=*), INTENT(IN) :: val
      READ(val, *) pnt
      NULLIFY(pnt)
      END SUBROUTINE assignINT

      SUBROUTINE assignINTarray(pnt,valArr)
      INTEGER, DIMENSION(:), POINTER, INTENT(INOUT) :: pnt
      TYPE(allocVal), DIMENSION(:), INTENT(IN) :: valArr
      INTEGER :: index
      DO index = 1,SIZE(valArr)
         READ(valArr(index)%val, *) pnt(index)
      END DO
      NULLIFY(pnt)
      END SUBROUTINE assignINTarray

      SUBROUTINE assignLOGICAL(pnt,val)
      LOGICAL, POINTER, INTENT(INOUT) :: pnt
      CHARACTER(LEN=*), INTENT(IN) :: val
      IF (val.EQ.".TRUE.") THEN
         pnt = .TRUE.
      ELSE IF (val.EQ.".FALSE.") THEN
         pnt = .FALSE.
      ELSE
         PRINT *, "Error in SUBROUTINE assignLOGICAL:"
         PRINT *, "Invalid LOGICAL value found"
         STOP 2
      END IF
      NULLIFY(pnt)
      END SUBROUTINE assignLOGICAL

      SUBROUTINE assignDOUBLE(pnt,val)
      REAL(8), POINTER, INTENT(INOUT) :: pnt
      CHARACTER(LEN=*), INTENT(IN) :: val
      READ(val, *) pnt
      NULLIFY(pnt)
      END SUBROUTINE assignDOUBLE
      
      END MODULE pointerAssignments
      
      
      
      
      
      PROGRAM FirstTimeStep

      USE pointerAssignments

      IMPLICIT NONE
      
      INTEGER :: i, j, k, rl, nb, npw, valLen, startInd, endInd, ios
      INTEGER :: ix, iy, iz, ixMax, iyMax, izMax, pw
      INTEGER :: endOfPhotcar, fileSize, ioErr
      INTEGER :: newLineCount, lineNum, strLen
      INTEGER :: commentInd, assignmentInd, numTotalTags
      INTEGER :: updateNewLineCount, wsInd
      INTEGER, ALLOCATABLE :: newLineIndices(:), tempNewLineIndices(:)
      REAL(8) :: recordLength, ispin, singleOrDouble, numKPoints
      REAL(8) :: numBands, enMax, eFermi, numPlaneWaves
      REAL(8) :: dt
      REAL(8) :: phaseE
      REAL(8) :: omega_0
      REAL(8) :: gxGrid, gyGrid, gzGrid
      REAL(8) :: gx, gy, gz, gMax
      REAL(8), DIMENSION(3,3) :: latticeVector
      REAL(8), DIMENSION(3) :: kPts
      REAL(8), DIMENSION(2,5) :: rk
      REAL(8), ALLOCATABLE :: bandEn(:,:), gweight(:,:)
      REAL(8), ALLOCATABLE :: occ(:,:), newOcc(:,:)
      REAL(8), ALLOCATABLE :: xMom(:), yMom(:), zMom(:), totMom(:)
      COMPLEX(4), ALLOCATABLE :: pwCoef(:,:) !band,PW
      COMPLEX(8), ALLOCATABLE :: bCoef(:), new_bCoef(:), d_bCoef(:)
      COMPLEX(8), ALLOCATABLE :: d2_bCoef(:)
      COMPLEX(8), ALLOCATABLE :: px(:,:), py(:,:), pz(:,:)
      CHARACTER(LEN=:), ALLOCATABLE :: photcarStr, updateStr
      CHARACTER(LEN=:), ALLOCATABLE :: tempKey, tempVal
      CHARACTER(LEN=72) :: tempStr, dipoleFmt, dataFileFmt
      CHARACTER(LEN=72), DIMENSION(2) :: tempStrArray
      CHARACTER(LEN=:), ALLOCATABLE :: tempSubStr1, tempSubStr2
      LOGICAL :: gsUpdated=.FALSE., esUpdated=.FALSE.
      
      REAL(8), PARAMETER :: pi = 3.1415926
      COMPLEX(4), PARAMETER :: im = (0,1)
      REAL(4), PARAMETER :: hbar = 1
      REAL(4), PARAMETER :: q = -1
      REAL(4), PARAMETER :: m = 1
      
      REAL(8) :: tempDouble
      COMPLEX(8) :: tempComplex
      
      LOGICAL, TARGET :: RWA
      REAL(8), TARGET :: E, timeStep, omega
      INTEGER, TARGET :: itScheme
      INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: bandArray

      TYPE :: tagStruct
         CHARACTER(LEN=:), ALLOCATABLE :: key, val
         TYPE(allocVal), DIMENSION(:), ALLOCATABLE :: valArr
         LOGICAL, POINTER :: lpnt => NULL()
         REAL(8), POINTER :: rpnt => NULL()
         INTEGER, POINTER :: ipnt => NULL()
         INTEGER, DIMENSION(:), POINTER :: iArrPnt => NULL()
         LOGICAL :: used=.FALSE., updated=.FALSE.
      END TYPE tagStruct
      
      TYPE :: photcar
         TYPE(tagStruct), ALLOCATABLE :: tags(:)
      END TYPE photcar
      
      TYPE(photcar) :: photcarParams
      
      
      
      
      
!!!!!!!!!! Read PHOTCAR tags to set parameters !!!!!!!!!!

      OPEN(UNIT=1,
     $     FILE='PHOTCAR',
     $     STATUS='old',
     $     ACCESS='stream',
     $     ACTION='read',
     $     FORM='unformatted',
     $     IOSTAT=ioErr)

      IF (ioErr.NE.0) THEN
         WRITE(*,*) "Failed to read PHOTCAR. Terminating program."
         STOP 2
      END IF

!!!!! Import PHOTCAR as a single string !!!!!
      
      INQUIRE(UNIT=1, SIZE=fileSize)
      ALLOCATE(CHARACTER(LEN=fileSize) :: photcarStr)
      READ(UNIT=1, POS=1) photcarStr
      CLOSE(1)

!!!!! Figure out # of lines and newline indices !!!!!
      
      ALLOCATE(tempNewLineIndices(50))
      newLineCount = 0
      tempNewLineIndices(1) = 0
      DO i=1,LEN(photcarStr)
         IF (photcarStr(i:i).EQ.NEW_LINE("")) THEN
            newLineCount = newLineCount + 1
            tempNewLineIndices(newLineCount+1) = i
         END IF
      END DO
      
      tempNewLineIndices(newLineCount+2) = LEN(photcarStr) + 1
      ALLOCATE(newLineIndices(newLineCount+2))
      DO i=1,newLineCount+2
         newLineIndices(i) = tempNewLineIndices(i)
      END DO
      
      DEALLOCATE(tempNewLineIndices)
      
      
!!!!! Default parameter values !!!!!
      
      E = 0.01
      timeStep = 0.01         !time step in femptoseconds
      omega = 0.2
      RWA = .FALSE.
      itScheme = 1
      ALLOCATE(bandArray(2))
      DO i=1,SIZE(bandArray)
         bandArray(i) = i
      END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      numTotalTags = 6 !!!!!!!! UPDATE THIS W/ NEW PHOTCAR TAGS !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
!!!!! Instantiate the PHOTCAR parameter data structure !!!!!
      
      ALLOCATE(photcarParams%tags(numTotalTags))
      
!      DO i=1,numTotalTags
!         NULLIFY(photcarParams%tags(i)%rpnt)
!         NULLIFY(photcarParams%tags(i)%lpnt)
!         NULLIFY(photcarParams%tags(i)%ipnt)
!         NULLIFY(photcarParams%tags(i)%iArrPnt)
!         ALLOCATE(photcarParams%tags(i)%valArr(SIZE(bandArray)))
!      END DO
      ALLOCATE(photcarParams%tags(6)%valArr(SIZE(bandArray)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! UPDATE THESE W/ NEW PHOTCAR TAGS !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      photcarParams%tags(1)%key = "E"
      photcarParams%tags(1)%rpnt => E
      photcarParams%tags(2)%key = "timeStep"
      photcarParams%tags(2)%rpnt => timeStep
      photcarParams%tags(3)%key = "omega"
      photcarParams%tags(3)%rpnt => omega
      photcarParams%tags(4)%key = "RWA"
      photcarParams%tags(4)%lpnt => RWA
      photcarParams%tags(5)%key = "itScheme"
      photcarParams%tags(5)%ipnt => itScheme
      photcarParams%tags(6)%key = "bands"
      photcarParams%tags(6)%iArrPnt => bandArray




!!!!! Set vals to defaults !!!!!

      DO i=1,numTotalTags
         IF (ASSOCIATED(photcarParams%tags(i)%rpnt)) THEN
            WRITE(tempStr,'(E11.5)') photcarParams%tags(i)%rpnt
            photcarParams%tags(i)%val = TRIM(ADJUSTL(tempStr))
         ELSE IF (ASSOCIATED(photcarParams%tags(i)%lpnt)) THEN
            IF (photcarParams%tags(i)%lpnt.EQV..TRUE.) THEN
               tempStr = ".TRUE."
            ELSE IF (photcarParams%tags(i)%lpnt.EQV..FALSE.) THEN
               tempStr = ".FALSE."
            ELSE
               WRITE(*,*) "Found non-LOGICAL value for tag"//
     $              photcarParams%tags(i)%key
               STOP 2
            END IF
            photcarParams%tags(i)%val = TRIM(ADJUSTL(tempStr))
         ELSE IF (ASSOCIATED(photcarParams%tags(i)%ipnt)) THEN
            WRITE(tempStr,'(I1)') photcarParams%tags(i)%ipnt
            photcarParams%tags(i)%val = TRIM(ADJUSTL(tempStr))
         ELSE
            WRITE(tempStrArray,'(I10.1)') photcarParams%tags(i)%iArrPnt
            DO j=1,SIZE(tempStrArray)
               photcarParams%tags(i)%valArr(j)%val =
     $              TRIM(ADJUSTL(tempStrArray(j)))
               !WRITE(*,*) 'valArr defaults:'
               !WRITE(*,*) photcarParams%tags(i)%valArr(j)%val
            END DO
         END IF
      END DO
      
      
      
      
      
!!!!! Read PHOTCAR string line by line !!!!!
!!!!! and save contents as key/value pairs !!!!!
      DO i=1,newLineCount+1
      ! Loop line-by-line !
         startInd = newLineIndices(i)+1
         endInd = newLineIndices(i+1)-1
         commentInd = INDEX(photcarStr(startInd:endInd), "#")
         assignmentInd = INDEX(photcarStr(startInd:endInd), "=")
         IF ((assignmentInd.NE.0 .AND.
     $        assignmentInd.LT.commentInd) .OR.
     $        (assignmentInd.GT.0 .AND. commentInd.EQ.0)) THEN
            ! Break line into key and val !
            ! if "=" is before "#" or "#" isn't in the line !
            tempKey = TRIM(ADJUSTL(
     $           photcarStr(startInd:startInd+assignmentInd-2)))
            IF (commentInd.EQ.0) THEN
               ! Save val where no "#" was found !
               tempVal = TRIM(ADJUSTL(
     $              photcarStr(startInd+assignmentInd:endInd)))
            ELSE
               ! Save val where "#" was found !
               tempVal = TRIM(ADJUSTL(
     $              photcarStr(startInd+assignmentInd:
     $              startInd+commentInd-2)))
            END IF
            DO k=1,numTotalTags
               ! Update the pointer and thereby the actual variable !
               IF (tempKey.EQ.photcarParams%tags(k)%key) THEN
                  ! Check for matching key !
                  photcarParams%tags(k)%used = .TRUE.
                  wsInd = INDEX(tempVal, " ")
                  IF (wsInd.EQ.0) THEN
                     photcarParams%tags(k)%val = tempVal
                     IF (ASSOCIATED(photcarParams%tags(k)%lpnt))
     $                    THEN
                        ! Assign value through logical pointer !
                        CALL assignLOGICAL(photcarParams%tags(k)%lpnt,
     $                       photcarParams%tags(k)%val)
                        EXIT
                     ELSE IF (ASSOCIATED(photcarParams%tags(k)%rpnt))
     $                       THEN
                        ! Assign value through double pointer !
                        CALL assignDOUBLE(photcarParams%tags(k)%rpnt,
     $                       photcarParams%tags(k)%val)
                        EXIT
                     ELSE
                        ! Assign value through integer pointer !
                        CALL assignINT(photcarParams%tags(k)%ipnt,
     $                       photcarParams%tags(k)%val)
                        EXIT
                     END IF
                  ELSE
                     !WRITE(*,*) "entering array assignment"
                     photcarParams%tags(k)%valArr(1)%val =
     $                    tempVal(1:wsInd-1)
                     photcarParams%tags(k)%valArr(2)%val =
     $                    tempVal(wsInd+1:LEN(tempVal))
                     ! Assign array of vals through int array pointer !
                     CALL assignINTarray(
     $                    photcarParams%tags(k)%iArrPnt,
     $                    photcarParams%tags(k)%valArr)
                  END IF
               END IF
            END DO
         END IF
      END DO
      
      
            
      
      
!!!!! Print values of tags !!!!!
      
      WRITE(*,'(T1,A,ES23.16E3)') "E = ", E
      WRITE(*,'(T1,A,ES23.16E3)') "timeStep = ", timeStep
      WRITE(*,'(T1,A,ES23.16E3)') "omega = ", omega
      WRITE(*,'(T1,A,L1)') "RWA = ", RWA
      WRITE(*,'(T1,A,I0.1)') "itScheme = ", itScheme
      WRITE(*,'(T1,A,2I5)') "bandArray = ", bandArray
      
      dt = timeStep*41.341374575751
      
      
!!!!!!!!!! Write tags into updateWAVECAR.f !!!!!!!!!!

      OPEN(UNIT=2,
     $     FILE='updateWAVECAR.f',
     $     ACCESS='stream',
     $     STATUS='old',
     $     ACTION='read',
     $     FORM='unformatted',
     $     IOSTAT=ioErr)
      
      IF (ioErr.NE.0) THEN
         WRITE(*,*) "Failed to read updateWAVECAR.f"
         WRITE(*,*) "Terminating program."
         STOP 2
      END IF

!!!!! Import updateWAVECAR.f as a single string !!!!!
      
      INQUIRE(UNIT=2, SIZE=fileSize)
      ALLOCATE(CHARACTER(LEN=fileSize) :: updateStr)
      READ(UNIT=2, POS=1) updateStr
      
      CLOSE(2)

!!!!! Open updateWAVECAR.f to hard-code the tags into it !!!!!
      
      OPEN(UNIT=2,
     $     FILE='updateWAVECAR.f',
     $     ACCESS='sequential',
     $     STATUS='old',
     $     FORM='formatted')
      
      ios = 0
      lineNum = 0
      DO WHILE (.NOT. ALL(photcarParams%tags%updated))
         ! make sure all tags get written !
         lineNum = lineNum + 1
         READ(2,FMT='(A)',IOSTAT=ios) tempStr
         IF (ios.NE.0) THEN
            ! Terminate program if it hits end of file or something !
            WRITE(*,*) "An error has occured while updating tags in"//
     $           " updateWAVECAR.f"
            WRITE(*,*) "Terminating program."
            STOP 2
         END IF
         DO i=1,numTotalTags
            ! Look through potential tags that might be on this line !
            IF (.NOT.photcarParams%tags(i)%updated) THEN
               ! Ignore tags that have already been set !
               IF (photcarParams%tags(i)%key.EQ."bands") THEN
                  IF (INDEX(tempStr,"      bandArray(1) =").NE.0) THEN
                     BACKSPACE 2
                     WRITE(2,FMT='(A)') "      bandArray(1) = "//
     $                    photcarParams%tags(i)%valArr(1)%val
                     gsUpdated = .TRUE.
                     k = 0
                     DO j=1,fileSize
                        ! Loop through string of updateWAVECAR.f !
                        IF (updateStr(j:j).EQ.NEW_LINE("")) THEN
                           ! Find newline chars in the string !
                           k = k+1
                           IF (k.EQ.lineNum) THEN
                             ! Found current file position in string !
                             ! Rewrite the rest of file from string !
                              WRITE(2,FMT='(A)')
     $                             updateStr(j+1:fileSize)
                              ! Return to previous position in the file !
                              CLOSE(2)                           
                              OPEN(UNIT=2,
     $                             FILE='updateWAVECAR.f',
     $                             ACCESS='sequential',
     $                             STATUS='old',
     $                             FORM='formatted')
                              DO k=1,lineNum
                                 READ(2,FMT='(A)',IOSTAT=ios) tempStr
                              END DO
                              EXIT
                           END IF
                        END IF
                     END DO
                  ELSE IF (INDEX(tempStr,"      bandArray(2) =").NE.0)
     $                    THEN
                     BACKSPACE 2
                     WRITE(2,FMT='(A)') "      bandArray(2) = "//
     $                    photcarParams%tags(i)%valArr(2)%val
                     esUpdated = .TRUE.
                     k = 0
                     DO j=1,fileSize
                        ! Loop through string of updateWAVECAR.f !
                        IF (updateStr(j:j).EQ.NEW_LINE("")) THEN
                           ! Find newline chars in the string !
                           k = k+1
                           IF (k.EQ.lineNum) THEN
                              ! Found current file position in string !
                              ! Rewrite the rest of file from string !
                              WRITE(2,FMT='(A)')
     $                             updateStr(j+1:fileSize)
                              ! Return to previous position in the file !
                              CLOSE(2)                           
                              OPEN(UNIT=2,
     $                             FILE='updateWAVECAR.f',
     $                             ACCESS='sequential',
     $                             STATUS='old',
     $                             FORM='formatted')
                              DO k=1,lineNum
                                 READ(2,FMT='(A)',IOSTAT=ios) tempStr
                              END DO
                              EXIT
                           END IF
                        END IF
                     END DO
                  END IF
                  IF (gsUpdated.AND.esUpdated) THEN
                     photcarParams%tags(i)%updated = .TRUE.
                  END IF
               END IF
               strLen = LEN("      "//photcarParams%tags(i)%key//" =")
               ALLOCATE(CHARACTER(LEN=strLen) :: tempSubStr1)
               ALLOCATE(CHARACTER(LEN=strLen-1) :: tempSubStr2)
               tempSubStr1 = "      "//photcarParams%tags(i)%key//
     $              " ="
               tempSubStr2 = "      "//photcarParams%tags(i)%key//
     $              "="
               IF (INDEX(tempStr,tempSubStr1).NE.0 .OR.
     $              INDEX(tempStr,tempSubStr2).NE.0) THEN
                  ! Tag found, write in value from PHOTCAR !
                  ! This truncates the file after the written line !
                  BACKSPACE 2
                  WRITE(2,FMT='(A)') "      "//
     $                 photcarParams%tags(i)%key//" = "//
     $                 photcarParams%tags(i)%val
                  photcarParams%tags(i)%updated = .TRUE.
                  k = 0
                  DO j=1,fileSize
                     ! Loop through string of updateWAVECAR.f !
                     IF (updateStr(j:j).EQ.NEW_LINE("")) THEN
                        ! Find newline chars in the string !
                        k = k+1
                        IF (k.EQ.lineNum) THEN
                           ! Found current file position in string !
                           ! Rewrite the rest of file from string !
                           WRITE(2,FMT='(A)')
     $                          updateStr(j+1:fileSize)
                           ! Return to previous position in the file !
                           CLOSE(2)                           
                           OPEN(UNIT=2,
     $                          FILE='updateWAVECAR.f',
     $                          ACCESS='sequential',
     $                          STATUS='old',
     $                          FORM='formatted')
                           DO k=1,lineNum
                              READ(2,FMT='(A)',IOSTAT=ios) tempStr
                           END DO
                           EXIT
                        END IF
                     END IF
                  END DO
               END IF
               DEALLOCATE(tempSubStr1, tempSubStr2)
            END IF
         END DO
      END DO
      
      CLOSE(2)




      !!!!! Nullify any unused pointers !!!!!
      
      DO i=1,numTotalTags
         IF (.NOT.photcarParams%tags(i)%used) THEN
            WRITE(*,*) "Unused PHOTCAR tag: ",
     $           photcarParams%tags(i)%key
            IF (ASSOCIATED(photcarParams%tags(i)%rpnt)) THEN
               NULLIFY(photcarParams%tags(i)%rpnt)
               WRITE(*,*) "Nullified REAL(8) pointer to ",
     $              photcarParams%tags(i)%key
            END IF
            IF (ASSOCIATED(photcarParams%tags(i)%lpnt)) THEN
               NULLIFY(photcarParams%tags(i)%lpnt)
               WRITE(*,*) "Nullified LOGICAL pointer to ",
     $              photcarParams%tags(i)%key
            END IF
            IF (ASSOCIATED(photcarParams%tags(i)%ipnt)) THEN
               NULLIFY(photcarParams%tags(i)%ipnt)
               WRITE(*,*) "Nullified INTEGER pointer to ",
     $              photcarParams%tags(i)%key
            END IF
         END IF
      END DO

      
      
      
      
!!!!!!!!!! Find record length of WAVECAR!!!!!!!!!!
      
      OPEN(UNIT=21,
     $     FILE='WAVECAR_OLD',
     $     FORM='unformatted',
     $     ACCESS='direct',
     $     STATUS='old',
     $     RECL=8)

      READ(21,REC=1) recordLength
      rl = INT(recordLength)

      CLOSE(21)
      
      
      
      

!!!!!!!!!! Open WAVECAR for data extraction !!!!!!!!!!

      OPEN(UNIT=21,
     $     FILE='WAVECAR_OLD',
     $     FORM='unformatted',
     $     ACCESS='direct',
     $     STATUS='old',
     $     RECL=rl)
      



      
!!!!!!!!!! Extract metadata from WAVECAR !!!!!!!!!!
      
      READ(21,REC=1) recordLength, ispin, singleOrDouble
      READ(21,REC=2) numKPoints, numBands, enMax,
     $     ((latticeVector(i,j),i=1,3),j=1,3),
     $     eFermi
      
      

      
      
!!!!!!!!!! Extract data from WAVECAR !!!!!!!!!!

      nb = INT(numBands)
      ALLOCATE(bandEn(nb,2))
      ALLOCATE(gweight(nb,2))
      ALLOCATE(occ(nb,2))
      ALLOCATE(newOcc(nb,2))
      ALLOCATE(bCoef(nb))
      ALLOCATE(new_bCoef(nb))
      ALLOCATE(px(nb,nb))
      ALLOCATE(py(nb,nb))
      ALLOCATE(pz(nb,nb))
      ALLOCATE(d_bCoef(nb))

!!!!! Spin-up !!!!!

      READ(21,REC=3) numPlaneWaves, (kPts(i),i=1,3),
     $     (bandEn(j,1),gweight(j,1),occ(j,1),j=1,nb)

      npw = INT(numPlaneWaves)
      ALLOCATE(pwCoef(nb,npw))
      ALLOCATE(xMom(npw))
      ALLOCATE(yMom(npw))
      ALLOCATE(zMom(npw))
      ALLOCATE(totMom(npw))

      WRITE(*,'(T1,A,F18.4)') 'numPlaneWaves = ', numPlaneWaves
      WRITE(*,'(T1,A,I0.1)') 'npw = ', npw
      
      bCoef(bandArray(1)) = SQRT(occ(bandArray(1),1))
      bCoef(bandArray(2)) = SQRT(occ(bandArray(2),1))

      DO i=4,3+nb
         READ(21,REC=i) (pwCoef(i-3,j),j=1,npw)
      END DO

!!!!! Spin-down !!!!!

!# I'll worry about this later #!
      
      
      CLOSE(21)






!      OPEN(UNIT=30,
!     $     FILE='pwInfo.txt')


!      WRITE(30,*) 'size of xMom=', SIZE(xMom)
!      WRITE(30,*) 'size of yMom=', SIZE(yMom)
!      WRITE(30,*) 'size of zMom=', SIZE(zMom)
!      WRITE(30,*) 'size of totMom=', SIZE(totMom)

      



!!!!!!!!!! Calculate momentum for each plane wave !!!!!!!!!!

      !!!!! Currently assuming rectangular simulation cell !!!!!
      gxGrid = 2*pi/(latticeVector(1,1)*1.8897259886)
      gyGrid = 2*pi/(latticeVector(2,2)*1.8897259886)
      gzGrid = 2*pi/(latticeVector(3,3)*1.8897259886)
!     0.682225431
!      WRITE(30,*) 'gxGrid=', gxGrid
!      WRITE(30,*) 'gyGrid=', gyGrid
!      WRITE(30,*) 'gzGrid=', gzGrid

      gMax = DSQRT(2*(enMax*0.036749405469679))
!      WRITE(30,*) 'gMax=', gMax

      izMax = INT(gMax/gzGrid)
!      WRITE(30,*) 'izMax=', izMax
      
      pw = 0
      DO i = 0, 2*izMax
!         WRITE(30,*), 'entering z loop'
         iz = i-(2*izMax+1)*INT(i/(izMax+1)) ! z=0,1,...,zMax,-zMax,...,-1
         gz = iz*gzGrid
         iyMax = INT(DSQRT(gMax**2-gz**2)/gyGrid)
!         WRITE(30,*) 'iz=', iz, 'iyMax=', iyMax
         DO j = 0, 2*iyMax
!            WRITE(30,*), 'entering y loop'
            iy = j-(2*iyMax+1)*INT(j/(iyMax+1))
            gy = iy*gyGrid
            ixMax = INT(DSQRT(gMax**2-gz**2-gy**2)/gxGrid)
!            WRITE(30,*) 'iy=', iy, 'ixMax=', ixMax
!            DO k = 0,2*ixMax
            DO k = 0,ixMax
               ix = k
               IF (ix==0 .AND. iy<0) CYCLE
               IF (ix==0 .AND. iy==0 .AND. iz<0) CYCLE
!               WRITE(30,*) 'entering x loop'
!               WRITE(30,*) 'size of yMom=', SIZE(yMom)
!               ix = k-(2*ixMax+1)*INT(k/(ixMax+1))
!               ix = k
!               WRITE(30,*) 'calculated ix', ix
               gx = ix*gxGrid
!               WRITE(30,*) 'calculated gx', gx
               pw = pw+1
!               WRITE(30,*) 'updated pw to ', pw
               xMom(pw) = gx
!               WRITE(30,*) 'wrote to xMom(pw)', xMom(pw)
               yMom(pw) = gy
!               WRITE(30,*) 'wrote to yMom(pw)', yMom(pw)
               zMom(pw) = gz
!               WRITE(30,*) 'wrote to zMom(pw)', zMom(pw)
               totMom(pw) = DSQRT(gx**2+gy**2+gz**2)
!               WRITE(30,*) 'pw#', pw, 'gx=', gx, 'gy=', gy, 'gz=', gz
!               WRITE(30,*) 'pw#', pw, 'total momentum =', totMom(pw)
            END DO
         END DO
      END DO



      !CLOSE(30)

      
!!!!!!!!!! Calculate new orbital occupations !!!!!!!!!!

      DO i=1,nb
         d_bCoef(i) = 0
         DO j=1,nb
            px(i,j)=(0,0)
            py(i,j)=(0,0)
            pz(i,j)=(0,0)
         END DO
      END DO
      
      DO i=1,npw
         px(bandArray(1),bandArray(2)) =
     $        px(bandArray(1),bandArray(2)) +
     $        CONJG(pwCoef(bandArray(1),i))*xMom(i)*
     $        pwCoef(bandArray(2),i) -
     $        pwCoef(bandArray(1),i)*xMom(i)*
     $        CONJG(pwCoef(bandArray(2),i))
         px(bandArray(2),bandArray(1)) =
     $        px(bandArray(2),bandArray(1)) +
     $        CONJG(pwCoef(bandArray(2),i))*xMom(i)*
     $        pwCoef(bandArray(1),i) -
     $        pwCoef(bandArray(2),i)*xMom(i)*
     $        CONJG(pwCoef(bandArray(1),i))
         py(bandArray(1),bandArray(2)) =
     $        py(bandArray(1),bandArray(2)) +
     $        CONJG(pwCoef(bandArray(1),i))*yMom(i)*
     $        pwCoef(bandArray(2),i) -
     $        pwCoef(bandArray(1),i)*yMom(i)*
     $        CONJG(pwCoef(bandArray(2),i))
         py(bandArray(2),bandArray(1)) =
     $        py(bandArray(2),bandArray(1)) +
     $        CONJG(pwCoef(bandArray(2),i))*yMom(i)*
     $        pwCoef(bandArray(1),i) -
     $        pwCoef(bandArray(2),i)*yMom(i)*
     $        CONJG(pwCoef(bandArray(1),i))
         pz(bandArray(1),bandArray(2)) =
     $        pz(bandArray(1),bandArray(2)) +
     $        CONJG(pwCoef(bandArray(1),i))*zMom(i)*
     $        pwCoef(bandArray(2),i) -
     $        pwCoef(bandArray(1),i)*zMom(i)*
     $        CONJG(pwCoef(bandArray(2),i))
         pz(bandArray(2),bandArray(1)) =
     $        pz(bandArray(2),bandArray(1)) +
     $        CONJG(pwCoef(bandArray(2),i))*zMom(i)*
     $        pwCoef(bandArray(1),i) -
     $        pwCoef(bandArray(2),i)*zMom(i)*
     $        CONJG(pwCoef(bandArray(1),i))
      END DO
      
      dipoleFmt =
     $     '(T1,A,I0.1,A,I0.1,A,"( ",ES23.15E3," , ",ES23.15E3," )")'
      WRITE(*,dipoleFmt)
     $     'px(',bandArray(1),',',bandArray(2),') = ',
     $     px(bandArray(1),bandArray(2))
      WRITE(*,dipoleFmt)
     $     'px(',bandArray(2),',',bandArray(1),') = ',
     $     px(bandArray(2),bandArray(1))
      WRITE(*,dipoleFmt) 'py(',bandArray(1),',',bandArray(2),') = ',
     $     py(bandArray(1),bandArray(2))
      WRITE(*,dipoleFmt) 'py(',bandArray(2),',',bandArray(1),') = ',
     $     py(bandArray(2),bandArray(1))
      WRITE(*,dipoleFmt) 'pz(',bandArray(1),',',bandArray(2),') = ',
     $     pz(bandArray(1),bandArray(2))
      WRITE(*,dipoleFmt) 'pz(',bandArray(2),',',bandArray(1),') = ',
     $     pz(bandArray(2),bandArray(1))
      
      omega_0 = 0.0367494*
     $     (bandEn(bandArray(2),1)-bandEn(bandArray(1),1))/hbar !eV to au
      WRITE(*,'(T1,A,F18.14)') 'omega_0 = ', omega_0

      phaseE = 1.0

      d_bCoef(bandArray(1)) = 0
      IF (RWA) THEN
         d_bCoef(bandArray(2)) = (q*E/(2*m*hbar*omega_0))*
     $        pz(bandArray(2),bandArray(1))
      ELSE
         d_bCoef(bandArray(2)) = (q*E/(m*hbar*omega_0))*
     $        pz(bandArray(2),bandArray(1))
      END IF
      WRITE(*,'(T1,A,I0.1,A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $     'd_bCoef(',bandArray(2),') = ', d_bCoef(bandArray(2))
      
      
      
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!! Iteration Scheme !!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      IF (itScheme.EQ.1) THEN   ! Euler
!         new_bCoef(1) = bCoef(1) + dt*d_bCoef(1)
         new_bCoef(bandArray(2)) = bCoef(bandArray(2)) +
     $        dt*d_bCoef(bandArray(2))
         new_bCoef(bandArray(1)) = SQRT((1,0) -
     $        CONJG(new_bCoef(bandArray(2)))*new_bCoef(bandArray(2)))
      ELSE IF (itScheme.EQ.2) THEN ! Verlet
         ALLOCATE(d2_bCoef(nb))
         IF (RWA) THEN
            PRINT *, "Solving for d^2/dt under RWA"
            d2_bCoef(bandArray(1)) = -((q**2)*(E**2)*
     $           pz(bandArray(1),bandArray(2))*
     $           pz(bandArray(2),bandArray(1)))/
     $           (4*m**2*hbar**2*omega_0**2)
            d2_bCoef(bandArray(2)) = (q*E*i*(omega-omega_0)*
     $           pz(bandArray(2),bandArray(1)))/
     $           (2*m*hbar*omega_0)
         ELSE
            PRINT *, "Solving for d^2/dt without RWA"
            d2_bCoef(bandArray(1)) = -((q**2)*(E**2)*
     $           pz(bandArray(1),bandArray(2))*
     $           pz(bandArray(2),bandArray(1)))/
     $           ((m**2)*(hbar**2)*(omega_0**2))
            d2_bCoef(bandArray(2)) = (q*E*im)/(m*hbar)
         END IF
         PRINT *, 'd_bCoef(',bandArray(2),') = ', d_bCoef(bandArray(2))
         PRINT *, "dt = ", dt
         PRINT *, 'd2_bCoef(',bandArray(1),') = ',
     $        d2_bCoef(bandArray(1))
         PRINT *, 'd2_bCoef(',bandArray(2),') = ',
     $        d2_bCoef(bandArray(2))
         tempComplex = d_bCoef(bandArray(2))*dt
         PRINT *, 'd_bCoef(',bandArray(2),')*dt = ', tempComplex
         tempComplex = 0.5*d2_bCoef(bandArray(2))*dt**2
         PRINT *, '0.5*d2_bCoef(',bandArray(2),')*dt**2 = ',
     $        tempComplex
         new_bCoef(bandArray(1)) = bCoef(bandArray(1)) +
     $        d_bCoef(bandArray(1))*dt +
     $        0.5*d2_bCoef(bandArray(1))*dt**2
         new_bCoef(bandArray(2)) = bCoef(bandArray(2)) +
     $        d_bCoef(bandArray(2))*dt +
     $        0.5*d2_bCoef(bandArray(2))*dt**2
         PRINT *, "Current end of Verlet code block"
      ELSE IF (itScheme.EQ.3) THEN ! Runge-Kutta 2nd Order
         rk(1,1) = dt*d_bCoef(bandArray(1))
         rk(2,1) = dt*d_bCoef(bandArray(2))
         IF (RWA) THEN
            rk(1,2) = dt*((-q*E*EXP(im*(omega-omega_0)*dt)*
     $           pz(bandArray(1),bandArray(2))*(rk(2,1)))
     $           /(2*m*hbar*omega_0))
            rk(2,2) = dt*((q*E*EXP(im*(omega_0-omega)*dt)*
     $           pz(bandArray(2),bandArray(1))*(1+rk(1,1)))/
     $           (2*m*hbar*omega_0))
         ELSE
            rk(1,2) = dt*((-q*E*DCOS(omega*dt)*
     $           pz(bandArray(1),bandArray(2))*
     $           EXP(-im*omega_0*dt)*rk(2,1))/
     $           (m*hbar*omega_0))
            rk(2,2) = dt*((q*E*DCOS(omega*dt)*
     $           pz(bandArray(2),bandArray(1))*
     $           EXP(im*omega_0*dt)*(1+rk(1,1)))/
     $           (m*hbar*omega_0))
         END IF
         rk(1,5) = (rk(1,1)+rk(1,2))/2
         rk(2,5) = (rk(2,1)+rk(2,2))/2
         new_bCoef(bandArray(1)) = bCoef(bandArray(1)) + rk(1,5)
         new_bCoef(bandArray(2)) = bCoef(bandArray(2)) + rk(2,5)
         PRINT *, "End of Runge-Kutta 2 block"
      ELSE IF (itScheme.EQ.4) THEN ! Runge-Kutta 4th Order
         rk(1,1) = dt*d_bCoef(bandArray(1))
         rk(2,1) = dt*d_bCoef(bandArray(2))
         IF (RWA) THEN
            rk(1,2) = dt*((-q*E*EXP(im*(omega-omega_0)*dt/2)*
     $           pz(bandArray(1),bandArray(2))*(rk(2,1)/2))
     $           /(2*m*hbar*omega_0))
            rk(2,2) = dt*((q*E*EXP(im*(omega_0-omega)*dt/2)*
     $           pz(bandArray(2),bandArray(1))*(1+rk(1,1)/2))/
     $           (2*m*hbar*omega_0))
            rk(1,3) = dt*((-q*E*EXP(im*(omega-omega_0)*dt/2)*
     $           pz(bandArray(1),bandArray(2))*(rk(2,2)/2))
     $           /(2*m*hbar*omega_0))
            rk(2,3) = dt*((q*E*EXP(im*(omega_0-omega)*dt/2)*
     $           pz(bandArray(2),bandArray(1))*(1+rk(1,2)/2))/
     $           (2*m*hbar*omega_0))
            rk(1,4) = dt*((-q*E*EXP(im*(omega-omega_0)*dt)*
     $           pz(bandArray(1),bandArray(2))*rk(2,3))/
     $           (2*m*hbar*omega_0))
            rk(2,4) = dt*((q*E*EXP(im*(omega_0-omega)*dt)*
     $           pz(bandArray(2),bandArray(1))*(1+rk(1,3)))/
     $           (2*m*hbar*omega_0))
         ELSE
            rk(1,2) = dt*((-q*E*DCOS(omega*dt/2)*
     $           pz(bandArray(1),bandArray(2))*
     $           EXP(-im*omega_0*dt/2)*(rk(2,1)/2))/
     $           (m*hbar*omega_0))
            rk(2,2) = dt*((q*E*DCOS(omega*dt/2)*
     $           pz(bandArray(2),bandArray(1))*
     $           EXP(im*omega_0*dt/2)*(1+rk(1,1)/2))/
     $           (m*hbar*omega_0))
            rk(1,3) = dt*((-q*E*DCOS(omega*dt/2)*
     $           pz(bandArray(1),bandArray(2))*
     $           EXP(-im*omega_0*dt/2)*(rk(2,2)/2))/
     $           (m*hbar*omega_0))
            rk(2,3) = dt*((q*E*DCOS(omega*dt/2)*
     $           pz(bandArray(2),bandArray(1))*
     $           EXP(im*omega_0*dt/2)*(1+rk(1,2)/2))/
     $           (m*hbar*omega_0))
            rk(1,4) = dt*((-q*E*DCOS(omega*dt)*
     $           pz(bandArray(1),bandArray(2))*
     $           EXP(-im*omega_0*dt)*rk(2,3))/
     $           (m*hbar*omega_0))
            rk(2,4) = dt*((q*E*DCOS(omega*dt)*
     $           pz(bandArray(2),bandArray(1))*
     $           EXP(im*omega_0*dt)*(1+rk(1,3)))/
     $           (m*hbar*omega_0))
         END IF
         rk(1,5) = (rk(1,1)+2*rk(1,2)+2*rk(1,3)+rk(1,4))/6
         rk(2,5) = (rk(2,1)+2*rk(2,2)+2*rk(2,3)+rk(2,4))/6
         new_bCoef(bandArray(1)) = bCoef(bandArray(1)) + rk(1,5)
         new_bCoef(bandArray(2)) = bCoef(bandArray(2)) + rk(2,5)
         PRINT *, "Reached end of Runge-Kutta 4 block"
      ELSE
         PRINT *, "Invalid iteration scheme provided"
         PRINT *, "Terminating program"
         STOP 2
      END IF

      WRITE(*,'(T1,A,I0.1,A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $     'new_bCoef(',bandArray(1),') = ', new_bCoef(bandArray(1))
      WRITE(*,'(T1,A,I0.1,A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $     'new_bCoef(',bandArray(2),') = ', new_bCoef(bandArray(2))
      WRITE(*,'(A)') "which means..."
      WRITE(*,'(T1,A,I0.1,A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $     'newOcc(',bandArray(1),') = ',
     $     CONJG(new_bCoef(bandArray(1)))*new_bCoef(bandArray(1))
      WRITE(*,'(T1,A,I0.1,A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $     'newOcc(',bandArray(2),') = ',
     $     CONJG(new_bCoef(bandArray(2)))*new_bCoef(bandArray(2))
      WRITE(*,'(T1,A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $     'total population = ',
     $     CONJG(new_bCoef(bandArray(1)))*new_bCoef(bandArray(1)) +
     $     CONJG(new_bCoef(bandArray(2)))*new_bCoef(bandArray(2))
      WRITE(*,'(T1,A)') "Reached end of time iteration loop"

      
      
      
!newOcc(1,1) = CONJG(new_bCoef(1))*new_bCoef(1)
      newOcc(bandArray(2),1) = CONJG(new_bCoef(bandArray(2)))*
     $     new_bCoef(bandArray(2))
      newOcc(bandArray(1),1) = CONJG(new_bCoef(bandArray(1)))*
     $     new_bCoef(bandArray(1))
      WRITE(*,'(T1,A,ES23.15E3,A,ES23.15E3)')
     $     'new Occ.s=', newOcc(bandArray(1),1), '     ',
     $     newOcc(bandArray(2),1)
      
      
      
      
!!!!!!!!!! Update WAVECAR occupations !!!!!!!!!!

      OPEN(UNIT=12,
     $     FILE='WAVECAR',
     $     FORM='unformatted',
     $     ACCESS='direct',
     $     STATUS='old',
     $     RECL=8)

      WRITE(*,'(A)') 'Opened WAVECAR'
      !# UPDATE THIS #!
      WRITE(12,REC=rl/4+4+3*bandArray(1)) newOcc(bandArray(1),1)
      WRITE(12,REC=rl/4+4+3*bandArray(2)) newOcc(bandArray(2),1)

      CLOSE(12)
      
      WRITE(*,'(A)') 'Updated WAVECAR'




!!!!!!!!!! Create data file !!!!!!!!!!

      OPEN(UNIT=50,
     $     FILE='dataFile.txt')

      dataFileFmt = '(A,6(A,I0.1,A),6(A,I0.1,A,I0.1,A),A)'
      WRITE(UNIT=50, FMT=dataFileFmt)
     $     'time     ',
     $     'new_bCoef(',bandArray(1),')     ',
     $     'new_bCoef(',bandArray(2),')     ',
     $     'bandEn(',bandArray(1),',1)     ',
     $     'bandEn(',bandArray(2),',1)     ',
     $     'newOcc(',bandArray(1),',1)     ',
     $     'newOcc(',bandArray(2),',1)     ',
     $     'px(',bandArray(1),',',bandArray(2),')     ',
     $     'px(',bandArray(2),',',bandArray(1),')     ',
     $     'py(',bandArray(1),',',bandArray(2),')     ',
     $     'py(',bandArray(2),',',bandArray(1),')     ',
     $     'pz(',bandArray(1),',',bandArray(2),')     ',
     $     'pz(',bandArray(2),',',bandArray(1),')     ',
     $     'phaseE     '
      
      WRITE(UNIT=50, FMT=*) timeStep,
     $     new_bCoef(bandArray(1)), new_bCoef(bandArray(2)),
     $     bandEn(bandArray(1),1), bandEn(bandArray(2),1),
     $     newOcc(bandArray(1),1), newOcc(bandArray(2),1),
     $     px(bandArray(1),bandArray(2)),
     $     px(bandArray(2),bandArray(1)),
     $     py(bandArray(1),bandArray(2)),
     $     py(bandArray(2),bandArray(1)),
     $     pz(bandArray(1),bandArray(2)),
     $     pz(bandArray(2),bandArray(1)),
     $     phaseE
      
      CLOSE(50)
      
      END PROGRAM FirstTimeStep
