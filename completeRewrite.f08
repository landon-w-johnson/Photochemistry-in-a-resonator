PROGRAM UnnamedJank
  
  IMPLICIT NONE

  ! Parameters !
  
  REAL(8), PARAMETER :: hbar = 1
  REAL(8), PARAMETER :: q = -1
  REAL(8), PARAMETER :: m = 1
  COMPLEX(8), PARAMETER :: im = (0,1)
  REAL(8), PARAMETER :: pi = 3.141592653589793
  REAL(8), PARAMETER :: fs_to_au = 41.34137333517
  REAL(8), PARAMETER :: au_to_fs = 0.024188843265864
  REAL(8), PARAMETER :: angstroms_to_au = 1.8897261259
  REAL(8), PARAMETER :: eV_to_au = 0.036749322175654


  
  ! PHOTCAR tags !
  
  REAL(8), TARGET :: EfieldAmp, timeStep, omega, excitationE
  LOGICAL, TARGET :: RWA
  INTEGER, TARGET :: intScheme, CHGint, maxSteps
  CHARACTER(LEN=:), ALLOCATABLE, TARGET :: VASPcmd


  
  ! INCAR tags !
  
  REAL(8) :: ispin_INCAR = 1
  LOGICAL :: ncl = .FALSE.
  INTEGER :: isym


  
  ! WAVECAR metadata !
  
  REAL(8) :: recordLength, ispin_WAVECAR, singleOrDouble
  REAL(8) :: numKpoints, numBands, enMax, eFermi, numPlaneWaves
  REAL(8), DIMENSION(3,3) :: latticeVector
  INTEGER :: rl, nb, nkp, npw



  ! As-of-yet uncategorized variables !
  
  INTEGER :: i, j, ioStatus, ispin, spinDownRec
  INTEGER :: izMax, ix, iy, iz
  REAL(8) :: dt
  REAL(8) :: gxGrid, gyGrid, gzGrid, gMax
  REAL(8), ALLOCATABLE :: bandEn(:,:), gweight(:,:), occ(:,:), newOcc(:,:)
  REAL(8), ALLOCATABLE :: bCoef(:,:), new_bCoef(:,:), d_bCoef(:,:)
  REAL(8), ALLOCATABLE :: px(:,:,:), py(:,:,:), pz(:,:,:)
  REAL(8), ALLOCATABLE :: xMom(:), yMom(:), zMom(:), totMom(:)
  REAL(8), ALLOCATABLE :: kPts(:,:), pwCoef(:,:,:)
  
  
  
  
  
  TYPE :: tagStruct
     CHARACTER(LEN=:), ALLOCATABLE :: key, val
     LOGICAL, POINTER :: lpnt => NULL()
     REAL(8), POINTER :: rpnt => NULL()
     INTEGER, POINTER :: ipnt => NULL()
     CHARACTER(LEN=:), POINTER :: spnt => NULL()
     LOGICAL :: used = .FALSE., updated = .FALSE.
  END TYPE tagStruct

  TYPE :: photcarTags
     TYPE(tagStruct), ALLOCATABLE :: tags (:)
  END TYPE photcarTags


  
!!!!!!!!!!!!!!!!!!!! BEGIN MAIN PROGRAM !!!!!!!!!!!!!!!!!!!!



  ! Parse PHOTCAR tags and set the corresponding variables!
  
  CALL parsePHOTCAR(EfieldAmp, timeStep, omega, RWA, intScheme, excitationE, CHGint, maxSteps, VASPcmd)

  WRITE(*,'(A)') 'Using the following values for PHOTCAR tags:'
  WRITE(*,'(T1,A,ES23.16E3)') '  EfieldAmp = ', EfieldAmp
  WRITE(*,'(T1,A,ES23.16E3)') '  timeStep = ', timeStep
  WRITE(*,'(T1,A,ES23.16E3)') '  omega = ', omega
  WRITE(*,'(T1,A,L1)') '  RWA = ', RWA
  WRITE(*,'(T1,A,I0.1)') '  intScheme = ', intScheme
  WRITE(*,'(T1,A,ES23.16E3)') '  excitationE = ', excitationE
  WRITE(*,'(T1,A,I0.1)') '  CHGint = ', CHGint
  WRITE(*,'(T1,A,I0.1)') '  maxSteps = ', maxSteps
  WRITE(*,'(T1,A,A)') '  VASPcmd = ', TRIM(VASPcmd)
  WRITE(*,*) ''

  dt = timeStep*fs_to_au

  

  ! Extract info from WAVECAR !
  
  rl = getRecordLength()
  WRITE(*,'(A,I0.1)') 'record length of WAVECAR_OLD = ', rl
  WRITE(*,*) ''
  
  CALL readMetadata(rl, recordLength, ispin_WAVECAR, singleOrDouble,&
       numKpoints, nkp, numBands, nb, enMax, latticeVector, eFermi)

  WRITE(*,'(A,I0.1)') 'number of bands in WAVECAR_OLD = ', nb
  WRITE(*,'(A,I0.1)') 'number of k points in WAVECAR_OLD = ', nkp
  WRITE(*,*) ''


  
  ! Determine how spins are being handled in INCAR !
  
  CALL parseINCAR(ncl, ispin_INCAR, isym)
  WRITE(*,'(A,L1)') 'using noncollinear spins: ', ncl
  WRITE(*,'(A,F0.1)') 'ISPIN from INCAR: ', ispin_INCAR
  WRITE(*,'(A,I0.1)') 'ISYM from INCAR: ', isym
  WRITE(*,*) ''

  IF (.NOT. ncl .AND. ispin_WAVECAR /= ispin_INCAR) THEN
     WRITE(*,'(A)') 'Error in spin handling!'
     WRITE(*,'(A)') 'Not using noncollinear spins and'
     WRITE(*,'(A,F0.1,A,F0.1,A)') 'ISPIN from WAVECAR_OLD (', ispin_WAVECAR,&
          ') does not match ISPIN from INCAR (', ispin_INCAR, ')!'
     STOP 2
  END IF

  
  
  ! Diverge behavior based on spins !

  IF (ncl) THEN
     WRITE(*,'(A)') 'Noncollinear spins not implemented yet! Terminating program!'
     STOP 2
  ELSE ! not using noncollinear spins
     ispin = INT(ispin_INCAR)
     SELECT CASE (ispin)
     CASE (1) ! closed-shell orbitals
        WRITE(*,'(A)') 'Closed-shell not implemented yet! Terminating program!'
        STOP 2
     CASE (2) ! polarized spins

        ! Instantiate necessary data structures !
        
        ALLOCATE(bandEn(nb,2))
        ALLOCATE(gweight(nb,2))
        ALLOCATE(occ(nb,2))
        ALLOCATE(newOcc(nb,2))
        ALLOCATE(bCoef(nb,2))
        ALLOCATE(new_bCoef(nb,2))
        ALLOCATE(px(nb,nb,2))
        ALLOCATE(py(nb,nb,2))
        ALLOCATE(pz(nb,nb,2))
        ALLOCATE(d_bCoef(nb,2))
        ALLOCATE(kPts(nkp,3))

        ! Get the rest of the data out of WAVECAR_OLD !
        OPEN(UNIT=21,&
         FILE='WAVECAR_OLD',&
         FORM='unformatted',&
         ACCESS='direct',&
         STATUS='old',&
         RECL=rl)
        
!!!!!!!!!! REMEMBER TO UPDATE THIS FOR K POINTS !!!!!!!!!!
        ! Spin up !
        READ(21,REC=3) numPlaneWaves, (kPts(1,i), i=1,3),&
             (bandEn(j,1), gweight(j,1), occ(j,1), j=1,nb)
        npw = INT(numPlaneWaves)
        ! Spin down !
        spinDownRec = 4 + nb
        READ(21,REC=spinDownRec) numPlaneWaves, (kPts(1,i), i=1,3),&
             (bandEn(j,2), gweight(j,2), occ(j,2), j=1,nb)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CLOSE(21)

        WRITE(*,*) 'kPts:'
        WRITE(*,*) kPts
        WRITE(*,*) ''
        WRITE(*,*) 'bandEn:'
        WRITE(*,*) bandEn
        WRITE(*,*) ''
        WRITE(*,*) 'occ:'
        WRITE(*,*) occ
        WRITE(*,*) ''
        WRITE(*,'(A,F18.4)') 'numPlaneWaves = ', numPlaneWaves
        WRITE(*,'(A,I0.1)') 'npw = ', npw
        WRITE(*,*) ''

        ALLOCATE(pwCoef(nb,npw,2))
        ALLOCATE(xMom(npw))
        ALLOCATE(yMom(npw))
        ALLOCATE(zMom(npw))
        ALLOCATE(totMom(npw))
        

        
        ! Calculate momentum for each plane wave !
        
        CALL calcMom(enMax, latticeVector, gxGrid, gyGrid, gzGrid, gMax, izMax,&
             isym, xMom, yMom, zMom, totMom)
        
        
        
     CASE DEFAULT ! undefined ISPIN
        WRITE(*,'(A,I0.1)') 'Error! Undefined ISPIN: ', ispin
        WRITE(*,'(A)') 'Terminating Program!'
        STOP 2
     END SELECT
  END IF
  
  








  




  



CONTAINS




  
  SUBROUTINE assignLOGICAL(pnt,val)
    IMPLICIT NONE
    
    LOGICAL, POINTER, INTENT(INOUT) :: pnt
    CHARACTER(LEN=*), INTENT(IN) :: val
    
    SELECT CASE(TRIM(ADJUSTL(val)))
    CASE('.TRUE.', 'T', 'TRUE')
       pnt = .TRUE.
    CASE('.FALSE.', 'F', 'FALSE')
       pnt = .FALSE.
    CASE DEFAULT
       WRITE(*,*) 'Error in SUBROUTINE assignLOGICAL.'
       WRITE(*,*) 'Invalid LOGICAL value: ', TRIM(val)
       STOP 2
    END SELECT
  END SUBROUTINE assignLOGICAL


  
  
  
  SUBROUTINE assignDOUBLE(pnt, val)
    IMPLICIT NONE
    
    REAL(8), POINTER, INTENT(INOUT) :: pnt
    CHARACTER(LEN=*), INTENT(IN) :: val
    
    READ(val, *, IOSTAT=ioStatus) pnt
    IF (ioStatus /= 0) THEN
       WRITE(*,*) 'Error in SUBROUTINE assignDOUBLE.'
       WRITE(*,*) 'Invaled REAL value: ', TRIM(val)
       STOP 2
    END IF
  END SUBROUTINE assignDOUBLE


  


  SUBROUTINE assignINT(pnt, val)
    IMPLICIT NONE
    
    INTEGER, INTENT(INOUT) :: pnt
    CHARACTER(LEN=*), INTENT(IN) :: val

    READ(val, *, IOSTAT=ioStatus) pnt
    IF (ioStatus /= 0) THEN
       WRITE(*,*) 'Error in SUBROUTINE assignINT.'
       WRITE(*,*) 'Invalid INTEGER value: ', TRIM(val)
       STOP 2
    END IF
  END SUBROUTINE assignINT



  

  SUBROUTINE assignSTRING(pnt, val)
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(INOUT) :: pnt
    CHARACTER(LEN=*), INTENT(IN) :: val

    READ(val, *, IOSTAT=ioStatus) pnt
    IF (ioStatus /= 0) THEN
       WRITE(*,*) 'Error in SUBROUTINE assignSTRING.'
       WRITE(*,*) 'Invalid STRING value: ', TRIM(val)
       STOP 2
    END IF
  END SUBROUTINE assignSTRING



  

  SUBROUTINE parsePHOTCAR(EfieldAmp, timeStep, omega, RWA, intScheme, excitationE, CHGint, maxSteps, VASPcmd)
    IMPLICIT NONE
    
    REAL(8), INTENT(OUT), TARGET :: EfieldAmp, timeStep, omega, excitationE
    LOGICAL, INTENT(OUT), TARGET :: RWA
    INTEGER, INTENT(OUT), TARGET :: intScheme, CHGint, maxSteps
    CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT), TARGET :: VASPcmd
    
    INTEGER :: ioStatus, fileSize, newLineCount, i, j, k, numTotalTags
    INTEGER :: startInd, endInd, commentInd, assignmentInd
    CHARACTER(LEN=:), ALLOCATABLE :: photcarStr, tmpKey, tmpVal
    CHARACTER(72) :: tmpStr
    !INTEGER, ALLOCATABLE :: newLineIndices(:), tmpNewLineIndices(:)
    INTEGER, DIMENSION(100) :: tmpNewLineIndices
    INTEGER, ALLOCATABLE :: newLineIndices(:)

    TYPE(photcarTags) :: photcarParams

    ! Read PHOTCAR !
    OPEN(UNIT=1,&
         FILE='PHOTCAR',&
         STATUS='old',&
         ACCESS='stream',&
         ACTION='read',&
         FORM='unformatted',&
         IOSTAT=ioStatus)
    IF (ioStatus /= 0) THEN
       WRITE(*,*) 'Failed to read PHOTCAR. Terminating program.'
       STOP 2
    END IF
    INQUIRE(UNIT=1, SIZE=fileSize)
    ALLOCATE(CHARACTER(LEN=fileSize) :: photcarStr)
    READ(UNIT=1, POS=1) photcarStr
    CLOSE(1)

    ! Find number and indices of newline characters (and beginning and end of file)!
    !ALLOCATE(tmpNewLineIndices(50))
    newLineCount = 0
    tmpNewLineIndices(1) = 0
    DO i = 1, LEN(photcarStr)
       IF (photcarStr(i:i) == NEW_LINE("")) THEN
          newLineCount = newLineCount + 1
          tmpNewLineIndices(newLineCount+1) = i
       END IF
    END DO
    tmpNewLineIndices(newLineCount+2) = LEN(photcarStr) + 1

    ! Save newline indices into an appropriately-sized array !
    ALLOCATE(newLineIndices(newLineCount+2))
    DO i = 1, newLineCount+2
       newLineIndices(i) = tmpNewLineIndices(i)
    END DO
    !DEALLOCATE(tmpNewLineIndices)

    ! Set default values for missing PHOTCAR tags !
    EfieldAmp = 0.01 ! maximum electric field amplitude in atomic units
    timeStep = 0.01 ! time step in femptoseconds
    omega = 0.2 ! angular frequency of E-field amplitude in Hz*radians
    RWA = .FALSE. ! whether or not to use the rotating wave approximation
    intScheme = 1 ! integration scheme, e.g. 1=Euler, 2=Verlet, 3=RK2, or 4=RK4
    excitationE = 1 ! energy of the photons in atomic units
    CHGint = 0 ! time steps between saved CHG files. 0=don't save any
    maxSteps = 0 ! maximum number of time steps to calculate. 0=no maximum
    VASPcmd = 'vasp_gam' ! VASP executable to run, e.g. 'vasp_gam', 'vasp_std', or 'vasp_ncl'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    numTotalTags = 9 ! UPDATE THIS W/ NEW PHOTCAR TAGS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Instantiate photcarParams key/val pairs !
    ALLOCATE(photcarParams%tags(numTotalTags))

    ! Set keys for PHOTCAR parameters and associate the pointers to the values
    photcarParams%tags(1)%key = 'EfieldAmp'
    photcarParams%tags(1)%rpnt => EfieldAmp
    photcarParams%tags(2)%key = 'timeStep'
    photcarParams%tags(2)%rpnt => timeStep
    photcarParams%tags(3)%key = 'omega'
    photcarParams%tags(3)%rpnt => omega
    photcarParams%tags(4)%key = 'RWA'
    photcarParams%tags(4)%lpnt => RWA
    photcarParams%tags(5)%key = 'intScheme'
    photcarParams%tags(5)%ipnt => intScheme
    photcarParams%tags(6)%key = 'excitationE'
    photcarParams%tags(6)%rpnt => excitationE
    photcarParams%tags(7)%key = 'CHGint'
    photcarParams%tags(7)%ipnt => CHGint
    photcarParams%tags(8)%key = 'maxSteps'
    photcarParams%tags(8)%ipnt => maxSteps
    photcarParams%tags(9)%key = 'VASPcmd'
    photcarParams%tags(9)%spnt => VASPcmd

    ! Set vals to defaults !
    DO i = 1, numTotalTags
       IF (ASSOCIATED(photcarParams%tags(i)%rpnt)) THEN
          WRITE(tmpStr,'(E11.5)') photcarParams%tags(i)%rpnt
          photcarParams%tags(i)%val = TRIM(ADJUSTL(tmpStr))
       ELSE IF (ASSOCIATED(photcarParams%tags(i)%lpnt)) THEN
          IF (photcarParams%tags(i)%lpnt.EQV..TRUE.) THEN
             tmpStr = '.TRUE.'
          ELSE IF (photcarParams%tags(i)%lpnt.EQV..FALSE.) THEN
             tmpStr = '.FALSE.'
          ELSE
             WRITE(*,*) 'Found non-logical value for tag '//&
                  photcarParams%tags(i)%key
             STOP 2
          END IF
          photcarParams%tags(i)%val = TRIM(ADJUSTL(tmpStr))
       ELSE IF (ASSOCIATED(photcarParams%tags(i)%ipnt)) THEN
          WRITE(tmpStr,'(I1)') photcarParams%tags(i)%ipnt
          photcarParams%tags(i)%val = TRIM(ADJUSTL(tmpStr))
       ELSE ! string
          WRITE(tmpStr,'(A)') photcarParams%tags(i)%spnt
          photcarParams%tags(i)%val = TRIM(ADJUSTL(tmpStr))
       END IF
    END DO
    
    ! Read PHOTCAR string line by line and save contents as key/val pairs !
    DO i = 1, newLineCount
       startInd = newLineIndices(i)+1
       endInd = newLineIndices(i+1)-1
       commentInd = INDEX(photcarStr(startInd:endInd), '#')
       assignmentInd = INDEX(photcarStr(startInd:endInd), '=')
       ! If '=' is there and before '#'
       IF ((assignmentInd /= 0 .AND. assignmentInd < commentInd) .OR. &
            ! Or if '=' is there and '#' is not
            (assignmentInd < 0 .AND. commentInd == 0)) THEN
          ! Extract key/val pair
          tmpKey = TRIM(ADJUSTL(photcarStr(startInd:startInd+assignmentInd-2)))
          IF (commentInd == 0) THEN
             ! Save everything to the right of '=' as the val
             tmpVal = TRIM(ADJUSTL(photcarStr(startIND+assignmentInd:endInd)))
          ELSE
             ! Save everything between '=' and '#' as the val
             tmpVal = TRIM(ADJUSTL(photcarStr(&
                  startInd+assignmentInd:startInd+commentInd-2)))
          END IF
          DO k = 1, numTotalTags
             ! Find the matching tag and update its pointer
             IF (tmpKey == photcarParams%tags(k)%key) THEN
                photcarParams%tags(k)%used = .TRUE.
                photcarParams%tags(k)%val = tmpVal
                IF (ASSOCIATED(photcarParams%tags(k)%lpnt)) THEN
                   CALL assignLOGICAL(&
                        photcarParams%tags(k)%lpnt,&
                        photcarParams%tags(k)%val)
                   EXIT
                ELSE IF (ASSOCIATED(photcarParams%tags(k)%rpnt)) THEN
                   CALL assignDOUBLE(&
                        photcarParams%tags(k)%rpnt,&
                        photcarParams%tags(k)%val)
                   EXIT
                ELSE IF (ASSOCIATED(photcarParams%tags(k)%ipnt)) THEN
                   CALL assignINT(&
                        photcarParams%tags(k)%ipnt,&
                        photcarParams%tags(k)%val)
                   EXIT
                ELSE ! string
                   CALL assignSTRING(&
                        photcarParams%tags(k)%spnt,&
                        photcarParams%tags(k)%val)
                   EXIT
                END IF
             END IF
          END DO
       END IF
    END DO

    ! Nullify any unused pointers !
    DO i = 1, numTotalTags
       IF (.NOT. photcarParams%tags(i)%used) THEN
          WRITE(*,'(A,A)') 'Unused PHOTCAR tag: ', TRIM(photcarParams%tags(i)%key)
          IF (ASSOCIATED(photcarParams%tags(i)%rpnt)) THEN
             NULLIFY(photcarParams%tags(i)%rpnt)
             WRITE(*,'(A,A)') 'Nullified REAL(8) POINTER to ',&
                  TRIM(photcarParams%tags(i)%key)
          END IF
          IF (ASSOCIATED(photcarParams%tags(i)%lpnt)) THEN
             NULLIFY(photcarParams%tags(i)%lpnt)
             WRITE(*,'(A,A)') 'Nullified LOGICAL POINTER to ',&
                  photcarParams%tags(i)%key
          END IF
          IF (ASSOCIATED(photcarParams%tags(i)%ipnt)) THEN
             NULLIFY(photcarParams%tags(i)%ipnt)
             WRITE(*,'(A,A)') 'Nullified INTEGER POINTER to ',&
                  photcarParams%tags(i)%key
          END IF
          IF (ASSOCIATED(photcarParams%tags(i)%spnt)) THEN
             NULLIFY(photcarParams%tags(i)%spnt)
             WRITE(*,'(A,A)') 'Nullified STRING POINTER to ',&
                  photcarParams%tags(i)%key
          END IF
          WRITE(*,*) ''
       END IF
    END DO
  END SUBROUTINE parsePHOTCAR


  


  INTEGER FUNCTION getRecordLength()
    IMPLICIT NONE
    
    REAL(8) :: recordLength
    LOGICAL :: fileExists
    
    INQUIRE(FILE='WAVECAR_OLD', EXIST=fileExists)
    IF (.NOT. fileExists) THEN
       CALL SYSTEM('cp WAVECAR WAVECAR_OLD')
    END IF

    OPEN(UNIT=21,&
         FILE='WAVECAR_OLD',&
         FORM='unformatted',&
         ACCESS='direct',&
         STATUS='old',&
         RECL=8)
    READ(21,REC=1) recordLength
    CLOSE(21)
    getRecordLength = INT(recordLength)
  END FUNCTION getRecordLength


  


  SUBROUTINE readMetadata(rl, recordLength, ispin_WAVECAR, singleOrDouble,&
       numKpoints, nkp, numBands, nb, enMax, latticeVector, eFermi)
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: rl
    REAL(8), INTENT(OUT) :: recordLength, ispin_WAVECAR, singleOrDouble
    REAL(8), INTENT(OUT) :: numKpoints, numBands, enMax, eFermi
    REAL(8), DIMENSION(3,3), INTENT(OUT) :: latticeVector
    INTEGER, INTENT(OUT) :: nb, nkp

    INTEGER :: i, j
    
    OPEN(UNIT=21,&
         FILE='WAVECAR_OLD',&
         FORM='unformatted',&
         ACCESS='direct',&
         STATUS='old',&
         RECL=rl)
    
    READ(21,REC=1) recordLength, ispin_WAVECAR, singleOrDouble
    READ(21,REC=2) numKpoints, numBands, enMax,&
         ((latticeVector(i,j), i=1,3), j=1,3),&
         eFermi
    nb = INT(numBands)
    nkp = INT(numKpoints)
    CLOSE(21)
  END SUBROUTINE readMetadata


  


  SUBROUTINE parseINCAR(ncl, ispin_INCAR, isym)
    IMPLICIT NONE
    
    LOGICAL, INTENT(OUT) :: ncl
    REAL(8), INTENT(OUT) :: ispin_INCAR
    INTEGER, INTENT(OUT) :: isym

    INTEGER :: fileSize, newLineCount, i, j, ioStatus
    CHARACTER(LEN=:), ALLOCATABLE :: incarStr, tmpKey, tmpVal
    INTEGER, DIMENSION(100) :: tmpNewLineIndices
    INTEGER, ALLOCATABLE :: newLineIndices(:)
    INTEGER :: startInd, endInd, commentInd, assignmentInd
    LOGICAL :: isymExists = .FALSE.

    ! Read INCAR !
    OPEN(UNIT=99,&
         FILE='INCAR',&
         STATUS='old',&
         ACCESS='stream',&
         ACTION='read',&
         FORM='unformatted',&
         IOSTAT=ioStatus)
    
    IF (ioStatus /= 0) THEN
       WRITE(*,'(A)') 'Failed to read INCAR. Terminating program.'
       STOP 2
    END IF

    INQUIRE(UNIT=99, SIZE=fileSize)
    ALLOCATE(CHARACTER(LEN=fileSize) :: incarStr)
    READ(UNIT=99, POS=1) incarStr
    
    CLOSE(99)

    ! Find number and indices of newlines !
    newLineCount = 0
    tmpNewLineIndices(1) = 0
    DO i = 1, LEN(incarStr)
       IF (incarStr(i:i) == NEW_LINE('')) THEN
          newLineCount = newLineCount + 1
          tmpNewLineIndices(newLineCount+1) = i
       END IF
    END DO

    tmpNewLineIndices(newLineCount+2) = LEN(incarStr) + 1
    ALLOCATE(newLineIndices(newLineCount+2))
    DO i = 1, newLineCount+2
       newLineIndices(i) = tmpNewLineIndices(i)
    END DO

    ! Read INCAR string line by line and search for LNONCOLLINEAR and ISPIN tags !
    DO i = 1, newLineCount+1
       startInd = newLineIndices(i) + 1
       endInd = newLineIndices(i+1) - 1
       commentInd = INDEX(incarStr(startInd:endInd), '#')
       assignmentInd = INDEX(incarStr(startInd:endInd), '=')
       ! If '=' exists and is to the left of '#'
       IF ((assignmentInd /= 0 .AND. assignmentInd < commentInd) .OR. &
            ! Or if '=' exists and '#' does not
            (assignmentInd > 0 .AND. commentInd == 0)) THEN
          tmpKey = TRIM(ADJUSTL(incarStr(startInd:startInd+assignmentInd-2)))
          SELECT CASE (tmpKey)
          CASE('LNONCOLLINEAR')
             IF (commentInd == 0) THEN
                tmpVal = TRIM(ADJUSTL(incarStr(startInd+assignmentInd:endInd)))
             ELSE
                tmpVal = TRIM(ADJUSTL(&
                     incarStr(startInd+assignmentInd:startInd+commentInd-2)))
             END IF
             READ(tmpVal,*) ncl
          CASE('ISPIN')
             IF (commentInd == 0) THEN
                tmpVal = TRIM(ADJUSTL(incarStr(startInd+assignmentInd:endInd)))
             ELSE
                tmpVal = TRIM(ADJUSTL(&
                     incarStr(startInd+assignmentInd:startInd+commentInd-2)))
             END IF
             READ(tmpVal,*) ispin_INCAR
          CASE('ISYM')
             IF (commentInd == 0) THEN
                tmpVal = TRIM(ADJUSTL(incarStr(startInd+assignmentInd:endInd)))
             ELSE
                tmpVal = TRIM(ADJUSTL(&
                     incarStr(startInd+assignmentInd:startInd+commentInd-2)))
             END IF
             READ(tmpVal,*) isym
             isymExists = .TRUE.
          END SELECT
       END IF
    END DO
    IF (.NOT. isymExists) THEN
       WRITE(*,'(A)') 'Error! ISYM not found in INCAR!'
       WRITE(*,'(A)') 'This program requires that ISYM be set explicitly!'
       WRITE(*,'(A)') 'Terminating program!'
       STOP 2
    END IF
  END SUBROUTINE parseINCAR
  
  
  
    
  
  SUBROUTINE calcMom(enMax, latticeVector, gxGrid, gyGrid, gzGrid, gMax, izMax,&
       isym, xMom, yMom, zMom, totMom)
    IMPLICIT NONE
    
    REAL(8), INTENT(IN) :: enMax
    REAL(8), DIMENSION(3,3), INTENT(IN) :: latticeVector
    INTEGER, INTENT(IN) :: isym
    REAL(8), INTENT(INOUT) :: gxGrid, gyGrid, gzGrid, gMax
    INTEGER, INTENT(INOUT) :: izMax
    REAL(8), DIMENSION(*), INTENT(OUT) :: xMom(:), yMom(:), zMom(:), totMom(:)

    REAL(8), DIMENSION(3) :: lenLatVec
    INTEGER :: i, pw, x, y, z, iyMax, ixMax
    REAL(8) :: gx, gy, gz
    
    DO i = 1, 3
       lenLatVec(i) = SQRT(latticeVector(i,1)**2&
            +latticeVector(i,2)**2&
            +latticeVector(i,3)**2)
    END DO
    
    gxGrid = 2*pi/(lenLatVec(1)*angstroms_to_au)
    gyGrid = 2*pi/(lenLatVec(2)*angstroms_to_au)
    gzGrid = 2*pi/(lenLatVec(3)*angstroms_to_au)
    WRITE(*,'(A,F18.4)') 'gxGrid = ', gxGrid
    WRITE(*,'(A,F18.4)') 'gyGrid = ', gyGrid
    WRITE(*,'(A,F18.4)') 'gzGrid = ', gzGrid
    WRITE(*,*) ''
    
    gMax = SQRT((2*m/hbar**2)*enMax*eV_to_au)
    WRITE(*, '(A,F18.4)') 'gMax = ', gMax
    WRITE(*,*) ''
    
    izMax = INT(gMax/gzGrid)
    WRITE(*,'(A,I0.1)') 'izMax = ', izMax
    WRITE(*,*) ''
    
    pw = 0
    SELECT CASE(isym)
    CASE(0)
       DO z = 0, 2*izMax
          iz = z - (2*izMax + 1)*INT(z/(izMax+1)) ! z=0,1,...,zMax,-zMax,...,-1
          gz = iz*gzGrid
          iyMax = INT(SQRT(gMax**2 - gz**2)/gyGrid)
          DO y = 0, 2*iyMax
             iy = y - (2*iyMax + 1)*INT(y/(iyMax+1))
             gy = iy*gyGrid
             ixMax = INT(SQRT(gMax**2 - gz**2 - gy**2)/gxGrid)
             DO x = 0, ixMax
                ix = x
                IF (ix==0 .AND. iy<0) CYCLE
                IF (ix==0 .AND. iy==0 .AND. iz<0) CYCLE
                gx = ix*gxGrid
                pw = pw + 1
                xMom(pw) = gx
                yMom(pw) = gy
                zMom(pw) = gz
                totMom(pw) = SQRT(gx**2 + gy**2 + gz**2)
             END DO
          END DO
       END DO
       
    CASE DEFAULT
       WRITE(*,'(A,I0.1,A)') 'ISYM=', isym, ' not implemented yet.'
       WRITE(*,'(A)') 'Terminating program!'
       STOP 2
    END SELECT
    
  END SUBROUTINE calcMom
  


  
  
END PROGRAM UnnamedJank
