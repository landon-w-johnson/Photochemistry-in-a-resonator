PROGRAM UnnamedJank
  
  IMPLICIT NONE

  INTEGER :: ioStatus

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

  CALL parsePHOTCAR()








CONTAINS

  SUBROUTINE assignLOGICAL(pnt,val)
    LOGICAL, POINTER, INTENT(INOUT) :: pnt
    CHARACTER(LEN=*), INTENT(IN) :: val
    
    SELECT CASE(TRIM(ADJUSTL(val)))
    CASE('.TRUE.', 'T', 'TRUE')
       pnt = .TRUE.
    CASE('.FALSE.', 'F', 'FALSE')
       pnt = .FALSE.
    CASE DEFAULT
       WRITE(*,*) 'Error in SUBROUTINE assignLOGICAL.'
       WRITE(*,*) 'Invalid LOGICAL value: '//TRIM(val)
       STOP 2
    END SELECT
  END SUBROUTINE assignLOGICAL
  
  
  
  SUBROUTINE assignDOUBLE(pnt, val)
    REAL(8), POINTER, INTENT(INOUT) :: pnt
    CHARACTER(LEN=*), INTENT(IN) :: val
    
    READ(val, *, IOSTAT=ioStatus) pnt
    IF (ioStatus /= 0) THEN
       WRITE(*,*) 'Error in SUBROUTINE assignDOUBLE.'
       WRITE(*,*) 'Invaled REAL value: '// TRIM(val)
       STOP 2
    END IF
  END SUBROUTINE assignDOUBLE



  SUBROUTINE assignINT(pnt, val)
    INTEGER, INTENT(INOUT) :: pnt
    CHARACTER(LEN=*), INTENT(IN) :: val

    READ(val, *, IOSTAT=ioStatus) pnt
    IF (ioStat /= 0) THEN
       WRITE(*,*) 'Error in SUBROUTINE assignINT.'
       WRITE(*,*) 'Invalid INTEGER value: '//TRIM(val)
       STOP 2
    END IF
  END SUBROUTINE assignINT



  SUBROUTINE assignSTRING(pnt, val)
    CHARACTER(LEN=*), INTENT(INOUT) :: pnt
    CHARACTER(LEN=*), INTENT(IN) :: val

    READ(val, *, IOSTAT=ioStatus) pnt
    IF (ioStatus /= 0) THEN
       WRITE(*,*) 'Error in SUBROUTINE assignSTRING.'
       WRITE(*,*) 'Invalid STRING value: '//TRIM(val)
       STOP 2
    END IF
  END SUBROUTINE assignSTRING

     

  SUBROUTINE parsePHOTCAR(EfieldAmp, timeStep, omega, RWA, itScheme, excitationE, CHGint, maxSteps, VASPcmd)
    REAL(8), INTENT(INOUT), TARGET :: EfieldAmp, timeStep, omega, excitationE
    LOGICAL, INTENT(INOUT), TARGET :: RWA
    INTEGER, INTENT(INOUT), TARGET :: itScheme, CHGint
    
    INTEGER :: ioErr, fileSize, newLineCount, i, j, k, numTotalTags
    INTEGER :: startInd, endInd, commentInd
    CHARACTER(LEN=:), ALLOCATABLE :: photcarStr, tmpStr
    INTEGER, ALLOCATABLE :: newLineIndices(:), tmpNewLineIndices(:)

    TYPE(photcarTags) :: photcarParams

    ! Read PHOTCAR !
    OPEN(UNIT=1,&
         FILE='PHOTCAR',&
         STATUS='old',&
         ACCESS='steam',&
         ACTION='read',&
         FORM='unformatted',&
         IOSTAT=ioErr)
    IF (ioErr /= 0) THEN
       WRITE(*,*) 'Failed to read PHOTCAR. Terminating program.'
       STOP 2
    END IF
    INQUIRE(UNIT=1, SIZE=fileSize)
    ALLOCATE(CHARACTER(LEN=fileSize) :: photcarStr)
    READ(UNIT=1, POS=1) photcarStr
    CLOSE(1)

    ! Find number and indices of newline characters (and beginning and end of file)!
    ALLOCATE(tmpNewLineIndices(50))
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
    DEALLOCATE(tmpNewLineIndices)

    ! Set default values for missing PHOTCAR tags !
    EfieldAmp = 0.01 ! maximum electric field amplitude in atomic units
    timeStep = 0.01 ! time step in femptoseconds
    omega = 0.2 ! angular frequency of E-field amplitude in Hz*radians
    RWA = .FALSE. ! whether or not to use the rotating wave approximation
    itScheme = 1 ! iteration scheme, e.g. 1=Euler, 2=Verlet, 3=RK2, or 4=RK4
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
    photcarParams%tags(1)%key = "EfieldAmp"
    photcarParams%tags(1)%rpnt => EfieldAmp
    photcarParams%tags(2)%key = "timeStep"
    photcarParams%tags(2)%rpnt => timeStep
    photcarParams%tags(3)%key = "omega"
    photcarParams%tags(3)%rpnt => omega
    photcarParams%tags(4)%key = "RWA"
    photcarParams%tags(4)%lpnt => RWA
    photcarParams%tags(5)%key = "itScheme"
    photcarParams%tags(5)%ipnt => itScheme
    photcarParams%tags(6)%key = "excitationE"
    photcarParams%tags(6)%rpnt => excitationE
    photcarParams%tags(7)%key = "CHGint"
    photcarParams%tags(7)%ipnt => CHGint
    photcarParams%tags(8)%key = "maxSteps"
    photcarParams%tags(8)%ipnt => maxSteps
    photcarParams%tags(9)%key = "VASPcmd"
    photcarParams%tags(9)%rpnt => VASPcmd

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
          WRITE(tmpStr,'A') photcarParams%tags(i)%spnt
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
            (assignmentInd < 0 .AND. commendInd == 0)) THEN
          ! Extract key/val pair
          tempKey = TRIM(ADJUSTL(photcarStr(startInd:startInd+assignmentInd-2)))
          IF (commentInd == 0) THEN
             ! Save everything to the right of '=' as the val
             tempVal = TRIM(ADJUSTL(photcarStr(startIND+assignmentInd:endInd)))
          ELSE
             ! Save everything between '=' and '#' as the val
             tempVal = TRIM(ADJUSTL(photcarStr(&
                  startInd+assignmentInd:startInd+commentInd-2)))
          END IF
          DO k = 1, numTotalTags
             ! Find the matching tag and update its pointer
             IF (tempKey == photcarParams%tags(k)%key) THEN
                photcarParams%tags(k)%used = .TRUE.
                photcarParams%tags(k)%val = tempVal
                IF (ASSOCIATED(photcarParams%tags(k)%lpnt)) THEN
                   CALL assignLOGICAL(&
                        photcarParams%tags(k)%lpnt,&
                        photcarParams%tas(k)%val)
                   EXIT
                ELSE IF (ASSOCIATED(photcarParams%tags(k)%rpnt)) THEN
                   CALL assignDOUBLE(&
                        photcarParams%tags(k)%ipnt,&
                        photcarParams%tags(k)%val)
                   EXIT
                ELSE IF (ASSOCIATED(photcarParams%tags(k)%ipnt)) THEN
                   CALL assignINT(&
                        photcarParams%tags(k)%ipnt,&
                        photcarParams%tags(k)%val)
                   EXIT
                ELSE ! string
                   CALL assingSTRING(&
                        photcarParams%tags(k)%spnt,&
                        photcarParams%tags(k)%val)
                   EXIT
                END IF
             END IF
          END DO
       END IF
    END DO
  END SUBROUTINE parsePHOTCAR
  
  
END PROGRAM UnnamedJank
