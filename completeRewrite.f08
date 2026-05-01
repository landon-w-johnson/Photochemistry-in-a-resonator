PROGRAM UnnamedJank
  
  IMPLICIT NONE

  ! Data structures for PHOTCAR tags
  
  TYPE :: tagStruct
     CHARACTER(LEN=:), ALLOCATABLE :: key, val
     LOGICAL, POINTER :: lpnt => NULL()
     REAL(8), POINTER :: rpnt => NULL()
     INTEGER, POINTER :: ipnt => NULL()
     CHARACTER(LEN=:), POINTER :: spnt => NULL()
     LOGICAL :: used = .FALSE., updated = .FALSE.
  END TYPE tagStruct

  TYPE :: photcarTags
     TYPE(tagStruct), ALLOCATABLE :: tags(:)
  END TYPE photcarTags



  ! Data structures for included transitions
  
  TYPE :: includedTransition
     INTEGER :: lb, ub ! lower band and upper band
     REAL(8) :: omega_0 ! resonant angular frequency of the transition
     COMPLEX(8), DIMENSION(2,3) :: tdm ! transition dipole moment (up/down,axis)
     ! NOTE: tdm is not literally the transition dipole moment !
     ! tdm(1,1) = <\psi_ub|p_x|\psi_lb>  and  tdm(2,1) = <\psi_lb|p_x|\psi_ub> !
  END TYPE INCLUDEDTRANSITION



  

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
  REAL(8), PARAMETER :: au_to_eV = 27.211386245988


  
  ! PHOTCAR tags !
  
  REAL(8), TARGET :: EfieldAmp, timeStep, omega, excitationE, resWidth
  LOGICAL, TARGET :: RWA
  INTEGER, TARGET :: intScheme, CHGint, maxSteps, polDir, POSint
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
  
  INTEGER :: i, j, ioStatus, ispin, spinDownRec, tInd
  INTEGER :: izMax, ix, iy, iz
  INTEGER :: numIncTrans, numIncTransDown, numIncBands, numIncBandsDown
  INTEGER :: vaspStat
  INTEGER(8) :: startTime, endTime, runTime
  REAL(8) :: dt, t, t_fs
  REAL(8) :: gxGrid, gyGrid, gzGrid, gMax
  REAL(8) :: lowerCutoff, upperCutoff, transE
  REAL(8), DIMENSION(3) :: bloch, blochDown
  REAL(8), ALLOCATABLE :: bandEn(:), bandEnDown(:)
  REAL(8), ALLOCATABLE :: gweight(:), gweightDown(:)
  REAL(8), ALLOCATABLE :: occ(:), occDown(:)
  COMPLEX(8), ALLOCATABLE :: bCoef(:), d_bCoef(:)
  COMPLEX(8), ALLOCATABLE :: bCoefDown(:), d_bCoefDown(:)
  COMPLEX(8), ALLOCATABLE :: pertHam(:,:), pertHamDown(:,:)
  COMPLEX(8), ALLOCATABLE :: interBandInProd(:,:), interBandInProdDown(:,:)
  REAL(8), ALLOCATABLE :: xMom(:), yMom(:), zMom(:), totMom(:)
  REAL(8), ALLOCATABLE :: kPts(:,:)
  COMPLEX(4), ALLOCATABLE :: pwCoefOld(:,:), pwCoefOldDown(:,:)
  COMPLEX(4), ALLOCATABLE :: pwCoefNew(:,:), pwCoefNewDown(:,:)
  COMPLEX(8), ALLOCATABLE :: d_pwCoef(:,:), d_pwCoefDown(:,:)
  INTEGER, ALLOCATABLE :: tmpTransBands(:,:)
  INTEGER, ALLOCATABLE :: tmpIncBands(:), incBands(:)
  INTEGER, ALLOCATABLE :: tmpTransBandsDown(:,:)
  INTEGER, ALLOCATABLE :: tmpIncBandsDown(:), incBandsDown(:)
  CHARACTER(LEN=72) :: dipoleFmt, occFileFmt, occFileHeadFmt
  CHARACTER(LEN=72) :: occFileFmtDown, occFileHeadFmtDown
  CHARACTER(LEN=5) :: numIncBandsStr, numIncBandsDownStr
  CHARACTER(LEN=10) :: tIndStr
  LOGICAL :: occFileExists, occFileDownExists
  
  TYPE(includedTransition), ALLOCATABLE :: trans(:), transDown(:)
  
  


  
  

  


  
!!!!!!!!!!!!!!!!!!!! BEGIN MAIN PROGRAM !!!!!!!!!!!!!!!!!!!!



  ! Parse PHOTCAR tags and set the corresponding variables!
  
  CALL parsePHOTCAR(EfieldAmp, timeStep, omega, RWA, intScheme, excitationE, CHGint, maxSteps, VASPcmd, POSint)

  WRITE(*,'(A)') 'Using the following values for PHOTCAR tags:'
  WRITE(*,'(T1,A,ES23.16E3)') '  EfieldAmp = ', EfieldAmp
  WRITE(*,'(T1,A,ES23.16E3)') '  timeStep = ', timeStep
  WRITE(*,'(T1,A,ES23.16E3)') '  omega = ', omega
  WRITE(*,'(T1,A,L1)') '  RWA = ', RWA
  WRITE(*,'(T1,A,I0.1)') '  intScheme = ', intScheme
  WRITE(*,'(T1,A,ES23.16E3)') '  excitationE = ', excitationE
  WRITE(*,'(A,F10.6)') '  resWidth = ', resWidth
  WRITE(*,'(A,I0.1)') '  polDir = ', polDir
  WRITE(*,'(T1,A,I0.1)') '  CHGint = ', CHGint
  WRITE(*,'(T1,A,I0.1)') '  POSint = ', POSint
  WRITE(*,'(T1,A,I0.1)') '  maxSteps = ', maxSteps
  WRITE(*,'(T1,A,A)') '  VASPcmd = ', TRIM(VASPcmd)
  WRITE(*,*) ''

  dt = timeStep*fs_to_au
  t = 0

  

  ! Extract info from WAVECAR !
  
  rl = getRecordLength()
  WRITE(*,'(A,I0.1)') 'record length of WAVECAR_OLD = ', rl
  WRITE(*,*) ''
  
  CALL readMetadata(rl, recordLength, ispin_WAVECAR, singleOrDouble,&
       numKpoints, nkp, numBands, nb, enMax, latticeVector, eFermi)

  WRITE(*,'(A,I0.1)') 'number of bands in WAVECAR_OLD = ', nb
  WRITE(*,'(A,I0.1)') 'number of k points in WAVECAR_OLD = ', nkp
  WRITE(*,*) ''



  ! Establish resonance cutoff !

  !resWidth = 1 ! (2*pi*hbar/dt)*au_to_eV ! (based on resonance curve)
  lowerCutoff = excitationE - resWidth
  upperCutoff = excitationE + resWidth
  WRITE(*,'(A,F10.4,A)') 'Using a resonance width of ', resWidth, ' eV'
  WRITE(*,'(A,F10.4,A,F10.4,A)')&
       'Only calculating transitions within energy range ',&
       lowerCutoff, ' to ', upperCutoff, ' eV'
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
        
        ALLOCATE(bandEn(nb))
        ALLOCATE(bandEnDown(nb))
        ALLOCATE(gweight(nb))
        ALLOCATE(gweightDown(nb))
        ALLOCATE(occ(nb))
        ALLOCATE(occDown(nb))
        ALLOCATE(kPts(nkp,3))
        ALLOCATE(pertHam(nb,nb))
        ALLOCATE(pertHamDown(nb,nb))
        ALLOCATE(tmpTransBands(nb*(nb-1),2))
        ALLOCATE(tmpTransBandsDown(nb*(nb-1),2))
        ALLOCATE(interBandInProd(nb,nb))
        ALLOCATE(interBandInProdDown(nb,nb))
        

        
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
             (bandEn(j), gweight(j), occ(j), j=1,nb)
        npw = INT(numPlaneWaves)
        ALLOCATE(pwCoefOld(nb,npw))
        ALLOCATE(pwCoefOldDown(nb,npw))
        ALLOCATE(pwCoefNew(nb,npw))
        ALLOCATE(pwCoefNewDown(nb,npw))
        ALLOCATE(d_pwCoef(nb,npw))
        ALLOCATE(d_pwCoefDown(nb,npw))
        DO i = 4, 3+nb
           READ(21,REC=i) (pwCoefNew(i-3,j), j=1,npw)
        END DO
        ! Spin down !
        spinDownRec = 4 + nb
        READ(21,REC=spinDownRec) numPlaneWaves, (kPts(1,i), i=1,3),&
             (bandEnDown(j), gweightDown(j), occDown(j), j=1,nb)
        DO i = 5+nb, 4+nb*2
           READ(21,REC=i) (pwCoefNewDown(i-4-nb,j), j=1,npw)
        END DO
        ! The plane waves are calculated later in the time iteration loop !
        ! This is needed again because of the rewriting of pwCoefOld with pwCoefNew !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CLOSE(21)

        WRITE(*,*) 'kPts:'
        WRITE(*,*) kPts
        WRITE(*,*) ''
        WRITE(*,*) 'bandEn:'
        WRITE(*,*) bandEn
        WRITE(*,*) ''
        WRITE(*,*) 'bandEnDown:'
        WRITE(*,*) bandEnDown
        WRITE(*,*) ''
        WRITE(*,*) 'occ:'
        WRITE(*,*) occ
        WRITE(*,*) ''
        WRITE(*,*) 'occDown:'
        WRITE(*,*) occDown
        WRITE(*,*) ''
        WRITE(*,'(A,F18.4)') 'numPlaneWaves = ', numPlaneWaves
        WRITE(*,'(A,I0.1)') 'npw = ', npw
        WRITE(*,*) ''
        
        
        
        ! Calculate momentum for each plane wave !

        ALLOCATE(xMom(npw))
        ALLOCATE(yMom(npw))
        ALLOCATE(zMom(npw))
        ALLOCATE(totMom(npw))

        CALL calcRecipGrid(enMax, latticeVector, gxGrid, gyGrid, gzGrid, gMax, izMax)
        CALL calcMom(gxGrid, gyGrid, gzGrid, gMax, izMax, isym, &
             xMom, yMom, zMom, totMom)


        
        ! Set text output format for transition dipole momenta !

        dipoleFmt = '(A,"(",ES23.15E3," , ",ES23.15E3," )")'
        
        
        
        ! Find near-resonant transitions to perform calculations on !

        DO i = 1, nb*(nb-1)
           tmpTransBands(i,1) = 0
           tmpTransBands(i,2) = 0
           tmpTransBandsDown(i,1) = 0
           tmpTransBandsDown(i,2) = 0
        END DO

        numIncTrans = 0
        DO i = 1, nb-1
           DO j = i+1, nb
              transE = bandEn(j) - bandEn(i)
              IF (transE >= lowerCutoff) THEN
                 IF (transE <= upperCutoff) THEN
                    numIncTrans = numIncTrans + 1
                    tmpTransBands(numIncTrans,1) = i
                    tmpTransBands(numIncTrans,2) = j
                 ELSE
                    EXIT
                 END IF
              END IF
           END DO
        END DO

        numIncTransDown = 0
        DO i = 1, nb-1
           DO j = i+1, nb
              transE = bandEnDown(j) - bandEnDown(i)
              IF (transE >= lowerCutoff) THEN
                 IF (transE <= upperCutoff) THEN
                    numIncTransDown = numIncTransDown + 1
                    tmpTransBandsDown(numIncTransDown,1) = i
                    tmpTransBandsDown(numIncTransDown,2) = j
                 ELSE
                    EXIT
                 END IF
              END IF
           END DO
        END DO

        IF (numIncTrans == 0 .AND. numIncTransDown == 0) THEN
           WRITE(*,'(A)') 'No transitions included in calculations!'
           WRITE(*,'(A)') 'Check excitationE in your PHOTCAR against the energy differences between bands in your previous OUTCAR.'
           WRITE(*,'(A)') 'Terminating program!'
           STOP 2
        END IF

        WRITE(*,*) 'tmpTransBands:'
        DO i = 1, nb*(nb-1)
           WRITE(*,*) tmpTransBands(i,1), tmpTransBands(i,2)
        END DO
        WRITE(*,*) ''
        WRITE(*,*) 'tmpTransBandsDown:'
        DO i = 1, nb*(nb-1)
           WRITE(*,*) tmpTransBandsDown(i,1), tmpTransBandsDown(i,2)
        END DO
        WRITE(*,*) ''
        
        ALLOCATE(tmpIncBands(2*numIncTrans))
        ALLOCATE(tmpIncBandsDown(2*numIncTransDown))
        ALLOCATE(trans(numIncTrans))
        ALLOCATE(transDown(numIncTransDown))
        
        DO i = 1, numIncTrans
           trans(i)%lb = tmpTransBands(i,1)
           trans(i)%ub = tmpTransBands(i,2)
           trans(i)%omega_0 = ((bandEn(tmpTransBands(i,2)) - bandEn(tmpTransBands(i,1)))/hbar)*eV_to_au
           tmpIncBands(i) = 0
           tmpIncBands(numIncTrans+i) = 0
        END DO

        DO i = 1, numIncTransDown
           transDown(i)%lb = tmpTransBandsDown(i,1)
           transDown(i)%ub = tmpTransBandsDown(i,2)
           transDown(i)%omega_0 = ((bandEnDown(tmpTransBandsDown(i,2)) - bandEnDown(tmpTransBandsDown(i,1)))/hbar)*eV_to_au
           tmpIncBandsDown(i) = 0
           tmpIncBandsDown(numIncTransDown+i) = 0
        END DO

        DEALLOCATE(tmpTransBands)
        DEALLOCATE(tmpTransBandsDown)

        WRITE(*,'(A)') 'SPIN-UP included transitions:'
        DO i = 1, numIncTrans
           !WRITE(*,'(A,I0.1,A,I0.1,A,F10.6,A,A,F10.6,A)') &
           WRITE(*,'(A,I0.1,A,I0.1,A,F10.6,A,A,F16.12,A)') &
                '  ', trans(i)%lb, ' -> ', trans(i)%ub, ' : ', &
                trans(i)%omega_0*au_to_eV, ' eV,', &
                ' omega_0 = ', trans(i)%omega_0, ' a.u.'
        END DO
        WRITE(*,*) ''
        WRITE(*,'(A)') 'SPIN-DOWN included transitions:'
        DO i = 1, numIncTransDown
           WRITE(*,'(A,I0.1,A,I0.1,A,F10.6,A,A,F10.6,A)') &
                '  ', transDown(i)%lb, ' -> ', transDown(i)%ub, ' : ', &
                transDown(i)%omega_0*au_to_eV, ' eV,', &
                ' omega_0 = ', transDown(i)%omega_0, ' a.u.'
        END DO
        WRITE(*,*) ''

        numIncBands = 0
        DO i = 1, numIncTrans
           IF (.NOT. ANY(tmpIncBands==trans(i)%lb) ) THEN
              numIncBands = numIncBands + 1
              tmpIncBands(numIncBands) = trans(i)%lb
           END IF
           IF (.NOT. ANY(tmpIncBands==trans(i)%ub) ) THEN
              numIncBands = numIncBands + 1
              tmpIncBands(numIncBands) = trans(i)%ub
           END IF
        END DO

        WRITE(*,*) 'tmpIncBands:'
        WRITE(*,*) tmpIncBands

        numIncBandsDown = 0
        DO i = 1, numIncTransDown
           IF (.NOT. ANY(tmpIncBandsDown==transDown(i)%lb) ) THEN
              numIncBandsDown = numIncBandsDown + 1
              tmpIncBandsDown(numIncBandsDown) = transDown(i)%lb
           END IF
           IF (.NOT. ANY(tmpIncBandsDown==transDown(i)%ub) ) THEN
              numIncBandsDown = numIncBandsDown + 1
              tmpIncBandsDown(numIncBandsDown) = transDown(i)%ub
           END IF
        END DO

        WRITE(*,*) 'tmpIncBandsDown:'
        WRITE(*,*) tmpIncBandsDown

        ALLOCATE(incBands(numIncBands))
        ALLOCATE(incBandsDown(numIncBandsDown))
        
        DO i = 1, numIncBands
           incBands(i) = tmpIncBands(i)
        END DO

        DO i = 1, numIncBandsDown
           incBandsDown(i) = tmpIncBandsDown(i)
        END DO

        DEALLOCATE(tmpIncBands)
        DEALLOCATE(tmpIncBandsDown)

        WRITE(*,'(A)') 'SPIN-UP bands included in calculations:'
        WRITE(*,*) incBands
        WRITE(*,*) ''
        WRITE(*,'(A)') 'SPIN-DOWN bands included in calculations:'
        WRITE(*,*) incBandsDown
        WRITE(*,*) ''
        
        
        
        ! Allocate data structures for the bands and initialize them
        
        ALLOCATE(bCoef(nb))
        ALLOCATE(d_bCoef(nb))
        ALLOCATE(bCoefDown(nb))
        ALLOCATE(d_bCoefDown(nb))

        DO i = 1, nb
           bCoef(i) = SQRT(occ(i))
           d_bCoef(i) = 0
           bCoefDown(i) = SQRT(occDown(i))
           d_bCoefDown(i) = 0
           DO j = 1, nb
              pertHam(i,j) = 0
              pertHamDown(i,j) = 0
           END DO
        END DO
        
        
        
        ! Check if this is a continuation job !
        
        occFileExists = fileExists('occupations.txt')
        occFileDownExists = fileExists('occupationsDOWN.txt')
        IF (occFileExists .OR. occFileDownExists) THEN
           WRITE(*,'(A)') 'Occupation file found'
           WRITE(*,'(A)') 'Skipping initialization'
           WRITE(*,*) ''
        ELSE
           WRITE(*,'(A)') 'No occupation files found'
           WRITE(*,'(A)') 'Initializing trajectory'
           WRITE(*,*) ''
           !CALL SYSTEM('cp WAVECAR WAVECAR_0')
           IF (numIncBands > 0) THEN
              WRITE(numIncBandsStr,'(I0.1)') numIncBands
              occFileHeadFmt = '(A10,'//numIncBandsStr//'(I24))'
              occFileFmt = '(F10.3,'//numIncBandsStr//'(ES24.15E3))'
              OPEN(UNIT=50, FILE='occupations.txt')
              WRITE(50, FMT=occFileHeadFmt) '     time', &
                   (incBands(i), i=1,numIncBands)
              IF (numIncBands == 2) THEN
                 OPEN(UNIT=60, FILE='bloch.txt')
                 WRITE(60, FMT='(A44)') ' time (fs)      x          y          z     '
              END IF
           END IF
           IF (numIncBandsDown > 0) THEN
              WRITE(numIncBandsDownStr,'(I0.1)') numIncBandsDown
              occFileHeadFmtDown = '(A10,'//numIncBandsDownStr//'(I24))'
              occFileFmtDown = '(F10.3,'//numIncBandsDownStr//'(ES24.15E3))'
              OPEN(UNIT=51, FILE='occupationsDOWN.txt')
              WRITE(51, FMT=occFileHeadFmtDown) '     time', &
                   (incBandsDown(i), i=1,numIncBandsDown)
              IF (numIncBandsDown == 2) THEN
                 OPEN(UNIT=61, FILE='blochDOWN.txt')
                 WRITE(61, FMT='(A44)') ' time (fs)      x          y          z     '
              END IF
           END IF
        END IF
        
        
        
        ! Update band expansion coefficients !
        
        DO tInd = 0, maxSteps
           startTime = getTime()
           t = tInd*dt
           t_fs = t*au_to_fs
           WRITE(*,*) ''
           WRITE(*,'(A,F10.3,A)') 't = ', t*au_to_fs, ' fs'
           WRITE(*,'(A)') '__________________________________________________'
           WRITE(*,*) ''
           ! Calculate transition dipole momenta !
           WRITE(*,'(A)') 'Preparing SPIN-UP data'
           CALL readBandEnsAndPWcoefs(rl, nb, npw, nkp, numIncTrans, numIncBands, &
                incBands, trans, pwCoefOld, pwCoefNew, interBandInProd)
           CALL calc_d_pwCoefs(dt, nb, npw, pwCoefOld, pwCoefNew, d_pwCoef)
           !CALL convertBandCoefs(nb, npw, pwCoefOld, pwCoefNew, bCoef)
           CALL calcTransDipMom(trans, numIncTrans, npw, xMom, yMom, zMom, pwCoefNew)
           CALL calc_pertHam(numIncTrans, EfieldAmp, omega, t, bandEn, &
                trans, polDir, pertHam)
           WRITE(*,'(A)') 'SPIN-UP perturbative Hamiltonian matrix elements:'
           WRITE(*,*) pertHam
           WRITE(*,'(A)') 'SPIN-UP transition dipole moments:'
           DO i = 1, numIncTrans
              WRITE(*,'(A,F10.6,A)') 'omega_0 = ', trans(i)%omega_0, ' a.u.'
              WRITE(*,'(A,I0.1,A,I0.1)') '  ', trans(i)%lb, ' -> ', trans(i)%ub
              WRITE(*,dipoleFmt) '    1: ', trans(i)%tdm(1,1)
              WRITE(*,dipoleFmt) '    2: ', trans(i)%tdm(1,2)
              WRITE(*,dipoleFmt) '    3: ', trans(i)%tdm(1,3)
              WRITE(*,'(A,I0.1,A,I0.1)') '  ', trans(i)%lb, ' <- ', trans(i)%ub
              WRITE(*,dipoleFmt) '    1: ', trans(i)%tdm(2,1)
              WRITE(*,dipoleFmt) '    2: ', trans(i)%tdm(2,2)
              WRITE(*,dipoleFmt) '    3: ', trans(i)%tdm(2,3)
           END DO
           WRITE(*,*) ''
           WRITE(*,'(A)') 'Preparing SPIN-DOWN data'
           CALL readBandEnsAndPWcoefsDOWN(rl, nb, npw, nkp, numIncTransDown, &
                numIncBandsDown, incBandsDown, transDown, &
                pwCoefOldDown, pwCoefNewDown, interBandInProdDown)
           !CALL convertBandCoefs(nb, npw, pwCoefOldDown, pwCoefNewDown, bCoefDown)
           CALL calc_d_pwCoefs(dt, nb, npw, pwCoefOldDown, pwCoefNewDown, d_pwCoefDown)
           CALL calcTransDipMom(transDown, numIncTransDown, npw, xMom, yMom, zMom, pwCoefNewDown)
           CALL calc_pertHam(numIncTransDown, EfieldAmp, omega, t, bandEnDown, &
                transDown, polDir, pertHamDown)
           WRITE(*,'(A)') 'SPIN-DOWN perturbative Hamiltonian matrix elements:'
           WRITE(*,*) pertHamDown
           WRITE(*,'(A)') 'SPIN-DOWN transition dipole moments:'
           DO i = 1, numIncTransDown
              WRITE(*,'(A,F10.6,A)') 'omega_0 = ', transDown(i)%omega_0, ' a.u.'
              WRITE(*,'(A,I0.1,A,I0.1)') '  ', transDown(i)%lb, ' -> ', transDown(i)%ub
              WRITE(*,dipoleFmt) '    1: ', transDown(i)%tdm(1,1)
              WRITE(*,dipoleFmt) '    2: ', transDown(i)%tdm(1,2)
              WRITE(*,dipoleFmt) '    3: ', transDown(i)%tdm(1,3)
              WRITE(*,'(A,I0.1,A,I0.1)') '  ', transDown(i)%lb, ' <- ', transDown(i)%ub
              WRITE(*,dipoleFmt) '    1: ', transDown(i)%tdm(2,1)
              WRITE(*,dipoleFmt) '    2: ', transDown(i)%tdm(2,2)
              WRITE(*,dipoleFmt) '    3: ', transDown(i)%tdm(2,3)
           END DO
           WRITE(*,*) ''
           ! Update band expansion coefficients
           DO i = 1, numIncBands
              d_bCoef(incBands(i)) = 0
           END DO
           DO i = 1, numIncBandsDown
              d_bCoefDown(incBandsDown(i)) = 0
           END DO
           !DO i = 1, numIncTrans
           !   WRITE(*,'(A,I0.1,A,I0.1)') &
           !        'calculating derivatives for SPIN-UP transition ', &
           !        trans(i)%lb, ' <-> ', trans(i)%ub
           !   ! Find derivatives of band expansion coefficients
           !   IF (occ(trans(i)%lb) == 0) THEN
           !      WRITE(*,'(A)') 'Lower band at 0 occupation'
           !      CALL calc_d_bCoefs_at0(t, omega, trans(i), &
           !           d_bCoef(trans(i)%ub), d_bCoef(trans(i)%lb), &
           !           bandEn(trans(i)%ub), bandEn(trans(i)%lb), &
           !           bCoef(trans(i)%ub), bCoef(trans(i)%lb), 2, polDir)
           !      WRITE(*,*) 'd_bCoef lower band =', d_bCoef(trans(i)%lb)
           !      WRITE(*,*) 'd_bCoef upper band =', d_bCoef(trans(i)%ub)
           !   ELSE IF (occ(trans(i)%ub) == 0) THEN
           !      WRITE(*,'(A)') 'Upper band at 0 occupation'
           !      CALL calc_d_bCoefs_at0(t, omega, trans(i), &
           !           d_bCoef(trans(i)%lb), d_bCoef(trans(i)%ub), &
           !           bandEn(trans(i)%lb), bandEn(trans(i)%ub), &
           !           bCoef(trans(i)%lb), bCoef(trans(i)%ub), 1, polDir)
           !      WRITE(*,*) 'd_bCoef lower band =', d_bCoef(trans(i)%lb)
           !      WRITE(*,*) 'd_bCoef upper band =', d_bCoef(trans(i)%ub)
           !   ELSE
           !      CALL calc_d_bCoefs(t, omega, trans(i), polDir, &
           !           bCoef(trans(i)%lb), bCoef(trans(i)%ub), &
           !           d_bCoef(trans(i)%lb), d_bCoef(trans(i)%ub))
           !      WRITE(*,*) 'd_bCoef lower band =', d_bCoef(trans(i)%lb)
           !      WRITE(*,*) 'd_bCoef upper band =', d_bCoef(trans(i)%ub)
           !   END IF
           !   WRITE(*,*) ''
           !END DO
           WRITE(*,'(A)') 'calculating SPIN-UP d_bCoefs'
           !DO i = 1, numIncBands
           !   CALL calc_d_bCoefs(nb, npw, t, omega, incBands(i), bCoef, bandEn, &
           !        pertHam, d_pwCoef, pwCoefNew, d_bCoef(incBands(i)))
           !END DO
           DO i = 1, numIncTrans
              CALL calc_d_bCoefs_2band(nb, npw, t, omega, bCoef, bandEn, &
                   pertHam, d_pwCoef, pwCoefNew, d_bCoef, &
                   trans(i), interBandInProd)
           END DO
           WRITE(*,*) ''
           !DO i = 1, numIncTransDown
           !   WRITE(*,'(A,I0.1,A,I0.1)') &
           !        'calculating derivatives for SPIN-DOWN transition ', &
           !        transDown(i)%lb, ' <-> ', transDown(i)%ub
           !   ! Find derivatives of band expansion coefficients !
           !   IF (occDown(transDown(i)%lb) == 0) THEN
           !      WRITE(*,'(A)') 'Lower band at 0 occupation (DOWN)'
           !      CALL calc_d_bCoefs_at0(t, omega, transDown(i), &
           !           d_bCoefDown(transDown(i)%ub), d_bCoefDown(transDown(i)%lb), &
           !           bandEnDown(transDown(i)%ub), bandEnDown(transDown(i)%lb), &
           !           bCoefDown(transDown(i)%ub), bCoefDown(transDown(i)%lb), 2, polDir)
           !      WRITE(*,*) 'd_bCoefDown lower band =', d_bCoefDown(transDown(i)%lb)
           !      WRITE(*,*) 'd_bCoefDown upper band =', d_bCoefDown(transDown(i)%ub)
           !   ELSE IF (occDown(transDown(i)%ub) == 0) THEN
           !      WRITE(*,'(A)') 'Upper band (DOWN) at 0 occupation'
           !      CALL calc_d_bCoefs_at0(t, omega, transDown(i), &
           !           d_bCoefDown(transDown(i)%lb), d_bCoefDown(transDown(i)%ub), &
           !           bandEnDown(transDown(i)%lb), bandEnDown(transDown(i)%ub), &
           !           bCoefDown(transDown(i)%lb), bCoefDown(transDown(i)%ub), 1, polDir)
           !      WRITE(*,*) 'd_bCoefDown lower band =', d_bCoefDown(transDown(i)%lb)
           !      WRITE(*,*) 'd_bCoefDown upper band =', d_bCoefDown(transDown(i)%ub)
           !   ELSE
           !      CALL calc_d_bCoefs(t, omega, transDown(i), polDir, &
           !           bCoefDown(transDown(i)%lb), bCoefDown(transDown(i)%ub), &
           !           d_bCoefDown(transDown(i)%lb), d_bCoefDown(transDown(i)%ub))
           !      WRITE(*,*) 'd_bCoefDown lower band =', d_bCoefDown(transDown(i)%lb)
           !      WRITE(*,*) 'd_bCoefDown upper band =', d_bCoefDown(transDown(i)%ub)
           !   END IF
           !   WRITE(*,*) ''
           !END DO
           WRITE(*,'(A)') 'Calculating SPIN-DOWN d_bCoefs'
           !DO i = 1, numIncBandsDown
           !   CALL calc_d_bCoefs(nb, npw, t, omega, incBandsDown(i), bCoefDown, &
           !        bandEnDown, pertHamDown, d_pwCoefDown, &
           !        pwCoefNewDown, d_bCoefDown(incBandsDown(i)))
           !END DO
           WRITE(*,*) ''
           ! Actually update band expansion coefficients !
           SELECT CASE (intScheme)
           CASE(1) ! Euler !
!!!!!!!!!!!!!!!!!!!! MODIFICATION !!!!!!!!!!!!!!!!!!!!
              IF (tInd==0) THEN
                 WRITE(*,*) 'special first step update...'
                 WRITE(*,*) ''
                 WRITE(*,*) 'bCoef upper before update:', bCoef(trans(1)%ub)
                 WRITE(*,*) 'd_bCoef upper before update:', d_bCoef(trans(1)%ub)
                 bCoef(trans(1)%ub) = bCoef(trans(1)%ub) + dt*d_bCoef(trans(1)%ub)
                 WRITE(*,'(A,I0.1,A,ES23.15E3,",",ES23.15E3,")")') &
                      'bCoef(', trans(1)%ub, ') = (', bCoef(trans(1)%ub)
                 WRITE(*,*) ''
                 WRITE(*,*) 'bCoef lower before update:', bCoef(trans(1)%lb)
                 bCoef(trans(1)%lb) = SQRT( (1,0) - bCoef(trans(1)%ub)*CONJG(bCoef(trans(1)%ub)) )
                 WRITE(*,'(A,I0.1,A,ES23.15E3,",",ES23.15E3,")")') &
                      'bCoef(', trans(1)%lb, ') = (', bCoef(trans(1)%lb)
                 WRITE(*,*) ''
              ELSE
!!!!!!!!!!!!!!!!!!!! END MODIFICATION !!!!!!!!!!!!!!!!!!!!
                 DO i = 1, numIncBands
                    WRITE(*,*) 'bCoef before update:', bCoef(incBands(i))
                    WRITE(*,*) 'd_bCoef before update:', d_bCoef(incBands(i))
                    bCoef(incBands(i)) = bCoef(incBands(i)) + dt*d_bCoef(incBands(i))
                    WRITE(*,'(A,I0.1,A,ES23.15E3,",",ES23.15E3,")")') &
                         'bCoef(', incBands(i), ') = (', bCoef(incBands(i))
                    WRITE(*,*) ''
                 END DO
                 DO i = 1, numIncBandsDown
                    WRITE(*,*) 'bCoefDown before update:', bCoefDown(incBandsDown(i))
                    WRITE(*,*) 'd_bCoefDown before update:', d_bCoefDown(incBandsDown(i))
                    bCoefDown(incBandsDown(i)) = bCoefDown(incBandsDown(i)) + &
                         dt*d_bCoefDown(incBandsDown(i))
                    WRITE(*,'(A,I0.1,A,ES23.15E3,",",ES23.15E3,")")') &
                         'bCoefDown(', incBandsDown(i), ') = (', bCoefDown(incBandsDown(i))
                    WRITE(*,*) ''
                 END DO
              END IF
              WRITE(*,*) ''
           CASE(2) ! Verlet !
              WRITE(*,'(A)') 'Verlet not implemented yet.'
              WRITE(*,'(A)') 'Terminating program!'
              STOP 2
           CASE(3) ! Runge-Kutta 2nd order !
              WRITE(*,'(A)') 'Runge-Kutta 2nd order not implemented yet.'
              WRITE(*,'(A)') 'Terminating program!'
              STOP 2
           CASE(4) ! Runge-Kutta 4th order !
              WRITE(*,'(A)') 'Runge-Kutta 4th order not implemented yet.'
              WRITE(*,'(A)') 'Terminating program!'
              STOP 2
           CASE DEFAULT
              WRITE(*,'(A,I0.1)') 'Error! Invalid intScheme: ', intScheme
              WRITE(*,'(A)') 'Terminating program!'
              STOP 2
           END SELECT
           ! Update occupations and overwrite WAVECAR !
           DO i = 1, numIncBands
              occ(incBands(i)) = bCoef(incBands(i))*CONJG(bCoef(incBands(i)))
              WRITE(*,'(A,I0.1,A,ES23.15E3)') &
                   'occ(', incBands(i), ') = ', occ(incBands(i))
           END DO
           DO i = 1, numIncBandsDown
              occDown(incBandsDown(i)) = &
                   bCoefDown(incBandsDown(i))*CONJG(bCoefDown(incBandsDown(i)))
              WRITE(*,'(A,I0.1,A,ES23.15E3)') &
                   'occDown(', incBandsDown(i), ') = ', occDown(incBandsDown(i))
           END DO
           CALL SYSTEM('cp WAVECAR WAVECAR_OLD')
           WRITE(tIndStr,'(I10)') tInd
           !CALL SYSTEM('cp WAVECAR WAVECAR_'//TRIM(ADJUSTL(tIndStr)))
           IF (POSint/=0 .AND. MODULO(tInd,POSint)==0) THEN
              CALL SYSTEM('cp POSCAR file-copies/POSCAR_'//TRIM(ADJUSTL(tIndStr)))
           END IF
           OPEN(UNIT=12,&
                FILE='WAVECAR',&
                FORM='unformatted',&
                ACCESS='direct',&
                STATUS='old',&
                RECL=8)
           DO i = 1, numIncBands
              WRITE(12,REC=rl/4+4+3*incBands(i)) occ(incBands(i))
           END DO
           DO i = 1, numIncBandsDown
              WRITE(12,REC=rl*3/8+rl*nb/8+4+3*incBandsDown(i)) occDown(incBandsDown(i))
           END DO
           CLOSE(12)
           endTime = getTime()
           runTime = endTime-startTime
           WRITE(*,*) ''
           WRITE(*,'(A,I0.1,A)') 'FORTRAN runtime: ', runTime, ' ms'
           WRITE(*,*) ''
           
           ! Run VASP !
           startTime = getTime()
           !CALL EXECUTE_COMMAND_LINE( &
           !     'mpiexec -machinefile $PBS_NODEFILE -np $NUM_PROC '// &
           !     VASPcmd, EXITSTAT=vaspStat) ! CCAST
           !CALL EXECUTE_COMMAND_LINE( &
           !     'srun -n 64 '//VASPcmd, EXITSTAT=vaspStat) ! NERSC
           CALL EXECUTE_COMMAND_LINE(VASPcmd, EXITSTAT=vaspStat) ! Trion
           !CALL SYSTEM('echo pretending to run VASP')
           endTime = getTime()
           runTime = endTime-startTime
           WRITE(*,*) ''
           WRITE(*,'(A,I0.1,A)') 'VASP runtime: ', runTime, ' ms'
           WRITE(*,*) ''
           IF (vaspStat /= 0) THEN
              WRITE(*,'(A)') 'VASP terminated with error. Terminating program.'
              STOP 2
           END IF
           CALL SYSTEM('cp CONTCAR POSCAR')
           
           ! Update occupations files !
           IF (numIncBands > 0) THEN
              WRITE(50, occFileFmt) t_fs, (occ(incBands(i)), i=1,numIncBands)
              ! Update Bloch sphere components if two-state system !
              IF (numIncBands == 2) THEN
                 bloch(1) = REALPART(CONJG(bCoef(incBands(2)))*bCoef(incBands(1)) &
                      + CONJG(bCoef(incBands(1)))*bCoef(incBands(2)))
                 WRITE(*,'(A)') 'complex part of bloch(1):'
                 WRITE(*,*) IMAGPART(CONJG(bCoef(incBands(2)))*bCoef(incBands(1)) &
                      + CONJG(bCoef(incBands(1)))*bCoef(incBands(2)))
                 bloch(2) = IMAGPART(CONJG(bCoef(incBands(2)))*bCoef(incBands(1)) &
                      - CONJG(bCoef(incBands(1)))*bCoef(incBands(2)))
                 WRITE(*,'(A)') 'real part of bloch(2):'
                 WRITE(*,*) REALPART(CONJG(bCoef(incBands(2)))*bCoef(incBands(1)) &
                      + CONJG(bCoef(incBands(1)))*bCoef(incBands(2)))
                 bloch(3) = occ(incBands(1)) - occ(incBands(2))
                 WRITE(60, '(F11.3,3(F11.6))') t_fs, (bloch(i), i=1,3)
              END IF
           END IF
           IF (numIncBandsDown > 0) THEN
              WRITE(51, occFileFmtDown) t_fs, (occDown(incBandsDown(i)), i=1,numIncBandsDown)
              IF (numIncBandsDown == 2) THEN
                 blochDown(1) = REALPART( &
                      CONJG(bCoefDown(incBandsDown(2)))*bCoefDown(incBandsDown(1)) &
                      + CONJG(bCoefDown(incBandsDown(1)))*bCoefDown(incBandsDown(2)))
                 WRITE(*,'(A)') 'complex part of blochDown(1):'
                 WRITE(*,*) IMAGPART( &
                      CONJG(bCoefDown(incBandsDown(2)))*bCoefDown(incBandsDown(1)) &
                      + CONJG(bCoefDown(incBandsDown(1)))*bCoefDown(incBandsDown(2)))
                 blochDown(2) = IMAGPART( &
                      CONJG(bCoefDown(incBandsDown(2)))*bCoefDown(incBandsDown(1)) &
                      - CONJG(bCoefDown(incBandsDown(1)))*bCoefDown(incBandsDown(2)))
                 WRITE(*,'(A)') 'real part of blochDown(2):'
                 WRITE(*,*) REALPART( &
                      CONJG(bCoefDown(incBandsDown(2)))*bCoefDown(incBandsDown(1)) &
                      + CONJG(bCoefDown(incBandsDown(1)))*bCoefDown(incBandsDown(2)))
                 blochDown(3) = occDown(incBandsDown(1)) - occDown(incBandsDown(2))
                 WRITE(61, '(F11.3,3(F11.6))') t_fs, (bloch(i), i=1,3)
              END IF
           END IF
           WRITE(*,'(A)') '__________________________________________________'
           WRITE(*,*) ''
           WRITE(*,*) ''
        END DO
        CLOSE(50)
        CLOSE(51)
        CLOSE(60)
        CLOSE(61)
        
        
        
     CASE DEFAULT ! undefined ISPIN !
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



  

  SUBROUTINE parsePHOTCAR(EfieldAmp, timeStep, omega, RWA, intScheme, excitationE, CHGint, maxSteps, VASPcmd, POSint)
    IMPLICIT NONE
    
    REAL(8), INTENT(OUT), TARGET :: EfieldAmp, timeStep, omega, excitationE
    LOGICAL, INTENT(OUT), TARGET :: RWA
    INTEGER, INTENT(OUT), TARGET :: intScheme, CHGint, maxSteps, POSint
    CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT), TARGET :: VASPcmd
    
    INTEGER :: ioStatus, fileSize, newLineCount, i, j, k, numTotalTags
    INTEGER :: startInd, endInd, commentInd, assignmentInd
    CHARACTER(LEN=:), ALLOCATABLE :: photcarStr, tmpKey, tmpVal
    CHARACTER(72) :: tmpStr
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
    excitationE = 1 ! energy of the photons in eV
    CHGint = 0 ! time steps between saved CHG files. 0=don't save any
    maxSteps = 0 ! maximum number of time steps to calculate. 0=no maximum
    VASPcmd = 'vasp_gam' ! VASP executable to run, e.g. 'vasp_gam', 'vasp_std', or 'vasp_ncl'
    polDir = 1
    resWidth = 0.1
    POSint = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    numTotalTags = 12 ! UPDATE THIS W/ NEW PHOTCAR TAGS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Instantiate photcarParams key/val pairs !
    ALLOCATE(photcarParams%tags(numTotalTags))

    ! Set keys for PHOTCAR parameters and associate the pointers to the values !
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
    photcarParams%tags(10)%key = 'polDir'
    photcarParams%tags(10)%ipnt => polDir
    photcarParams%tags(11)%key = 'resWidth'
    photcarParams%tags(11)%rpnt => resWidth
    photcarParams%tags(12)%key = 'POSint'
    photcarParams%tags(12)%ipnt => POSint

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

    ! Create a subdirectory to store copies of files !
    IF (CHGint/=0 .OR. POSint/=0) THEN
       CALL SYSTEM('mkdir -p file-copies')
    END IF
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


  


  SUBROUTINE readMetadata(rl, recordLength, ispin_WAVECAR, singleOrDouble, &
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




  SUBROUTINE calcRecipGrid(enMax, latticeVector, gxGrid, gyGrid, gzGrid, gMax, izMax)
    IMPLICIT NONE

    REAL(8), INTENT(IN) :: enMax
    REAL(8), DIMENSION(3,3), INTENT(IN) :: latticeVector
    REAL(8), INTENT(OUT) :: gxGrid, gyGrid, gzGrid, gMax
    INTEGER, INTENT(OUT) :: izMax

    REAL(8), DIMENSION(3) :: lenLatVec
    INTEGER :: i

    DO i = 1, 3
       lenLatVec(i) = SQRT(latticeVector(i,1)**2&
            +latticeVector(i,2)**2&
            +latticeVector(i,3)**2)
       WRITE(*,*) 'lenLatVec', i, ':', lenLatVec(i)
    END DO
    WRITE(*,*) ''
    
    gxGrid = 2*pi/(lenLatVec(1)*angstroms_to_au)
    gyGrid = 2*pi/(lenLatVec(2)*angstroms_to_au)
    gzGrid = 2*pi/(lenLatVec(3)*angstroms_to_au)
    WRITE(*,'(A,F18.14)') 'gxGrid = ', gxGrid
    WRITE(*,'(A,F18.14)') 'gyGrid = ', gyGrid
    WRITE(*,'(A,F18.14)') 'gzGrid = ', gzGrid
    WRITE(*,*) ''
    
    gMax = SQRT((2*m/hbar**2)*enMax*eV_to_au)
    WRITE(*, '(A,F18.14)') 'gMax = ', gMax
    WRITE(*,*) ''
    
    izMax = INT(gMax/gzGrid)
    WRITE(*,'(A,I0.1)') 'izMax = ', izMax
    WRITE(*,*) ''
  END SUBROUTINE calcRecipGrid
  
  
  
    
  
  SUBROUTINE calcMom(gxGrid, gyGrid, gzGrid, gMax, izMax, isym, &
       xMom, yMom, zMom, totMom)
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: isym
    REAL(8), INTENT(IN) :: gxGrid, gyGrid, gzGrid, gMax
    INTEGER, INTENT(IN) :: izMax
    REAL(8), DIMENSION(*), INTENT(OUT) :: xMom(:), yMom(:), zMom(:), totMom(:)
    
    INTEGER :: pw, x, y, z, iyMax, ixMax
    REAL(8) :: gx, gy, gz
    
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
                !WRITE(*,*) 'pw#', pw, 'gx=', gx, 'gy=', gy, 'gz=', gz
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





  SUBROUTINE calcTransDipMom(trans, numIncTrans, npw, xMom, yMom, zMom, pwCoef)
    IMPLICIT NONE

    TYPE(includedTransition), DIMENSION(*), INTENT(INOUT) :: trans

    INTEGER, INTENT(IN) :: numIncTrans, npw
    REAL(8), DIMENSION(*), INTENT(IN) :: xMom(:), yMom(:), zMom(:)
    COMPLEX(4), DIMENSION(*), INTENT(IN) :: pwCoef(:,:)

    INTEGER :: i, j

    DO i = 1, numIncTrans
       trans(i)%tdm(1,1) = 0
       trans(i)%tdm(1,2) = 0
       trans(i)%tdm(1,3) = 0
       trans(i)%tdm(2,1) = 0
       trans(i)%tdm(2,2) = 0
       trans(i)%tdm(2,3) = 0
    END DO

    DO i = 1, numIncTrans
       DO j = 1, npw
          trans(i)%tdm(1,1) = trans(i)%tdm(1,1) + & ! Upward x transition
               CONJG(pwCoef(trans(i)%ub,j))*xMom(j)*pwCoef(trans(i)%lb,j) - &
               pwCoef(trans(i)%ub,j)*xMom(j)*CONJG(pwCoef(trans(i)%lb,j))

          trans(i)%tdm(2,1) = trans(i)%tdm(2,1) + & ! Downward x transition
               CONJG(pwCoef(trans(i)%lb,j))*xMom(j)*pwCoef(trans(i)%ub,j) - &
               pwCoef(trans(i)%lb,j)*xMom(j)*CONJG(pwCoef(trans(i)%ub,j))

          trans(i)%tdm(1,2) = trans(i)%tdm(1,2) + &
               CONJG(pwCoef(trans(i)%ub,j))*yMom(j)*pwCoef(trans(i)%lb,j) - &
               pwCoef(trans(i)%ub,j)*yMom(j)*CONJG(pwCoef(trans(i)%lb,j))
          
          trans(i)%tdm(2,2) = trans(i)%tdm(2,2) + &
               CONJG(pwCoef(trans(i)%lb,j))*yMom(j)*pwCoef(trans(i)%ub,j) - &
               pwCoef(trans(i)%lb,j)*yMom(j)*CONJG(pwCoef(trans(i)%ub,j))

          trans(i)%tdm(1,3) = trans(i)%tdm(1,3) + &
               CONJG(pwCoef(trans(i)%ub,j))*zMom(j)*pwCoef(trans(i)%lb,j) - &
               pwCoef(trans(i)%ub,j)*zMom(j)*CONJG(pwCoef(trans(i)%lb,j))
          
          trans(i)%tdm(2,3) = trans(i)%tdm(2,3) + &
               CONJG(pwCoef(trans(i)%lb,j))*zMom(j)*pwCoef(trans(i)%ub,j) - &
               pwCoef(trans(i)%lb,j)*zMom(j)*CONJG(pwCoef(trans(i)%ub,j))
       END DO
    END DO
  END SUBROUTINE calcTransDipMom





  LOGICAL FUNCTION fileExists(fileName)
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN) :: fileName

    LOGICAL :: result

    INQUIRE(file=TRIM(fileName), EXIST=result)
    fileExists = result
  END FUNCTION fileExists





  SUBROUTINE calc_d_bCoefs_at0(t, omega, transition, d_bCoefNOT0, d_bCoef0, &
       bandEnNOT0, bandEn0, bCoefNOT0, bCoef0, upOrDown, polDir)
    IMPLICIT NONE

    REAL(8), INTENT(IN) :: t, omega
    TYPE(includedTransition), INTENT(INOUT) :: transition
    COMPLEX(8), INTENT(INOUT) :: d_bCoefNOT0, d_bCoef0
    REAL(8), INTENT(IN) :: bandEnNOT0, bandEn0
    COMPLEX(8), INTENT(IN) :: bCoefNOT0, bCoef0
    INTEGER, INTENT(IN) :: upOrDown, polDir

    REAL(4) :: l, bSq

    !WRITE(*,*) 'd_bCoef0 = ', d_bCoef0
    !WRITE(*,*) 'q = ', q
    !WRITE(*,*) 'EfieldAmp = ', EfieldAmp
    !WRITE(*,*) 'm = ', m
    !WRITE(*,*) 'bandEn0 = ', bandEn0
    !WRITE(*,*) 'bandEnNOT0 = ', bandEnNOT0
    !WRITE(*,*) 'omega = ', omega
    !WRITE(*,*) 't = ', t
    !WRITE(*,*) 'upOrDown = ', upOrDown
    !WRITE(*,*) 'polDir = ', polDir
    !WRITE(*,*) 'transition%tdm(upOrDown,polDir) = ', transition%tdm(upOrDown,polDir)
    !WRITE(*,*) 'im = ', im
    !WRITE(*,*) 'hbar = ', hbar
    !WRITE(*,*) 'bCoefNOT0 = ', bCoefNOT0
    !WRITE(*,*) 'dt = ', dt
    !WRITE(*,*) 'transition%omega_0 = ', transition%omega_0

    WRITE(*,*) 'omega_0-omega =', transition%omega_0 - omega
    WRITE(*,*) 'sin = ', SIN(dt*(transition%omega_0-omega)/2)
    WRITE(*,*) 'sin^2 = ', SIN(dt*(transition%omega_0-omega)/2)**2
    WRITE(*,*) 'sin^2/stuff^2 = ', (SIN(dt*(transition%omega_0-omega)/2)/(dt*(transition%omega_0-omega)/2))**2
    WRITE(*,*) ''

    d_bCoef0 = d_bCoef0 + (q*EfieldAmp/(m*(bandEn0-bandEnNOT0)))* &
         COS(omega*t)*transition%tdm(upOrDown,polDir)* &
         EXP(im*((bandEn0-bandEnNOT0)/hbar)*t)*bCoefNOT0* &
         (SIN(dt*(transition%omega_0-omega)/2)/(dt*(transition%omega_0-omega)/2))**2
         !SIN(dt*(transition%omega_0-omega)/2)**2

    l = REALPART(bCoefNOT0)
    bSq = REALPART(d_bCoef0*CONJG(d_bCoef0))
    WRITE(*,*) 'bCoefNOT0 = ', bCoefNOT0
    WRITE(*,*) 'l = ', l
    WRITE(*,*) '(dt**2)*|d_bCoef0|^2 = ', (dt**2)*bSq
    WRITE(*,*) ''

    IF (l > 0) THEN
       d_bCoefNOT0 = d_bCoefNOT0 + (-l + SQRT(l**2-(dt**2)*bSq))/dt
    ELSE IF (l < 0) THEN
       d_bCoefNOT0 = d_bCoefNOT0 + (-l - SQRT(l**2-(dt**2)*bSq))/dt
    END IF
  END SUBROUTINE calc_d_bCoefs_at0





!  SUBROUTINE calc_d_bCoefs(t, omega, transition, polDir, bCoefLow, bCoefHigh, d_bCoefLow, d_bCoefHigh)
!    IMPLICIT NONE
!
!    REAL(8), INTENT(IN) :: t, omega
!    TYPE(includedTransition), INTENT(IN) :: transition
!    INTEGER, INTENT(IN) :: polDir
!    COMPLEX(8), INTENT(IN) :: bCoefLow, bCoefHigh
!    COMPLEX(8), INTENT(OUT) :: d_bCoefLow, d_bCoefHigh
!    
!    d_bCoefLow = d_bCoefLow - (q*EfieldAmp/(m*hbar*transition%omega_0))* &
!         COS(omega*t)*transition%tdm(2,polDir)*EXP(-im*transition%omega_0*t)* &
!         bCoefHigh*(SIN(dt*(transition%omega_0-omega)/2)/(dt*(transition%omega_0-omega)/2))**2
!
!    d_bCoefHigh = d_bCoefHigh + (q*EfieldAmp/(m*hbar*transition%omega_0))* &
!         COS(omega*t)*transition%tdm(1,polDir)*EXP(im*transition%omega_0*t)* &
!         bCoefLow*(SIN(dt*(transition%omega_0-omega)/2)/(dt*(transition%omega_0-omega)/2))**2
  !  END SUBROUTINE calc_d_bCoefs





  SUBROUTINE calc_d_bCoefs(nb, npw, t, omega, band, bCoef, bandEn, &
       pertHam, d_pwCoef, pwCoef, d_bCoef)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nb, npw, band
    REAL(8), INTENT(IN) :: t, omega
    COMPLEX(8), DIMENSION(*), INTENT(IN) :: bCoef(:), pertHam(:,:)
    REAL(8), DIMENSION(*), INTENT(IN) :: bandEn(:)
    COMPLEX(8), DIMENSION(*), INTENT(IN) :: d_pwCoef(:,:)
    COMPLEX(4), DIMENSION(*), INTENT(IN) :: pwCoef(:,:)
    COMPLEX(8), INTENT(OUT) :: d_bCoef

    INTEGER :: i, j
    COMPLEX(8) :: be, pwProd

    WRITE(*,'(A,I0.1)') 'band ', band

    d_bCoef = 0

    DO i = 1, nb
       IF (bCoef(i)==0) THEN
          CYCLE
       END IF
       be = bCoef(i)*EXP(im*(bandEn(band)-bandEn(i))*ev_to_au*t/hbar)
       pwProd = 0
       DO j = 1, npw
          pwProd = pwProd + d_pwCoef(i,j)*CONJG(pwCoef(band,j))
       END DO
       d_bCoef = d_bCoef - be*((im/hbar)*pertHam(i,band) + pwProd)
       WRITE(*,'(A,I0.1,A)') 'be(', i, '):'
       WRITE(*,*) be
       WRITE(*,'(A,I0.1,A)') 'pwProd(', i, '):'
       WRITE(*,*) pwProd
       WRITE(*,'(A,I0.1,A)') 'contribution to d_bCoef from band ', i, ':'
       WRITE(*,*) -be*((im/hbar)*pertHam(i,band) + pwProd)
       WRITE(*,*) ''
    END DO
    WRITE(*,'(A,I0.1,A)') 'd_bCoef(', band, '):'
    WRITE(*,*) d_bCoef
    WRITE(*,*) ''
  END SUBROUTINE calc_d_bCoefs





  SUBROUTINE calc_d_bCoefs_2band(nb, npw, t, omega, bCoef, bandEn, &
       pertHam, d_pwCoef, pwCoef, d_bCoef, trans, interBandInProd)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nb, npw
    REAL(8), INTENT(IN) :: t, omega
    COMPLEX(8), DIMENSION(*), INTENT(IN) :: bCoef(:), pertHam(:,:)
    REAL(8), DIMENSION(*), INTENT(IN) :: bandEn(:)
    COMPLEX(8), DIMENSION(*), INTENT(IN) :: d_pwCoef(:,:)
    COMPLEX(4), DIMENSION(*), INTENT(IN) :: pwCoef(:,:)
    COMPLEX(8), DIMENSION(*), INTENT(OUT) :: d_bCoef
    TYPE(includedTransition), INTENT(IN) :: trans
    COMPLEX(8), DIMENSION(*), INTENT(IN) :: interBandInProd(:,:)

    INTEGER :: i, j
    COMPLEX(8) :: be, pwProd_lb, pwProd_ub, pwDiff_lb, pwDiff_ub
    COMPLEX(8) :: s11s22_minus_s12s21
    COMPLEX(8) :: db_lb, db_ub

    WRITE(*,'(A,I0.1,A,I0.1)') 'transition: ', trans%lb, ' <-> ', trans%ub

    s11s22_minus_s12s21 = interBandInProd(trans%lb,trans%lb)* &
         interBandInProd(trans%ub,trans%ub) - &
         interBandInProd(trans%lb,trans%ub)* &
         interBandInProd(trans%ub,trans%lb)

    pwDiff_lb = 0
    pwDiff_ub = 0
    DO j = 1,npw
       pwDiff_lb = CONJG(pwCoef(trans%lb,j))*interBandInProd(trans%ub,trans%ub) - &
            CONJG(pwCoef(trans%ub,j))*interBandInProd(trans%ub,trans%lb)
       pwDiff_ub = CONJG(pwCoef(trans%ub,j))*interBandInProd(trans%lb,trans%lb) - &
            CONJG(pwCoef(trans%lb,j))*interBandInProd(trans%lb,trans%ub)
    END DO
    
    DO i = 1, nb
       d_bCoef(i) = 0
       IF (bCoef(i)==0) THEN
          CYCLE
       END IF
       be = bCoef(i)*EXP(im*(bandEn(trans%lb)-bandEn(i))*ev_to_au*t/hbar)
       pwProd_lb = 0
       pwProd_ub = 0
       DO j = 1, npw
          pwProd_lb = pwProd_lb + d_pwCoef(i,j)*pwDiff_lb
          pwProd_ub = pwProd_ub + d_pwCoef(i,j)*pwDiff_ub
       END DO
       !WRITE(*,'(A,I0.1,A)') 'be(', i, '):'
       !WRITE(*,*) be
       !WRITE(*,'(A,I0.1,A)') 'pwProd(', i, '):'
       !WRITE(*,*) pwProd
       !WRITE(*,'(A,I0.1,A)') 'contribution to d_bCoef from band ', i, ':'
       !WRITE(*,*) -be*((im/hbar)*pertHam(i,band) + pwProd)
    END DO

    db_lb = 0
    db_ub = 0
    db_lb = db_lb - be*pwProd_lb
    WRITE(*,*) 'db_lb after subtracting be*pwProd:'
    WRITE(*,*) db_lb
    WRITE(*,*) ''
    db_lb = db_lb - &
         (im/hbar)*interBandInProd(trans%ub,trans%ub)* &
         EXP(im*(bandEn(trans%lb)-bandEn(trans%ub))*ev_to_au*t/hbar)* &
         bCoef(trans%ub)*pertHam(trans%ub,trans%lb)
    WRITE(*,*) 'db_lb after subtracting stuff*bCoef(upper):'
    WRITE(*,*) db_lb
    WRITE(*,*) ''
    db_lb = db_lb + &
         (im/hbar)*bCoef(trans%lb)*pertHam(trans%lb,trans%ub)* &
         interBandInProd(trans%ub,trans%lb)
    WRITE(*,*) 'db_lb after adding stuff*bCoef(lower):'
    WRITE(*,*) db_lb
    WRITE(*,*) ''
    db_lb = db_lb/s11s22_minus_s12s21
    WRITE(*,*) 'db_lb after dividing by nonorthonormality:'
    WRITE(*,*) db_lb
    WRITE(*,*) ''
    d_bCoef(trans%lb) = d_bCoef(trans%lb) + db_lb!* &
         !(SIN(dt*(trans%omega_0-omega)/2)/(dt*(trans%omega_0-omega)/2))**2

    db_ub = db_ub - be*pwProd_ub
    WRITE(*,*) 'im/hbar =', im/hbar
    WRITE(*,*) 'exp =', EXP(im*(bandEn(trans%ub)-bandEn(trans%lb))*ev_to_au*t/hbar)
    WRITE(*,*) 'bCoef(lower) =', bCoef(trans%lb)
    WRITE(*,*) ''
    WRITE(*,*) 'db_ub after subtracting be*pwProd:'
    WRITE(*,*) db_ub
    WRITE(*,*) ''
    db_ub = db_ub - &
         (im/hbar)*&!interBandInProd(trans%lb,trans%lb)* &
         EXP(im*(bandEn(trans%ub)-bandEn(trans%lb))*ev_to_au*t/hbar)* &
         bCoef(trans%lb)*pertHam(trans%lb,trans%ub)
    WRITE(*,*) 'db_ub after subtracting stuff*bCoef(lower):'
    WRITE(*,*) db_ub
    WRITE(*,*) ''
    db_ub = db_ub + &
         (im/hbar)*bCoef(trans%ub)*pertHam(trans%ub,trans%lb)* &
         interBandInProd(trans%lb,trans%ub)
    WRITE(*,*) 'db_ub after adding stuff*bCoef(upper):'
    WRITE(*,*) db_ub
    WRITE(*,*) ''
    db_ub = db_ub/s11s22_minus_s12s21
    WRITE(*,*) 'db_ub after dividing by nonorthonormality:'
    WRITE(*,*) db_ub
    WRITE(*,*) ''
    d_bCoef(trans%ub) = d_bCoef(trans%ub) + db_ub!* &
         !(SIN(dt*(trans%omega_0-omega)/2)/(dt*(trans%omega_0-omega)/2))**2
    WRITE(*,*) ''
    WRITE(*,'(A,I0.1,A)') 'd_bCoef(', trans%lb, '):'
    WRITE(*,*) d_bCoef(trans%lb)
    WRITE(*,'(A,I0.1,A)') 'd_bCoef(', trans%ub, '):'
    WRITE(*,*) d_bCoef(trans%ub)
    WRITE(*,*) ''
  END SUBROUTINE calc_d_bCoefs_2band





  SUBROUTINE calc_d_pwCoefs(dt, nb, npw, pwCoefOld, pwCoefNew, d_pwCoef)
    IMPLICIT NONE

    REAL(8), INTENT(IN) :: dt
    INTEGER, INTENT(IN) :: nb, npw
    COMPLEX(4), DIMENSION(*), INTENT(IN) :: pwCoefOld(:,:), pwCoefNew(:,:)
    COMPLEX(8), DIMENSION(*), INTENT(OUT) :: d_pwCoef(:,:)

    INTEGER :: i, j

    DO i = 1, nb
       DO j = 1, npw
          d_pwCoef(i,j) = (pwCoefNew(i,j) - pwCoefOld(i,j))/dt
          !WRITE(*,'(A,I0.1,A,I0.1,A)') 'd_pwCoef(', i, ',', j, '):'
          !WRITE(*,*) d_pwCoef(i,j)
       END DO
    END DO
  END SUBROUTINE calc_d_pwCoefs





  SUBROUTINE calc_pertHam(numIncTrans, EfieldAmp, omega, t, bandEn, &
       trans, polDir, pertHam)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: numIncTrans
    REAL(8), INTENT(IN) :: EfieldAmp, omega, t
    REAL(8), DIMENSION(*), INTENT(IN) :: bandEn
    TYPE(includedTransition), DIMENSION(*), INTENT(IN) :: trans
    INTEGER, INTENT(IN) :: polDir
    COMPLEX(8), DIMENSION(*), INTENT(OUT) :: pertHam(:,:)

    INTEGER :: i

    DO i = 1, numIncTrans
       WRITE(*,*) 'q =', q
       WRITE(*,*) 'im =', im
       WRITE(*,*) 'EfieldAmp =', EfieldAmp
       WRITE(*,*) 'm =', m
       WRITE(*,*) 'omega_0 =', trans(i)%omega_0
       WRITE(*,*) 'cos(wt) =', COS(omega*t)
       WRITE(*,*) 'tdm =', trans(i)%tdm(1,polDir)
       pertHam(trans(i)%lb,trans(i)%ub) = (q*im*EfieldAmp/(m*trans(i)%omega_0))* &
            COS(omega*t)*trans(i)%tdm(1,polDir) ! H'_if
       pertHam(trans(i)%ub,trans(i)%lb) = -(q*im*EfieldAmp/(m*trans(i)%omega_0))* &
            COS(omega*t)*trans(i)%tdm(2,polDir) ! H'_fi
       WRITE(*,'(A,I0.1,A,I0.1,A)') &
            'pertHam(', trans(i)%lb, ',', trans(i)%ub, '):'
       WRITE(*,*) pertHam(trans(i)%lb,trans(i)%ub)
       WRITE(*,'(A,I0.1,A,I0.1,A)') &
            'pertHam(', trans(i)%ub, ',', trans(i)%lb, '):'
       WRITE(*,*) pertHam(trans(i)%ub,trans(i)%lb)
       WRITE(*,*) ''
    END DO
  END SUBROUTINE calc_pertHam





  INTEGER(8) FUNCTION getTime()
    IMPLICIT NONE

    CALL SYSTEM('date +%s%3N > tmpTime.txt')
    OPEN(UNIT=98,&
         FILE='tmpTime.txt',&
         ACTION='read')
    READ(98,*) getTime
    CLOSE(98)
  END FUNCTION getTime




  SUBROUTINE readBandEnsAndPWcoefs(rl, nb, npw, nkp, numIncTrans, numIncBands, &
       incBands, trans, pwCoefOld, pwCoefNew, interBandInProd)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: rl, nb, npw, nkp, numIncTrans, numIncBands
    INTEGER, DIMENSION(*), INTENT(IN) :: incBands(:)
    TYPE(includedTransition), DIMENSION(*), INTENT(INOUT) :: trans(:)
    COMPLEX(4), DIMENSION(*), INTENT(OUT) :: pwCoefOld(:,:), pwCoefNew(:,:)
    COMPLEX(8), DIMENSION(*), INTENT(INOUT) :: interBandInProd(:,:)
    
    INTEGER :: i, j, k
    REAL(8) :: tmpNumPW
    REAL(8), DIMENSION(nkp,3) :: tmpKpts
    REAL(8), DIMENSION(nb) :: tmpBandEn, tmpGweight, tmpOcc
    COMPLEX(8), DIMENSION(nb) :: inProd
    
    REAL(8), DIMENSION(numIncBands) :: absInProd
    INTEGER, DIMENSION(numIncBands) :: bandMap

    OPEN(UNIT=21,&
         FILE='WAVECAR_OLD',&
         FORM='unformatted',&
         ACCESS='direct',&
         STATUS='old',&
         RECL=rl)
    READ(21,REC=3) tmpNumPW, (tmpKpts(1,i), i=1,3),&
         (tmpBandEn(j), tmpGweight(j), tmpOcc(j), j=1,nb)
    !DO i = 4, 3+nb
    !   READ(21,REC=i) (pwCoefOld(i-3,j), j=1,npw)
    !END DO
    ! This should be equivalent but faster to do in memory
    CLOSE(21)

    DO i = 1, nb
       DO j = 1, npw
          pwCoefOld(i,j) = pwCoefNew(i,j)
       END DO
    END DO
    
    OPEN(UNIT=12,&
         FILE='WAVECAR',&
         FORM='unformatted',&
         ACCESS='direct',&
         STATUS='old',&
         RECL=rl)
    DO i = 4, 3+nb
       READ(12,REC=i) (pwCoefNew(i-3,j), j=1,npw)
    END DO

    ! Check for flipped indices of degenerate bands !
    DO i = 1, numIncBands
       DO j = 1, numIncBands
          interBandInProd(i,j) = 0
          DO k = 1, npw
             interBandInProd(i,j) = interBandInProd(i,j) + &
                  pwCoefOld(incBands(i),k)*CONJG(pwCoefNew(incBands(j),k))
          END DO
          WRITE(*,'(A,I0.1,A,I0.1,A,"(",ES23.15E3,",",ES23.15E3,")")') &
               'interBandInProd(', incBands(i), ',', incBands(j), ') = ', interBandInProd(i,j)
          absInProd(j) = ABS(interBandInProd(i,j))
          WRITE(*,'(A,ES23.15E3)') '  Absolute value = ', absInProd(j)
          WRITE(*,*) ''
       END DO
       bandMap(i) = MAXLOC(absInProd, DIM=1)
    END DO
    WRITE(*,'(A)') 'bandMap:'
    WRITE(*,*) bandMap
    DO i = 1, numIncBands
       IF (i /= bandMap(i)) THEN
          WRITE(12,REC=3+incBands(bandMap(i))) (pwCoefNew(incBands(i),j), j=1,npw)
       END IF
    END DO
    DO i = 1, numIncBands
       READ(12,REC=3+incBands(i)) (pwCoefNew(incBands(i),j), j=1,npw)
    END DO
    WRITE(*,*) ''

    ! Check for flipped signs on plane wave expansion coefficients !
    DO i = 1, nb
       inProd(i) = 0
       DO j = 1, npw
          inProd(i) = inProd(i) + pwCoefOld(i,j)*CONJG(pwCoefNew(i,j))
       END DO
       WRITE(*,'(A,I0.1,A,ES23.15E3,",",ES23.15E3,A)') &
            'inner product ', i, ' = (', inProd(i), ')'
       IF (REALPART(inProd(i)) < -0.5) THEN
          WRITE(*,'(A)') '  flipping signs of plane waves'
          inProd(i) = 0
          DO j = 1, npw
             pwCoefNew(i,j) = -(1,0)*pwCoefNew(i,j)
             inProd(i) = inProd(i) + pwCoefOld(i,j)*CONJG(pwCoefNew(i,j))
          END DO
          WRITE(12,REC=i+3) (pwCoefNew(i,j), j=1,npw)
          WRITE(*,'(A,I0.1,A,ES23.15E3,",",ES23.15E3,A)') &
               '  new inner product ', i, ' = (', inProd(i), ')'
       END IF
    END DO
    WRITE(*,*) ''

    CLOSE(12)

    DO i = 1, numIncTrans
       trans(i)%omega_0 = ((tmpBandEn(trans(i)%ub) - tmpBandEn(trans(i)%lb))/hbar)*eV_to_au
    END DO
  END SUBROUTINE readBandEnsAndPWcoefs




  SUBROUTINE readBandEnsAndPWcoefsDOWN(rl, nb, npw, nkp, numIncTrans, numIncBands, &
       incBands, trans, pwCoefOld, pwCoefNew, interBandInProd)
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: rl, nb, npw, nkp, numIncTrans, numIncBands
    INTEGER, DIMENSION(*), INTENT(IN) :: incBands(:)
    TYPE(includedTransition), DIMENSION(*), INTENT(INOUT) :: trans(:)
    COMPLEX(4), DIMENSION(*), INTENT(OUT) :: pwCoefOld(:,:), pwCoefNew(:,:)
    COMPLEX(8), DIMENSION(*), INTENT(INOUT) :: interBandInProd(:,:)
    
    INTEGER :: i, j, k
    REAL(8) :: tmpNumPW
    REAL(8), DIMENSION(nkp,3) :: tmpKpts
    REAL(8), DIMENSION(nb) :: tmpBandEn, tmpGweight, tmpOcc
    COMPLEX(8), DIMENSION(nb) :: inProd
    

    OPEN(UNIT=21,&
         FILE='WAVECAR_OLD',&
         FORM='unformatted',&
         ACCESS='direct',&
         STATUS='old',&
         RECL=rl)
    READ(21,REC=spinDownRec) tmpNumPW, (tmpKpts(1,i), i=1,3),&
         (tmpBandEn(j), tmpGweight(j), tmpOcc(j), j=1,nb)
    !DO i = 4, 3+nb
    !   READ(21,REC=i) (pwCoefOld(i-3,j), j=1,npw)
    !END DO
    ! This should be equivalent but faster if done in memory
    CLOSE(21)

    DO i = 1, nb
       DO j = 1, npw
          pwCoefOld(i,j) = pwCoefNew(i,j)
       END DO
    END DO
    
    OPEN(UNIT=12,&
         FILE='WAVECAR',&
         FORM='unformatted',&
         ACCESS='direct',&
         STATUS='old',&
         RECL=rl)
    spinDownRec = 4 + nb
    DO i = 5+nb, 4+nb*2
       READ(12,REC=i) (pwCoefNew(i-4-nb,j), j=1,npw)
    END DO

    ! Check for flipped indices of degenerate bands !
    DO i = 1, numIncBands
       DO j = 1, numIncBands
          interBandInProd(i,j) = 0
          DO k = 1, npw
             interBandInProd(i,j) = interBandInProd(i,j) + &
                  pwCoefOld(incBands(i),k)*CONJG(pwCoefNew(incBands(j),k))
          END DO
          WRITE(*,'(A,I0.1,A,I0.1,A,"(",ES23.15E3,",",ES23.15E3,")")') &
               'interBandInProdDown(', incBands(i), ',', incBands(j), ') = ', interBandInProd(i,j)
       END DO
    END DO
    
    ! Check for flipped signs on plane wave expansion coefficients !
    DO i = 1, nb
       inProd(i) = 0
       DO j = 1, npw
          inProd(i) = inProd(i) + pwCoefOld(i,j)*CONJG(pwCoefNew(i,j))
       END DO
       WRITE(*,'(A,I0.1,A,ES23.15E3,",",ES23.15E3,A)') &
            'inner product ', i, ' = (', inProd(i), ')'
       IF (REALPART(inProd(i)) < -0.5) THEN
          WRITE(*,'(A)') '  flipping signs of plane waves'
          inProd(i) = 0
          DO j = 1, npw
             pwCoefNew(i,j) = -(1,0)*pwCoefNew(i,j)
             inProd(i) = inProd(i) + pwCoefOld(i,j)*CONJG(pwCoefNew(i,j))
          END DO
          WRITE(12,REC=i+3) (pwCoefNew(i,j), j=1,npw)
          WRITE(*,'(A,I0.1,A,ES23.15E3,",",ES23.15E3,A)') &
               '  new inner product ', i, ' = (', inProd(i), ')'
       END IF
    END DO
    WRITE(*,*) ''
    
    CLOSE(12)
    
    DO i = 1, numIncTrans
       trans(i)%omega_0 = ((tmpBandEn(trans(i)%ub) - tmpBandEn(trans(i)%lb))/hbar)*eV_to_au
    END DO
  END SUBROUTINE readBandEnsAndPWcoefsDOWN





  SUBROUTINE convertBandCoefs(nb, npw, pwCoefOld, pwCoefNew, bCoef)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nb, npw
    COMPLEX(4), DIMENSION(*), INTENT(IN) :: pwCoefOld(:,:), pwCoefNew(:,:)
    COMPLEX(8), DIMENSION(*), INTENT(INOUT) :: bCoef

    INTEGER :: i, j, k
    COMPLEX(8), DIMENSION(nb) :: bCoefOld

    DO i = 1, nb
       bCoefOld(i) = bCoef(i)
       bCoef(i) = 0
    END DO

    DO i = 1, nb
       DO j = 1, nb
          DO k = 1, npw
             bCoef(i) = bCoef(i) + bCoefOld(j)*CONJG(pwCoefNew(i,k))*pwCoefOld(j,k)
          END DO
       END DO
       WRITE(*,'(A,I0.1,A,"(",ES23.15E3,",",ES23.15E3,")")') &
            'bCoefOld(', i, ') = ', bCoefOld(i)
       WRITE(*,'(A,I0.1,A,"(",ES23.15E3,",",ES23.15E3,")")') &
            'bCoefNew(', i, ') = ', bCoef(i)
    END DO
  END SUBROUTINE convertBandCoefs

  


    
END PROGRAM UnnamedJank
