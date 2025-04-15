      PROGRAM UpdateWAVECAR

      IMPLICIT NONE
      
      INTEGER :: i, j, k, rl, rl_OLD, nb, npw, ios
      INTEGER :: ix, iy, iz, ixMax, iyMax, izMax, pw
      REAL(8) :: recordLength, recordLength_OLD, ispin, ispin_OLD
      REAL(8) :: singleOrDouble, singleOrDouble_OLD
      REAL(8) :: numKPoints, numKPoints_OLD
      REAL(8) :: numBands, numBands_OLD, enMax, enMax_OLD
      REAL(8) :: eFermi, eFermi_OLD, numPlaneWaves, numPlaneWaves_OLD
      REAL(8) :: dt, timeStep, t, time, omega, E
      REAL(8) :: omega_0, flag
      REAL(8) :: phaseE
      REAL(8) :: gxGrid, gyGrid, gzGrid
      REAL(8) :: gx, gy, gz, gMax
      REAL(8), DIMENSION(3,3) :: latticeVector, latticeVector_OLD
      REAL(8), DIMENSION(3) :: kPts, kPts_OLD
      REAL(8), DIMENSION(2,5) :: rk
      REAL(8), allocatable :: bandEn(:,:), bandEn_OLD(:,:)
      REAL(8), allocatable :: gweight(:,:), gweight_OLD(:,:)
      REAL(8), allocatable :: occ(:,:), occ_OLD(:,:), newOcc(:,:)
      REAL(8), allocatable :: xMom(:), yMom(:), zMom(:), totMom(:)
      COMPLEX(4), allocatable :: pwCoef(:,:), pwCoef_OLD(:,:) !band,PW
      COMPLEX(8), allocatable :: bCoef(:), new_bCoef(:), prev_bCoef(:)
      COMPLEX(8), allocatable :: d_bCoef(:), d2_bCoef(:)
      COMPLEX(8), allocatable :: px(:,:), py(:,:), pz(:,:)
      COMPLEX(8), allocatable :: innerProduct(:)
      COMPLEX(8), DIMENSION(10) :: temp, oldTemp
      LOGICAL :: RWA
      INTEGER :: itScheme
      INTEGER, DIMENSION(2) :: bandArray
      CHARACTER(72) :: dipoleFmt


      REAL(8), PARAMETER :: pi = 3.1415926
      COMPLEX(4), PARAMETER :: im = (0,1)
      REAL(4), PARAMETER :: hbar = 1
      REAL(4), PARAMETER :: q = -1
      REAL(4), PARAMETER :: m = 1




!!!!!!!!!! Placeholder for parameters !!!!!!!!!!

      E = 0.000001
      omega = 0.38
      RWA = .FALSE.
      timeStep = 0.01
      itScheme = 1
      bandArray(1) = 1
      bandArray(2) = 2





!!!!!!!!!! Read in the time info and band coefficients !!!!!!!!!!

      OPEN(UNIT=5,
     $     FILE='dataFile.txt',
     $     ACCESS='sequential',
     $     STATUS='old',
     $     ACTION='read',
     $     FORM='formatted')

!      READ(5,*) timeStep, temp(1), temp(2)
!      PRINT*, 'timeStep=', timeStep
!      time = timeStep
      
      flag = 1
      temp(1) = 1
      temp(2) = 0
      oldTemp(3) = 1
      oldTemp(4) = 0
      READ(5,*)                 ! Skip the header line
      DO WHILE (flag.GE.0)
         oldTemp(1) = oldTemp(3)
         oldTemp(2) = oldTemp(4)
         oldTemp(3) = temp(1)
         oldTemp(4) = temp(2)
         READ(5,*,IOSTAT=ios) flag, temp(1), temp(2)
!         PRINT*, 'flag=', flag, 'IOSTAT=', ios
         IF (ios.EQ.0) THEN
            time = flag
         ELSE
            flag = -1
         END IF
      END DO
      time = time + timeStep
      WRITE(*,'(A,F18.6)') 'time = ', time
      WRITE(*,'(A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $     'temp(1) = ', temp(1)
      WRITE(*,'(A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $     'temp(2) = ', temp(2)
      WRITE(*,'(A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $     'oldTemp(1) = ', oldTemp(1)
      WRITE(*,'(A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $     'oldTemp(2) = ', oldTemp(2)
      
      CLOSE(5)

      

      
!!!!!!!!!! Find record lengths !!!!!!!!!!
      
      OPEN(UNIT=12,
     $     FILE='WAVECAR',
     $     FORM='unformatted',
     $     ACCESS='direct',
     $     STATUS='old',
     $     RECL=8)

      OPEN(UNIT=21,
     $     FILE='WAVECAR_OLD',
     $     FORM='unformatted',
     $     ACCESS='direct',
     $     STATUS='old',
     $     RECL=8)
      
      READ(12,REC=1) recordLength
      rl = INT(recordLength)

      READ(21,REC=1) recordLength_OLD
      rl_OLD = INT(recordLength_OLD)

      CLOSE(12)

      CLOSE(21)

      IF (rl/=rl_OLD) THEN
         PRINT*, "Mismatching WAVECAR record lengths between time steps"
         PRINT*, "Terminating updateWAVECAR.exe"
         STOP 2
      END IF
      
      
      
      

!!!!!!!!!! Open WAVECARS for data extraction !!!!!!!!!!

      OPEN(UNIT=12,
     $     FILE='WAVECAR',
     $     FORM='unformatted',
     $     ACCESS='direct',
     $     STATUS='old',
     $     RECL=rl)

      OPEN(UNIT=21,
     $     FILE='WAVECAR_OLD',
     $     FORM='unformatted',
     $     ACCESS='direct',
     $     STATUS='old',
     $     RECL=rl_OLD)
         
      



!!!!!!!!!! Extract metadata from WAVECARS !!!!!!!!!!
      
      READ(12,REC=1) recordLength, ispin, singleOrDouble
      READ(12,REC=2) numKPoints, numBands, enMax,
     $     ((latticeVector(i,j),i=1,3),j=1,3),
     $     eFermi


      READ(21,REC=1) recordLength_OLD, ispin_OLD, singleOrDouble_OLD
      READ(21,REC=2) numKPoints_OLD, numBands_OLD, enMax_OLD,
     $     ((latticeVector_OLD(i,j),i=1,3),j=1,3),
     $     eFermi_OLD
      
      

      
      
!!!!!!!!!! Extract data from WAVECARS !!!!!!!!!!

      IF (INT(numBands)==INT(numBands_OLD)) THEN
         nb = INT(numBands)
         ALLOCATE(bandEn(nb,2))
         ALLOCATE(bandEn_OLD(nb,2))
         ALLOCATE(gweight(nb,2))
         ALLOCATE(gweight_OLD(nb,2))
         ALLOCATE(occ(nb,2))
         ALLOCATE(occ_OLD(nb,2))
         ALLOCATE(newOcc(nb,2))
         ALLOCATE(bCoef(nb))
         ALLOCATE(new_bCoef(nb))
         ALLOCATE(prev_bCoef(nb))
         ALLOCATE(px(nb,nb))
         ALLOCATE(py(nb,nb))
         ALLOCATE(pz(nb,nb))
         ALLOCATE(d_bCoef(nb))
         ALLOCATE(innerProduct(nb))
      ELSE
         CLOSE(12)
         CLOSE(21)
         PRINT*, "Mismatching number of bands between time steps!"
         PRINT*, "Terminating updateWAVECAR.exe"
         STOP 2
      END IF

!!!!! Spin-up !!!!!

      READ(12,REC=3) numPlaneWaves, (kPts(i),i=1,3),
     $     (bandEn(j,1),gweight(j,1),occ(j,1),j=1,nb)

      READ(21,REC=3) numPlaneWaves_OLD, (kPts_OLD(i),i=1,3),
     $     (bandEn_OLD(j,1),gweight_OLD(j,1),occ_OLD(j,1),j=1,nb)

      IF (INT(numPlaneWaves)==INT(numPlaneWaves_OLD)) THEN
         npw = INT(numPlaneWaves)
         ALLOCATE(pwCoef(nb,npw))
         ALLOCATE(pwCoef_OLD(nb,npw))
         ALLOCATE(xMom(npw))
         ALLOCATE(yMom(npw))
         ALLOCATE(zMom(npw))
         ALLOCATE(totMom(npw))
      ELSE
         CLOSE(12)
         CLOSE(21)
         PRINT*, "Mismatching number of PWs between time steps!"
         PRINT*, "Terminating updateWAVECAR.exe"
         STOP 2
      END IF


      PRINT*, "Assigning bCoef variables:"
      bCoef(bandArray(1)) = temp(1)
      bCoef(bandArray(2)) = temp(2)
      prev_bCoef(bandArray(1)) = oldTemp(1)
      prev_bCoef(bandArray(2)) = oldTemp(2)
      WRITE(*,'(A,I0.1,A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $     'bCoef(',bandArray(1),') = ', bCoef(bandArray(1))
      WRITE(*,'(A,I0.1,A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $     'bCoef(',bandArray(2),') = ', bCoef(bandArray(2))
      WRITE(*,'(A,I0.1,A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $     'prev_bCoef(',bandArray(1),') = ', prev_bCoef(bandArray(1))
      WRITE(*,'(A,I0.1,A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $     'prev_bCoef(',bandArray(2),') = ', prev_bCoef(2)
      WRITE(*,'(A,I0.1,A,ES23.15E3)')
     $     'occ(',bandArray(2),',1) = ', occ(bandArray(2),1)

      DO i=4,3+nb
         READ(12,REC=i) (pwCoef(i-3,j),j=1,npw)
         READ(21,REC=i) (pwCoef_OLD(i-3,j),j=1,npw)
      END DO

      DO i=1,nb
         innerProduct(i) = 0
         DO j=1,npw
            innerProduct(i) = innerProduct(i)
     $           + pwCoef_OLD(i,j)*CONJG(pwCoef(i,j))
         END DO
         WRITE(*,'(A,I0.1,A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $        'innerProduct ', i, ' = ', innerProduct(i)
         IF (REALPART(innerProduct(i))<-0.5) THEN
            PRINT*, "Flipping signs..."
            innerProduct(i)=0
            DO j=1,npw
               pwCoef(i,j) = -(1,0)*pwCoef(i,j)
               innerProduct(i) = innerProduct(i)
     $              + pwCoef_OLD(i,j)*CONJG(pwCoef(i,j))
            END DO
            WRITE(12,REC=i+3) (pwCoef(i,j),j=1,npw)
            PRINT*, "New innerProduct", i, "=", innerProduct(i)
         END IF
      END DO

      temp(3)=0
      temp(4)=0
      DO i=1,npw
         temp(3)=temp(3)+
     $        pwCoef(bandArray(1),i)*CONJG(pwCoef(bandArray(1),i))
         temp(4)=temp(4)+
     $        pwCoef(bandArray(2),i)*CONJG(pwCoef(bandArray(2),i))
      END DO
      WRITE(*,'(A,I0.1,A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $     'Band ',bandArray(1),' Occ. = ',
     $     bCoef(bandArray(1))*CONJG(bCoef(bandArray(1)))
      WRITE(*,'(A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $     'by PW = ', temp(3)
      WRITE(*,'(A,I0.1,A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $     'Band ',bandArray(2),' Occ. = ',
     $     bCoef(bandArray(2))*CONJG(bCoef(bandArray(2)))
      WRITE(*,'(A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $     'by PW = ', temp(4)

!!!!! Spin-down !!!!!

!# I'll worry about this later #!
      
      
      CLOSE(12)
      CLOSE(21)





!!!!!!!!!! Calculate momentum for each plane wave !!!!!!!!!!

      !!!!! Currently assuming rectangular simulation cell !!!!!
      gxGrid = 2*pi/(latticeVector(1,1)*1.8897259886)
      gyGrid = 2*pi/(latticeVector(2,2)*1.8897259886)
      gzGrid = 2*pi/(latticeVector(3,3)*1.8897259886)

      gMax = DSQRT(2*(enMax*0.036749405469679))

      izMax = INT(gMax/gzGrid)
      
      pw = 0
      DO i = 0, 2*izMax
         iz = i-(2*izMax+1)*INT(i/(izMax+1)) ! z=0,1,...,zMax,-zMax,...,-1
         gz = iz*gzGrid
         iyMax = INT(DSQRT(gMax**2-gz**2)/gyGrid)
         DO j = 0, 2*iyMax
            iy = j-(2*iyMax+1)*INT(j/(iyMax+1))
            gy = iy*gyGrid
            ixMax = INT(DSQRT(gMax**2-gz**2-gy**2)/gxGrid)
            DO k = 0,ixMax
               ix = k
               IF (ix==0 .AND. iy<0) CYCLE
               IF (ix==0 .AND. iy==0 .AND. iz<0) CYCLE
               gx = ix*gxGrid
               pw = pw+1
               xMom(pw) = gx
               yMom(pw) = gy
               zMom(pw) = gz
               totMom(pw) = DSQRT(gx**2+gy**2+gz**2)
               !PRINT*, 'pw#', pw, 'gx=', gx, 'gy=', gy, 'gz=', gz
               !PRINT*, 'pw#', pw, 'total momentum =', totMom(pw)
            END DO
         END DO
      END DO      




!!!!!!!!!! Calculate new orbital occupations !!!!!!!!!!

! fs to au !
      dt = timeStep*41.341374575751
      t = time*41.341374575751

      
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
      WRITE(*,'(A,F18.12)') 'omega_0 = ', omega_0

      phaseE = DCOS(omega*t)
      WRITE(*,'(A,F18.12)') 'cos(omega*t) = ', phaseE
      
      WRITE(*,'(A,I0.1,A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $     'bCoef(',bandArray(1),') = ', bCoef(bandArray(1))
      WRITE(*,'(A,I0.1,A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $     'bCoef(',bandArray(2),') = ', bCoef(bandArray(2))
      WRITE(*,'(A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $     'exp = ', EXP(-im*omega_0*t)
      WRITE(*,'(A,I0.1,A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $     'exp*bCoef(',bandArray(2),') = ',
     $     EXP(-im*omega_0*t)*bCoef(bandArray(2))

      IF (RWA) THEN
         d_bCoef(bandArray(1)) = -((q*E)/(2*m*hbar*omega_0))*
     $        EXP(im*(omega-omega_0)*t)*
     $        pz(bandArray(1),bandArray(2))*bCoef(bandArray(2))
         d_bCoef(bandArray(2)) = ((q*E)/(2*m*hbar*omega_0))*
     $        EXP(im*(omega_0-omega)*t)*
     $        pz(bandArray(2),bandArray(1))*bCoef(bandArray(1))
      ELSE
         d_bCoef(bandArray(1)) = -(q*E/(m*hbar*omega_0))*
     $        DCOS(omega*t)*pz(bandArray(1),bandArray(2))*
     $        EXP(-im*omega_0*t)*bCoef(bandArray(2))
         d_bCoef(bandArray(2)) = (q*E/(m*hbar*omega_0))*
     $        DCOS(omega*t)*pz(bandArray(2),bandArray(1))*
     $        EXP(im*omega_0*t)*bCoef(bandArray(1))
      END IF
      WRITE(*,'(A,I0.1,A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $     'd_bCoef(',bandArray(1),') = ', d_bCoef(bandArray(1))
      WRITE(*,'(A,I0.1,A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $     'd_bCoef(',bandArray(2),') = ', d_bCoef(bandArray(2))

      IF (itScheme.EQ.1) THEN
         new_bCoef(bandArray(1)) = bCoef(bandArray(1)) +
     $        dt*d_bCoef(bandArray(1))
         new_bCoef(bandArray(2)) = bCoef(bandArray(2)) +
     $        dt*d_bCoef(bandArray(2))
      ELSE IF (itScheme.EQ.2) THEN
         ALLOCATE(d2_bCoef(nb))
         IF (RWA) THEN
            d2_bCoef(bandArray(1)) =
     $           -(q*E*pz(bandArray(1),bandArray(2))
     $           *EXP(im*(omega-omega_0)*t)*
     $           im*(omega-omega_0)*bCoef(bandArray(2)))/
     $           (2*m*hbar*omega_0) -
     $           ((q**2)*(E**2)*pz(bandArray(1),bandArray(2))
     $           *pz(bandArray(2),bandArray(1))*bCoef(bandArray(1)))/
     $           (4*(m**2)*(hbar**2)*(omega_0**2))
            d2_bCoef(bandArray(2)) =
     $           (q*E*im*(omega_0-omega)*
     $           pz(bandArray(2),bandArray(1))*
     $           EXP(im*(omega_0-omega)*t)*bCoef(bandArray(1)))/
     $           (2*m*hbar*omega_0) -
     $           ((q**2)*(E**2)*pz(bandArray(1),bandArray(2))*
     $           pz(bandArray(2),bandArray(1))*bCoef(bandArray(2)))/
     $           (4*(m**2)*(hbar**2)*(omega_0**2))
         ELSE
            d2_bCoef(bandArray(1)) =
     $           (q*E*pz(bandArray(1),bandArray(2))*
     $           EXP(-im*omega_0*t)*bCoef(bandArray(2))*
     $           (omega*DSIN(omega*t)+im*omega_0*DCOS(omega*t)))/
     $           (m*hbar*omega_0) -
     $           ((q**2)*(E**2)*pz(bandArray(1),bandArray(2))*
     $           pz(bandArray(2),bandArray(1))*
     $           ((DCOS(omega*t))**2)*bCoef(bandArray(1)))/
     $           ((m**2)*(hbar**2)*(omega_0**2))
            d2_bCoef(bandArray(2)) =
     $           -(q*E*pz(bandArray(2),bandArray(1))*
     $           EXP(im*omega_0*t)*bCoef(bandArray(1))*
     $           (omega*DSIN(omega*t)-im*omega_0*DCOS(omega*t)))/
     $           (hbar*m*omega_0) -
     $           ((q**2)*(E**2)*pz(bandArray(1),bandArray(2))*
     $           pz(bandArray(2),bandArray(1))*
     $           ((DCOS(omega*t))**2)*bCoef(bandArray(2)))/
     $           ((m**2)*(hbar**2)*(omega_0**2))
         END IF
         new_bCoef(bandArray(1)) = 2*bCoef(bandArray(1)) -
     $        prev_bCoef(bandArray(1)) + d2_bCoef(bandArray(1))*dt**2
         new_bCoef(bandArray(2)) = 2*bCoef(bandArray(2)) -
     $        prev_bCoef(bandArray(2)) + d2_bCoef(bandArray(2))*dt**2
         WRITE(*,'(A,I0.1,A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $        'new_bCoef(',bandArray(1),') = ', new_bCoef(bandArray(1))
         WRITE(*,'(A,I0.1,A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $        'new_bCoef(',bandArray(2),') = ', new_bCoef(bandArray(2))
         PRINT *, 'End of Verlet code block'
      ELSE IF (itScheme.EQ.3) THEN
         rk(1,1) = dt*d_bCoef(bandArray(1))
         rk(2,1) = dt*d_bCoef(bandArray(2))
         IF (RWA) THEN
            rk(1,2) = dt*((-q*E*EXP(im*(omega-omega_0)*(t+dt))*
     $           pz(bandArray(1),bandArray(2))*
     $           (bCoef(bandArray(2))+rk(2,1)))/
     $           (2*m*hbar*omega_0))
            rk(2,2) = dt*((q*E*EXP(im*(omega_0-omega)*(t+dt))*
     $           pz(bandArray(2),bandArray(1))*
     $           (bCoef(bandArray(1))+rk(1,1)))/
     $           (2*m*hbar*omega_0))
         ELSE
            rk(1,2) = dt*((-q*E*DCOS(omega*(t+dt))*
     $           pz(bandArray(1),bandArray(2))*
     $           EXP(-im*omega_0*(t+dt))*
     $           (bCoef(bandArray(2))+rk(2,1)))/
     $           (m*hbar*omega_0))
            rk(2,2) = dt*((q*E*DCOS(omega*(t+dt))*
     $           pz(bandArray(2),bandArray(1))*
     $           EXP(im*omega_0*(t+dt))*
     $           (bCoef(bandArray(1))+rk(1,1)))/
     $           (m*hbar*omega_0))
         END IF
         rk(1,5) = (rk(1,1)+rk(1,2))/2
         rk(2,5) = (rk(2,1)+rk(2,2))/2
         new_bCoef(bandArray(1)) = bCoef(bandArray(1)) + rk(1,5)
         new_bCoef(bandArray(2)) = bCoef(bandArray(2)) + rk(2,5)
         PRINT *, 'End of Runge-Kutta 2 block'
      ELSE IF (itScheme.EQ.4) THEN
         rk(1,1) = dt*d_bCoef(bandArray(1))
         rk(2,1) = dt*d_bCoef(bandArray(2))
         IF (RWA) THEN
            rk(1,2) = dt*((-q*E*EXP(im*(omega-omega_0)*(t+dt/2))*
     $           pz(bandArray(1),bandArray(2))*
     $           (bCoef(bandArray(2))+rk(2,1)/2))/
     $           (2*m*hbar*omega_0))
            rk(2,2) = dt*((q*E*EXP(im*(omega_0-omega)*(t+dt/2))*
     $           pz(bandArray(2),bandArray(1))*
     $           (bCoef(bandArray(1))+rk(1,1)/2))/
     $           (2*m*hbar*omega_0))
            rk(1,3) = dt*((-q*E*EXP(im*(omega-omega_0)*(t+dt/2))*
     $           pz(bandArray(1),bandArray(2))*
     $           (bCoef(bandArray(2))+rk(2,2)/2))/
     $           (2*m*hbar*omega_0))
            rk(2,3) = dt*((q*E*EXP(im*(omega_0-omega)*(t+dt/2))*
     $           pz(bandArray(2),bandArray(1))*
     $           (bCoef(bandArray(1))+rk(1,2)/2))/
     $           (2*m*hbar*omega_0))
            rk(1,4) = dt*((-q*E*EXP(im*(omega-omega_0)*(t+dt))*
     $           pz(bandArray(1),bandArray(2))*
     $           (bCoef(bandArray(2))+rk(2,3)))/
     $           (2*m*hbar*omega_0))
            rk(2,4) = dt*((q*E*EXP(im*(omega_0-omega)*(t+dt))*
     $           pz(bandArray(2),bandArray(1))*
     $           (bCoef(bandArray(1))+rk(1,3)))/
     $           (2*m*hbar*omega_0))
         ELSE
            rk(1,2) = dt*((-q*E*DCOS(omega*(t+dt/2))*
     $           pz(bandArray(1),bandArray(2))*
     $           EXP(-im*omega_0*(t+dt/2))*
     $           (bCoef(bandArray(2))+rk(2,1)/2))/
     $           (m*hbar*omega_0))
            rk(2,2) = dt*((q*E*DCOS(omega*(t+dt/2))*
     $           pz(bandArray(2),bandArray(1))*
     $           EXP(im*omega_0*(t+dt/2))*
     $           (bCoef(bandArray(1))+rk(1,1)/2))/
     $           (m*hbar*omega_0))
            rk(1,3) = dt*((-q*E*DCOS(omega*(t+dt/2))*
     $           pz(bandArray(1),bandArray(2))*
     $           EXP(-im*omega_0*(t+dt/2))*
     $           (bCoef(bandArray(2))+rk(2,2)/2))/
     $           (m*hbar*omega_0))
            rk(2,3) = dt*((q*E*DCOS(omega*(t+dt/2))*
     $           pz(bandArray(2),bandArray(1))*
     $           EXP(im*omega_0*(t+dt/2))*
     $           (bCoef(bandArray(1))+rk(1,2)/2))/
     $           (m*hbar*omega_0))
            rk(1,4) = dt*((-q*E*DCOS(omega*(t+dt))*
     $           pz(bandArray(1),bandArray(2))*
     $           EXP(-im*omega_0*(t+dt))*
     $           (bCoef(bandArray(2))+rk(2,3)))/
     $           (m*hbar*omega_0))
            rk(2,4) = dt*((q*E*DCOS(omega*(t+dt))*
     $           pz(bandArray(2),bandArray(1))*
     $           EXP(im*omega_0*(t+dt))*
     $           (bCoef(bandArray(1))+rk(1,3)))/
     $           (m*hbar*omega_0))
         END IF
         rk(1,5) = (rk(1,1)+2*rk(1,2)+2*rk(1,3)+rk(1,4))/6
         rk(2,5) = (rk(2,1)+2*rk(2,2)+2*rk(2,3)+rk(2,4))/6
         new_bCoef(bandArray(1)) = bCoef(bandArray(1)) + rk(1,5)
         new_bCoef(bandArray(2)) = bCoef(bandArray(2)) + rk(2,5)
         PRINT*, 'Reached end of Runge-Kutta 4 block'
      ELSE
         PRINT*, 'Invalid iteration scheme detected in updateWAVECAR.f'
         PRINT*, 'Terminating program'
         STOP 2
      END IF
         
      WRITE(*,'(A,I0.1,A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $     'new_bCoef(',bandArray(1),') = ', new_bCoef(bandArray(1))
      WRITE(*,'(A,I0.1,A,"( ",ES23.15E3," , ",ES23.15E3," )")')
     $     'new_bCoef(',bandArray(2),') = ', new_bCoef(bandArray(2))
      
      newOcc(bandArray(1),1) =
     $     CONJG(new_bCoef(bandArray(1)))*new_bCoef(bandArray(1))
      newOcc(bandArray(2),1) =
     $     CONJG(new_bCoef(bandArray(2)))*new_bCoef(bandArray(2))
      WRITE(*,'(A,ES23.15E3,A,ES23.15E3)') 'new Occ.s = ',
     $     newOcc(bandArray(1),1), '     ', newOcc(bandArray(2),1)

      


!!!!!!!!!! Update WAVECAR occupations !!!!!!!!!!

      OPEN(UNIT=12,
     $     FILE='WAVECAR',
     $     FORM='unformatted',
     $     ACCESS='direct',
     $     STATUS='old',
     $     RECL=8)

      !# UPDATE THIS #!
      WRITE(12,REC=rl/4+4+3*bandArray(1)) newOcc(bandArray(1),1)
      WRITE(12,REC=rl/4+4+3*bandArray(2)) newOcc(bandArray(2),1)

      CLOSE(12)




!!!!!!!!!! Create data file !!!!!!!!!!

      OPEN(UNIT=50,
     $     FILE='dataFile.txt',
     $     ACCESS='append')

      WRITE(UNIT=50, FMT=*) time,
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
      
      END PROGRAM UpdateWAVECAR





