      PROGRAM test

      IMPLICIT NONE

      INTEGER :: i, sz, rl
      REAL(8) :: temp, rlFloat
      COMPLEX(4) :: compTemp
      
      OPEN(UNIT=1, FILE='WAVECAR', FORM='unformatted',
     $     ACCESS='direct', STATUS='old', RECL=8)

      INQUIRE(UNIT=1, SIZE=sz)
      PRINT*, 'size=', sz

      OPEN(UNIT=20, FILE='testDump')
      
      DO i=1,sz/8
         READ(1,REC=i) compTemp
         WRITE(20,*) 'i=', i, compTemp
      END DO

      

      CLOSE(1)
      CLOSE(20)
      
      END PROGRAM test
