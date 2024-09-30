!Use gfortran -o test test.f fDSS_HESSIAN17.f to compile

      PROGRAM TestDSS17
      IMPLICIT NONE
    
      INTEGER :: IS, IL, IC, IO
      REAL*8 :: x, Q2
      REAL*8 :: U, UB, D, DB, S, SB, C, B, GL
      REAL*8 :: xmin, xmax, dx
      INTEGER :: FINI17
      INTEGER :: outFile

      COMMON / FRAGINI17 / FINI17
      FINI17 = 0

      IS = 0      ! Best Fit
      IL = 1      ! 68% 
      IC = 1      ! K+
      IO = 1      ! NLO 
      Q2 = 10.0D0 

      xmin = 0.06D0
      xmax = 0.99D0
      dx = 0.01D0  

      OPEN(UNIT=outFile, FILE='nFF_u.csv', STATUS='REPLACE')

      WRITE(outFile, '(A)') 'X,U'

      x = xmin
      DO WHILE (x <= xmax)
        CALL fDSSH(IS, IL, IC, IO, x, Q2, U, UB, D, DB, S, SB, C, B, GL)

        WRITE(outFile, '(F6.4,1X,F12.6)') x, U+UB

        x = x + dx
      END DO

      CLOSE(outFile)

      PRINT *, 'Finished nFF_u.csv'

      END PROGRAM TestDSS17

