C------------
C E772 1-16
C CMS5 17-24
C ATLAS1 25-31
C E886 32-39
C RHIC 40-43
C------------
      Program test
        integer position(1), I, J, K, D
        real chi2(40000), NG1, NG2
        
        CALL BILINEAR(chi2, 17, 31)
        
        print *, 'global minimum of chi2 is ', MINVAL(chi2)
        position = MINLOC(chi2)
        print *, 'number ', position, ' is where minimum is.'
        K = floor(position(1)/200.0 - 0.000001)
        I = K+1
        J = position(1) - 200*K
        NG1 = (I-1)*10/199.0
        NG2 = (J-1)*0.1/199.0
        print *, 'best NG1 is ', NG1
        print *, 'best NG2 is ', NG2

      End program test


      SUBROUTINE BILINEAR(chi2, A, B)
      INTEGER M, A, B
      REAL tmp, chi2(40000)
      
      INTEGER I, J, K
     
      character(len=70) :: filename, I_str 
      CHARACTER*72 line, dummy, filename_patterns(31)
      real, dimension(7) :: data_array(7), data_array1(4)

      tmp = 0

C------------------------------------------------------------

      filename_patterns(1) = 'E772_C1.dat'
      filename_patterns(2) = 'E772_C2.dat'
      filename_patterns(3) = 'E772_C3.dat'
      filename_patterns(4) = 'E772_C4.dat'
      filename_patterns(5) = 'E772_Ca1.dat'
      filename_patterns(6) = 'E772_Ca2.dat'
      filename_patterns(7) = 'E772_Ca3.dat'
      filename_patterns(8) = 'E772_Ca4.dat'
      filename_patterns(9) = 'E772_Fe1.dat'
      filename_patterns(10) = 'E772_Fe2.dat'
      filename_patterns(11) = 'E772_Fe3.dat'
      filename_patterns(12) = 'E772_Fe4.dat'
      filename_patterns(13) = 'E772_W1.dat'
      filename_patterns(14) = 'E772_W2.dat'
      filename_patterns(15) = 'E772_W3.dat'
      filename_patterns(16) = 'E772_W4.dat'

      filename_patterns(17) = 'CMS5_1.dat'
      filename_patterns(18) = 'CMS5_2.dat'
      filename_patterns(19) = 'CMS5_3.dat'
      filename_patterns(20) = 'CMS5_4.dat'
      filename_patterns(21) = 'CMS5_5.dat'
      filename_patterns(22) = 'CMS5_6.dat'
      filename_patterns(23) = 'CMS5_7.dat'
      filename_patterns(24) = 'CMS5_8.dat'

      filename_patterns(25) = 'ATLAS1_1.dat'
      filename_patterns(26) = 'ATLAS1_2.dat'
      filename_patterns(27) = 'ATLAS1_3.dat'
      filename_patterns(28) = 'ATLAS1_4.dat'
      filename_patterns(29) = 'ATLAS1_5.dat'
      filename_patterns(30) = 'ATLAS1_6.dat'
      filename_patterns(31) = 'ATLAS1_7.dat'

      filename_patterns(32) = 'E866_Fe1.dat'
      filename_patterns(33) = 'E866_Fe2.dat'
      filename_patterns(34) = 'E866_Fe3.dat'
      filename_patterns(35) = 'E866_Fe4.dat'
      filename_patterns(36) = 'E866_W1.dat'
      filename_patterns(37) = 'E866_W2.dat'
      filename_patterns(38) = 'E866_W3.dat'
      filename_patterns(39) = 'E866_W4.dat'

      filename_patterns(40) = 'RHIC_1.dat'
      filename_patterns(41) = 'RHIC_2.dat'
      filename_patterns(42) = 'RHIC_3.dat'
      filename_patterns(43) = 'RHIC_4.dat'

C------------------------filename-----------------------------
      OPEN(UNIT=10, FILE='chi2.dat', STATUS='unknown')

      DO I = 1,200

      write(I_str, '(I0)') I+99

      DO J = 1,200 
      DO M = A,B

      filename = '../log_grid/rep_' // trim(I_str) //
     &              trim('/grid_data/') // trim(filename_patterns(M))

      OPEN (UNIT = 1, FILE = filename)
      READ(1,'(A)') LINE    ! skip first row

      DO K = 1,J

        IF ( ((M .GT. 0) .AND. (M .LT. 17)) .OR.
     &       ((M .GT. 31) .AND. (M .LT. 44)) ) THEN
          read(1,*) data_array, dummy
        ELSEIF ( (M .GT. 16) .AND. (M .LT. 32)) THEN
          read(1,*) data_array1
        ENDIF

      ENDDO

        IF ( ((M .GT. 0) .AND. (M .LT. 17)) .OR.
     &       ((M .GT. 31) .AND. (M .LT. 44)) ) THEN
          tmp = tmp + data_array(7)
        ELSEIF ( (M .GT. 16) .AND. (M .LT. 32)) THEN
          tmp = tmp + data_array1(4)
        ENDIF      

      CLOSE(1)

      ENDDO !(M)

      chi2(200*(I-1)+J) = tmp
      tmp = 0      

      ENDDO !(J)

      ENDDO !(I)

      CLOSE(10)
      RETURN
      END

