      Program test
        REAL NG1, NG2, dy
        integer M
        NG1 = 0.0502
        NG2 = 0.00653
        M = 1

        CALL BILINEAR(NG1, NG2, dy, 1)

        print *, 'dy = ', dy

      End program test


      SUBROUTINE BILINEAR(XVAL, YVAL, dy, M)
      INTEGER M, N
      REAL*8 X(200), Y(200), F(200,200), dy
      REAL*8 XVAL, YVAL 
      INTEGER I, J, k
      REAL*8 A, B
      REAL*8 g00, g01, g10, g11
      REAL*8 g0, g1
      character(len=70) :: filename, filename1, I_str, I_str1
      CHARACTER*72 line, dummy
      real*8, dimension(7) :: data_array, data_array1
      CHARACTER(len=20) :: filename_patterns(63)

      N = 200

      DO I = 1, N
        X(I) = (I-1)*10.0/199
      ENDDO
      DO J = 1, N
        Y(J) = (J-1)*0.1/199
      ENDDO

      I = 1
      DO WHILE (XVAL > X(I+1) .AND. I < N-1)
        I = I + 1
      END DO
      J = 1
      DO WHILE (YVAL > Y(J+1) .AND. J < N-1)
        J = J + 1
      END DO

      !print *, 'I,J =', I, J
      write(I_str, '(I0)') I+99
      write(I_str1, '(I0)') I+100

C------------------------filename-----------------------------
      
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
      filename_patterns(36) = 'E866_Fe5.dat'
      filename_patterns(37) = 'E866_Fe6.dat'
      filename_patterns(38) = 'E866_Fe7.dat'
      filename_patterns(39) = 'E866_Fe8.dat'
      filename_patterns(40) = 'E866_Fe9.dat'
      filename_patterns(41) = 'E866_Fe10.dat'
      filename_patterns(42) = 'E866_Fe11.dat'
      filename_patterns(43) = 'E866_Fe12.dat'
      filename_patterns(44) = 'E866_Fe13.dat'
      filename_patterns(45) = 'E866_Fe14.dat'

      filename_patterns(46) = 'E866_W1.dat'
      filename_patterns(47) = 'E866_W2.dat'
      filename_patterns(48) = 'E866_W3.dat'
      filename_patterns(49) = 'E866_W4.dat'
      filename_patterns(50) = 'E866_W5.dat'
      filename_patterns(51) = 'E866_W6.dat'
      filename_patterns(52) = 'E866_W7.dat'
      filename_patterns(53) = 'E866_W8.dat'
      filename_patterns(54) = 'E866_W9.dat'
      filename_patterns(55) = 'E866_W10.dat'
      filename_patterns(56) = 'E866_W11.dat'
      filename_patterns(57) = 'E866_W12.dat'
      filename_patterns(58) = 'E866_W13.dat'
      filename_patterns(59) = 'E866_W14.dat'

      filename_patterns(60) = 'RHIC_1.dat'
      filename_patterns(61) = 'RHIC_2.dat'
      filename_patterns(62) = 'RHIC_3.dat'
      filename_patterns(63) = 'RHIC_4.dat'
 
      filename = '../log_grid/rep_' // trim(I_str) //
     &              trim('/grid_data/') // trim(filename_patterns(M))
      filename1 = '../log_grid/rep_' // trim(I_str1) // 
     &              trim('/grid_data/') // trim(filename_patterns(M))

C---------------------------------Openfile----------------

      OPEN (UNIT = 1, FILE = filename)
      READ(1,'(A)') LINE    ! skip first row

        DO K = 1,J

          IF ( ((M .GT. 0) .AND. (M .LT. 17)) .OR. 
     &         (M .GT. 31) ) THEN
            read(1,*) data_array, dummy
          ELSEIF ( (M .GT. 16) .AND. (M .LT. 32)) THEN
            read(1,*) data_array1
          ENDIF

        ENDDO

          IF ( ((M .GT. 0) .AND. (M .LT. 17)) .OR.
     &         (M .GT. 31) ) THEN
          F(I,J) = data_array(6)
          read(1,*) data_array, dummy
          F(I,J+1) = data_array(6)
        ELSEIF ( (M .GT. 16) .AND. (M .LT. 32)) THEN
          F(I,J) = data_array1(3)
          read(1,*) data_array1
          F(I,J+1) = data_array1(3)
        ENDIF

      CLOSE(1)

      OPEN (UNIT = 2, FILE = filename1)
      READ(2,'(A)') LINE    ! skip first row
        DO K = 1,J
          IF ( ((M .GT. 0) .AND. (M .LT. 17)) .OR.
     &         (M .GT. 31) ) THEN
            read(2,*) data_array, dummy
          ELSEIF ( (M .GT. 16) .AND. (M .LT. 32)) THEN
            read(2,*) data_array1
          ENDIF
        ENDDO

          IF ( ((M .GT. 0) .AND. (M .LT. 17)) .OR.
     &         (M .GT. 31) ) THEN
            F(I,J) = data_array(6)
            read(2,*) data_array, dummy
            F(I,J+1) = data_array(6)
          ELSEIF ( (M .GT. 16) .AND. (M .LT. 32)) THEN
            F(I,J) = data_array1(3)
          read(2,*) data_array1
            F(I,J+1) = data_array1(3)
          ENDIF

      CLOSE(2)

C-------------------------interpolation--------------------

      A = (XVAL - X(I)) / (X(I+1) - X(I))
      B = (YVAL - Y(J)) / (Y(J+1) - Y(J))

      !print *, 'X(I) = ', X(1)
      !print *, 'X(I+1) = ', X(2)
      !print *, 'X(I) = ', Y(2)

      g00 = F(I,J)
      g01 = F(I,J+1)
      g10 = F(I+1,J)
      g11 = F(I+1,J+1)

      !print *, 'g00 = ', g00
      !print *, 'g10 = ', g10
      !print *, 'g01 = ', g01
      !print *, 'g11 = ', g11

      g0 = g00*(1-A)+g10*A
      g1 = g01*(1-A)+g11*A

      !print *, 'g0 = ', g0
      !print *, 'g1 = ', g1
      
      dy = g0*(1-B)+g1*B

      RETURN
      END

