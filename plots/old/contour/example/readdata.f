c-----------------------------------------
c     read data file which has 15 columns
c-----------------------------------------
      SUBROUTINE READDATA15(FILEEXP,NLINES,
     >                    X1,X2,X3,X4 ,X5 ,X6,
     >                    X7,X8,X9,X10,X11,X12,
     >                    X13,X14,X15,
     >                    NTOT)
      CHARACTER*72 FILEEXP,LINE
      REAL*8 X1(300) ,X2(300) ,X3(300)
      REAL*8 X4(300) ,X5(300) ,X6(300)
      REAL*8 X7(300) ,X8(300) ,X9(300)
      REAL*8 X10(300),X11(300),X12(300)
      REAL*8 X13(300),X14(300),X15(300)
      INTEGER NLINES,NTOT,I,NT

      OPEN(unit = 91,FILE=FILEEXP,STATUS='unknown')

      DO I = 1,NLINES,1
         READ(91,'(A)') LINE    ! skip NLINES..
      ENDDO

      NT = 0
 10   NT =  NT +1

      READ(91,*,ERR = 99, end = 11)
     +     X1(NT),X2(NT),X3(NT),X4(NT) ,X5(NT) ,X6(NT),
     +     X7(NT),X8(NT),X9(NT),X10(NT),X11(NT),X12(NT),
     +     X13(NT),X14(NT),X15(NT)

      goto 10
 11   if(NT.GT.NTOT) goto 101

      CLOSE(91)

      GOTO 101
99    PRINT *, 'ERROR IN READING FILE    ', FILEEXP
101   CONTINUE
      RETURN
      END


      SUBROUTINE READDATA
      IMPLICIT NONE
      REAL*8 X1(300),X2(300),X3(300),X4(300),X5(300)
      REAL*8 X6(300),X7(300),X8(300),X9(300),X10(300)
      REAL*8 X11(300),X12(300),X13(300),X14(300),X15(300)
      CHARACTER*72 FILEEXP
      INTEGER I
      INTEGER NLINES ! NUMBER OF LINES SKIPPED
      INTEGER NTOT ! NUMBER OF ROWS
      integer usenc
      common /nckern/ usenc
      real*8 ANA(101), BNA(101), GDA(101), ADA(101), GPA(101), APA(101)
      COMMON /param_dat/ ANA, BNA, GDA, ADA, GPA, APA
      NTOT = 101
      NLINES = 1

      FILEEXP='../../collected_results.txt'

      CALL READDATA9(FILEEXP,NLINES,
     >                    X1, X2, X3, X4, X5, X6, X7, X8, X9, NTOT)

      DO I=1,NTOT
         ANA(I)  =X4(I)
         BNA(I)  =X5(I)
         GDA(I)  =X6(I)
         ADA(I)  =X7(I)
         GPA(I)  =X8(I)
         APA(I)  =X9(I)

      END DO

      RETURN
      END

c-----------------------------------------
c     read data file which has 9 columns
c-----------------------------------------
      SUBROUTINE READDATA9(FILEEXP,NLINES,
     >                    X1,X2,X3,X4,X5,X6,X7,X8,X9,NTOT)
      CHARACTER*72 FILEEXP,LINE
      REAL*8 X1(300),X2(300),X3(300),X4(300),X5(300),X6(300)
      REAL*8 X7(300),X8(300),X9(300)
      INTEGER NLINES,NTOT,I,NT

      OPEN(unit = 91,FILE=FILEEXP,STATUS='unknown')

      DO I = 1,NLINES,1
         READ(91,'(A)') LINE    ! skip NLINES..
      ENDDO

      NT = 0
 10   NT =  NT +1

      READ(91,*,ERR = 99, end = 11)
     +     X1(NT),X2(NT),X3(NT),X4(NT),X5(NT),X6(NT),
     +     X7(NT),X8(NT),X9(NT)

      goto 10
 11   if(NT.GT.NTOT) goto 101

      CLOSE(91)

      GOTO 101
99    PRINT *, 'ERROR IN READING FILE    ', FILEEXP
101   CONTINUE
      RETURN
      END
c-----------------------------------------
c     read data file which has 8 columns
c-----------------------------------------
      SUBROUTINE READDATA8(FILEEXP,NLINES,
     >                    X1,X2,X3,X4,X5,X6,X7,X8,NTOT)
      CHARACTER*72 FILEEXP,LINE
      REAL*8 X1(300),X2(300),X3(300),X4(300),X5(300),X6(300)
      REAL*8 X7(300),X8(300)
      INTEGER NLINES,NTOT,I,NT

      OPEN(unit = 91,FILE=FILEEXP,STATUS='unknown')

      DO I = 1,NLINES,1
         READ(91,'(A)') LINE    ! skip NLINES..
      ENDDO

      NT = 0
 10   NT =  NT +1

      READ(91,*,ERR = 99, end = 11)
     +     X1(NT),X2(NT),X3(NT),X4(NT),X5(NT),X6(NT),
     +     X7(NT),X8(NT)

      goto 10
 11   if(NT.GT.NTOT) goto 101

      CLOSE(91)

      GOTO 101
99    PRINT *, 'ERROR IN READING FILE    ', FILEEXP
101   CONTINUE
      RETURN
      END

c-----------------------------------------
c     read data file which has 7 columns
c-----------------------------------------
      SUBROUTINE READDATA7(FILEEXP,NLINES,
     >                    X1,X2,X3,X4,X5,X6,X7,NTOT)
      CHARACTER*72 FILEEXP,LINE
      REAL*8 X1(300),X2(300),X3(300),X4(300),X5(300),X6(300)
      REAL*8 X7(300)
      INTEGER NLINES,NTOT,I,NT

      OPEN(unit = 91,FILE=FILEEXP,STATUS='unknown')

      DO I = 1,NLINES,1
         READ(91,'(A)') LINE    ! skip NLINES..
      ENDDO

      NT = 0
 10   NT =  NT +1

      READ(91,*,ERR = 99, end = 11)
     +     X1(NT),X2(NT),X3(NT),X4(NT),X5(NT),X6(NT),
     +     X7(NT)

      goto 10
 11   if(NT.GT.NTOT) goto 101

      CLOSE(91)

      GOTO 101
99    PRINT *, 'ERROR IN READING FILE    ', FILEEXP
101   CONTINUE
      RETURN
      END


c-----------------------------------------
c     read data file which has 6 columns
c-----------------------------------------
      SUBROUTINE READDATA6(FILEEXP,NLINES,
     >                    X1,X2,X3,X4,X5,X6,NTOT)
      CHARACTER*72 FILEEXP,LINE
      REAL*8 X1(300),X2(300),X3(300),X4(300),X5(300),X6(300)
      INTEGER NLINES,NTOT,I,NT

      OPEN(unit = 91,FILE=FILEEXP,STATUS='unknown')

      DO I = 1,NLINES,1
         READ(91,'(A)') LINE    ! skip NLINES..
      ENDDO

      NT = 0
 10   NT =  NT +1

      READ(91,*,ERR = 99, end = 11)
     +     X1(NT),X2(NT),X3(NT),X4(NT),X5(NT),X6(NT)

      goto 10
 11   if(NT.GT.NTOT) goto 101

      CLOSE(91)

      GOTO 101
99    PRINT *, 'ERROR IN READING FILE    ', FILEEXP
101   CONTINUE
      RETURN
      END
c-----------------------------------------
c     read data file which has 5 columns
c-----------------------------------------
      SUBROUTINE READDATA5(FILEEXP,NLINES,
     >                    X1,X2,X3,X4,X5,NTOT)
      CHARACTER*72 FILEEXP,LINE
      REAL*8 X1(300),X2(300),X3(300),X4(300),X5(300)
      INTEGER NLINES,NTOT,I,NT

      OPEN(unit = 91,FILE=FILEEXP,STATUS='unknown')

      DO I = 1,NLINES,1
         READ(91,'(A)') LINE    ! skip NLINES..
      ENDDO

      NT = 0
 10   NT =  NT +1

      READ(91,*,ERR = 99, end = 11)
     +     X1(NT),X2(NT),X3(NT),X4(NT),X5(NT)

      goto 10
 11   if(NT.GT.NTOT) goto 101

      CLOSE(91)

      GOTO 101
99    PRINT *, 'ERROR IN READING FILE    ', FILEEXP
101   CONTINUE
      RETURN
      END

c-----------------------------------------
c     read data file which has 4 columns
c-----------------------------------------
      SUBROUTINE READDATA4(FILEEXP,NLINES,
     >                    X1,X2,X3,X4,NTOT)
      CHARACTER*72 FILEEXP,LINE
      REAL*8 X1(300),X2(300),X3(300),X4(300)
      INTEGER NLINES,NTOT,I,NT

      OPEN(unit = 91,FILE=FILEEXP,STATUS='unknown')

      DO I = 1,NLINES,1
         READ(91,'(A)') LINE    ! skip NLINES..
      ENDDO

      NT = 0
 10   NT =  NT +1

      READ(91,*,ERR = 99, end = 11)
     +     X1(NT),X2(NT),X3(NT),X4(NT)

      goto 10
 11   if(NT.GT.NTOT) goto 101

      CLOSE(91)

      GOTO 101
99    PRINT *, 'ERROR IN READING FILE    ', FILEEXP
101   CONTINUE
      RETURN
      END

c-----------------------------------------
c     read data file which has 3 columns
c-----------------------------------------
      SUBROUTINE READDATA3(FILEEXP,NLINES,
     >                    X1,X2,X3,NTOT)
      CHARACTER*72 FILEEXP,LINE
      REAL*8 X1(300),X2(300),X3(300)
      INTEGER NLINES,NTOT,I,NT

      OPEN(unit = 91,FILE=FILEEXP,STATUS='unknown')

      DO I = 1,NLINES,1
         READ(91,'(A)') LINE    ! skip NLINES..
      ENDDO

      NT = 0
 10   NT =  NT +1

      READ(91,*,ERR = 99, end = 11)
     +     X1(NT),X2(NT),X3(NT)

      goto 10
 11   if(NT.GT.NTOT) goto 101

      CLOSE(91)

      GOTO 101
99    PRINT *, 'ERROR IN READING FILE    ', FILEEXP
101   CONTINUE
      RETURN
      END
