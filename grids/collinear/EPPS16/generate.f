      PROGRAM MASTER
      IMPLICIT NONE
      CHARACTER*72 EPPS16_HE_DIR, EPPS16_HE_GRIDS
      CHARACTER*72 EPPS16_NE_DIR, EPPS16_NE_GRIDS
      CHARACTER*72 EPPS16_KR_DIR, EPPS16_KR_GRIDS
      CHARACTER*72 EPPS16_XE_DIR, EPPS16_XE_GRIDS
      LOGICAL DIRECTORY_EXISTS
      INTEGER I_HE, I_NE, I_KR, I_XE
      INTEGER NUM_XVALUES
      INTEGER NUM_QVALUES
      INTEGER I, J
      CHARACTER*72 MKDIR
      REAL*8 X(81)
      REAL*8 Q(31)
      INTEGER PARTICLE_ID(-5:21)
      INTEGER A,Z
      REAL*8 U,D,UB,DB,SS,SB,GL
      INTEGER NLOOPS
      character*100 outdir

C----- Set nloops
      nloops = 1


C----- This program generates LHAPDF Grids from the EPPS16 Driver

C-----CHOOSE GRIDS TO EXPORT
      I_HE = 1
      I_NE = 1
      I_KR = 1
      I_XE = 1

C-----SET UP LHAPDF PROTON PDF GRIDS (CT14NLO)
      CALL InitPDFsetByNameM(1,"CT14nlo")
      call InitPDFM(1,0)

*-----Set up EPPS16's ordered list of x-values
      NUM_XVALUES = 81
      OPEN(UNIT = 1, FILE = 'GRID_XVALUES.csv',STATUS = 'unknown')
      DO I = 1, NUM_XVALUES, 1
        READ(1,*) X(I)
      ENDDO
      CLOSE(1)

*-----Set up EPPS16's ordered list of Q-values
      NUM_QVALUES = 31
      OPEN(UNIT = 2, FILE = 'GRID_QVALUES.csv',STATUS = 'unknown')
      DO I = 1, NUM_QVALUES, 1
        READ(2,*) Q(I)
      ENDDO
      CLOSE(2)

C-----PARTICLE ID NUMBERS: Based on Standard PDG ID Numbers:
C------https://pdg.lbl.gov/2006/reviews/pdf-files/montecarlo-web.pdf
      PARTICLE_ID( 1) =  1  ! d
      PARTICLE_ID(-1) = -1  ! dbar
      PARTICLE_ID( 2) =  2  ! u
      PARTICLE_ID(-2) = -2  ! ubar
      PARTICLE_ID( 3) =  3  ! s
      PARTICLE_ID(-3) = -3  ! sbar
      PARTICLE_ID(21) = 21  ! g

*-----Set Grid Directory Names
      EPPS16_HE_DIR = 'EPPS16nlo_CT14nlo_He4/'
      EPPS16_NE_DIR = 'EPPS16nlo_CT14nlo_Ne20/'
      EPPS16_KR_DIR = 'EPPS16nlo_CT14nlo_Kr84/'
      EPPS16_XE_DIR = 'EPPS16nlo_CT14nlo_Xe131/'

*-----GRID FILE NAMES:
      EPPS16_HE_GRIDS = 'EPPS16nlo_CT14nlo_He4_0000.dat'
      EPPS16_NE_GRIDS = 'EPPS16nlo_CT14nlo_Ne20_0000.dat'
      EPPS16_KR_GRIDS = 'EPPS16nlo_CT14nlo_Kr84_0000.dat'
      EPPS16_XE_GRIDS = 'EPPS16nlo_CT14nlo_Xe131_0000.dat'


*-----Create Directories if they do not exist
      MKDIR = 'mkdir'

C------HE
      IF (I_HE.eq.1) THEN
        INQUIRE(FILE=EPPS16_HE_DIR, EXIST=DIRECTORY_EXISTS)
        IF (.not.DIRECTORY_EXISTS) THEN
            CALL SYSTEM(MKDIR//EPPS16_HE_DIR)
        ENDIF
*-------Generate LHAPDF Grids
        outdir = trim(EPPS16_HE_DIR)//trim(EPPS16_HE_GRIDS)
        OPEN(UNIT = 4, FILE = outdir)
        WRITE(4,*) 'PdfType: central'
        WRITE(4,*) 'Format: lhagrid1'
        WRITE(4,*) '---'
        WRITE(4,100) X
        WRITE(4,101) Q
        WRITE(4,102) PARTICLE_ID(-3:-1),PARTICLE_ID(1:3),PARTICLE_ID(21)

C-------HE (A = 4, Z = 2)
        A = 4
        Z = 2

        DO I = 1, NUM_XVALUES
          DO J = 1, NUM_QVALUES
            print *, X(I),Q(J)
            CALL PDF_A_EPPS16(X(I),Q(J),A,Z,U,D,UB,DB,SS,SB,GL,nloops)
            WRITE(4,103) X(I)*SB, X(I)*UB, X(I)*DB, X(I)*D, X(I)*U,
     +                   X(I)*SS, X(I)*GL
          ENDDO
        ENDDO


        WRITE(4,*) '---'
        CLOSE(4)
      ENDIF

C------NE
      IF (I_NE.eq.1) THEN
        INQUIRE(FILE=EPPS16_NE_DIR, EXIST=DIRECTORY_EXISTS)
        IF (.not.DIRECTORY_EXISTS) THEN
            CALL SYSTEM(MKDIR//EPPS16_NE_DIR)
        ENDIF
*-------Generate LHAPDF Grids

        outdir = trim(EPPS16_NE_DIR)//trim(EPPS16_NE_GRIDS)
        OPEN(UNIT = 20, FILE = outdir)
        WRITE(20,*) 'PdfType: central'
        WRITE(20,*) 'Format: lhagrid1'
        WRITE(20,*) '---'
        WRITE(20,100) X
        WRITE(20,101) Q
        WRITE(20,102) PARTICLE_ID(-3:-1),PARTICLE_ID(1:3),
     +                PARTICLE_ID(21)

C-------NE (A = 20, Z = 10)
        A = 20
        Z = 10

        DO I = 1, NUM_XVALUES
          DO J = 1, NUM_QVALUES
            CALL PDF_A_EPPS16(X(I),Q(J),A,Z,U,D,UB,DB,SS,SB,GL,nloops)
            WRITE(20,103) X(I)*SB, X(I)*UB, X(I)*DB, X(I)*D, X(I)*U,
     +                   X(I)*SS, X(I)*GL
          ENDDO
        ENDDO

        WRITE(20,*) '---'
        CLOSE(20)
      ENDIF


C------KR
      IF (I_KR.eq.1) THEN
        INQUIRE(FILE=EPPS16_KR_DIR, EXIST=DIRECTORY_EXISTS)
        IF (.not.DIRECTORY_EXISTS) THEN
            CALL SYSTEM(MKDIR//EPPS16_KR_DIR)
        ENDIF
*-------Generate LHAPDF Grids
        outdir = trim(EPPS16_KR_DIR)//trim(EPPS16_KR_GRIDS)
        OPEN(UNIT = 84, FILE = outdir)
        WRITE(84,*) 'PdfType: central'
        WRITE(84,*) 'Format: lhagrid1'
        WRITE(84,*) '---'
        WRITE(84,100) X
        WRITE(84,101) Q
        WRITE(84,102) PARTICLE_ID(-3:-1),PARTICLE_ID(1:3),
     +                PARTICLE_ID(21)

C-------KR (A = 84, Z = 36)
        A = 84
        Z = 36

        DO I = 1, NUM_XVALUES
          DO J = 1, NUM_QVALUES
            CALL PDF_A_EPPS16(X(I),Q(J),A,Z,U,D,UB,DB,SS,SB,GL,nloops)
            WRITE(84,103) X(I)*SB, X(I)*UB, X(I)*DB, X(I)*D, X(I)*U,
     +                   X(I)*SS, X(I)*GL
          ENDDO
        ENDDO

        WRITE(84,*) '---'
        CLOSE(84)
      ENDIF


C------XE
      IF (I_XE.eq.1) THEN
        INQUIRE(FILE=EPPS16_XE_DIR, EXIST=DIRECTORY_EXISTS)
        IF (.not.DIRECTORY_EXISTS) THEN
            CALL SYSTEM(MKDIR//EPPS16_XE_DIR)
        ENDIF
*-------Generate LHAPDF Grids
        outdir = trim(EPPS16_XE_DIR)//trim(EPPS16_XE_GRIDS)
        OPEN(UNIT = 131, FILE = outdir)
        WRITE(131,*) 'PdfType: central'
        WRITE(131,*) 'Format: lhagrid1'
        WRITE(131,*) '---'
        WRITE(131,100) X
        WRITE(131,101) Q
        WRITE(131,102) PARTICLE_ID(-3:-1),PARTICLE_ID(1:3),
     +                 PARTICLE_ID(21)

C-------XE (A = 131, Z = 54)
        A = 131
        Z = 54

        DO I = 1, NUM_XVALUES
          DO J = 1, NUM_QVALUES
            CALL PDF_A_EPPS16(X(I),Q(J),A,Z,U,D,UB,DB,SS,SB,GL,nloops)
            WRITE(131,103) X(I)*SB, X(I)*UB, X(I)*DB, X(I)*D, X(I)*U,
     +                   X(I)*SS, X(I)*GL
          ENDDO
        ENDDO

        WRITE(131,*) '---'
        CLOSE(131)
      ENDIF

*------FORMATTING
  100 FORMAT(81(E12.6E2,' '))
  101 FORMAT(31(E12.6E2,' '))
  102 FORMAT(' ',3(I2,' '),3(I1,' '),I2)
  103 FORMAT(7(E12.6E2,'  '))
  202 FORMAT(' ',5(I2,' '),5(I1,' '),I2)



      END PROGRAM MASTER
