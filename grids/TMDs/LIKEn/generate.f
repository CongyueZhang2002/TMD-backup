      program master
      implicit none
      CHARACTER*72 LIKEn_VC_DIR, LIKEn_VC_GRIDS
      CHARACTER*72 LIKEn_HE_DIR, LIKEn_HE_GRIDS
      CHARACTER*72 LIKEn_NE_DIR, LIKEn_NE_GRIDS
      CHARACTER*72 LIKEn_KR_DIR, LIKEn_KR_GRIDS
      CHARACTER*72 LIKEn_XE_DIR, LIKEn_XE_GRIDS
      CHARACTER*72 LIKEn_CC_DIR, LIKEn_CC_GRIDS
      CHARACTER*72 LIKEn_FE_DIR, LIKEn_FE_GRIDS
      CHARACTER*72 LIKEn_AU_DIR, LIKEn_AU_GRIDS
      CHARACTER*72 LIKEn_PB_DIR, LIKEn_PB_GRIDS
      LOGICAL DIRECTORY_EXISTS
      INTEGER I_VC
      INTEGER I_HE, I_NE, I_KR, I_XE
      INTEGER I_CC, I_FE, I_AU, I_PB
      CHARACTER*72 MKDIR
      integer NB, NZ, NQ
      PARAMETER (NZ=47, NQ=24, NB=51)
      real*8 Q, BB(NB), zh
      real*8 Q2(NQ), QQ(NQ)
      real*8 ZHH(NZ)
      integer I,J
      integer A,Z
      INTEGER PARTICLE_ID(-3:21)
      real*8 b
      real*8 u,d,s,c,bq,ub,db,sb,cb,bqb,g
      integer nloops,IH,IC
      common /meson/ IH,IC
      integer IIREAD, FININUCLEAR
      COMMON / FRAGININUCLEAR / FININUCLEAR
      integer FINI14
      COMMON / FRAGINI14 / FINI14

C-----CHOOSE GRIDS TO EXPORT
      I_VC = 0
      I_HE = 0
      I_NE = 0
      I_KR = 0
      I_XE = 0
      I_CC = 0
      I_FE = 0
      I_AU = 1
      I_PB = 0

      nloops = 2
      IH = 1
      IC = 1

C.....Q**2 VALUES OF THE GRID.
       DATA Q2 / 1.d0, 1.25D0, 1.5D0, 2.5D0,
     1           4.0D0, 6.4D0, 1.0D1, 1.5D1, 2.5D1, 4.0D1, 6.4D1,
     2           1.0D2, 1.8D2, 3.2D2, 5.8D2, 1.0D3, 1.8D3,
     3           3.2D3, 5.8D3, 1.0D4, 1.8D4, 3.2D4, 5.8D4, 1.0D5/
C......Z VALUES OF THE GRID.
       DATA ZHH /0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
     4        0.095, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
     5        0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475,  0.5,
     6        0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7,
     7        0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9,
     8        0.925, 0.95, 0.975, 1.0/

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
      LIKEn_VC_DIR = 'LIKEnVC/'
      LIKEn_HE_DIR = 'LIKEnHE/'
      LIKEn_NE_DIR = 'LIKEnNE/'
      LIKEn_KR_DIR = 'LIKEnKR/'
      LIKEn_XE_DIR = 'LIKEnXE/'
      LIKEn_CC_DIR = 'LIKEnCC/'
      LIKEn_FE_DIR = 'LIKEnFE/'
      LIKEn_AU_DIR = 'LIKEnAU/'
      LIKEn_PB_DIR = 'LIKEnPB/'

*-----GRID FILE NAMES:
      LIKEn_VC_GRIDS = 'LIKEnVC_0000.dat'
      LIKEn_HE_GRIDS = 'LIKEnHE_0000.dat'
      LIKEn_NE_GRIDS = 'LIKEnNE_0000.dat'
      LIKEn_KR_GRIDS = 'LIKEnKR_0000.dat'
      LIKEn_XE_GRIDS = 'LIKEnXE_0000.dat'
      LIKEn_CC_GRIDS = 'LIKEnCC_0000.dat'
      LIKEn_FE_GRIDS = 'LIKEnFE_0000.dat'
      LIKEn_AU_GRIDS = 'LIKEnAU_0000.dat'
      LIKEn_PB_GRIDS = 'LIKEnPB_0000.dat'

C-------Vacuum

      IF (I_VC.eq.1) THEN

      FINI14 = 0

      DO J = 1, NQ
      QQ(J) = sqrt(Q2(J))
      ENDDO

        INQUIRE(FILE=LIKEn_VC_DIR, EXIST=DIRECTORY_EXISTS)
        IF (.not.DIRECTORY_EXISTS) THEN
            CALL SYSTEM(MKDIR//LIKEn_VC_DIR)
        ENDIF
*-------Generate LHAPDF Grids
        OPEN(UNIT = 4, FILE = trim(LIKEn_VC_DIR)//LIKEn_VC_GRIDS)
        WRITE(4,*) 'PdfType: central'
        WRITE(4,*) 'Format: lhagrid1'
        WRITE(4,*) '---'
        WRITE(4,100) ZHH
        WRITE(4,101) QQ
        WRITE(4,102) PARTICLE_ID(-3:-1),PARTICLE_ID(1:3),PARTICLE_ID(21)

        DO I = 1, NZ
           zh = ZHH(I)
          DO J = 1, NQ
            Q = QQ(J)
            print *, "VC",zh,Q
            call FF(zh,Q,u,ub,d,db,s,sb,g,nloops)
            WRITE(4,103) zh*SB,zh*UB,zh*DB,zh*D,zh*U,zh*S,zh*G
          ENDDO
        ENDDO


        WRITE(4,*) '---'
        CLOSE(4)
      ENDIF

C-------HE

      IF (I_HE.eq.1) THEN
      A = 4
      Z = 2

      FININUCLEAR = 0

      DO J = 1, NQ
      QQ(J) = sqrt(Q2(J))
      ENDDO

        INQUIRE(FILE=LIKEn_HE_DIR, EXIST=DIRECTORY_EXISTS)
        IF (.not.DIRECTORY_EXISTS) THEN
            CALL SYSTEM(MKDIR//LIKEn_HE_DIR)
        ENDIF
*-------Generate LHAPDF Grids
        OPEN(UNIT = 4, FILE = trim(LIKEn_HE_DIR)//LIKEn_HE_GRIDS)
        WRITE(4,*) 'PdfType: central'
        WRITE(4,*) 'Format: lhagrid1'
        WRITE(4,*) '---'
        WRITE(4,100) ZHH
        WRITE(4,101) QQ
        WRITE(4,102) PARTICLE_ID(-3:-1),PARTICLE_ID(1:3),PARTICLE_ID(21)

        DO I = 1, NZ
           zh = ZHH(I)
          DO J = 1, NQ
            Q = QQ(J)
            print *, "HE",zh,Q
            call LIKEnFF(zh,Q,A,u,ub,d,db,s,sb,g,nloops)
            WRITE(4,103) zh*SB,zh*UB,zh*DB,zh*D,zh*U,zh*S,zh*G
          ENDDO
        ENDDO


        WRITE(4,*) '---'
        CLOSE(4)
      ENDIF

C-------NE

      IF (I_NE.eq.1) THEN
      A = 20
      Z = 10

      FININUCLEAR = 0

      DO J = 1, NQ
      QQ(J) = sqrt(Q2(J))
      ENDDO

        INQUIRE(FILE=LIKEn_NE_DIR, EXIST=DIRECTORY_EXISTS)
        IF (.not.DIRECTORY_EXISTS) THEN
            CALL SYSTEM(MKDIR//LIKEn_NE_DIR)
        ENDIF
*-------Generate LHAPDF Grids
        OPEN(UNIT = 4, FILE = trim(LIKEn_NE_DIR)//LIKEn_NE_GRIDS)
        WRITE(4,*) 'PdfType: central'
        WRITE(4,*) 'Format: lhagrid1'
        WRITE(4,*) '---'
        WRITE(4,100) ZHH
        WRITE(4,101) QQ
        WRITE(4,102) PARTICLE_ID(-3:-1),PARTICLE_ID(1:3),PARTICLE_ID(21)

        DO I = 1, NZ
          zh = ZHH(I)
          DO J = 1, NQ
            Q = QQ(J)
            print *, "NE",zh,Q
            call LIKEnFF(zh,Q,A,u,ub,d,db,s,sb,g,nloops)
            WRITE(4,103) zh*SB,zh*UB,zh*DB,zh*D,zh*U,zh*S,zh*G
          ENDDO
        ENDDO


        WRITE(4,*) '---'
        CLOSE(4)
      ENDIF

C-------KR

      IF (I_KR.eq.1) THEN
      A = 84
      Z = 36

      FININUCLEAR = 0

      DO J = 1, NQ
      QQ(J) = sqrt(Q2(J))
      ENDDO

        INQUIRE(FILE=LIKEn_KR_DIR, EXIST=DIRECTORY_EXISTS)
        IF (.not.DIRECTORY_EXISTS) THEN
            CALL SYSTEM(MKDIR//LIKEn_KR_DIR)
        ENDIF
*-------Generate LHAPDF Grids
        OPEN(UNIT = 4, FILE = trim(LIKEn_KR_DIR)//LIKEn_KR_GRIDS)
        WRITE(4,*) 'PdfType: central'
        WRITE(4,*) 'Format: lhagrid1'
        WRITE(4,*) '---'
        WRITE(4,100) ZHH
        WRITE(4,101) QQ
        WRITE(4,102) PARTICLE_ID(-3:-1),PARTICLE_ID(1:3),PARTICLE_ID(21)

        DO I = 1, NZ
          zh = ZHH(I)
          DO J = 1, NQ
            Q = QQ(J)
            print *, "KR",zh,Q
            call LIKEnFF(zh,Q,A,u,ub,d,db,s,sb,g,nloops)
            WRITE(4,103) zh*SB,zh*UB,zh*DB,zh*D,zh*U,zh*S,zh*G
          ENDDO
        ENDDO


        WRITE(4,*) '---'
        CLOSE(4)
      ENDIF


C-------XE

      IF (I_XE.eq.1) THEN
      A = 131
      Z = 54

      FININUCLEAR = 0

      DO J = 1, NQ
      QQ(J) = sqrt(Q2(J))
      ENDDO
        INQUIRE(FILE=LIKEn_XE_DIR, EXIST=DIRECTORY_EXISTS)
        IF (.not.DIRECTORY_EXISTS) THEN
            CALL SYSTEM(MKDIR//LIKEn_XE_DIR)
        ENDIF
*-------Generate LHAPDF Grids
        OPEN(UNIT = 4, FILE = trim(LIKEn_XE_DIR)//LIKEn_XE_GRIDS)
        WRITE(4,*) 'PdfType: central'
        WRITE(4,*) 'Format: lhagrid1'
        WRITE(4,*) '---'
        WRITE(4,100) ZHH
        WRITE(4,101) QQ
        WRITE(4,102) PARTICLE_ID(-3:-1),PARTICLE_ID(1:3),PARTICLE_ID(21)

        DO I = 1, NZ
          zh = ZHH(I)
          DO J = 1, NQ
            Q = QQ(J)
            print *, "XE",zh,Q
            call LIKEnFF(zh,Q,A,u,ub,d,db,s,sb,g,nloops)
            WRITE(4,103) zh*SB,zh*UB,zh*DB,zh*D,zh*U,zh*S,zh*G
          ENDDO
        ENDDO


        WRITE(4,*) '---'
        CLOSE(4)
      ENDIF

C-------CC

      IF (I_CC.eq.1) THEN
      A = 12
      Z = 6

      FININUCLEAR = 0

      DO J = 1, NQ
      QQ(J) = sqrt(Q2(J))
      ENDDO

        INQUIRE(FILE=LIKEn_CC_DIR, EXIST=DIRECTORY_EXISTS)
        IF (.not.DIRECTORY_EXISTS) THEN
            CALL SYSTEM(MKDIR//LIKEn_CC_DIR)
        ENDIF
*-------Generate LHAPDF Grids
        OPEN(UNIT = 4, FILE = trim(LIKEn_CC_DIR)//LIKEn_CC_GRIDS)
        WRITE(4,*) 'PdfType: central'
        WRITE(4,*) 'Format: lhagrid1'
        WRITE(4,*) '---'
        WRITE(4,100) ZHH
        WRITE(4,101) QQ
        WRITE(4,102) PARTICLE_ID(-3:-1),PARTICLE_ID(1:3),PARTICLE_ID(21)

        DO I = 1, NZ
          zh = ZHH(I)
          DO J = 1, NQ
            Q = QQ(J)
            print *, "CC",zh,Q
            call LIKEnFF(zh,Q,A,u,ub,d,db,s,sb,g,nloops)
            WRITE(4,103) zh*SB,zh*UB,zh*DB,zh*D,zh*U,zh*S,zh*G
          ENDDO
        ENDDO


        WRITE(4,*) '---'
        CLOSE(4)
      ENDIF

C-------FE

      IF (I_FE.eq.1) THEN
      A = 56
      Z = 12

      FININUCLEAR = 0

      DO J = 1, NQ
      QQ(J) = sqrt(Q2(J))
      ENDDO

        INQUIRE(FILE=LIKEn_FE_DIR, EXIST=DIRECTORY_EXISTS)
        IF (.not.DIRECTORY_EXISTS) THEN
            CALL SYSTEM(MKDIR//LIKEn_FE_DIR)
        ENDIF
*-------Generate LHAPDF Grids
        OPEN(UNIT = 4, FILE = trim(LIKEn_FE_DIR)//LIKEn_FE_GRIDS)
        WRITE(4,*) 'PdfType: central'
        WRITE(4,*) 'Format: lhagrid1'
        WRITE(4,*) '---'
        WRITE(4,100) ZHH
        WRITE(4,101) QQ
        WRITE(4,102) PARTICLE_ID(-3:-1),PARTICLE_ID(1:3),PARTICLE_ID(21)

        DO I = 1, NZ
          zh = ZHH(I)
          DO J = 1, NQ
            Q = QQ(J)
            print *, "FE",zh,Q
            call LIKEnFF(zh,Q,A,u,ub,d,db,s,sb,g,nloops)
            WRITE(4,103) zh*SB,zh*UB,zh*DB,zh*D,zh*U,zh*S,zh*G
          ENDDO
        ENDDO


        WRITE(4,*) '---'
        CLOSE(4)
      ENDIF

C-------AU

      IF (I_AU.eq.1) THEN
      A = 197
      Z = 79

      FININUCLEAR = 0

      DO J = 1, NQ
      QQ(J) = sqrt(Q2(J))
      ENDDO

        INQUIRE(FILE=LIKEn_AU_DIR, EXIST=DIRECTORY_EXISTS)
        IF (.not.DIRECTORY_EXISTS) THEN
            CALL SYSTEM(MKDIR//LIKEn_AU_DIR)
        ENDIF
*-------Generate LHAPDF Grids
        OPEN(UNIT = 4, FILE = trim(LIKEn_AU_DIR)//LIKEn_AU_GRIDS)
        WRITE(4,*) 'PdfType: central'
        WRITE(4,*) 'Format: lhagrid1'
        WRITE(4,*) '---'
        WRITE(4,100) ZHH
        WRITE(4,101) QQ
        WRITE(4,102) PARTICLE_ID(-3:-1),PARTICLE_ID(1:3),PARTICLE_ID(21)

        DO I = 1, NZ
          zh = ZHH(I)
          DO J = 1, NQ
            Q = QQ(J)
            print *, "AU",zh,Q
            call LIKEnFF(zh,Q,A,u,ub,d,db,s,sb,g,nloops)
            WRITE(4,103) zh*SB,zh*UB,zh*DB,zh*D,zh*U,zh*S,zh*G
          ENDDO
        ENDDO


        WRITE(4,*) '---'
        CLOSE(4)
      ENDIF

C-------PB

      IF (I_PB.eq.1) THEN
      A = 208
      Z = 82

      FININUCLEAR = 0

      DO J = 1, NQ
      QQ(J) = sqrt(Q2(J))
      ENDDO

        INQUIRE(FILE=LIKEn_PB_DIR, EXIST=DIRECTORY_EXISTS)
        IF (.not.DIRECTORY_EXISTS) THEN
            CALL SYSTEM(MKDIR//LIKEn_PB_DIR)
        ENDIF
*-------Generate LHAPDF Grids
        OPEN(UNIT = 4, FILE = trim(LIKEn_PB_DIR)//LIKEn_PB_GRIDS)
        WRITE(4,*) 'PdfType: central'
        WRITE(4,*) 'Format: lhagrid1'
        WRITE(4,*) '---'
        WRITE(4,100) ZHH
        WRITE(4,101) QQ
        WRITE(4,102) PARTICLE_ID(-3:-1),PARTICLE_ID(1:3),PARTICLE_ID(21)

        DO I = 1, NZ
          zh = ZHH(I)
          DO J = 1, NQ
            Q = QQ(J)
            print *, "PB",zh,Q
            call LIKEnFF(zh,Q,A,u,ub,d,db,s,sb,g,nloops)
            WRITE(4,103) zh*SB,zh*UB,zh*DB,zh*D,zh*U,zh*S,zh*G
          ENDDO
        ENDDO


        WRITE(4,*) '---'
        CLOSE(4)
      ENDIF

*------FORMATTING
  100 FORMAT(47(E12.6E2,' '))
  101 FORMAT(51(E12.6E2,' '))
  102 FORMAT(' ',3(I2,' '),3(I1,' '),I2)
  103 FORMAT(7(E12.6E2,'  '))

      return
      end
