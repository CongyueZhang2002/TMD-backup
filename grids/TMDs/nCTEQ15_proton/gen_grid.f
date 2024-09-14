      program master
      implicit none
      CHARACTER*72 nCTEQ15_HE_DIR   ,nCTEQ15_HE_GRIDS
      CHARACTER*72 nCTEQ15_NE_DIR   ,nCTEQ15_NE_GRIDS
      CHARACTER*72 nCTEQ15_CA_DIR   ,nCTEQ15_CA_GRIDS
      CHARACTER*72 nCTEQ15_KR_DIR   ,nCTEQ15_KR_GRIDS
      CHARACTER*72 nCTEQ15_XE_DIR   ,nCTEQ15_XE_GRIDS
      CHARACTER*72 nCTEQ15_BE_DIR   ,nCTEQ15_BE_GRIDS
      CHARACTER*72 nCTEQ15_FE_DIR   ,nCTEQ15_FE_GRIDS
      CHARACTER*72 nCTEQ15_WW_DIR   ,nCTEQ15_WW_GRIDS
      CHARACTER*72 nCTEQ15_JCC_DIR  ,nCTEQ15_JCC_GRIDS
      CHARACTER*72 nCTEQ15_JFE_DIR  ,nCTEQ15_JFE_GRIDS
      CHARACTER*72 nCTEQ15_JPB_DIR  ,nCTEQ15_JPB_GRIDS
      CHARACTER*72 nCTEQ15_LPB_DIR  ,nCTEQ15_LPB_GRIDS
      CHARACTER*72 nCTEQ15_AU_DIR   ,nCTEQ15_AU_GRIDS
      CHARACTER*72 nCTEQ15_PR_DIR   ,nCTEQ15_PR_GRIDS
      LOGICAL DIRECTORY_EXISTS
      INTEGER I_HE, I_NE, I_KR, I_XE
      INTEGER I_BE, I_FE, I_WW
      INTEGER I_JCC,I_JFE,I_JPB
      INTEGER I_PB, I_AU ,I_PR, I_CA
      CHARACTER*72 MKDIR
      integer NXB, NQQ
      parameter (NXB = 240,NQQ = 37)
      real*8 x, Q, XB(NXB), QQ(NQQ)
      integer I,J
      integer A,Z
      INTEGER PARTICLE_ID(-3:21)
      real*8 u,d,s,c,bq,ub,db,sb,cb,bqb,g
      integer nloops, TIC, TIH
      COMMON / ttarget / TIH,TIC
      integer IIREAD

C-----CHOOSE GRIDS TO EXPORT
      I_PR = 0
      I_HE = 1
      I_NE = 1
      I_KR = 1
      I_XE = 1
      I_BE = 1
      I_FE = 1
      I_WW = 1
      I_JCC= 1
      I_JFE= 1
      I_JPB= 1
      I_PB = 1
      I_AU = 1
      I_CA = 1

      nloops = 2

      call SetLHAPARM('SILENT') ! To not show the calls, although they are called
C-----INITIALIZE PDF SETS USING LHAPDF
      CALL InitPDFsetByNameM(1,"CT14nlo")
      CALL InitPDFM(1,0)

      DATA XB /
     * 1.00000000e-07, 1.06976578e-07, 1.14439883e-07, 1.22423871e-07,
     * 1.30964869e-07, 1.40101735e-07, 1.49876043e-07, 1.60332262e-07,
     * 1.71517968e-07, 1.83484054e-07, 1.96284962e-07, 2.09978937e-07,
     * 2.24628282e-07, 2.40299650e-07, 2.57064343e-07, 2.74998639e-07,
     * 2.94184134e-07, 3.14708121e-07, 3.36663980e-07, 3.60151606e-07,
     * 3.85277865e-07, 4.12157077e-07, 4.40911539e-07, 4.71672078e-07,
     * 5.04578650e-07, 5.39780975e-07, 5.77439217e-07, 6.17724717e-07,
     * 6.60820766e-07, 7.06923445e-07, 7.56242513e-07, 8.09002364e-07,
     * 8.65443048e-07, 9.25821361e-07, 9.90412013e-07, 1.05950888e-06,
     * 1.13342635e-06, 1.21250073e-06, 1.29709179e-06, 1.38758442e-06,
     * 1.48439033e-06, 1.58794999e-06, 1.69873456e-06, 1.81724811e-06,
     * 1.94402985e-06, 2.07965661e-06, 2.22474549e-06, 2.37995660e-06,
     * 2.54599614e-06, 2.72361955e-06, 2.91363501e-06, 3.11690704e-06,
     * 3.33436050e-06, 3.56698477e-06, 3.81583826e-06, 4.08205321e-06,
     * 4.36684085e-06, 4.67149692e-06, 4.99740756e-06, 5.34605562e-06,
     * 5.71902738e-06, 6.11801981e-06, 6.54484825e-06, 7.00145472e-06,
     * 7.48991670e-06, 8.01245660e-06, 8.57145192e-06, 9.16944598e-06,
     * 9.80915956e-06, 1.04935033e-05, 1.12255907e-05, 1.20087529e-05,
     * 1.28465529e-05, 1.37428028e-05, 1.47015802e-05, 1.57272474e-05,
     * 1.68244712e-05, 1.79982436e-05, 1.92539052e-05, 2.05971690e-05,
     * 2.20341466e-05, 2.35713761e-05, 2.52158516e-05, 2.69750553e-05,
     * 2.88569911e-05, 3.08702217e-05, 3.30239069e-05, 3.53278457e-05,
     * 3.77925205e-05, 4.04291453e-05, 4.32497164e-05, 4.62670667e-05,
     * 4.94949249e-05, 5.29479771e-05, 5.66419342e-05, 6.05936032e-05,
     * 6.48209634e-05, 6.93432487e-05, 7.41810348e-05, 7.93563328e-05,
     * 8.48926895e-05, 9.08152945e-05, 9.71510947e-05, 1.03928917e-04,
     * 1.11179599e-04, 1.18936131e-04, 1.27233804e-04, 1.36110370e-04,
     * 1.45606216e-04, 1.55764548e-04, 1.66631584e-04, 1.78256767e-04,
     * 1.90692990e-04, 2.03996836e-04, 2.18228835e-04, 2.33453741e-04,
     * 2.49740824e-04, 2.67164188e-04, 2.85803107e-04, 3.05742385e-04,
     * 3.27072742e-04, 3.49891228e-04, 3.74301664e-04, 4.00415112e-04,
     * 4.28350387e-04, 4.58234587e-04, 4.90203682e-04, 5.24403126e-04,
     * 5.60988521e-04, 6.00126325e-04, 6.41994608e-04, 6.86783865e-04,
     * 7.34697880e-04, 7.85954653e-04, 8.40787396e-04, 8.99445587e-04,
     * 9.62196113e-04, 1.02932448e-03, 1.10113611e-03, 1.17795773e-03,
     * 1.26013888e-03, 1.34805345e-03, 1.44210146e-03, 1.54271080e-03,
     * 1.65033922e-03, 1.76547643e-03, 1.88864628e-03, 2.02040917e-03,
     * 2.16136460e-03, 2.31215389e-03, 2.47346312e-03, 2.64602621e-03,
     * 2.83062831e-03, 3.02810931e-03, 3.23936773e-03, 3.46536475e-03,
     * 3.70712864e-03, 3.96575938e-03, 4.24243369e-03, 4.53841040e-03,
     * 4.85503616e-03, 5.19375156e-03, 5.55609771e-03, 5.94372322e-03,
     * 6.35839173e-03, 6.80198991e-03, 7.27653606e-03, 7.78418931e-03,
     * 8.32725937e-03, 8.90821715e-03, 9.52970590e-03, 1.01945533e-02,
     * 1.09057843e-02, 1.16666349e-02, 1.24805668e-02, 1.33512833e-02,
     * 1.42827461e-02, 1.52791931e-02, 1.63451579e-02, 1.74854907e-02,
     * 1.87053797e-02, 2.00103751e-02, 2.14064146e-02, 2.28998499e-02,
     * 2.44974759e-02, 2.62065615e-02, 2.80348828e-02, 2.99907584e-02,
     * 3.20830871e-02, 3.43213888e-02, 3.67158474e-02, 3.92773573e-02,
     * 4.20175729e-02, 4.49489618e-02, 4.80848614e-02, 5.14395394e-02,
     * 5.50282592e-02, 5.88673488e-02, 6.29742755e-02, 6.73677252e-02,
     * 7.20676873e-02, 7.70955460e-02, 8.24741772e-02, 8.82280528e-02,
     * 9.43833521e-02, 1.00968081e-01, 1.08012198e-01, 1.15547754e-01,
     * 1.23609033e-01, 1.32232714e-01, 1.41458033e-01, 1.51326964e-01,
     * 1.61884408e-01, 1.73178400e-01, 1.85260327e-01, 1.98185159e-01,
     * 2.12011702e-01, 2.26802865e-01, 2.42625944e-01, 2.59552933e-01,
     * 2.77660847e-01, 2.97032074e-01, 3.17754749e-01, 3.39923158e-01,
     * 3.63638164e-01, 3.89007665e-01, 4.16147090e-01, 4.45179918e-01,
     * 4.76238244e-01, 5.09463378e-01, 5.45006490e-01, 5.83029294e-01,
     * 6.23704790e-01, 6.67218044e-01, 7.13767033e-01, 7.63563550e-01,
     * 8.16834159e-01, 8.73821234e-01, 9.34784058e-01, 1.00000000e+00/

      DATA QQ /
     * 1.30000000e+00, 1.66682256e+00, 2.13715189e+00, 2.74019460e+00,
     * 3.51339860e+00, 4.50477851e+00, 5.77589727e+00, 7.40568914e+00,
     * 9.49536134e+00, 1.21746789e+01, 1.56100226e+01, 2.00147214e+01,
     * 2.56622994e+01, 3.29034613e+01, 4.21878705e+01, 5.40920726e+01,
     * 6.93552977e+01, 8.89253654e+01, 1.14017543e+02, 1.46190009e+02,
     * 1.87440620e+02, 2.40330965e+02, 3.08145442e+02, 3.95095211e+02,
     * 5.06579702e+02, 6.49521905e+02, 8.32798281e+02, 1.06778997e+03,
     * 1.36908955e+03, 1.75540719e+03, 2.25073255e+03, 2.88582446e+03,
     * 3.70012101e+03, 4.74418860e+03, 6.08286200e+03, 7.79927048e+03,
     * 1.00000000e+04
     * /

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
      nCTEQ15_HE_DIR = 'nCTEQ15HE_proton/'
      nCTEQ15_CA_DIR = 'nCTEQ15CA_proton/'
      nCTEQ15_NE_DIR = 'nCTEQ15NE_proton/'
      nCTEQ15_KR_DIR = 'nCTEQ15KR_proton/'
      nCTEQ15_XE_DIR = 'nCTEQ15XE_proton/'
      nCTEQ15_BE_DIR = 'nCTEQ15BE_proton/'
      nCTEQ15_FE_DIR = 'nCTEQ15FE_proton/'
      nCTEQ15_WW_DIR = 'nCTEQ15WW_proton/'
      nCTEQ15_JCC_DIR= 'nCTEQ15JCC_proton/'
      nCTEQ15_JFE_DIR= 'nCTEQ15JFE_proton/'
      nCTEQ15_JPB_DIR= 'nCTEQ15JPB_proton/'
      nCTEQ15_LPB_DIR= 'nCTEQ15LPB_proton/'
      nCTEQ15_AU_DIR = 'nCTEQ15AU_proton/'

*-----GRID FILE NAMES:
      nCTEQ15_HE_GRIDS   = 'nCTEQ15HE_proton_0000.dat'
      nCTEQ15_CA_GRIDS   = 'nCTEQ15CA_proton_0000.dat'
      nCTEQ15_NE_GRIDS   = 'nCTEQ15NE_proton_0000.dat'
      nCTEQ15_KR_GRIDS   = 'nCTEQ15KR_proton_0000.dat'
      nCTEQ15_XE_GRIDS   = 'nCTEQ15XE_proton_0000.dat'
      nCTEQ15_BE_GRIDS   = 'nCTEQ15BE_proton_0000.dat'
      nCTEQ15_FE_GRIDS   = 'nCTEQ15FE_proton_0000.dat'
      nCTEQ15_WW_GRIDS   = 'nCTEQ15WW_proton_0000.dat'
      nCTEQ15_JCC_GRIDS  = 'nCTEQ15JCC_proton_0000.dat'
      nCTEQ15_JFE_GRIDS  = 'nCTEQ15JFE_proton_0000.dat'
      nCTEQ15_JPB_GRIDS  = 'nCTEQ15JPB_proton_0000.dat'
      nCTEQ15_LPB_GRIDS  = 'nCTEQ15LPB_proton_0000.dat'
      nCTEQ15_AU_GRIDS   = 'nCTEQ15AU_proton_0000.dat'

C-----INITIALIZE PDF SETS USING LHAPDF
            print *, "1"
            CALL InitPDFsetByNameM(1,"CT14nlo")
            CALL InitPDFM(1,0)
            print *, "2"
            CALL InitPDFsetByNameM(2,"nCTEQ15_4")  !He
            CALL InitPDFM(2,0)
            print *, "3"
            CALL InitPDFsetByNameM(3,"nCTEQ15_20") !Ne
            CALL InitPDFM(3,0)
            print *, "4"
            CALL InitPDFsetByNameM(4,"nCTEQ15_84") !Kr
            CALL InitPDFM(4,0)
            print *, "5"
            CALL InitPDFsetByNameM(5,"nCTEQ15_131") !Xe
            CALL InitPDFM(5,0)
            print *, "6"
            CALL InitPDFsetByNameM(6,"LIKEnHE")
            CALL InitPDFM(6,0)
            print *, "7"
            CALL InitPDFsetByNameM(7,"LIKEnNE")
            CALL InitPDFM(7,0)
            print *, "8"
            CALL InitPDFsetByNameM(8,"LIKEnKR")
            CALL InitPDFM(8,0)
            print *, "9"
            CALL InitPDFsetByNameM(9,"LIKEnXE")
            CALL InitPDFM(9,0)
            print *, "10"
            CALL InitPDFsetByNameM(10,"nCTEQ15_9") ! Be
            CALL InitPDFM(10,0)
            print *, "11"
            CALL InitPDFsetByNameM(11,"nCTEQ15_56") !Fe
            CALL InitPDFM(11,0)
            print *, "12"
            CALL InitPDFsetByNameM(12,"nCTEQ15_184") ! W
            CALL InitPDFM(12,0)
            print *, "13"
            CALL InitPDFsetByNameM(13,"nCTEQ15_12") ! C
            CALL InitPDFM(13,0)
            print *, "14"
            CALL InitPDFsetByNameM(14,"nCTEQ15_56") ! Fe
            CALL InitPDFM(14,0)
            print *, "15"
            CALL InitPDFsetByNameM(15,"nCTEQ15_208") ! Pb
            CALL InitPDFM(15,0)
            print *, "16"
      !      CALL InitPDFsetByNameM(16,"LIKEnCC")
      !      CALL InitPDFM(16,0)
            print *, "17"
      !      CALL InitPDFsetByNameM(17,"LIKEnFE")
      !      CALL InitPDFM(17,0)
            print *, "18"
      !      CALL InitPDFsetByNameM(18,"LIKEnPB")
      !      CALL InitPDFM(18,0)
            print *, "19"
            CALL InitPDFsetByNameM(19,"nCTEQ15_208") ! Pb
            CALL InitPDFM(19,0)
            print *, "20"
            CALL InitPDFsetByNameM(20,"nCTEQ15_197") ! Au
            CALL InitPDFM(20,0)
            print *, "21"
            CALL InitPDFsetByNameM(21,"CT14nlo") ! PR
            CALL InitPDFM(21,0)
            print *, "22"
            CALL InitPDFsetByNameM(22,"nCTEQ15_40") ! Ca
            CALL InitPDFM(22,0)
            print *, "23"
            CALL InitPDFsetByNameM(23,"LIKEnVC")
            CALL InitPDFM(23,0)
            CALL InitPDFsetByNameM(24,"LIKEnAU")
            CALL InitPDFM(24,0)
C-------CA

      IF (I_CA.eq.1) THEN
      A = 40
      Z = 20

        INQUIRE(FILE=nCTEQ15_CA_DIR, EXIST=DIRECTORY_EXISTS)
        IF (.not.DIRECTORY_EXISTS) THEN
            CALL SYSTEM(MKDIR//nCTEQ15_CA_DIR)
        ENDIF
*-------Generate LHAPDF Grids
        OPEN(UNIT = 4, FILE = trim(nCTEQ15_CA_DIR)//nCTEQ15_CA_GRIDS)
        WRITE(4,*) 'PdfType: central'
        WRITE(4,*) 'Format: lhagrid1'
        WRITE(4,*) '---'
        WRITE(4,100) XB
        WRITE(4,101) QQ
        WRITE(4,102) PARTICLE_ID(-3:-1),PARTICLE_ID(1:3),PARTICLE_ID(21)

        DO I = 1, NXB
          x = XB(I)
          DO J = 1, NQQ
            Q = QQ(J)
            print *, "CA",x, Q
            call PDF_A_nCTEQ15(x,Q,A,Z,U,D,UB,DB,S,SB,G,nloops)
            WRITE(4,103) x*SB,x*UB,x*DB,x*D,x*U,x*S,x*G
          ENDDO
        ENDDO

        WRITE(4,*) '---'
	       CLOSE(4)
      ENDIF


C-------HE

      IF (I_HE.eq.1) THEN
      A = 4
      Z = 2

        INQUIRE(FILE=nCTEQ15_HE_DIR, EXIST=DIRECTORY_EXISTS)
        IF (.not.DIRECTORY_EXISTS) THEN
            CALL SYSTEM(MKDIR//nCTEQ15_HE_DIR)
        ENDIF
*-------Generate LHAPDF Grids
        OPEN(UNIT = 4, FILE = trim(nCTEQ15_HE_DIR)//nCTEQ15_HE_GRIDS)
        WRITE(4,*) 'PdfType: central'
        WRITE(4,*) 'Format: lhagrid1'
        WRITE(4,*) '---'
        WRITE(4,100) XB
        WRITE(4,101) QQ
        WRITE(4,102) PARTICLE_ID(-3:-1),PARTICLE_ID(1:3),PARTICLE_ID(21)

        DO I = 1, NXB
          x = XB(I)
          DO J = 1, NQQ
            Q = QQ(J)
            print *, "HE",x, Q
            call PDF_A_nCTEQ15(x,Q,A,Z,U,D,UB,DB,S,SB,G,nloops)
            WRITE(4,103) x*SB,x*UB,x*DB,x*D,x*U,x*S,x*G
          ENDDO
        ENDDO


        WRITE(4,*) '---'
        CLOSE(4)
      ENDIF

C-------NE

      IF (I_NE.eq.1) THEN
      A = 20
      Z = 10

        INQUIRE(FILE=nCTEQ15_NE_DIR, EXIST=DIRECTORY_EXISTS)
        IF (.not.DIRECTORY_EXISTS) THEN
            CALL SYSTEM(MKDIR//nCTEQ15_NE_DIR)
        ENDIF
*-------Generate LHAPDF Grids
        OPEN(UNIT = 4, FILE = trim(nCTEQ15_NE_DIR)//nCTEQ15_NE_GRIDS)
        WRITE(4,*) 'PdfType: central'
        WRITE(4,*) 'Format: lhagrid1'
        WRITE(4,*) '---'
        WRITE(4,100) XB
        WRITE(4,101) QQ
        WRITE(4,102) PARTICLE_ID(-3:-1),PARTICLE_ID(1:3),PARTICLE_ID(21)

        DO I = 1, NXB
          x = XB(I)
          DO J = 1, NQQ
            Q = QQ(J)
       	    print *, "NE",x, Q
            call PDF_A_nCTEQ15(x,Q,A,Z,U,D,UB,DB,S,SB,G,nloops)
            WRITE(4,103) x*SB,x*UB,x*DB,x*D,x*U,x*S,x*G
          ENDDO
        ENDDO


        WRITE(4,*) '---'
        CLOSE(4)
      ENDIF

C-------KR

      IF (I_KR.eq.1) THEN
      A = 84
      Z = 36

        INQUIRE(FILE=nCTEQ15_KR_DIR, EXIST=DIRECTORY_EXISTS)
        IF (.not.DIRECTORY_EXISTS) THEN
            CALL SYSTEM(MKDIR//nCTEQ15_KR_DIR)
        ENDIF
*-------Generate LHAPDF Grids
        OPEN(UNIT = 4, FILE = trim(nCTEQ15_KR_DIR)//nCTEQ15_KR_GRIDS)
        WRITE(4,*) 'PdfType: central'
        WRITE(4,*) 'Format: lhagrid1'
        WRITE(4,*) '---'
        WRITE(4,100) XB
        WRITE(4,101) QQ
        WRITE(4,102) PARTICLE_ID(-3:-1),PARTICLE_ID(1:3),PARTICLE_ID(21)

        DO I = 1, NXB
          x = XB(I)
          DO J = 1, NQQ
            Q = QQ(J)
       	    print *, "KR",x, Q
            call PDF_A_nCTEQ15(x,Q,A,Z,U,D,UB,DB,S,SB,G,nloops)
            WRITE(4,103) x*SB,x*UB,x*DB,x*D,x*U,x*S,x*G
          ENDDO
        ENDDO


        WRITE(4,*) '---'
        CLOSE(4)
      ENDIF

C-------XE

      IF (I_XE.eq.1) THEN
      A = 131
      Z = 54

        INQUIRE(FILE=nCTEQ15_XE_DIR, EXIST=DIRECTORY_EXISTS)
        IF (.not.DIRECTORY_EXISTS) THEN
            CALL SYSTEM(MKDIR//nCTEQ15_XE_DIR)
        ENDIF
*-------Generate LHAPDF Grids
        OPEN(UNIT = 4, FILE = trim(nCTEQ15_XE_DIR)//nCTEQ15_XE_GRIDS)
        WRITE(4,*) 'PdfType: central'
        WRITE(4,*) 'Format: lhagrid1'
        WRITE(4,*) '---'
        WRITE(4,100) XB
        WRITE(4,101) QQ
        WRITE(4,102) PARTICLE_ID(-3:-1),PARTICLE_ID(1:3),PARTICLE_ID(21)

        DO I = 1, NXB
          x = XB(I)
          DO J = 1, NQQ
            Q = QQ(J)
       	    print *, "XE",x, Q
            call PDF_A_nCTEQ15(x,Q,A,Z,U,D,UB,DB,S,SB,G,nloops)
            WRITE(4,103) x*SB,x*UB,x*DB,x*D,x*U,x*S,x*G
          ENDDO
        ENDDO


        WRITE(4,*) '---'
        CLOSE(4)
      ENDIF

C-------BE

      IF (I_BE.eq.1) THEN
      A = 9
      Z = 4

        INQUIRE(FILE=nCTEQ15_BE_DIR, EXIST=DIRECTORY_EXISTS)
        IF (.not.DIRECTORY_EXISTS) THEN
            CALL SYSTEM(MKDIR//nCTEQ15_HE_DIR)
        ENDIF
*-------Generate LHAPDF Grids
        OPEN(UNIT = 4, FILE = trim(nCTEQ15_BE_DIR)//nCTEQ15_BE_GRIDS)
        WRITE(4,*) 'PdfType: central'
        WRITE(4,*) 'Format: lhagrid1'
        WRITE(4,*) '---'
        WRITE(4,100) XB
        WRITE(4,101) QQ
        WRITE(4,102) PARTICLE_ID(-3:-1),PARTICLE_ID(1:3),PARTICLE_ID(21)

        DO I = 1, NXB
          x = XB(I)
          DO J = 1, NQQ
            Q = QQ(J)
       	    print *, "BE",x, Q
            call PDF_A_nCTEQ15(x,Q,A,Z,U,D,UB,DB,S,SB,G,nloops)
            WRITE(4,103) x*SB,x*UB,x*DB,x*D,x*U,x*S,x*G
          ENDDO
        ENDDO


        WRITE(4,*) '---'
        CLOSE(4)
      ENDIF

C-------FE

      IF (I_FE.eq.1) THEN
      A = 56
      Z = 26

        INQUIRE(FILE=nCTEQ15_FE_DIR, EXIST=DIRECTORY_EXISTS)
        IF (.not.DIRECTORY_EXISTS) THEN
            CALL SYSTEM(MKDIR//nCTEQ15_HE_DIR)
        ENDIF
*-------Generate LHAPDF Grids
        OPEN(UNIT = 4, FILE = trim(nCTEQ15_FE_DIR)//nCTEQ15_FE_GRIDS)
        WRITE(4,*) 'PdfType: central'
        WRITE(4,*) 'Format: lhagrid1'
        WRITE(4,*) '---'
        WRITE(4,100) XB
        WRITE(4,101) QQ
        WRITE(4,102) PARTICLE_ID(-3:-1),PARTICLE_ID(1:3),PARTICLE_ID(21)

        DO I = 1, NXB
          x = XB(I)
          DO J = 1, NQQ
            Q = QQ(J)
       	    print *, "FE",x, Q
            call PDF_A_nCTEQ15(x,Q,A,Z,U,D,UB,DB,S,SB,G,nloops)
            WRITE(4,103) x*SB,x*UB,x*DB,x*D,x*U,x*S,x*G
          ENDDO
        ENDDO


        WRITE(4,*) '---'
        CLOSE(4)
      ENDIF

C-------WW

      IF (I_WW.eq.1) THEN
      A = 184
      Z = 74

        INQUIRE(FILE=nCTEQ15_WW_DIR, EXIST=DIRECTORY_EXISTS)
        IF (.not.DIRECTORY_EXISTS) THEN
            CALL SYSTEM(MKDIR//nCTEQ15_HE_DIR)
        ENDIF
*-------Generate LHAPDF Grids
        OPEN(UNIT = 4, FILE = trim(nCTEQ15_WW_DIR)//nCTEQ15_WW_GRIDS)
        WRITE(4,*) 'PdfType: central'
        WRITE(4,*) 'Format: lhagrid1'
        WRITE(4,*) '---'
        WRITE(4,100) XB
        WRITE(4,101) QQ
        WRITE(4,102) PARTICLE_ID(-3:-1),PARTICLE_ID(1:3),PARTICLE_ID(21)

        DO I = 1, NXB
          x = XB(I)
          DO J = 1, NQQ
            Q = QQ(J)
       	    print *, "WW",x, Q
            call PDF_A_nCTEQ15(x,Q,A,Z,U,D,UB,DB,S,SB,G,nloops)
            WRITE(4,103) x*SB,x*UB,x*DB,x*D,x*U,x*S,x*G
          ENDDO
        ENDDO


        WRITE(4,*) '---'
        CLOSE(4)
      ENDIF

C-------JLAB C

      IF (I_JCC.eq.1) THEN
      A = 12
      Z = 6

        INQUIRE(FILE=nCTEQ15_JCC_DIR, EXIST=DIRECTORY_EXISTS)
        IF (.not.DIRECTORY_EXISTS) THEN
            CALL SYSTEM(MKDIR//nCTEQ15_JCC_DIR)
        ENDIF
*-------Generate LHAPDF Grids
        OPEN(UNIT = 4, FILE = trim(nCTEQ15_JCC_DIR)//nCTEQ15_JCC_GRIDS)
        WRITE(4,*) 'PdfType: central'
        WRITE(4,*) 'Format: lhagrid1'
        WRITE(4,*) '---'
        WRITE(4,100) XB
        WRITE(4,101) QQ
        WRITE(4,102) PARTICLE_ID(-3:-1),PARTICLE_ID(1:3),PARTICLE_ID(21)

        DO I = 1, NXB
          x = XB(I)
          DO J = 1, NQQ
            Q = QQ(J)
            print *, "CC", x, Q
            call PDF_A_nCTEQ15(x,Q,A,Z,U,D,UB,DB,S,SB,G,nloops)
            WRITE(4,103) x*SB,x*UB,x*DB,x*D,x*U,x*S,x*G
          ENDDO
        ENDDO


        WRITE(4,*) '---'
        CLOSE(4)
      ENDIF

C-------JLAB FE

      IF (I_JFE.eq.1) THEN
      A = 56
      Z = 26

        INQUIRE(FILE=nCTEQ15_JFE_DIR, EXIST=DIRECTORY_EXISTS)
        IF (.not.DIRECTORY_EXISTS) THEN
            CALL SYSTEM(MKDIR//nCTEQ15_JFE_DIR)
        ENDIF
*-------Generate LHAPDF Grids
        OPEN(UNIT = 4, FILE = trim(nCTEQ15_JFE_DIR)//nCTEQ15_JFE_GRIDS)
        WRITE(4,*) 'PdfType: central'
        WRITE(4,*) 'Format: lhagrid1'
        WRITE(4,*) '---'
        WRITE(4,100) XB
        WRITE(4,101) QQ
        WRITE(4,102) PARTICLE_ID(-3:-1),PARTICLE_ID(1:3),PARTICLE_ID(21)

        DO I = 1, NXB
          x = XB(I)
          DO J = 1, NQQ
            Q = QQ(J)
            print *, "FE", x, Q
            call PDF_A_nCTEQ15(x,Q,A,Z,U,D,UB,DB,S,SB,G,nloops)
            WRITE(4,103) x*SB,x*UB,x*DB,x*D,x*U,x*S,x*G
          ENDDO
        ENDDO


        WRITE(4,*) '---'
        CLOSE(4)
      ENDIF

C-------JLAB PB

      IF (I_JPB.eq.1) THEN
      A = 208
      Z = 82

        INQUIRE(FILE=nCTEQ15_JPB_DIR, EXIST=DIRECTORY_EXISTS)
        IF (.not.DIRECTORY_EXISTS) THEN
            CALL SYSTEM(MKDIR//nCTEQ15_JPB_DIR)
        ENDIF
*-------Generate LHAPDF Grids
        OPEN(UNIT = 4, FILE = trim(nCTEQ15_JPB_DIR)//nCTEQ15_JPB_GRIDS)
        WRITE(4,*) 'PdfType: central'
        WRITE(4,*) 'Format: lhagrid1'
        WRITE(4,*) '---'
        WRITE(4,100) XB
        WRITE(4,101) QQ
        WRITE(4,102) PARTICLE_ID(-3:-1),PARTICLE_ID(1:3),PARTICLE_ID(21)

        DO I = 1, NXB
          x = XB(I)
          DO J = 1, NQQ
            Q = QQ(J)
            print *, "PB", x, Q
            call PDF_A_nCTEQ15(x,Q,A,Z,U,D,UB,DB,S,SB,G,nloops)
            WRITE(4,103) x*SB,x*UB,x*DB,x*D,x*U,x*S,x*G
          ENDDO
        ENDDO


        WRITE(4,*) '---'
        CLOSE(4)
      ENDIF

C-------CMS PB

      IF (I_PB.eq.1) THEN
      A = 209
      Z = 82

        INQUIRE(FILE=nCTEQ15_LPB_DIR, EXIST=DIRECTORY_EXISTS)
        IF (.not.DIRECTORY_EXISTS) THEN
            CALL SYSTEM(MKDIR//nCTEQ15_LPB_DIR)
        ENDIF
*-------Generate LHAPDF Grids
        OPEN(UNIT = 4, FILE = trim(nCTEQ15_LPB_DIR)//nCTEQ15_LPB_GRIDS)
        WRITE(4,*) 'PdfType: central'
        WRITE(4,*) 'Format: lhagrid1'
        WRITE(4,*) '---'
        WRITE(4,100) XB
        WRITE(4,101) QQ
        WRITE(4,102) PARTICLE_ID(-3:-1),PARTICLE_ID(1:3),PARTICLE_ID(21)

        DO I = 1, NXB
          x = XB(I)
          DO J = 1, NQQ
            Q = QQ(J)
            print *, "PB", x, Q
            call PDF_A_nCTEQ15(x,Q,A,Z,U,D,UB,DB,S,SB,G,nloops)
            WRITE(4,103) x*SB,x*UB,x*DB,x*D,x*U,x*S,x*G
          ENDDO
        ENDDO


        WRITE(4,*) '---'
        CLOSE(4)
      ENDIF

      IF (I_AU.eq.1) THEN
      A = 197
      Z = 79

        INQUIRE(FILE=nCTEQ15_AU_DIR, EXIST=DIRECTORY_EXISTS)
        IF (.not.DIRECTORY_EXISTS) THEN
            CALL SYSTEM(MKDIR//nCTEQ15_AU_DIR)
        ENDIF
*-------Generate LHAPDF Grids
        OPEN(UNIT = 4, FILE = trim(nCTEQ15_AU_DIR)//nCTEQ15_AU_GRIDS)
        WRITE(4,*) 'PdfType: central'
        WRITE(4,*) 'Format: lhagrid1'
        WRITE(4,*) '---'
        WRITE(4,100) XB
        WRITE(4,101) QQ
        WRITE(4,102) PARTICLE_ID(-3:-1),PARTICLE_ID(1:3),PARTICLE_ID(21)

        DO I = 1, NXB
          x = XB(I)
          DO J = 1, NQQ
            Q = QQ(J)
            print *, "AU", x, Q
            call PDF_A_nCTEQ15(x,Q,A,Z,U,D,UB,DB,S,SB,G,nloops)
            WRITE(4,103) x*SB,x*UB,x*DB,x*D,x*U,x*S,x*G
          ENDDO
        ENDDO


        WRITE(4,*) '---'
        CLOSE(4)
      ENDIF

*------FORMATTING
  100 FORMAT(240(E12.6E2,' '))
!  100 FORMAT(105(E12.6E2,' '))
  104 FORMAT(125(E12.6E2,' '))
  101 FORMAT(51(E12.6E2,' '))
  102 FORMAT(' ',3(I2,' '),3(I1,' '),I2)
  103 FORMAT(7(E12.6E2,'  '))

      return
      end
