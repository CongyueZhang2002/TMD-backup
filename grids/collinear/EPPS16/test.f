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
      REAL*8 Xgen
      REAL*8 X(998)
      REAL*8 Q
      INTEGER PARTICLE_ID(-5:21)
      INTEGER A,Z
      REAL*8 U,D,UB,DB,SS,SB,GL
      REAL*8 Ut,Dt,UBt,DBt,SSt,SBt,GLt
      INTEGER NLOOPS
      REAL*8 pdf(-6:6)
C----- Set nloops
      nloops = 1


C----- This program tests the generated LHAPDF grids

C-----SET UP LHAPDF PROTON PDF GRIDS (CT14NLO)
      CALL InitPDFsetByNameM(1,"CT14nlo")
      CALL InitPDFM(1,0)

C---- SET UP HE PDF GRIDS (EPPS16)
      CALL InitPDFsetByNameM(4,"EPPS16nlo_CT14nlo_He4")
      CALL InitPDFM(4,0)

C---- SET UP NE PDF GRIDS (EPPS16)
      CALL InitPDFsetByNameM(20,"EPPS16nlo_CT14nlo_Ne20")
      CALL InitPDFM(20,0)

C---- SET UP KR PDF GRIDS (EPPS16)
      CALL InitPDFsetByNameM(84,"EPPS16nlo_CT14nlo_Kr84")
      CALL InitPDFM(84,0)

C---- SET UP XE PDF GRIDS (EPPS16)
      CALL InitPDFsetByNameM(131,"EPPS16nlo_CT14nlo_Xe131")
      CALL InitPDFM(131,0)

C-----PARTICLE ID NUMBERS: Based on Standard PDG ID Numbers:
C------https://pdg.lbl.gov/2006/reviews/pdf-files/montecarlo-web.pdf
      PARTICLE_ID( 1) =  1  ! d
      PARTICLE_ID(-1) = -1  ! dbar
      PARTICLE_ID( 2) =  2  ! u
      PARTICLE_ID(-2) = -2  ! ubar
      PARTICLE_ID( 3) =  3  ! s
      PARTICLE_ID(-3) = -3  ! sbar
      !PARTICLE_ID( 4) =  4  ! c
      !PARTICLE_ID(-4) = -4  ! cbar
      !PARTICLE_ID( 5) =  5  ! b
      !PARTICLE_ID(-5) = -5  ! bbar
      PARTICLE_ID(21) = 21  ! g

      Xgen = 0d0
      DO I = 1, 998, 1
        Xgen = Xgen + 0.001d0
        X(I) = Xgen
      ENDDO

C-----TEST 1: Q = 2
      Q = 2d0

      OPEN(UNIT = 4, FILE = 'plot_data/He_Q2.dat')
      WRITE(4,*) 'x ', 'ut ', 'dt ', 'st ', 'ubt ', 'dbt ', 'sbt ',
     +           'glt ','u ', 'd ', 's ', 'ub ', 'db ', 'sb ', 'gl '
      A = 4
      Z  = 2
      DO I = 1, 998, 1
        CALL evolvePDFM(4,X(I),Q,pdf)
        Ut=pdf(2)
        Dt=pdf(1)
        SSt=pdf(3)
        UBt=pdf(-2)
        DBt=pdf(-1)
        SBt=pdf(-3)
        GLt=pdf(0)
        CALL PDF_A_EPPS16(x(I),Q,A,Z,U,D,UB,DB,SS,SB,GL,1)
        WRITE(4,*) x(I), Ut, Dt, SSt, UBt, DBt, SBt, GLt,
     +   U*x(I),D*x(I),SS*x(I),UB*x(I),DB*x(I),SB*x(I),GL*x(I)

        print *, UT/x(I)/U, UT/U
      ENDDO
      CLOSE(4)


      OPEN(UNIT = 20, FILE = 'plot_data/Ne_Q2.dat')
      WRITE(20,*) 'x ', 'ut ', 'dt ', 'st ', 'ubt ', 'dbt ', 'sbt ',
     +           'glt ','u ', 'd ', 's ', 'ub ', 'db ', 'sb ', 'gl '
      A = 20
      Z  = 10
      DO I = 1, 998, 1
        CALL evolvePDFM(20,X(I),Q,pdf)
        Ut=pdf(2)
        Dt=pdf(1)
        SSt=pdf(3)
        UBt=pdf(-2)
        DBt=pdf(-1)
        SBt=pdf(-3)
        GLt=pdf(0)
        CALL PDF_A_EPPS16(x(I),Q,A,Z,U,D,UB,DB,SS,SB,GL,1)
        WRITE(20,*) x(I), Ut, Dt, SSt, UBt, DBt, SBt, GLt,
     +   U*x(I),D*x(I),SS*x(I),UB*x(I),DB*x(I),SB*x(I),GL*x(I)
      ENDDO
      CLOSE(20)

      OPEN(UNIT = 84, FILE = 'plot_data/Kr_Q2.dat')
      WRITE(84,*) 'x ', 'ut ', 'dt ', 'st ', 'ubt ', 'dbt ', 'sbt ',
     +           'glt ','u ', 'd ', 's ', 'ub ', 'db ', 'sb ', 'gl '
      A = 84
      Z  = 36
      DO I = 1, 998, 1
        CALL evolvePDFM(84,X(I),Q,pdf)
        Ut=pdf(2)
        Dt=pdf(1)
        SSt=pdf(3)
        UBt=pdf(-2)
        DBt=pdf(-1)
        SBt=pdf(-3)
        GLt=pdf(0)
        CALL PDF_A_EPPS16(x(I),Q,A,Z,U,D,UB,DB,SS,SB,GL,1)
        WRITE(84,*) x(I), Ut, Dt, SSt, UBt, DBt, SBt, GLt,
     +   U*x(I),D*x(I),SS*x(I),UB*x(I),DB*x(I),SB*x(I),GL*x(I)
      ENDDO
      CLOSE(84)


      OPEN(UNIT = 131, FILE = 'plot_data/Xe_Q2.dat')
      WRITE(131,*) 'x ', 'ut ', 'dt ', 'st ', 'ubt ', 'dbt ', 'sbt ',
     +           'glt ','u ', 'd ', 's ', 'ub ', 'db ', 'sb ', 'gl '
      A = 131
      Z  = 54
      DO I = 1, 998, 1
        CALL evolvePDFM(131,X(I),Q,pdf)
        Ut=pdf(2)
        Dt=pdf(1)
        SSt=pdf(3)
        UBt=pdf(-2)
        DBt=pdf(-1)
        SBt=pdf(-3)
        GLt=pdf(0)
        CALL PDF_A_EPPS16(x(I),Q,A,Z,U,D,UB,DB,SS,SB,GL,1)
        WRITE(131,*) x(I), Ut, Dt, SSt, UBt, DBt, SBt, GLt,
     +   U*x(I),D*x(I),SS*x(I),UB*x(I),DB*x(I),SB*x(I),GL*x(I)
      ENDDO
      CLOSE(131)

C-----TEST 2: Q = 5
      Q = 5d0

      OPEN(UNIT = 4, FILE = 'plot_data/He_Q5.dat')
      WRITE(4,*) 'x ', 'ut ', 'dt ', 'st ', 'ubt ', 'dbt ', 'sbt ',
     +           'glt ','u ', 'd ', 's ', 'ub ', 'db ', 'sb ', 'gl '
      A = 4
      Z  = 2
      DO I = 1, 998, 1
        CALL evolvePDFM(4,X(I),Q,pdf)
        Ut=pdf(2)
        Dt=pdf(1)
        SSt=pdf(3)
        UBt=pdf(-2)
        DBt=pdf(-1)
        SBt=pdf(-3)
        GLt=pdf(0)
        CALL PDF_A_EPPS16(x(I),Q,A,Z,U,D,UB,DB,SS,SB,GL,1)
        WRITE(4,*) x(I), Ut, Dt, SSt, UBt, DBt, SBt, GLt,
     +   U*x(I),D*x(I),SS*x(I),UB*x(I),DB*x(I),SB*x(I),GL*x(I)
      ENDDO
      CLOSE(4)


      OPEN(UNIT = 20, FILE = 'plot_data/Ne_Q5.dat')
      WRITE(20,*) 'x ', 'ut ', 'dt ', 'st ', 'ubt ', 'dbt ', 'sbt ',
     +           'glt ','u ', 'd ', 's ', 'ub ', 'db ', 'sb ', 'gl '
      A = 20
      Z  = 10
      DO I = 1, 998, 1
        CALL evolvePDFM(20,X(I),Q,pdf)
        Ut=pdf(2)
        Dt=pdf(1)
        SSt=pdf(3)
        UBt=pdf(-2)
        DBt=pdf(-1)
        SBt=pdf(-3)
        GLt=pdf(0)
        CALL PDF_A_EPPS16(x(I),Q,A,Z,U,D,UB,DB,SS,SB,GL,1)
        WRITE(20,*) x(I), Ut, Dt, SSt, UBt, DBt, SBt, GLt,
     +   U*x(I),D*x(I),SS*x(I),UB*x(I),DB*x(I),SB*x(I),GL*x(I)
      ENDDO
      CLOSE(20)

      OPEN(UNIT = 84, FILE = 'plot_data/Kr_Q5.dat')
      WRITE(84,*) 'x ', 'ut ', 'dt ', 'st ', 'ubt ', 'dbt ', 'sbt ',
     +           'glt ','u ', 'd ', 's ', 'ub ', 'db ', 'sb ', 'gl '
      A = 84
      Z  = 36
      DO I = 1, 998, 1
        CALL evolvePDFM(84,X(I),Q,pdf)
        Ut=pdf(2)
        Dt=pdf(1)
        SSt=pdf(3)
        UBt=pdf(-2)
        DBt=pdf(-1)
        SBt=pdf(-3)
        GLt=pdf(0)
        CALL PDF_A_EPPS16(x(I),Q,A,Z,U,D,UB,DB,SS,SB,GL,1)
        WRITE(84,*) x(I), Ut, Dt, SSt, UBt, DBt, SBt, GLt,
     +   U*x(I),D*x(I),SS*x(I),UB*x(I),DB*x(I),SB*x(I),GL*x(I)
      ENDDO
      CLOSE(84)


      OPEN(UNIT = 131, FILE = 'plot_data/Xe_Q5.dat')
      WRITE(131,*) 'x ', 'ut ', 'dt ', 'st ', 'ubt ', 'dbt ', 'sbt ',
     +           'glt ','u ', 'd ', 's ', 'ub ', 'db ', 'sb ', 'gl '
      A = 131
      Z  = 54
      DO I = 1, 998, 1
        CALL evolvePDFM(131,X(I),Q,pdf)
        Ut=pdf(2)
        Dt=pdf(1)
        SSt=pdf(3)
        UBt=pdf(-2)
        DBt=pdf(-1)
        SBt=pdf(-3)
        GLt=pdf(0)
        CALL PDF_A_EPPS16(x(I),Q,A,Z,U,D,UB,DB,SS,SB,GL,1)
        WRITE(131,*) x(I), Ut, Dt, SSt, UBt, DBt, SBt, GLt,
     +   U*x(I),D*x(I),SS*x(I),UB*x(I),DB*x(I),SB*x(I),GL*x(I)
      ENDDO
      CLOSE(131)

C-----TEST 3: Q = 10
      Q = 10d0

      OPEN(UNIT = 4, FILE = 'plot_data/He_Q10.dat')
      WRITE(4,*) 'x ', 'ut ', 'dt ', 'st ', 'ubt ', 'dbt ', 'sbt ',
     +           'glt ','u ', 'd ', 's ', 'ub ', 'db ', 'sb ', 'gl '
      A = 4
      Z  = 2
      DO I = 1, 998, 1
        CALL evolvePDFM(4,X(I),Q,pdf)
        Ut=pdf(2)
        Dt=pdf(1)
        SSt=pdf(3)
        UBt=pdf(-2)
        DBt=pdf(-1)
        SBt=pdf(-3)
        GLt=pdf(0)
        CALL PDF_A_EPPS16(x(I),Q,A,Z,U,D,UB,DB,SS,SB,GL,1)
        WRITE(4,*) x(I), Ut, Dt, SSt, UBt, DBt, SBt, GLt,
     +   U*x(I),D*x(I),SS*x(I),UB*x(I),DB*x(I),SB*x(I),GL*x(I)
      ENDDO
      CLOSE(4)


      OPEN(UNIT = 20, FILE = 'plot_data/Ne_Q10.dat')
      WRITE(20,*) 'x ', 'ut ', 'dt ', 'st ', 'ubt ', 'dbt ', 'sbt ',
     +           'glt ','u ', 'd ', 's ', 'ub ', 'db ', 'sb ', 'gl '
      A = 20
      Z  = 10
      DO I = 1, 998, 1
        CALL evolvePDFM(20,X(I),Q,pdf)
        Ut=pdf(2)
        Dt=pdf(1)
        SSt=pdf(3)
        UBt=pdf(-2)
        DBt=pdf(-1)
        SBt=pdf(-3)
        GLt=pdf(0)
        CALL PDF_A_EPPS16(x(I),Q,A,Z,U,D,UB,DB,SS,SB,GL,1)
        WRITE(20,*) x(I), Ut, Dt, SSt, UBt, DBt, SBt, GLt,
     +   U*x(I),D*x(I),SS*x(I),UB*x(I),DB*x(I),SB*x(I),GL*x(I)
      ENDDO
      CLOSE(20)

      OPEN(UNIT = 84, FILE = 'plot_data/Kr_Q10.dat')
      WRITE(84,*) 'x ', 'ut ', 'dt ', 'st ', 'ubt ', 'dbt ', 'sbt ',
     +           'glt ','u ', 'd ', 's ', 'ub ', 'db ', 'sb ', 'gl '
      A = 84
      Z  = 36
      DO I = 1, 998, 1
        CALL evolvePDFM(84,X(I),Q,pdf)
        Ut=pdf(2)
        Dt=pdf(1)
        SSt=pdf(3)
        UBt=pdf(-2)
        DBt=pdf(-1)
        SBt=pdf(-3)
        GLt=pdf(0)
        CALL PDF_A_EPPS16(x(I),Q,A,Z,U,D,UB,DB,SS,SB,GL,1)
        WRITE(84,*) x(I), Ut, Dt, SSt, UBt, DBt, SBt, GLt,
     +   U*x(I),D*x(I),SS*x(I),UB*x(I),DB*x(I),SB*x(I),GL*x(I)
      ENDDO
      CLOSE(84)


      OPEN(UNIT = 131, FILE = 'plot_data/Xe_Q10.dat')
      WRITE(131,*) 'x ', 'ut ', 'dt ', 'st ', 'ubt ', 'dbt ', 'sbt ',
     +           'glt ','u ', 'd ', 's ', 'ub ', 'db ', 'sb ', 'gl '
      A = 131
      Z  = 54
      DO I = 1, 998, 1
        CALL evolvePDFM(131,X(I),Q,pdf)
        Ut=pdf(2)
        Dt=pdf(1)
        SSt=pdf(3)
        UBt=pdf(-2)
        DBt=pdf(-1)
        SBt=pdf(-3)
        GLt=pdf(0)
        CALL PDF_A_EPPS16(x(I),Q,A,Z,U,D,UB,DB,SS,SB,GL,1)
        WRITE(131,*) x(I), Ut, Dt, SSt, UBt, DBt, SBt, GLt,
     +   U*x(I),D*x(I),SS*x(I),UB*x(I),DB*x(I),SB*x(I),GL*x(I)
      ENDDO
      CLOSE(131)

C-----TEST 4: Q = 50
      Q = 50d0

      OPEN(UNIT = 4, FILE = 'plot_data/He_Q50.dat')
      WRITE(4,*) 'x ', 'ut ', 'dt ', 'st ', 'ubt ', 'dbt ', 'sbt ',
     +           'glt ','u ', 'd ', 's ', 'ub ', 'db ', 'sb ', 'gl '
      A = 4
      Z  = 2
      DO I = 1, 998, 1
        CALL evolvePDFM(4,X(I),Q,pdf)
        Ut=pdf(2)
        Dt=pdf(1)
        SSt=pdf(3)
        UBt=pdf(-2)
        DBt=pdf(-1)
        SBt=pdf(-3)
        GLt=pdf(0)
        CALL PDF_A_EPPS16(x(I),Q,A,Z,U,D,UB,DB,SS,SB,GL,1)
        WRITE(4,*) x(I), Ut, Dt, SSt, UBt, DBt, SBt, GLt,
     +   U*x(I),D*x(I),SS*x(I),UB*x(I),DB*x(I),SB*x(I),GL*x(I)
      ENDDO
      CLOSE(4)


      OPEN(UNIT = 20, FILE = 'plot_data/Ne_Q50.dat')
      WRITE(20,*) 'x ', 'ut ', 'dt ', 'st ', 'ubt ', 'dbt ', 'sbt ',
     +           'glt ','u ', 'd ', 's ', 'ub ', 'db ', 'sb ', 'gl '
      A = 20
      Z  = 10
      DO I = 1, 998, 1
        CALL evolvePDFM(20,X(I),Q,pdf)
        Ut=pdf(2)
        Dt=pdf(1)
        SSt=pdf(3)
        UBt=pdf(-2)
        DBt=pdf(-1)
        SBt=pdf(-3)
        GLt=pdf(0)
        CALL PDF_A_EPPS16(x(I),Q,A,Z,U,D,UB,DB,SS,SB,GL,1)
        WRITE(20,*) x(I), Ut, Dt, SSt, UBt, DBt, SBt, GLt,
     +   U*x(I),D*x(I),SS*x(I),UB*x(I),DB*x(I),SB*x(I),GL*x(I)
      ENDDO
      CLOSE(20)

      OPEN(UNIT = 84, FILE = 'plot_data/Kr_Q50.dat')
      WRITE(84,*) 'x ', 'ut ', 'dt ', 'st ', 'ubt ', 'dbt ', 'sbt ',
     +           'glt ','u ', 'd ', 's ', 'ub ', 'db ', 'sb ', 'gl '
      A = 84
      Z  = 36
      DO I = 1, 998, 1
        CALL evolvePDFM(84,X(I),Q,pdf)
        Ut=pdf(2)
        Dt=pdf(1)
        SSt=pdf(3)
        UBt=pdf(-2)
        DBt=pdf(-1)
        SBt=pdf(-3)
        GLt=pdf(0)
        CALL PDF_A_EPPS16(x(I),Q,A,Z,U,D,UB,DB,SS,SB,GL,1)
        WRITE(84,*) x(I), Ut, Dt, SSt, UBt, DBt, SBt, GLt,
     +   U*x(I),D*x(I),SS*x(I),UB*x(I),DB*x(I),SB*x(I),GL*x(I)
      ENDDO
      CLOSE(84)


      OPEN(UNIT = 131, FILE = 'plot_data/Xe_Q50.dat')
      WRITE(131,*) 'x ', 'ut ', 'dt ', 'st ', 'ubt ', 'dbt ', 'sbt ',
     +           'glt ','u ', 'd ', 's ', 'ub ', 'db ', 'sb ', 'gl '
      A = 131
      Z  = 54
      DO I = 1, 998, 1
        CALL evolvePDFM(131,X(I),Q,pdf)
        Ut=pdf(2)
        Dt=pdf(1)
        SSt=pdf(3)
        UBt=pdf(-2)
        DBt=pdf(-1)
        SBt=pdf(-3)
        GLt=pdf(0)
        CALL PDF_A_EPPS16(x(I),Q,A,Z,U,D,UB,DB,SS,SB,GL,1)
        WRITE(131,*) x(I), Ut, Dt, SSt, UBt, DBt, SBt, GLt,
     +   U*x(I),D*x(I),SS*x(I),UB*x(I),DB*x(I),SB*x(I),GL*x(I)
      ENDDO
      CLOSE(131)



      END PROGRAM MASTER
