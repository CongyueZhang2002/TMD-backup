      program example
      implicit none
      real*8 x,Q, z
      integer Evolve3D_TMDPDF, Evolve3D_TMDFF
      integer plot3D_TMDPDF, plot3D_TMDFF
      real*8 up,ubp,dp,dbp,sp,sbp,glp
      real*8 u1,ub1,d1,db1,s1,sb1,gl1
      real*8 u2,ub2,d2,db2,s2,sb2,gl2
      integer i,j,k0, k1,k2,k3,k4,k5,k6,set
      integer usenc
      Integer nucleon, Nucleus1, Nucleus2
      common /nckern/ usenc
      integer nloops,hop,nll
      common /scheme/ nloops,hop,nll
      integer REPLICA_NUM
      COMMON /rep_num/ REPLICA_NUM
      integer EIGENSET
      real*8 kt

C-----perturbative orders:
      NLOOPS = 2
      HOP = 0
      NLL = 3

C-----plot3D = 1 : generate data for TMD 3D plots (vs. kt and x) (slow)
C-----plot3D = 0 : do not generate data for TMD 3D plots (fast)
      plot3D_TMDPDF = 0
      plot3D_TMDFF = 0

C-----Evolveplot3D = 1 : generate data for 3D TMD Evolution plots (vs. kt and Q) (slow)
C-----plot3D = 0 : do not generate data for 3D TMD Evolution plots (fast)
      Evolve3D_TMDPDF = 0
      Evolve3D_TMDFF = 0

C-----AVAILABLE NUCLEI FOR OUTPUTTING nTMDPDFs:
C----- HE (A=4)  , BE (A=9),  NE (A=20),  FE (A=56), KR(A=84)
C----- XE (A=131), W (A=184)

C-----Set Nuclei by MASS NUMBER:
      nucleon = 1
      Nucleus1 = 56
      Nucleus2 = 197

      SET = 1
      call init_nTMDs(set)

      IF(REPLICA_NUM.le.40) then
      EIGENSET = 0
      ELSE
      EIGENSET = REPLICA_NUM-40
      ENDIF
      IF(REPLICA_NUM.gt.96) then
      EIGENSET = 0
      ENDIF

      print*, "REPLICA_NUM,  EIGENSET:", REPLICA_NUM, EIGENSET

C-----INITIALIZE PDF SETS USING LHAPDF
      print *, "1"
      CALL InitPDFsetByNameM(1,"CT14nlo")
      CALL InitPDFM(1,EIGENSET)
      !CALL InitPDFsetByNameM(2,"EPPS16nlo_CT14nlo_He4")
      !CALL InitPDFM(2,REPLICA_NUM)
      CALL InitPDFsetByNameM(3,"EPPS16nlo_CT14nlo_Be9")
      CALL InitPDFM(3,REPLICA_NUM)
      CALL InitPDFsetByNameM(4,"EPPS16nlo_CT14nlo_C12")
      CALL InitPDFM(4,REPLICA_NUM)
      CALL InitPDFsetByNameM(5,"EPPS16nlo_CT14nlo_Fe56")
      CALL InitPDFM(5,REPLICA_NUM)
      CALL InitPDFsetByNameM(6,"EPPS16nlo_CT14nlo_W184")
      CALL InitPDFM(6,REPLICA_NUM)
      CALL InitPDFsetByNameM(7,"EPPS16nlo_CT14nlo_Au197")
      CALL InitPDFM(7,REPLICA_NUM)

C-----Get TMDs:
      DO J = 1,1

      IF(J.EQ.1) THEN
      Q = sqrt(2.4d0) ! SCALE 1: Q = Q0
      ELSEIF(J.EQ.2) THEN
      Q = 10d0        ! SCALE 2: Q = 10 GeV
      ENDIF
C-----NUMBERING SYSTEM FOR PARTON TYPE FOR PDF_KT
C-----1: u   2:  d    3: s   4: ub   5: db   6: sb
C---------------------------------------------------
C-----Getting TMDPDFs vs. KT for u,d quark for HE and XE
      k0 = j
      WRITE(k0, *) 'kt up dp u1 d1 u2 d2'
      WRITE(k0, *) 0d0,0d0,0d0,0d0,0d0,0d0,0d0

      DO i = 1, 100
C-----KT RANGE:
      kt = 0.01*i  ! from pt = 0.01 to pt = 1.0 GeV
C-----Set x-value for nTMDPDFs:
      x = 0.004

      CALL PDF_kt(kt,x,Q,1,nucleon,up)
      CALL PDF_kt(kt,x,Q,-1,nucleon,dp)

      CALL PDF_kt(kt,x,Q,1,Nucleus1,u1)
      CALL PDF_kt(kt,x,Q,-1,nucleus1,d1)

      CALL PDF_kt(kt,x,Q,1,Nucleus2,u2)
      CALL PDF_kt(kt,x,Q,-1,nucleus2,d2)

      WRITE(k0,*) kt,up,dp, u1,d1,u2,d2
      ENDDO
      CLOSE(k0)


C-----Getting TMDPDFs vs. KT for u,d quark for HE and XE
      k1 = 100+j
      WRITE(k1, *) 'kt up dp u1 d1 u2 d2'
      WRITE(k1, *) 0d0,0d0,0d0,0d0,0d0,0d0,0d0

      DO i = 1, 100
C-----KT RANGE:
      kt = 0.01*i  ! from pt = 0.01 to pt = 1.0 GeV
C-----Set x-value for nTMDPDFs:
      x = 0.04d0

      CALL PDF_kt(kt,x,Q,1,nucleon,up)
      CALL PDF_kt(kt,x,Q,-1,nucleon,dp)

      CALL PDF_kt(kt,x,Q,1,Nucleus1,u1)
      CALL PDF_kt(kt,x,Q,-1,nucleus1,d1)

      CALL PDF_kt(kt,x,Q,1,Nucleus2,u2)
      CALL PDF_kt(kt,x,Q,-1,nucleus2,d2)

      WRITE(k1,*) kt,up,dp, u1,d1,u2,d2
      ENDDO
      CLOSE(k1)

C-----Getting TMDPDFs vs. KT for u,d quark for HE and XE
      k0 = 1000*j
      WRITE(k0, *) 'x up dp u1 d1 u2 d2'
      WRITE(k0, *) 0d0,0d0,0d0,0d0,0d0,0d0,0d0

      DO i = 1, 150
C-----KT RANGE:
      kt = 0.01 !
C-----Set x-value for nTMDPDFs:
      x = 0.004*i !from x = 0.01 to x = 0.6

      CALL PDF_kt(kt,x,Q,1,nucleon,up)
      CALL PDF_kt(kt,x,Q,-1,nucleon,dp)

      CALL PDF_kt(kt,x,Q,1,Nucleus1,u1)
      CALL PDF_kt(kt,x,Q,-1,nucleus1,d1)

      CALL PDF_kt(kt,x,Q,1,Nucleus2,u2)
      CALL PDF_kt(kt,x,Q,-1,nucleus2,d2)

      WRITE(k0,*) x,up,dp, u1,d1,u2,d2
      ENDDO
      CLOSE(k0)

      k3 = 300+j
      WRITE(k3, *) 'kt up dp u1 d1 u2 d2'
      WRITE(k3, *) 0d0,0d0,0d0,0d0,0d0,0d0,0d0

      DO i = 1, 100
C-----KT RANGE:
      kt = 0.01*i  ! from pt = 0.01 to pt = 1.0 GeV
C-----Set x-value for nTMDPDFs:
      x = 0.40d0

      CALL PDF_kt(kt,x,Q,1,nucleon,up)
      CALL PDF_kt(kt,x,Q,-1,nucleon,dp)

      CALL PDF_kt(kt,x,Q,1,Nucleus1,u1)
      CALL PDF_kt(kt,x,Q,-1,nucleus1,d1)

      CALL PDF_kt(kt,x,Q,1,Nucleus2,u2)
      CALL PDF_kt(kt,x,Q,-1,nucleus2,d2)

      WRITE(k3,*) kt,up,dp, u1,d1,u2,d2
      ENDDO
      CLOSE(k3)

      k5 = 500+j
      WRITE(k5, *) 'kt up dp u1 d1 u2 d2'
      WRITE(k5, *) 0d0,0d0,0d0,0d0,0d0,0d0,0d0

      DO i = 1, 100
C-----KT RANGE:
      kt = 0.01*i  ! from pt = 0.01 to pt = 1.0 GeV
C-----Set x-value for nTMDPDFs:
      x = 0.60d0

      CALL PDF_kt(kt,x,Q,1,nucleon,up)
      CALL PDF_kt(kt,x,Q,-1,nucleon,dp)

      CALL PDF_kt(kt,x,Q,1,Nucleus1,u1)
      CALL PDF_kt(kt,x,Q,-1,nucleus1,d1)

      CALL PDF_kt(kt,x,Q,1,Nucleus2,u2)
      CALL PDF_kt(kt,x,Q,-1,nucleus2,d2)

      WRITE(k5,*) kt,up,dp, u1,d1,u2,d2
      ENDDO
      CLOSE(k5)


      ENDDO






      END PROGRAM

      subroutine init_nTMDs(set)
      implicit none
      integer set, I
      real*8 Ans ,Bns, GDs, ADs, GPs, APs
      common /fitp/
     >       Ans ,Bns ,GDs, Ads, GPs, APs
      real*8 ANA(101), BNA(101), GDA(101), ADA(101), GPA(101), APA(101)
      COMMON /param_dat/ ANA, BNA, GDA, ADA, GPA, APA

      CALL READDATA

      ANs =       ANA(set)
      BNs =       BNA(set)
      GPs =       GDA(set) ! 0d0
      ADs =       ADA(set) ! 0d0
      GDs =       0d0
      APs =       0d0

      !DO I=1,101
         !PRINT*, I, ANA(I), BNA(I), GBA(I)
      !END DO

      end
