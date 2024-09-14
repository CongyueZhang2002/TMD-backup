      PROGRAM MASTER
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL FCNG1
      DIMENSION NPRM(6),VSTRT(6),STP(6),BL(6),BUP(6),ARGLIS(6)
      CHARACTER*10 PNAM(6)
      INTEGER IRD,IWR,ISAV
      INTEGER NFIT
      INTEGER DAY, HOUR, MINUTE, SECOND
      CHARACTER NAME*50
      COMMON /PARAMS/ NFIT
      INTEGER IFLAG
      COMMON  / IEND / IFLAG

C     INITALIZE THE PARAMETRIZATION
      REAL*8 aN, bN, g2A, A2, g2B, B2
      COMMON /FITP/ aN, bN, g2A, A2, g2B, B2

      data NPRM /1, 2, 3, 4, 5, 6/
      data PNAM /'aN', 'bN', 'g2A', 'A2', 'g2B', 'B2'/

C-----STARTING STEP SIZE OR APPROXIMATE PARAMETER ERROR.
      data STP /0.01, 0.01, 0.01, 0.01, 0.01, 0.01/

C-----LOWER BOUND (LIMIT) ON PARAMETER VALUE, IF ANY
      data BL   /-0.086, -0.034, -10.0, -10.0, -10.0, -10.0/

C-----UPPER BOUND (LIMIT) ON PARAMETER VALUE, IF ANY
      data BUP  /0.5,  0.5,  10.0,  10.0,  10.0,  10.0/

C-----Local Variable Declarations
      include "tools/data-inc.f"
      real*8 rts,y,Q,Qmin,Qmax,pt,fuu, Rds, sigma, x2, ds,y2
      real*8 dcorr,fN
      real*8 fuuA,fuuB
      real*8 R_dy
      real*8 xfmin, xfmax, ymin, ymax, Qbar, const
      integer nloops,hop,nll, pre, prepdf,preff
      INTEGER REPLICA_NUM
      COMMON /REPLICA_INFO/ REPLICA_NUM
      INTEGER EIGENSET
      common /scheme/ nloops,hop,nll,pre
      real*8 Q0
      character*2 proc
      integer IH,IC
      real*8 AA,ZZ
      common /meson/ IH,IC
      REAL*8 ANMAX  ,ANMIN
      REAL*8 BNMAX  ,BNMIN
      REAL*8 GAMAX  ,GAMIN
      REAL*8 A2MAX  ,A2MIN
      REAL*8 GBMAX  ,GBMIN
      REAL*8 B2MAX  ,B2MIN
      integer seed,initseed
      REAL*8 XX(6)

      INTEGER NUM_EXP
      COMMON /NNUM_EXP/ NUM_EXP
      real*8 r0,r01000

      REAL*8 r1, r2 , r3, r4
      REAL*8 r5, r6

      REAL*8 ANBEST,BNBEST,GABEST,A2BEST,GBBEST,B2BEST


 1115 CONTINUE
C      INCLUDE "best-params-rand.f"

C      INCLUDE "job_output/best-params.f"

      INCLUDE "best-params.f"


C-----RHIC Ratios
      I_RHIC_Ratio_pAu1 = 0
      I_RHIC_Ratio_pAu2 = 0

C-----ATLAS 5 TEV
      I_ATLAS5_Y1 = 0
      I_ATLAS5_Y2 = 0
      I_ATLAS5_Y3 = 0

C-----CMS 5 TEV
      I_CMS5 = 0

C-----E866
      I_E866_800q  = 0

C-----E772
      I_E772_800   = 0

C-----SIDIS (HERMES)
      I_PIP_HE_Z   = 0
      I_PIP_NE_Z   = 0
      I_PIP_KR_Z   = 0
      I_PIP_XE_Z   = 0
      I_PI0_HE_Z   = 0
      I_PI0_NE_Z   = 0
      I_PI0_KR_Z   = 0
      I_PI0_XE_Z   = 0
      I_PIM_HE_Z   = 0
      I_PIM_NE_Z   = 0
      I_PIM_KR_Z   = 0
      I_PIM_XE_Z   = 0


      I_PIP_HE_PT2 = 0
      I_PIP_NE_PT2 = 0
      I_PIP_KR_PT2 = 0
      I_PIP_XE_PT2 = 0
      I_PI0_HE_PT2 = 0
      I_PI0_NE_PT2 = 0
      I_PI0_KR_PT2 = 0
      I_PI0_XE_PT2 = 0
      I_PIM_HE_PT2 = 0
      I_PIM_NE_PT2 = 0
      I_PIM_KR_PT2 = 0
      I_PIM_XE_PT2 = 0


      I_JLAB_PIP_12 = 0
      I_JLAB_PIM_12 = 0

      I_EIC_PRED = 1
      I_JLAB_PRED = 1



c-----FIT PARAMETERS
      nloops = 2 ! LO: 1 NLO: 2
      nll = 3
      pre = 1

!     CALL SETLHAPARM('SILENT') ! TO NOT SHOW THE CALLS, ALTHOUGH THEY ARE CALLED
      CALL SETLHAPARM('SILENT')

      CALL READDATA


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
      CALL InitPDFM(1,0)
      print *, "2"
      CALL InitPDFsetByNameM(2,"EPPS16HE")
      CALL InitPDFM(2,0)
      print *, "3"
      CALL InitPDFsetByNameM(3,"EPPS16NE")
      CALL InitPDFM(3,0)
      print *, "4"
      CALL InitPDFsetByNameM(4,"EPPS16KR")
      CALL InitPDFM(4,0)
      print *, "5"
      CALL InitPDFsetByNameM(5,"EPPS16XE")
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
      CALL InitPDFsetByNameM(10,"EPPS16BE")
      CALL InitPDFM(10,0)
      print *, "11"
      CALL InitPDFsetByNameM(11,"EPPS16FE")
      CALL InitPDFM(11,0)
      print *, "12"
      CALL InitPDFsetByNameM(12,"EPPS16WW")
      CALL InitPDFM(12,0)
      print *, "13"
      CALL InitPDFsetByNameM(13,"EPPS16JCC")
      CALL InitPDFM(13,0)
      print *, "14"
      CALL InitPDFsetByNameM(14,"EPPS16JFE")
      CALL InitPDFM(14,0)
      print *, "15"
      CALL InitPDFsetByNameM(15,"EPPS16JPB")
      CALL InitPDFM(15,0)
      print *, "16"
      CALL InitPDFsetByNameM(16,"LIKEnCC")
      CALL InitPDFM(16,0)
      print *, "17"
      CALL InitPDFsetByNameM(17,"LIKEnFE")
      CALL InitPDFM(17,0)
      print *, "18"
      CALL InitPDFsetByNameM(18,"LIKEnPB")
      CALL InitPDFM(18,0)
      print *, "20"
      CALL InitPDFsetByNameM(20,"EPPS16AU")
      CALL InitPDFM(20,0)
      print *, "21"
      CALL InitPDFsetByNameM(21,"EPPS16PR")
      CALL InitPDFM(21,0)
      print *, "22"
      CALL InitPDFsetByNameM(22,"EPPS16CA")
      CALL InitPDFM(22,0)
      print *, "23"
      CALL InitPDFsetByNameM(23,"LIKEnVC")
      CALL InitPDFM(23,0)
      CALL InitPDFsetByNameM(24,"LIKEnAU")
      CALL InitPDFM(24,0)


*.....CALCULATE THE STARTING TIME OF THE FIT
      CALL TIMESTAMP( )         ! WRITE OUT THE STARTING TIME
      CALL CPU_TIME(TIME1)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     START FIT
*.....INITIALIZATION:
      IRD = 5         ! UNIT NUMBER FOR INPUT TO MINUIT (5 KEYBOARD)
      IWR = 6         ! UNIT NUMBER FOR OUTPUT FROM MINUIT (6 SCREEN)
      ISAV = 7        ! UNIT NUMBER FOR USE OF THE SAVE COMMAND
      CALL MNINIT (IRD,IWR,ISAV)
*.....DEFINITON OF THE PARAMETERS :
      DO 11 I = 1, 2
         CALL MNPARM (NPRM(I),PNAM(I),VSTRT(I),
     >                        STP(I),BL(I),BUP(I),IERFLG)
         IF (IERFLG .NE. 0) THEN
            WRITE (6,*) 'UNABLE TO DEFINE PARAMETER NO.', I
            STOP
         END IF
 11   CONTINUE
*.....OUTPUT PRINT LEVEL (FROM -1 TO 3)
      ARGLIS(1) = 3.
      CALL MNEXCM (FCNG1,'SET PRINT',ARGLIS, 1, IERFLG, DUM)
*.....FIRST CALL :
!      ARGLIS(1) = 1.            !   IFLAG = 1
!      CALL MNEXCM (FCNG1,'CALL FCN', ARGLIS, 1, IERFLG, DUM)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     NUMBER OF PARAMETERS
      NFIT = 2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*.....SIMPLEX FIT :
!      ARGLIS(1) = 2000.
!      CALL MNEXCM (FCNG1, 'SIMPLEX', ARGLIS, 1, IERFLG, DUM)
*.....MINIMIZE FIT :
!     IT IS EQUIVALENT TO MIGRAD, EXCEPT THAT IF MIGRAD FAILS,
*     IT REVERTS TO SIMPLEX AND THEN CALLS MIGRAD AGAIN
C      ARGLIS(1) = 4000.
C      CALL MNEXCM (FCNG1, 'MINIMIZE', ARGLIS, 1, IERFLG, DUM)
*     CHECK THAT RESULT HAS CONVERGED
C      XX(1)  = AN
C      XX(2)  = BN
C      XX(3)  = G2A
C      XX(4)  = A2
C      XX(5)  = G2B
C      XX(6)  = B2
C      IF (CHISQUARE(XX)/NUM_EXP.GT.20D0) THEN
C      GOTO 1115
C      ELSE IF (ISNAN(CHISQUARE(XX))) THEN
C      GOTO 1115
C      ENDIF
*.....LAST CALL :
*     THIS CREATES THE FILE WITH THE RESULTS
      ARGLIS(1) = 3.            !   IFLAG = 3
      CALL MNEXCM (FCNG1, 'CALL FCN', ARGLIS, 1, IERFLG, DUM)
*.....STOP
      CALL MNEXCM (FCNG1, 'STOP', ARGLIS, 1, IERFLG, DUM)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C
*.....CALCULATE AND PRINT THE TOTAL TIME OF THE FIT
      CALL CPU_TIME(TIME2)
      TIMETOTAL = TIME2-TIME1
        DAY=INT(TIMETOTAL/86400.)
        HOUR=INT((TIMETOTAL-86400.*DAY)/3600.)
        MINUTE=INT((TIMETOTAL-86400.*DAY-3600.*HOUR)/60.)
        SECOND=INT(TIMETOTAL-86400.*DAY-3600.*HOUR-60.*MINUTE)
 127  FORMAT(A40,I2,A1,I2,A1,I2,A1,I2,A1)
C      WRITE(*,*) 'ELAPSED CPU TIME FOR THE FIT (TOTAL SECS):',
C     <              TIMETOTAL
      WRITE(*,127) 'ELAPSED CPU TIME FOR THE FIT (DDHHMMSS):',
     <              DAY,'D',HOUR,'H',MINUTE,'M',SECOND,'S'

       ANBEST  = AN
       BNBEST  = BN
       G2ABEST = G2A
       A2BEST  = A2
       G2BBEST = G2B
       B2BEST  = B2

       OPEN(UNIT=999,FILE='./best-params.f')
       WRITE(999, *) "      REAL*8 ANBEST"
       WRITE(999, *) "      REAL*8 BNBEST"
       WRITE(999, *) "      REAL*8 G2ABEST"
       WRITE(999, *) "      REAL*8 A2BEST"
       WRITE(999, *) "      REAL*8 G2BBEST"
       WRITE(999, *) "      REAL*8 B2BEST"
       WRITE(999, *) "      ANBEST =    ",  ANBEST
       WRITE(999, *) "      BNBEST =    ",  BNBEST
       WRITE(999, *) "      G2ABEST =    ", G2ABEST
       WRITE(999, *) "      A2BEST =    ",  A2BEST
       WRITE(999, *) "      G2BBEST =    ", G2BBEST
       WRITE(999, *) "      B2BEST =    ",  B2BEST
       WRITE(999, *) "      VSTRT(1)  =  ANBEST"
       WRITE(999, *) "      VSTRT(2)  =  BNBEST"
       WRITE(999, *) "      VSTRT(3)  =  G2ABEST"
       WRITE(999, *) "      VSTRT(4)  =  A2BEST"
       WRITE(999, *) "      VSTRT(5)  =  G2BBEST"
       WRITE(999, *) "      VSTRT(6)  =  B2BEST"


       CLOSE(999)

       RETURN
       END

      SUBROUTINE FCNG1 (NPAR, G, F, X, IIFLAG, dum)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION X(11), G(11)
      integer iflag,iiflag
      COMMON  / iend / iflag

      iflag = iiflag
      F = chisquare(X)

      RETURN
      END

      function chisquare(XX)
      implicit none
      integer nfit
      common /params/ nfit
      real*8 xx(nfit)
      real*8 chisquare

      real*8 CHI2_E866_800q, CHI2_E772_800
      REAL*8 CHI2_HERMES
      REAL*8 CHI2_PIP_HE_Z,CHI2_PIP_NE_Z,CHI2_PIP_KR_Z,CHI2_PIP_XE_Z
      REAL*8 CHI2_PIM_HE_Z,CHI2_PIM_NE_Z,CHI2_PIM_KR_Z,CHI2_PIM_XE_Z
      REAL*8 CHI2_PI0_HE_Z,CHI2_PI0_NE_Z,CHI2_PI0_KR_Z,CHI2_PI0_XE_Z
      REAL*8 CHI2_PIP_HE_PT2, CHI2_PIM_HE_PT2, CHI2_PI0_HE_PT2
      REAL*8 CHI2_PIP_NE_PT2, CHI2_PIM_NE_PT2, CHI2_PI0_NE_PT2
      REAL*8 CHI2_PIP_KR_PT2, CHI2_PIM_KR_PT2, CHI2_PI0_KR_PT2
      REAL*8 CHI2_PIP_XE_PT2, CHI2_PIM_XE_PT2, CHI2_PI0_XE_PT2
      REAL*8 CHI2_PIP_Z,CHI2_PI0_Z,CHI2_PIM_Z
      REAL*8 CHI2_PIP_PT2,CHI2_PI0_PT2,CHI2_PIM_PT2


C-----RHIC (DY)
      real*8 CHI2_RHIC_Ratio_pAu1, CHI2_RHIC_Ratio_pAu2

C-----LHC (DY)
      real*8 CHI2_ATLAS5_Y1,CHI2_ATLAS5_Y2,CHI2_ATLAS5_Y3
      real*8 CHI2_CMS5

C-----RHIC (DY)
      integer NUM_RHIC_Ratio_pAu1, NUM_RHIC_Ratio_pAu2

C-----LHC (DY)
      integer NUM_ATLAS5_Y1,NUM_ATLAS5_Y2,NUM_ATLAS5_Y3
      integer NUM_CMS5

      integer num_e772_800

C-----HERMES
      integer NUM_E866_800q
      integer NUM_HERMES
      INTEGER NUM_PIP_HE_Z,NUM_PIP_NE_Z,NUM_PIP_KR_Z,NUM_PIP_XE_Z
      INTEGER NUM_PIM_HE_Z,NUM_PIM_NE_Z,NUM_PIM_KR_Z,NUM_PIM_XE_Z
      INTEGER NUM_PI0_HE_Z,NUM_PI0_NE_Z,NUM_PI0_KR_Z,NUM_PI0_XE_Z
      INTEGER NUM_PIP_HE_PT2, NUM_PIM_HE_PT2, NUM_PI0_HE_PT2
      INTEGER NUM_PIP_NE_PT2, NUM_PIM_NE_PT2, NUM_PI0_NE_PT2
      INTEGER NUM_PIP_KR_PT2, NUM_PIM_KR_PT2, NUM_PI0_KR_PT2
      INTEGER NUM_PIP_XE_PT2, NUM_PIM_XE_PT2, NUM_PI0_XE_PT2
      INTEGER NUM_PIP_Z   ,NUM_PI0_Z   ,NUM_PIM_Z
      INTEGER NUM_PIP_PT2   ,NUM_PI0_PT2   ,NUM_PIM_PT2



C----- JLAB (12 GeV)
      integer NUM_JLAB12

      REAL*8 CMS5_NORM
      REAL*8 ATLAS5_Y1_NORM,ATLAS5_Y2_NORM,ATLAS5_Y3_NORM
      REAL*8 CX_STORE(100)

      integer IT,IH,IC
      common /meson/ IH,IC

      REAL*8 qTdQcut,phtcut
      REAL*8 pt2cut,zcut, pTdQdZcut

      real*8 norm,normsumnum,normsumden,ptmin,ptmax

      include "../fit/tools/data-inc.f"

      REAL*8 aN, bN, g2A, A2, g2B, B2
      COMMON /FITP/ aN, bN, g2A, A2, g2B, B2

      real*8 Sep_hermes
      real*8 Sep_JLAB
      data Sep_hermes/52.7d0/
      data Sep_JLAB/23.4d0/
      COMMON /NNUM_EXP/ NUM_EXP

      integer iflag,num_exp
      common /iend/ iflag

      !LHC variables
      real*8  RTS_LHC

      integer cutflag
      real*8 ptcut,etamin,etamax,cutfac
      common /cuts/ ptcut,etamin,etamax,cutflag

      !HERMES variables
      real*8 MULT, NU, Q2, PT2, STAT, SYS, VAL, Z
      !JLAB variables
      real*8 PTLOW, PTHIGH, nucl_target
      real*8 nulow, nuhigh
      real*8 Q2low, Q2high
      real*8 Qs


      real*8 CHI2
      real*8 xb,zh,pht,tmp
      real*8 sigma
      real*8 dcorr,fN
      real*8 rts,y,Q,Qmin,Qmax,pt,fuu,Rds
      real*8 fuuA,fuuB, DIS, R_a
      real*8 R_dy
      real*8 xfmin, xfmax, ymin, ymax, Qbar, const, ybar
      real*8 pi
      REAL*8 AN_REP, BN_REP
      COMMON /REP_PARAMS/ AN_REP, BN_REP
      data pi/3.1415926535d0/
      real*8 M
      data M/0.938d0/
      real*8 CX_LHC,CX_RHIC
      integer i

      aN    = AN_REP !
      bN    = BN_REP !
      g2A   = xx( 3)
      a2    = xx( 4)
      g2B   = xx( 5)
      b2    = xx( 6)

      print*, "aN", aN
      print*, "bN", bN
      print*, "gA", g2A
      print*, "gB", g2B


C-----CUTS TO DRELLYAN
      qTdQcut = 0.3d0
      phtcut = 1d0

C-----CUTS TO SIDIS
      pt2cut    = 0.30d0
      zcut      = 0.7d0
      pTdQdzcut = 1000

C-----RHIC (Ratios) Num of Data
      NUM_RHIC_Ratio_pAu1 = 0
      NUM_RHIC_Ratio_pAu2 = 0

C-----LHC NUM OF DATA
      NUM_ATLAS5_Y1 = 0
      NUM_ATLAS5_Y2 = 0
      NUM_ATLAS5_Y3 = 0

      NUM_CMS5 = 0

c-----FNAL Num of Data
      num_E772_800 = 0
      num_E866_800q = 0

C-----HERMES NUM OF DATA
      NUM_HERMES     = 0

      NUM_PIP_HE_Z   = 0
      NUM_PIP_NE_Z   = 0
      NUM_PIP_KR_Z   = 0
      NUM_PIP_XE_Z   = 0

      NUM_PIM_HE_Z   = 0
      NUM_PIM_NE_Z   = 0
      NUM_PIM_KR_Z   = 0
      NUM_PIM_XE_Z   = 0

      NUM_PI0_HE_Z   = 0
      NUM_PI0_NE_Z   = 0
      NUM_PI0_KR_Z   = 0
      NUM_PI0_XE_Z   = 0

      NUM_PIP_HE_PT2   = 0
      NUM_PIP_NE_PT2   = 0
      NUM_PIP_KR_PT2   = 0
      NUM_PIP_XE_PT2   = 0

      NUM_PIM_HE_PT2   = 0
      NUM_PIM_NE_PT2   = 0
      NUM_PIM_KR_PT2   = 0
      NUM_PIM_XE_PT2   = 0

      NUM_PI0_HE_PT2   = 0
      NUM_PI0_NE_PT2   = 0
      NUM_PI0_KR_PT2   = 0
      NUM_PI0_XE_PT2   = 0

      NUM_PIP_Z = 0
      NUM_PI0_Z = 0
      NUM_PIM_Z = 0

      NUM_PIP_PT2 = 0
      NUM_PI0_PT2 = 0
      NUM_PIM_PT2 = 0

      NUM_JLAB12 = 0


      CHI2 = 0d0

C-----RHIC CHI2
      CHI2_RHIC_Ratio_pAu1 = 0
      CHI2_RHIC_Ratio_pAu2 = 0

C-----CMS 5 TEV CHI2
      CHI2_CMS5  = 0d0

C-----ATLAS AT 5 TEV CHI2
      CHI2_ATLAS5_Y1  = 0d0
      CHI2_ATLAS5_Y2  = 0d0
      CHI2_ATLAS5_Y3  = 0d0

C-----DRELLYAN CHI2
      CHI2_E772_800 = 0d0
      CHI2_E866_800q = 0d0

C-----HERMES CHI2
      CHI2_HERMES     = 0d0

      CHI2_PIP_HE_Z   = 0d0
      CHI2_PIP_NE_Z   = 0d0
      CHI2_PIP_KR_Z   = 0d0
      CHI2_PIP_XE_Z   = 0d0

      CHI2_PIM_HE_Z   = 0d0
      CHI2_PIM_NE_Z   = 0d0
      CHI2_PIM_KR_Z   = 0d0
      CHI2_PIM_XE_Z   = 0d0

      CHI2_PI0_HE_Z   = 0d0
      CHI2_PI0_NE_Z   = 0d0
      CHI2_PI0_KR_Z   = 0d0
      CHI2_PI0_XE_Z   = 0d0

      CHI2_PIP_Z = 0d0
      CHI2_PI0_Z = 0d0
      CHI2_PIM_Z = 0d0


      CHI2_PIP_HE_PT2   = 0d0
      CHI2_PIP_NE_PT2   = 0d0
      CHI2_PIP_KR_PT2   = 0d0
      CHI2_PIP_XE_PT2   = 0d0

      CHI2_PIM_HE_PT2   = 0d0
      CHI2_PIM_NE_PT2   = 0d0
      CHI2_PIM_KR_PT2   = 0d0
      CHI2_PIM_XE_PT2   = 0d0

      CHI2_PI0_HE_PT2   = 0d0
      CHI2_PI0_NE_PT2   = 0d0
      CHI2_PI0_KR_PT2   = 0d0

      CHI2_PIP_PT2 = 0d0
      CHI2_PI0_PT2 = 0d0
      CHI2_PIM_PT2 = 0d0




      IF(I_EIC_PRED.EQ.1) THEN
        open(unit=  5,file='plot_data/EIC_Qlo.dat')
        WRITE(5,*) "PhT R"
        DO I = 1,10
          pht = 0.07*i
          rts = 2d0*2d0*10d0*110d0
          xb = 0.05d0
          z  = 0.4d0
          Qs = 4d0
          IH = 1
          IC = 1
          IT = 197
          CALL DISUU(rts,Qs,xb,z,pht,fuua,IT,IH,IC)
          IT = 3
          CALL DISUU(rts,Qs,xb,z,pht,fuu,IT,IH,IC)
          DIS = EIC_PRED_DIS(1)
          R_a = fuua/(fuu*DIS)
          WRITE(5,*) Pht,R_a
        ENDDO
      ENDIF


      IF(I_EIC_PRED.EQ.1) THEN
        open(unit=  5,file='plot_data/EIC_Qmed.dat')
        WRITE(5,*) "PhT R"
        DO I = 1,10
          pht = 0.07*i
          rts = 2d0*2d0*10d0*110d0
          xb = 0.05d0
          z  = 0.4d0
          Qs = 25d0
          IH = 1
          IC = 1
          IT = 197
          CALL DISUU(rts,Qs,xb,z,pht,fuua,IT,IH,IC)
          IT = 3
          CALL DISUU(rts,Qs,xb,z,pht,fuu,IT,IH,IC)
          DIS = EIC_PRED_DIS(2)
          R_a = fuua/(fuu*DIS)
          WRITE(5,*) Pht,R_a
        ENDDO
      ENDIF


      IF(I_EIC_PRED.EQ.1) THEN
        open(unit=  5,file='plot_data/EIC_Qhi.dat')
        WRITE(5,*) "PhT R"
        DO I = 1,10
          pht = 0.07*i
          rts = 2d0*2d0*10d0*110d0
          xb = 0.05d0
          z  = 0.4d0
          Qs = 100d0
          IH = 1
          IC = 1
          IT = 197
          CALL DISUU(rts,Qs,xb,z,pht,fuua,IT,IH,IC)
          IT = 3
          CALL DISUU(rts,Qs,xb,z,pht,fuu,IT,IH,IC)
          DIS = EIC_PRED_DIS(3)
          R_a = fuua/(fuu*DIS)
          WRITE(5,*) Pht,R_a
        ENDDO
      ENDIF


      IF(I_JLAB_PRED.EQ.1) THEN
        open(unit=  5,file='plot_data/JLAB_pred.dat')
        WRITE(5,*) "PhT R"
        DO I = 1,10
          pht = 0.07*i
          rts = Sep_JLAB
          xb = 0.4d0
          z  = 0.4d0
          Qs = 2.5d0
          IH = 1
          IC = 1
          IT = 208
          CALL DISUU(rts,Qs,xb,z,pht,fuua,IT,IH,IC)
          IT = 3
          CALL DISUU(rts,Qs,xb,z,pht,fuu,IT,IH,IC)
          DIS = JLAB_PRED_DIS
          R_a = fuua/(fuu*DIS)
          WRITE(5,*) Pht,R_a
        ENDDO
      ENDIF



C------------- JLAB (Miguel's data)
      if(I_JLAB_PIP_12.eq.1) then
      open(unit=  5,file='plot_data/JLAB_piplus.dat')
      WRITE(5,*) 'target  ', 'Z ' , 'PTMIN ', 'PTMAX ', 'MULT-RATIO ',
     *       'ERR ', 'fuu ', 'fuua ',  'DIS'
      do I=1,180
         rts = Sep_JLAB
         nucl_target = JLAB_PIP_TARGET(I)
         z = JLAB_PIP_Z(I)-0.05d0
         MULT = JLAB_PIP_MULT_RATIO(I)
         sigma = JLAB_PIP_ERR(I)
         Q2LOW  = 1d0
         Q2HIGH = 4.1d0
         NULOW  = 2.2d0
         NUHIGH = 4.2d0
         PTLOW =  JLAB_PIP_PTLOW(I)
         PTHIGH = JLAB_PIP_PTHIGH(I)
         Qs = 0.5d0*(Q2LOW+Q2HIGH)
         Nu = 0.5d0*(NULOW+NUHIGH)
         xb  = Qs/(2*M*Nu)
         pht = 0.5d0*(PTLOW+PTHIGH)
         IH = 1
         IC = 1
         IF(nucl_target.eq.0) THEN
         IT = 12
         ELSEIF(nucl_target.eq.1) THEN
         IT = 56
         ELSEIF(nucl_target.eq.2) THEN
         IT = 208
         ENDIF
         ! Computing SIDIS CX for e+A
         CALL DISUU(rts,Qs,xb,z,pht,fuua,IT,IH,IC)
         ! Computing SIDIS CX for d+A
         IT = 3
         CALL DISUU(rts,Qs,xb,z,pht,fuu,IT,IH,IC)
         DIS = JLAB_PIP_DIS(I)
         R_a = fuua/(fuu*DIS)
         tmp = (R_a - MULT)**2d0/sigma**2d0
         CHI2 = CHI2 + tmp
         CHI2_HERMES = CHI2_HERMES + tmp
         NUM_JLAB12 = NUM_JLAB12 + 1
         WRITE(5,*) nucl_target, z, ptlow, pthigh, MULT, sigma,
     *               fuu, fuua, DIS
C         if((nucl_target.eq.2) .AND. (z.lt.0.7) .AND. (pt.lt.0.7) )THEN
C            CALL DISUU(rts,Qs,xb,z,pht,fuua,IT,IH,IC)
C            ! Computing SIDIS CX for d+A
C            IT = 3
c            CALL DISUU(rts,Qs,xb,z,pht,fuu,IT,IH,IC)
C            DIS = JLAB_DIS(I)
c            R_a = fuua/(fuu*DIS)
c            tmp = (R_a - MULT)**2d0/sigma**2d0
c            CHI2 = CHI2 + tmp
c            CHI2_HERMES = CHI2_HERMES + tmp
c            NUM_JLAB12 = NUM_JLAB12 + 1
C            WRITE(5,*) nucl_target, z, ptlow, pthigh, MULT, sigma,
c     *               fuu, fuua, DIS
C         ENDIF
      enddo
      close(5)
      endif

C
C------------- JLAB (Miguel's data)
      if(I_JLAB_PIM_12.eq.1) then
      open(unit=  5,file='plot_data/JLAB_piminus.dat')
      WRITE(5,*) 'target  ', 'Z ' , 'PTMIN ', 'PTMAX ', 'MULT-RATIO ',
     *       'ERR ', 'fuu ', 'fuua ',  'DIS'
      do I=1,180
         rts = Sep_JLAB
         nucl_target = JLAB_PIM_TARGET(I)
         z = JLAB_PIM_Z(I)-0.05d0
         MULT = JLAB_PIM_MULT_RATIO(I)
         sigma = JLAB_PIM_ERR(I)
         Q2LOW  = 1d0
         Q2HIGH = 4.1d0
         NULOW  = 2.2d0
         NUHIGH = 4.2d0
         PTLOW =  JLAB_PIM_PTLOW(I)
         PTHIGH = JLAB_PIM_PTHIGH(I)
         Qs = 0.5d0*(Q2LOW+Q2HIGH)
         Nu = 0.5d0*(NULOW+NUHIGH)
         xb  = Qs/(2*M*Nu)
         pht = 0.5d0*(PTLOW+PTHIGH)
         IH = 1
         IC = -1
         IF(nucl_target.eq.0) THEN
         IT = 12
         ELSEIF(nucl_target.eq.1) THEN
         IT = 56
         ELSEIF(nucl_target.eq.2) THEN
         IT = 208
         ENDIF
         ! Computing SIDIS CX for e+A
         CALL DISUU(rts,Qs,xb,z,pht,fuua,IT,IH,IC)
         ! Computing SIDIS CX for d+A
         IT = 3
         CALL DISUU(rts,Qs,xb,z,pht,fuu,IT,IH,IC)
         DIS = JLAB_PIM_DIS(I)
         R_a = fuua/(fuu*DIS)
         tmp = (R_a - MULT)**2d0/sigma**2d0
         CHI2 = CHI2 + tmp
         CHI2_HERMES = CHI2_HERMES + tmp
         NUM_JLAB12 = NUM_JLAB12 + 1
         WRITE(5,*) nucl_target, z, ptlow, pthigh, MULT, sigma,
     *               fuu, fuua, DIS
C         if((nucl_target.eq.2) .AND. (z.lt.0.7) .AND. (pt.lt.0.7) )THEN
C            CALL DISUU(rts,Qs,xb,z,pht,fuua,IT,IH,IC)
C            ! Computing SIDIS CX for d+A
C            IT = 3
c            CALL DISUU(rts,Qs,xb,z,pht,fuu,IT,IH,IC)
C            DIS = JLAB_DIS(I)
c            R_a = fuua/(fuu*DIS)
c            tmp = (R_a - MULT)**2d0/sigma**2d0
c            CHI2 = CHI2 + tmp
c            CHI2_HERMES = CHI2_HERMES + tmp
c            NUM_JLAB12 = NUM_JLAB12 + 1
C            WRITE(5,*) nucl_target, z, ptlow, pthigh, MULT, sigma,
c     *               fuu, fuua, DIS
C         ENDIF
      enddo
      close(5)
      endif

C-------------DY (RHIC) p + Au -> muons (Au-going)
      if(I_RHIC_Ratio_pAu1.eq.1) then
      open(unit=  5,file='plot_data/RHIC_Ratio_pAu1.dat')
      write(5 ,*) 'pt ','CX ', 'ERR ', 'FUU'
      do I=1,5
         rts       = 200d0
         pt        = RHIC_Ratio_pAu1_pT(I)
         CX_RHIC   = RHIC_Ratio_pAu1_CX(I)
         sigma     = RHIC_Ratio_pAu1_ERR(I)
         fuuB      = RHIC_pp_FUU(I)
         Qmin      = 4.8d0
         Qmax      = 8.2d0
         Qbar     = (Qmin+Qmax)/2d0
         ymin     = -2.2d0
         ymax     = -1.2d0
         if (pt/Qbar.lt.qTdQcut) then
           IT = 197
           call DY_overyu(1,rts,pt,ymin,ymax,Qmin,Qmax,fuuA,IT)
           R_dy = fuuA/fuuB
           write(5,*) pt, CX_RHIC, sigma, R_dy
           tmp = (R_dy - CX_RHIC)**2d0/sigma**2d0
           CHI2_RHIC_Ratio_pAu1 = CHI2_RHIC_Ratio_pAu1 + tmp
           CHI2 = CHI2 + tmp
           num_RHIC_Ratio_pAu1 = num_RHIC_Ratio_pAu1 + 1
        endif
      enddo
      close(5)
      endif

C-------------DY (RHIC) p + Au -> muons (p-going)
      if(I_RHIC_Ratio_pAu2.eq.1) then
      open(unit=  5,file='plot_data/RHIC_Ratio_pAu2.dat')
      write(5 ,*) 'pt ','CX ', 'ERR ', 'FUU'
      do I=1,5
         rts       = 200d0
         pt        = RHIC_Ratio_pAu2_pT(I)
         CX_RHIC   = RHIC_Ratio_pAu2_CX(I)
         sigma     = RHIC_Ratio_pAu2_ERR(I)
         fuuB      = RHIC_pp_FUU(I)
         Qmin      = 4.8d0 ! GeV
         Qmax      = 8.2d0 ! GeV
         Qbar     = (Qmin+Qmax)/2d0
         ymin     = 1.2d0
         ymax     = 2.2d0
         if (pt/Qbar.lt.qTdQcut) then
           ! Compute X section for p+Au (A = 195 OR 197 available)
           IT = 197
           call DY_overyu(1,rts,pt,ymin,ymax,Qmin,Qmax,fuuA,IT)
           R_dy = fuuA/fuuB
           write(5,*) pt, CX_RHIC, sigma, R_dy
           tmp = (R_dy - CX_RHIC)**2d0/sigma**2d0
           CHI2_RHIC_Ratio_pAu2 = CHI2_RHIC_Ratio_pAu2 + tmp
           CHI2 = CHI2 + tmp
           num_RHIC_Ratio_pAu2 = num_RHIC_Ratio_pAu2 + 1
        endif
      enddo
      close(5)
      endif

      normsumnum = 0d0
      normsumden = 0d0
      cutflag = 1
      etamin = -2.4d0
      etamax =  2.4d0
      ptcut = 20d0
      if(I_CMS5.eq.1) then
      open(unit=  5,file='plot_data/CMS5.dat')
      write(5 ,*) 'pt ','CX ', 'ERR ', 'FUU'
      do I=1,13
         rts      = 5020D0
         pt       = CMS5_PT(I)
         CX_LHC   = CMS5_CX(I)*1d3
         sigma    = CMS5_ER(I)*1d3
         ymin     = -2.8D0
         ymax     =  2D0
         QBAR = 91.2D0
         dcorr = 0.035
         IT = 208
         if (pt/Qbar.lt.qTdQcut) then
         call DY_overyu(1,rts,pt,ymin,ymax,60d0,120d0,fuuA,IT)
         FUUA = 2d0*pi*pt*FUUA*208
         CX_STORE(I) = FUUA
         normsumnum = normsumnum+CX_LHC*CX_LHC/sigma/sigma
         normsumden = normsumden+CX_LHC*FUUA  /sigma/sigma
         endif
      enddo
      normsumnum = normsumnum+1./dcorr/dcorr
      normsumden = normsumden+1./dcorr/dcorr
      close(5)
      endif
      norm =  normsumnum/normsumden
      fN = norm
      CMS5_NORM = norm

      print *, "CMS5_NORM", CMS5_NORM

C----------------
      if(I_CMS5.eq.1) then
      open(unit=  5,file='plot_data/CMS5.dat')
      write(5 ,*) 'pt ','CX ', 'ERR ', 'FUU'
      do I=1,13
         rts      = 5020D0
         pt       = CMS5_PT(I)
         CX_LHC   = CMS5_CX(I)*1d3
         sigma    = CMS5_ER(I)*1d3
         ymin     = -2.8D0
         ymax     =  2D0
         QBAR = 91.2D0
         IT = 208
         if (pt/Qbar.lt.qTdQcut) then
         fuua = norm*CX_STORE(I)
         tmp = (fuuA - CX_LHC)**2d0/sigma**2d0
         CHI2_CMS5 = CHI2_CMS5 + tmp
         CHI2 = CHI2 + tmp
         num_CMS5 = num_CMS5 + 1
         write(5,*) pt, CX_LHC, sigma, fuuA
         endif
      enddo
      CHI2_CMS5      = CHI2_CMS5+(1.-fN)**2./(dcorr**2.)
      chi2           = chi2     +(1.-fN)**2./(dcorr**2.)
      close(5)
      endif
      cutflag = 0

      normsumnum = 0d0
      normsumden = 0d0
      if(I_ATLAS5_Y1.eq.1) then
      open(unit=  5,file='plot_data/ATLAS5_Y1.dat')
      write(5 ,*) 'pt ','CX ', 'ERR ', 'FUU'
      do I=1,14
         rts      = 5020D0
         pt       = ATLAS5_Y1_PT(I)
         CX_LHC   = ATLAS5_Y1_CX(I)*1d3
         sigma    = ATLAS5_Y1_ER(I)*1d3
         dcorr    = 0.027
         ymin     = -3D0
         ymax     =  2D0
         QBAR = 91.2D0
         IT = 208
         if (pt/Qbar.lt.qTdQcut) then
         call DY_overyu(1,rts,pt,ymin,ymax,66d0,116d0,fuuA,IT)
         FUUA = 2d0*pi*pt*FUUA/PT*208
         CX_STORE(I) = FUUA
         normsumnum = normsumnum+CX_LHC*CX_LHC/sigma/sigma
         normsumden = normsumden+CX_LHC*FUUA  /sigma/sigma
         endif
      enddo
      normsumnum = normsumnum+1./dcorr/dcorr
      normsumden = normsumden+1./dcorr/dcorr
      close(5)
      endif
      norm =  normsumnum/normsumden
      fN = norm
      ATLAS5_Y1_NORM = norm

      print *, "ATLAS5_Y1_NORM", ATLAS5_Y1_NORM

      if(I_ATLAS5_Y1.eq.1) then
      open(unit=  5,file='plot_data/ATLAS5_Y1.dat')
      write(5 ,*) 'pt ','CX ', 'ERR ', 'FUU'
      do I=1,14
         rts      = 5020D0
         pt       = ATLAS5_Y1_PT(I)
         CX_LHC   = ATLAS5_Y1_CX(I)*1d3
         sigma    = ATLAS5_Y1_ER(I)*1d3
         dcorr    = 0.027
         ymin     = -3D0
         ymax     =  2D0
         QBAR = 91.2D0
         IT = 208
         if (pt/Qbar.lt.qTdQcut) then
         fuua = norm*CX_STORE(I)
         tmp = (fuuA - CX_LHC)**2d0/sigma**2d0
         CHI2_ATLAS5_Y1 = CHI2_ATLAS5_Y1 + tmp
         CHI2 = CHI2 + tmp
         num_ATLAS5_Y1 = num_ATLAS5_Y1 + 1
         write(5,*) pt, CX_LHC, sigma, fuuA
         endif
      enddo
      close(5)
      chi2_atlas5_Y1 = CHI2_ATLAS5_Y1+(1.-fN)**2./(dcorr**2.)
      chi2           = chi2          +(1.-fN)**2./(dcorr**2.)
      endif



      normsumnum = 0d0
      normsumden = 0d0
      if(I_ATLAS5_Y2.eq.1) then
      open(unit=  5,file='plot_data/ATLAS5_Y2.dat')
      write(5 ,*) 'pt ','CX ', 'ERR ', 'FUU'
      do I=1,8
         rts      = 5020D0
         pt       = ATLAS5_Y2_PT(I)
         CX_LHC   = ATLAS5_Y2_CX(I)*1d3
         sigma    = ATLAS5_Y2_ER(I)*1d3
         dcorr    = 0.027
         ymin     = -2D0
         ymax     =  0D0
         QBAR = 91.2D0
         IT = 208
         if (pt/Qbar.lt.qTdQcut) then
         call DY_overyu(1,rts,pt,ymin,ymax,66d0,116d0,fuuA,IT)
         FUUA = 2d0*pi*pt*FUUA/PT*208
         CX_STORE(I) = FUUA
         normsumnum = normsumnum+CX_LHC*CX_LHC/sigma/sigma
         normsumden = normsumden+CX_LHC*FUUA  /sigma/sigma
         endif
      enddo
      normsumnum = normsumnum+1./dcorr/dcorr
      normsumden = normsumden+1./dcorr/dcorr
      close(5)
      endif
      norm =  normsumnum/normsumden
      fN = norm
      ATLAS5_Y2_NORM = norm

      print *, "ATLAS5_Y2_NORM", ATLAS5_Y2_NORM


      if(I_ATLAS5_Y2.eq.1) then
      open(unit=  5,file='plot_data/ATLAS5_Y2.dat')
      write(5 ,*) 'pt ','CX ', 'ERR ', 'FUU'
      do I=1,8
         rts      = 5020D0
         pt       = ATLAS5_y2_PT(I)
         CX_LHC   = ATLAS5_y2_CX(I)*1d3
         sigma    = ATLAS5_y2_ER(I)*1d3
         ymin     = -2D0
         ymax     =  0D0
         QBAR = 91.2D0
         IT = 208
         if (pt/Qbar.lt.qTdQcut) then
         fuua =	norm*CX_STORE(I)
         tmp = (fuuA - CX_LHC)**2d0/sigma**2d0
         CHI2_ATLAS5_Y2 = CHI2_ATLAS5_Y2 + tmp
         CHI2 = CHI2 + tmp
         num_ATLAS5_Y2 = num_ATLAS5_Y2 + 1
         write(5,*) pt, CX_LHC, sigma, fuuA
         endif
      enddo
      close(5)
      chi2_atlas5_Y2 = CHI2_ATLAS5_Y2+(1.-fN)**2./(dcorr**2.)
      chi2           = chi2          +(1.-fN)**2./(dcorr**2.)
      endif

      normsumnum = 0d0
      normsumden = 0d0
      if(I_ATLAS5_Y3.eq.1) then
      open(unit=  5,file='plot_data/ATLAS5_Y3.dat')
      write(5 ,*) 'pt ','CX ', 'ERR ', 'FUU'
      do I=1,8
         rts      = 5020D0
         pt       = ATLAS5_Y3_PT(I)
         CX_LHC   = ATLAS5_Y3_CX(I)*1d3
         sigma    = ATLAS5_Y3_ER(I)*1d3
         dcorr    = 0.027
         ymin     =  0D0
         ymax     =  2D0
         QBAR = 91.2D0
         IT = 208
         if (pt/Qbar.lt.qTdQcut) then
         call DY_overyu(1,rts,pt,ymin,ymax,66d0,116d0,fuuA,IT)
         FUUA = 2d0*pi*pt*FUUA/PT*208
         CX_STORE(I) = FUUA
         normsumnum = normsumnum+CX_LHC*CX_LHC/sigma/sigma
         normsumden = normsumden+CX_LHC*FUUA  /sigma/sigma
         endif
      enddo
      normsumnum = normsumnum+1./dcorr/dcorr
      normsumden = normsumden+1./dcorr/dcorr
      close(5)
      endif
      norm =  normsumnum/normsumden
      fN = norm
      ATLAS5_Y3_NORM = norm

      print *, "ATLAS5_Y3_NORM", ATLAS5_Y3_NORM


      if(I_ATLAS5_y3.eq.1) then
      open(unit=  5,file='plot_data/ATLAS5_Y3.dat')
      write(5 ,*) 'pt ','CX ', 'ERR ', 'FUU'
      do I=1,8
         rts      = 5020D0
         pt       = ATLAS5_y3_PT(I)
         CX_LHC   = ATLAS5_y3_CX(I)*1d3
         sigma    = ATLAS5_y3_ER(I)*1d3
         ymin     =  0D0
         ymax     =  2D0
         QBAR = 91.2D0
         IT = 208
         if (pt/Qbar.lt.qTdQcut) then
         fuua =	norm*CX_STORE(I)
         tmp = (fuuA - CX_LHC)**2d0/sigma**2d0
         CHI2_ATLAS5_Y3 = CHI2_ATLAS5_Y3 + tmp
         CHI2 = CHI2 + tmp
         num_ATLAS5_Y3 = num_ATLAS5_Y3 + 1
         write(5,*) pt, CX_LHC, sigma, fuuA
         endif
      enddo
      close(5)
      chi2_atlas5_Y3 = CHI2_ATLAS5_Y3+(1.-fN)**2./(dcorr**2.)
      chi2           = chi2          +(1.-fN)**2./(dcorr**2.)
      endif

C------E772
C------- C/D
      if(I_E772_800.eq.1) then
      open(unit=  5,file='plot_data/E772_800_CD.dat')
      write(5 ,*) 'pt ','DY-RATIO ', 'error ', 'R_dy '
      do I=1,7
         rts   = dsqrt(1600d0*0.938d0)
         pt    = E772_800_C_pT(I)
         Rds   = E772_800_C_R(I)
         sigma = E772_800_C_Err(I)
         Qmin = 4d0
         Qmax = 9d0
         Qbar = (Qmin+Qmax)/2d0
         xfmin = 0.05d0
         xfmax = 0.3d0
         const = rts/(Qbar)
         ymin = -dlog(const*xfmin)
         ymax = -dlog(const*xfmax)
         ybar = (ymin+ymax)/2d0
         if (pt/Qbar.lt.qTdQcut) then
           ! Compute X section for C (A = 12)
           !call DY_overQ(1,rts,xf,ybar,pt,Qmin,Qmax,fuuA,12)
           call DY_overxbQ(4,rts,pt,xfmin,xfmax,Qmin,Qmax,fuuA,12)
           ! Compute X section for Deuteron
           !call DY_overQ(1,rts,xf,ybar,pt,Qmin,Qmax,fuuB,3)
           call DY_overxbQ(4,rts,pt,xfmin,xfmax,Qmin,Qmax,fuuB,3)
           !print *, fuuA,fuuB
           R_dy = fuuA/fuuB
           write(5,*) pt, Rds, sigma, R_dy
           tmp = (R_dy - Rds)**2d0/sigma**2d0
           CHI2_E772_800 = CHI2_E772_800 + tmp
           CHI2 = CHI2 + tmp
           num_E772_800 = num_E772_800 + 1
        endif
      enddo
      close(5)
      endif

C--------- Ca/D
      if(I_E772_800.eq.1) then
      open(unit=  66,file='plot_data/E772_800_CaD.dat')
      write(66 ,*) 'pt ','DY-RATIO ', 'error ', 'R_dy '
      do I=1,7
         rts   = dsqrt(1600d0*0.938d0)
         pt    = E772_800_Ca_pT(I)
         Rds   = E772_800_Ca_R(I)
         sigma = E772_800_Ca_Err(I)
         Qmin = 4d0
         Qmax = 9d0
         Qbar = (Qmin+Qmax)/2d0
         xfmin = 0.05d0
         xfmax = 0.3d0
         const = rts/(Qbar)
         ymin = -dlog(const*xfmin)
         ymax = -dlog(const*xfmax)
         ybar = (ymin+ymax)/2d0
         if (pt/Qbar.lt.qTdQcut) then
           ! Compute X section for Ca (A = 40)
           !call DY_overQ(1,rts,xf,ybar,pt,Qmin,Qmax,fuuA,40)
           call DY_overxbQ(4,rts,pt,xfmin,xfmax,Qmin,Qmax,fuuA,40)
           ! Compute X section for Deuteron
           !call DY_overQ(1,rts,xf,ybar,pt,Qmin,Qmax,fuuB,3)
           call DY_overxbQ(4,rts,pt,xfmin,xfmax,Qmin,Qmax,fuuB,3)
           !print *, fuuA,fuuB
           R_dy = fuuA/fuuB
           write(66,*) pt, Rds, sigma, R_dy
           tmp = (R_dy - Rds)**2d0/sigma**2d0
           CHI2_E772_800 = CHI2_E772_800 + tmp
           CHI2 = CHI2 + tmp
           num_E772_800 = num_E772_800 + 1
         endif
      enddo
      close(66)
      endif

C--------- Fe/D
      if(I_E772_800.eq.1) then
      open(unit=  7,file='plot_data/E772_800_FeD.dat')
      write(7 ,*) 'pt ','DY-RATIO ', 'error ', 'R_dy '
      do I=1,7
         rts   = dsqrt(1600d0*0.938d0)
         pt    = E772_800_Fe_pT(I)
         Rds   = E772_800_Fe_R(I)
         sigma = E772_800_Fe_Err(I)
         Qmin = 4d0
         Qmax = 9d0
         Qbar = (Qmin+Qmax)/2d0
         xfmin = 0.05d0
         xfmax = 0.3d0
         const = rts/(Qbar)
         ymin = -dlog(const*xfmin)
         ymax = -dlog(const*xfmax)
         ybar = (ymin+ymax)/2d0
         if (pt/Qbar.lt.qTdQcut) then
           ! Compute X section for Fe (A = 56)
           !call DY_overQ(1,rts,xf,ybar,pt,Qmin,Qmax,fuuA,56)
           call DY_overxbQ(4,rts,pt,xfmin,xfmax,Qmin,Qmax,fuuA,56)
           ! Compute X section for Deuteron
           !call DY_overQ(1,rts,xf,ybar,pt,Qmin,Qmax,fuuB,3)
           call DY_overxbQ(4,rts,pt,xfmin,xfmax,Qmin,Qmax,fuuB,3)
           !print *, fuuA,fuuB
           R_dy = fuuA/fuuB
           write(7,*) pt, Rds, sigma, R_dy
           tmp = (R_dy - Rds)**2d0/sigma**2d0
           CHI2_E772_800 = CHI2_E772_800 + tmp
           CHI2 = CHI2 + tmp
           num_E772_800 = num_E772_800 + 1
         endif
      enddo
      close(7)
      endif

C--------- W/D
      if(I_E772_800.eq.1) then
      open(unit=  8,file='plot_data/E772_800_WD.dat')
      write(8 ,*) 'pt ','DY-RATIO ', 'error ', 'R_dy '
      do I=1,7
         rts   = dsqrt(1600d0*0.938d0)
         pt    = E772_800_W_pT(I)
         Rds   = E772_800_W_R(I)
         sigma = E772_800_W_Err(I)
         Qmin = 4d0
         Qmax = 9d0
         Qbar = (Qmin+Qmax)/2d0
         xfmin = 0.05d0
         xfmax = 0.3d0
         const = rts/(Qbar)
         ymin = -dlog(const*xfmin)
         ymax = -dlog(const*xfmax)
         ybar = (ymin+ymax)/2d0
         if (pt/Qbar.lt.qTdQcut) then
           ! Compute X section for W (A = 184)
           !call DY_overQ(1,rts,xf,ybar,pt,Qmin,Qmax,fuuA,184)
           call DY_overxbQ(4,rts,pt,xfmin,xfmax,Qmin,Qmax,fuuA,184)
           ! Compute X section for Deuteron
           !call DY_overQ(1,rts,xf,ybar,pt,Qmin,Qmax,fuuB,3)
           call DY_overxbQ(4,rts,pt,xfmin,xfmax,Qmin,Qmax,fuuB,3)
           !print *, fuuA,fuuB
           R_dy = fuuA/fuuB
           write(8,*) pt, Rds, sigma, R_dy
           tmp = (R_dy - Rds)**2d0/sigma**2d0
           CHI2_E772_800 = CHI2_E772_800 + tmp
           CHI2 = CHI2 + tmp
           num_E772_800 = num_E772_800 + 1
         endif
      enddo
      close(8)
      endif

C-----E866 (Q binned data)
      if(I_E866_800q.eq.1) then
      open(unit=  40,file='plot_data/E866_800_Qbin_FeBe.dat')
      write(40 ,*) 'pt ','DY-RATIO ','Q ', 'error ', 'R_dy '
      do I=1,32
         rts   = dsqrt(1600d0*0.938d0)
         pt    = E866_800q_pT(I)
         Rds   = E866_800q_R_FeBe(I)
         Q     = E866_800q_Q(I)
         sigma = E866_800q_Err_FeBe(I)
         xfmin = 0.13d0
         xfmax = 0.93d0
         ymin = asinh(rts/2d0/Q*xfmin)
         ymax = asinh(rts/2d0/Q*xfmax)
         if (pt/Q.lt.qTdQcut) then
           !Compute X section for Be (A = 9)
           call DY_overxF(3,rts,pt,xfmin,xfmax,Q,fuuB,9)
           !Compute X section for Fe (A = 56)
           call DY_overxF(3,rts,pt,xfmin,xfmax,Q,fuuA,56)
           R_dy = fuuA/fuuB
           write(40,*) pt, Rds, Q, sigma, R_dy
           tmp = (R_dy - Rds)**2d0/sigma**2d0
           CHI2_E866_800q = CHI2_E866_800q + tmp
           CHI2 = CHI2 + tmp
           num_E866_800q = num_E866_800q + 1
         endif
      enddo
      close(40)
      endif

C------W/Be
      if(I_E866_800q.eq.1) then
      open(unit=  50,file='plot_data/E866_800_Qbin_WBe.dat')
      write(50 ,*) 'pt ','DY-RATIO ','x2 ', 'error ', 'R_dy '
      do I=1,32
         rts   = dsqrt(1600d0*0.938d0)
         pt    = E866_800q_pT(I)
         Rds   = E866_800q_R_WBe(I)
         Q     = E866_800q_Q(I)
         sigma = E866_800q_Err_WBe(I)
         xfmin = 0.13d0
         xfmax = 0.93d0
         ymin = asinh(rts/2d0/Q*xfmin)
         ymax = asinh(rts/2d0/Q*xfmax)
         if (pt/Q.lt.qTdQcut) then
           !Compute X section for Be (A = 9)
           call DY_overxF(3,rts,pt,xfmin,xfmax,Q,fuuB,9)
           !Compute X section for W (A = 184)
           call DY_overxF(3,rts,pt,xfmin,xfmax,Q,fuuA,184)
           R_dy = fuuA/fuuB
           write(50,*) pt, Rds, Q, sigma , R_dy
           tmp = (R_dy - Rds)**2d0/sigma**2d0
           CHI2_E866_800q = CHI2_E866_800q + tmp
           CHI2 = CHI2 + tmp
           num_E866_800q = num_E866_800q + 1
         endif
      enddo
      close(50)
      endif

      IF(I_PIP_HE_Z.EQ.1) THEN
      OPEN(UNIT = 5,FILE ='plot_data/HERMES/pip_he_z.dat')
      WRITE(5,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
      DO I=1,9
          VAL  = PIP_HE_Z_VALUE(I)
          MULT = PIP_HE_Z_MULT(I)
          STAT = PIP_HE_Z_STAT(I)
          SYS  = PIP_HE_Z_SYS(I)
          Q2   = PIP_HE_Z_Q2(I)
          Nu   = PIP_HE_Z_NU(I)
          Z    = VAL
          PT2  = PIP_HE_Z_PT2(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          rts = Sep_hermes
          IT = 4
          IH = 1
          IC = 1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(rts,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(rts,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PIP_HE_Z_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PIP_HE_Z = CHI2_PIP_HE_Z + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PIP_HE_Z = NUM_PIP_HE_Z + 1
            WRITE(5,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(5)
      ENDIF

      IF(I_PIP_NE_Z.EQ.1) THEN
      OPEN(UNIT = 66,FILE ='plot_data/HERMES/pip_ne_z.dat')
      WRITE(66,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
      DO I=1,9
          VAL  = PIP_NE_Z_VALUE(I)
          MULT = PIP_NE_Z_MULT(I)
          STAT = PIP_NE_Z_STAT(I)
          SYS  = PIP_NE_Z_SYS(I)
          Q2   = PIP_NE_Z_Q2(I)
          Nu   = PIP_NE_Z_NU(I)
          Z    = VAL
          PT2  = PIP_NE_Z_PT2(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          rts = Sep_hermes
          IT = 20
          IH = 1
          IC = 1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(rts,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(rts,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PIP_NE_Z_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PIP_NE_Z = CHI2_PIP_NE_Z + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PIP_NE_Z = NUM_PIP_NE_Z + 1
            WRITE(66,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(66)
      ENDIF

      IF(I_PIP_KR_Z.EQ.1) THEN
      OPEN(UNIT = 7,FILE ='plot_data/HERMES/pip_kr_z.dat')
      WRITE(7,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
      DO I=1,9
          VAL  = PIP_KR_Z_VALUE(I)
          MULT = PIP_KR_Z_MULT(I)
          STAT = PIP_KR_Z_STAT(I)
          SYS  = PIP_KR_Z_SYS(I)
          Q2   = PIP_KR_Z_Q2(I)
          Nu   = PIP_KR_Z_NU(I)
          Z    = VAL
          PT2  = PIP_KR_Z_PT2(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          rts = Sep_hermes
          IT = 84
          IH = 1
          IC = 1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(rts,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(rts,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PIP_KR_Z_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PIP_KR_Z = CHI2_PIP_KR_Z + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PIP_KR_Z = NUM_PIP_KR_Z + 1
            WRITE(7,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(7)
      ENDIF

      IF(I_PIP_XE_Z.EQ.1) THEN
      OPEN(UNIT = 8,FILE ='plot_data/HERMES/pip_xe_z.dat')
      WRITE(8,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
      DO I=1,9
          VAL  = PIP_XE_Z_VALUE(I)
          MULT = PIP_XE_Z_MULT(I)
          STAT = PIP_XE_Z_STAT(I)
          SYS  = PIP_XE_Z_SYS(I)
          Q2   = PIP_XE_Z_Q2(I)
          Nu   = PIP_XE_Z_NU(I)
          Z    = VAL
          PT2  = PIP_XE_Z_PT2(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          rts = Sep_hermes
          IT = 131
          IH = 1
          IC = 1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(rts,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(rts,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PIP_XE_Z_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PIP_XE_Z = CHI2_PIP_XE_Z + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PIP_XE_Z = NUM_PIP_XE_Z + 1
            WRITE(8,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(8)
      ENDIF

      IF(I_PIM_HE_Z.EQ.1) THEN
      OPEN(UNIT = 41,FILE ='plot_data/HERMES/pim_he_z.dat')
      WRITE(41,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
      DO I=1,9
          VAL  = PIM_HE_Z_VALUE(I)
          MULT = PIM_HE_Z_MULT(I)
          STAT = PIM_HE_Z_STAT(I)
          SYS  = PIM_HE_Z_SYS(I)
          Q2   = PIM_HE_Z_Q2(I)
          Nu   = PIM_HE_Z_NU(I)
          Z    = VAL
          PT2  = PIM_HE_Z_PT2(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          rts = Sep_hermes
          IT = 4
          IH = 1
          IC = -1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(rts,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(rts,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PIM_HE_Z_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PIM_HE_Z = CHI2_PIM_HE_Z + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PIM_HE_Z = NUM_PIM_HE_Z + 1
            WRITE(41,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(41)
      ENDIF

      IF(I_PIM_NE_Z.EQ.1) THEN
      OPEN(UNIT = 42,FILE ='plot_data/HERMES/pim_ne_z.dat')
      WRITE(42,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
      DO I=1,9
          VAL  = PIM_NE_Z_VALUE(I)
          MULT = PIM_NE_Z_MULT(I)
          STAT = PIM_NE_Z_STAT(I)
          SYS  = PIM_NE_Z_SYS(I)
          Q2   = PIM_NE_Z_Q2(I)
          Nu   = PIM_NE_Z_NU(I)
          Z    = VAL
          PT2  = PIM_NE_Z_PT2(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          rts = Sep_hermes
          IT = 20
          IH = 1
          IC = -1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(rts,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(rts,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PIM_NE_Z_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PIM_NE_Z = CHI2_PIM_NE_Z + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PIM_NE_Z = NUM_PIM_NE_Z + 1
            WRITE(42,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(42)
      ENDIF

      IF(I_PIM_KR_Z.EQ.1) THEN
      OPEN(UNIT = 43,FILE ='plot_data/HERMES/pim_kr_z.dat')
      WRITE(43,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
      DO I=1,9
          VAL  = PIM_KR_Z_VALUE(I)
          MULT = PIM_KR_Z_MULT(I)
          STAT = PIM_KR_Z_STAT(I)
          SYS  = PIM_KR_Z_SYS(I)
          Q2   = PIM_KR_Z_Q2(I)
          Nu   = PIM_KR_Z_NU(I)
          Z    = VAL
          PT2  = PIM_KR_Z_PT2(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          rts = Sep_hermes
          IT = 84
          IH = 1
          IC = -1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(rts,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(rts,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PIM_KR_Z_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PIM_KR_Z = CHI2_PIM_KR_Z + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PIM_KR_Z = NUM_PIM_KR_Z + 1
            WRITE(43,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(43)
      ENDIF

      IF(I_PIM_XE_Z.EQ.1) THEN
      OPEN(UNIT = 44,FILE ='plot_data/HERMES/pim_xe_z.dat')
      WRITE(44,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
      DO I=1,9
          VAL  = PIM_XE_Z_VALUE(I)
          MULT = PIM_XE_Z_MULT(I)
          STAT = PIM_XE_Z_STAT(I)
          SYS  = PIM_XE_Z_SYS(I)
          Q2   = PIM_XE_Z_Q2(I)
          Nu   = PIM_XE_Z_NU(I)
          Z    = VAL
          PT2  = PIM_XE_Z_PT2(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          rts = Sep_hermes
          IT = 131
          IH = 1
          IC = -1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(rts,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(rts,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PIM_XE_Z_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PIM_XE_Z = CHI2_PIM_XE_Z + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PIM_XE_Z = NUM_PIM_XE_Z + 1
            WRITE(44,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(44)
      ENDIF

C------------------
C------ PT2 DEPENDET DATA
C------- PI + DATA
      IF(I_PIP_HE_PT2.EQ.1) THEN
      OPEN(UNIT = 73,FILE ='plot_data/HERMES/pip_he_pt2.dat')
      WRITE(73,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
      DO I=1,9
          VAL  = PIP_HE_PT2_VALUE(I)
          MULT = PIP_HE_PT2_MULT(I)
          STAT = PIP_HE_PT2_STAT(I)
          SYS  = PIP_HE_PT2_SYS(I)
          Q2   = PIP_HE_PT2_Q2(I)
          Nu   = PIP_HE_PT2_NU(I)
          PT2  = VAL
          Z    = PIP_HE_PT2_Z(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          rts = Sep_hermes
          IT = 4
          IH = 1
          IC = 1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(rts,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(rts,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PIP_HE_PT2_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PIP_HE_PT2 = CHI2_PIP_HE_PT2 + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PIP_HE_PT2 = NUM_PIP_HE_PT2 + 1
            WRITE(73,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(73)
      ENDIF

      IF(I_PIP_NE_PT2.EQ.1) THEN
      OPEN(UNIT = 74,FILE ='plot_data/HERMES/pip_ne_pt2.dat')
      WRITE(74,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
      DO I=1,9
          VAL  = PIP_NE_PT2_VALUE(I)
          MULT = PIP_NE_PT2_MULT(I)
          STAT = PIP_NE_PT2_STAT(I)
          SYS  = PIP_NE_PT2_SYS(I)
          Q2   = PIP_NE_PT2_Q2(I)
          Nu   = PIP_NE_PT2_NU(I)
          PT2  = VAL
          Z    = PIP_NE_PT2_Z(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          rts = Sep_hermes
          IT = 20
          IH = 1
          IC = 1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(rts,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(rts,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PIP_NE_PT2_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PIP_NE_PT2 = CHI2_PIP_NE_PT2 + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PIP_NE_PT2 = NUM_PIP_NE_PT2 + 1
            WRITE(74,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(74)
      ENDIF

      IF(I_PIP_KR_PT2.EQ.1) THEN
      OPEN(UNIT = 75,FILE ='plot_data/HERMES/pip_kr_pt2.dat')
      WRITE(75,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
      DO I=1,9
          VAL  = PIP_KR_PT2_VALUE(I)
          MULT = PIP_KR_PT2_MULT(I)
          STAT = PIP_KR_PT2_STAT(I)
          SYS  = PIP_KR_PT2_SYS(I)
          Q2   = PIP_KR_PT2_Q2(I)
          Nu   = PIP_KR_PT2_NU(I)
          PT2  = VAL
          Z    = PIP_KR_PT2_Z(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          rts = Sep_hermes
          IT = 84
          IH = 1
          IC = 1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(rts,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(rts,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PIP_KR_PT2_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PIP_KR_PT2 = CHI2_PIP_KR_PT2 + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PIP_KR_PT2 = NUM_PIP_KR_PT2 + 1
            WRITE(75,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(75)
      ENDIF

      IF(I_PIP_XE_PT2.EQ.1) THEN
      OPEN(UNIT = 76,FILE ='plot_data/HERMES/pip_xe_pt2.dat')
      WRITE(76,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
      DO I=1,9
          VAL  = PIP_XE_PT2_VALUE(I)
          MULT = PIP_XE_PT2_MULT(I)
          STAT = PIP_XE_PT2_STAT(I)
          SYS  = PIP_XE_PT2_SYS(I)
          Q2   = PIP_XE_PT2_Q2(I)
          Nu   = PIP_XE_PT2_NU(I)
          PT2  = VAL
          Z    = PIP_XE_PT2_Z(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          rts = Sep_hermes
          IT = 131
          IH = 1
          IC = 1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(rts,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(rts,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PIP_XE_PT2_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PIP_XE_PT2 = CHI2_PIP_XE_PT2 + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PIP_XE_PT2 = NUM_PIP_XE_PT2 + 1
            WRITE(76,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(76)
      ENDIF

C------- NEGATIVELY CHARGED HADRONS
C------ PI - DATA
      IF(I_PIM_HE_PT2.EQ.1) THEN
      OPEN(UNIT = 85,FILE ='plot_data/HERMES/pim_he_pt2.dat')
      WRITE(85,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
      DO I=1,9
          VAL  = PIM_HE_PT2_VALUE(I)
          MULT = PIM_HE_PT2_MULT(I)
          STAT = PIM_HE_PT2_STAT(I)
          SYS  = PIM_HE_PT2_SYS(I)
          Q2   = PIM_HE_PT2_Q2(I)
          Nu   = PIM_HE_PT2_NU(I)
          PT2  = VAL
          Z    = PIM_HE_PT2_Z(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          rts = Sep_hermes
          IT = 4
          IH = 1
          IC = -1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(rts,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(rts,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PIM_HE_PT2_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PIM_HE_PT2 = CHI2_PIM_HE_PT2 + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PIM_HE_PT2 = NUM_PIM_HE_PT2 + 1
            WRITE(85,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(85)
      ENDIF

      IF(I_PIM_NE_PT2.EQ.1) THEN
      OPEN(UNIT = 86,FILE ='plot_data/HERMES/pim_ne_pt2.dat')
      WRITE(86,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
      DO I=1,9
          VAL  = PIM_NE_PT2_VALUE(I)
          MULT = PIM_NE_PT2_MULT(I)
          STAT = PIM_NE_PT2_STAT(I)
          SYS  = PIM_NE_PT2_SYS(I)
          Q2   = PIM_NE_PT2_Q2(I)
          Nu   = PIM_NE_PT2_NU(I)
          PT2  = VAL
          Z    = PIM_NE_PT2_Z(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          rts = Sep_hermes
          IT = 20
          IH = 1
          IC = -1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(rts,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(rts,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PIM_NE_PT2_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PIM_NE_PT2 = CHI2_PIM_NE_PT2 + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PIM_NE_PT2 = NUM_PIM_NE_PT2 + 1
            WRITE(86,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(86)
      ENDIF

      IF(I_PIM_KR_PT2.EQ.1) THEN
      OPEN(UNIT = 87,FILE ='plot_data/HERMES/pim_kr_pt2.dat')
      WRITE(87,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
      DO I=1,9
          VAL  = PIM_KR_PT2_VALUE(I)
          MULT = PIM_KR_PT2_MULT(I)
          STAT = PIM_KR_PT2_STAT(I)
          SYS  = PIM_KR_PT2_SYS(I)
          Q2   = PIM_KR_PT2_Q2(I)
          Nu   = PIM_KR_PT2_NU(I)
          PT2  = VAL
          Z    = PIM_KR_PT2_Z(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          rts = Sep_hermes
          IT = 84
          IH = 1
          IC = -1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(rts,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(rts,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PIM_KR_PT2_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PIM_KR_PT2 = CHI2_PIM_KR_PT2 + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PIM_KR_PT2 = NUM_PIM_KR_PT2 + 1
            WRITE(87,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(87)
      ENDIF

      IF(I_PIM_XE_PT2.EQ.1) THEN
      OPEN(UNIT = 88,FILE ='plot_data/HERMES/pim_xe_pt2.dat')
      WRITE(88,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
      DO I=1,9
          VAL  = PIM_XE_PT2_VALUE(I)
          MULT = PIM_XE_PT2_MULT(I)
          STAT = PIM_XE_PT2_STAT(I)
          SYS  = PIM_XE_PT2_SYS(I)
          Q2   = PIM_XE_PT2_Q2(I)
          Nu   = PIM_XE_PT2_NU(I)
          PT2  = VAL
          Z    = PIM_XE_PT2_Z(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          rts = Sep_hermes
          IT = 131
          IH = 1
          IC = -1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(rts,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(rts,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PIM_XE_PT2_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PIM_XE_PT2 = CHI2_PIM_XE_PT2 + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PIM_XE_PT2 = NUM_PIM_XE_PT2 + 1
            WRITE(88,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(88)
      ENDIF


C------- PI0 DATA
      IF(I_PI0_HE_PT2.EQ.1) THEN
      OPEN(UNIT = 203,FILE ='plot_data/HERMES/pi0_he_pt2.dat')
      WRITE(203,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
      DO I=1,9
          VAL  = PI0_HE_PT2_VALUE(I)
          MULT = PI0_HE_PT2_MULT(I)
          STAT = PI0_HE_PT2_STAT(I)
          SYS  = PI0_HE_PT2_SYS(I)
          Q2   = PI0_HE_PT2_Q2(I)
          Nu   = PI0_HE_PT2_NU(I)
          PT2  = VAL
          Z    = PI0_HE_PT2_Z(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          rts = Sep_hermes
          IT = 4
          IH = 1
          IC = 0
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(rts,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(rts,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PI0_HE_PT2_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PI0_HE_PT2 = CHI2_PI0_HE_PT2 + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PI0_HE_PT2 = NUM_PI0_HE_PT2 + 1
            WRITE(203,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(203)
      ENDIF

      IF(I_PI0_NE_PT2.EQ.1) THEN
      OPEN(UNIT = 204,FILE ='plot_data/HERMES/pi0_ne_pt2.dat')
      WRITE(204,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
      DO I=1,9
          VAL  = PI0_NE_PT2_VALUE(I)
          MULT = PI0_NE_PT2_MULT(I)
          STAT = PI0_NE_PT2_STAT(I)
          SYS  = PI0_NE_PT2_SYS(I)
          Q2   = PI0_NE_PT2_Q2(I)
          Nu   = PI0_NE_PT2_NU(I)
          PT2  = VAL
          Z    = PI0_NE_PT2_Z(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          rts = Sep_hermes
          IT = 20
          IH = 1
          IC = 0
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(rts,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(rts,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PI0_NE_PT2_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PI0_NE_PT2 = CHI2_PI0_NE_PT2 + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PI0_NE_PT2 = NUM_PI0_NE_PT2 + 1
            WRITE(204,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(204)
      ENDIF

      IF(I_PI0_KR_PT2.EQ.1) THEN
      OPEN(UNIT = 205,FILE ='plot_data/HERMES/pi0_kr_pt2.dat')
      WRITE(205,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
      DO I=1,9
          VAL  = PI0_KR_PT2_VALUE(I)
          MULT = PI0_KR_PT2_MULT(I)
          STAT = PI0_KR_PT2_STAT(I)
          SYS  = PI0_KR_PT2_SYS(I)
          Q2   = PI0_KR_PT2_Q2(I)
          Nu   = PI0_KR_PT2_NU(I)
          PT2  = VAL
          Z    = PI0_KR_PT2_Z(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          rts = Sep_hermes
          IT = 84
          IH = 1
          IC = 0
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(rts,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(rts,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PI0_KR_PT2_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PI0_KR_PT2 = CHI2_PI0_KR_PT2 + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PI0_KR_PT2 = NUM_PI0_KR_PT2 + 1
            WRITE(205,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(205)
      ENDIF

      IF(I_PI0_XE_PT2.EQ.1) THEN
      OPEN(UNIT = 206,FILE ='plot_data/HERMES/pi0_xe_pt2.dat')
      WRITE(206,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
      DO I=1,9
          VAL  = PI0_XE_PT2_VALUE(I)
          MULT = PI0_XE_PT2_MULT(I)
          STAT = PI0_XE_PT2_STAT(I)
          SYS  = PI0_XE_PT2_SYS(I)
          Q2   = PI0_XE_PT2_Q2(I)
          Nu   = PI0_XE_PT2_NU(I)
          PT2  = VAL
          Z    = PI0_XE_PT2_Z(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          rts = Sep_hermes
          IT = 131
          IH = 1
          IC = 0
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(rts,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(rts,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PI0_XE_PT2_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PI0_XE_PT2 = CHI2_PI0_XE_PT2 + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PI0_XE_PT2 = NUM_PI0_XE_PT2 + 1
            WRITE(206,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(206)
      ENDIF

      chisquare=CHI2

      CHI2_PIP_Z  =  CHI2_PIP_Z   + CHI2_PIP_HE_Z  + CHI2_PIP_NE_Z
     >                            + CHI2_PIP_KR_Z  + CHI2_PIP_XE_Z
      CHI2_PI0_Z  =  CHI2_PI0_Z   + CHI2_PI0_HE_Z   + CHI2_PI0_NE_Z
     >                            + CHI2_PI0_KR_Z   + CHI2_PI0_XE_Z
      CHI2_PIM_Z  =  CHI2_PIM_Z   + CHI2_PIM_HE_Z   + CHI2_PIM_NE_Z
     >                            + CHI2_PIM_KR_Z   + CHI2_PIM_XE_Z

      NUM_PIP_Z  =  NUM_PIP_Z   + NUM_PIP_HE_Z  + NUM_PIP_NE_Z
     >                            + NUM_PIP_KR_Z  + NUM_PIP_XE_Z
      NUM_PI0_Z  =  NUM_PI0_Z   + NUM_PI0_HE_Z   + NUM_PI0_NE_Z
     >                            + NUM_PI0_KR_Z   + NUM_PI0_XE_Z
      NUM_PIM_Z  =  NUM_PIM_Z   + NUM_PIM_HE_Z   + NUM_PIM_NE_Z
     >                            + NUM_PIM_KR_Z   + NUM_PIM_XE_Z

      CHI2_PIP_PT2 = CHI2_PIP_PT2 + CHI2_PIP_HE_PT2 + CHI2_PIP_NE_PT2
     >                            + CHI2_PIP_KR_PT2 + CHI2_PIP_XE_PT2
      CHI2_PI0_PT2 = CHI2_PI0_PT2 + CHI2_PI0_HE_PT2 + CHI2_PI0_NE_PT2
     >                            + CHI2_PI0_KR_PT2 + CHI2_PI0_XE_PT2
      CHI2_PIM_PT2 = CHI2_PIM_PT2 + CHI2_PIM_HE_PT2 + CHI2_PIM_NE_PT2
     >                            + CHI2_PIM_KR_PT2 + CHI2_PIM_XE_PT2

      NUM_PIP_PT2 = NUM_PIP_PT2 + NUM_PIP_HE_PT2 + NUM_PIP_NE_PT2
     >                            + NUM_PIP_KR_PT2 + NUM_PIP_XE_PT2

      NUM_PI0_PT2 = NUM_PI0_PT2 + NUM_PI0_HE_PT2 + NUM_PI0_NE_PT2
     >                            + NUM_PI0_KR_PT2 + NUM_PI0_XE_PT2

      NUM_PIM_PT2 = NUM_PIM_PT2 + NUM_PIM_HE_PT2 + NUM_PIM_NE_PT2
     >                            + NUM_PIM_KR_PT2 + NUM_PIM_XE_PT2

      NUM_HERMES = NUM_HERMES
     >       +NUM_PIP_HE_Z
     >       +NUM_PIP_NE_Z
     >       +NUM_PIP_KR_Z
     >       +NUM_PIP_XE_Z
     >       +NUM_PIM_HE_Z
     >       +NUM_PIM_NE_Z
     >       +NUM_PIM_KR_Z
     >       +NUM_PIM_XE_Z
     >       +NUM_PI0_HE_Z
     >       +NUM_PI0_NE_Z
     >       +NUM_PI0_KR_Z
     >       +NUM_PI0_XE_Z
     >       +NUM_PIP_HE_PT2
     >       +NUM_PIP_NE_PT2
     >       +NUM_PIP_KR_PT2
     >       +NUM_PIP_XE_PT2
     >       +NUM_PIM_HE_PT2
     >       +NUM_PIM_NE_PT2
     >       +NUM_PIM_KR_PT2
     >       +NUM_PIM_XE_PT2
     >       +NUM_PI0_HE_PT2
     >       +NUM_PI0_NE_PT2
     >       +NUM_PI0_KR_PT2
     >       +NUM_PI0_XE_PT2

      num_exp = 0
     >       +num_E772_800
     >       +num_E866_800q
     >       +num_HERMES
     >       +NUM_RHIC_Ratio_pAu1
     >       +NUM_RHIC_Ratio_pAu2
     >       +NUM_ATLAS5_Y1
     >       +NUM_ATLAS5_Y2
     >       +NUM_ATLAS5_Y3
     >       +NUM_CMS5
     >       +NUM_JLAB12

c-----format to write out the final result-------f
 100  format(7(1pE12.4))
 101  format(6(1pE10.2))
 110  format(A35,F10.3,I2,I2)
 112  format(A35,I5)
 113  format(A35,13(1pe12.4,','))
 114  format(A35,F10.3,',',I5,',',I5,',',F10.3)
 115  format(A25,F10.3,I2,',',A25,F10.3,I2)
! 116  format(3(A25,F10.3,I2,I2,','),A25,F10.3,I2,I2,',',A25,F10.3,I2,I2)
 116  format(3(A25,F10.3,I2,I2,','))
 117  format(4(A25,F10.3,I2,I2,','))
 1117 format(4(A25,F10.3,I2,','))
c-----write out the partial result on the screen-------
c         write(* ,  *) '---------------------------------'
c         write(* ,114) 'chi2,points,params,chi2/dof:',
c     >   chisquare,num_exp,nfit,chisquare/(num_exp-nfit)
c         write(* ,  *) '---------------------------------'
c      write(*,*) 'aN,aT,bT=',xx( 1),xx( 2),xx( 3)
c
c         WRITE(*,  *) '-------------------------'
c         write(* ,  *) '---------------------------------'


c         write(* ,  *) '---------------------------------'
c         write(* ,114) 'chi2,points,params,chi2/dof:',
c     >   chisquare,num_exp,nfit,chisquare/(num_exp-nfit)
c         write(* ,  *) '---------------------------------'
c      write(*,*) 'aN , bN, g2A, a2, g2B, B2',
c     ^            xx(1), xx(2), xx(3), xx(4), xx(5), xx(6)
c         WRITE(*,  *) '-------------------------'
c         write(* ,  *) '---------------------------------'

         write(6,110) 'chi2:',chisquare
         write(6,112) 'total num of data in fit:',num_exp
         write(6,112) 'number of parameters:',nfit
         write(6,110) 'chi2/dof:',chisquare/(num_exp-nfit)
         WRITE(6,  *) '-------------------------'
         write(6,113) 'fitted parameters:'
         write(6,*) 'aN, bN, g2A, a2, g2B, B2'
         write(6,*) xx(1), xx(2), xx(3), xx(4), xx(5), xx(6)
         WRITE(6,  *) '-------------------------'
      write(6,*) 'E866 (Q bins):' ,CHI2_E866_800q, NUM_E866_800q
      write(6,*) 'E772 (800 GeV):',CHI2_E772_800,  NUM_E772_800
      write(6,*) 'RHIC(Ratio) Au going:',
     *                      CHI2_RHIC_Ratio_pAu1, NUM_RHIC_Ratio_pAu1
      write(6,*) 'RHIC(Ratio) p going:',
     *                      CHI2_RHIC_Ratio_pAu2, NUM_RHIC_Ratio_pAu2
      write(6,*) 'ATLAS5 Y1:',CHI2_ATLAS5_Y1, NUM_ATLAS5_Y1
      write(6,*) 'ATLAS5 Y2:',CHI2_ATLAS5_Y2, NUM_ATLAS5_Y2
      write(6,*) 'ATLAS5 Y3:',CHI2_ATLAS5_Y3, NUM_ATLAS5_Y3
      write(6,*) 'CMS5:',CHI2_CMS5, NUM_CMS5
      write(6,*) 'HERMES (SIDIS):',  CHI2_HERMES, NUM_HERMES
      write(6,*) 'HERMES z'
      write(6,*) 'HERMES (PIP Nu-He):',CHI2_PIP_HE_Z  ,NUM_PIP_HE_Z
      write(6,*) 'HERMES (PIP Nu-NE):',CHI2_PIP_NE_Z  ,NUM_PIP_NE_Z
      write(6,*) 'HERMES (PIP Nu-KR):',CHI2_PIP_KR_Z  ,NUM_PIP_KR_Z
      write(6,*) 'HERMES (PIP Nu-XE):',CHI2_PIP_XE_Z  ,NUM_PIP_XE_Z
      write(6,*) 'HERMES (PI0 Nu-He):',CHI2_PI0_HE_Z  ,NUM_PI0_HE_Z
      write(6,*) 'HERMES (PI0 Nu-NE):',CHI2_PI0_NE_Z  ,NUM_PI0_NE_Z
      write(6,*) 'HERMES (PI0 Nu-KR):',CHI2_PI0_KR_Z  ,NUM_PI0_KR_Z
      write(6,*) 'HERMES (PI0 Nu-XE):',CHI2_PI0_XE_Z  ,NUM_PI0_XE_Z
      write(6,*) 'HERMES (PIM Nu-He):',CHI2_PIM_HE_Z  ,NUM_PIM_HE_Z
      write(6,*) 'HERMES (PIM Nu-NE):',CHI2_PIM_NE_Z  ,NUM_PIM_NE_Z
      write(6,*) 'HERMES (PIM Nu-KR):',CHI2_PIM_KR_Z  ,NUM_PIM_KR_Z
      write(6,*) 'HERMES (PIM Nu-XE):',CHI2_PIM_XE_Z  ,NUM_PIM_XE_Z
      write(6,*) 'HERMES PT2'
      write(6,*) 'HERMES (PIP Nu-He):',CHI2_PIP_HE_PT2 ,NUM_PIP_HE_PT2
      write(6,*) 'HERMES (PIP Nu-NE):',CHI2_PIP_NE_PT2  ,NUM_PIP_NE_PT2
      write(6,*) 'HERMES (PIP Nu-KR):',CHI2_PIP_KR_PT2  ,NUM_PIP_KR_PT2
      write(6,*) 'HERMES (PIP Nu-XE):',CHI2_PIP_XE_PT2  ,NUM_PIP_XE_PT2
      write(6,*) 'HERMES (PI0 Nu-He):',CHI2_PI0_HE_PT2  ,NUM_PI0_HE_PT2
      write(6,*) 'HERMES (PI0 Nu-NE):',CHI2_PI0_NE_PT2  ,NUM_PI0_NE_PT2
      write(6,*) 'HERMES (PI0 Nu-KR):',CHI2_PI0_KR_PT2  ,NUM_PI0_KR_PT2
      write(6,*) 'HERMES (PI0 Nu-XE):',CHI2_PI0_XE_PT2  ,NUM_PI0_XE_PT2
      write(6,*) 'HERMES (PIM Nu-He):',CHI2_PIM_HE_PT2  ,NUM_PIM_HE_PT2
      write(6,*) 'HERMES (PIM Nu-NE):',CHI2_PIM_NE_PT2  ,NUM_PIM_NE_PT2
      write(6,*) 'HERMES (PIM Nu-KR):',CHI2_PIM_KR_PT2  ,NUM_PIM_KR_PT2
      write(6,*) 'HERMES (PIM Nu-XE):',CHI2_PIM_XE_PT2  ,NUM_PIM_XE_PT2

c-----write out the final result on a file-------
      if(iflag.eq.3) then

         open(unit=82,file='result.txt',status='unknown')
         write(82,110) 'chi2:',chisquare
         write(82,112) 'total num of data in fit:',num_exp
         write(82,112) 'number of parameters:',nfit
         write(82,110) 'chi2/dof:',chisquare/(num_exp-nfit)
         WRITE(82,  *) '-------------------------'
      write(82,113) 'fitted parameters:'
      write(82,*) 'aN, bN, g2A, a2, g2B, B2'
      write(82,*) xx(1), xx(2), xx(3), xx(4), xx(5), xx(6)
         WRITE(82,  *) '-------------------------'

         write(82,*) 'RHIC(Ratio) Au going:',
     *                      CHI2_RHIC_Ratio_pAu1, NUM_RHIC_Ratio_pAu1
         write(82,*) 'RHIC(Ratio) p going:',
     *                      CHI2_RHIC_Ratio_pAu2, NUM_RHIC_Ratio_pAu2
         write(82,*) 'ATLAS5 Y1:',CHI2_ATLAS5_Y1, NUM_ATLAS5_Y1
         write(82,*) 'NORMALIZATION:', ATLAS5_Y1_NORM
         write(82,*) 'ATLAS5 Y2:',CHI2_ATLAS5_Y2, NUM_ATLAS5_Y2
         write(82,*) 'NORMALIZATION:', ATLAS5_Y2_NORM
         write(82,*) 'ATLAS5 Y3:',CHI2_ATLAS5_Y3, NUM_ATLAS5_Y3
         write(82,*) 'NORMALIZATION:', ATLAS5_Y3_NORM
         write(82,*) 'CMS5:',CHI2_CMS5, NUM_CMS5
         write(82,*) 'NORMALIZATION:', CMS5_NORM

         write(82,*) 'E866 (Q bins):' ,CHI2_E866_800q, NUM_E866_800q
         write(82,*) 'E772 (800 GeV):',CHI2_E772_800,  NUM_E772_800
         write(82,*) 'HERMES (SIDIS):',  CHI2_HERMES, NUM_HERMES
         write(82,*) 'HERMES z'
         write(82,*) 'HERMES (PIP Z-He):',CHI2_PIP_HE_Z, NUM_PIP_HE_Z
         write(82,*) 'HERMES (PIP Z-Ne):',CHI2_PIP_NE_Z, NUM_PIP_NE_Z
         write(82,*) 'HERMES (PIP Z-Kr):',CHI2_PIP_KR_Z, NUM_PIP_KR_Z
         write(82,*) 'HERMES (PIP Z-XE):',CHI2_PIP_XE_Z, NUM_PIP_XE_Z
         write(82,*) 'HERMES (PI0 Z-He):',CHI2_PI0_HE_Z, NUM_PI0_HE_Z
         write(82,*) 'HERMES (PI0 Z-Ne):',CHI2_PI0_NE_Z, NUM_PI0_NE_Z
         write(82,*) 'HERMES (PI0 Z-Kr):',CHI2_PI0_KR_Z, NUM_PI0_KR_Z
         write(82,*) 'HERMES (PI0 Z-XE):',CHI2_PI0_XE_Z, NUM_PI0_XE_Z
         write(82,*) 'HERMES (PIM Z-He):',CHI2_PIM_HE_Z, NUM_PIM_HE_Z
         write(82,*) 'HERMES (PIM Z-Ne):',CHI2_PIM_NE_Z, NUM_PIM_NE_Z
         write(82,*) 'HERMES (PIM Z-Kr):',CHI2_PIM_KR_Z, NUM_PIM_KR_Z
         write(82,*) 'HERMES (PIM Z-XE):',CHI2_PIM_XE_Z, NUM_PIM_XE_Z
         write(82,*) 'HERMES PT2'
         write(82,*) 'HERMES (PIP Nu-He):',CHI2_PIP_HE_PT2,
     *      NUM_PIP_HE_PT2
         write(82,*) 'HERMES (PIP Nu-NE):',CHI2_PIP_NE_PT2  ,
     *      NUM_PIP_NE_PT2
         write(82,*) 'HERMES (PIP Nu-KR):',CHI2_PIP_KR_PT2  ,
     *      NUM_PIP_KR_PT2
         write(82,*) 'HERMES (PIP Nu-XE):',CHI2_PIP_XE_PT2  ,
     *      NUM_PIP_XE_PT2
         write(82,*) 'HERMES (PI0 Nu-He):',CHI2_PI0_HE_PT2  ,
     *      NUM_PI0_HE_PT2
         write(82,*) 'HERMES (PI0 Nu-NE):',CHI2_PI0_NE_PT2  ,
     *      NUM_PI0_NE_PT2
         write(82,*) 'HERMES (PI0 Nu-KR):',CHI2_PI0_KR_PT2  ,
     *      NUM_PI0_KR_PT2
         write(82,*) 'HERMES (PI0 Nu-XE):',CHI2_PI0_XE_PT2  ,
     *      NUM_PI0_XE_PT2
         write(82,*) 'HERMES (PIM Nu-He):',CHI2_PIM_HE_PT2  ,
     *      NUM_PIM_HE_PT2
         write(82,*) 'HERMES (PIM Nu-NE):',CHI2_PIM_NE_PT2  ,
     *      NUM_PIM_NE_PT2
         write(82,*) 'HERMES (PIM Nu-KR):',CHI2_PIM_KR_PT2
     *      ,NUM_PIM_KR_PT2
         write(82,*) 'HERMES (PIM Nu-XE):',CHI2_PIM_XE_PT2  ,
     *      NUM_PIM_XE_PT2
         close(82)

         close(82)

      endif

      return
      end
