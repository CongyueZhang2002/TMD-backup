      PROGRAM MASTER
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL FCNG1
      DIMENSION NPRM(12),VSTRT(12),STP(12),BL(12),BUP(12),ARGLIS(12) !2023
      CHARACTER*10 PNAM(12) !2023
      INTEGER IRD,IWR,ISAV
      INTEGER NFIT,NFIT1
      INTEGER DAY, HOUR, MINUTE, SECOND
      COMMON /PARAMS/ NFIT
      COMMON /new/ NFIT1 !2023
      INTEGER IFLAG
      COMMON  / IEND / IFLAG
      real chi2

C     INITALIZE THE PARAMETRIZATION

C-----Parameter Names
      REAL*8 gamma, g3f, g3D, Nq1, gq1, dq1,
     &       Nq2, gq2, dq2, p_10, p_11, p_12 !2023 
      COMMON /FITP/ gamma, g3f, g3D, Nq1, gq1, dq1,
     &              Nq2, gq2, dq2, p_10, p_11, p_12 !2023 
      data NPRM /1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/ !2023
      data PNAM /'gamma', 'g3f', 'g3D', 
     &           'Nq1', 'gq1', 'dq1',
     &           'Nq2', 'gq2', 'dq2',  
     &           'p_10', 'p_11', 'p_12'/ !2023

C-----STARTING STEP SIZE 
      data STP /0.1, 0.005, 0.005, 
     &          0.2, 0.2, 0.2,
     &          0.1, 0.1, 0.1, 
     &          0.2, 0.1, 100.0/ !2023

C-----LOWER BOUND (LIMIT) ON PARAMETER VALUE, IF ANY
      data BL   /0.0, 0.0, 0.0, 
     &           -5.0, -5.0, -5.0,
     &           -2.0, -2.0, -2.0, 
     &           -5.0, -2.0, -100.0/ !2023

C-----UPPER BOUND (LIMIT) ON PARAMETER VALUE, IF ANY
      data BUP  /100.0, 0.1, 0.1,  
     &           5.0, 5.0, 5.0,
     &           2.0, 2.0, 2.0,  
     &           5.0, 2.0, 100.0/ !2023
     
C-----Local Variable Declarations
      include "tools/data-inc.f"
      integer nloops,hop,nll,prepdf,preff, pre
      integer collFF, LIKEn, newFF, newAPFEL
      common /scheme/ nloops,hop,nll,pre
      common /collinear/ collFF, LIKEn, newFF, newAPFEL
      REAL*8 gammaMAX ,gammaMIN !2023
      REAL*8 g3fMAX   ,g3fMIN
      REAL*8 g3DMAX   ,g3DMIN
      REAL*8 Nq1MAX   ,Nq1MIN
      REAL*8 gq1MAX   ,gq1MIN
      REAL*8 dq1MAX   ,dq1MIN
      REAL*8 Nq2MAX   ,Nq2MIN
      REAL*8 gq2MAX   ,gq2MIN
      REAL*8 dq2MAX   ,dq2MIN
      REAL*8 p_10MAX  ,p_10MIN
      REAL*8 p_11MAX  ,p_11MIN
      REAL*8 p_12MAX  ,p_12MIN
      integer seed,initseed
      REAL*8 XX(12) !2023

      INTEGER NUM_EXP
      COMMON /NNUM_EXP/ NUM_EXP
      real*8 r0,r01000

      REAL*8 r1, r2, r3, r4, r5, r6 !2023
      REAL*8 r7, r8, r9, r10, r11, r12 

      REAL*8 gammaBEST,g3fBEST,g3DBEST,Nq1BEST,gq1BEST,dq1BEST,
     &       Nq2BEST,gq2BEST,dq2BEST,p_10BEST,p_11BEST,p_12BEST !2023

C-------------Actual Number of Parameters Used !2023
      NFIT1 = 10
C-------------

 1115 CONTINUE

C      INCLUDE "best-params-rand.f"

C-----CHOOSE EXP DATA SETS
C-----SETS TO INCLUDE IN INITALIZATION (DO NOT CHANGE)
C-----These are not sets which enter the fit
C-----These are only sets to include when testing the code initially
C-----To check that the parameter values chosen (randomly) do not cause the program to break
C-----To choose the sets for the fit, go to line 263.

C-----RHIC Ratios
      I_RHIC_Ratio_pAu1 = 1 !1
      I_RHIC_Ratio_pAu2 = 1 !1

C-----ATLAS 5 TEV
      I_ATLAS5_Y1 = 1 !1
      I_ATLAS5_Y2 = 0
      I_ATLAS5_Y3 = 0

C-----CMS 5 TEV
      I_CMS5 = 1 !1

C-----E866
      I_E866_800q  = 1 !1 !1

C-----E772
      I_E772_800   = 1

C-----SIDIS (HERMES)
      I_PIP_HE_Z   = 0
      I_PIP_NE_Z   = 1
      I_PIP_KR_Z   = 1
      I_PIP_XE_Z   = 1
      I_PI0_HE_Z   = 0
      I_PI0_NE_Z   = 1
      I_PI0_KR_Z   = 1
      I_PI0_XE_Z   = 1
      I_PIM_HE_Z   = 0
      I_PIM_NE_Z   = 1
      I_PIM_KR_Z   = 1
      I_PIM_XE_Z   = 1


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

C------ J-LAB 2022
      I_pip_2022_C  = 0
      I_pip_2022_Fe = 0
      I_pip_2022_Pb = 0
      I_pim_2022_C  = 0
      I_pim_2022_Fe = 0
      I_pim_2022_Pb = 0
C----------------------------------------------------------------------
      LIKEn = 0 
      newFF = 1 
      newAPFEL = 2 

c-----FIT PARAMETERS
      nloops = 2 ! LO: 1 NLO: 2
      nll = 3
      pre = 1
      collFF = newFF

!     CALL SETLHAPARM('SILENT')!TO NOT SHOW THE CALLS, ALTHOUGH THEY ARE CALLED
      CALL SETLHAPARM('SILENT')

C------TO NOT SHOW THE APFEL CALLS, ALTHOUGH THEY ARE CALLED
      CALL EnableWelcomeMessage(.False.)      

C-----INITIALIZE PDF SETS USING LHAPDF
      print *, "1"
      CALL InitPDFsetByNameM(1,"CT18ANLO")
      CALL InitPDFM(1,0)
      print *, "2"
      CALL InitPDFsetByNameM(2,"EPPS21HE")
      CALL InitPDFM(2,0)
      print *, "3"
      CALL InitPDFsetByNameM(3,"EPPS21NE")
      CALL InitPDFM(3,0)
      print *, "4"
      CALL InitPDFsetByNameM(4,"EPPS21KR")
      CALL InitPDFM(4,0)
      print *, "5"
      CALL InitPDFsetByNameM(5,"EPPS21XE")
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
      CALL InitPDFsetByNameM(10,"EPPS21BE")
      CALL InitPDFM(10,0)
      print *, "11"
      CALL InitPDFsetByNameM(11,"EPPS21FE")
      CALL InitPDFM(11,0)
      print *, "12"
      CALL InitPDFsetByNameM(12,"EPPS21WW")
      CALL InitPDFM(12,0)
      print *, "13"
      CALL InitPDFsetByNameM(13,"EPPS21JCC")
      CALL InitPDFM(13,0)
      print *, "14"
      CALL InitPDFsetByNameM(14,"EPPS21JFE")
      CALL InitPDFM(14,0)
      print *, "15"
      CALL InitPDFsetByNameM(15,"EPPS21JPB")
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
      print *, "20"
      CALL InitPDFsetByNameM(20,"EPPS21AU")
      CALL InitPDFM(20,0)
      print *, "21"
      CALL InitPDFsetByNameM(21,"EPPS21PR")
      CALL InitPDFM(21,0)
      print *, "22"
      CALL InitPDFsetByNameM(22,"EPPS21CA")
      CALL InitPDFM(22,0)
      print *, "23"
      CALL InitPDFsetByNameM(23,"LIKEnVC")
      CALL InitPDFM(23,0)
      CALL InitPDFsetByNameM(24,"LIKEnAU")
      CALL InitPDFM(24,0)

      CALL READDATA

*.....CALCULATE THE STARTING TIME OF THE FIT
      CALL TIMESTAMP( )         ! WRITE OUT THE STARTING TIME
      CALL CPU_TIME(TIME1)

      print*, "Making sure chi2 is not NaN"

      INCLUDE "best-params-rand.f" ! 2023
      XX(1)  = gammaBEST
      XX(2)  = g3fBEST
      XX(3)  = g3DBEST
      XX(4)  = Nq1BEST
      XX(5)  = gq1BEST
      XX(6)  = dq1BEST
      XX(7)  = Nq2BEST
      XX(8)  = gq2BEST
      XX(9)  = dq2BEST
      XX(10) = p_10BEST
      XX(11) = p_11BEST
      XX(12) = p_12BEST

      chi2 = 0
      chi2 = chisquare_calculate(XX)
      print*, "Initializing"
      print*, "chi2/dof = ", chi2/100d0 ! 2023, Need to change here if include different dat sets
      print*, "gamma,g3f,g3D,Nq1,gq1,dq1,Nq2,gq2,dq2,p_10,p_11,p_12" !2023
      print*, gamma,g3f,g3D,Nq1,gq1,dq1,Nq2,gq2,dq2,p_10,p_11,p_12 
      do while( (chi2/100d0 .GT. 100d0) .OR. (isNAN(chi2)) )
          chi2 = 0
          INCLUDE "best-params-rand.f"
          XX(1)  = gammaBEST
          XX(2)  = g3fBEST
          XX(3)  = g3DBEST
          XX(4)  = Nq1BEST
          XX(5)  = gq1BEST
          XX(6)  = dq1BEST
          XX(7)  = Nq2BEST
          XX(8)  = gq2BEST
          XX(9)  = dq2BEST
          XX(10) = p_10BEST
          XX(11) = p_11BEST
          XX(12) = p_12BEST
          chi2 = chisquare_calculate(XX)
          print*, "chi2/dof = ", chi2/100d0
          print*, "gamma,g3f,g3D,Nq1,gq1,dq1,Nq2,gq2,dq2,p_10,p_11,p_12" !2023
          print*, gamma,g3f,g3D,Nq1,gq1,dq1,Nq2,gq2,dq2,p_10,p_11,p_12 
      end do

      print*, "It's not NaN"

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     START FIT
*.....INITIALIZATION:
      IRD = 5         ! UNIT NUMBER FOR INPUT TO MINUIT (5 KEYBOARD)
      IWR = 6         ! UNIT NUMBER FOR OUTPUT FROM MINUIT (6 SCREEN)
      ISAV = 7        ! UNIT NUMBER FOR USE OF THE SAVE COMMAND
      CALL MNINIT (IRD,IWR,ISAV)
*.....DEFINITON OF THE PARAMETERS :
      DO 11 I = 1, 12 ! 2023
         CALL MNPARM (NPRM(I),PNAM(I),VSTRT(I),
     >                        STP(I),BL(I),BUP(I),IERFLG)
         IF (IERFLG .NE. 0) THEN
            WRITE (6,*) 'UNABLE TO DEFINE PARAMETER NO.', I
            STOP
         END IF
 11   CONTINUE

C-----Parameters 10-12 Haven't Been Activated
      CALL MNFIXP(1)
      CALL MNFIXP(11)

*.....OUTPUT PRINT LEVEL (FROM -1 TO 3)
      ARGLIS(1) = 3.
      CALL MNEXCM (FCNG1,'SET PRINT',ARGLIS, 1, IERFLG, DUM)
*.....FIRST CALL :
      ARGLIS(1) = 1.            !   IFLAG = 1
      CALL MNEXCM (FCNG1,'CALL FCN', ARGLIS, 1, IERFLG, DUM)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     NUMBER OF PARAMETERS
      NFIT = 12 !2023
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*.....SIMPLEX FIT :
      ARGLIS(1) = 2000.
      CALL MNEXCM (FCNG1, 'SIMPLEX', ARGLIS, 2, IERFLG, DUM)

C---- Choose sets to include in fit.
C-----RHIC Ratios
      I_RHIC_Ratio_pAu1 = 1 !1
      I_RHIC_Ratio_pAu2 = 1 !1

C-----ATLAS 5 TEV
      I_ATLAS5_Y1 = 1 !1
      I_ATLAS5_Y2 = 0
      I_ATLAS5_Y3 = 0

C-----CMS 5 TEV
      I_CMS5 = 1 !1

C-----E866
      I_E866_800q  = 1 !1 !1

C-----E772
      I_E772_800   = 1

C-----SIDIS (HERMES)
      I_PIP_HE_Z   = 0
      I_PIP_NE_Z   = 1
      I_PIP_KR_Z   = 1
      I_PIP_XE_Z   = 1
      I_PI0_HE_Z   = 0
      I_PI0_NE_Z   = 1
      I_PI0_KR_Z   = 1
      I_PI0_XE_Z   = 1
      I_PIM_HE_Z   = 0
      I_PIM_NE_Z   = 1
      I_PIM_KR_Z   = 1
      I_PIM_XE_Z   = 1

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

C------ J-LAB 2022
      I_pip_2022_C = 0
      I_pip_2022_Fe = 0
      I_pip_2022_Pb = 0
      I_pim_2022_C = 0
      I_pim_2022_Fe = 0
      I_pim_2022_Pb = 0

*.....MINIMIZE FIT :
!     IT IS EQUIVALENT TO MIGRAD, EXCEPT THAT IF MIGRAD FAILS,
*     IT REVEs TO SIMPLEX AND THEN CALLS MIGRAD AGAIN
      ARGLIS(1) = 4000.
      CALL MNEXCM (FCNG1, 'MINIMIZE', ARGLIS, 1, IERFLG, DUM)
*     CHECK THAT RESULT HAS CONVERGED
      XX(1)  = gammaBEST !2023
      XX(2)  = g3fBEST
      XX(3)  = g3DBEST
      XX(4)  = Nq1BEST
      XX(5)  = gq1BEST
      XX(6)  = dq1BEST
      XX(7)  = Nq2BEST
      XX(8)  = gq2BEST
      XX(9)  = dq2BEST
      XX(10) = p_10BEST
      XX(11) = p_11BEST
      XX(12) = p_12BEST
      IF (CHISQUARE(XX)/NUM_EXP.GT.20D0) THEN
      GOTO 1115
      ELSE IF (ISNAN(CHISQUARE(XX))) THEN
      GOTO 1115
      ENDIF

C---- This call creates the files for plotting as well as result.txt 
*.....LAST CALL :
*     THIS CREATES THE FILE WITH THE RESULTS
      ARGLIS(1) = 3.            !   IFLAG = 3
      CALL MNEXCM (FCNG1, 'CALL FCN', ARGLIS, 1, IERFLG, DUM)
*.....STOP
      CALL MNEXCM (FCNG1, 'STOP', ARGLIS, 1, IERFLG, DUM)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

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

      gammaBEST = gamma!2023
      g3fBEST   = g3f
      g3DBEST   = g3D
      Nq1BEST   = Nq1
      gq1BEST   = gq1
      dq1BEST   = dq1
      Nq2BEST   = Nq2
      gq2BEST   = gq2
      dq2BEST   = dq2
      p_10BEST  = p_10
      p_11BEST  = p_11
      p_12BEST  = p_12

      OPEN(UNIT=999,FILE='./best-params.f') !2023
      WRITE(999, *) "gammaBEST = ", gammaBEST
      WRITE(999, *) "g3fBEST =   ", g3fBEST
      WRITE(999, *) "g3DBEST =   ", g3DBEST
      WRITE(999, *) "Nq1BEST =   ", Nq1BEST
      WRITE(999, *) "gq1BEST =   ", gq1BEST
      WRITE(999, *) "dq1BEST =   ", dq1BEST
      WRITE(999, *) "Nq2BEST =   ", Nq2BEST
      WRITE(999, *) "gq2BEST =   ", gq2BEST
      WRITE(999, *) "dq2BEST =   ", dq2BEST
      WRITE(999, *) "p_10BEST =  ", p_10BEST
      WRITE(999, *) "p_11BEST =  ", p_11BEST
      WRITE(999, *) "p_12BEST =  ", p_12BEST

      CLOSE(999)

      RETURN
      END

C-----

      SUBROUTINE FCNG1 (NPAR, G, F, X, IIFLAG, dum)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION X(12), G(12) !2023
      integer iflag,iiflag
      COMMON  / iend / iflag

      iflag = iiflag
      F = chisquare(X)

      RETURN
      END

C-----

      function chisquare(XX)
      IMPLICIT NONE
      integer nfit,nfit1
      common /params/ nfit
      common /new/ nfit1
      real*8 xx(nfit)
      real*8 chisquare
      real*8 ktw_value
      real*8  ptw_value

      real*8 CHI2_E866_800q, CHI2_E772_800
      REAL*8 CHI2_HERMES
      REAL*8 CHI2_PIP_HE_NU, CHI2_PIP_HE_Z, CHI2_PIP_HE_Q2
      REAL*8 CHI2_PIP_NE_NU, CHI2_PIP_NE_Z, CHI2_PIP_NE_Q2
      REAL*8 CHI2_PIP_KR_NU, CHI2_PIP_KR_Z, CHI2_PIP_KR_Q2
      REAL*8 CHI2_PIP_XE_NU, CHI2_PIP_XE_Z, CHI2_PIP_XE_Q2
      REAL*8 CHI2_PIM_HE_NU, CHI2_PIM_HE_Z, CHI2_PIM_HE_Q2
      REAL*8 CHI2_PIM_NE_NU, CHI2_PIM_NE_Z, CHI2_PIM_NE_Q2
      REAL*8 CHI2_PIM_KR_NU, CHI2_PIM_KR_Z, CHI2_PIM_KR_Q2
      REAL*8 CHI2_PIM_XE_NU, CHI2_PIM_XE_Z, CHI2_PIM_XE_Q2
      REAL*8 CHI2_PI0_HE_NU, CHI2_PI0_HE_Z, CHI2_PI0_HE_Q2
      REAL*8 CHI2_PI0_NE_NU, CHI2_PI0_NE_Z, CHI2_PI0_NE_Q2
      REAL*8 CHI2_PI0_KR_NU, CHI2_PI0_KR_Z, CHI2_PI0_KR_Q2
      REAL*8 CHI2_PI0_XE_NU, CHI2_PI0_XE_Z, CHI2_PI0_XE_Q2
      REAL*8 CHI2_PIP_HE_PT2, CHI2_PIM_HE_PT2, CHI2_PI0_HE_PT2
      REAL*8 CHI2_PIP_NE_PT2, CHI2_PIM_NE_PT2, CHI2_PI0_NE_PT2
      REAL*8 CHI2_PIP_KR_PT2, CHI2_PIM_KR_PT2, CHI2_PI0_KR_PT2
      REAL*8 CHI2_PIP_XE_PT2, CHI2_PIM_XE_PT2, CHI2_PI0_XE_PT2
      REAL*8 CHI2_KP_HE_NU, CHI2_KP_HE_Z, CHI2_KP_HE_Q2
      REAL*8 CHI2_KP_NE_NU, CHI2_KP_NE_Z, CHI2_KP_NE_Q2
      REAL*8 CHI2_KP_KR_NU, CHI2_KP_KR_Z, CHI2_KP_KR_Q2
      REAL*8 CHI2_KP_XE_NU, CHI2_KP_XE_Z, CHI2_KP_XE_Q2
      REAL*8 CHI2_KM_HE_NU, CHI2_KM_HE_Z, CHI2_KM_HE_Q2
      REAL*8 CHI2_KM_NE_NU, CHI2_KM_NE_Z, CHI2_KM_NE_Q2
      REAL*8 CHI2_KM_KR_NU, CHI2_KM_KR_Z, CHI2_KM_KR_Q2
      REAL*8 CHI2_KM_XE_NU, CHI2_KM_XE_Z, CHI2_KM_XE_Q2
      REAL*8 CHI2_KP_HE_PT2, CHI2_KM_HE_PT2
      REAL*8 CHI2_KP_NE_PT2, CHI2_KM_NE_PT2
      REAL*8 CHI2_KP_KR_PT2, CHI2_KM_KR_PT2
      REAL*8 CHI2_KP_XE_PT2, CHI2_KM_XE_PT2
      REAL*8 CHI2_PIP_NU
      REAL*8 CHI2_PIP_Z
      REAL*8 CHI2_PIP_Q2
      REAL*8 CHI2_PIP_PT2
      REAL*8 CHI2_PI0_NU
      REAL*8 CHI2_PI0_Z
      REAL*8 CHI2_PI0_Q2
      REAL*8 CHI2_PI0_PT2
      REAL*8 CHI2_PIM_NU
      REAL*8 CHI2_PIM_Z
      REAL*8 CHI2_PIM_Q2
      REAL*8 CHI2_PIM_PT2
      REAL*8 CHI2_PIP_2022_C !2023
      REAL*8 CHI2_PIP_2022_Fe
      REAL*8 CHI2_PIP_2022_Pb
      REAL*8 CHI2_PIM_2022_C
      REAL*8 CHI2_PIM_2022_Fe
      REAL*8 CHI2_PIM_2022_Pb

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
      INTEGER NUM_PIP_HE_NU, NUM_PIP_HE_Z, NUM_PIP_HE_Q2
      INTEGER NUM_PIP_NE_NU, NUM_PIP_NE_Z, NUM_PIP_NE_Q2
      INTEGER NUM_PIP_KR_NU, NUM_PIP_KR_Z, NUM_PIP_KR_Q2
      INTEGER NUM_PIP_XE_NU, NUM_PIP_XE_Z, NUM_PIP_XE_Q2
      INTEGER NUM_PIM_HE_NU, NUM_PIM_HE_Z, NUM_PIM_HE_Q2
      INTEGER NUM_PIM_NE_NU, NUM_PIM_NE_Z, NUM_PIM_NE_Q2
      INTEGER NUM_PIM_KR_NU, NUM_PIM_KR_Z, NUM_PIM_KR_Q2
      INTEGER NUM_PIM_XE_NU, NUM_PIM_XE_Z, NUM_PIM_XE_Q2
      INTEGER NUM_PI0_HE_NU, NUM_PI0_HE_Z, NUM_PI0_HE_Q2
      INTEGER NUM_PI0_NE_NU, NUM_PI0_NE_Z, NUM_PI0_NE_Q2
      INTEGER NUM_PI0_KR_NU, NUM_PI0_KR_Z, NUM_PI0_KR_Q2
      INTEGER NUM_PI0_XE_NU, NUM_PI0_XE_Z, NUM_PI0_XE_Q2
      INTEGER NUM_PIP_HE_PT2, NUM_PIM_HE_PT2, NUM_PI0_HE_PT2
      INTEGER NUM_PIP_NE_PT2, NUM_PIM_NE_PT2, NUM_PI0_NE_PT2
      INTEGER NUM_PIP_KR_PT2, NUM_PIM_KR_PT2, NUM_PI0_KR_PT2
      INTEGER NUM_PIP_XE_PT2, NUM_PIM_XE_PT2, NUM_PI0_XE_PT2
      INTEGER NUM_KP_HE_NU, NUM_KP_HE_Z, NUM_KP_HE_Q2
      INTEGER NUM_KP_NE_NU, NUM_KP_NE_Z, NUM_KP_NE_Q2
      INTEGER NUM_KP_KR_NU, NUM_KP_KR_Z, NUM_KP_KR_Q2
      INTEGER NUM_KP_XE_NU, NUM_KP_XE_Z, NUM_KP_XE_Q2
      INTEGER NUM_KM_HE_NU, NUM_KM_HE_Z, NUM_KM_HE_Q2
      INTEGER NUM_KM_NE_NU, NUM_KM_NE_Z, NUM_KM_NE_Q2
      INTEGER NUM_KM_KR_NU, NUM_KM_KR_Z, NUM_KM_KR_Q2
      INTEGER NUM_KM_XE_NU, NUM_KM_XE_Z, NUM_KM_XE_Q2
      INTEGER NUM_KP_HE_PT2, NUM_KM_HE_PT2
      INTEGER NUM_KP_NE_PT2, NUM_KM_NE_PT2
      INTEGER NUM_KP_KR_PT2, NUM_KM_KR_PT2
      INTEGER NUM_KP_XE_PT2, NUM_KM_XE_PT2
      INTEGER NUM_PIP_NU
      INTEGER NUM_PIP_Z
      INTEGER NUM_PIP_Q2
      INTEGER NUM_PIP_PT2
      INTEGER NUM_PI0_NU
      INTEGER NUM_PI0_Z
      INTEGER NUM_PI0_Q2
      INTEGER NUM_PI0_PT2
      INTEGER NUM_PIM_NU
      INTEGER NUM_PIM_Z
      INTEGER NUM_PIM_Q2
      INTEGER NUM_PIM_PT2
C-------JLAB2022
      INTEGER NUM_PIP_2022_C !2023
      INTEGER NUM_PIP_2022_Fe
      INTEGER NUM_PIP_2022_Pb
      INTEGER NUM_PIM_2022_C
      INTEGER NUM_PIM_2022_Fe
      INTEGER NUM_PIM_2022_Pb

      REAL*8 CMS5_NORM
      REAL*8 ATLAS5_Y1_NORM,ATLAS5_Y2_NORM,ATLAS5_Y3_NORM
      REAL*8 CX_STORE(100)

      integer IT,IH,IC,j
      common /meson/ IH,IC

      REAL*8 qTdQcut,phtcut,zcut_up,zcut_low, zcut_1, zcut_2
      REAL*8 pt2cut,zcut, pTdQdZcut, pTdQdzcut_JLAB2022

      real*8 norm,normsumnum,normsumden,ptmin,ptmax

      include "../fit/tools/data-inc.f"

      REAL*8 gamma, g3f, g3D, Nq1, gq1, dq1,
     &       Nq2, gq2, dq2, p_10, p_11, p_12 !2023 
      COMMON /FITP/ gamma, g3f, g3D, Nq1, gq1, dq1,
     &              Nq2, gq2, dq2, p_10, p_11, p_12 !2023 

      real*8 Sep_hermes
      data Sep_hermes/52.7d0/
      real*8 Sep_JLAB2022 !2023
      data Sep_JLAB2022/10.29d0/

      COMMON /NNUM_EXP/ NUM_EXP

      integer iflag,num_exp
      common /iend/ iflag

      !LHC variables
      real*8  s_LHC

      integer cutflag
      real*8 ptcut,etamin,etamax,cutfac
      common /cuts/ ptcut,etamin,etamax,cutflag

      !HERMES variables
      real*8 MULT, NU, Q2, PT2, STAT, SYS, VAL, Z

      real*8 CHI2
      real*8 xb,zh,pht,tmp
      real*8 sigma
      real*8 dcorr,fN
      real*8 s,y,Q,Qmin,Qmax,pt,fuu,Rds
      real*8 fuuA,fuuB, DIS, R_a
      real*8 R_dy
      real*8 xfmin, xfmax, ymin, ymax, Qbar, const, ybar
      real*8 xbmin, xbmax
      real*8 pi
      data pi/3.1415926535d0/
      real*8 M
      data M/0.938d0/
      real*8 CX_LHC,CX_RHIC
      integer i,bool
      real*8 test_pdf, test_ff

      double precision H_AA

*     Input COMMON BLOCKS !2023
      COMMON /HMASS/ H_AA 

C-------QCDNUM initialization
      DOUBLE PRECISION :: array(47), def(-6:6,12), qq(2),wt(2)
      DOUBLE PRECISION :: pdf(-6:6) 
      Double precision as0, r20, xmin, eps 
      integer :: iord, nfin, itype, iset, jset, iosp, nx
      integer :: nxin, nqin, lun, idbug, iqc, iqb, iq0, nq 
      double precision :: q2c, q2b, q0 
      double precision Qf 
      double precision Qf2  

C------------- X grid and mu2 grid paramaters
      data idbug/0/                                     !debug printpout
      data xmin/1.D-4/, nxin/100/                       !x grid, 
      data qq/1.D0,1.D5/, wt/1.D0,1.D0/, nqin/60/       !mu2 grid

C------ Z-array for output file 
      DATA array/0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
     4        0.095, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
     5        0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475,  0.5,
     6        0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7,
     7        0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9,
     8        0.925, 0.95, 0.975, 1.0/

C--------------------------------------------------------------   
      external func                       !input parton dists
      integer, external :: iqfrmq
      data def  /                         !flavor composition
C--   tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
C--   -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6  
     + 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.,         !d
     + 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.,         !u
     + 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0.,         !s
     + 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,         !dbar
     + 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,         !ubar
     + 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,         !sbar
     + 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,         !c
     + 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,         !cbar
     + 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,         !b
     + 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,         !bbar
     + 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,         !t
     + 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.  /       !tbar 
                                       !pdfout

C--   Weight files
      character*26 fnam(3)

      data fnam /'weights/unpolarised.wgt',
     +           'weights/polarised.wgt  ',
     +           'weights/timelike.wgt   '/

      gamma = XX(1) !2023
      g3f   = XX(2)
      g3D   = XX(3)
      Nq1   = XX(4)
      gq1   = XX(5)
      dq1   = XX(6)
      Nq2   = XX(7)
      gq2   = XX(8)
      dq2   = XX(9)
      p_10  = XX(10)
      p_11  = XX(11)
      p_12  = XX(12)

C-----CUTS TO DRELLYAN
      qTdQcut = 0.3d0
      phtcut = 1d0

C-----CUTS TO SIDIS
      pt2cut    = 0.3d0
      zcut      = 0.7d0
      pTdQdzcut = 1000d0
      zcut_up = 0.7d0
      zcut_low = 0.0d0
      zcut_1 = 0.4d0
      zcut_2 = 0.4d0
      pTdQdzcut_JLAB2022 = 0.5d0

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

C-----JLAB2022
      NUM_PIP_2022_C  = 0!2023
      NUM_PIP_2022_Fe = 0
      NUM_PIP_2022_Pb = 0
      NUM_PIM_2022_C  = 0 
      NUM_PIM_2022_Fe = 0 
      NUM_PIM_2022_Pb = 0

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
      CHI2_PI0_XE_PT2   = 0d0

      CHI2_PIP_Z = 0d0
      CHI2_PI0_Z = 0d0
      CHI2_PIM_Z = 0d0

      CHI2_PIP_PT2 = 0d0
      CHI2_PI0_PT2 = 0d0
      CHI2_PIM_PT2 = 0d0

C-----JLAB2022
      CHI2_PIP_2022_C = 0 !2023
      CHI2_PIP_2022_Fe = 0
      CHI2_PIP_2022_Pb = 0
      CHI2_PIM_2022_C = 0
      CHI2_PIM_2022_Fe = 0
      CHI2_PIM_2022_Pb = 0

C-----Make Sure Positive Normalization      
      Call Check(bool)
      if (((1 + Nq1*(1-208**Nq2)).le.0d0) .OR. 
     *      (bool .NE. 1)) then
      print*, "-------------------------------------------"
      chi2 = 10D5
      chisquare = chi2
      print*, "chi2", chi2
      print*, "-------------------------------------------"
      else

C-------------DY (RHIC) p + Au -> muons (Au-going)
      if(I_RHIC_Ratio_pAu1.eq.1) then
      open(unit=  5,file='plot_data/RHIC_Ratio_pAu1.dat')
      write(5 ,*) 'pt ','CX ', 'ERR ', 'FUU'
      do I=1,5
         s       = 200d0
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
           !call DY_overyu(1,s,pt,ymin,ymax,Qmin,Qmax,fuuA,IT)
           !R_dy = fuuA/fuuB
           Call BILINEAR(gamma, g3f, R_dy, I+59)
           tmp = (R_dy - CX_RHIC)**2d0/sigma**2d0
           CHI2_RHIC_Ratio_pAu1 = CHI2_RHIC_Ratio_pAu1 + tmp
           CHI2 = CHI2 + tmp
           num_RHIC_Ratio_pAu1 = num_RHIC_Ratio_pAu1 + 1
           write(5,*) pt, CX_RHIC, sigma, R_dy
        endif
      enddo
      close(5)
      endif

C-------------DY (RHIC) p + Au -> muons (p-going)
      if(I_RHIC_Ratio_pAu2.eq.1) then
      open(unit=  5,file='plot_data/RHIC_Ratio_pAu2.dat')
      write(5 ,*) 'pt ','CX ', 'ERR ', 'FUU'
      do I=1,5
         s       = 200d0
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
           !call DY_overyu(1,s,pt,ymin,ymax,Qmin,Qmax,fuuA,IT)
           !R_dy = fuuA/fuuB
           Call BILINEAR(gamma, g3f, R_dy, I+61)
           tmp = (R_dy - CX_RHIC)**2d0/sigma**2d0
           CHI2_RHIC_Ratio_pAu2 = CHI2_RHIC_Ratio_pAu2 + tmp
           CHI2 = CHI2 + tmp
           num_RHIC_Ratio_pAu2 = num_RHIC_Ratio_pAu2 + 1
           write(5,*) pt, CX_RHIC, sigma, R_dy
        endif
      enddo
      close(5)
      endif

C-----------------------------------------------------------------

      normsumnum = 0d0
      normsumden = 0d0
      cutflag = 1
      etamin = -2.4d0
      etamax =  2.4d0
      ptcut = 20d0
      if(I_CMS5.eq.1) then
      do I=1,13
         s      = 5020D0
         pt       = CMS5_PT(I)
         CX_LHC   = CMS5_CX(I)*1d3
         sigma    = CMS5_ER(I)*1d3
         ymin     = -2.8D0
         ymax     =  2D0
         QBAR = 91.2D0
         dcorr = 0.035
         IT = 208
         if (pt/Qbar.lt.qTdQcut) then
         !call DY_overyu(1,s,pt,ymin,ymax,60d0,120d0,fuuA,IT)
         !FUUA = 2d0*pi*pt*FUUA*208
         Call BILINEAR(gamma, g3f, FUUA, I+16)
         !print *, I
         CX_STORE(I) = FUUA
         normsumnum = normsumnum+CX_LHC*CX_LHC/sigma/sigma
         normsumden = normsumden+CX_LHC*FUUA  /sigma/sigma
         endif
      enddo
      normsumnum = normsumnum+1./dcorr/dcorr
      normsumden = normsumden+1./dcorr/dcorr
      endif
      norm =  normsumnum/normsumden
      fN = norm
      CMS5_NORM = norm

      print *, "CMS5_NORM", CMS5_NORM

      if(I_CMS5.eq.1) then
      open(unit=  5,file='plot_data/CMS5.dat')
      write(5 ,*) 'pt ','CX ', 'ERR ', 'FUU'
      do I=1,13
         s      = 5020D0
         pt       = CMS5_ PT(I)
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

C--------------------------------------------------------------------

      normsumnum = 0d0
      normsumden = 0d0
      if(I_ATLAS5_Y1.eq.1) then
      do I=1,14
         s      = 5020D0
         pt       = ATLAS5_Y1_PT(I)
         CX_LHC   = ATLAS5_Y1_CX(I)*1d3
         sigma    = ATLAS5_Y1_ER(I)*1d3
         dcorr    = 0.027
         ymin     = -3D0
         ymax     =  2D0
         QBAR = 91.2D0
         IT = 208
         if (pt/Qbar.lt.qTdQcut) then
         !call DY_overyu(1,s,pt,ymin,ymax,66d0,116d0,fuuA,IT)
         !FUUA = 2d0*pi*pt*FUUA/PT*208
         Call BILINEAR(gamma, g3f, FUUA, I+24)
         CX_STORE(I) = FUUA
         normsumnum = normsumnum+CX_LHC*CX_LHC/sigma/sigma
         normsumden = normsumden+CX_LHC*FUUA  /sigma/sigma
         endif
      enddo
      normsumnum = normsumnum+1./dcorr/dcorr
      normsumden = normsumden+1./dcorr/dcorr
      endif
      norm =  normsumnum/normsumden
      fN = norm
      ATLAS5_Y1_NORM = norm

      print *, "ATLAS5_Y1_NORM", ATLAS5_Y1_NORM

      if(I_ATLAS5_Y1.eq.1) then
      open(unit=  5,file='plot_data/ATLAS5_Y1.dat')
      write(5 ,*) 'pt ','CX ', 'ERR ', 'FUU'
      do I=1,14
         s      = 5020D0
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

C---------------------------------------------------------------------

      normsumnum = 0d0
      normsumden = 0d0
      if(I_ATLAS5_Y2.eq.1) then
      open(unit=  5,file='plot_data/ATLAS5_Y2.dat')
      write(5 ,*) 'pt ','CX ', 'ERR ', 'FUU'
      do I=1,8
         s      = 5020D0
         pt       = ATLAS5_Y2_PT(I)
         CX_LHC   = ATLAS5_Y2_CX(I)*1d3
         sigma    = ATLAS5_Y2_ER(I)*1d3
         dcorr    = 0.027
         ymin     = -2D0
         ymax     =  0D0
         QBAR = 91.2D0
         IT = 208
         if (pt/Qbar.lt.qTdQcut) then
         call DY_overyu(1,s,pt,ymin,ymax,66d0,116d0,fuuA,IT)
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
         s      = 5020D0
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
         s      = 5020D0
         pt       = ATLAS5_Y3_PT(I)
         CX_LHC   = ATLAS5_Y3_CX(I)*1d3
         sigma    = ATLAS5_Y3_ER(I)*1d3
         dcorr    = 0.027
         ymin     =  0D0
         ymax     =  2D0
         QBAR = 91.2D0
         IT = 208
         if (pt/Qbar.lt.qTdQcut) then
         call DY_overyu(1,s,pt,ymin,ymax,66d0,116d0,fuuA,IT)
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
         s      = 5020D0
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
      s   = dsqrt(1600d0*0.938d0)
      pt    = E772_800_C_pT(I)
      Rds   = E772_800_C_R(I)
      sigma = E772_800_C_Err(I)
      Qmin = 4d0
      Qmax = 9d0
      Qbar = (Qmin+Qmax)/2d0
      xbmin = 0.05d0
      xbmax = 0.3d0
      if (pt/Qbar.lt.qTdQcut) then
      !call DY_overxbQ(4,s,pt,xbmin,xbmax,Qmin,Qmax,fuuB,3)
      !call DY_overxbQ(4,s,pt,xbmin,xbmax,Qmin,Qmax,fuuA,12)
      !R_dy = fuuA/fuuB
      Call BILINEAR(gamma, g3f, R_dy, I)
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
      s   = dsqrt(1600d0*0.938d0)
      pt    = E772_800_Ca_pT(I)
      Rds   = E772_800_Ca_R(I)
      sigma = E772_800_Ca_Err(I)
      Qmin = 4d0
      Qmax = 9d0
      Qbar = (Qmin+Qmax)/2d0
      xbmin = 0.05d0
      xbmax = 0.3d0
      if (pt/Qbar.lt.qTdQcut) then
      !call DY_overxbQ(4,s,pt,xbmin,xbmax,Qmin,Qmax,fuuB,3)
      !call DY_overxbQ(4,s,pt,xbmin,xbmax,Qmin,Qmax,fuuA,40)
      !R_dy = fuuA/fuuB
      Call BILINEAR(gamma, g3f, R_dy, I+4)
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
      s   = dsqrt(1600d0*0.938d0)
      pt    = E772_800_Fe_pT(I)
      Rds   = E772_800_Fe_R(I)
      sigma = E772_800_Fe_Err(I)
      Qmin = 4d0
      Qmax = 9d0
      Qbar = (Qmin+Qmax)/2d0
      xfmin = 0.05d0
      xfmax = 0.3d0
      if (pt/Qbar.lt.qTdQcut) then
      !call DY_overxbQ(4,s,pt,xfmin,xfmax,Qmin,Qmax,fuuB,3)
      !call DY_overxbQ(4,s,pt,xfmin,xfmax,Qmin,Qmax,fuuA,56)
      !R_dy = fuuA/fuuB
      Call BILINEAR(gamma, g3f, R_dy, I+8)
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
      s   = dsqrt(1600d0*0.938d0)
      pt    = E772_800_W_pT(I)
      Rds   = E772_800_W_R(I)
      sigma = E772_800_W_Err(I)
      Qmin = 4d0
      Qmax = 9d0
      Qbar = (Qmin+Qmax)/2d0
      xfmin = 0.05d0
      xfmax = 0.3d0
      if (pt/Qbar.lt.qTdQcut) then
      !call DY_overxbQ(4,s,pt,xfmin,xfmax,Qmin,Qmax,fuuB,3)
      !call DY_overxbQ(4,s,pt,xfmin,xfmax,Qmin,Qmax,fuuA,184)
      !R_dy = fuuA/fuuB
      Call BILINEAR(gamma, g3f, R_dy, I+12)
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
      j = 0
      do I=1,32
      s   = dsqrt(1600d0*0.938d0)
      pt    = E866_800q_pT(I)
      Rds   = E866_800q_R_FeBe(I)
      Q     = E866_800q_Q(I)
      sigma = E866_800q_Err_FeBe(I)
      xfmin = 0.13d0
      xfmax = 0.93d0
      if (pt/Q.lt.qTdQcut) then
      !call DY_overxF(2,s,pt,xfmin,xfmax,Q,fuuB,9)
      !call DY_overxF(2,s,pt,xfmin,xfmax,Q,fuuA,56)
      !R_dy = fuua/fuub
      !print *, 'Fe_I = ', i
      j = j+1
      Call BILINEAR(gamma, g3f, R_dy, j+31) 
      tmp = (R_dy - Rds)**2d0/sigma**2d0
      CHI2_E866_800q = CHI2_E866_800q + tmp
      CHI2 = CHI2 + tmp
      num_E866_800q = num_E866_800q + 1
      write(40,*) pt, Rds, Q, sigma, R_dy
      endif
      enddo
      close(40)
      endif

C------W/Be
      if(I_E866_800q.eq.1) then
      open(unit=  50,file='plot_data/E866_800_Qbin_WBe.dat')
      write(50 ,*) 'pt ','DY-RATIO ','x2 ', 'error ', 'R_dy '
      j = 0
      do I=1,32
      s   = dsqrt(1600d0*0.938d0)
      pt    = E866_800q_pT(I)
      Rds   = E866_800q_R_WBe(I)
      Q     = E866_800q_Q(I)
      sigma = E866_800q_Err_WBe(I)
      xfmin = 0.13d0
      xfmax = 0.93d0
      if (pt/Q.lt.qTdQcut) then
      !call DY_overxF(2,s,pt,xfmin,xfmax,Q,fuuB,9)
      !call DY_overxF(2,s,pt,xfmin,xfmax,Q,fuuA,184)
      !R_dy = fuuA/fuuB
      j = j + 1
      Call BILINEAR(gamma, g3f, R_dy, j+45)
      !print *, 'W_I = ', I
      tmp = (R_dy - Rds)**2d0/sigma**2d0
      CHI2_E866_800q = CHI2_E866_800q + tmp
      CHI2 = CHI2 + tmp
      num_E866_800q = num_E866_800q + 1
      write(50,*) pt, Rds, Q, sigma , R_dy
      endif
      enddo
      close(50)
      endif

C-----z-dependent HERMES
      IF(I_PIP_HE_Z.EQ.1) THEN
      OPEN(UNIT = 5,FILE ='plot_data/HERMES/pip_he_z.dat')
      WRITE(5,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 4
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 4
          IH = 1
          IC = 1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.((z.LT.zcut_1) .OR. (z.GT.zcut_2))
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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

C----- Set-up -----------------------------------------------------------
      H_AA = 20
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 20
          IH = 1
          IC = 1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.((z.LT.zcut_1) .OR. (z.GT.zcut_2))
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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

C----- Set-up -----------------------------------------------------------
      H_AA = 84
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 84
          IH = 1
          IC = 1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.((z.LT.zcut_1) .OR. (z.GT.zcut_2))
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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

C----- Set-up -----------------------------------------------------------
      H_AA = 131
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 131
          IH = 1
          IC = 1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.((z.LT.zcut_1) .OR. (z.GT.zcut_2))
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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

C----- Set-up -----------------------------------------------------------
      H_AA = 4
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 4
          IH = 1
          IC = -1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.((z.LT.zcut_1) .OR. (z.GT.zcut_2))
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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

C----- Set-up -----------------------------------------------------------
      H_AA = 20
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 20
          IH = 1
          IC = -1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.((z.LT.zcut_1) .OR. (z.GT.zcut_2))
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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

C----- Set-up -----------------------------------------------------------
      H_AA = 84
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 84
          IH = 1
          IC = -1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.((z.LT.zcut_1) .OR. (z.GT.zcut_2))
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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

C----- Set-up -----------------------------------------------------------
      H_AA = 131
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 131
          IH = 1
          IC = -1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.((z.LT.zcut_1) .OR. (z.GT.zcut_2))
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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

      IF(I_PI0_HE_Z.EQ.1) THEN
      OPEN(UNIT = 105,FILE ='plot_data/HERMES/pi0_he_z.dat')
      WRITE(105,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 4
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

      DO I=1,9
          VAL  = PI0_HE_Z_VALUE(I)
          MULT = PI0_HE_Z_MULT(I)
          STAT = PI0_HE_Z_STAT(I)
          SYS  = PI0_HE_Z_SYS(I)
          Q2   = PI0_HE_Z_Q2(I)
          Nu   = PI0_HE_Z_NU(I)
          Z    = VAL
          PT2  = PI0_HE_Z_PT2(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          s = Sep_hermes
          IT = 4
          IH = 1
          IC = 0
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.((z.LT.zcut_1) .OR. (z.GT.zcut_2))
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PI0_HE_Z_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PI0_HE_Z = CHI2_PI0_HE_Z + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PI0_HE_Z = NUM_PI0_HE_Z + 1
            WRITE(105,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(105)
      ENDIF

      IF(I_PI0_NE_Z.EQ.1) THEN
      OPEN(UNIT = 106,FILE ='plot_data/HERMES/pi0_ne_z.dat')
      WRITE(106,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 20
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

      DO I=1,9
          VAL  = PI0_NE_Z_VALUE(I)
          MULT = PI0_NE_Z_MULT(I)
          STAT = PI0_NE_Z_STAT(I)
          SYS  = PI0_NE_Z_SYS(I)
          Q2   = PI0_NE_Z_Q2(I)
          Nu   = PI0_NE_Z_NU(I)
          Z    = VAL
          PT2  = PI0_NE_Z_PT2(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          s = Sep_hermes
          IT = 20
          IH = 1
          IC = 0
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.((z.LT.zcut_1) .OR. (z.GT.zcut_2))
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PI0_NE_Z_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PI0_NE_Z = CHI2_PI0_NE_Z + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PI0_NE_Z = NUM_PI0_NE_Z + 1
            WRITE(106,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(106)
      ENDIF

      IF(I_PI0_KR_Z.EQ.1) THEN
      OPEN(UNIT = 107,FILE ='plot_data/HERMES/pi0_kr_z.dat')
      WRITE(107,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 84
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

      DO I=1,9
          VAL  = PI0_KR_Z_VALUE(I)
          MULT = PI0_KR_Z_MULT(I)
          STAT = PI0_KR_Z_STAT(I)
          SYS  = PI0_KR_Z_SYS(I)
          Q2   = PI0_KR_Z_Q2(I)
          Nu   = PI0_KR_Z_NU(I)
          Z    = VAL
          PT2  = PI0_KR_Z_PT2(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          s = Sep_hermes
          IT = 84
          IH = 1
          IC = 0
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.((z.LT.zcut_1) .OR. (z.GT.zcut_2))
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PI0_KR_Z_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PI0_KR_Z = CHI2_PI0_KR_Z + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PI0_KR_Z = NUM_PI0_KR_Z + 1
            WRITE(107,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(107)
      ENDIF

      IF(I_PI0_XE_Z.EQ.1) THEN
      OPEN(UNIT = 108,FILE ='plot_data/HERMES/pi0_xe_z.dat')
      WRITE(108,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 131
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

      DO I=1,9
          VAL  = PI0_XE_Z_VALUE(I)
          MULT = PI0_XE_Z_MULT(I)
          STAT = PI0_XE_Z_STAT(I)
          SYS  = PI0_XE_Z_SYS(I)
          Q2   = PI0_XE_Z_Q2(I)
          Nu   = PI0_XE_Z_NU(I)
          Z    = VAL
          PT2  = PI0_XE_Z_PT2(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          s = Sep_hermes
          IT = 131
          IH = 1
          IC = 0
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.((z.LT.zcut_1) .OR. (z.GT.zcut_2))
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PI0_XE_Z_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PI0_XE_Z = CHI2_PI0_XE_Z + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PI0_XE_Z = NUM_PI0_XE_Z + 1
            WRITE(108,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(108)
      ENDIF


C------------------
C------ PT2 DEPENDET DATA
C------------------
C------ PT2 DEPENDET DATA

C------- PI + DATA
      IF(I_PIP_HE_PT2.EQ.1) THEN
      OPEN(UNIT = 73,
     *STATUS='UNKNOWN',FILE ='plot_data/HERMES/pip_he_pt2.dat')
      WRITE(73,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 4
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 4
          IH = 1
          IC = 1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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
      OPEN(UNIT = 74,
     *STATUS='UNKNOWN',FILE ='plot_data/HERMES/pip_ne_pt2.dat')
      WRITE(74,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 20
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 20
          IH = 1
          IC = 1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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
      OPEN(UNIT = 75,
     *STATUS='UNKNOWN',FILE ='plot_data/HERMES/pip_kr_pt2.dat')
      WRITE(75,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 84
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 84
          IH = 1
          IC = 1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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
      OPEN(UNIT = 76,
     *STATUS='UNKNOWN',FILE ='plot_data/HERMES/pip_xe_pt2.dat')
      WRITE(76,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 131
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 131
          IH = 1
          IC = 1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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
      OPEN(UNIT = 85,
     *STATUS='UNKNOWN',FILE ='plot_data/HERMES/pim_he_pt2.dat')
      WRITE(85,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 4
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 4
          IH = 1
          IC = -1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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
      OPEN(UNIT = 86,
     *STATUS='UNKNOWN',FILE ='plot_data/HERMES/pim_ne_pt2.dat')
      WRITE(86,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
      
C----- Set-up -----------------------------------------------------------
      H_AA = 20
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 20
          IH = 1
          IC = -1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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
      OPEN(UNIT = 87,
     *STATUS='UNKNOWN',FILE ='plot_data/HERMES/pim_kr_pt2.dat')
      WRITE(87,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 84
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 84
          IH = 1
          IC = -1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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
      OPEN(UNIT = 88,
     *STATUS='UNKNOWN',FILE ='plot_data/HERMES/pim_xe_pt2.dat')
      WRITE(88,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 131
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 131
          IH = 1
          IC = -1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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
      OPEN(UNIT = 203,
     *STATUS='UNKNOWN',FILE ='plot_data/HERMES/pi0_he_pt2.dat')
      WRITE(203,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 4
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 4
          IH = 1
          IC = 0
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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
      OPEN(UNIT = 204,
     *STATUS='UNKNOWN',FILE ='plot_data/HERMES/pi0_ne_pt2.dat')
      WRITE(204,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 20
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 20
          IH = 1
          IC = 0
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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
      OPEN(UNIT = 205,
     *STATUS='UNKNOWN',FILE ='plot_data/HERMES/pi0_kr_pt2.dat')
      WRITE(205,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
C----- Set-up -----------------------------------------------------------
      H_AA = 84
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 84
          IH = 1
          IC = 0
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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
      OPEN(UNIT = 206,
     *STATUS='UNKNOWN',FILE ='plot_data/HERMES/pi0_xe_pt2.dat')
      WRITE(206,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 131
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 131
          IH = 1
          IC = 0
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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

c--------------
c JLAB2022
c--------------

      IF(I_pip_2022_C.EQ.1) THEN
      OPEN(UNIT = 207,
     *STATUS='UNKNOWN',FILE ='plot_data/JLAB2022/pi+C_pre.dat')
      WRITE(207,*)  'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 12
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

      DO I=1,48
            MULT = pip2022_c(I)
            STAT = pip2022_cstat(I)
            SYS  = pip2022_csys(I)
            Q2   = pip2022_qc(I)
            xb   = pip2022_xc(I)
            PT2  = 0.5d0*(pip2022_ptlow(I)+pip2022_ptup(I))
            Z    = 0.5d0*(pip2022_zlow(I)+pip2022_zup(I))
            Q    = dsqrt(Q2)
            pht  = sqrt(PT2)
            s  = Sep_JLAB2022
            IT = 12
            IH = 1
            IC = 1
            IF(  (PT2.LT.pt2cut).AND.(z<zcut_up).AND.(z>zcut_low)
     >      .AND.( (pht/(Q*Z)).LT.pTdQdzcut_JLAB2022) ) THEN
                  CALL DISUU(s,Q2,xb,Z,pht,fuua,IT,IH,IC)
                  IT = 3
                  CALL DISUU(s,Q2,xb,Z,pht,fuu,IT,IH,IC)
                  DIS = pip2022_cdis(I)
                  R_a = fuua/(fuu*DIS)
                  sigma = sqrt(STAT**2 + SYS**2)
                  tmp = (R_a - MULT)**2d0/sigma**2d0
                  CHI2_PIP_2022_C = CHI2_PIP_2022_C + tmp
                  CHI2 = CHI2 + tmp
                  NUM_PIP_2022_C = NUM_PIP_2022_C +1
     
                  WRITE(207,*) Z,Q2,PT2,MULT,STAT,SYS,fuu,fuua,DIS
            ENDIF
      ENDDO
      CLOSE(207)
      ENDIF

      IF(I_pip_2022_Fe.EQ.1) THEN
      OPEN(UNIT = 208,
     *STATUS='UNKNOWN',FILE ='plot_data/JLAB2022/pi+Fe_pre.dat')
      WRITE(208,*)  'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 56
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

      DO I=1,48
            MULT = pip2022_fe(I)
            STAT = pip2022_festat(I)
            SYS  = pip2022_fesys(I)
            Q2   = pip2022_qfe(I)
            xb   = pip2022_xfe(I)
            PT2  = 0.5d0*(pip2022_ptlow(I)+pip2022_ptup(I))
            Z    = 0.5d0*(pip2022_zlow(I)+pip2022_zup(I))
            Q    = dsqrt(Q2)
            pht  = sqrt(PT2)
            s  = Sep_JLAB2022
            IT = 56
            IH = 1
            IC = 1
            IF(  (PT2.LT.pt2cut).AND.(z > zcut_low).AND.(z < zcut_up)
     >      .AND.( (pht/(Q*Z)).LT.pTdQdzcut_JLAB2022) ) THEN
                  CALL DISUU(s,Q2,xb,Z,pht,fuua,IT,IH,IC)
                  IT = 3
                  CALL DISUU(s,Q2,xb,Z,pht,fuu,IT,IH,IC)
                  DIS = pip2022_fedis(I)
                  R_a = fuua/(fuu*DIS)
                  sigma = sqrt(STAT**2 + SYS**2)
                  tmp = (R_a - MULT)**2d0/sigma**2d0
                  CHI2_PIP_2022_Fe = CHI2_PIP_2022_Fe + tmp
                  CHI2 = CHI2 + tmp
                  NUM_PIP_2022_Fe = NUM_PIP_2022_Fe +1
              
                  WRITE(208,*) Z,Q2,PT2,MULT,STAT,SYS,fuu,fuua,DIS
            ENDIF
      ENDDO
      CLOSE(208)
      ENDIF

      IF(I_pip_2022_Pb.EQ.1) THEN
      OPEN(UNIT = 209,
     *STATUS='UNKNOWN',FILE ='plot_data/JLAB2022/pi+Pb_pre.dat')
      WRITE(209,*) 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 208
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

      DO I=1,48
            MULT = pip2022_pb(I)
            STAT = pip2022_pbstat(I)
            SYS  = pip2022_pbsys(I)
            Q2   = pip2022_qpb(I)
            xb   = pip2022_xpb(I)
            PT2  = 0.5d0*(pip2022_ptlow(I)+pip2022_ptup(I))
            Z    = 0.5d0*(pip2022_zlow(I)+pip2022_zup(I))
            Q    = dsqrt(Q2)
            pht  = sqrt(PT2)
            s  = Sep_JLAB2022
            IT = 208
            IH = 1
            IC = 1
            IF(  (PT2.LT.pt2cut).AND.(z > zcut_low).AND.(z < zcut_up)
     >      .AND.( (pht/(Q*Z)).LT.pTdQdzcut_JLAB2022) ) THEN
                  CALL DISUU(s,Q2,xb,Z,pht,fuua,IT,IH,IC)
                  IT = 3
                  CALL DISUU(s,Q2,xb,Z,pht,fuu,IT,IH,IC)
                  DIS = pip2022_pbdis(I)
                  R_a = fuua/(fuu*DIS)
                  sigma = sqrt(STAT**2 + SYS**2)
                  tmp = (R_a - MULT)**2d0/sigma**2d0
                  CHI2_PIP_2022_Pb = CHI2_PIP_2022_Pb + tmp
                  CHI2 = CHI2 + tmp
                  NUM_PIP_2022_Pb = NUM_PIP_2022_Pb +1
              
                  WRITE(209,*) Z,Q2,PT2,MULT,STAT,SYS,fuu,fuua,DIS
            ENDIF
      ENDDO
      CLOSE(209)
      ENDIF

      IF(I_pim_2022_C.EQ.1) THEN
      OPEN(UNIT = 210,
     *STATUS='UNKNOWN',FILE ='plot_data/JLAB2022/pi-C_pre.dat')
      WRITE(210,*) 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 12
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

      DO I=1,40
            MULT = pim2022_c(I)
            STAT = pim2022_cstat(I)
            SYS  = pim2022_csys(I)
            Q2   = pim2022_qc(I)
            xb   = pim2022_xc(I)
            PT2  = 0.5d0*(pim2022_ptlow(I)+pim2022_ptup(I))
            Z    = 0.5d0*(pim2022_zlow(I)+pim2022_zup(I))
            Q    = dsqrt(Q2)
            pht  = sqrt(PT2)
            s  = Sep_JLAB2022
            IT = 12
            IH = 1
            IC = -1
            IF(  (PT2.LT.pt2cut).AND.(z > zcut_low).AND.(z < zcut_up)
     >      .AND.( (pht/(Q*Z)).LT.pTdQdzcut_JLAB2022) ) THEN
                  CALL DISUU(s,Q2,xb,Z,pht,fuua,IT,IH,IC)
                  IT = 3
                  CALL DISUU(s,Q2,xb,Z,pht,fuu,IT,IH,IC)
                  DIS = pim2022_cdis(I)
                  R_a = fuua/(fuu*DIS)
                  sigma = sqrt(STAT**2 + SYS**2)
                  tmp = (R_a - MULT)**2d0/sigma**2d0
                  CHI2_PIM_2022_C = CHI2_PIM_2022_C + tmp
                  CHI2 = CHI2 + tmp
                  NUM_PIM_2022_C = NUM_PIM_2022_C +1
              
                  WRITE(210,*) Z,Q2,PT2,MULT,STAT,SYS,fuu,fuua,DIS
            ENDIF
      ENDDO
      CLOSE(210)
      ENDIF

      IF(I_pim_2022_Fe.EQ.1) THEN
      OPEN(UNIT = 211,
     *STATUS='UNKNOWN',FILE ='plot_data/JLAB2022/pi-Fe_pre.dat')
      WRITE(211,*) 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 56
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

      DO I=1,40
            MULT = pim2022_fe(I)
            STAT = pim2022_festat(I)
            SYS  = pim2022_fesys(I)
            Q2   = pim2022_qfe(I)
            xb   = pim2022_xfe(I)
            PT2  = 0.5d0*(pim2022_ptlow(I)+pim2022_ptup(I))
            Z    = 0.5d0*(pim2022_zlow(I)+pim2022_zup(I))
            Q    = dsqrt(Q2)
            pht  = sqrt(PT2)
            s  = Sep_JLAB2022
            IT = 56
            IH = 1
            IC = -1
            IF(  (PT2.LT.pt2cut).AND.(z > zcut_low).AND.(z < zcut_up)
     >      .AND.( (pht/(Q*Z)).LT.pTdQdzcut_JLAB2022) ) THEN
                  CALL DISUU(s,Q2,xb,Z,pht,fuua,IT,IH,IC)
                  IT = 3
                  CALL DISUU(s,Q2,xb,Z,pht,fuu,IT,IH,IC)
                  DIS = pim2022_fedis(I)
                  R_a = fuua/(fuu*DIS)
                  sigma = sqrt(STAT**2 + SYS**2)
                  tmp = (R_a - MULT)**2d0/sigma**2d0
                  CHI2_PIM_2022_Fe = CHI2_PIM_2022_Fe + tmp
                  CHI2 = CHI2 + tmp
                  NUM_PIM_2022_Fe = NUM_PIM_2022_Fe +1
            
                  WRITE(211,*) Z,Q2,PT2,MULT,STAT,SYS,fuu,fuua,DIS
            ENDIF
      ENDDO
      CLOSE(211)
      ENDIF

      IF(I_pim_2022_Pb.EQ.1) THEN
      OPEN(UNIT = 212,
     *STATUS='UNKNOWN',FILE ='plot_data/JLAB2022/pi-Pb_pre.dat')
      WRITE(212,*) 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 208
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

      DO I=1,40
            MULT = pim2022_pb(I)
            STAT = pim2022_pbstat(I)
            SYS  = pim2022_pbsys(I)
            Q2   = pim2022_qpb(I)
            xb   = pim2022_xpb(I)
            PT2  = 0.5d0*(pim2022_ptlow(I)+pim2022_ptup(I))
            Z    = 0.5d0*(pim2022_zlow(I)+pim2022_zup(I))
            Q    = dsqrt(Q2)
            pht  = sqrt(PT2)
            s  = Sep_JLAB2022
            IT = 208
            IH = 1
            IC = -1
            IF(  (PT2.LT.pt2cut).AND.(z > zcut_low).AND.(z < zcut_up)
     >      .AND.( (pht/(Q*Z)).LT.pTdQdzcut_JLAB2022) ) THEN
                  CALL DISUU(s,Q2,xb,Z,pht,fuua,IT,IH,IC)
                  IT = 3
                  CALL DISUU(s,Q2,xb,Z,pht,fuu,IT,IH,IC)
                  DIS = pim2022_pbdis(I)
                  R_a = fuua/(fuu*DIS)
                  sigma = sqrt(STAT**2 + SYS**2)
                  tmp = (R_a - MULT)**2d0/sigma**2d0
                  CHI2_PIM_2022_Pb = CHI2_PIM_2022_Pb + tmp
                  CHI2 = CHI2 + tmp
                  NUM_PIM_2022_Pb = NUM_PIM_2022_Pb +1
          
                  WRITE(212,*) Z,Q2,PT2,MULT,STAT,SYS,fuu,fuua,DIS
            ENDIF
      ENDDO
      CLOSE(212)
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
     >       +NUM_PIP_2022_C  !2023
     >       +NUM_PIP_2022_Fe
     >       +NUM_PIP_2022_Pb
     >       +NUM_PIM_2022_C
     >       +NUM_PIM_2022_Fe
     >       +NUM_PIM_2022_Pb

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

      write(6,110) 'chi2:',chisquare
      write(6,112) 'total num of data in fit:',num_exp
      write(6,112) 'number of parameters:',nfit1
      write(6,110) 'chi2/dof:',chisquare/(num_exp-nfit1)
      WRITE(6,  *) '-------------------------'
      write(6,113) 'fitted parameters:'
      write(6,*) 'gamma, g3f, g3D, Nq1, gq1, dq1, Nq2, dq2, dq2, p_10, 
     & p_11, p_12' !2023
      write(6,*) xx(1), xx(2), xx(3), xx(4), xx(5), xx(6), 
     &           xx(7), xx(8), xx(9), xx(10), xx(11), xx(12) !2023
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

      write(6,*) 'JLAB2022 (PIP-C):',CHI2_PIP_2022_C  ,NUM_PIP_2022_C !2023
      write(6,*) 'JLAB2022 (PIP-Fe):',CHI2_PIP_2022_Fe  ,NUM_PIP_2022_Fe
      write(6,*) 'JLAB2022 (PIP-Pb):',CHI2_PIP_2022_Pb  ,NUM_PIP_2022_Pb
      write(6,*) 'JLAB2022 (PIM-C):',CHI2_PIM_2022_C  ,NUM_PIM_2022_C
      write(6,*) 'JLAB2022 (PIM-Fe):',CHI2_PIM_2022_Fe  ,NUM_PIM_2022_Fe
      write(6,*) 'JLAB2022 (PIM-Pb):',CHI2_PIM_2022_Pb  ,NUM_PIM_2022_Pb

c-----write out the final result on a file-------
      if(iflag.eq.3) then

         open(unit=82,file='result.txt',status='unknown')
         write(82,110) 'chi2:',chisquare
         write(82,112) 'total num of data in fit:',num_exp
         write(82,112) 'number of parameters:',nfit1
         write(82,110) 'chi2/dof:',chisquare/(num_exp-nfit1)
         WRITE(82,  *) '-------------------------'
      write(82,113) 'fitted parameters:'
      write(82,*) 'gamma, g3f, g3D, Nq1, gq1, dq1, Nq2, dq2, dq2, p_10, 
     & p_11, p_12' !2023
      write(82,*) xx(1), xx(2), xx(3), xx(4), xx(5), xx(6), 
     &           xx(7), xx(8), xx(9), xx(10), xx(11), xx(12) !2023
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
         write(82,*) 'JLAB2022 (PIP-C):',CHI2_PIP_2022_C  ,
     *      NUM_PIP_2022_C !2023
         write(82,*) 'JLAB2022 (PIP-Fe):',CHI2_PIP_2022_Fe  ,
     *      NUM_PIP_2022_Fe
         write(82,*) 'JLAB2022 (PIP-Pb):',CHI2_PIP_2022_Pb  ,
     *      NUM_PIP_2022_Pb
         write(82,*) 'JLAB2022 (PIM-C):',CHI2_PIM_2022_C  ,
     *      NUM_PIM_2022_C
         write(82,*) 'JLAB2022 (PIM-Fe):',CHI2_PIM_2022_Fe  ,
     *      NUM_PIM_2022_Fe
         write(82,*) 'JLAB2022 (PIM-Pb):',CHI2_PIM_2022_Pb  ,
     *      NUM_PIM_2022_Pb
         write(82,*)
         close(82)

         open(unit=1000,file='chi2.dat',status='unknown')
         write(1000,*) 'chi2/dof ', 'RHIC-Au ',
     *                 'RHIC-p ',
     *                 'ATLAS5-Y1 ', 'CMS5 ', 'E866 ', 'E772 ',
     *                 'HERMES ',
     *                 'PIP-C ', 'PIP-Fe ', 'PIP-Pb ',
     *                 'PIM-C ', 'PIM-Fe ', 'PIM-Pb '
         write(1000,*) chisquare/(num_exp-nfit1),
     *                 CHI2_RHIC_Ratio_pAu1,
     *                 CHI2_RHIC_Ratio_pAu2,
     *                 CHI2_ATLAS5_Y1, CHI2_CMS5,
     *                 CHI2_E866_800q, CHI2_E772_800,
     *                 CHI2_HERMES,
     *                 CHI2_PIP_2022_C, CHI2_PIP_2022_Fe,
     *                 CHI2_PIP_2022_Pb, CHI2_PIM_2022_C,
     *                 CHI2_PIM_2022_Fe, CHI2_PIM_2022_Pb
         close(1000)

      endif

      endif

      return
      end


      function chisquare_calculate(XX)
      IMPLICIT NONE
      integer nfit
      common /params/ nfit
      real*8 xx(nfit)
      real*8 chisquare
      real*8 chisquare_calculate
      real*8 CHI2_E866_800q, CHI2_E772_800
      REAL*8 CHI2_HERMES
      REAL*8 CHI2_PIP_HE_NU, CHI2_PIP_HE_Z, CHI2_PIP_HE_Q2
      REAL*8 CHI2_PIP_NE_NU, CHI2_PIP_NE_Z, CHI2_PIP_NE_Q2
      REAL*8 CHI2_PIP_KR_NU, CHI2_PIP_KR_Z, CHI2_PIP_KR_Q2
      REAL*8 CHI2_PIP_XE_NU, CHI2_PIP_XE_Z, CHI2_PIP_XE_Q2
      REAL*8 CHI2_PIM_HE_NU, CHI2_PIM_HE_Z, CHI2_PIM_HE_Q2
      REAL*8 CHI2_PIM_NE_NU, CHI2_PIM_NE_Z, CHI2_PIM_NE_Q2
      REAL*8 CHI2_PIM_KR_NU, CHI2_PIM_KR_Z, CHI2_PIM_KR_Q2
      REAL*8 CHI2_PIM_XE_NU, CHI2_PIM_XE_Z, CHI2_PIM_XE_Q2
      REAL*8 CHI2_PI0_HE_NU, CHI2_PI0_HE_Z, CHI2_PI0_HE_Q2
      REAL*8 CHI2_PI0_NE_NU, CHI2_PI0_NE_Z, CHI2_PI0_NE_Q2
      REAL*8 CHI2_PI0_KR_NU, CHI2_PI0_KR_Z, CHI2_PI0_KR_Q2
      REAL*8 CHI2_PI0_XE_NU, CHI2_PI0_XE_Z, CHI2_PI0_XE_Q2
      REAL*8 CHI2_PIP_HE_PT2, CHI2_PIM_HE_PT2, CHI2_PI0_HE_PT2
      REAL*8 CHI2_PIP_NE_PT2, CHI2_PIM_NE_PT2, CHI2_PI0_NE_PT2
      REAL*8 CHI2_PIP_KR_PT2, CHI2_PIM_KR_PT2, CHI2_PI0_KR_PT2
      REAL*8 CHI2_PIP_XE_PT2, CHI2_PIM_XE_PT2, CHI2_PI0_XE_PT2
      REAL*8 CHI2_KP_HE_NU, CHI2_KP_HE_Z, CHI2_KP_HE_Q2
      REAL*8 CHI2_KP_NE_NU, CHI2_KP_NE_Z, CHI2_KP_NE_Q2
      REAL*8 CHI2_KP_KR_NU, CHI2_KP_KR_Z, CHI2_KP_KR_Q2
      REAL*8 CHI2_KP_XE_NU, CHI2_KP_XE_Z, CHI2_KP_XE_Q2
      REAL*8 CHI2_KM_HE_NU, CHI2_KM_HE_Z, CHI2_KM_HE_Q2
      REAL*8 CHI2_KM_NE_NU, CHI2_KM_NE_Z, CHI2_KM_NE_Q2
      REAL*8 CHI2_KM_KR_NU, CHI2_KM_KR_Z, CHI2_KM_KR_Q2
      REAL*8 CHI2_KM_XE_NU, CHI2_KM_XE_Z, CHI2_KM_XE_Q2
      REAL*8 CHI2_KP_HE_PT2, CHI2_KM_HE_PT2
      REAL*8 CHI2_KP_NE_PT2, CHI2_KM_NE_PT2
      REAL*8 CHI2_KP_KR_PT2, CHI2_KM_KR_PT2
      REAL*8 CHI2_KP_XE_PT2, CHI2_KM_XE_PT2
      REAL*8 CHI2_PIP_NU
      REAL*8 CHI2_PIP_Z
      REAL*8 CHI2_PIP_Q2
      REAL*8 CHI2_PIP_PT2
      REAL*8 CHI2_PI0_NU
      REAL*8 CHI2_PI0_Z
      REAL*8 CHI2_PI0_Q2
      REAL*8 CHI2_PI0_PT2
      REAL*8 CHI2_PIM_NU
      REAL*8 CHI2_PIM_Z
      REAL*8 CHI2_PIM_Q2
      REAL*8 CHI2_PIM_PT2
      REAL*8 CHI2_PIP_2022_C !2023
      REAL*8 CHI2_PIP_2022_Fe
      REAL*8 CHI2_PIP_2022_Pb
      REAL*8 CHI2_PIM_2022_C
      REAL*8 CHI2_PIM_2022_Fe
      REAL*8 CHI2_PIM_2022_Pb

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
      INTEGER NUM_PIP_HE_NU, NUM_PIP_HE_Z, NUM_PIP_HE_Q2
      INTEGER NUM_PIP_NE_NU, NUM_PIP_NE_Z, NUM_PIP_NE_Q2
      INTEGER NUM_PIP_KR_NU, NUM_PIP_KR_Z, NUM_PIP_KR_Q2
      INTEGER NUM_PIP_XE_NU, NUM_PIP_XE_Z, NUM_PIP_XE_Q2
      INTEGER NUM_PIM_HE_NU, NUM_PIM_HE_Z, NUM_PIM_HE_Q2
      INTEGER NUM_PIM_NE_NU, NUM_PIM_NE_Z, NUM_PIM_NE_Q2
      INTEGER NUM_PIM_KR_NU, NUM_PIM_KR_Z, NUM_PIM_KR_Q2
      INTEGER NUM_PIM_XE_NU, NUM_PIM_XE_Z, NUM_PIM_XE_Q2
      INTEGER NUM_PI0_HE_NU, NUM_PI0_HE_Z, NUM_PI0_HE_Q2
      INTEGER NUM_PI0_NE_NU, NUM_PI0_NE_Z, NUM_PI0_NE_Q2
      INTEGER NUM_PI0_KR_NU, NUM_PI0_KR_Z, NUM_PI0_KR_Q2
      INTEGER NUM_PI0_XE_NU, NUM_PI0_XE_Z, NUM_PI0_XE_Q2
      INTEGER NUM_PIP_HE_PT2, NUM_PIM_HE_PT2, NUM_PI0_HE_PT2
      INTEGER NUM_PIP_NE_PT2, NUM_PIM_NE_PT2, NUM_PI0_NE_PT2
      INTEGER NUM_PIP_KR_PT2, NUM_PIM_KR_PT2, NUM_PI0_KR_PT2
      INTEGER NUM_PIP_XE_PT2, NUM_PIM_XE_PT2, NUM_PI0_XE_PT2
      INTEGER NUM_KP_HE_NU, NUM_KP_HE_Z, NUM_KP_HE_Q2
      INTEGER NUM_KP_NE_NU, NUM_KP_NE_Z, NUM_KP_NE_Q2
      INTEGER NUM_KP_KR_NU, NUM_KP_KR_Z, NUM_KP_KR_Q2
      INTEGER NUM_KP_XE_NU, NUM_KP_XE_Z, NUM_KP_XE_Q2
      INTEGER NUM_KM_HE_NU, NUM_KM_HE_Z, NUM_KM_HE_Q2
      INTEGER NUM_KM_NE_NU, NUM_KM_NE_Z, NUM_KM_NE_Q2
      INTEGER NUM_KM_KR_NU, NUM_KM_KR_Z, NUM_KM_KR_Q2
      INTEGER NUM_KM_XE_NU, NUM_KM_XE_Z, NUM_KM_XE_Q2
      INTEGER NUM_KP_HE_PT2, NUM_KM_HE_PT2
      INTEGER NUM_KP_NE_PT2, NUM_KM_NE_PT2
      INTEGER NUM_KP_KR_PT2, NUM_KM_KR_PT2
      INTEGER NUM_KP_XE_PT2, NUM_KM_XE_PT2
      INTEGER NUM_PIP_NU
      INTEGER NUM_PIP_Z
      INTEGER NUM_PIP_Q2
      INTEGER NUM_PIP_PT2
      INTEGER NUM_PI0_NU
      INTEGER NUM_PI0_Z
      INTEGER NUM_PI0_Q2
      INTEGER NUM_PI0_PT2
      INTEGER NUM_PIM_NU
      INTEGER NUM_PIM_Z
      INTEGER NUM_PIM_Q2
      INTEGER NUM_PIM_PT2

C-------JLAB2022
      INTEGER NUM_PIP_2022_C !2023
      INTEGER NUM_PIP_2022_Fe
      INTEGER NUM_PIP_2022_Pb
      INTEGER NUM_PIM_2022_C
      INTEGER NUM_PIM_2022_Fe
      INTEGER NUM_PIM_2022_Pb

      REAL*8 CMS5_NORM
      REAL*8 ATLAS5_Y1_NORM,ATLAS5_Y2_NORM,ATLAS5_Y3_NORM
      REAL*8 CX_STORE(100)

      integer IT,IH,IC
      common /meson/ IH,IC

      REAL*8 qTdQcut,phtcut,zcut_low,zcut_up,zcut_1,zcut_2
      REAL*8 pt2cut,zcut,pTdQdZcut,pTdQdZcut_JLAB2022

      real*8 norm,normsumnum,normsumden,ptmin,ptmax

      include "../fit/tools/data-inc.f"

      REAL*8 gamma, g3f, g3D, Nq1, gq1, dq1,
     &       Nq2, gq2, dq2, p_10, p_11, p_12 !2023 
      COMMON /FITP/ gamma, g3f, g3D, Nq1, gq1, dq1,
     &              Nq2, gq2, dq2, p_10, p_11, p_12 !2023

      real*8 Sep_hermes
      data Sep_hermes/52.7d0/
      real*8 Sep_JLAB2022 !2023
      data Sep_JLAB2022/10.29d0/

      COMMON /NNUM_EXP/ NUM_EXP

      integer iflag,num_exp
      common /iend/ iflag

      !LHC variables
      real*8  s_LHC

      integer cutflag
      real*8 ptcut,etamin,etamax,cutfac
      common /cuts/ ptcut,etamin,etamax,cutflag

      !HERMES variables
      real*8 MULT, NU, Q2, PT2, STAT, SYS, VAL, Z

      real*8 CHI2
      real*8 xb,zh,pht,tmp
      real*8 sigma
      real*8 dcorr,fN
      real*8 s,y,Q,Qmin,Qmax,pt,fuu,Rds
      real*8 fuuA,fuuB, DIS, R_a
      real*8 R_dy
      real*8 xfmin, xfmax, ymin, ymax, Qbar, const, ybar
      real*8 xbmin, xbmax
      real*8 pi
      data pi/3.1415926535d0/
      real*8 M
      data M/0.938d0/
      real*8 CX_LHC,CX_RHIC
      integer i, j, bool

      double precision H_AA

*     Input COMMON BLOCKS !2023
      COMMON /HMASS/ H_AA 

C-------QCDNUM initialization
      DOUBLE PRECISION :: array(47), def(-6:6,12), qq(2),wt(2)
      DOUBLE PRECISION :: pdf(-6:6) 
      Double precision as0, r20, xmin, eps 
      integer :: iord, nfin, itype, iset, jset, iosp, nx
      integer :: nxin, nqin, lun, idbug, iqc, iqb, iq0, nq 
      double precision :: q2c, q2b, q0 
      double precision Qf 
      double precision Qf2  

C------------- X grid and mu2 grid paramaters
      data idbug/0/                                     !debug printpout
      data xmin/1.D-4/, nxin/100/                       !x grid, 
      data qq/1.D0,1.D5/, wt/1.D0,1.D0/, nqin/60/       !mu2 grid

C------ Z-array for output file 
      DATA array/0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
     4        0.095, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
     5        0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475,  0.5,
     6        0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7,
     7        0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9,
     8        0.925, 0.95, 0.975, 1.0/

C--------------------------------------------------------------   
      external func                       !input parton dists
      integer, external :: iqfrmq
      data def  /                         !flavor composition
C--   tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
C--   -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6  
     + 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.,         !d
     + 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.,         !u
     + 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0.,         !s
     + 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,         !dbar
     + 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,         !ubar
     + 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,         !sbar
     + 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,         !c
     + 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,         !cbar
     + 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,         !b
     + 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,         !bbar
     + 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,         !t
     + 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.  /       !tbar 
                                       !pdfout

C--   Weight files
      character*26 fnam(3)

      data fnam /'weights/unpolarised.wgt',
     +           'weights/polarised.wgt  ',
     +           'weights/timelike.wgt   '/

      gamma = XX(1) !2023
      g3f   = XX(2)
      g3D   = XX(3)
      Nq1   = XX(4)
      gq1   = XX(5)
      dq1   = XX(6)
      Nq2   = XX(7)
      gq2   = XX(8)
      dq2   = XX(9)
      p_10  = XX(10)
      p_11  = XX(11)
      p_12  = XX(12)

C-----CUTS TO DRELLYAN
      qTdQcut = 0.3d0
      phtcut = 1d0

C-----CUTS TO SIDIS
      pt2cut    = 0.30d0
      zcut      = 0.7d0
      pTdQdzcut = 1000d0
      zcut_low = 0.0d0
      zcut_up = 0.7d0
      pTdQdZcut_JLAB2022 = 0.5
      zcut_1 = 0.4d0
      zcut_2 = 0.4d0

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

C-----JLAB2022
      NUM_PIP_2022_C  = 0!2023
      NUM_PIP_2022_Fe = 0
      NUM_PIP_2022_Pb = 0
      NUM_PIM_2022_C  = 0
      NUM_PIM_2022_Fe = 0
      NUM_PIM_2022_Pb = 0
C------------------

      CHI2 = 0

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
      CHI2_PI0_XE_PT2   = 0d0

      CHI2_PIP_Z = 0d0
      CHI2_PI0_Z = 0d0
      CHI2_PIM_Z = 0d0

      CHI2_PIP_PT2 = 0d0
      CHI2_PI0_PT2 = 0d0
      CHI2_PIM_PT2 = 0d0

C-----JLAB2022
      CHI2_PIP_2022_C = 0 !2023
      CHI2_PIP_2022_Fe = 0
      CHI2_PIP_2022_Pb = 0
      CHI2_PIM_2022_C = 0
      CHI2_PIM_2022_Fe = 0
      CHI2_PIM_2022_Pb = 0

C-----Make Sure Positive Normalization      
      Call Check(bool)
      if (((1 + Nq1*(1-208**Nq2)).le.0d0) .OR. 
     *      (bool .NE. 1)) then
      print*, "-------------------------------------------"
      chi2 = 10D5
      chisquare = chi2
      print*, "chi2", chi2
      print*, "-------------------------------------------"
      else

C-------------DY (RHIC) p + Au -> muons (Au-going)
      if(I_RHIC_Ratio_pAu1.eq.1) then
      do I=1,5
         s       = 200d0
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
           !call DY_overyu(1,s,pt,ymin,ymax,Qmin,Qmax,fuuA,IT)
           !R_dy = fuuA/fuuB
           Call BILINEAR(gamma, g3f, R_dy, I+59)
           tmp = (R_dy - CX_RHIC)**2d0/sigma**2d0
           CHI2_RHIC_Ratio_pAu1 = CHI2_RHIC_Ratio_pAu1 + tmp
           CHI2 = CHI2 + tmp
           num_RHIC_Ratio_pAu1 = num_RHIC_Ratio_pAu1 + 1
        endif
      enddo
      endif

C-------------DY (RHIC) p + Au -> muons (p-going)
      if(I_RHIC_Ratio_pAu2.eq.1) then
      do I=1,5
         s       = 200d0
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
           !call DY_overyu(1,s,pt,ymin,ymax,Qmin,Qmax,fuuA,IT)
           !R_dy = fuuA/fuuB
           Call BILINEAR(gamma, g3f, R_dy, I+61)
           tmp = (R_dy - CX_RHIC)**2d0/sigma**2d0
           CHI2_RHIC_Ratio_pAu2 = CHI2_RHIC_Ratio_pAu2 + tmp
           CHI2 = CHI2 + tmp
           num_RHIC_Ratio_pAu2 = num_RHIC_Ratio_pAu2 + 1
        endif
      enddo
      close(5)
      endif

C-----------------------------------------------------------------

      normsumnum = 0d0
      normsumden = 0d0
      cutflag = 1
      etamin = -2.4d0
      etamax =  2.4d0
      ptcut = 20d0
      if(I_CMS5.eq.1) then
      do I=1,13
         s      = 5020D0
         pt       = CMS5_PT(I)
         CX_LHC   = CMS5_CX(I)*1d3
         sigma    = CMS5_ER(I)*1d3
         ymin     = -2.8D0
         ymax     =  2D0
         QBAR = 91.2D0
         dcorr = 0.035
         IT = 208
         if (pt/Qbar.lt.qTdQcut) then
         !call DY_overyu(1,s,pt,ymin,ymax,60d0,120d0,fuuA,IT)
         !FUUA = 2d0*pi*pt*FUUA*208
         Call BILINEAR(gamma, g3f, FUUA, I+16)
         !print *, I
         CX_STORE(I) = FUUA
         normsumnum = normsumnum+CX_LHC*CX_LHC/sigma/sigma
         normsumden = normsumden+CX_LHC*FUUA  /sigma/sigma
         endif
      enddo
      normsumnum = normsumnum+1./dcorr/dcorr
      normsumden = normsumden+1./dcorr/dcorr
      endif
      norm =  normsumnum/normsumden
      fN = norm
      CMS5_NORM = norm

      print *, "CMS5_NORM", CMS5_NORM

      if(I_CMS5.eq.1) then
      open(unit=  5,file='plot_data/CMS5.dat')
      write(5 ,*) 'pt ','CX ', 'ERR ', 'FUU'
      do I=1,13
         s      = 5020D0
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

C--------------------------------------------------------------------

      normsumnum = 0d0
      normsumden = 0d0
      if(I_ATLAS5_Y1.eq.1) then
      do I=1,14
         s      = 5020D0
         pt       = ATLAS5_Y1_PT(I)
         CX_LHC   = ATLAS5_Y1_CX(I)*1d3
         sigma    = ATLAS5_Y1_ER(I)*1d3
         dcorr    = 0.027
         ymin     = -3D0
         ymax     =  2D0
         QBAR = 91.2D0
         IT = 208
         if (pt/Qbar.lt.qTdQcut) then
         !call DY_overyu(1,s,pt,ymin,ymax,66d0,116d0,fuuA,IT)
         !FUUA = 2d0*pi*pt*FUUA/PT*208
         Call BILINEAR(gamma, g3f, FUUA, I+24)
         CX_STORE(I) = FUUA
         normsumnum = normsumnum+CX_LHC*CX_LHC/sigma/sigma
         normsumden = normsumden+CX_LHC*FUUA  /sigma/sigma
         endif
      enddo
      normsumnum = normsumnum+1./dcorr/dcorr
      normsumden = normsumden+1./dcorr/dcorr
      endif
      norm =  normsumnum/normsumden
      fN = norm
      ATLAS5_Y1_NORM = norm

      print *, "ATLAS5_Y1_NORM", ATLAS5_Y1_NORM

      if(I_ATLAS5_Y1.eq.1) then
      open(unit=  5,file='plot_data/ATLAS5_Y1.dat')
      write(5 ,*) 'pt ','CX ', 'ERR ', 'FUU'
      do I=1,14
         s      = 5020D0
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

C---------------------------------------------------------------------

      normsumnum = 0d0
      normsumden = 0d0
      if(I_ATLAS5_Y2.eq.1) then
      open(unit=  5,file='plot_data/ATLAS5_Y2.dat')
      write(5 ,*) 'pt ','CX ', 'ERR ', 'FUU'
      do I=1,8
         s      = 5020D0
         pt       = ATLAS5_Y2_PT(I)
         CX_LHC   = ATLAS5_Y2_CX(I)*1d3
         sigma    = ATLAS5_Y2_ER(I)*1d3
         dcorr    = 0.027
         ymin     = -2D0
         ymax     =  0D0
         QBAR = 91.2D0
         IT = 208
         if (pt/Qbar.lt.qTdQcut) then
         call DY_overyu(1,s,pt,ymin,ymax,66d0,116d0,fuuA,IT)
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
         s      = 5020D0
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
         s      = 5020D0
         pt       = ATLAS5_Y3_PT(I)
         CX_LHC   = ATLAS5_Y3_CX(I)*1d3
         sigma    = ATLAS5_Y3_ER(I)*1d3
         dcorr    = 0.027
         ymin     =  0D0
         ymax     =  2D0
         QBAR = 91.2D0
         IT = 208
         if (pt/Qbar.lt.qTdQcut) then
         call DY_overyu(1,s,pt,ymin,ymax,66d0,116d0,fuuA,IT)
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
         s      = 5020D0
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
      s   = dsqrt(1600d0*0.938d0)
      pt    = E772_800_C_pT(I)
      Rds   = E772_800_C_R(I)
      sigma = E772_800_C_Err(I)
      Qmin = 4d0
      Qmax = 9d0
      Qbar = (Qmin+Qmax)/2d0
      xbmin = 0.05d0
      xbmax = 0.3d0
      if (pt/Qbar.lt.qTdQcut) then
      !call DY_overxbQ(4,s,pt,xbmin,xbmax,Qmin,Qmax,fuuB,3)
      !call DY_overxbQ(4,s,pt,xbmin,xbmax,Qmin,Qmax,fuuA,12)
      !R_dy = fuuA/fuuB
      Call BILINEAR(gamma, g3f, R_dy, I)
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
      s   = dsqrt(1600d0*0.938d0)
      pt    = E772_800_Ca_pT(I)
      Rds   = E772_800_Ca_R(I)
      sigma = E772_800_Ca_Err(I)
      Qmin = 4d0
      Qmax = 9d0
      Qbar = (Qmin+Qmax)/2d0
      xbmin = 0.05d0
      xbmax = 0.3d0
      if (pt/Qbar.lt.qTdQcut) then
      !call DY_overxbQ(4,s,pt,xbmin,xbmax,Qmin,Qmax,fuuB,3)
      !call DY_overxbQ(4,s,pt,xbmin,xbmax,Qmin,Qmax,fuuA,40)
      !R_dy = fuuA/fuuB
      Call BILINEAR(gamma, g3f, R_dy, I+4)
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
      s   = dsqrt(1600d0*0.938d0)
      pt    = E772_800_Fe_pT(I)
      Rds   = E772_800_Fe_R(I)
      sigma = E772_800_Fe_Err(I)
      Qmin = 4d0
      Qmax = 9d0
      Qbar = (Qmin+Qmax)/2d0
      xfmin = 0.05d0
      xfmax = 0.3d0
      if (pt/Qbar.lt.qTdQcut) then
      !call DY_overxbQ(4,s,pt,xfmin,xfmax,Qmin,Qmax,fuuB,3)
      !call DY_overxbQ(4,s,pt,xfmin,xfmax,Qmin,Qmax,fuuA,56)
      !R_dy = fuuA/fuuB
      Call BILINEAR(gamma, g3f, R_dy, I+8)
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
      s   = dsqrt(1600d0*0.938d0)
      pt    = E772_800_W_pT(I)
      Rds   = E772_800_W_R(I)
      sigma = E772_800_W_Err(I)
      Qmin = 4d0
      Qmax = 9d0
      Qbar = (Qmin+Qmax)/2d0
      xfmin = 0.05d0
      xfmax = 0.3d0
      if (pt/Qbar.lt.qTdQcut) then
      !call DY_overxbQ(4,s,pt,xfmin,xfmax,Qmin,Qmax,fuuB,3)
      !call DY_overxbQ(4,s,pt,xfmin,xfmax,Qmin,Qmax,fuuA,184)
      !R_dy = fuuA/fuuB
      Call BILINEAR(gamma, g3f, R_dy, I+12)
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
      j = 0
      do I=1,32
      s   = dsqrt(1600d0*0.938d0)
      pt    = E866_800q_pT(I)
      Rds   = E866_800q_R_FeBe(I)
      Q     = E866_800q_Q(I)
      sigma = E866_800q_Err_FeBe(I)
      xfmin = 0.13d0
      xfmax = 0.93d0
      if (pt/Q.lt.qTdQcut) then
      !call DY_overxF(2,s,pt,xfmin,xfmax,Q,fuuB,9)
      !call DY_overxF(2,s,pt,xfmin,xfmax,Q,fuuA,56)
      !R_dy = fuua/fuub
      !print *, 'Fe_I = ', i
      j = j+1
      Call BILINEAR(gamma, g3f, R_dy, j+31) 
      tmp = (R_dy - Rds)**2d0/sigma**2d0
      CHI2_E866_800q = CHI2_E866_800q + tmp
      CHI2 = CHI2 + tmp
      num_E866_800q = num_E866_800q + 1
      endif
      enddo
      endif

C------W/Be
      if(I_E866_800q.eq.1) then
      j = 0
      do I=1,32
      s   = dsqrt(1600d0*0.938d0)
      pt    = E866_800q_pT(I)
      Rds   = E866_800q_R_WBe(I)
      Q     = E866_800q_Q(I)
      sigma = E866_800q_Err_WBe(I)
      xfmin = 0.13d0
      xfmax = 0.93d0
      if (pt/Q.lt.qTdQcut) then
      !call DY_overxF(2,s,pt,xfmin,xfmax,Q,fuuB,9)
      !call DY_overxF(2,s,pt,xfmin,xfmax,Q,fuuA,184)
      !R_dy = fuuA/fuuB
      j = j + 1
      Call BILINEAR(gamma, g3f, R_dy, j+45)
      !print *, 'W_I = ', I
      tmp = (R_dy - Rds)**2d0/sigma**2d0
      CHI2_E866_800q = CHI2_E866_800q + tmp
      CHI2 = CHI2 + tmp
      num_E866_800q = num_E866_800q + 1
      endif
      enddo
      endif

C-----z-dependent HERMES
      IF(I_PIP_HE_Z.EQ.1) THEN
      OPEN(UNIT = 5,FILE ='plot_data/HERMES/pip_he_z.dat')
      WRITE(5,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 4
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 4
          IH = 1
          IC = 1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.((z.LT.zcut_1) .OR. (z.GT.zcut_2))
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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

C----- Set-up -----------------------------------------------------------
      H_AA = 20
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 20
          IH = 1
          IC = 1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.((z.LT.zcut_1) .OR. (z.GT.zcut_2))
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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

C----- Set-up -----------------------------------------------------------
      H_AA = 84
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 84
          IH = 1
          IC = 1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.((z.LT.zcut_1) .OR. (z.GT.zcut_2))
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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

C----- Set-up -----------------------------------------------------------
      H_AA = 131
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 131
          IH = 1
          IC = 1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.((z.LT.zcut_1) .OR. (z.GT.zcut_2))
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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

C----- Set-up -----------------------------------------------------------
      H_AA = 4
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 4
          IH = 1
          IC = -1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.((z.LT.zcut_1) .OR. (z.GT.zcut_2))
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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

C----- Set-up -----------------------------------------------------------
      H_AA = 20
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 20
          IH = 1
          IC = -1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.((z.LT.zcut_1) .OR. (z.GT.zcut_2))
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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

C----- Set-up -----------------------------------------------------------
      H_AA = 84
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 84
          IH = 1
          IC = -1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.((z.LT.zcut_1) .OR. (z.GT.zcut_2))
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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

C----- Set-up -----------------------------------------------------------
      H_AA = 131
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 131
          IH = 1
          IC = -1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.((z.LT.zcut_1) .OR. (z.GT.zcut_2))
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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

      IF(I_PI0_HE_Z.EQ.1) THEN
      OPEN(UNIT = 105,FILE ='plot_data/HERMES/pi0_he_z.dat')
      WRITE(105,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 4
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

      DO I=1,9
          VAL  = PI0_HE_Z_VALUE(I)
          MULT = PI0_HE_Z_MULT(I)
          STAT = PI0_HE_Z_STAT(I)
          SYS  = PI0_HE_Z_SYS(I)
          Q2   = PI0_HE_Z_Q2(I)
          Nu   = PI0_HE_Z_NU(I)
          Z    = VAL
          PT2  = PI0_HE_Z_PT2(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          s = Sep_hermes
          IT = 4
          IH = 1
          IC = 0
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.((z.LT.zcut_1) .OR. (z.GT.zcut_2))
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PI0_HE_Z_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PI0_HE_Z = CHI2_PI0_HE_Z + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PI0_HE_Z = NUM_PI0_HE_Z + 1
            WRITE(105,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(105)
      ENDIF

      IF(I_PI0_NE_Z.EQ.1) THEN
      OPEN(UNIT = 106,FILE ='plot_data/HERMES/pi0_ne_z.dat')
      WRITE(106,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 20
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

      DO I=1,9
          VAL  = PI0_NE_Z_VALUE(I)
          MULT = PI0_NE_Z_MULT(I)
          STAT = PI0_NE_Z_STAT(I)
          SYS  = PI0_NE_Z_SYS(I)
          Q2   = PI0_NE_Z_Q2(I)
          Nu   = PI0_NE_Z_NU(I)
          Z    = VAL
          PT2  = PI0_NE_Z_PT2(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          s = Sep_hermes
          IT = 20
          IH = 1
          IC = 0
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.((z.LT.zcut_1) .OR. (z.GT.zcut_2))
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PI0_NE_Z_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PI0_NE_Z = CHI2_PI0_NE_Z + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PI0_NE_Z = NUM_PI0_NE_Z + 1
            WRITE(106,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(106)
      ENDIF

      IF(I_PI0_KR_Z.EQ.1) THEN
      OPEN(UNIT = 107,FILE ='plot_data/HERMES/pi0_kr_z.dat')
      WRITE(107,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 84
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

      DO I=1,9
          VAL  = PI0_KR_Z_VALUE(I)
          MULT = PI0_KR_Z_MULT(I)
          STAT = PI0_KR_Z_STAT(I)
          SYS  = PI0_KR_Z_SYS(I)
          Q2   = PI0_KR_Z_Q2(I)
          Nu   = PI0_KR_Z_NU(I)
          Z    = VAL
          PT2  = PI0_KR_Z_PT2(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          s = Sep_hermes
          IT = 84
          IH = 1
          IC = 0
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.((z.LT.zcut_1) .OR. (z.GT.zcut_2))
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PI0_KR_Z_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PI0_KR_Z = CHI2_PI0_KR_Z + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PI0_KR_Z = NUM_PI0_KR_Z + 1
            WRITE(107,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(107)
      ENDIF

      IF(I_PI0_XE_Z.EQ.1) THEN
      OPEN(UNIT = 108,FILE ='plot_data/HERMES/pi0_xe_z.dat')
      WRITE(108,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 131
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

      DO I=1,9
          VAL  = PI0_XE_Z_VALUE(I)
          MULT = PI0_XE_Z_MULT(I)
          STAT = PI0_XE_Z_STAT(I)
          SYS  = PI0_XE_Z_SYS(I)
          Q2   = PI0_XE_Z_Q2(I)
          Nu   = PI0_XE_Z_NU(I)
          Z    = VAL
          PT2  = PI0_XE_Z_PT2(I)
          Q    = dsqrt(Q2)
          xb   = Q2/(2*M*Nu)
          pht = sqrt(PT2)
          s = Sep_hermes
          IT = 131
          IH = 1
          IC = 0
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.((z.LT.zcut_1) .OR. (z.GT.zcut_2))
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
            DIS = PI0_XE_Z_DIS(I)
            R_a = fuua/(fuu*DIS)
            sigma = sqrt(STAT**2 + SYS**2)
            tmp = (R_a - MULT)**2d0/sigma**2d0
            CHI2_PI0_XE_Z = CHI2_PI0_XE_Z + tmp
            CHI2 = CHI2 + tmp
            CHI2_HERMES = CHI2_HERMES + tmp
            NUM_PI0_XE_Z = NUM_PI0_XE_Z + 1
            WRITE(108,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
          ENDIF
      ENDDO
      CLOSE(108)
      ENDIF


C------------------
C------ PT2 DEPENDET DATA

C------- PI + DATA
      IF(I_PIP_HE_PT2.EQ.1) THEN
      OPEN(UNIT = 73,
     *STATUS='UNKNOWN',FILE ='plot_data/HERMES/pip_he_pt2.dat')
      WRITE(73,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 4
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 4
          IH = 1
          IC = 1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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
      OPEN(UNIT = 74,
     *STATUS='UNKNOWN',FILE ='plot_data/HERMES/pip_ne_pt2.dat')
      WRITE(74,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 20
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 20
          IH = 1
          IC = 1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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
      OPEN(UNIT = 75,
     *STATUS='UNKNOWN',FILE ='plot_data/HERMES/pip_kr_pt2.dat')
      WRITE(75,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 84
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 84
          IH = 1
          IC = 1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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
      OPEN(UNIT = 76,
     *STATUS='UNKNOWN',FILE ='plot_data/HERMES/pip_xe_pt2.dat')
      WRITE(76,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 131
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 131
          IH = 1
          IC = 1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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
      OPEN(UNIT = 85,
     *STATUS='UNKNOWN',FILE ='plot_data/HERMES/pim_he_pt2.dat')
      WRITE(85,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 4
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 4
          IH = 1
          IC = -1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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
      OPEN(UNIT = 86,
     *STATUS='UNKNOWN',FILE ='plot_data/HERMES/pim_ne_pt2.dat')
      WRITE(86,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
      
C----- Set-up -----------------------------------------------------------
      H_AA = 20
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 20
          IH = 1
          IC = -1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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
      OPEN(UNIT = 87,
     *STATUS='UNKNOWN',FILE ='plot_data/HERMES/pim_kr_pt2.dat')
      WRITE(87,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 84
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 84
          IH = 1
          IC = -1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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
      OPEN(UNIT = 88,
     *STATUS='UNKNOWN',FILE ='plot_data/HERMES/pim_xe_pt2.dat')
      WRITE(88,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 131
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 131
          IH = 1
          IC = -1
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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
      OPEN(UNIT = 203,
     *STATUS='UNKNOWN',FILE ='plot_data/HERMES/pi0_he_pt2.dat')
      WRITE(203,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 4
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 4
          IH = 1
          IC = 0
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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
      OPEN(UNIT = 204,
     *STATUS='UNKNOWN',FILE ='plot_data/HERMES/pi0_ne_pt2.dat')
      WRITE(204,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 20
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 20
          IH = 1
          IC = 0
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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
      OPEN(UNIT = 205,
     *STATUS='UNKNOWN',FILE ='plot_data/HERMES/pi0_kr_pt2.dat')
      WRITE(205,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
C----- Set-up -----------------------------------------------------------
      H_AA = 84
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 84
          IH = 1
          IC = 0
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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
      OPEN(UNIT = 206,
     *STATUS='UNKNOWN',FILE ='plot_data/HERMES/pi0_xe_pt2.dat')
      WRITE(206,*)  'Nu','Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 131
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

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
          s = Sep_hermes
          IT = 131
          IH = 1
          IC = 0
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(s,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(s,Q2,xb,z,pht,fuu,IT,IH,IC)
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

c--------------
c JLAB2022
c--------------

      IF(I_pip_2022_C.EQ.1) THEN
      OPEN(UNIT = 207,
     *STATUS='UNKNOWN',FILE ='plot_data/JLAB2022/pi+C_pre.dat')
      WRITE(207,*)  'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 12
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

      DO I=1,48
            MULT = pip2022_c(I)
            STAT = pip2022_cstat(I)
            SYS  = pip2022_csys(I)
            Q2   = pip2022_qc(I)
            xb   = pip2022_xc(I)
            PT2  = 0.5d0*(pip2022_ptlow(I)+pip2022_ptup(I))
            Z    = 0.5d0*(pip2022_zlow(I)+pip2022_zup(I))
            Q    = dsqrt(Q2)
            pht  = sqrt(PT2)
            s  = Sep_JLAB2022
            IT = 12
            IH = 1
            IC = 1
            IF(  (PT2.LT.pt2cut).AND.(z<zcut_up).AND.(z>zcut_low)
     >      .AND.( (pht/(Q*Z)).LT.pTdQdzcut_JLAB2022) ) THEN
                  CALL DISUU(s,Q2,xb,Z,pht,fuua,IT,IH,IC)
                  IT = 3
                  CALL DISUU(s,Q2,xb,Z,pht,fuu,IT,IH,IC)
                  DIS = pip2022_cdis(I)
                  R_a = fuua/(fuu*DIS)
                  sigma = sqrt(STAT**2 + SYS**2)
                  tmp = (R_a - MULT)**2d0/sigma**2d0
                  CHI2_PIP_2022_C = CHI2_PIP_2022_C + tmp
                  CHI2 = CHI2 + tmp
                  NUM_PIP_2022_C = NUM_PIP_2022_C +1
     
                  WRITE(207,*) Z,Q2,PT2,MULT,STAT,SYS,fuu,fuua,DIS
            ENDIF
      ENDDO
      CLOSE(207)
      ENDIF

      IF(I_pip_2022_Fe.EQ.1) THEN
      OPEN(UNIT = 208,
     *STATUS='UNKNOWN',FILE ='plot_data/JLAB2022/pi+Fe_pre.dat')
      WRITE(208,*)  'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 56
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

      DO I=1,48
            MULT = pip2022_fe(I)
            STAT = pip2022_festat(I)
            SYS  = pip2022_fesys(I)
            Q2   = pip2022_qfe(I)
            xb   = pip2022_xfe(I)
            PT2  = 0.5d0*(pip2022_ptlow(I)+pip2022_ptup(I))
            Z    = 0.5d0*(pip2022_zlow(I)+pip2022_zup(I))
            Q    = dsqrt(Q2)
            pht  = sqrt(PT2)
            s  = Sep_JLAB2022
            IT = 56
            IH = 1
            IC = 1
            IF(  (PT2.LT.pt2cut).AND.(z > zcut_low).AND.(z < zcut_up)
     >      .AND.( (pht/(Q*Z)).LT.pTdQdzcut_JLAB2022) ) THEN
                  CALL DISUU(s,Q2,xb,Z,pht,fuua,IT,IH,IC)
                  IT = 3
                  CALL DISUU(s,Q2,xb,Z,pht,fuu,IT,IH,IC)
                  DIS = pip2022_fedis(I)
                  R_a = fuua/(fuu*DIS)
                  sigma = sqrt(STAT**2 + SYS**2)
                  tmp = (R_a - MULT)**2d0/sigma**2d0
                  CHI2_PIP_2022_Fe = CHI2_PIP_2022_Fe + tmp
                  CHI2 = CHI2 + tmp
                  NUM_PIP_2022_Fe = NUM_PIP_2022_Fe +1
              
                  WRITE(208,*) Z,Q2,PT2,MULT,STAT,SYS,fuu,fuua,DIS
            ENDIF
      ENDDO
      CLOSE(208)
      ENDIF

      IF(I_pip_2022_Pb.EQ.1) THEN
      OPEN(UNIT = 209,
     *STATUS='UNKNOWN',FILE ='plot_data/JLAB2022/pi+Pb_pre.dat')
      WRITE(209,*) 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 208
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

      DO I=1,48
            MULT = pip2022_pb(I)
            STAT = pip2022_pbstat(I)
            SYS  = pip2022_pbsys(I)
            Q2   = pip2022_qpb(I)
            xb   = pip2022_xpb(I)
            PT2  = 0.5d0*(pip2022_ptlow(I)+pip2022_ptup(I))
            Z    = 0.5d0*(pip2022_zlow(I)+pip2022_zup(I))
            Q    = dsqrt(Q2)
            pht  = sqrt(PT2)
            s  = Sep_JLAB2022
            IT = 208
            IH = 1
            IC = 1
            IF(  (PT2.LT.pt2cut).AND.(z > zcut_low).AND.(z < zcut_up)
     >      .AND.( (pht/(Q*Z)).LT.pTdQdzcut_JLAB2022) ) THEN
                  CALL DISUU(s,Q2,xb,Z,pht,fuua,IT,IH,IC)
                  IT = 3
                  CALL DISUU(s,Q2,xb,Z,pht,fuu,IT,IH,IC)
                  DIS = pip2022_pbdis(I)
                  R_a = fuua/(fuu*DIS)
                  sigma = sqrt(STAT**2 + SYS**2)
                  tmp = (R_a - MULT)**2d0/sigma**2d0
                  CHI2_PIP_2022_Pb = CHI2_PIP_2022_Pb + tmp
                  CHI2 = CHI2 + tmp
                  NUM_PIP_2022_Pb = NUM_PIP_2022_Pb +1
              
                  WRITE(209,*) Z,Q2,PT2,MULT,STAT,SYS,fuu,fuua,DIS
            ENDIF
      ENDDO
      CLOSE(209)
      ENDIF

      IF(I_pim_2022_C.EQ.1) THEN
      OPEN(UNIT = 210,
     *STATUS='UNKNOWN',FILE ='plot_data/JLAB2022/pi-C_pre.dat')
      WRITE(210,*) 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 12
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

      DO I=1,40
            MULT = pim2022_c(I)
            STAT = pim2022_cstat(I)
            SYS  = pim2022_csys(I)
            Q2   = pim2022_qc(I)
            xb   = pim2022_xc(I)
            PT2  = 0.5d0*(pim2022_ptlow(I)+pim2022_ptup(I))
            Z    = 0.5d0*(pim2022_zlow(I)+pim2022_zup(I))
            Q    = dsqrt(Q2)
            pht  = sqrt(PT2)
            s  = Sep_JLAB2022
            IT = 12
            IH = 1
            IC = -1
            IF(  (PT2.LT.pt2cut).AND.(z > zcut_low).AND.(z < zcut_up)
     >      .AND.( (pht/(Q*Z)).LT.pTdQdzcut_JLAB2022) ) THEN
                  CALL DISUU(s,Q2,xb,Z,pht,fuua,IT,IH,IC)
                  IT = 3
                  CALL DISUU(s,Q2,xb,Z,pht,fuu,IT,IH,IC)
                  DIS = pim2022_cdis(I)
                  R_a = fuua/(fuu*DIS)
                  sigma = sqrt(STAT**2 + SYS**2)
                  tmp = (R_a - MULT)**2d0/sigma**2d0
                  CHI2_PIM_2022_C = CHI2_PIM_2022_C + tmp
                  CHI2 = CHI2 + tmp
                  NUM_PIM_2022_C = NUM_PIM_2022_C +1
              
                  WRITE(210,*) Z,Q2,PT2,MULT,STAT,SYS,fuu,fuua,DIS
            ENDIF
      ENDDO
      CLOSE(210)
      ENDIF

      IF(I_pim_2022_Fe.EQ.1) THEN
      OPEN(UNIT = 211,
     *STATUS='UNKNOWN',FILE ='plot_data/JLAB2022/pi-Fe_pre.dat')
      WRITE(211,*) 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 56
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

      DO I=1,40
            MULT = pim2022_fe(I)
            STAT = pim2022_festat(I)
            SYS  = pim2022_fesys(I)
            Q2   = pim2022_qfe(I)
            xb   = pim2022_xfe(I)
            PT2  = 0.5d0*(pim2022_ptlow(I)+pim2022_ptup(I))
            Z    = 0.5d0*(pim2022_zlow(I)+pim2022_zup(I))
            Q    = dsqrt(Q2)
            pht  = sqrt(PT2)
            s  = Sep_JLAB2022
            IT = 56
            IH = 1
            IC = -1
            IF(  (PT2.LT.pt2cut).AND.(z > zcut_low).AND.(z < zcut_up)
     >      .AND.( (pht/(Q*Z)).LT.pTdQdzcut_JLAB2022) ) THEN
                  CALL DISUU(s,Q2,xb,Z,pht,fuua,IT,IH,IC)
                  IT = 3
                  CALL DISUU(s,Q2,xb,Z,pht,fuu,IT,IH,IC)
                  DIS = pim2022_fedis(I)
                  R_a = fuua/(fuu*DIS)
                  sigma = sqrt(STAT**2 + SYS**2)
                  tmp = (R_a - MULT)**2d0/sigma**2d0
                  CHI2_PIM_2022_Fe = CHI2_PIM_2022_Fe + tmp
                  CHI2 = CHI2 + tmp
                  NUM_PIM_2022_Fe = NUM_PIM_2022_Fe +1
            
                  WRITE(211,*) Z,Q2,PT2,MULT,STAT,SYS,fuu,fuua,DIS
            ENDIF
      ENDDO
      CLOSE(211)
      ENDIF

      IF(I_pim_2022_Pb.EQ.1) THEN
      OPEN(UNIT = 212,
     *STATUS='UNKNOWN',FILE ='plot_data/JLAB2022/pi-Pb_pre.dat')
      WRITE(212,*) 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'

C----- Set-up -----------------------------------------------------------
      H_AA = 208
      lun = -6    ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 2 ! 2 for linear interpolarion, 3 for spline       
      call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
      call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
      nfin = 1 
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout

      call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

      DO I=1,40
            MULT = pim2022_pb(I)
            STAT = pim2022_pbstat(I)
            SYS  = pim2022_pbsys(I)
            Q2   = pim2022_qpb(I)
            xb   = pim2022_xpb(I)
            PT2  = 0.5d0*(pim2022_ptlow(I)+pim2022_ptup(I))
            Z    = 0.5d0*(pim2022_zlow(I)+pim2022_zup(I))
            Q    = dsqrt(Q2)
            pht  = sqrt(PT2)
            s  = Sep_JLAB2022
            IT = 208
            IH = 1
            IC = -1
            IF(  (PT2.LT.pt2cut).AND.(z > zcut_low).AND.(z < zcut_up)
     >      .AND.( (pht/(Q*Z)).LT.pTdQdzcut_JLAB2022) ) THEN
                  CALL DISUU(s,Q2,xb,Z,pht,fuua,IT,IH,IC)
                  IT = 3
                  CALL DISUU(s,Q2,xb,Z,pht,fuu,IT,IH,IC)
                  DIS = pim2022_pbdis(I)
                  R_a = fuua/(fuu*DIS)
                  sigma = sqrt(STAT**2 + SYS**2)
                  tmp = (R_a - MULT)**2d0/sigma**2d0
                  CHI2_PIM_2022_Pb = CHI2_PIM_2022_Pb + tmp
                  CHI2 = CHI2 + tmp
                  NUM_PIM_2022_Pb = NUM_PIM_2022_Pb +1
          
                  WRITE(212,*) Z,Q2,PT2,MULT,STAT,SYS,fuu,fuua,DIS
            ENDIF
      ENDDO
      CLOSE(212)
      ENDIF

      endif !(from normalization check)

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
     >       +NUM_PIP_2022_C !2023
     >       +NUM_PIP_2022_Fe
     >       +NUM_PIP_2022_Pb
     >       +NUM_PIM_2022_C
     >       +NUM_PIM_2022_Fe
     >       +NUM_PIM_2022_Pb

      chisquare_calculate = chi2

      return
      end

      SUBROUTINE BILINEAR(XVAL, YVAL, dy, M)
      INTEGER M, N, L
      REAL*8 X(100), Y(100), F(100,100), dy
      REAL*8 XVAL, YVAL 
      INTEGER I, J, k
      REAL*8 A, B
      REAL*8 g00, g01, g10, g11
      REAL*8 g0, g1
      character(len=70) :: filename, filename1, I_str, I_str1
      CHARACTER*72 line, dummy
      real*8, dimension(7) :: data_array
      real*8, dimension(4) :: data_array1
      CHARACTER(len=20) :: filename_patterns(63)

      N = 100

      DO I = 1, N
        X(I) = 0.0 + (I-1)*4.0/(N-1)
      ENDDO
      DO J = 1, N
        Y(J) = (J-1)*0.5/(N-1)
      ENDDO

      I = 1
      DO WHILE (XVAL > X(I+1) .AND. I < (N-1))
        I = I + 1
      END DO
      J = 1
      DO WHILE (YVAL > Y(J+1) .AND. J < (N-1))
        J = J + 1
      END DO

      !print *, 'I,J,NG2,Yvalue =', I, J, NG2, YVAL
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
      filename_patterns(35) = 'E866_Fe10.dat'
      filename_patterns(36) = 'E866_Fe11.dat'
      filename_patterns(37) = 'E866_Fe12.dat'
      filename_patterns(38) = 'E866_Fe18.dat'
      filename_patterns(39) = 'E866_Fe19.dat'
      filename_patterns(40) = 'E866_Fe20.dat'
      filename_patterns(41) = 'E866_Fe21.dat'
      filename_patterns(42) = 'E866_Fe26.dat'
      filename_patterns(43) = 'E866_Fe27.dat'
      filename_patterns(44) = 'E866_Fe28.dat'
      filename_patterns(45) = 'E866_Fe29.dat'

      filename_patterns(46) = 'E866_W1.dat'
      filename_patterns(47) = 'E866_W2.dat'
      filename_patterns(48) = 'E866_W3.dat'
      filename_patterns(49) = 'E866_W10.dat'
      filename_patterns(50) = 'E866_W11.dat'
      filename_patterns(51) = 'E866_W12.dat'
      filename_patterns(52) = 'E866_W18.dat'
      filename_patterns(53) = 'E866_W19.dat'
      filename_patterns(54) = 'E866_W20.dat'
      filename_patterns(55) = 'E866_W21.dat'
      filename_patterns(56) = 'E866_W26.dat'
      filename_patterns(57) = 'E866_W27.dat'
      filename_patterns(58) = 'E866_W28.dat'
      filename_patterns(59) = 'E866_W29.dat'

      filename_patterns(60) = 'RHIC_1.dat'
      filename_patterns(61) = 'RHIC_2.dat'
      filename_patterns(62) = 'RHIC_3.dat'
      filename_patterns(63) = 'RHIC_4.dat'
 
      filename = '../pwr_grid/rep_' // trim(I_str) //
     &              trim('/grid_data/') // trim(filename_patterns(M))
      filename1 = '../pwr_grid/rep_' // trim(I_str1) // 
     &              trim('/grid_data/') // trim(filename_patterns(M))

C---------------------------------Openfile----------------
      IF ( ((M .GT. 0) .AND. (M .LT. 17)) .OR. 
     &     (M .GT. 31) ) THEN

        OPEN (UNIT = 1, FILE = filename)
        READ(1,'(A)') LINE    ! skip first row
        DO K = 1,J
          read(1,*) data_array, dummy
        ENDDO
        F(I,J) = data_array(6)
        read(1,*) data_array, dummy
        F(I,J+1) = data_array(6)
        CLOSE(1)

        OPEN (UNIT = 2, FILE = filename1)
        READ(2,'(A)') LINE    ! skip first row
        DO K = 1,J
         read(2,*) data_array, dummy
        ENDDO
        F(I+1,J) = data_array(6)
        read(2,*) data_array, dummy
        F(I+1,J+1) = data_array(6)
        CLOSE(2)

      ELSEIF ( (M .GT. 16) .AND. (M .LT. 32)) THEN

        OPEN (UNIT = 1, FILE = filename)
        READ(1,'(A)') LINE    ! skip first row
        DO K = 1,J
          read(1,*) data_array1
          !print *, 'YVAL, FILE = ', yval, data_array1(2)
        ENDDO
        F(I,J) = data_array1(3)
        read(1,*) data_array1
        F(I,J+1) = data_array1(3)
        CLOSE(1)

        OPEN (UNIT = 2, FILE = filename1)
        READ(2,'(A)') LINE    ! skip first row
        DO K = 1,J
          read(2,*) data_array1
        ENDDO
        F(I+1,J) = data_array1(3)
        read(2,*) data_array1
        F(I+1,J+1) = data_array1(3)
        CLOSE(2)

      ENDIF
      CLOSE(1)
      CLOSE(2)

C      OPEN (UNIT = 1, FILE = filename)
C      READ(1,'(A)') LINE    ! skip first row
C        DO K = 1,J
C        read(1,*) data_array, dummy
C        ENDDO
C      F(I,J) = data_array(6)
C      read(1,*) data_array, dummy
C      F(I,J+1) = data_array(6)
C      CLOSE(1)
C
C      OPEN (UNIT = 2, FILE = filename1)
C      READ(2,'(A)') LINE    ! skip first row
C        DO K = 1,J
C        read(2,*) data_array, dummy
C        ENDDO
C      F(I+1,J) = data_array(6)
C      read(2,*) data_array, dummy
C      F(I+1,J+1) = data_array(6)
C      CLOSE(2)

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

      Subroutine Check(bool)
      implicit none
      integer i, bool, j
      double precision Ni(0:6),ai(0:6),bi(0:6),gi(0:6),di(0:6)
      double precision Nq1,Nq2,Ng1,Ng2 
      double precision aq1,bq1,gq1,dq1,ag1,bg1,gg1,dg1
      double precision aq2,bq2,gq2,dq2,ag2,bg2,gg2,dg2
      REAL*8 gamma, g3f, g3D, p_10, p_11, p_12, A !2023

      real*8, dimension(6) :: A_array = (/1,20,56,84,131,208/)

      COMMON /FITP/ gamma, g3f, g3D, Nq1, gq1, dq1,
     &              Nq2, gq2, dq2, p_10, p_11, p_12 !2023 

      bool = 1
         
      DO j=1,6
        A = A_array(j)
        DO i=1,6

*     utot
      Ni(1)  = 0.387d0 
      ai(1)  = -0.388d0
      bi(1)  = 0.910d0
      gi(1)  = 7.15d0
      di(1)  = 3.96d0

*     dtot
      Ni(2)  = 0.388d0 
      ai(2)  = ai(1)
      bi(2)  = bi(1)
      gi(2)  = gi(1)
      di(2)  = di(1)

*     ubar = d
      Ni(3) = 0.105d0 
      ai(3) = 1.649d0
      bi(3) = 3.286d0
      gi(3) = 49.95d0
      di(3) =  8.67d0

*     stot
      Ni(4)  = 0.273d0 
      ai(4)  = 1.449d0
      bi(4)  = bi(3)
      gi(4)  = gi(3)
      di(4)  = di(3)

*     ctot
      Ni(5)  = 0.306d0 
      ai(5)  = 1.345d0
      bi(5)  = 5.519d0
      gi(5)  = 19.78d0
      di(5)  = 10.22d0

*     btot
      Ni(6)  =  0.372d0 
      ai(6)  = -0.127d0 
      bi(6)  =  4.490d0 
      gi(6)  =  24.49d0 
      di(6)  =  12.80d0

          ai(i) = ai(i) 
          bi(i) = bi(i) + p_10 *(1d0-A**p_11)
          di(i) = di(i) + dq1 *(1d0-A**dq2)         

          IF ((ai(i) .LT. (-2)) .OR. 
     *        (bi(i) .LT. (-1)) .OR.
     *        ((bi(i)+di(i)) .LT. (-1))) THEN
            bool = 0
          ENDIF
        ENDDO
      ENDDO 

      RETURN
      END

