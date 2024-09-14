      PROGRAM MASTER
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL FCNG1
      DIMENSION NPRM(12),VSTRT(12),STP(12),BL(12),BUP(12),ARGLIS(12) !2023
      CHARACTER*10 PNAM(12) !2023
      INTEGER IRD,IWR,ISAV
      INTEGER NFIT
      INTEGER DAY, HOUR, MINUTE, SECOND
      COMMON /PARAMS/ NFIT
      INTEGER IFLAG
      COMMON  / IEND / IFLAG

      INTEGER i,j,k

C     INITALIZE THE PARAMETRIZATION
      REAL*8 gamma, g3f, g3D, Nq1, gq1, dq1,
     &       Nq2, gq2, dq2, p_10, p_11, p_12 !2023 
      COMMON /FITP/ gamma, g3f, g3D, Nq1, gq1, dq1,
     &              Nq2, gq2, dq2, p_10, p_11, p_12 !2023

      CHARACTER(len=20) :: E772_files(16),CMS5_files(8),ATLAS1_files(7),      
     &                     E866_files(28),RHIC_files(4)
      Integer E866_sets(14)

C------- specify values of parameters in line 217.

      data NPRM /1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/ !2023
      data PNAM /'gamma', 'g3f', 'g3D', 
     &           'Nq1', 'gq1', 'dq1',
     &           'Nq2', 'gq2', 'dq2',  
     &           'p_10', 'p_11', 'p_12'/ !2023

C-----STARTING STEP SIZE 
      data STP /0.1, 0.025, 0.005, 
     &          0.5, 0.5, 0.5,
     &          0.2, 0.2, 0.2, 
     &          100.0, 100.0, 100.0/ !2023

C-----LOWER BOUND (LIMIT) ON PARAMETER VALUE, IF ANY
      data BL   /0.0, 0.0, 0.0, 
     &           -5.0, -5.0, -5.0,
     &           -2.0, -2.0, -2.0, 
     &           -100.0, -100.0, -100.0/ !2023

C-----UPPER BOUND (LIMIT) ON PARAMETER VALUE, IF ANY
      data BUP  /4.0, 0.5, 0.1,  
     &           5.0, 5.0, 5.0,
     &           2.0, 2.0, 2.0,  
     &           100.0, 100.0, 100.0/ !2023

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

 1115 CONTINUE

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
      I_E866_800q  = 0 !1 !1

C-----E772
      I_E772_800   = 1

      LIKEn = 0 
      newFF = 1 
      newAPFEL = 2

c-----FIT PARAMETERS
      nloops = 2 ! LO: 1 NLO: 2
      nll = 3
      pre = 1
      collFF = newFF

!     CALL SETLHAPARM('SILENT') ! TO NOT SHOW THE CALLS, ALTHOUGH THEY ARE CALLED
      CALL SETLHAPARM('SILENT')

      ! TO NOT SHOW THE APFEL CALLS, ALTHOUGH THEY ARE CALLED
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
      CALL InitPDFsetByNameM(16,"LIKEnCC")
      CALL InitPDFM(16,0)
      print *, "17"
      CALL InitPDFsetByNameM(17,"LIKEnFE")
      CALL InitPDFM(17,0)
      print *, "18"
      CALL InitPDFsetByNameM(18,"LIKEnPB")
      CALL InitPDFM(18,0)
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

      print*, "start initializaition"

C---------------Open files E772-----------------------------------------
     
      E772_files(1) = 'E772_C1.dat'
      E772_files(2) = 'E772_C2.dat' 
      E772_files(3) = 'E772_C3.dat' 
      E772_files(4) = 'E772_C4.dat' 
                   
      E772_files(5) = 'E772_Ca1.dat' 
      E772_files(6) = 'E772_Ca2.dat' 
      E772_files(7) = 'E772_Ca3.dat' 
      E772_files(8) = 'E772_Ca4.dat'
                   
      E772_files(9) = 'E772_Fe1.dat' 
      E772_files(10) = 'E772_Fe2.dat' 
      E772_files(11) = 'E772_Fe3.dat' 
      E772_files(12) = 'E772_Fe4.dat' 
                   
      E772_files(13) = 'E772_W1.dat' 
      E772_files(14) = 'E772_W2.dat' 
      E772_files(15) = 'E772_W3.dat' 
      E772_files(16) = 'E772_W4.dat'

      if(I_E772_800.eq.1) then
        DO i = 1,16
          open(unit=i+100, file='grid_data/' // trim(E772_files(i)), 
     &         status='unknown')
          write(i+100, *) 'gamma ', 'g3f ', 'pt ','DY-RATIO ', 
     &                'error ', 'R_dy ', 'tmp ', 'element '
        END DO
      END IF

C---------------Open files CMS5-----------------------------------------

      CMS5_files(1) = 'CMS5_1.dat'
      CMS5_files(2) = 'CMS5_2.dat'
      CMS5_files(3) = 'CMS5_3.dat'
      CMS5_files(4) = 'CMS5_4.dat'
      CMS5_files(5) = 'CMS5_5.dat'
      CMS5_files(6) = 'CMS5_6.dat'
      CMS5_files(7) = 'CMS5_7.dat'
      CMS5_files(8) = 'CMS5_8.dat'

      if(I_CMS5.eq.1) then
        DO i = 1,8
          open(unit=i+200, file='grid_data/' // trim(CMS5_files(i)),
     &         status='unknown')
          write(i+200, *) 'gamma ', 'g3f ', 'fuua ', 'tmp ' 
        END DO
      END IF

C---------------Open files ATLAS1---------------------------------------

      ATLAS1_files(1) = 'ATLAS1_1.dat'
      ATLAS1_files(2) = 'ATLAS1_2.dat'
      ATLAS1_files(3) = 'ATLAS1_3.dat'
      ATLAS1_files(4) = 'ATLAS1_4.dat'
      ATLAS1_files(5) = 'ATLAS1_5.dat'
      ATLAS1_files(6) = 'ATLAS1_6.dat'
      ATLAS1_files(7) = 'ATLAS1_7.dat'

      if(I_ATLAS5_Y1.eq.1) then
        DO i = 1,7
          open(unit=i+300, file='grid_data/' // trim(ATLAS1_files(i)),
     &         status='unknown')
          write(i+300, *) 'gamma ', 'g3f ', 'fuua ', 'tmp '
        END DO
      END IF

C---------------Open files E866-----------------------------------------

      E866_files(1) = 'E866_Fe1.dat'
      E866_files(2) = 'E866_Fe2.dat'
      E866_files(3) = 'E866_Fe3.dat'
      E866_files(4) = 'E866_Fe10.dat'
      E866_files(5) = 'E866_Fe11.dat'
      E866_files(6) = 'E866_Fe12.dat'
      E866_files(7) = 'E866_Fe18.dat'
      E866_files(8) = 'E866_Fe19.dat'
      E866_files(9) = 'E866_Fe20.dat'
      E866_files(10) = 'E866_Fe21.dat'
      E866_files(11) = 'E866_Fe26.dat'
      E866_files(12) = 'E866_Fe27.dat'
      E866_files(13) = 'E866_Fe28.dat'
      E866_files(14) = 'E866_Fe29.dat'

      E866_files(15) = 'E866_W1.dat'
      E866_files(16) = 'E866_W2.dat'
      E866_files(17) = 'E866_W3.dat'
      E866_files(18) = 'E866_W10.dat'
      E866_files(19) = 'E866_W11.dat'
      E866_files(20) = 'E866_W12.dat'
      E866_files(21) = 'E866_W18.dat'
      E866_files(22) = 'E866_W19.dat'
      E866_files(23) = 'E866_W20.dat'
      E866_files(24) = 'E866_W21.dat'
      E866_files(25) = 'E866_W26.dat'
      E866_files(26) = 'E866_W27.dat'
      E866_files(27) = 'E866_W28.dat'
      E866_files(28) = 'E866_W29.dat'

      DATA E866_sets /1,2,3,10,11,12,18,19,20,21,26,27,28,29/

      if(I_E866_800q.eq.1) then
        DO I = 1,14
          open(unit=(E866_sets(i)+400), 
     &         file='grid_data/' // trim(E866_files(i)),
     &         status='unknown')
          write(E866_sets(i)+400, *) 'gamma ', 'g3f ', 'pt '
     &        ,'DY-RATIO ', 'error ', 'R_dy ', 'tmp ', 'element '
        ENDDO

        DO I = 1,14
          open(unit=(E866_sets(i)+500), 
     &         file='grid_data/' // trim(E866_files(i+14)),
     &         status='unknown')
          write(E866_sets(i)+500, *) 'gamma ', 'g3f ', 'pt '
     &        ,'DY-RATIO ', 'error ', 'R_dy ', 'tmp ', 'element '
        ENDDO
      END IF

C---------------Open files RHIC-----------------------------------------

      RHIC_files(1) = 'RHIC_1.dat'
      RHIC_files(2) = 'RHIC_2.dat'
      RHIC_files(3) = 'RHIC_3.dat'
      RHIC_files(4) = 'RHIC_4.dat'

      if(I_RHIC_Ratio_pAu1.eq.1) then
        DO i = 1,4
          open(unit=i+600, file='grid_data/' // trim(RHIC_files(i)),
     &         status='unknown')
          write(i+600, *) 'gamma ', 'g3f ', 'pt ', 'CX ',
     &                'error ', 'R_dy ', 'tmp ', 'element '
        END DO
      END IF

C---------------Parameter assignment------------------------------------

      j = 100
      j = j - 100
      DO k = 0,99

      gammaBEST = j*(BUP(1)-BL(1))/99 + BL(1)
      g3fBEST = k*(BUP(2)-BL(2))/99 + BL(2)

C-----------------------------------------------------------------------

      VSTRT(1)  =  gammaBEST
      VSTRT(2)  =  g3fBEST
      VSTRT(3)  =  0.0
      VSTRT(4)  =  0.0
      VSTRT(5)  =  0.0
      VSTRT(6)  =  0.0
      VSTRT(7)  =  0.0
      VSTRT(8)  =  0.0
      VSTRT(9)  =  0.0
      VSTRT(10)  =  0.0
      VSTRT(11)  =  0.0
      VSTRT(12)  =  0.0

C-----RHIC Ratios
      I_RHIC_Ratio_pAu1 = 1 !1
      I_RHIC_Ratio_pAu2 = 1 !1

C-----ATLAS 5 TEV
      I_ATLAS5_Y1 = 1!1
      I_ATLAS5_Y2 = 0
      I_ATLAS5_Y3 = 0

C-----CMS 5 TEV
      I_CMS5 = 1 !1

C-----E866
      I_E866_800q  = 0 !1

C-----E772
      I_E772_800   = 1

C     NUMBER OF PARAMETERS
      NFIT = 12 !2023

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
*.....OUTPUT PRINT LEVEL (FROM -1 TO 3)
      ARGLIS(1) = 3.
      CALL MNEXCM (FCNG1,'SET PRINT',ARGLIS, 1, IERFLG, DUM)

*     THIS CREATES THE FILE WITH THE RESULTS
      ARGLIS(1) = 3.            !   IFLAG = 3
      CALL MNEXCM (FCNG1, 'CALL FCN', ARGLIS, 1, IERFLG, DUM)
*.....STOP
      CALL MNEXCM (FCNG1, 'STOP', ARGLIS, 1, IERFLG, DUM)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C
*.....CALCULATE AND PRINT THE TOTAL TIME OF THE PREDICTION
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

      ENDDO !2023

C-----------------------------CLOSE----------------------------
      
      DO I = 1,30
      close(I+100)
      close(I+200)
      close(I+300)
      close(I+400)
      close(I+500)
      close(I+600)
      ENDDO

C-----------------------------CLOSE----------------------------

      RETURN
      END

      SUBROUTINE FCNG1 (NPAR, G, F, X, IIFLAG, dum)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION X(12), G(12)
      integer iflag,iiflag
      COMMON  / iend / iflag

      iflag = iiflag
      F = chisquare(X)

      RETURN
      END

      function chisquare(XX) 
      integer nfit
      common /params/ nfit
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

C---------------

      REAL*8 CMS5_NORM
      REAL*8 ATLAS5_Y1_NORM,ATLAS5_Y2_NORM,ATLAS5_Y3_NORM
      REAL*8 CX_STORE(100)

      integer IT,IH,IC
      common /meson/ IH,IC

      REAL*8 qTdQcut,phtcut,zcut
      REAL*8 pt2cut,zcut_low,zcut_up, pTdQdZcut,qcut

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

      integer iflag,num_exp
      common /iend/ iflag

      !LHC variables
      real*8  RTS_LHC

      integer cutflag
      real*8 ptcut,etamin,etamax,cutfac
      common /cuts/ ptcut,etamin,etamax,cutflag,qflg

      !HERMES variables
      real*8 MULT, NU, Q2, PT2, STAT, SYS, VAL, Z
	
      !JLAB2022 variables
      real*8 zlow, zup, bin, ptlow, ptup, c, cstat
      real*8 csys, fe, festat, fesys, pb,pbstat,pbsys
      real*8 cerr, feerr, pberr

      real*8 CHI2
      real*8 xb,zh,pht,tmp
      real*8 sigma
      real*8 dcorr,fN
      real*8 rts,y,Q,Qmin,Qmax,pt,fuu,Rds
      real*8 fuuA,fuuB, DIS, R_a
      real*8 R_dy
      real*8 xfmin, xfmax, xbbar, ymin, ymax, Qbar, const, ybar, yy, xf
      real*8 xbmin, xbmax
      real*8 pi
      data pi/3.1415926535d0/
      real*8 M
      data M/0.938d0/
      real*8 CX_LHC,CX_RHIC
      integer i
      real*8 test_pdf, test_ff

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
      pt2cut   = 0.3d0
      zcut_up  = 0.7d0
      zcut_low = 0.0
      pTdQdzcut = 0.5d0 !2023
      qcut = 10d0
      zcut = 0.7d0

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
      NUM_PIP_2022_C = 0 !2023
      NUM_PIP_2022_Fe = 0
      NUM_PIP_2022_Pb = 0
      NUM_PIM_2022_C = 0
      NUM_PIM_2022_Fe = 0
      NUM_PIM_2022_Pb = 0
C-------------

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

C-------------

      ktw_value = 0.424d0
      ptw_value = 0.168d0

C-------------DY (RHIC) p + Au -> muons (Au-going)
      if(I_RHIC_Ratio_pAu1.eq.1) then
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
           !print *, 'RHIC ', i
           IT = 197
           call DY_overyu(1,rts,pt,ymin,ymax,Qmin,Qmax,fuuA,IT)
           R_dy = fuuA/fuuB
           tmp = (R_dy - CX_RHIC)**2d0/sigma**2d0
           write(I+600,*)  gamma, g3f,pt, CX_RHIC, 
     &                     sigma, R_dy, tmp, 'Au'
           CHI2_RHIC_Ratio_pAu1 = CHI2_RHIC_Ratio_pAu1 + tmp
           CHI2 = CHI2 + tmp
           num_RHIC_Ratio_pAu1 = num_RHIC_Ratio_pAu1 + 1
        endif
      enddo
      endif

C-------------DY (RHIC) p + Au -> muons (p-going)
      if(I_RHIC_Ratio_pAu2.eq.1) then
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
           !print *, 'RHIC',i
           ! Compute X section for p+Au (A = 195 OR 197 available)
           IT = 197
           call DY_overyu(1,rts,pt,ymin,ymax,Qmin,Qmax,fuuA,IT)
           R_dy = fuuA/fuuB
           tmp = (R_dy - CX_RHIC)**2d0/sigma**2d0
           write(I+602,*)  gamma, g3f,pt, CX_RHIC,
     &                     sigma, R_dy, tmp, 'Au'
           CHI2_RHIC_Ratio_pAu2 = CHI2_RHIC_Ratio_pAu2 + tmp
           CHI2 = CHI2 + tmp
           num_RHIC_Ratio_pAu2 = num_RHIC_Ratio_pAu2 + 1
        endif
      enddo
      endif

C--------------------------------------------------------------

      normsumnum = 0d0
      normsumden = 0d0
      cutflag = 1
      etamin = -2.4d0
      etamax =  2.4d0
      ptcut = 20d0
      if(I_CMS5.eq.1) then
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
           !print *, 'CMS',i
         call DY_overyu(1,rts,pt,ymin,ymax,60d0,120d0,fuuA,IT)
         FUUA = 2d0*pi*pt*FUUA*208
         tmp = 0
         write(I+200,*) gamma, g3f, fuua, tmp
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
         endif
      enddo
      CHI2_CMS5      = CHI2_CMS5+(1.-fN)**2./(dcorr**2.)
      chi2           = chi2     +(1.-fN)**2./(dcorr**2.)
      endif
      cutflag = 0

C--------------------------------------------------------------

      normsumnum = 0d0
      normsumden = 0d0
      if(I_ATLAS5_Y1.eq.1) then
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
           !print *, 'ATLAS',i
         call DY_overyu(1,rts,pt,ymin,ymax,66d0,116d0,fuuA,IT)
         FUUA = 2d0*pi*pt*FUUA/PT*208
         tmp = 0
         write(I+300,*) gamma, g3f, fuua, tmp
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
         endif
      enddo
      chi2_atlas5_Y1 = CHI2_ATLAS5_Y1+(1.-fN)**2./(dcorr**2.)
      chi2           = chi2          +(1.-fN)**2./(dcorr**2.)
      endif

C------------------------------------------------------------------------

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
             fuua =     norm*CX_STORE(I)
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
             fuua =     norm*CX_STORE(I)
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

C------E772 DY_overQ(Ixxfy,rtss,xxf,yy,qtt,QQmin,QQmax,fuu,ITT)
C---------DY_overxbQ(Ixxfy,rtss,qqt,xbmin,xbmax,QQmin,QQmax,fuu,ITT)
C------- C/D

      if(I_E772_800.eq.1) then
      do I=1,7
      rts   = dsqrt(1600d0*0.938d0)
      pt    = E772_800_C_pT(I)
      Rds   = E772_800_C_R(I)
      sigma = E772_800_C_Err(I)
      Qmin = 4d0
      Qmax = 9d0
      Qbar = (Qmin+Qmax)/2d0
      xbmin = 0.05d0
      xbmax = 0.3d0
      xbbar = (xbmin+xbmax)/2d0
      xf = 0.26
      yy = ASINH((rts**2)*xf/(2*Qbar))
      if (pt/Qbar.lt.qTdQcut) then
           !print *, 'E772',i
      !call DY_overQ(4,rts,xbbar,yy,pt,Qmin,Qmax,fuuB,3)
      !call DY_overQ(4,rts,xbbar,yy,pt,Qmin,Qmax,fuuA,12)
      call DY_overxbQ(4,rts,pt,xbmin,xbmax,Qmin,Qmax,fuuB,3)
      call DY_overxbQ(4,rts,pt,xbmin,xbmax,Qmin,Qmax,fuuA,12)
      R_dy = fuuA/fuuB
      tmp = (R_dy - Rds)**2d0/sigma**2d0
      CHI2_E772_800 = CHI2_E772_800 + tmp
      CHI2 = CHI2 + tmp
      write(I+100,*) gamma, g3f, pt, Rds, sigma, R_dy, tmp, 'C'    

      num_E772_800 = num_E772_800 + 1
      endif
      enddo
      endif

C--------- Ca/D
      if(I_E772_800.eq.1) then
      do I=1,7
      rts   = dsqrt(1600d0*0.938d0)
      pt    = E772_800_Ca_pT(I)
      Rds   = E772_800_Ca_R(I)
      sigma = E772_800_Ca_Err(I)
      Qmin = 4d0
      Qmax = 9d0
      Qbar = (Qmin+Qmax)/2d0
      xbmin = 0.05d0
      xbmax = 0.3d0
      xbbar = (xbmin+xbmax)/2d0
      xf = 0.26
      yy = ASINH((rts**2)*xf/(2*Qbar))
      print *, xbbar
      if (pt/Qbar.lt.qTdQcut) then
           !print *, 'E772',i
      !call DY_overQ(4,rts,xbbar,yy,pt,Qmin,Qmax,fuuB,3)
      !call DY_overQ(4,rts,xbbar,yy,pt,Qmin,Qmax,fuuA,40)
      call DY_overxbQ(4,rts,pt,xbmin,xbmax,Qmin,Qmax,fuuB,3)
      call DY_overxbQ(4,rts,pt,xbmin,xbmax,Qmin,Qmax,fuuA,40)
      R_dy = fuuA/fuuB
      tmp = (R_dy - Rds)**2d0/sigma**2d0
      CHI2_E772_800 = CHI2_E772_800 + tmp
      CHI2 = CHI2 + tmp
      write(I+104,*) gamma, g3f, pt, Rds, sigma, R_dy, tmp, 'Ca'     

      num_E772_800 = num_E772_800 + 1
      endif
      enddo
      endif

C--------- Fe/D
      if(I_E772_800.eq.1) then
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
      xbbar = (xbmin+xbmax)/2d0
      xf = 0.26
      yy = ASINH((rts**2)*xf/(2*Qbar))
      if (pt/Qbar.lt.qTdQcut) then
           !print *, 'E772',i
      !call DY_overQ(4,rts,xbbar,yy,pt,Qmin,Qmax,fuuB,3)
      !call DY_overQ(4,rts,xbbar,yy,pt,Qmin,Qmax,fuuA,56)
      call DY_overxbQ(4,rts,pt,xbmin,xbmax,Qmin,Qmax,fuuB,3)
      call DY_overxbQ(4,rts,pt,xbmin,xbmax,Qmin,Qmax,fuuA,56)
      R_dy = fuuA/fuuB
      tmp = (R_dy - Rds)**2d0/sigma**2d0
      CHI2_E772_800 = CHI2_E772_800 + tmp
      CHI2 = CHI2 + tmp
      write(I+108,*)  gamma, g3f, pt, Rds, sigma, R_dy, tmp, 'Fe'

      num_E772_800 = num_E772_800 + 1
      endif
      enddo
      endif

C--------- W/D
      if(I_E772_800.eq.1) then
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
      xbbar = (xbmin+xbmax)/2d0
      xf = 0.26
      yy = ASINH((rts**2)*xf/(2*Qbar))
      if (pt/Qbar.lt.qTdQcut) then
           !print *, 'E772',i
      !call DY_overQ(4,rts,xbbar,yy,pt,Qmin,Qmax,fuuB,3)
      !call DY_overQ(4,rts,xbbar,yy,pt,Qmin,Qmax,fuuA,184)
      call DY_overxbQ(4,rts,pt,xbmin,xbmax,Qmin,Qmax,fuuB,3)
      call DY_overxbQ(4,rts,pt,xbmin,xbmax,Qmin,Qmax,fuuA,184)
      R_dy = fuuA/fuuB
   
      tmp = (R_dy - Rds)**2d0/sigma**2d0
      CHI2_E772_800 = CHI2_E772_800 + tmp
   
      CHI2 = CHI2 + tmp
      num_E772_800 = num_E772_800 + 1
      write(I+112,*)  gamma, g3f,pt, Rds, sigma, R_dy, tmp, 'W'
 
      endif
      enddo
      endif

C-----E866 (Q binned data)
      if(I_E866_800q.eq.1) then
      do I=1,32
      rts   = dsqrt(1600d0*0.938d0)
      pt    = E866_800q_pT(I)
      Rds   = E866_800q_R_FeBe(I)
      Q     = E866_800q_Q(I)
      sigma = E866_800q_Err_FeBe(I)
      xfmin = 0.13d0
      xfmax = 0.93d0
      if (pt/Q.lt.qTdQcut) then
           !print *, 'E866',i
      call DY_overxF(2,rts,pt,xfmin,xfmax,Q,fuuB,9)
      call DY_overxF(2,rts,pt,xfmin,xfmax,Q,fuuA,56)
      R_dy = fuuA/fuuB
      tmp = (R_dy - Rds)**2d0/sigma**2d0
      write(I+400,*)  gamma, g3f,pt, Rds, sigma, R_dy, tmp, 'Fe'
      CHI2_E866_800q = CHI2_E866_800q + tmp
      CHI2 = CHI2 + tmp
      num_E866_800q = num_E866_800q + 1
      endif
      enddo
      endif

C------W/Be
      if(I_E866_800q.eq.1) then
      do I=1,32
      rts   = dsqrt(1600d0*0.938d0)
      pt    = E866_800q_pT(I)
      Rds   = E866_800q_R_WBe(I)
      Q     = E866_800q_Q(I)
      sigma = E866_800q_Err_WBe(I)
      xfmin = 0.13d0
      xfmax = 0.93d0
      if (pt/Q.lt.qTdQcut) then
           !print *, 'E866',i
      call DY_overxF(2,rts,pt,xfmin,xfmax,Q,fuuB,9)
      call DY_overxF(2,rts,pt,xfmin,xfmax,Q,fuuA,184)
      R_dy = fuuA/fuuB
      tmp = (R_dy - Rds)**2d0/sigma**2d0
      write(I+500,*)  gamma, g3f,pt, Rds, sigma, R_dy, tmp, 'W'
      CHI2_E866_800q = CHI2_E866_800q + tmp
      CHI2 = CHI2 + tmp
      num_E866_800q = num_E866_800q + 1
      endif
      enddo
      endif

      chisquare=CHI2

      num_exp = 0
     >       +num_E772_800
     >       +num_E866_800q
     >       +NUM_RHIC_Ratio_pAu1
     >       +NUM_RHIC_Ratio_pAu2
     >       +NUM_ATLAS5_Y1
     >       +NUM_ATLAS5_Y2
     >       +NUM_ATLAS5_Y3
     >       +NUM_CMS5

 100  format(7(1pE12.4))
 101  format(6(1pE10.2))
 110  format(A35,F10.3,I2,I2)
 112  format(A35,I5)
 113  format(A35,13(1pe12.4,','))
 114  format(A35,F10.3,',',I5,',',I5,',',F10.3)
 115  format(A25,F10.3,I2,',',A25,F10.3,I2)

 116  format(3(A25,F10.3,I2,I2,','))
 117  format(4(A25,F10.3,I2,I2,','))
 1117 format(4(A25,F10.3,I2,','))

         write(6,110) 'chi2:',chisquare
         write(6,112) 'total num of data in fit:',num_exp
         write(6,112) 'number of parameters:',nfit
         write(6,110) 'chi2/dof:',chisquare/(num_exp-nfit)
         WRITE(6,  *) '-------------------------'
         write(6,113) 'fitted parameters:'
      write(6,*) 'gamma, g3f' !2023
      write(6,*) xx(1), xx(2) !2023
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

      return
      end
