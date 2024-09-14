      SUBROUTINE READDATA
      IMPLICIT NONE
      REAL*8 X1(3000),X2(3000),X3(3000),
     $       X4(3000),X5(3000),X6(3000),
     $       X7(3000),X8(3000),X9(3000),X10(3000),
     $       X11(3000), X12(3000), X13(3000), X14(3000), X15(3000),
     $       X16(3000), X17(3000), X18(3000), X19(3000), X20(3000)

      CHARACTER*72 FILEEXP
      INTEGER NLINES,NTOT,NTOT2,I
      include "./data-inc.f"
      CHARACTER*8 rt

      rt = 'expdata/'

C-----CMS AT 8TEV
      NLINES = 1
      NTOT = 2
      FILEEXP = rt//'DRELLYAN/CMS8/DY_no.dat'
      CALL READDATA2(FILEEXP,NLINES,X1,X2,NTOT)
      DO I =1, NTOT
        CMS8_DY_no_PT(I) = X1(I)
        CMS8_DY_no_CX(I) = X2(I)
      ENDDO

      NLINES = 1
      NTOT = 2
      FILEEXP = rt//'DRELLYAN/CMS8/DY_fid.dat'
      CALL READDATA2(FILEEXP,NLINES,X1,X2,NTOT)
      DO I =1, NTOT
        CMS8_DY_fid_PT(I) = X1(I)
        CMS8_DY_fid_CX(I) = X2(I)
      ENDDO
      NLINES = 1
      NTOT = 9
      FILEEXP = rt//'DRELLYAN/CMS8/Z_no.dat'
      CALL READDATA2(FILEEXP,NLINES,X1,X2,NTOT)
      DO I =1, NTOT
        CMS8_ZZ_no_PT(I) = X1(I)
        CMS8_ZZ_no_CX(I) = X2(I)
      ENDDO

      NLINES = 1
      NTOT = 9
      FILEEXP = rt//'DRELLYAN/CMS8/Z_fid.dat'
      CALL READDATA2(FILEEXP,NLINES,X1,X2,NTOT)
      DO I =1, NTOT
        CMS8_ZZ_fid_PT(I) = X1(I)
        CMS8_ZZ_fid_CX(I) = X2(I)
      ENDDO

C-----ATLAS AT 3TEV
      NLINES = 1
      NTOT = 12
      FILEEXP = rt//'DRELLYAN/ATLAS3/1.dat'
      CALL READDATA8(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     $                              X6,X7,X8,NTOT)
      DO I =1, NTOT
        ATLAS3_PT(I) = X1(I)
        ATLAS3_CX(I) = X4(I)
        ATLAS3_ER(I) = DSQRT(X5(I)**2D0+X7(I)**2D0)
      ENDDO

C-----E906 at 120 GeV
      NLINES = 1
      NTOT = 5
      FILEEXP = rt//'DRELLYAN/E906/E906_C.csv'
      CALL READDATA3(FILEEXP, NLINES, X1, X2, X3, NTOT)
     $
      DO I = 1, NTOT
        E906_C_pT(I)  = X1(I)
        E906_C_Ra(I)  = X2(I)
        E906_C_err(I) = X3(I)
      ENDDO

      NLINES = 1
      NTOT = 5
      FILEEXP = rt//'DRELLYAN/E906/E906_FE.csv'
      CALL READDATA3(FILEEXP, NLINES, X1, X2, X3, NTOT)
     $
      DO I = 1, NTOT
        E906_Fe_pT(I)  = X1(I)
        E906_Fe_Ra(I)  = X2(I)
        E906_Fe_err(I) = X3(I)
      ENDDO


      NLINES = 1
      NTOT = 5
      FILEEXP = rt//'DRELLYAN/E906/E906_W.csv'
      CALL READDATA3(FILEEXP, NLINES, X1, X2, X3, NTOT)
     $
      DO I = 1, NTOT
        E906_W_pT(I)  = X1(I)
        E906_W_Ra(I)  = X2(I)
        E906_W_err(I) = X3(I)
      ENDDO



C-----CMS AT 5TEV
      NLINES = 1
      NTOT = 13
      FILEEXP = rt//'DRELLYAN/CMS5/1.dat'
      CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     $                            X6,X7,X8,X9,X10,NTOT)
      DO I =1, NTOT
        CMS5_PT(I) = X1(I)
        CMS5_CX(I) = X4(I)
        CMS5_ER(I) = DSQRT(X5(I)**2D0+X7(I)**2D0)
      ENDDO

C-----ATLAS AT 5TEV
      NLINES = 1
      NTOT = 14
      FILEEXP = rt//'DRELLYAN/ATLAS5/1.dat'
      CALL READDATA8(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     $                              X6,X7,X8,NTOT)
      DO I =1, NTOT
        ATLAS5_Y1_PT(I) = X1(I)
        ATLAS5_Y1_CX(I) = X4(I)
        ATLAS5_Y1_ER(I) = DSQRT(X5(I)**2D0+X7(I)**2D0)
      ENDDO

      NLINES = 1
      NTOT = 8
      FILEEXP = rt//'DRELLYAN/ATLAS5/2.dat'
      CALL READDATA8(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     $                              X6,X7,X8,NTOT)
      DO I =1, NTOT
        ATLAS5_Y2_PT(I) = X1(I)
        ATLAS5_Y2_CX(I) = X4(I)
        ATLAS5_Y2_ER(I) = DSQRT(X5(I)**2D0+X7(I)**2D0)
      ENDDO

      NLINES = 1
      NTOT = 8
      FILEEXP = rt//'DRELLYAN/ATLAS5/3.dat'
      CALL READDATA8(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     $                              X6,X7,X8,NTOT)
      DO I =1, NTOT
        ATLAS5_Y3_PT(I) = X1(I)
        ATLAS5_Y3_CX(I) = X4(I)
        ATLAS5_Y3_ER(I) = DSQRT(X5(I)**2D0+X7(I)**2D0)
      ENDDO

      FILEEXP = rt//'DRELLYAN/RHIC/RHIC_pAu1.csv'
      CALL READDATA5(FILEEXP, NLINES, X1,X2,X3,X4,X5,NTOT)
      DO I =1, NTOT
        RHIC_pAu1_pT(I)     = X1(I)
        RHIC_pAu1_CX(I)     = X2(I)
        RHIC_pAu1_err(I)    = X3(I)
      ENDDO

      FILEEXP = rt//'DRELLYAN/RHIC/RHIC_pAu2.csv'
      CALL READDATA5(FILEEXP, NLINES, X1,X2,X3,X4,X5,NTOT)
      DO I =1, NTOT
        RHIC_pAu2_pT(I)     = X1(I)
        RHIC_pAu2_CX(I)     = X2(I)
        RHIC_pAu2_err(I)    = X3(I)
      ENDDO

      FILEEXP = rt//'DRELLYAN/RHIC2/Ratio_Augoing.csv'
      CALL READDATA3(FILEEXP, 1, X1,X2,X3,NTOT)
      DO I =1, 5
        RHIC_Ratio_pAu1_pT(I)     = X1(I)
        RHIC_Ratio_pAu1_CX(I)     = X2(I)
        RHIC_Ratio_pAu1_err(I)    = X3(I)
      ENDDO

      FILEEXP = rt//'DRELLYAN/RHIC2/Ratio_pgoing.csv'
      CALL READDATA3(FILEEXP, 1, X1,X2,X3,NTOT)
      DO I =1, 5
        RHIC_Ratio_pAu2_pT(I)     = X1(I)
        RHIC_Ratio_pAu2_CX(I)     = X2(I)
        RHIC_Ratio_pAu2_err(I)    = X3(I)
      ENDDO

C-----RHIC
C-----RHIC RTS = 200 GEV
      NLINES = 1
      NTOT = 5
      FILEEXP = rt//'DRELLYAN/RHIC2/RHIC_pp.dat'
      CALL READDATA4(FILEEXP, NLINES, X1,X2,X3,X4, NTOT)
      DO I =1, NTOT
        RHIC_pp_pT(I)     = X1(I)
        RHIC_pp_CX(I)     = X2(I)
        RHIC_pp_err(I)    = X3(I)
        RHIC_pp_FUU(I)    = X4(I)
      ENDDO

      FILEEXP = rt//'DRELLYAN/RHIC/RHIC_pAu1.csv'
      CALL READDATA5(FILEEXP, NLINES, X1,X2,X3,X4,X5,NTOT)
      DO I =1, NTOT
        RHIC_pAu1_pT(I)     = X1(I)
        RHIC_pAu1_CX(I)     = X2(I)
        RHIC_pAu1_err(I)    = X3(I)
      ENDDO

      FILEEXP = rt//'DRELLYAN/RHIC/RHIC_pAu2.csv'
      CALL READDATA5(FILEEXP, NLINES, X1,X2,X3,X4,X5,NTOT)
      DO I =1, NTOT
        RHIC_pAu2_pT(I)     = X1(I)
        RHIC_pAu2_CX(I)     = X2(I)
        RHIC_pAu2_err(I)    = X3(I)
      ENDDO

C-----DRELL YAN
c-----E288: DY 200 GeV
      NLINES = 1
      NTOT = 119
      FILEEXP=rt//'DRELLYAN/E288/E288_200.csv'
      CALL READDATA6(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     $                X6,NTOT)
      DO I=1,NTOT
         E288_200_pT(I)   =X1(I)
         E288_200_ds(I)   =X4(I)
         E288_200_sigma(I)=X5(I)
         E288_200_Q(I)    =X6(I)

      ENDDO

c-----E288: DY 300 GeV
      NLINES = 1
      NTOT = 184
      FILEEXP=rt//'DRELLYAN/E288/E288_300.csv'
      CALL READDATA6(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     $                X6,NTOT)
      DO I=1,NTOT
         E288_300_pT(I)   =X1(I)
         E288_300_ds(I)   =X4(I)
         E288_300_sigma(I)=X5(I)
         E288_300_Q(I)    =X6(I)

      ENDDO

c-----E288: DY 400 GeV
      NLINES = 1
      NTOT = 225
      FILEEXP=rt//'DRELLYAN/E288/E288_400.csv'
      CALL READDATA6(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     $                X6,NTOT)
      DO I=1,NTOT
         E288_400_pT(I)   =X1(I)
         E288_400_ds(I)   =X4(I)
         E288_400_sigma(I)=X5(I)
         E288_400_Q(I)    =X6(I)

      ENDDO

C-----E866: DY 800 GeV (x2 bins)
      NLINES = 1
      NTOT = 40
      FILEEXP=rt//'DRELLYAN/E866/FNAL866_x2bins.csv'
      CALL READDATA8(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     $                X6,X7,X8, NTOT)
      DO I=1,NTOT
         E866_800_pT(I)       =X1(I)
         E866_800_R_FeBe(I)   =X2(I)
         E866_800_Err_FeBe(I) =X3(I)
         E866_800_R_WBe(I)    =X4(I)
         E866_800_Err_WBe(I)  =X5(I)
         E866_800_x2(I)       =X6(I)
      ENDDO

C-----E866: DY 800 GeV (Q bins)
      NLINES = 1
      NTOT = 32
      FILEEXP=rt//'DRELLYAN/E866/FNAL866_Qbins.csv'
      CALL READDATA8(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     $                X6,X7,X8, NTOT)
      DO I=1,NTOT
         E866_800q_pT(I)       =X1(I)
         E866_800q_R_FeBe(I)   =X2(I)
         E866_800q_Err_FeBe(I) =X3(I)
         E866_800q_R_WBe(I)    =X4(I)
         E866_800q_Err_WBe(I)  =X5(I)
         E866_800q_Q(I)        =X6(I)
      ENDDO

C-----E772: C/D 800 GeV
      NLINES = 1
      NTOT = 7
      FILEEXP=rt//'DRELLYAN/E772/E772_800_CD.csv'
      CALL READDATA3(FILEEXP,NLINES,X1,X2,X3,NTOT)
      DO I=1, NTOT
         E772_800_C_pT(I)    =X1(I)
         E772_800_C_R(I)     =X2(I)
         E772_800_C_Err(I)   =X3(I)
      ENDDO

C-----E772: Ca/D 800 GeV
      NLINES = 1
      NTOT = 7
      FILEEXP=rt//'DRELLYAN/E772/E772_800_CaD.csv'
      CALL READDATA3(FILEEXP,NLINES,X1,X2,X3,NTOT)
      DO I=1, NTOT
         E772_800_Ca_pT(I)    =X1(I)
         E772_800_Ca_R(I)     =X2(I)
         E772_800_Ca_Err(I)   =X3(I)
      ENDDO

C-----E772: Fe/D 800 GeV
      NLINES = 1
      NTOT = 7
      FILEEXP=rt//'DRELLYAN/E772/E772_800_FeD.csv'
      CALL READDATA3(FILEEXP,NLINES,X1,X2,X3,NTOT)
      DO I=1, NTOT
         E772_800_Fe_pT(I)    =X1(I)
         E772_800_Fe_R(I)     =X2(I)
         E772_800_Fe_Err(I)   =X3(I)
      ENDDO

C-----E772: W/D 800 GeV
      NLINES = 1
      NTOT = 7
      FILEEXP=rt//'DRELLYAN/E772/E772_800_WD.csv'
      CALL READDATA3(FILEEXP,NLINES,X1,X2,X3,NTOT)
      DO I=1, NTOT
         E772_800_W_pT(I)    =X1(I)
         E772_800_W_R(I)     =X2(I)
         E772_800_W_Err(I)   =X3(I)
      ENDDO

C-----E605: Cu: RTS= 38.8 GeV
      NLINES = 1
      NTOT = 74
      FILEEXP=rt//'DRELLYAN/E605/E605_Cu.csv'
      CALL READDATA6(FILEEXP,NLINES,X1,X2,X3,X4,X5,X6,NTOT)
      DO I=1, NTOT
         E605_pT(I)    =X1(I)
         E605_ds(I)    =X2(I)
         E605_err(I)   =X3(I)
         E605_xf(I)    =X4(I)
         E605_Qlo(I)   =X5(I)
         E605_Qhi(I)   =X6(I)
      ENDDO

C-----SIDIS (HERMES)
C-------------------
C-------------------
C-----------IMPORTATION OF ALL HERMES (HERA-DESY) SIDIS DATA
      NTOT = 9 ! 9 ROWS
      NTOT2 = 8! 8 ROWS
      NLINES = 1 !SKIP THE FIRST LINE
C-----------------------------------------
C------PI+ NU DEPENDENT DATA------------------------
C------- HE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pip_he_nu.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PIP_HE_NU_VALUE(I)   =X1(I)
               PIP_HE_NU_Z(I)       =X2(I)
               PIP_HE_NU_Q2(I)      =X3(I)
               PIP_HE_NU_PT2(I)     =X4(I)
               PIP_HE_NU_MULT(I)    =X5(I)
               PIP_HE_NU_STAT(I)    =X6(I)
               PIP_HE_NU_SYS(I)     =X7(I)
               PIP_HE_NU_DIS(I)     =X10(I)
            END DO

C------- NE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pip_ne_nu.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PIP_NE_NU_VALUE(I)   =X1(I)
               PIP_NE_NU_Z(I)       =X2(I)
               PIP_NE_NU_Q2(I)      =X3(I)
               PIP_NE_NU_PT2(I)     =X4(I)
               PIP_NE_NU_MULT(I)    =X5(I)
               PIP_NE_NU_STAT(I)    =X6(I)
               PIP_NE_NU_SYS(I)     =X7(I)
               PIP_NE_NU_DIS(I)     =X10(I)
            END DO
C------- KR
            FILEEXP = RT//'SIDIS/HERMES_DIS/pip_kr_nu.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PIP_KR_NU_VALUE(I)   =X1(I)
               PIP_KR_NU_Z(I)       =X2(I)
               PIP_KR_NU_Q2(I)      =X3(I)
               PIP_KR_NU_PT2(I)     =X4(I)
               PIP_KR_NU_MULT(I)    =X5(I)
               PIP_KR_NU_STAT(I)    =X6(I)
               PIP_KR_NU_SYS(I)     =X7(I)
               PIP_KR_NU_DIS(I)     =X10(I)
            END DO
C------- XE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pip_xe_nu.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PIP_XE_NU_VALUE(I)   =X1(I)
               PIP_XE_NU_Z(I)       =X2(I)
               PIP_XE_NU_Q2(I)      =X3(I)
               PIP_XE_NU_PT2(I)     =X4(I)
               PIP_XE_NU_MULT(I)    =X5(I)
               PIP_XE_NU_STAT(I)    =X6(I)
               PIP_XE_NU_SYS(I)     =X7(I)
               PIP_XE_NU_DIS(I)     =X10(I)
            END DO
C-----------------------------------------
C------PI+ Z DEPENDENT DATA------------------------
C------- HE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pip_he_z.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PIP_HE_Z_VALUE(I)   =X2(I)
               PIP_HE_Z_NU(I)      =X1(I)
               PIP_HE_Z_Q2(I)      =X3(I)
               PIP_HE_Z_PT2(I)     =X4(I)
               PIP_HE_Z_MULT(I)    =X5(I)
               PIP_HE_Z_STAT(I)    =X6(I)
               PIP_HE_Z_SYS(I)     =X7(I)
               PIP_HE_Z_DIS(I)     =X10(I)
            END DO
C------- NE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pip_ne_z.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PIP_NE_Z_VALUE(I)   =X2(I)
               PIP_NE_Z_NU(I)      =X1(I)
               PIP_NE_Z_Q2(I)      =X3(I)
               PIP_NE_Z_PT2(I)     =X4(I)
               PIP_NE_Z_MULT(I)    =X5(I)
               PIP_NE_Z_STAT(I)    =X6(I)
               PIP_NE_Z_SYS(I)     =X7(I)
               PIP_NE_Z_DIS(I)     =X10(I)
            END DO
C------- KR
            FILEEXP = RT//'SIDIS/HERMES_DIS/pip_kr_z.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PIP_KR_Z_VALUE(I)   =X2(I)
               PIP_KR_Z_NU(I)      =X1(I)
               PIP_KR_Z_Q2(I)      =X3(I)
               PIP_KR_Z_PT2(I)     =X4(I)
               PIP_KR_Z_MULT(I)    =X5(I)
               PIP_KR_Z_STAT(I)    =X6(I)
               PIP_KR_Z_SYS(I)     =X7(I)
               PIP_KR_Z_DIS(I)     =X10(I)
            END DO
C------- XE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pip_xe_z.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PIP_XE_Z_VALUE(I)   =X2(I)
               PIP_XE_Z_NU(I)      =X1(I)
               PIP_XE_Z_Q2(I)      =X3(I)
               PIP_XE_Z_PT2(I)     =X4(I)
               PIP_XE_Z_MULT(I)    =X5(I)
               PIP_XE_Z_STAT(I)    =X6(I)
               PIP_XE_Z_SYS(I)     =X7(I)
               PIP_XE_Z_DIS(I)     =X10(I)
            END DO
C-----------------------------------------
C------PI+ Q2 DEPENDENT DATA------------------------
C------- HE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pip_he_q2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT2)
            DO I=1,NTOT2
               PIP_HE_Q2_VALUE(I)   =X3(I)
               PIP_HE_Q2_NU(I)      =X1(I)
               PIP_HE_Q2_Z(I)       =X2(I)
               PIP_HE_Q2_PT2(I)     =X4(I)
               PIP_HE_Q2_MULT(I)    =X5(I)
               PIP_HE_Q2_STAT(I)    =X6(I)
               PIP_HE_Q2_SYS(I)     =X7(I)
               PIP_HE_Q2_DIS(I)     =X10(I)
            END DO
C------- NE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pip_ne_q2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT2)
            DO I=1,NTOT2
               PIP_NE_Q2_VALUE(I)   =X3(I)
               PIP_NE_Q2_NU(I)      =X1(I)
               PIP_NE_Q2_Z(I)       =X2(I)
               PIP_NE_Q2_PT2(I)     =X4(I)
               PIP_NE_Q2_MULT(I)    =X5(I)
               PIP_NE_Q2_STAT(I)    =X6(I)
               PIP_NE_Q2_SYS(I)     =X7(I)
               PIP_NE_Q2_DIS(I)     =X10(I)
            END DO
C------- KR
            FILEEXP = RT//'SIDIS/HERMES_DIS/pip_kr_q2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT2)
            DO I=1,NTOT2
               PIP_KR_Q2_VALUE(I)   =X3(I)
               PIP_KR_Q2_NU(I)      =X1(I)
               PIP_KR_Q2_Z(I)       =X2(I)
               PIP_KR_Q2_PT2(I)     =X4(I)
               PIP_KR_Q2_MULT(I)    =X5(I)
               PIP_KR_Q2_STAT(I)    =X6(I)
               PIP_KR_Q2_SYS(I)     =X7(I)
               PIP_KR_Q2_DIS(I)     =X10(I)
            END DO

C------- XE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pip_xe_q2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT2)
            DO I=1,NTOT2
               PIP_XE_Q2_VALUE(I)   =X3(I)
               PIP_XE_Q2_NU(I)      =X1(I)
               PIP_XE_Q2_Z(I)       =X2(I)
               PIP_XE_Q2_PT2(I)     =X4(I)
               PIP_XE_Q2_MULT(I)    =X5(I)
               PIP_XE_Q2_STAT(I)    =X6(I)
               PIP_XE_Q2_SYS(I)     =X7(I)
               PIP_XE_Q2_DIS(I)     =X10(I)
            END DO
C-----------------------------------------
C-----------------------------------------
C-----------------------------------------
C------PI- NU DEPENDENT DATA------------------------
C------- HE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pim_he_nu.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PIM_HE_NU_VALUE(I)   =X1(I)
               PIM_HE_NU_Z(I)       =X2(I)
               PIM_HE_NU_Q2(I)      =X3(I)
               PIM_HE_NU_PT2(I)     =X4(I)
               PIM_HE_NU_MULT(I)    =X5(I)
               PIM_HE_NU_STAT(I)    =X6(I)
               PIM_HE_NU_SYS(I)     =X7(I)
               PIM_HE_NU_DIS(I)     =X10(I)
            END DO
C------- NE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pim_ne_nu.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PIM_NE_NU_VALUE(I)   =X1(I)
               PIM_NE_NU_Z(I)       =X2(I)
               PIM_NE_NU_Q2(I)      =X3(I)
               PIM_NE_NU_PT2(I)     =X4(I)
               PIM_NE_NU_MULT(I)    =X5(I)
               PIM_NE_NU_STAT(I)    =X6(I)
               PIM_NE_NU_SYS(I)     =X7(I)
               PIM_NE_NU_DIS(I)     =X10(I)
            END DO

C------- KR
            FILEEXP = RT//'SIDIS/HERMES_DIS/pim_kr_nu.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PIM_KR_NU_VALUE(I)   =X1(I)
               PIM_KR_NU_Z(I)       =X2(I)
               PIM_KR_NU_Q2(I)      =X3(I)
               PIM_KR_NU_PT2(I)     =X4(I)
               PIM_KR_NU_MULT(I)    =X5(I)
               PIM_KR_NU_STAT(I)    =X6(I)
               PIM_KR_NU_SYS(I)     =X7(I)
               PIM_KR_NU_DIS(I)     =X10(I)
            END DO

C------- XE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pim_xe_nu.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PIM_XE_NU_VALUE(I)   =X1(I)
               PIM_XE_NU_Z(I)       =X2(I)
               PIM_XE_NU_Q2(I)      =X3(I)
               PIM_XE_NU_PT2(I)     =X4(I)
               PIM_XE_NU_MULT(I)    =X5(I)
               PIM_XE_NU_STAT(I)    =X6(I)
               PIM_XE_NU_SYS(I)     =X7(I)
               PIM_XE_NU_DIS(I)     =X10(I)
            END DO

C-----------------------------------------
C------PI- Z DEPENDENT DATA------------------------
C------- HE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pim_he_z.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PIM_HE_Z_VALUE(I)   =X2(I)
               PIM_HE_Z_NU(I)      =X1(I)
               PIM_HE_Z_Q2(I)      =X3(I)
               PIM_HE_Z_PT2(I)     =X4(I)
               PIM_HE_Z_MULT(I)    =X5(I)
               PIM_HE_Z_STAT(I)    =X6(I)
               PIM_HE_Z_SYS(I)     =X7(I)
               PIM_HE_Z_DIS(I)     =X10(I)
            END DO

C------- NE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pim_ne_z.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PIM_NE_Z_VALUE(I)   =X2(I)
               PIM_NE_Z_NU(I)      =X1(I)
               PIM_NE_Z_Q2(I)      =X3(I)
               PIM_NE_Z_PT2(I)     =X4(I)
               PIM_NE_Z_MULT(I)    =X5(I)
               PIM_NE_Z_STAT(I)    =X6(I)
               PIM_NE_Z_SYS(I)     =X7(I)
               PIM_NE_Z_DIS(I)     =X10(I)
            END DO

C------- KR
            FILEEXP = RT//'SIDIS/HERMES_DIS/pim_kr_z.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PIM_KR_Z_VALUE(I)   =X2(I)
               PIM_KR_Z_NU(I)      =X1(I)
               PIM_KR_Z_Q2(I)      =X3(I)
               PIM_KR_Z_PT2(I)     =X4(I)
               PIM_KR_Z_MULT(I)    =X5(I)
               PIM_KR_Z_STAT(I)    =X6(I)
               PIM_KR_Z_SYS(I)     =X7(I)
               PIM_KR_Z_DIS(I)     =X10(I)
            END DO

C------- XE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pim_xe_z.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PIM_XE_Z_VALUE(I)   =X2(I)
               PIM_XE_Z_NU(I)      =X1(I)
               PIM_XE_Z_Q2(I)      =X3(I)
               PIM_XE_Z_PT2(I)     =X4(I)
               PIM_XE_Z_MULT(I)    =X5(I)
               PIM_XE_Z_STAT(I)    =X6(I)
               PIM_XE_Z_SYS(I)     =X7(I)
               PIM_XE_Z_DIS(I)     =X10(I)
            END DO

C-----------------------------------------
C------PI- Q2 DEPENDENT DATA------------------------
C------- HE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pim_he_q2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT2)
            DO I=1,NTOT2
               PIM_HE_Q2_VALUE(I)   =X3(I)
               PIM_HE_Q2_NU(I)      =X1(I)
               PIM_HE_Q2_Z(I)       =X2(I)
               PIM_HE_Q2_PT2(I)     =X4(I)
               PIM_HE_Q2_MULT(I)    =X5(I)
               PIM_HE_Q2_STAT(I)    =X6(I)
               PIM_HE_Q2_SYS(I)     =X7(I)
               PIM_HE_Q2_DIS(I)     =X10(I)
            END DO

C------- NE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pim_ne_q2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT2)
            DO I=1,NTOT2
               PIM_NE_Q2_VALUE(I)   =X3(I)
               PIM_NE_Q2_NU(I)      =X1(I)
               PIM_NE_Q2_Z(I)       =X2(I)
               PIM_NE_Q2_PT2(I)     =X4(I)
               PIM_NE_Q2_MULT(I)    =X5(I)
               PIM_NE_Q2_STAT(I)    =X6(I)
               PIM_NE_Q2_SYS(I)     =X7(I)
               PIM_NE_Q2_DIS(I)     =X10(I)
            END DO

C------- KR
            FILEEXP = RT//'SIDIS/HERMES_DIS/pim_kr_q2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT2)
            DO I=1,NTOT2
               PIM_KR_Q2_VALUE(I)   =X3(I)
               PIM_KR_Q2_NU(I)      =X1(I)
               PIM_KR_Q2_Z(I)       =X2(I)
               PIM_KR_Q2_PT2(I)     =X4(I)
               PIM_KR_Q2_MULT(I)    =X5(I)
               PIM_KR_Q2_STAT(I)    =X6(I)
               PIM_KR_Q2_SYS(I)     =X7(I)
               PIM_KR_Q2_DIS(I)     =X10(I)
            END DO

C------- XE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pim_xe_q2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT2)
            DO I=1,NTOT2
               PIM_XE_Q2_VALUE(I)   =X3(I)
               PIM_XE_Q2_NU(I)      =X1(I)
               PIM_XE_Q2_Z(I)       =X2(I)
               PIM_XE_Q2_PT2(I)     =X4(I)
               PIM_XE_Q2_MULT(I)    =X5(I)
               PIM_XE_Q2_STAT(I)    =X6(I)
               PIM_XE_Q2_SYS(I)     =X7(I)
               PIM_XE_Q2_DIS(I)     =X10(I)
            END DO

C-----------------------------------------
C------PI+ PT2 DATA-----------------------
C------- HE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pip_he_pt2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PIP_HE_PT2_VALUE(I)   =X4(I)
               PIP_HE_PT2_NU(I)      =X1(I)
               PIP_HE_PT2_Z(I)       =X2(I)
               PIP_HE_PT2_Q2(I)      =X3(I)
               PIP_HE_PT2_MULT(I)    =X5(I)
               PIP_HE_PT2_STAT(I)    =X6(I)
               PIP_HE_PT2_SYS(I)     =X7(I)
               PIP_HE_PT2_DIS(I)     =X10(I)
            END DO


C------- NE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pip_ne_pt2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PIP_NE_PT2_VALUE(I)   =X4(I)
               PIP_NE_PT2_NU(I)      =X1(I)
               PIP_NE_PT2_Z(I)       =X2(I)
               PIP_NE_PT2_Q2(I)      =X3(I)
               PIP_NE_PT2_MULT(I)    =X5(I)
               PIP_NE_PT2_STAT(I)    =X6(I)
               PIP_NE_PT2_SYS(I)     =X7(I)
               PIP_NE_PT2_DIS(I)     =X10(I)
            END DO

C------- KR
            FILEEXP = RT//'SIDIS/HERMES_DIS/pip_kr_pt2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PIP_KR_PT2_VALUE(I)   =X4(I)
               PIP_KR_PT2_NU(I)      =X1(I)
               PIP_KR_PT2_Z(I)       =X2(I)
               PIP_KR_PT2_Q2(I)      =X3(I)
               PIP_KR_PT2_MULT(I)    =X5(I)
               PIP_KR_PT2_STAT(I)    =X6(I)
               PIP_KR_PT2_SYS(I)     =X7(I)
               PIP_KR_PT2_DIS(I)     =X10(I)
            END DO

C------- XE

            FILEEXP = RT//'SIDIS/HERMES_DIS/pip_xe_pt2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PIP_XE_PT2_VALUE(I)   =X4(I)
               PIP_XE_PT2_NU(I)      =X1(I)
               PIP_XE_PT2_Z(I)       =X2(I)
               PIP_XE_PT2_Q2(I)      =X3(I)
               PIP_XE_PT2_MULT(I)    =X5(I)
               PIP_XE_PT2_STAT(I)    =X6(I)
               PIP_XE_PT2_SYS(I)     =X7(I)
               PIP_XE_PT2_DIS(I)     =X10(I)
            END DO
C-----------------------------------------
C------PI- PT2 DATA-----------------------
C------- HE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pim_he_pt2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PIM_HE_PT2_VALUE(I)   =X4(I)
               PIM_HE_PT2_NU(I)      =X1(I)
               PIM_HE_PT2_Z(I)       =X2(I)
               PIM_HE_PT2_Q2(I)      =X3(I)
               PIM_HE_PT2_MULT(I)    =X5(I)
               PIM_HE_PT2_STAT(I)    =X6(I)
               PIM_HE_PT2_SYS(I)     =X7(I)
               PIM_HE_PT2_DIS(I)     =X10(I)
            END DO

C------- NE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pim_ne_pt2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PIM_NE_PT2_VALUE(I)   =X4(I)
               PIM_NE_PT2_NU(I)      =X1(I)
               PIM_NE_PT2_Z(I)       =X2(I)
               PIM_NE_PT2_Q2(I)      =X3(I)
               PIM_NE_PT2_MULT(I)    =X5(I)
               PIM_NE_PT2_STAT(I)    =X6(I)
               PIM_NE_PT2_SYS(I)     =X7(I)
               PIM_NE_PT2_DIS(I)     =X10(I)
            END DO
C------- KR
            FILEEXP = RT//'SIDIS/HERMES_DIS/pim_kr_pt2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PIM_KR_PT2_VALUE(I)   =X4(I)
               PIM_KR_PT2_NU(I)      =X1(I)
               PIM_KR_PT2_Z(I)       =X2(I)
               PIM_KR_PT2_Q2(I)      =X3(I)
               PIM_KR_PT2_MULT(I)    =X5(I)
               PIM_KR_PT2_STAT(I)    =X6(I)
               PIM_KR_PT2_SYS(I)     =X7(I)
               PIM_KR_PT2_DIS(I)     =X10(I)
            END DO
C------- XE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pim_xe_pt2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PIM_XE_PT2_VALUE(I)   =X4(I)
               PIM_XE_PT2_NU(I)      =X1(I)
               PIM_XE_PT2_Z(I)       =X2(I)
               PIM_XE_PT2_Q2(I)      =X3(I)
               PIM_XE_PT2_MULT(I)    =X5(I)
               PIM_XE_PT2_STAT(I)    =X6(I)
               PIM_XE_PT2_SYS(I)     =X7(I)
               PIM_XE_PT2_DIS(I)     =X10(I)
            END DO

C-----------------------------------------
C------PI0 NU DEPENDENT DATA------------------------
C------- HE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pi0_he_nu.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PI0_HE_NU_VALUE(I)   =X1(I)
               PI0_HE_NU_Z(I)       =X2(I)
               PI0_HE_NU_Q2(I)      =X3(I)
               PI0_HE_NU_PT2(I)     =X4(I)
               PI0_HE_NU_MULT(I)    =X5(I)
               PI0_HE_NU_STAT(I)    =X6(I)
               PI0_HE_NU_SYS(I)     =X7(I)
               PI0_HE_NU_DIS(I)     =X10(I)
            END DO

C------- NE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pi0_ne_nu.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PI0_NE_NU_VALUE(I)   =X1(I)
               PI0_NE_NU_Z(I)       =X2(I)
               PI0_NE_NU_Q2(I)      =X3(I)
               PI0_NE_NU_PT2(I)     =X4(I)
               PI0_NE_NU_MULT(I)    =X5(I)
               PI0_NE_NU_STAT(I)    =X6(I)
               PI0_NE_NU_SYS(I)     =X7(I)
               PI0_NE_NU_DIS(I)     =X10(I)
            END DO
C------- KR
            FILEEXP = RT//'SIDIS/HERMES_DIS/pi0_kr_nu.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PI0_KR_NU_VALUE(I)   =X1(I)
               PI0_KR_NU_Z(I)       =X2(I)
               PI0_KR_NU_Q2(I)      =X3(I)
               PI0_KR_NU_PT2(I)     =X4(I)
               PI0_KR_NU_MULT(I)    =X5(I)
               PI0_KR_NU_STAT(I)    =X6(I)
               PI0_KR_NU_SYS(I)     =X7(I)
               PI0_KR_NU_DIS(I)     =X10(I)
            END DO
C------- XE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pi0_xe_nu.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PI0_XE_NU_VALUE(I)   =X1(I)
               PI0_XE_NU_Z(I)       =X2(I)
               PI0_XE_NU_Q2(I)      =X3(I)
               PI0_XE_NU_PT2(I)     =X4(I)
               PI0_XE_NU_MULT(I)    =X5(I)
               PI0_XE_NU_STAT(I)    =X6(I)
               PI0_XE_NU_SYS(I)     =X7(I)
               PI0_XE_NU_DIS(I)     =X10(I)
            END DO
C-----------------------------------------
C------PI0 Z DEPENDENT DATA------------------------
C------- HE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pi0_he_z.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PI0_HE_Z_VALUE(I)   =X2(I)
               PI0_HE_Z_NU(I)      =X1(I)
               PI0_HE_Z_Q2(I)      =X3(I)
               PI0_HE_Z_PT2(I)     =X4(I)
               PI0_HE_Z_MULT(I)    =X5(I)
               PI0_HE_Z_STAT(I)    =X6(I)
               PI0_HE_Z_SYS(I)     =X7(I)
               PI0_HE_Z_DIS(I)     =X10(I)
            END DO

C------- NE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pi0_ne_z.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PI0_NE_Z_VALUE(I)   =X2(I)
               PI0_NE_Z_NU(I)      =X1(I)
               PI0_NE_Z_Q2(I)      =X3(I)
               PI0_NE_Z_PT2(I)     =X4(I)
               PI0_NE_Z_MULT(I)    =X5(I)
               PI0_NE_Z_STAT(I)    =X6(I)
               PI0_NE_Z_SYS(I)     =X7(I)
               PI0_NE_Z_DIS(I)     =X10(I)
            END DO

C------- KR
            FILEEXP = RT//'SIDIS/HERMES_DIS/pi0_kr_z.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PI0_KR_Z_VALUE(I)   =X2(I)
               PI0_KR_Z_NU(I)      =X1(I)
               PI0_KR_Z_Q2(I)      =X3(I)
               PI0_KR_Z_PT2(I)     =X4(I)
               PI0_KR_Z_MULT(I)    =X5(I)
               PI0_KR_Z_STAT(I)    =X6(I)
               PI0_KR_Z_SYS(I)     =X7(I)
               PI0_KR_Z_DIS(I)     =X10(I)
            END DO

C------- XE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pi0_xe_z.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PI0_XE_Z_VALUE(I)   =X2(I)
               PI0_XE_Z_NU(I)      =X1(I)
               PI0_XE_Z_Q2(I)      =X3(I)
               PI0_XE_Z_PT2(I)     =X4(I)
               PI0_XE_Z_MULT(I)    =X5(I)
               PI0_XE_Z_STAT(I)    =X6(I)
               PI0_XE_Z_SYS(I)     =X7(I)
               PI0_XE_Z_DIS(I)     =X10(I)
            END DO

C-----------------------------------------
C------PI0 Q2 DEPENDENT DATA------------------------
C------- HE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pi0_he_q2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT2)
            DO I=1,NTOT2
               PI0_HE_Q2_VALUE(I)   =X3(I)
               PI0_HE_Q2_NU(I)      =X1(I)
               PI0_HE_Q2_Z(I)       =X2(I)
               PI0_HE_Q2_PT2(I)     =X4(I)
               PI0_HE_Q2_MULT(I)    =X5(I)
               PI0_HE_Q2_STAT(I)    =X6(I)
               PI0_HE_Q2_SYS(I)     =X7(I)
               PI0_HE_Q2_DIS(I)     =X10(I)
            END DO

C------- NE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pi0_ne_q2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT2)
            DO I=1,NTOT2
               PI0_NE_Q2_VALUE(I)   =X3(I)
               PI0_NE_Q2_NU(I)      =X1(I)
               PI0_NE_Q2_Z(I)       =X2(I)
               PI0_NE_Q2_PT2(I)     =X4(I)
               PI0_NE_Q2_MULT(I)    =X5(I)
               PI0_NE_Q2_STAT(I)    =X6(I)
               PI0_NE_Q2_SYS(I)     =X7(I)
               PI0_NE_Q2_DIS(I)     =X10(I)
            END DO

C------- KR
            FILEEXP = RT//'SIDIS/HERMES_DIS/pi0_kr_q2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT2)
            DO I=1,NTOT2
               PI0_KR_Q2_VALUE(I)   =X3(I)
               PI0_KR_Q2_NU(I)      =X1(I)
               PI0_KR_Q2_Z(I)       =X2(I)
               PI0_KR_Q2_PT2(I)     =X4(I)
               PI0_KR_Q2_MULT(I)    =X5(I)
               PI0_KR_Q2_STAT(I)    =X6(I)
               PI0_KR_Q2_SYS(I)     =X7(I)
               PI0_KR_Q2_DIS(I)     =X10(I)
            END DO

C------- XE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pi0_xe_q2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT2)
            DO I=1,NTOT2
               PI0_XE_Q2_VALUE(I)   =X3(I)
               PI0_XE_Q2_NU(I)      =X1(I)
               PI0_XE_Q2_Z(I)       =X2(I)
               PI0_XE_Q2_PT2(I)     =X4(I)
               PI0_XE_Q2_MULT(I)    =X5(I)
               PI0_XE_Q2_STAT(I)    =X6(I)
               PI0_XE_Q2_SYS(I)     =X7(I)
               PI0_XE_Q2_DIS(I)     =X10(I)
            END DO

C-----------------------------------------
C------PI0 PT2 DATA-----------------------
C------- HE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pi0_he_pt2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PI0_HE_PT2_VALUE(I)   =X4(I)
               PI0_HE_PT2_NU(I)      =X1(I)
               PI0_HE_PT2_Z(I)       =X2(I)
               PI0_HE_PT2_Q2(I)      =X3(I)
               PI0_HE_PT2_MULT(I)    =X5(I)
               PI0_HE_PT2_STAT(I)    =X6(I)
               PI0_HE_PT2_SYS(I)     =X7(I)
               PI0_HE_PT2_DIS(I)     =X10(I)
            END DO

C------- NE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pi0_ne_pt2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PI0_NE_PT2_VALUE(I)   =X4(I)
               PI0_NE_PT2_NU(I)      =X1(I)
               PI0_NE_PT2_Z(I)       =X2(I)
               PI0_NE_PT2_Q2(I)      =X3(I)
               PI0_NE_PT2_MULT(I)    =X5(I)
               PI0_NE_PT2_STAT(I)    =X6(I)
               PI0_NE_PT2_SYS(I)     =X7(I)
               PI0_NE_PT2_DIS(I)     =X10(I)
            END DO


C------- KR
            FILEEXP = RT//'SIDIS/HERMES_DIS/pi0_kr_pt2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PI0_KR_PT2_VALUE(I)   =X4(I)
               PI0_KR_PT2_NU(I)      =X1(I)
               PI0_KR_PT2_Z(I)       =X2(I)
               PI0_KR_PT2_Q2(I)      =X3(I)
               PI0_KR_PT2_MULT(I)    =X5(I)
               PI0_KR_PT2_STAT(I)    =X6(I)
               PI0_KR_PT2_SYS(I)     =X7(I)
               PI0_KR_PT2_DIS(I)     =X10(I)
            END DO


C------- XE
            FILEEXP = RT//'SIDIS/HERMES_DIS/pi0_xe_pt2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               PI0_XE_PT2_VALUE(I)   =X4(I)
               PI0_XE_PT2_NU(I)      =X1(I)
               PI0_XE_PT2_Z(I)       =X2(I)
               PI0_XE_PT2_Q2(I)      =X3(I)
               PI0_XE_PT2_MULT(I)    =X5(I)
               PI0_XE_PT2_STAT(I)    =X6(I)
               PI0_XE_PT2_SYS(I)     =X7(I)
               PI0_XE_PT2_DIS(I)     =X10(I)
            END DO


C--------------------------------------------------
C------K- NU DEPENDENT DATA------------------------
C------- HE
            FILEEXP = RT//'SIDIS/HERMES_DIS/km_he_nu.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               KM_HE_NU_VALUE(I)   =X1(I)
               KM_HE_NU_Z(I)       =X2(I)
               KM_HE_NU_Q2(I)      =X3(I)
               KM_HE_NU_PT2(I)     =X4(I)
               KM_HE_NU_MULT(I)    =X5(I)
               KM_HE_NU_STAT(I)    =X6(I)
               KM_HE_NU_SYS(I)     =X7(I)
               KM_HE_NU_DIS(I)     =X10(I)
            END DO

C------- NE
            FILEEXP = RT//'SIDIS/HERMES_DIS/km_ne_nu.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               KM_NE_NU_VALUE(I)   =X1(I)
               KM_NE_NU_Z(I)       =X2(I)
               KM_NE_NU_Q2(I)      =X3(I)
               KM_NE_NU_PT2(I)     =X4(I)
               KM_NE_NU_MULT(I)    =X5(I)
               KM_NE_NU_STAT(I)    =X6(I)
               KM_NE_NU_SYS(I)     =X7(I)
               KM_NE_NU_DIS(I)     =X10(I)
            END DO


C------- KR
            FILEEXP = RT//'SIDIS/HERMES_DIS/km_kr_nu.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               KM_KR_NU_VALUE(I)   =X1(I)
               KM_KR_NU_Z(I)       =X2(I)
               KM_KR_NU_Q2(I)      =X3(I)
               KM_KR_NU_PT2(I)     =X4(I)
               KM_KR_NU_MULT(I)    =X5(I)
               KM_KR_NU_STAT(I)    =X6(I)
               KM_KR_NU_SYS(I)     =X7(I)
               KM_KR_NU_DIS(I)     =X10(I)
            END DO


C------- XE
            FILEEXP = RT//'SIDIS/HERMES_DIS/km_xe_nu.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               KM_XE_NU_VALUE(I)   =X1(I)
               KM_XE_NU_Z(I)       =X2(I)
               KM_XE_NU_Q2(I)      =X3(I)
               KM_XE_NU_PT2(I)     =X4(I)
               KM_XE_NU_MULT(I)    =X5(I)
               KM_XE_NU_STAT(I)    =X6(I)
               KM_XE_NU_SYS(I)     =X7(I)
               KM_XE_NU_DIS(I)     =X10(I)
            END DO


C-----------------------------------------
C------K- Z DEPENDENT DATA------------------------
C------- HE
            FILEEXP = RT//'SIDIS/HERMES_DIS/km_he_z.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               KM_HE_Z_VALUE(I)   =X2(I)
               KM_HE_Z_NU(I)      =X1(I)
               KM_HE_Z_Q2(I)      =X3(I)
               KM_HE_Z_PT2(I)     =X4(I)
               KM_HE_Z_MULT(I)    =X5(I)
               KM_HE_Z_STAT(I)    =X6(I)
               KM_HE_Z_SYS(I)     =X7(I)
               KM_HE_Z_DIS(I)     =X10(I)
            END DO


C------- NE
            FILEEXP = RT//'SIDIS/HERMES_DIS/km_ne_z.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               KM_NE_Z_VALUE(I)   =X2(I)
               KM_NE_Z_NU(I)      =X1(I)
               KM_NE_Z_Q2(I)      =X3(I)
               KM_NE_Z_PT2(I)     =X4(I)
               KM_NE_Z_MULT(I)    =X5(I)
               KM_NE_Z_STAT(I)    =X6(I)
               KM_NE_Z_SYS(I)     =X7(I)
               KM_NE_Z_DIS(I)     =X10(I)
            END DO


C------- KR
            FILEEXP = RT//'SIDIS/HERMES_DIS/km_kr_z.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               KM_KR_Z_VALUE(I)   =X2(I)
               KM_KR_Z_NU(I)      =X1(I)
               KM_KR_Z_Q2(I)      =X3(I)
               KM_KR_Z_PT2(I)     =X4(I)
               KM_KR_Z_MULT(I)    =X5(I)
               KM_KR_Z_STAT(I)    =X6(I)
               KM_KR_Z_SYS(I)     =X7(I)
               KM_KR_Z_DIS(I)     =X10(I)
            END DO


C------- XE
            FILEEXP = RT//'SIDIS/HERMES_DIS/km_xe_z.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               KM_XE_Z_VALUE(I)   =X2(I)
               KM_XE_Z_NU(I)      =X1(I)
               KM_XE_Z_Q2(I)      =X3(I)
               KM_XE_Z_PT2(I)     =X4(I)
               KM_XE_Z_MULT(I)    =X5(I)
               KM_XE_Z_STAT(I)    =X6(I)
               KM_XE_Z_SYS(I)     =X7(I)
               KM_XE_Z_DIS(I)     =X10(I)
            END DO


C-----------------------------------------
C------K- Q2 DEPENDENT DATA------------------------
C------- HE
            FILEEXP = RT//'SIDIS/HERMES_DIS/km_he_q2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT2)
            DO I=1,NTOT2
               KM_HE_Q2_VALUE(I)   =X3(I)
               KM_HE_Q2_NU(I)      =X1(I)
               KM_HE_Q2_Z(I)       =X2(I)
               KM_HE_Q2_PT2(I)     =X4(I)
               KM_HE_Q2_MULT(I)    =X5(I)
               KM_HE_Q2_STAT(I)    =X6(I)
               KM_HE_Q2_SYS(I)     =X7(I)
               KM_HE_Q2_DIS(I)     =X10(I)
            END DO



C------- NE
            FILEEXP = RT//'SIDIS/HERMES_DIS/km_ne_q2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT2)
            DO I=1,NTOT2
               KM_NE_Q2_VALUE(I)   =X3(I)
               KM_NE_Q2_NU(I)      =X1(I)
               KM_NE_Q2_Z(I)       =X2(I)
               KM_NE_Q2_PT2(I)     =X4(I)
               KM_NE_Q2_MULT(I)    =X5(I)
               KM_NE_Q2_STAT(I)    =X6(I)
               KM_NE_Q2_SYS(I)     =X7(I)
               KM_NE_Q2_DIS(I)     =X10(I)
            END DO


C------- KR
            FILEEXP = RT//'SIDIS/HERMES_DIS/km_kr_q2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT2)
            DO I=1,NTOT2
               KM_KR_Q2_VALUE(I)   =X3(I)
               KM_KR_Q2_NU(I)      =X1(I)
               KM_KR_Q2_Z(I)       =X2(I)
               KM_KR_Q2_PT2(I)     =X4(I)
               KM_KR_Q2_MULT(I)    =X5(I)
               KM_KR_Q2_STAT(I)    =X6(I)
               KM_KR_Q2_SYS(I)     =X7(I)
               KM_KR_Q2_DIS(I)     =X10(I)
            END DO

C------- XE
            FILEEXP = RT//'SIDIS/HERMES_DIS/km_xe_q2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT2)
            DO I=1,NTOT2
               KM_XE_Q2_VALUE(I)   =X3(I)
               KM_XE_Q2_NU(I)      =X1(I)
               KM_XE_Q2_Z(I)       =X2(I)
               KM_XE_Q2_PT2(I)     =X4(I)
               KM_XE_Q2_MULT(I)    =X5(I)
               KM_XE_Q2_STAT(I)    =X6(I)
               KM_XE_Q2_SYS(I)     =X7(I)
               KM_XE_Q2_DIS(I)     =X10(I)
            END DO

C------K- PT2 DATA-----------------------
C------- HE
            FILEEXP = RT//'SIDIS/HERMES_DIS/km_he_pt2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               KM_HE_PT2_VALUE(I)   =X4(I)
               KM_HE_PT2_NU(I)      =X1(I)
               KM_HE_PT2_Z(I)       =X2(I)
               KM_HE_PT2_Q2(I)      =X3(I)
               KM_HE_PT2_MULT(I)    =X5(I)
               KM_HE_PT2_STAT(I)    =X6(I)
               KM_HE_PT2_SYS(I)     =X7(I)
               KM_HE_PT2_DIS(I)     =X10(I)
            END DO

C------- NE
            FILEEXP = RT//'SIDIS/HERMES_DIS/km_ne_pt2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               KM_NE_PT2_VALUE(I)   =X4(I)
               KM_NE_PT2_NU(I)      =X1(I)
               KM_NE_PT2_Z(I)       =X2(I)
               KM_NE_PT2_Q2(I)      =X3(I)
               KM_NE_PT2_MULT(I)    =X5(I)
               KM_NE_PT2_STAT(I)    =X6(I)
               KM_NE_PT2_SYS(I)     =X7(I)
               KM_NE_PT2_DIS(I)     =X10(I)
            END DO
C------- KR
            FILEEXP = RT//'SIDIS/HERMES_DIS/km_kr_pt2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               KM_KR_PT2_VALUE(I)   =X4(I)
               KM_KR_PT2_NU(I)      =X1(I)
               KM_KR_PT2_Z(I)       =X2(I)
               KM_KR_PT2_Q2(I)      =X3(I)
               KM_KR_PT2_MULT(I)    =X5(I)
               KM_KR_PT2_STAT(I)    =X6(I)
               KM_KR_PT2_SYS(I)     =X7(I)
               KM_KR_PT2_DIS(I)     =X10(I)
            END DO
C------- XE
            FILEEXP = RT//'SIDIS/HERMES_DIS/km_xe_pt2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               KM_XE_PT2_VALUE(I)   =X4(I)
               KM_XE_PT2_NU(I)      =X1(I)
               KM_XE_PT2_Z(I)       =X2(I)
               KM_XE_PT2_Q2(I)      =X3(I)
               KM_XE_PT2_MULT(I)    =X5(I)
               KM_XE_PT2_STAT(I)    =X6(I)
               KM_XE_PT2_SYS(I)     =X7(I)
               KM_XE_PT2_DIS(I)     =X10(I)
            END DO


C-------------------------------------------------
C------K+ NU DEPENDENT DATA------------------------
C------- HE
            FILEEXP = RT//'SIDIS/HERMES_DIS/kp_he_nu.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               KP_HE_NU_VALUE(I)   =X1(I)
               KP_HE_NU_Z(I)       =X2(I)
               KP_HE_NU_Q2(I)      =X3(I)
               KP_HE_NU_PT2(I)     =X4(I)
               KP_HE_NU_MULT(I)    =X5(I)
               KP_HE_NU_STAT(I)    =X6(I)
               KP_HE_NU_SYS(I)     =X7(I)
               KP_HE_NU_DIS(I)     =X10(I)
            END DO

C------- NE
            FILEEXP = RT//'SIDIS/HERMES_DIS/kp_ne_nu.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               KP_NE_NU_VALUE(I)   =X1(I)
               KP_NE_NU_Z(I)       =X2(I)
               KP_NE_NU_Q2(I)      =X3(I)
               KP_NE_NU_PT2(I)     =X4(I)
               KP_NE_NU_MULT(I)    =X5(I)
               KP_NE_NU_STAT(I)    =X6(I)
               KP_NE_NU_SYS(I)     =X7(I)
               KP_NE_NU_DIS(I)     =X10(I)
            END DO
C------- KR
            FILEEXP = RT//'SIDIS/HERMES_DIS/kp_kr_nu.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               KP_KR_NU_VALUE(I)   =X1(I)
               KP_KR_NU_Z(I)       =X2(I)
               KP_KR_NU_Q2(I)      =X3(I)
               KP_KR_NU_PT2(I)     =X4(I)
               KP_KR_NU_MULT(I)    =X5(I)
               KP_KR_NU_STAT(I)    =X6(I)
               KP_KR_NU_SYS(I)     =X7(I)
               KP_KR_NU_DIS(I)     =X10(I)
            END DO
C------- XE
            FILEEXP = RT//'SIDIS/HERMES_DIS/kp_xe_nu.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               KP_XE_NU_VALUE(I)   =X1(I)
               KP_XE_NU_Z(I)       =X2(I)
               KP_XE_NU_Q2(I)      =X3(I)
               KP_XE_NU_PT2(I)     =X4(I)
               KP_XE_NU_MULT(I)    =X5(I)
               KP_XE_NU_STAT(I)    =X6(I)
               KP_XE_NU_SYS(I)     =X7(I)
               KP_XE_NU_DIS(I)     =X10(I)
            END DO
C-----------------------------------------
C------K+ Z DEPENDENT DATA------------------------
C------- HE
            FILEEXP = RT//'SIDIS/HERMES_DIS/kp_he_z.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               KP_HE_Z_VALUE(I)   =X2(I)
               KP_HE_Z_NU(I)      =X1(I)
               KP_HE_Z_Q2(I)      =X3(I)
               KP_HE_Z_PT2(I)     =X4(I)
               KP_HE_Z_MULT(I)    =X5(I)
               KP_HE_Z_STAT(I)    =X6(I)
               KP_HE_Z_SYS(I)     =X7(I)
               KP_HE_Z_DIS(I)     =X10(I)
            END DO
C------- NE
            FILEEXP = RT//'SIDIS/HERMES_DIS/kp_ne_z.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               KP_NE_Z_VALUE(I)   =X2(I)
               KP_NE_Z_NU(I)      =X1(I)
               KP_NE_Z_Q2(I)      =X3(I)
               KP_NE_Z_PT2(I)     =X4(I)
               KP_NE_Z_MULT(I)    =X5(I)
               KP_NE_Z_STAT(I)    =X6(I)
               KP_NE_Z_SYS(I)     =X7(I)
               KP_NE_Z_DIS(I)     =X10(I)
            END DO
C------- KR
            FILEEXP = RT//'SIDIS/HERMES_DIS/kp_kr_z.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               KP_KR_Z_VALUE(I)   =X2(I)
               KP_KR_Z_NU(I)      =X1(I)
               KP_KR_Z_Q2(I)      =X3(I)
               KP_KR_Z_PT2(I)     =X4(I)
               KP_KR_Z_MULT(I)    =X5(I)
               KP_KR_Z_STAT(I)    =X6(I)
               KP_KR_Z_SYS(I)     =X7(I)
               KP_KR_Z_DIS(I)     =X10(I)
            END DO
C------- XE
            FILEEXP = RT//'SIDIS/HERMES_DIS/kp_xe_z.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               KP_XE_Z_VALUE(I)   =X2(I)
               KP_XE_Z_NU(I)      =X1(I)
               KP_XE_Z_Q2(I)      =X3(I)
               KP_XE_Z_PT2(I)     =X4(I)
               KP_XE_Z_MULT(I)    =X5(I)
               KP_XE_Z_STAT(I)    =X6(I)
               KP_XE_Z_SYS(I)     =X7(I)
               KP_XE_Z_DIS(I)     =X10(I)
            END DO
C-----------------------------------------
C------K+ Q2 DEPENDENT DATA------------------------
C------- HE
            FILEEXP = RT//'SIDIS/HERMES_DIS/kp_he_q2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT2)
            DO I=1,NTOT2
               KP_HE_Q2_VALUE(I)   =X3(I)
               KP_HE_Q2_NU(I)      =X1(I)
               KP_HE_Q2_Z(I)       =X2(I)
               KP_HE_Q2_PT2(I)     =X4(I)
               KP_HE_Q2_MULT(I)    =X5(I)
               KP_HE_Q2_STAT(I)    =X6(I)
               KP_HE_Q2_SYS(I)     =X7(I)
               KP_HE_Q2_DIS(I)     =X10(I)
            END DO
C------- NE
            FILEEXP = RT//'SIDIS/HERMES_DIS/kp_ne_q2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT2)
            DO I=1,NTOT2
               KP_NE_Q2_VALUE(I)   =X3(I)
               KP_NE_Q2_NU(I)      =X1(I)
               KP_NE_Q2_Z(I)       =X2(I)
               KP_NE_Q2_PT2(I)     =X4(I)
               KP_NE_Q2_MULT(I)    =X5(I)
               KP_NE_Q2_STAT(I)    =X6(I)
               KP_NE_Q2_SYS(I)     =X7(I)
               KP_NE_Q2_DIS(I)     =X10(I)
            END DO
C------- KR
            FILEEXP = RT//'SIDIS/HERMES_DIS/kp_kr_q2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT2)
            DO I=1,NTOT2
               KP_KR_Q2_VALUE(I)   =X3(I)
               KP_KR_Q2_NU(I)      =X1(I)
               KP_KR_Q2_Z(I)       =X2(I)
               KP_KR_Q2_PT2(I)     =X4(I)
               KP_KR_Q2_MULT(I)    =X5(I)
               KP_KR_Q2_STAT(I)    =X6(I)
               KP_KR_Q2_SYS(I)     =X7(I)
               KP_KR_Q2_DIS(I)     =X10(I)
            END DO

C------- XE
            FILEEXP = RT//'SIDIS/HERMES_DIS/kp_xe_q2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT2)
            DO I=1,NTOT2
               KP_XE_Q2_VALUE(I)   =X3(I)
               KP_XE_Q2_NU(I)      =X1(I)
               KP_XE_Q2_Z(I)       =X2(I)
               KP_XE_Q2_PT2(I)     =X4(I)
               KP_XE_Q2_MULT(I)    =X5(I)
               KP_XE_Q2_STAT(I)    =X6(I)
               KP_XE_Q2_SYS(I)     =X7(I)
               KP_XE_Q2_DIS(I)     =X10(I)
            END DO

C------K+ PT2 DATA-----------------------
C------- HE
            FILEEXP = RT//'SIDIS/HERMES_DIS/kp_he_pt2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               KP_HE_PT2_VALUE(I)   =X4(I)
               KP_HE_PT2_NU(I)      =X1(I)
               KP_HE_PT2_Z(I)       =X2(I)
               KP_HE_PT2_Q2(I)      =X3(I)
               KP_HE_PT2_MULT(I)    =X5(I)
               KP_HE_PT2_STAT(I)    =X6(I)
               KP_HE_PT2_SYS(I)     =X7(I)
               KP_HE_PT2_DIS(I)     =X10(I)
            END DO
C------- NE
            FILEEXP = RT//'SIDIS/HERMES_DIS/kp_ne_pt2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               KP_NE_PT2_VALUE(I)   =X4(I)
               KP_NE_PT2_NU(I)      =X1(I)
               KP_NE_PT2_Z(I)       =X2(I)
               KP_NE_PT2_Q2(I)      =X3(I)
               KP_NE_PT2_MULT(I)    =X5(I)
               KP_NE_PT2_STAT(I)    =X6(I)
               KP_NE_PT2_SYS(I)     =X7(I)
               KP_NE_PT2_DIS(I)     =X10(I)
            END DO

C------- KR
            FILEEXP = RT//'SIDIS/HERMES_DIS/kp_kr_pt2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               KP_KR_PT2_VALUE(I)   =X4(I)
               KP_KR_PT2_NU(I)      =X1(I)
               KP_KR_PT2_Z(I)       =X2(I)
               KP_KR_PT2_Q2(I)      =X3(I)
               KP_KR_PT2_MULT(I)    =X5(I)
               KP_KR_PT2_STAT(I)    =X6(I)
               KP_KR_PT2_SYS(I)     =X7(I)
               KP_KR_PT2_DIS(I)     =X10(I)
            END DO
C------- XE
            FILEEXP = RT//'SIDIS/HERMES_DIS/kp_xe_pt2.dat'
            CALL READDATA10(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, NTOT)
            DO I=1,NTOT
               KP_XE_PT2_VALUE(I)   =X4(I)
               KP_XE_PT2_NU(I)      =X1(I)
               KP_XE_PT2_Z(I)       =X2(I)
               KP_XE_PT2_Q2(I)      =X3(I)
               KP_XE_PT2_MULT(I)    =X5(I)
               KP_XE_PT2_STAT(I)    =X6(I)
               KP_XE_PT2_SYS(I)     =X7(I)
               KP_XE_PT2_DIS(I)     =X10(I)
            END DO

C-----------------------------
C------- JLAB
      NTOT = 22

            FILEEXP = RT//'SIDIS/JLAB_DIS/Clas_bin1.dat'
            CALL READDATA19(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,  X7,  X8,  X9,  X10, X11, X12, X13,
     *                      X14, X15, X16, X17, X18, X19, NTOT)
            DO I=1,NTOT
               CLAS1_PTLOW(I)       =X1(I)
               CLAS1_PTHIGH(I)      =X2(I)
               CLAS1_NULOW(I)       =X3(I)
               CLAS1_NUHIGH(I)      =X4(I)
               CLAS1_Q2LOW(I)       =X5(I)
               CLAS1_Q2HIGH(I)      =X6(I)
               CLAS1_RA_C(I)        =X7(I)
               CLAS1_RA_Cerr(I)     =X8(I)
               CLAS1_RA_FE(I)       =X9(I)
               CLAS1_RA_FEerr(I)    =X10(I)
               CLAS1_RA_PB(I)       =X11(I)
               CLAS1_RA_PBerr(I)    =X12(I)
               CLAS1_C_DIS(I)       =X17(I)
               CLAS1_FE_DIS(I)      =X18(I)
               CLAS1_PB_DIS(I)      =X19(I)
            END DO

            FILEEXP = RT//'SIDIS/JLAB_DIS/Clas_bin2.dat'
            CALL READDATA19(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, X11, X12, X13,
     *                     X14, X15, X16, X17, X18, X19, NTOT)
            DO I=1,NTOT
               CLAS2_PTLOW(I)       =X1(I)
               CLAS2_PTHIGH(I)      =X2(I)
               CLAS2_NULOW(I)       =X3(I)
               CLAS2_NUHIGH(I)      =X4(I)
               CLAS2_Q2LOW(I)       =X5(I)
               CLAS2_Q2HIGH(I)      =X6(I)
               CLAS2_RA_C(I)        =X7(I)
               CLAS2_RA_Cerr(I)     =X8(I)
               CLAS2_RA_FE(I)       =X9(I)
               CLAS2_RA_FEerr(I)    =X10(I)
               CLAS2_RA_PB(I)       =X11(I)
               CLAS2_RA_PBerr(I)    =X12(I)
               CLAS2_C_DIS(I)       =X17(I)
               CLAS2_FE_DIS(I)      =X18(I)
               CLAS2_PB_DIS(I)      =X19(I)
            END DO

            FILEEXP = RT//'SIDIS/JLAB_DIS/Clas_bin3.dat'
            CALL READDATA19(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, X11, X12, X13,
     *                     X14, X15, X16, X17, X18, X19, NTOT)
            DO I=1,NTOT
               CLAS3_PTLOW(I)       =X1(I)
               CLAS3_PTHIGH(I)      =X2(I)
               CLAS3_NULOW(I)       =X3(I)
               CLAS3_NUHIGH(I)      =X4(I)
               CLAS3_Q2LOW(I)       =X5(I)
               CLAS3_Q2HIGH(I)      =X6(I)
               CLAS3_RA_C(I)        =X7(I)
               CLAS3_RA_Cerr(I)     =X8(I)
               CLAS3_RA_FE(I)       =X9(I)
               CLAS3_RA_FEerr(I)    =X10(I)
               CLAS3_RA_PB(I)       =X11(I)
               CLAS3_RA_PBerr(I)    =X12(I)
               CLAS3_C_DIS(I)       =X17(I)
               CLAS3_FE_DIS(I)      =X18(I)
               CLAS3_PB_DIS(I)      =X19(I)
            END DO

            FILEEXP = RT//'SIDIS/JLAB_DIS/Clas_bin4.dat'
            CALL READDATA19(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, X11, X12, X13,
     *                     X14, X15, X16, X17, X18, X19, NTOT)
            DO I=1,NTOT
               CLAS4_PTLOW(I)       =X1(I)
               CLAS4_PTHIGH(I)      =X2(I)
               CLAS4_NULOW(I)       =X3(I)
               CLAS4_NUHIGH(I)      =X4(I)
               CLAS4_Q2LOW(I)       =X5(I)
               CLAS4_Q2HIGH(I)      =X6(I)
               CLAS4_RA_C(I)        =X7(I)
               CLAS4_RA_Cerr(I)     =X8(I)
               CLAS4_RA_FE(I)       =X9(I)
               CLAS4_RA_FEerr(I)    =X10(I)
               CLAS4_RA_PB(I)       =X11(I)
               CLAS4_RA_PBerr(I)    =X12(I)
               CLAS4_C_DIS(I)       =X17(I)
               CLAS4_FE_DIS(I)      =X18(I)
               CLAS4_PB_DIS(I)      =X19(I)
            END DO

            FILEEXP = RT//'SIDIS/JLAB_DIS/Clas_bin5.dat'
            CALL READDATA19(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, X11, X12, X13,
     *                     X14, X15, X16, X17, X18, X19, NTOT)
            DO I=1,NTOT
               CLAS5_PTLOW(I)       =X1(I)
               CLAS5_PTHIGH(I)      =X2(I)
               CLAS5_NULOW(I)       =X3(I)
               CLAS5_NUHIGH(I)      =X4(I)
               CLAS5_Q2LOW(I)       =X5(I)
               CLAS5_Q2HIGH(I)      =X6(I)
               CLAS5_RA_C(I)        =X7(I)
               CLAS5_RA_Cerr(I)     =X8(I)
               CLAS5_RA_FE(I)       =X9(I)
               CLAS5_RA_FEerr(I)    =X10(I)
               CLAS5_RA_PB(I)       =X11(I)
               CLAS5_RA_PBerr(I)    =X12(I)
               CLAS5_C_DIS(I)       =X17(I)
               CLAS5_FE_DIS(I)      =X18(I)
               CLAS5_PB_DIS(I)      =X19(I)
            END DO

            FILEEXP = RT//'SIDIS/JLAB_DIS/Clas_bin6.dat'
            CALL READDATA19(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, X11, X12, X13,
     *                     X14, X15, X16, X17, X18, X19, NTOT)
            DO I=1,NTOT
               CLAS6_PTLOW(I)       =X1(I)
               CLAS6_PTHIGH(I)      =X2(I)
               CLAS6_NULOW(I)       =X3(I)
               CLAS6_NUHIGH(I)      =X4(I)
               CLAS6_Q2LOW(I)       =X5(I)
               CLAS6_Q2HIGH(I)      =X6(I)
               CLAS6_RA_C(I)        =X7(I)
               CLAS6_RA_Cerr(I)     =X8(I)
               CLAS6_RA_FE(I)       =X9(I)
               CLAS6_RA_FEerr(I)    =X10(I)
               CLAS6_RA_PB(I)       =X11(I)
               CLAS6_RA_PBerr(I)    =X12(I)
               CLAS6_C_DIS(I)       =X17(I)
               CLAS6_FE_DIS(I)      =X18(I)
               CLAS6_PB_DIS(I)      =X19(I)
            END DO


            FILEEXP = RT//'SIDIS/JLAB_DIS/Clas_bin7.dat'
            CALL READDATA19(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, X11, X12, X13,
     *                     X14, X15, X16, X17, X18, X19, NTOT)
            DO I=1,NTOT
               CLAS7_PTLOW(I)       =X1(I)
               CLAS7_PTHIGH(I)      =X2(I)
               CLAS7_NULOW(I)       =X3(I)
               CLAS7_NUHIGH(I)      =X4(I)
               CLAS7_Q2LOW(I)       =X5(I)
               CLAS7_Q2HIGH(I)      =X6(I)
               CLAS7_RA_C(I)        =X7(I)
               CLAS7_RA_Cerr(I)     =X8(I)
               CLAS7_RA_FE(I)       =X9(I)
               CLAS7_RA_FEerr(I)    =X10(I)
               CLAS7_RA_PB(I)       =X11(I)
               CLAS7_RA_PBerr(I)    =X12(I)
               CLAS7_C_DIS(I)       =X17(I)
               CLAS7_FE_DIS(I)      =X18(I)
               CLAS7_PB_DIS(I)      =X19(I)
            END DO

            FILEEXP = RT//'SIDIS/JLAB_DIS/Clas_bin8.dat'
            CALL READDATA19(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, X11, X12, X13,
     *                     X14, X15, X16, X17, X18, X19, NTOT)
            DO I=1,NTOT
               CLAS8_PTLOW(I)       =X1(I)
               CLAS8_PTHIGH(I)      =X2(I)
               CLAS8_NULOW(I)       =X3(I)
               CLAS8_NUHIGH(I)      =X4(I)
               CLAS8_Q2LOW(I)       =X5(I)
               CLAS8_Q2HIGH(I)      =X6(I)
               CLAS8_RA_C(I)        =X7(I)
               CLAS8_RA_Cerr(I)     =X8(I)
               CLAS8_RA_FE(I)       =X9(I)
               CLAS8_RA_FEerr(I)    =X10(I)
               CLAS8_RA_PB(I)       =X11(I)
               CLAS8_RA_PBerr(I)    =X12(I)
               CLAS8_C_DIS(I)       =X17(I)
               CLAS8_FE_DIS(I)      =X18(I)
               CLAS8_PB_DIS(I)      =X19(I)
            END DO

            FILEEXP = RT//'SIDIS/JLAB_DIS/Clas_bin9.dat'
            CALL READDATA19(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, X11, X12, X13,
     *                     X14, X15, X16, X17, X18, X19, NTOT)
            DO I=1,NTOT
               CLAS9_PTLOW(I)       =X1(I)
               CLAS9_PTHIGH(I)      =X2(I)
               CLAS9_NULOW(I)       =X3(I)
               CLAS9_NUHIGH(I)      =X4(I)
               CLAS9_Q2LOW(I)       =X5(I)
               CLAS9_Q2HIGH(I)      =X6(I)
               CLAS9_RA_C(I)        =X7(I)
               CLAS9_RA_Cerr(I)     =X8(I)
               CLAS9_RA_FE(I)       =X9(I)
               CLAS9_RA_FEerr(I)    =X10(I)
               CLAS9_RA_PB(I)       =X11(I)
               CLAS9_RA_PBerr(I)    =X12(I)
               CLAS9_C_DIS(I)       =X17(I)
               CLAS9_FE_DIS(I)      =X18(I)
               CLAS9_PB_DIS(I)      =X19(I)
            END DO

C--------JLAB2022 DATA
C--------pi+
      NLINES=1
      NTOT=48
      FILEEXP = RT//'JLAB2022/pi+.csv'
            CALL READDATA20(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, X11, X12, X13,
     *                     X14,X15,X16,X17,X18,X19,X20,NTOT)
            DO I=1,NTOT
               pip2022_zlow(I)   =X1(I)
               pip2022_zup(I)    =X2(I)
               pip2022_bin(I)    =X3(I)
               pip2022_ptlow(I)  =X4(I)
               pip2022_ptup(I)   =X5(I)
               pip2022_c(I)      =X6(I)
               pip2022_cstat(I)  =X7(I)
               pip2022_csys(I)   =X8(I)
               pip2022_fe(I)     =X9(I)
               pip2022_festat(I) =X10(I)
               pip2022_fesys(I)  =X11(I)
               pip2022_pb(I)     =X12(I)
               pip2022_pbstat(I) =X13(I)
               pip2022_pbsys(I)  =X14(I)
               pip2022_qc(I)     =x15(I)
               pip2022_qfe(I)    =x16(I)
               pip2022_qpb(I)    =x17(I)
               pip2022_xc(I)     =x18(I)
               pip2022_xfe(I)    =x19(I)
               pip2022_xpb(I)    =x20(I)
            ENDDO

C--------C+ 2023
      NLINES=1
      NTOT=48
      FILEEXP = 'expdata/JLAB2022/pi+C.dat'
            CALL READDATA8(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, NTOT)
            DO I=1,NTOT
               pip2022_cfuu(I)  =X6(I)
               pip2022_cfuua(I) =X7(I)
               pip2022_cdis(I)  =X8(I)
            ENDDO

C--------Fe+ 2023
      NLINES=1
      NTOT=48
      FILEEXP = 'expdata/JLAB2022/pi+Fe.dat'
            CALL READDATA8(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, NTOT)
            DO I=1,NTOT
               pip2022_fefuu(I)  =X6(I)
               pip2022_fefuua(I) =X7(I)
               pip2022_fedis(I)  =X8(I)
            ENDDO

C--------Pb+ 2023
      NLINES=1
      NTOT=48
      FILEEXP = 'expdata/JLAB2022/pi+Pb.dat'
            CALL READDATA8(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, NTOT)
            DO I=1,NTOT
               pip2022_pbfuu(I)  =X6(I)
               pip2022_pbfuua(I) =X7(I)
               pip2022_pbdis(I)  =X8(I)
            ENDDO

C--------pi-
      NLINES=1
      NTOT=40
      FILEEXP = RT//'JLAB2022/pi-.csv'
            CALL READDATA20(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, X11, X12, X13,
     *                     X14,X15,X16,X17,X18,X19,X20,NTOT)
            DO I=1,NTOT
               pim2022_zlow(I)   =X1(I)
               pim2022_zup(I)    =X2(I)
               pim2022_bin(I)    =X3(I)
               pim2022_ptlow(I)  =X4(I)
               pim2022_ptup(I)   =X5(I)
               pim2022_c(I)      =X6(I)
               pim2022_cstat(I)  =X7(I)
               pim2022_csys(I)   =X8(I)
               pim2022_fe(I)     =X9(I)
               pim2022_festat(I) =X10(I)
               pim2022_fesys(I)  =X11(I)
               pim2022_pb(I)     =X12(I)
               pim2022_pbstat(I) =X13(I)
               pim2022_pbsys(I)  =X14(I)
               pim2022_qc(I)     =x15(I)
               pim2022_qfe(I)    =x16(I)
               pim2022_qpb(I)    =x17(I)
               pim2022_xc(I)     =x18(I)
               pim2022_xfe(I)    =x19(I)
               pim2022_xpb(I)    =x20(I)
            ENDDO

C--------C- 2023
      NLINES=1
      NTOT=40
      FILEEXP = 'expdata/JLAB2022/pi-C.dat'
            CALL READDATA8(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, NTOT)
            DO I=1,NTOT
               pim2022_cfuu(I)  =X6(I)
               pim2022_cfuua(I) =X7(I)
               pim2022_cdis(I)  =X8(I)
            ENDDO

C--------Fe- 2023
      NLINES=1
      NTOT=40
      FILEEXP = 'expdata/JLAB2022/pi-Fe.dat'
            CALL READDATA8(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, NTOT)
            DO I=1,NTOT
               pim2022_fefuu(I)  =X6(I)
               pim2022_fefuua(I) =X7(I)
               pim2022_fedis(I)  =X8(I)
            ENDDO

C--------Pb- 2023
      NLINES=1
      NTOT=40
      FILEEXP = 'expdata/JLAB2022/pi-Pb.dat'
            CALL READDATA8(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, NTOT)
            DO I=1,NTOT
               pim2022_pbfuu(I)  =X6(I)
               pim2022_pbfuua(I) =X7(I)
               pim2022_pbdis(I)  =X8(I)
            ENDDO

C-----predictions
      NLINES=1
      NTOT=30
      FILEEXP = 'expdata/PRED/EIC.dat'
            CALL READDATA3(FILEEXP,NLINES,X1,X2,X3,
     *                     NTOT)
            DO I=1,NTOT
              EIC_PRED_dis(I) = X3(I)
            ENDDO
    
      NLINES=1
      NTOT=30
      FILEEXP = 'expdata/PRED/EICC.dat'
            CALL READDATA3(FILEEXP,NLINES,X1,X2,X3,
     *                     NTOT)
            DO I=1,NTOT
              EICC_PRED_dis(I) = X3(I)
            ENDDO

      NLINES=1
      NTOT=30
      FILEEXP = 'expdata/PRED/JLAB.dat'
            CALL READDATA3(FILEEXP,NLINES,X1,X2,X3,
     *                     NTOT)
            DO I=1,NTOT
              JLAB_PRED_dis(I) = X3(I)
            ENDDO

C-----triple
      NLINES=1
      NTOT=270
      FILEEXP = 'expdata/PRED/DIS_x.dat'
            CALL READDATA5(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                     NTOT)
            DO I=1,NTOT
              x_PRED_dis(I) = X5(I)
            ENDDO
    
      NLINES=1
      NTOT=270
      FILEEXP = 'expdata/PRED/DIS_z.dat'
            CALL READDATA5(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                     NTOT)
            DO I=1,NTOT
              z_PRED_dis(I) = X5(I)
            ENDDO

      NLINES=1
      NTOT=270
      FILEEXP = 'expdata/PRED/DIS_pht.dat'
            CALL READDATA5(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                     NTOT)
            DO I=1,NTOT
              pht_PRED_dis(I) = X5(I)
            ENDDO

      return
      end
c-----------------------------------------
c     read data file which has 20 columns
c-----------------------------------------
      SUBROUTINE READDATA20(FILEEXP,NLINES,X1,
     $           X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,
     $           X16, X17, X18, X19, X20, NTOT)
      CHARACTER*72 FILEEXP,LINE
      REAL*8 X1(100),X2(100),X3(100),X4(100),X5(100),X6(100)
      REAL*8 X7(100),X8(100),X9(100),X10(100),X11(100),X12(100),X13(100)
      REAL*8 X14(100), X15(100), X16(100), X17(100), X18(100), X19(100)
      REAL*8 X20(100)
      INTEGER NLINES,NTOT,I,NT

      OPEN(unit = 91,FILE=FILEEXP,STATUS='unknown')

      DO I = 1,NLINES,1
         READ(91,'(A)') LINE    ! skip NLINES..
      ENDDO

      NT = 0
 10   NT =  NT +1

      READ(91,*,ERR = 99, end = 11)
     +     X1(NT),X2(NT),X3(NT),X4(NT),X5(NT),X6(NT),
     +     X7(NT),X8(NT),X9(NT),X10(NT),X11(NT),X12(NT),X13(NT),
     +     X14(NT), X15(NT), X16(NT), X17(NT), X18(NT), X19(NT),X20(NT)

      goto 10
 11   if(NT.GT.NTOT) goto 101

      CLOSE(91)

      GOTO 101
99    PRINT *, 'ERROR IN READING FILE    ', FILEEXP
101   CONTINUE
      RETURN
      END

c-----------------------------------------
c     read data file which has 19 columns
c-----------------------------------------
      SUBROUTINE READDATA19(FILEEXP,NLINES,X1,
     $           X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,
     $           X16, X17, X18, X19, NTOT)
      CHARACTER*72 FILEEXP,LINE
      REAL*8 X1(100),X2(100),X3(100),X4(100),X5(100),X6(100)
      REAL*8 X7(100),X8(100),X9(100),X10(100),X11(100),X12(100),X13(100)
      REAL*8 X14(100), X15(100), X16(100), X17(100), X18(100), X19(100)
      INTEGER NLINES,NTOT,I,NT

      OPEN(unit = 91,FILE=FILEEXP,STATUS='unknown')

      DO I = 1,NLINES,1
         READ(91,'(A)') LINE    ! skip NLINES..
      ENDDO

      NT = 0
 10   NT =  NT +1

      READ(91,*,ERR = 99, end = 11)
     +     X1(NT),X2(NT),X3(NT),X4(NT),X5(NT),X6(NT),
     +     X7(NT),X8(NT),X9(NT),X10(NT),X11(NT),X12(NT),X13(NT),
     +     X14(NT), X15(NT), X16(NT), X17(NT), X18(NT), X19(NT)

      goto 10
 11   if(NT.GT.NTOT) goto 101

      CLOSE(91)

      GOTO 101
99    PRINT *, 'ERROR IN READING FILE    ', FILEEXP
101   CONTINUE
      RETURN
      END

c-----------------------------------------
c     read data file which has 14 columns
c-----------------------------------------
      SUBROUTINE READDATA14(FILEEXP,NLINES,X1,
     $           X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12, X13,X14,NTOT)
      CHARACTER*72 FILEEXP,LINE
      REAL*8 X1(100),X2(100),X3(100),X4(100),X5(100),X6(100)
      REAL*8 X7(100),X8(100),X9(100),X10(100),X11(100),X12(100),X13(100)
      REAL*8 X14(100)
      INTEGER NLINES,NTOT,I,NT

      OPEN(unit = 91,FILE=FILEEXP,STATUS='unknown')

      DO I = 1,NLINES,1
       READ(91,'(A)') LINE    ! skip NLINES..
      ENDDO

      NT = 0
   10   NT =  NT +1

      READ(91,*,ERR = 99, end = 11)
     +     X1(NT),X2(NT),X3(NT),X4(NT),X5(NT),X6(NT),
     +     X7(NT),X8(NT),X9(NT),X10(NT),X11(NT),X12(NT),X13(NT),
     +     x14(NT)

      goto 10
   11 if(NT.GT.NTOT) goto 101

      CLOSE(91)

      GOTO 101
99    PRINT *, 'ERROR IN READING FILE    ', FILEEXP
101   CONTINUE
      RETURN
      END

c-----------------------------------------
c     read data file which has 13 columns
c-----------------------------------------
      SUBROUTINE READDATA13(FILEEXP,NLINES,X1,
     $           X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12, X13, NTOT)
      CHARACTER*72 FILEEXP,LINE
      REAL*8 X1(100),X2(100),X3(100),X4(100),X5(100),X6(100)
      REAL*8 X7(100),X8(100),X9(100),X10(100),X11(100),X12(100),X13(100)
      INTEGER NLINES,NTOT,I,NT

      OPEN(unit = 91,FILE=FILEEXP,STATUS='unknown')

      DO I = 1,NLINES,1
         READ(91,'(A)') LINE    ! skip NLINES..
      ENDDO

      NT = 0
 10   NT =  NT +1

      READ(91,*,ERR = 99, end = 11)
     +     X1(NT),X2(NT),X3(NT),X4(NT),X5(NT),X6(NT),
     +     X7(NT),X8(NT),X9(NT),X10(NT),X11(NT),X12(NT),X13(NT)

      goto 10
 11   if(NT.GT.NTOT) goto 101

      CLOSE(91)

      GOTO 101
99    PRINT *, 'ERROR IN READING FILE    ', FILEEXP
101   CONTINUE
      RETURN
      END

c-----------------------------------------
c     read data file which has 10 columns
c-----------------------------------------
      SUBROUTINE READDATA10(FILEEXP,NLINES,X1,
     $           X2,X3,X4,X5,X6,X7,X8,X9,X10,NTOT)
      CHARACTER*72 FILEEXP,LINE
      REAL*8 X1(100),X2(100),X3(100),X4(100),X5(100),X6(100)
      REAL*8 X7(100),X8(100),X9(100),X10(100)
      INTEGER NLINES,NTOT,I,NT

      OPEN(unit = 91,FILE=FILEEXP,STATUS='unknown')

      DO I = 1,NLINES,1
         READ(91,'(A)') LINE    ! skip NLINES..
      ENDDO

      NT = 0
 10   NT =  NT +1

      READ(91,*,ERR = 99, end = 11)
     +     X1(NT),X2(NT),X3(NT),X4(NT),X5(NT),X6(NT),
     +     X7(NT),X8(NT),X9(NT),X10(NT)

      goto 10
 11   if(NT.GT.NTOT) goto 101

      CLOSE(91)

      GOTO 101
99    PRINT *, 'ERROR IN READING FILE    ', FILEEXP
101   CONTINUE
      RETURN
      END

c-----------------------------------------
c     read data file which has 9 columns
c-----------------------------------------
      SUBROUTINE READDATA9(FILEEXP,NLINES,
     >                    X1,X2,X3,X4,X5,X6,X7,X8,X9,NTOT)
      CHARACTER*72 FILEEXP,LINE
      REAL*8 X1(100),X2(100),X3(100),X4(100),X5(100),X6(100)
      REAL*8 X7(100),X8(100),X9(100)
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
      REAL*8 X1(100),X2(100),X3(100),X4(100),X5(100),X6(100)
      REAL*8 X7(100),X8(100)
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
c     read data file which has 6 columns
c-----------------------------------------
      SUBROUTINE READDATA7(FILEEXP,NLINES,
     >                    X1,X2,X3,X4,X5,X6,X7, NTOT)
      CHARACTER*72 FILEEXP,LINE
      REAL*8 X1(100),X2(100),X3(100),X4(100),X5(100),X6(100),X7(100)
      INTEGER NLINES,NTOT,I,NT

      OPEN(unit = 91,FILE=FILEEXP,STATUS='unknown')

      DO I = 1,NLINES,1
         READ(91,'(A)') LINE    ! skip NLINES..
      ENDDO

      NT = 0
 10   NT =  NT +1

      READ(91,*,ERR = 99, end = 11)
     +     X1(NT),X2(NT),X3(NT),X4(NT),X5(NT),X6(NT),X7(NT)

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
      REAL*8 X1(100),X2(100),X3(100),X4(100),X5(100),X6(100)
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
     >                    X1,X2,X3,X4,X5, NTOT)
      CHARACTER*72 FILEEXP,LINE
      REAL*8 X1(100),X2(100),X3(100),X4(100),X5(100)
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
c     read data file which has 5 columns
c-----------------------------------------
      SUBROUTINE READDATA4(FILEEXP,NLINES,
     >                    X1,X2,X3,X4, NTOT)
      CHARACTER*72 FILEEXP,LINE
      REAL*8 X1(100),X2(100),X3(100),X4(100)
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
      REAL*8 X1(100),X2(100),X3(100)
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

c-----------------------------------------
c     read data file which has 2 columns
c-----------------------------------------
      SUBROUTINE READDATA2(FILEEXP,NLINES,
     >                    X1,X2,NTOT)
      CHARACTER*72 FILEEXP,LINE
      REAL*8 X1(100),X2(100)
      INTEGER NLINES,NTOT,I,NT

      OPEN(unit = 91,FILE=FILEEXP,STATUS='unknown')

      DO I = 1,NLINES,1
         READ(91,'(A)') LINE    ! skip NLINES..
      ENDDO

      NT = 0
 10   NT =  NT +1

      READ(91,*,ERR = 99, end = 11)
     +     X1(NT),X2(NT)

      goto 10
 11   if(NT.GT.NTOT) goto 101

      CLOSE(91)

      GOTO 101
99    PRINT *, 'ERROR IN READING FILE    ', FILEEXP
101   CONTINUE
      RETURN
      END
