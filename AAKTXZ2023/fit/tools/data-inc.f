C-----------------------------------------
C     DECLARATION OF ARRAYS USED IN READDATA.F
C-----------------------------------------
      INTEGER I_EIC_PRED
      INTEGER I_JLAB_PRED
      INTEGER I_EICC_PRED

      INTEGER NUM_EIC_z1 
      INTEGER NUM_EIC_z5 
      INTEGER NUM_EICC_z3
      INTEGER NUM_EICC_z7 
      INTEGER NUM_JLAB 
      
      INTEGER I_x_PRED
      INTEGER I_z_PRED
      INTEGER I_pht_PRED

      REAL*8      EIC_PRED_DIS, JLAB_PRED_DIS, EICC_PRED_DIS

      COMMON /PRED_DIS/ EIC_PRED_DIS, JLAB_PRED_DIS, EICC_PRED_DIS

C----- triple
      REAL*8 :: x_PRED_Qs(270),x_PRED_xb(270),
     *          x_PRED_zh(270),x_PRED_pht(270),x_PRED_dis(270) !2023
      REAL*8 :: z_PRED_Qs(270),z_PRED_xb(270),
     *          z_PRED_zh(270),z_PRED_pht(270),z_PRED_dis(270)
      REAL*8 :: pht_PRED_Qs(270),pht_PRED_xb(270),
     *          pht_PRED_zh(270),pht_PRED_pht(270),pht_PRED_dis(270)  
      COMMON /x_PRED/ x_PRED_Qs,x_PRED_xb,
     *          x_PRED_zh,x_PRED_pht,x_PRED_dis !2023
      COMMON /z_PRED/ z_PRED_Qs,z_PRED_xb,
     *          z_PRED_zh,z_PRED_pht,z_PRED_dis
      COMMON /pht_PRED/ pht_PRED_Qs,pht_PRED_xb,
     *          pht_PRED_zh,pht_PRED_pht,pht_PRED_dis  

C----- JLAB 12 GeV (preliminary data from Miguel)

      REAL*8 ::       JLAB_PIP_Z(180),
     *                JLAB_PIP_PTLOW(180),      JLAB_PIP_PTHIGH(180),
     *                JLAB_PIP_MULT_RATIO(180), JLAB_PIP_ERR(180),
     *                JLAB_PIP_FUU(180),        JLAB_PIP_FUUA(180),
     *                JLAB_PIP_DIS(180)

      INTEGER ::      JLAB_PIP_TARGET(180)


      COMMON /JLAB_PIP_12/ JLAB_PIP_TARGET,     JLAB_PIP_Z,
     *                 JLAB_PIP_PTLOW,      JLAB_PIP_PTHIGH,
     *                 JLAB_PIP_MULT_RATIO, JLAB_PIP_ERR,
     *                 JLAB_PIP_FUU,        JLAB_PIP_FUUA,
     *                 JLAB_PIP_DIS


      INTEGER I_JLAB_PIP_12


      REAL*8 ::       JLAB_PIM_Z(180),
     *                JLAB_PIM_PTLOW(180),      JLAB_PIM_PTHIGH(180),
     *                JLAB_PIM_MULT_RATIO(180), JLAB_PIM_ERR(180),
     *                JLAB_PIM_FUU(180),        JLAB_PIM_FUUA(180),
     *                JLAB_PIM_DIS(180)

      INTEGER ::      JLAB_PIM_TARGET(180)


      COMMON /JLAB_PIM_12/ JLAB_PIM_TARGET,     JLAB_PIM_Z,
     *                 JLAB_PIM_PTLOW,      JLAB_PIM_PTHIGH,
     *                 JLAB_PIM_MULT_RATIO, JLAB_PIM_ERR,
     *                 JLAB_PIM_FUU,        JLAB_PIM_FUUA,
     *                 JLAB_PIM_DIS


      INTEGER I_JLAB_PIM_12

C--------------- LHC (CMS) p + pB --> muons
      real*8 CMS8_DY_no_pt(2)
      real*8 CMS8_DY_no_CX(2)

      COMMON /CMS8_DY_no/ CMS8_DY_no_pt,CMS8_DY_no_CX

      INTEGER I_CMS8_DY_no

C--------------- LHC (CMS) p + pB --> muons
      real*8 CMS8_DY_fid_pt(2)
      real*8 CMS8_DY_fid_CX(2)

      COMMON /CMS8_DY_fid/ CMS8_DY_fid_pt,CMS8_DY_fid_CX

      INTEGER I_CMS8_DY_fid

C--------------- LHC (CMS) p + pB --> muons
      real*8 CMS8_ZZ_no_pt(9)
      real*8 CMS8_ZZ_no_CX(9)

      COMMON /CMS8_ZZ_no/ CMS8_ZZ_no_pt,CMS8_ZZ_no_CX

      INTEGER I_CMS8_ZZ_no

C--------------- LHC (CMS) p + pB --> muons
      real*8 CMS8_ZZ_fid_pt(9)
      real*8 CMS8_ZZ_fid_CX(9)

      COMMON /CMS8_ZZ_fid/ CMS8_ZZ_fid_pt,CMS8_ZZ_fid_CX

      INTEGER I_CMS8_ZZ_fid

C--------------
C-----ATLAS 5TEV

      real*8 ATLAS3_pT(12)
      real*8 ATLAS3_CX(12)
      real*8 ATLAS3_ER(12)

      integer I_ATLAS3

      common /ATLAS3/
     $  ATLAS3_pT, ATLAS3_CX, ATLAS3_ER

C--------------
C-----DRELL YAN
C-----E906 at 120 GeV
      real*8 E906_C_pT(5)
      real*8 E906_C_Ra(5)
      real*8 E906_C_err(5)

      real*8 E906_Fe_pT(5)
      real*8 E906_Fe_Ra(5)
      real*8 E906_Fe_err(5)

      real*8 E906_W_pT(5)
      real*8 E906_W_Ra(5)
      real*8 E906_W_err(5)

      integer I_E906_C
      integer I_E906_Fe
      integer I_E906_W

      common /E906/
     $   E906_C_pT, E906_C_Ra, E906_C_Err,
     $   E906_Fe_pT, E906_Fe_Ra, E906_Fe_Err,
     $   E906_W_pT, E906_W_Ra, E906_W_Err


C--------------
C-----CMS 5TEV
      real*8 CMS5_pT(13)
      real*8 CMS5_CX(13)
      real*8 CMS5_ER(13)

      integer I_CMS5

      common /CMS5/
     $  CMS5_pT, CMS5_CX, CMS5_ER

C--------------
C-----ATLAS 5TEV
      real*8 ATLAS5_Y1_pT(14)
      real*8 ATLAS5_Y1_CX(14)
      real*8 ATLAS5_Y1_ER(14)

      real*8 ATLAS5_Y2_pT(8)
      real*8 ATLAS5_Y2_CX(8)
      real*8 ATLAS5_Y2_ER(8)

      real*8 ATLAS5_Y3_pT(8)
      real*8 ATLAS5_Y3_CX(8)
      real*8 ATLAS5_Y3_ER(8)

      integer I_ATLAS5_Y1,I_ATLAS5_Y2,I_ATLAS5_Y3
      integer I_ATLAS5_Y_RAT

      common /ATLAS5/
     $  ATLAS5_Y1_pT, ATLAS5_Y1_CX, ATLAS5_Y1_ER,
     $  ATLAS5_Y2_pT, ATLAS5_Y2_CX, ATLAS5_Y2_ER,
     $  ATLAS5_Y3_pT, ATLAS5_Y3_CX, ATLAS5_Y3_ER

C--------------
C-----DRELL YAN
C-----RHIC
      real*8 RHIC_pp_pT(5)
      real*8 RHIC_pp_CX(5)
      real*8 RHIC_pp_err(5)
      real*8 RHIC_pp_FUU(5)

      real*8 RHIC_pAu1_pT(5)
      real*8 RHIC_pAu1_CX(5)
      real*8 RHIC_pAu1_err(5)

      real*8 RHIC_pAu2_pT(5)
      real*8 RHIC_pAu2_CX(5)
      real*8 RHIC_pAu2_err(5)

      integer I_RHIC_pp
      integer I_RHIC_pAu1
      integer I_RHIC_pAu2

      common /RHIC/ RHIC_pp_pT,   RHIC_pp_CX,   RHIC_pp_err,
     $              RHIC_pAu1_pT, RHIC_pAu1_CX, RHIC_pAu1_err,
     $              RHIC_pAu2_pT, RHIC_pAu2_CX, RHIC_pAu2_err,
     $              RHIC_pp_FUU

C-----RHIC (Ratios)

      real*8 RHIC_Ratio_pAu1_pT(5)
      real*8 RHIC_Ratio_pAu1_CX(5)
      real*8 RHIC_Ratio_pAu1_err(5)

      real*8 RHIC_Ratio_pAu2_pT(5)
      real*8 RHIC_Ratio_pAu2_CX(5)
      real*8 RHIC_Ratio_pAu2_err(5)

      integer I_RHIC_Ratio_pAu1
      integer I_RHIC_Ratio_pAu2

      common /RHIC_Ratio/ RHIC_Ratio_pAu1_pT, RHIC_Ratio_pAu1_CX,
     $  RHIC_Ratio_pAu1_err, RHIC_Ratio_pAu2_pT,
     $  RHIC_Ratio_pAu2_CX, RHIC_Ratio_pAu2_err

C-------------
c-----E288
      real*8 E288_200_pT(119)
      real*8 E288_200_Q(119)
      real*8 E288_200_ds(119)
      real*8 E288_200_sigma(119)
      integer I_E288_200
      common /E288_200/ E288_200_pT,
     $                  E288_200_Q,
     $                  E288_200_ds,
     $ 	       	        E288_200_sigma

      real*8 E288_300_pT(184)
      real*8 E288_300_Q(184)
      real*8 E288_300_ds(184)
      real*8 E288_300_sigma(184)
      integer I_E288_300
      common /E288_300/ E288_300_pT,
     $                  E288_300_Q,
     $                  E288_300_ds,
     $ 	       	        E288_300_sigma

      real*8 E288_400_pT(225)
      real*8 E288_400_Q(225)
      real*8 E288_400_ds(225)
      real*8 E288_400_sigma(225)
      integer I_E288_400
      common /E288_400/ E288_400_pT,
     $                  E288_400_Q,
     $                  E288_400_ds,
     $ 	       	        E288_400_sigma


C-----E866
      real*8 E866_800_pT(40)
      real*8 E866_800_R_FeBe(40)
      real*8 E866_800_Err_FeBe(40)
      real*8 E866_800_R_WBe(40)
      real*8 E866_800_Err_WBe(40)
      real*8 E866_800_x2(40)
      integer I_E866_800

      common /E866_800/ E866_800_pT,
     $                  E866_800_R_FeBe,
     $                  E866_800_Err_FeBe,
     $ 	       	        E866_800_R_WBe,
     $                  E866_800_Err_WBe,
     $                  E866_800_x2

C-----E866
      real*8 E866_800q_pT(32)
      real*8 E866_800q_R_FeBe(32)
      real*8 E866_800q_Err_FeBe(32)
      real*8 E866_800q_R_WBe(32)
      real*8 E866_800q_Err_WBe(32)
      real*8 E866_800q_Q(32)

      integer I_E866_800q
      common /E866_800q/ E866_800q_pT,
     $                   E866_800q_R_FeBe,
     $                   E866_800q_Err_FeBe,
     $ 	       	         E866_800q_R_WBe,
     $                   E866_800q_Err_WBe,
     $                   E866_800q_Q

C-----E772
      real*8 E772_800_C_pT(7)
      real*8 E772_800_C_R(7)
      real*8 E772_800_C_Err(7)
      real*8 E772_800_Ca_pT(7)
      real*8 E772_800_Ca_R(7)
      real*8 E772_800_Ca_Err(7)
      real*8 E772_800_Fe_pT(7)
      real*8 E772_800_Fe_R(7)
      real*8 E772_800_Fe_Err(7)
      real*8 E772_800_W_pT(7)
      real*8 E772_800_W_R(7)
      real*8 E772_800_W_Err(7)

      integer I_E772_800

      common /E772_800/ E772_800_C_pT,
     $                  E772_800_C_R,
     $                  E772_800_C_Err,
     $ 	       	        E772_800_Ca_pT,
     $                  E772_800_Ca_R,
     $                  E772_800_Ca_Err,
     $ 	       	        E772_800_Fe_pT,
     $                  E772_800_Fe_R,
     $                  E772_800_Fe_Err,
     $ 	       	        E772_800_W_pT,
     $                  E772_800_W_R,
     $                  E772_800_W_Err

C-----E605: RTS : 38.8 GeV
      real*8 E605_pT(74)
      real*8 E605_ds(74)
      real*8 E605_err(74)
      real*8 E605_xf(74)
      real*8 E605_Qlo(74)
      real*8 E605_Qhi(74)

      integer I_E605_800

      common /E605_38pt8/ E605_pT, E605_ds, E605_err,
     $                    E605_xf, E605_Qlo, E605_Qhi

C-----SIDIS (HERMES DATA)
C-----------------------------------------
C------PI+ NU DEPENDENT DATA------------------------
C------- HE
      REAL*8 ::          PIP_HE_NU_VALUE(9), PIP_HE_NU_Z(9),
     *                   PIP_HE_NU_Q2(9),    PIP_HE_NU_PT2(9),
     *                   PIP_HE_NU_MULT(9),  PIP_HE_NU_FUU(9),
     *                   PIP_HE_NU_STAT(9),  PIP_HE_NU_SYS(9),
     *                   PIP_HE_NU_FUUA(9),  PIP_HE_NU_DIS(9)

      COMMON /PIP_HE_NU/ PIP_HE_NU_VALUE,    PIP_HE_NU_Z,
     *                   PIP_HE_NU_Q2,       PIP_HE_NU_PT2,
     *                   PIP_HE_NU_MULT,     PIP_HE_NU_FUU,
     *                   PIP_HE_NU_STAT,     PIP_HE_NU_SYS,
     *                   PIP_HE_NU_FUUA,     PIP_HE_NU_DIS


      INTEGER I_PIP_HE_NU
C------- NE
      REAL*8 ::          PIP_NE_NU_VALUE(9), PIP_NE_NU_Z(9),
     *                   PIP_NE_NU_Q2(9),    PIP_NE_NU_PT2(9),
     *                   PIP_NE_NU_MULT(9),  PIP_NE_NU_FUU(9),
     *                   PIP_NE_NU_STAT(9),  PIP_NE_NU_SYS(9),
     *                   PIP_NE_NU_FUUA(9),  PIP_NE_NU_DIS(9)

      COMMON /PIP_NE_NU/ PIP_NE_NU_VALUE,    PIP_NE_NU_Z,
     *                   PIP_NE_NU_Q2,       PIP_NE_NU_PT2,
     *                   PIP_NE_NU_MULT,     PIP_NE_NU_FUU,
     *                   PIP_NE_NU_STAT,     PIP_NE_NU_SYS,
     *                   PIP_NE_NU_FUUA,     PIP_NE_NU_DIS


      INTEGER I_PIP_NE_NU
C------- KR
      REAL*8 ::          PIP_KR_NU_VALUE(9), PIP_KR_NU_Z(9),
     *                   PIP_KR_NU_Q2(9),    PIP_KR_NU_PT2(9),
     *                   PIP_KR_NU_MULT(9),  PIP_KR_NU_FUU(9),
     *                   PIP_KR_NU_STAT(9),  PIP_KR_NU_SYS(9),
     *                   PIP_KR_NU_FUUA(9),  PIP_KR_NU_DIS(9)

      COMMON /PIP_KR_NU/ PIP_KR_NU_VALUE,    PIP_KR_NU_Z,
     *                   PIP_KR_NU_Q2,       PIP_KR_NU_PT2,
     *                   PIP_KR_NU_MULT,     PIP_KR_NU_FUU,
     *                   PIP_KR_NU_STAT,     PIP_KR_NU_SYS,
     *                   PIP_KR_NU_FUUA,     PIP_KR_NU_DIS

      INTEGER I_PIP_KR_NU
C------- XE
      REAL*8 ::          PIP_XE_NU_VALUE(9), PIP_XE_NU_Z(9),
     *                   PIP_XE_NU_Q2(9),    PIP_XE_NU_PT2(9),
     *                   PIP_XE_NU_MULT(9),  PIP_XE_NU_FUU(9),
     *                   PIP_XE_NU_STAT(9),  PIP_XE_NU_SYS(9),
     *                   PIP_XE_NU_FUUA(9),  PIP_XE_NU_DIS(9)

      COMMON /PIP_XE_NU/ PIP_XE_NU_VALUE,    PIP_XE_NU_Z,
     *                   PIP_XE_NU_Q2,       PIP_XE_NU_PT2,
     *                   PIP_XE_NU_MULT,     PIP_XE_NU_FUU,
     *                   PIP_XE_NU_STAT,     PIP_XE_NU_SYS,
     *                   PIP_XE_NU_FUUA,     PIP_XE_NU_DIS

      INTEGER I_PIP_XE_NU

C------PI+ Z DEPENDENT DATA------------------------
C------- HE
      REAL*8 ::          PIP_HE_Z_VALUE(9), PIP_HE_Z_NU(9),
     *                   PIP_HE_Z_Q2(9),    PIP_HE_Z_PT2(9),
     *                   PIP_HE_Z_MULT(9),  PIP_HE_Z_FUU(9),
     *                   PIP_HE_Z_STAT(9),  PIP_HE_Z_SYS(9),
     *                   PIP_HE_Z_FUUA(9),  PIP_HE_Z_DIS(9)

      COMMON /PIP_HE_Z/  PIP_HE_Z_VALUE,    PIP_HE_Z_NU,
     *                   PIP_HE_Z_Q2,       PIP_HE_Z_PT2,
     *                   PIP_HE_Z_MULT,     PIP_HE_Z_FUU,
     *                   PIP_HE_Z_STAT,     PIP_HE_Z_SYS,
     *                   PIP_HE_Z_FUUA,     PIP_HE_Z_DIS


      INTEGER I_PIP_HE_Z
C------- NE
      REAL*8 ::          PIP_NE_Z_VALUE(9), PIP_NE_Z_NU(9),
     *                   PIP_NE_Z_Q2(9),    PIP_NE_Z_PT2(9),
     *                   PIP_NE_Z_MULT(9),  PIP_NE_Z_FUU(9),
     *                   PIP_NE_Z_STAT(9),  PIP_NE_Z_SYS(9),
     *                   PIP_NE_Z_FUUA(9),  PIP_NE_Z_DIS(9)

      COMMON /PIP_NE_Z/  PIP_NE_Z_VALUE,    PIP_NE_Z_NU,
     *                   PIP_NE_Z_Q2,       PIP_NE_Z_PT2,
     *                   PIP_NE_Z_MULT,     PIP_NE_Z_FUU,
     *                   PIP_NE_Z_STAT,     PIP_NE_Z_SYS,
     *                   PIP_NE_Z_FUUA,     PIP_NE_Z_DIS

      INTEGER I_PIP_NE_Z
C------- KR
      REAL*8 ::          PIP_KR_Z_VALUE(9), PIP_KR_Z_NU(9),
     *                   PIP_KR_Z_Q2(9),    PIP_KR_Z_PT2(9),
     *                   PIP_KR_Z_MULT(9),  PIP_KR_Z_FUU(9),
     *                   PIP_KR_Z_STAT(9),  PIP_KR_Z_SYS(9),
     *                   PIP_KR_Z_FUUA(9),  PIP_KR_Z_DIS(9)

      COMMON /PIP_KR_Z/  PIP_KR_Z_VALUE,    PIP_KR_Z_NU,
     *                   PIP_KR_Z_Q2,       PIP_KR_Z_PT2,
     *                   PIP_KR_Z_MULT,     PIP_KR_Z_FUU,
     *                   PIP_KR_Z_STAT,     PIP_KR_Z_SYS,
     *                   PIP_KR_Z_FUUA,     PIP_KR_Z_DIS

      INTEGER I_PIP_KR_Z
C------- XE
      REAL*8 ::          PIP_XE_Z_VALUE(9), PIP_XE_Z_NU(9),
     *                   PIP_XE_Z_Q2(9),    PIP_XE_Z_PT2(9),
     *                   PIP_XE_Z_MULT(9),  PIP_XE_Z_FUU(9),
     *                   PIP_XE_Z_STAT(9),  PIP_XE_Z_SYS(9),
     *                   PIP_XE_Z_FUUA(9),  PIP_XE_Z_DIS(9)

      COMMON /PIP_XE_NU/ PIP_XE_Z_VALUE,    PIP_XE_Z_NU,
     *                   PIP_XE_Z_Q2,       PIP_XE_Z_PT2,
     *                   PIP_XE_Z_MULT,     PIP_XE_Z_FUU,
     *                   PIP_XE_Z_STAT,     PIP_XE_Z_SYS,
     *                   PIP_XE_Z_FUUA,     PIP_XE_Z_DIS

      INTEGER I_PIP_XE_Z

C------PI+ Q2 DEPENDENT DATA------------------------
C------- HE
      REAL*8 ::          PIP_HE_Q2_VALUE(9), PIP_HE_Q2_NU(9),
     *                   PIP_HE_Q2_Z(9),     PIP_HE_Q2_PT2(9),
     *                   PIP_HE_Q2_MULT(9),  PIP_HE_Q2_FUU(9),
     *                   PIP_HE_Q2_STAT(9),  PIP_HE_Q2_SYS(9),
     *                   PIP_HE_Q2_FUUA(9),  PIP_HE_Q2_DIS(9)

      COMMON /PIP_HE_Q2/ PIP_HE_Q2_VALUE,    PIP_HE_Q2_NU,
     *                   PIP_HE_Q2_Z,        PIP_HE_Q2_PT2,
     *                   PIP_HE_Q2_MULT,     PIP_HE_Q2_FUU,
     *                   PIP_HE_Q2_STAT,     PIP_HE_Q2_SYS,
     *                   PIP_HE_Q2_FUUA,     PIP_HE_Q2_DIS

      INTEGER I_PIP_HE_Q2
C------- NE
      REAL*8 ::          PIP_NE_Q2_VALUE(9), PIP_NE_Q2_NU(9),
     *                   PIP_NE_Q2_Z(9),     PIP_NE_Q2_PT2(9),
     *                   PIP_NE_Q2_MULT(9),  PIP_NE_Q2_FUU(9),
     *                   PIP_NE_Q2_STAT(9),  PIP_NE_Q2_SYS(9),
     *                   PIP_NE_Q2_FUUA(9),  PIP_NE_Q2_DIS(9)

      COMMON /PIP_NE_Q2/ PIP_NE_Q2_VALUE,    PIP_NE_Q2_NU,
     *                   PIP_NE_Q2_Z,        PIP_NE_Q2_PT2,
     *                   PIP_NE_Q2_MULT,     PIP_NE_Q2_FUU,
     *                   PIP_NE_Q2_STAT,     PIP_NE_Q2_SYS,
     *                   PIP_NE_Q2_FUUA,     PIP_NE_Q2_DIS


      INTEGER I_PIP_NE_Q2
C------- KR
      REAL*8 ::          PIP_KR_Q2_VALUE(9), PIP_KR_Q2_NU(9),
     *                   PIP_KR_Q2_Z(9),     PIP_KR_Q2_PT2(9),
     *                   PIP_KR_Q2_MULT(9),  PIP_KR_Q2_FUU(9),
     *                   PIP_KR_Q2_STAT(9),  PIP_KR_Q2_SYS(9),
     *                   PIP_KR_Q2_FUUA(9),  PIP_KR_Q2_DIS(9)

      COMMON /PIP_KR_Q2/ PIP_KR_Q2_VALUE,    PIP_KR_Q2_NU,
     *                   PIP_KR_Q2_Z,        PIP_KR_Q2_PT2,
     *                   PIP_KR_Q2_MULT,     PIP_KR_Q2_FUU,
     *                   PIP_KR_Q2_STAT,     PIP_KR_Q2_SYS,
     *                   PIP_KR_Q2_FUUA,     PIP_KR_Q2_DIS

      INTEGER I_PIP_KR_Q2
C------- XE
      REAL*8 ::          PIP_XE_Q2_VALUE(9), PIP_XE_Q2_NU(9),
     *                   PIP_XE_Q2_Z(9),     PIP_XE_Q2_PT2(9),
     *                   PIP_XE_Q2_MULT(9),  PIP_XE_Q2_FUU(9),
     *                   PIP_XE_Q2_STAT(9),  PIP_XE_Q2_SYS(9),
     *                   PIP_XE_Q2_FUUA(9),  PIP_XE_Q2_DIS(9)

      COMMON /PIP_XE_Q2/ PIP_XE_Q2_VALUE,    PIP_XE_Q2_NU,
     *                   PIP_XE_Q2_Z,        PIP_XE_Q2_PT2,
     *                   PIP_XE_Q2_MULT,     PIP_XE_Q2_FUU,
     *                   PIP_XE_Q2_STAT,     PIP_XE_Q2_SYS,
     *                   PIP_XE_Q2_FUUA,     PIP_XE_Q2_DIS

      INTEGER I_PIP_XE_Q2
C-----------------------------------------
C------K+ NU DEPENDENT DATA------------------------
C------- HE
      REAL*8 ::          KP_HE_NU_VALUE(9), KP_HE_NU_Z(9),
     *                   KP_HE_NU_Q2(9),    KP_HE_NU_PT2(9),
     *                   KP_HE_NU_MULT(9),  KP_HE_NU_FUU(9),
     *                   KP_HE_NU_STAT(9),  KP_HE_NU_SYS(9),
     *                   KP_HE_NU_FUUA(9),  KP_HE_NU_DIS(9)

      COMMON /KP_HE_NU/  KP_HE_NU_VALUE,    KP_HE_NU_Z,
     *                   KP_HE_NU_Q2,       KP_HE_NU_PT2,
     *                   KP_HE_NU_MULT,     KP_HE_NU_FUU,
     *                   KP_HE_NU_STAT,     KP_HE_NU_SYS,
     *                   KP_HE_NU_FUUA,     KP_HE_NU_DIS


      INTEGER I_KP_HE_NU
C------- NE
      REAL*8 ::          KP_NE_NU_VALUE(9), KP_NE_NU_Z(9),
     *                   KP_NE_NU_Q2(9),    KP_NE_NU_PT2(9),
     *                   KP_NE_NU_MULT(9),  KP_NE_NU_FUU(9),
     *                   KP_NE_NU_STAT(9),  KP_NE_NU_SYS(9),
     *                   KP_NE_NU_FUUA(9),  KP_NE_NU_DIS(9)

      COMMON /KP_NE_NU/  KP_NE_NU_VALUE,    KP_NE_NU_Z,
     *                   KP_NE_NU_Q2,       KP_NE_NU_PT2,
     *                   KP_NE_NU_MULT,     KP_NE_NU_FUU,
     *                   KP_NE_NU_STAT,     KP_NE_NU_SYS,
     *                   KP_NE_NU_FUUA,     KP_NE_NU_DIS

      INTEGER I_KP_NE_NU
C------- KR
      REAL*8 ::          KP_KR_NU_VALUE(9), KP_KR_NU_Z(9),
     *                   KP_KR_NU_Q2(9),    KP_KR_NU_PT2(9),
     *                   KP_KR_NU_MULT(9),  KP_KR_NU_FUU(9),
     *                   KP_KR_NU_STAT(9),  KP_KR_NU_SYS(9),
     *                   KP_KR_NU_FUUA(9),  KP_KR_NU_DIS(9)

      COMMON /KP_KR_NU/  KP_KR_NU_VALUE,    KP_KR_NU_Z,
     *                   KP_KR_NU_Q2,       KP_KR_NU_PT2,
     *                   KP_KR_NU_MULT,     KP_KR_NU_FUU,
     *                   KP_KR_NU_STAT,     KP_KR_NU_SYS,
     *                   KP_KR_NU_FUUA,     KP_KR_NU_DIS

      INTEGER I_KP_KR_NU
C------- XE
      REAL*8 ::          KP_XE_NU_VALUE(9), KP_XE_NU_Z(9),
     *                   KP_XE_NU_Q2(9),    KP_XE_NU_PT2(9),
     *                   KP_XE_NU_MULT(9),  KP_XE_NU_FUU(9),
     *                   KP_XE_NU_STAT(9),  KP_XE_NU_SYS(9),
     *                   KP_XE_NU_FUUA(9),  KP_XE_NU_DIS(9)

      COMMON /KP_XE_NU/  KP_XE_NU_VALUE,    KP_XE_NU_Z,
     *                   KP_XE_NU_Q2,       KP_XE_NU_PT2,
     *                   KP_XE_NU_MULT,     KP_XE_NU_FUU,
     *                   KP_XE_NU_STAT,     KP_XE_NU_SYS,
     *                   KP_XE_NU_FUUA,     KP_XE_NU_DIS

      INTEGER I_KP_XE_NU

C------K+ Z DEPENDENT DATA------------------------
C------- HE
      REAL*8 ::          KP_HE_Z_VALUE(9), KP_HE_Z_NU(9),
     *                   KP_HE_Z_Q2(9),    KP_HE_Z_PT2(9),
     *                   KP_HE_Z_MULT(9),  KP_HE_Z_FUU(9),
     *                   KP_HE_Z_STAT(9),  KP_HE_Z_SYS(9),
     *                   KP_HE_Z_FUUA(9),  KP_HE_Z_DIS(9)

      COMMON /KP_HE_Z/   KP_HE_Z_VALUE,    KP_HE_Z_NU,
     *                   KP_HE_Z_Q2,       KP_HE_Z_PT2,
     *                   KP_HE_Z_MULT,     KP_HE_Z_FUU,
     *                   KP_HE_Z_STAT,     KP_HE_Z_SYS,
     *                   KP_HE_Z_FUUA,     KP_HE_Z_DIS


      INTEGER I_KP_HE_Z
C------- NE
      REAL*8 ::          KP_NE_Z_VALUE(9), KP_NE_Z_NU(9),
     *                   KP_NE_Z_Q2(9),    KP_NE_Z_PT2(9),
     *                   KP_NE_Z_MULT(9),  KP_NE_Z_FUU(9),
     *                   KP_NE_Z_STAT(9),  KP_NE_Z_SYS(9),
     *                   KP_NE_Z_FUUA(9),  KP_NE_Z_DIS(9)

      COMMON /KP_NE_Z/   KP_NE_Z_VALUE,    KP_NE_Z_NU,
     *                   KP_NE_Z_Q2,       KP_NE_Z_PT2,
     *                   KP_NE_Z_MULT,     KP_NE_Z_FUU,
     *                   KP_NE_Z_STAT,     KP_NE_Z_SYS,
     *                   KP_NE_Z_FUUA,     KP_NE_Z_DIS

      INTEGER I_KP_NE_Z
C------- K
      REAL*8 ::          KP_KR_Z_VALUE(9), KP_KR_Z_NU(9),
     *                   KP_KR_Z_Q2(9),    KP_KR_Z_PT2(9),
     *                   KP_KR_Z_MULT(9),  KP_KR_Z_FUU(9),
     *                   KP_KR_Z_STAT(9),  KP_KR_Z_SYS(9),
     *                   KP_KR_Z_FUUA(9),  KP_KR_Z_DIS(9)

      COMMON /KP_KR_Z/   KP_KR_Z_VALUE,    KP_KR_Z_NU,
     *                   KP_KR_Z_Q2,       KP_KR_Z_PT2,
     *                   KP_KR_Z_MULT,     KP_KR_Z_FUU,
     *                   KP_KR_Z_STAT,     KP_KR_Z_SYS,
     *                   KP_KR_Z_FUUA,     KP_KR_Z_DIS

      INTEGER I_KP_KR_Z
C------- XE
      REAL*8 ::          KP_XE_Z_VALUE(9), KP_XE_Z_NU(9),
     *                   KP_XE_Z_Q2(9),    KP_XE_Z_PT2(9),
     *                   KP_XE_Z_MULT(9),  KP_XE_Z_FUU(9),
     *                   KP_XE_Z_STAT(9),  KP_XE_Z_SYS(9),
     *                   KP_XE_Z_FUUA(9),  KP_XE_Z_DIS(9)

      COMMON /KP_XE_NU/  KP_XE_Z_VALUE,    KP_XE_Z_NU,
     *                   KP_XE_Z_Q2,       KP_XE_Z_PT2,
     *                   KP_XE_Z_MULT,     KP_XE_Z_FUU,
     *                   KP_XE_Z_STAT,     KP_XE_Z_SYS,
     *                   KP_XE_Z_FUUA,     KP_XE_Z_DIS

      INTEGER I_KP_XE_Z

C------K+ Q2 DEPENDENT DATA------------------------
C------- HE
      REAL*8 ::          KP_HE_Q2_VALUE(9), KP_HE_Q2_NU(9),
     *                   KP_HE_Q2_Z(9),     KP_HE_Q2_PT2(9),
     *                   KP_HE_Q2_MULT(9),  KP_HE_Q2_FUU(9),
     *                   KP_HE_Q2_STAT(9),  KP_HE_Q2_SYS(9),
     *                   KP_HE_Q2_FUUA(9),  KP_HE_Q2_DIS(9)

      COMMON /KP_HE_Q2/  KP_HE_Q2_VALUE,    KP_HE_Q2_NU,
     *                   KP_HE_Q2_Z,        KP_HE_Q2_PT2,
     *                   KP_HE_Q2_MULT,     KP_HE_Q2_FUU,
     *                   KP_HE_Q2_STAT,     KP_HE_Q2_SYS,
     *                   KP_HE_Q2_FUUA,     KP_HE_Q2_DIS

      INTEGER I_KP_HE_Q2
C------- NE
      REAL*8 ::          KP_NE_Q2_VALUE(9), KP_NE_Q2_NU(9),
     *                   KP_NE_Q2_Z(9),     KP_NE_Q2_PT2(9),
     *                   KP_NE_Q2_MULT(9),  KP_NE_Q2_FUU(9),
     *                   KP_NE_Q2_STAT(9),  KP_NE_Q2_SYS(9),
     *                   KP_NE_Q2_FUUA(9),  KP_NE_Q2_DIS(9)

      COMMON /KP_NE_Q2/  KP_NE_Q2_VALUE,    KP_NE_Q2_NU,
     *                   KP_NE_Q2_Z,        KP_NE_Q2_PT2,
     *                   KP_NE_Q2_MULT,     KP_NE_Q2_FUU,
     *                   KP_NE_Q2_STAT,     KP_NE_Q2_SYS,
     *                   KP_NE_Q2_FUUA,     KP_NE_Q2_DIS


      INTEGER I_KP_NE_Q2
C------- KR
      REAL*8 ::          KP_KR_Q2_VALUE(9), KP_KR_Q2_NU(9),
     *                   KP_KR_Q2_Z(9),     KP_KR_Q2_PT2(9),
     *                   KP_KR_Q2_MULT(9),  KP_KR_Q2_FUU(9),
     *                   KP_KR_Q2_STAT(9),  KP_KR_Q2_SYS(9),
     *                   KP_KR_Q2_FUUA(9),  KP_KR_Q2_DIS(9)

      COMMON /KP_KR_Q2/  KP_KR_Q2_VALUE,    KP_KR_Q2_NU,
     *                   KP_KR_Q2_Z,        KP_KR_Q2_PT2,
     *                   KP_KR_Q2_MULT,     KP_KR_Q2_FUU,
     *                   KP_KR_Q2_STAT,     KP_KR_Q2_SYS,
     *                   KP_KR_Q2_FUUA,     KP_KR_Q2_DIS
      INTEGER I_KP_KR_Q2
C------- XE
      REAL*8 ::          KP_XE_Q2_VALUE(9), KP_XE_Q2_NU(9),
     *                   KP_XE_Q2_Z(9),     KP_XE_Q2_PT2(9),
     *                   KP_XE_Q2_MULT(9),  KP_XE_Q2_FUU(9),
     *                   KP_XE_Q2_STAT(9),  KP_XE_Q2_SYS(9),
     *                   KP_XE_Q2_FUUA(9),  KP_XE_Q2_DIS(9)

      COMMON /KP_XE_Q2/  KP_XE_Q2_VALUE,    KP_XE_Q2_NU,
     *                   KP_XE_Q2_Z,        KP_XE_Q2_PT2,
     *                   KP_XE_Q2_MULT,     KP_XE_Q2_FUU,
     *                   KP_XE_Q2_STAT,     KP_XE_Q2_SYS,
     *                   KP_XE_Q2_FUUA,     KP_XE_Q2_DIS
      INTEGER I_KP_XE_Q2

C-----------------------------------------
C------PI- NU DEPENDENT DATA------------------------
C------- HE
      REAL*8 ::          PIM_HE_NU_VALUE(9), PIM_HE_NU_Z(9),
     *                   PIM_HE_NU_Q2(9),    PIM_HE_NU_PT2(9),
     *                   PIM_HE_NU_MULT(9),  PIM_HE_NU_FUU(9),
     *                   PIM_HE_NU_STAT(9),  PIM_HE_NU_SYS(9),
     *                   PIM_HE_NU_FUUA(9),  PIM_HE_NU_DIS(9)

      COMMON /PIM_HE_NU/ PIM_HE_NU_VALUE,    PIM_HE_NU_Z,
     *                   PIM_HE_NU_Q2,       PIM_HE_NU_PT2,
     *                   PIM_HE_NU_MULT,     PIM_HE_NU_FUU,
     *                   PIM_HE_NU_STAT,     PIM_HE_NU_SYS,
     *                   PIM_HE_NU_FUUA,     PIM_HE_NU_DIS

      INTEGER I_PIM_HE_NU
C------- NE
      REAL*8 ::          PIM_NE_NU_VALUE(9), PIM_NE_NU_Z(9),
     *                   PIM_NE_NU_Q2(9),    PIM_NE_NU_PT2(9),
     *                   PIM_NE_NU_MULT(9),  PIM_NE_NU_FUU(9),
     *                   PIM_NE_NU_STAT(9),  PIM_NE_NU_SYS(9),
     *                   PIM_NE_NU_FUUA(9),  PIM_NE_NU_DIS(9)

      COMMON /PIM_NE_NU/ PIM_NE_NU_VALUE,    PIM_NE_NU_Z,
     *                   PIM_NE_NU_Q2,       PIM_NE_NU_PT2,
     *                   PIM_NE_NU_MULT,     PIM_NE_NU_FUU,
     *                   PIM_NE_NU_STAT,     PIM_NE_NU_SYS,
     *                   PIM_NE_NU_FUUA,     PIM_NE_NU_DIS

      INTEGER I_PIM_NE_NU
C------- KR
      REAL*8 ::          PIM_KR_NU_VALUE(9), PIM_KR_NU_Z(9),
     *                   PIM_KR_NU_Q2(9),    PIM_KR_NU_PT2(9),
     *                   PIM_KR_NU_MULT(9),  PIM_KR_NU_FUU(9),
     *                   PIM_KR_NU_STAT(9),  PIM_KR_NU_SYS(9),
     *                   PIM_KR_NU_FUUA(9),  PIM_KR_NU_DIS(9)

      COMMON /PIM_KR_NU/ PIM_KR_NU_VALUE,    PIM_KR_NU_Z,
     *                   PIM_KR_NU_Q2,       PIM_KR_NU_PT2,
     *                   PIM_KR_NU_MULT,     PIM_KR_NU_FUU,
     *                   PIM_KR_NU_STAT,     PIM_KR_NU_SYS,
     *                   PIM_KR_NU_FUUA,     PIM_KR_NU_DIS

      INTEGER I_PIM_KR_NU
C------- XE
      REAL*8 ::          PIM_XE_NU_VALUE(9), PIM_XE_NU_Z(9),
     *                   PIM_XE_NU_Q2(9),    PIM_XE_NU_PT2(9),
     *                   PIM_XE_NU_MULT(9),  PIM_XE_NU_FUU(9),
     *                   PIM_XE_NU_STAT(9),  PIM_XE_NU_SYS(9),
     *                   PIM_XE_NU_FUUA(9),  PIM_XE_NU_DIS(9)

      COMMON /PIM_XE_NU/ PIM_XE_NU_VALUE,    PIM_XE_NU_Z,
     *                   PIM_XE_NU_Q2,       PIM_XE_NU_PT2,
     *                   PIM_XE_NU_MULT,     PIM_XE_NU_FUU,
     *                   PIM_XE_NU_STAT,     PIM_XE_NU_SYS,
     *                   PIM_XE_NU_FUUA,     PIM_XE_NU_DIS

      INTEGER I_PIM_XE_NU

C------PI- Z DEPENDENT DATA------------------------
C------- HE
      REAL*8 ::          PIM_HE_Z_VALUE(9), PIM_HE_Z_NU(9),
     *                   PIM_HE_Z_Q2(9),    PIM_HE_Z_PT2(9),
     *                   PIM_HE_Z_MULT(9),  PIM_HE_Z_FUU(9),
     *                   PIM_HE_Z_STAT(9),  PIM_HE_Z_SYS(9),
     *                   PIM_HE_Z_FUUA(9),  PIM_HE_Z_DIS(9)

      COMMON /PIM_HE_Z/  PIM_HE_Z_VALUE,    PIM_HE_Z_NU,
     *                   PIM_HE_Z_Q2,       PIM_HE_Z_PT2,
     *                   PIM_HE_Z_MULT,     PIM_HE_Z_FUU,
     *                   PIM_HE_Z_STAT,     PIM_HE_Z_SYS,
     *                   PIM_HE_Z_FUUA,     PIM_HE_Z_DIS

      INTEGER I_PIM_HE_Z
C------- NE
      REAL*8 ::          PIM_NE_Z_VALUE(9), PIM_NE_Z_NU(9),
     *                   PIM_NE_Z_Q2(9),    PIM_NE_Z_PT2(9),
     *                   PIM_NE_Z_MULT(9),  PIM_NE_Z_FUU(9),
     *                   PIM_NE_Z_STAT(9),  PIM_NE_Z_SYS(9),
     *                   PIM_NE_Z_FUUA(9),  PIM_NE_Z_DIS(9)

      COMMON /PIM_NE_Z/  PIM_NE_Z_VALUE,    PIM_NE_Z_NU,
     *                   PIM_NE_Z_Q2,       PIM_NE_Z_PT2,
     *                   PIM_NE_Z_MULT,     PIM_NE_Z_FUU,
     *                   PIM_NE_Z_STAT,     PIM_NE_Z_SYS,
     *                   PIM_NE_Z_FUUA,     PIM_NE_Z_DIS

      INTEGER I_PIM_NE_Z
C------- KR
      REAL*8 ::          PIM_KR_Z_VALUE(9), PIM_KR_Z_NU(9),
     *                   PIM_KR_Z_Q2(9),    PIM_KR_Z_PT2(9),
     *                   PIM_KR_Z_MULT(9),  PIM_KR_Z_FUU(9),
     *                   PIM_KR_Z_STAT(9),  PIM_KR_Z_SYS(9),
     *                   PIM_KR_Z_FUUA(9),  PIM_KR_Z_DIS(9)

      COMMON /PIM_KR_Z/  PIM_KR_Z_VALUE,    PIM_KR_Z_NU,
     *                   PIM_KR_Z_Q2,       PIM_KR_Z_PT2,
     *                   PIM_KR_Z_MULT,     PIM_KR_Z_FUU,
     *                   PIM_KR_Z_STAT,     PIM_KR_Z_SYS,
     *                   PIM_KR_Z_FUUA,     PIM_KR_Z_DIS

      INTEGER I_PIM_KR_Z
C------- XE
      REAL*8 ::          PIM_XE_Z_VALUE(9), PIM_XE_Z_NU(9),
     *                   PIM_XE_Z_Q2(9),    PIM_XE_Z_PT2(9),
     *                   PIM_XE_Z_MULT(9),  PIM_XE_Z_FUU(9),
     *                   PIM_XE_Z_STAT(9),  PIM_XE_Z_SYS(9),
     *                   PIM_XE_Z_FUUA(9),  PIM_XE_Z_DIS(9)

      COMMON /PIM_XE_NU/ PIM_XE_Z_VALUE,    PIM_XE_Z_NU,
     *                   PIM_XE_Z_Q2,       PIM_XE_Z_PT2,
     *                   PIM_XE_Z_MULT,     PIM_XE_Z_FUU,
     *                   PIM_XE_Z_STAT,     PIM_XE_Z_SYS,
     *                   PIM_XE_Z_FUUA,     PIM_XE_Z_DIS

      INTEGER I_PIM_XE_Z

C------PI- Q2 DEPENDENT DATA------------------------
C------- HE
      REAL*8 ::          PIM_HE_Q2_VALUE(9), PIM_HE_Q2_NU(9),
     *                   PIM_HE_Q2_Z(9),     PIM_HE_Q2_PT2(9),
     *                   PIM_HE_Q2_MULT(9),  PIM_HE_Q2_FUU(9),
     *                   PIM_HE_Q2_STAT(9),  PIM_HE_Q2_SYS(9),
     *                   PIM_HE_Q2_FUUA(9),  PIM_HE_Q2_DIS(9)

      COMMON /PIM_HE_Q2/ PIM_HE_Q2_VALUE,    PIM_HE_Q2_NU,
     *                   PIM_HE_Q2_Z,        PIM_HE_Q2_PT2,
     *                   PIM_HE_Q2_MULT,     PIM_HE_Q2_FUU,
     *                   PIM_HE_Q2_STAT,     PIM_HE_Q2_SYS,
     *                   PIM_HE_Q2_FUUA,     PIM_HE_Q2_DIS


      INTEGER I_PIM_HE_Q2
C------- NE
      REAL*8 ::          PIM_NE_Q2_VALUE(9), PIM_NE_Q2_NU(9),
     *                   PIM_NE_Q2_Z(9),     PIM_NE_Q2_PT2(9),
     *                   PIM_NE_Q2_MULT(9),  PIM_NE_Q2_FUU(9),
     *                   PIM_NE_Q2_STAT(9),  PIM_NE_Q2_SYS(9),
     *                   PIM_NE_Q2_FUUA(9),  PIM_NE_Q2_DIS(9)

      COMMON /PIM_NE_Q2/ PIM_NE_Q2_VALUE,    PIM_NE_Q2_NU,
     *                   PIM_NE_Q2_Z,        PIM_NE_Q2_PT2,
     *                   PIM_NE_Q2_MULT,     PIM_NE_Q2_FUU,
     *                   PIM_NE_Q2_STAT,     PIM_NE_Q2_SYS,
     *                   PIM_NE_Q2_FUUA,     PIM_NE_Q2_DIS

      INTEGER I_PIM_NE_Q2
C------- KR
      REAL*8 ::          PIM_KR_Q2_VALUE(9), PIM_KR_Q2_NU(9),
     *                   PIM_KR_Q2_Z(9),     PIM_KR_Q2_PT2(9),
     *                   PIM_KR_Q2_MULT(9),  PIM_KR_Q2_FUU(9),
     *                   PIM_KR_Q2_STAT(9),  PIM_KR_Q2_SYS(9),
     *                   PIM_KR_Q2_FUUA(9),  PIM_KR_Q2_DIS(9)

      COMMON /PIM_KR_Q2/ PIM_KR_Q2_VALUE,    PIM_KR_Q2_NU,
     *                   PIM_KR_Q2_Z,        PIM_KR_Q2_PT2,
     *                   PIM_KR_Q2_MULT,     PIM_KR_Q2_FUU,
     *                   PIM_KR_Q2_STAT,     PIM_KR_Q2_SYS,
     *                   PIM_KR_Q2_FUUA,     PIM_KR_Q2_DIS

      INTEGER I_PIM_KR_Q2
C------- XE
      REAL*8 ::          PIM_XE_Q2_VALUE(9), PIM_XE_Q2_NU(9),
     *                   PIM_XE_Q2_Z(9),     PIM_XE_Q2_PT2(9),
     *                   PIM_XE_Q2_MULT(9),  PIM_XE_Q2_FUU(9),
     *                   PIM_XE_Q2_STAT(9),  PIM_XE_Q2_SYS(9),
     *                   PIM_XE_Q2_FUUA(9),  PIM_XE_Q2_DIS(9)

      COMMON /PIM_XE_Q2/ PIM_XE_Q2_VALUE,    PIM_XE_Q2_NU,
     *                   PIM_XE_Q2_Z,        PIM_XE_Q2_PT2,
     *                   PIM_XE_Q2_MULT,     PIM_XE_Q2_FUU,
     *                   PIM_XE_Q2_STAT,     PIM_XE_Q2_SYS,
     *                   PIM_XE_Q2_FUUA,     PIM_XE_Q2_DIS

      INTEGER I_PIM_XE_Q2
C-----------------------------------------
C------K- NU DEPENDENT DATA------------------------
C------- HE
      REAL*8 ::          KM_HE_NU_VALUE(9), KM_HE_NU_Z(9),
     *                   KM_HE_NU_Q2(9),    KM_HE_NU_PT2(9),
     *                   KM_HE_NU_MULT(9),  KM_HE_NU_FUU(9),
     *                   KM_HE_NU_STAT(9),  KM_HE_NU_SYS(9),
     *                   KM_HE_NU_FUUA(9),  KM_HE_NU_DIS(9)

      COMMON /KM_HE_NU/  KM_HE_NU_VALUE,    KM_HE_NU_Z,
     *                   KM_HE_NU_Q2,       KM_HE_NU_PT2,
     *                   KM_HE_NU_MULT,     KM_HE_NU_FUU,
     *                   KM_HE_NU_STAT,     KM_HE_NU_SYS,
     *                   KM_HE_NU_FUUA,     KM_HE_NU_DIS


      INTEGER I_KM_HE_NU
C------- NE
      REAL*8 ::          KM_NE_NU_VALUE(9), KM_NE_NU_Z(9),
     *                   KM_NE_NU_Q2(9),    KM_NE_NU_PT2(9),
     *                   KM_NE_NU_MULT(9),  KM_NE_NU_FUU(9),
     *                   KM_NE_NU_STAT(9),  KM_NE_NU_SYS(9),
     *                   KM_NE_NU_FUUA(9),  KM_NE_NU_DIS(9)

      COMMON /KM_NE_NU/  KM_NE_NU_VALUE,    KM_NE_NU_Z,
     *                   KM_NE_NU_Q2,       KM_NE_NU_PT2,
     *                   KM_NE_NU_MULT,     KM_NE_NU_FUU,
     *                   KM_NE_NU_STAT,     KM_NE_NU_SYS,
     *                   KM_NE_NU_FUUA,     KM_NE_NU_DIS

      INTEGER I_KM_NE_NU
C------- KR
      REAL*8 ::          KM_KR_NU_VALUE(9), KM_KR_NU_Z(9),
     *                   KM_KR_NU_Q2(9),    KM_KR_NU_PT2(9),
     *                   KM_KR_NU_MULT(9),  KM_KR_NU_FUU(9),
     *                   KM_KR_NU_STAT(9),  KM_KR_NU_SYS(9),
     *                   KM_KR_NU_FUUA(9),  KM_KR_NU_DIS(9)

      COMMON /KM_KR_NU/  KM_KR_NU_VALUE,    KM_KR_NU_Z,
     *                   KM_KR_NU_Q2,       KM_KR_NU_PT2,
     *                   KM_KR_NU_MULT,     KM_KR_NU_FUU,
     *                   KM_KR_NU_STAT,     KM_KR_NU_SYS,
     *                   KM_KR_NU_FUUA,     KM_KR_NU_DIS

      INTEGER I_KM_KR_NU
C------- XE
      REAL*8 ::          KM_XE_NU_VALUE(9), KM_XE_NU_Z(9),
     *                   KM_XE_NU_Q2(9),    KM_XE_NU_PT2(9),
     *                   KM_XE_NU_MULT(9),  KM_XE_NU_FUU(9),
     *                   KM_XE_NU_STAT(9),  KM_XE_NU_SYS(9),
     *                   KM_XE_NU_FUUA(9),  KM_XE_NU_DIS(9)

      COMMON /KM_XE_NU/  KM_XE_NU_VALUE,    KM_XE_NU_Z,
     *                   KM_XE_NU_Q2,       KM_XE_NU_PT2,
     *                   KM_XE_NU_MULT,     KM_XE_NU_FUU,
     *                   KM_XE_NU_STAT,     KM_XE_NU_SYS,
     *                   KM_XE_NU_FUUA,     KM_XE_NU_DIS

      INTEGER I_KM_XE_NU

C------K- Z DEPENDENT DATA------------------------
C------- HE
      REAL*8 ::          KM_HE_Z_VALUE(9), KM_HE_Z_NU(9),
     *                   KM_HE_Z_Q2(9),    KM_HE_Z_PT2(9),
     *                   KM_HE_Z_MULT(9),  KM_HE_Z_FUU(9),
     *                   KM_HE_Z_STAT(9),  KM_HE_Z_SYS(9),
     *                   KM_HE_Z_FUUA(9),  KM_HE_Z_DIS(9)

      COMMON /KM_HE_Z/   KM_HE_Z_VALUE,    KM_HE_Z_NU,
     *                   KM_HE_Z_Q2,       KM_HE_Z_PT2,
     *                   KM_HE_Z_MULT,     KM_HE_Z_FUU,
     *                   KM_HE_Z_STAT,     KM_HE_Z_SYS,
     *                   KM_HE_Z_FUUA,     KM_HE_Z_DIS


      INTEGER I_KM_HE_Z
C------- NE
      REAL*8 ::          KM_NE_Z_VALUE(9), KM_NE_Z_NU(9),
     *                   KM_NE_Z_Q2(9),    KM_NE_Z_PT2(9),
     *                   KM_NE_Z_MULT(9),  KM_NE_Z_FUU(9),
     *                   KM_NE_Z_STAT(9),  KM_NE_Z_SYS(9),
     *                   KM_NE_Z_FUUA(9),  KM_NE_Z_DIS(9)

      COMMON /KM_NE_Z/   KM_NE_Z_VALUE,    KM_NE_Z_NU,
     *                   KM_NE_Z_Q2,       KM_NE_Z_PT2,
     *                   KM_NE_Z_MULT,     KM_NE_Z_FUU,
     *                   KM_NE_Z_STAT,     KM_NE_Z_SYS,
     *                   KM_NE_Z_FUUA,     KM_NE_Z_DIS

      INTEGER I_KM_NE_Z
C------- KR
      REAL*8 ::          KM_KR_Z_VALUE(9), KM_KR_Z_NU(9),
     *                   KM_KR_Z_Q2(9),    KM_KR_Z_PT2(9),
     *                   KM_KR_Z_MULT(9),  KM_KR_Z_FUU(9),
     *                   KM_KR_Z_STAT(9),  KM_KR_Z_SYS(9),
     *                   KM_KR_Z_FUUA(9),  KM_KR_Z_DIS(9)

      COMMON /KM_KR_Z/   KM_KR_Z_VALUE,    KM_KR_Z_NU,
     *                   KM_KR_Z_Q2,       KM_KR_Z_PT2,
     *                   KM_KR_Z_MULT,     KM_KR_Z_FUU,
     *                   KM_KR_Z_STAT,     KM_KR_Z_SYS,
     *                   KM_KR_Z_FUUA,     KM_KR_Z_DIS

      INTEGER I_KM_KR_Z
C------- XE
      REAL*8 ::          KM_XE_Z_VALUE(9), KM_XE_Z_NU(9),
     *                   KM_XE_Z_Q2(9),    KM_XE_Z_PT2(9),
     *                   KM_XE_Z_MULT(9),  KM_XE_Z_FUU(9),
     *                   KM_XE_Z_STAT(9),  KM_XE_Z_SYS(9),
     *                   KM_XE_Z_FUUA(9),  KM_XE_Z_DIS(9)

      COMMON /KM_XE_NU/  KM_XE_Z_VALUE,    KM_XE_Z_NU,
     *                   KM_XE_Z_Q2,       KM_XE_Z_PT2,
     *                   KM_XE_Z_MULT,     KM_XE_Z_FUU,
     *                   KM_XE_Z_STAT,     KM_XE_Z_SYS,
     *                   KM_XE_Z_FUUA,     KM_XE_Z_DIS

      INTEGER I_KM_XE_Z

C------K- Q2 DEPENDENT DATA------------------------
C------- HE
      REAL*8 ::          KM_HE_Q2_VALUE(9), KM_HE_Q2_NU(9),
     *                   KM_HE_Q2_Z(9),     KM_HE_Q2_PT2(9),
     *                   KM_HE_Q2_MULT(9),  KM_HE_Q2_FUU(9),
     *                   KM_HE_Q2_STAT(9),  KM_HE_Q2_SYS(9),
     *                   KM_HE_Q2_FUUA(9),  KM_HE_Q2_DIS(9)

      COMMON /KM_HE_Q2/  KM_HE_Q2_VALUE,    KM_HE_Q2_NU,
     *                   KM_HE_Q2_Z,        KM_HE_Q2_PT2,
     *                   KM_HE_Q2_MULT,     KM_HE_Q2_FUU,
     *                   KM_HE_Q2_STAT,     KM_HE_Q2_SYS,
     *                   KM_HE_Q2_FUUA,     KM_HE_Q2_DIS

      INTEGER I_KM_HE_Q2
C------- NE
      REAL*8 ::          KM_NE_Q2_VALUE(9), KM_NE_Q2_NU(9),
     *                   KM_NE_Q2_Z(9),     KM_NE_Q2_PT2(9),
     *                   KM_NE_Q2_MULT(9),  KM_NE_Q2_FUU(9),
     *                   KM_NE_Q2_STAT(9),  KM_NE_Q2_SYS(9),
     *                   KM_NE_Q2_FUUA(9),  KM_NE_Q2_DIS(9)

      COMMON /KM_NE_Q2/  KM_NE_Q2_VALUE,    KM_NE_Q2_NU,
     *                   KM_NE_Q2_Z,        KM_NE_Q2_PT2,
     *                   KM_NE_Q2_MULT,     KM_NE_Q2_FUU,
     *                   KM_NE_Q2_STAT,     KM_NE_Q2_SYS,
     *                   KM_NE_Q2_FUUA,     KM_NE_Q2_DIS

      INTEGER I_KM_NE_Q2
C------- KR
      REAL*8 ::          KM_KR_Q2_VALUE(9), KM_KR_Q2_NU(9),
     *                   KM_KR_Q2_Z(9),     KM_KR_Q2_PT2(9),
     *                   KM_KR_Q2_MULT(9),  KM_KR_Q2_FUU(9),
     *                   KM_KR_Q2_STAT(9),  KM_KR_Q2_SYS(9),
     *                   KM_KR_Q2_FUUA(9),  KM_KR_Q2_DIS(9)

      COMMON /KM_KR_Q2/  KM_KR_Q2_VALUE,    KM_KR_Q2_NU,
     *                   KM_KR_Q2_Z,        KM_KR_Q2_PT2,
     *                   KM_KR_Q2_MULT,     KM_KR_Q2_FUU,
     *                   KM_KR_Q2_STAT,     KM_KR_Q2_SYS,
     *                   KM_KR_Q2_FUUA,     KM_KR_Q2_DIS

      INTEGER I_KM_KR_Q2
C------- X
      REAL*8 ::          KM_XE_Q2_VALUE(9), KM_XE_Q2_NU(9),
     *                   KM_XE_Q2_Z(9),     KM_XE_Q2_PT2(9),
     *                   KM_XE_Q2_MULT(9),  KM_XE_Q2_FUU(9),
     *                   KM_XE_Q2_STAT(9),  KM_XE_Q2_SYS(9),
     *                   KM_XE_Q2_FUUA(9),  KM_XE_Q2_DIS(9)

      COMMON /KM_XE_Q2/  KM_XE_Q2_VALUE,    KM_XE_Q2_NU,
     *                   KM_XE_Q2_Z,        KM_XE_Q2_PT2,
     *                   KM_XE_Q2_MULT,     KM_XE_Q2_FUU,
     *                   KM_XE_Q2_STAT,     KM_XE_Q2_SYS,
     *                   KM_XE_Q2_FUUA,     KM_XE_Q2_DIS

      INTEGER I_KM_XE_Q2

C - - - - - - - - - - - - - - - - - - - - -
C-----------------------------------------
C------PI+ PT2 DATA------------------------
C------- HE
      REAL*8 ::           PIP_HE_PT2_VALUE(9), PIP_HE_PT2_NU(9),
     *                    PIP_HE_PT2_Z(9),     PIP_HE_PT2_Q2(9),
     *                    PIP_HE_PT2_MULT(9),  PIP_HE_PT2_FUU(9),
     *                    PIP_HE_PT2_STAT(9),  PIP_HE_PT2_SYS(9),
     *                    PIP_HE_PT2_FUUA(9),  PIP_HE_PT2_DIS(9)

      COMMON /PIP_HE_PT2/ PIP_HE_PT2_VALUE,    PIP_HE_PT2_NU,
     *                    PIP_HE_PT2_Z,        PIP_HE_PT2_Q2,
     *                    PIP_HE_PT2_MULT,     PIP_HE_PT2_FUU,
     *                    PIP_HE_PT2_STAT,     PIP_HE_PT2_SYS,
     *                    PIP_HE_PT2_FUUA,     PIP_HE_PT2_DIS

      INTEGER I_PIP_HE_PT2

C-------NE
      REAL*8 ::           PIP_NE_PT2_VALUE(9), PIP_NE_PT2_NU(9),
     *                    PIP_NE_PT2_Z(9),     PIP_NE_PT2_Q2(9),
     *                    PIP_NE_PT2_MULT(9),  PIP_NE_PT2_FUU(9),
     *                    PIP_NE_PT2_STAT(9),  PIP_NE_PT2_SYS(9),
     *                    PIP_NE_PT2_FUUA(9),  PIP_NE_PT2_DIS(9)

      COMMON /PIP_NE_PT2/ PIP_NE_PT2_VALUE,    PIP_NE_PT2_NU,
     *                    PIP_NE_PT2_Z,        PIP_NE_PT2_Q2,
     *                    PIP_NE_PT2_MULT,     PIP_NE_PT2_FUU,
     *                    PIP_NE_PT2_STAT,     PIP_NE_PT2_SYS,
     *                    PIP_NE_PT2_FUUA,     PIP_NE_PT2_DIS

      INTEGER I_PIP_NE_PT2
C------- KR
      REAL*8 ::           PIP_KR_PT2_VALUE(9), PIP_KR_PT2_NU(9),
     *                    PIP_KR_PT2_Z(9),     PIP_KR_PT2_Q2(9),
     *                    PIP_KR_PT2_MULT(9),  PIP_KR_PT2_FUU(9),
     *                    PIP_KR_PT2_STAT(9),  PIP_KR_PT2_SYS(9),
     *                    PIP_KR_PT2_FUUA(9),  PIP_KR_PT2_DIS(9)

      COMMON /PIP_KR_PT2/ PIP_KR_PT2_VALUE,    PIP_KR_PT2_NU,
     *                    PIP_KR_PT2_Z,        PIP_KR_PT2_Q2,
     *                    PIP_KR_PT2_MULT,     PIP_KR_PT2_FUU,
     *                    PIP_KR_PT2_STAT,     PIP_KR_PT2_SYS,
     *                    PIP_KR_PT2_FUUA,     PIP_KR_PT2_DIS

      INTEGER I_PIP_KR_PT2
C------- XE
      REAL*8 ::           PIP_XE_PT2_VALUE(9), PIP_XE_PT2_NU(9),
     *                    PIP_XE_PT2_Z(9),     PIP_XE_PT2_Q2(9),
     *                    PIP_XE_PT2_MULT(9),  PIP_XE_PT2_FUU(9),
     *                    PIP_XE_PT2_STAT(9),  PIP_XE_PT2_SYS(9),
     *                    PIP_XE_PT2_FUUA(9),  PIP_XE_PT2_DIS(9)

      COMMON /PIP_XE_PT2/ PIP_XE_PT2_VALUE,    PIP_XE_PT2_NU,
     *                    PIP_XE_PT2_Z,        PIP_XE_PT2_Q2,
     *                    PIP_XE_PT2_MULT,     PIP_XE_PT2_FUU,
     *                    PIP_XE_PT2_STAT,     PIP_XE_PT2_SYS,
     *                    PIP_XE_PT2_FUUA,     PIP_XE_PT2_DIS

      INTEGER I_PIP_XE_PT2
C------K+ PT2 DATA------------------------
C------- HE
      REAL*8 ::           KP_HE_PT2_VALUE(9), KP_HE_PT2_NU(9),
     *                    KP_HE_PT2_Z(9),     KP_HE_PT2_Q2(9),
     *                    KP_HE_PT2_MULT(9),  KP_HE_PT2_FUU(9),
     *                    KP_HE_PT2_STAT(9),  KP_HE_PT2_SYS(9),
     *                    KP_HE_PT2_FUUA(9),  KP_HE_PT2_DIS(9)

      COMMON /KP_HE_PT2/  KP_HE_PT2_VALUE,    KP_HE_PT2_NU,
     *                    KP_HE_PT2_Z,        KP_HE_PT2_Q2,
     *                    KP_HE_PT2_MULT,     KP_HE_PT2_FUU,
     *                    KP_HE_PT2_STAT,     KP_HE_PT2_SYS,
     *                    KP_HE_PT2_FUUA,     KP_HE_PT2_DIS

      INTEGER I_KP_HE_PT2

C-------NE
      REAL*8 ::           KP_NE_PT2_VALUE(9), KP_NE_PT2_NU(9),
     *                    KP_NE_PT2_Z(9),     KP_NE_PT2_Q2(9),
     *                    KP_NE_PT2_MULT(9),  KP_NE_PT2_FUU(9),
     *                    KP_NE_PT2_STAT(9),  KP_NE_PT2_SYS(9),
     *                    KP_NE_PT2_FUUA(9),  KP_NE_PT2_DIS(9)

      COMMON /KP_NE_PT2/  KP_NE_PT2_VALUE,    KP_NE_PT2_NU,
     *                    KP_NE_PT2_Z,        KP_NE_PT2_Q2,
     *                    KP_NE_PT2_MULT,     KP_NE_PT2_FUU,
     *                    KP_NE_PT2_STAT,     KP_NE_PT2_SYS,
     *                    KP_NE_PT2_FUUA,     KP_NE_PT2_DIS

      INTEGER I_KP_NE_PT2
C------- KR
      REAL*8 ::           KP_KR_PT2_VALUE(9), KP_KR_PT2_NU(9),
     *                    KP_KR_PT2_Z(9),     KP_KR_PT2_Q2(9),
     *                    KP_KR_PT2_MULT(9),  KP_KR_PT2_FUU(9),
     *                    KP_KR_PT2_STAT(9),  KP_KR_PT2_SYS(9),
     *                    KP_KR_PT2_FUUA(9),  KP_KR_PT2_DIS(9)

      COMMON /KP_KR_PT2/  KP_KR_PT2_VALUE,    KP_KR_PT2_NU,
     *                    KP_KR_PT2_Z,        KP_KR_PT2_Q2,
     *                    KP_KR_PT2_MULT,     KP_KR_PT2_FUU,
     *                    KP_KR_PT2_STAT,     KP_KR_PT2_SYS,
     *                    KP_KR_PT2_FUUA,     KP_KR_PT2_DIS

      INTEGER I_KP_KR_PT2
C------- XE
      REAL*8 ::           KP_XE_PT2_VALUE(9), KP_XE_PT2_NU(9),
     *                    KP_XE_PT2_Z(9),     KP_XE_PT2_Q2(9),
     *                    KP_XE_PT2_MULT(9),  KP_XE_PT2_FUU(9),
     *                    KP_XE_PT2_STAT(9),  KP_XE_PT2_SYS(9),
     *                    KP_XE_PT2_FUUA(9),  KP_XE_PT2_DIS(9)

      COMMON /KP_XE_PT2/  KP_XE_PT2_VALUE,    KP_XE_PT2_NU,
     *                    KP_XE_PT2_Z,        KP_XE_PT2_Q2,
     *                    KP_XE_PT2_MULT,     KP_XE_PT2_FUU,
     *                    KP_XE_PT2_STAT,     KP_XE_PT2_SYS,
     *                    KP_XE_PT2_FUUA,     KP_XE_PT2_DIS

      INTEGER I_KP_XE_PT2
C-----------------------------------------
C-----------------------------------------
C------PI- PT2 DATA------------------------
C------- HE
      REAL*8 ::           PIM_HE_PT2_VALUE(9), PIM_HE_PT2_NU(9),
     *                    PIM_HE_PT2_Z(9),     PIM_HE_PT2_Q2(9),
     *                    PIM_HE_PT2_MULT(9),  PIM_HE_PT2_FUU(9),
     *                    PIM_HE_PT2_STAT(9),  PIM_HE_PT2_SYS(9),
     *                    PIM_HE_PT2_FUUA(9),  PIM_HE_PT2_DIS(9)

      COMMON /PIM_HE_PT2/ PIM_HE_PT2_VALUE,    PIM_HE_PT2_NU,
     *                    PIM_HE_PT2_Z,        PIM_HE_PT2_Q2,
     *                    PIM_HE_PT2_MULT,     PIM_HE_PT2_FUU,
     *                    PIM_HE_PT2_STAT,     PIM_HE_PT2_SYS,
     *                    PIM_HE_PT2_FUUA,     PIM_HE_PT2_DIS

      INTEGER I_PIM_HE_PT2

C-------NE
      REAL*8 ::           PIM_NE_PT2_VALUE(9), PIM_NE_PT2_NU(9),
     *                    PIM_NE_PT2_Z(9),     PIM_NE_PT2_Q2(9),
     *                    PIM_NE_PT2_MULT(9),  PIM_NE_PT2_FUU(9),
     *                    PIM_NE_PT2_STAT(9),  PIM_NE_PT2_SYS(9),
     *                    PIM_NE_PT2_FUUA(9),  PIM_NE_PT2_DIS(9)

      COMMON /PIM_NE_PT2/ PIM_NE_PT2_VALUE,    PIM_NE_PT2_NU,
     *                    PIM_NE_PT2_Z,        PIM_NE_PT2_Q2,
     *                    PIM_NE_PT2_MULT,     PIM_NE_PT2_FUU,
     *                    PIM_NE_PT2_STAT,     PIM_NE_PT2_SYS,
     *                    PIM_NE_PT2_FUUA,     PIM_NE_PT2_DIS

      INTEGER I_PIM_NE_PT2
C------- KR
      REAL*8 ::           PIM_KR_PT2_VALUE(9), PIM_KR_PT2_NU(9),
     *                    PIM_KR_PT2_Z(9),     PIM_KR_PT2_Q2(9),
     *                    PIM_KR_PT2_MULT(9),  PIM_KR_PT2_FUU(9),
     *                    PIM_KR_PT2_STAT(9),  PIM_KR_PT2_SYS(9),
     *                    PIM_KR_PT2_FUUA(9),  PIM_KR_PT2_DIS(9)

      COMMON /PIM_KR_PT2/ PIM_KR_PT2_VALUE,    PIM_KR_PT2_NU,
     *                    PIM_KR_PT2_Z,        PIM_KR_PT2_Q2,
     *                    PIM_KR_PT2_MULT,     PIM_KR_PT2_FUU,
     *                    PIM_KR_PT2_STAT,     PIM_KR_PT2_SYS,
     *                    PIM_KR_PT2_FUUA,     PIM_KR_PT2_DIS

      INTEGER I_PIM_KR_PT2
C------- XE
      REAL*8 ::           PIM_XE_PT2_VALUE(9), PIM_XE_PT2_NU(9),
     *                    PIM_XE_PT2_Z(9),     PIM_XE_PT2_Q2(9),
     *                    PIM_XE_PT2_MULT(9),  PIM_XE_PT2_FUU(9),
     *                    PIM_XE_PT2_STAT(9),  PIM_XE_PT2_SYS(9),
     *                    PIM_XE_PT2_FUUA(9),  PIM_XE_PT2_DIS(9)

      COMMON /PIM_XE_PT2/ PIM_XE_PT2_VALUE,    PIM_XE_PT2_NU,
     *                    PIM_XE_PT2_Z,        PIM_XE_PT2_Q2,
     *                    PIM_XE_PT2_MULT,     PIM_XE_PT2_FUU,
     *                    PIM_XE_PT2_STAT,     PIM_XE_PT2_SYS,
     *                    PIM_XE_PT2_FUUA,     PIM_XE_PT2_DIS

      INTEGER I_PIM_XE_PT2
C------K- PT2 DATA------------------------
C------- HE
      REAL*8 ::           KM_HE_PT2_VALUE(9), KM_HE_PT2_NU(9),
     *                    KM_HE_PT2_Z(9),     KM_HE_PT2_Q2(9),
     *                    KM_HE_PT2_MULT(9),  KM_HE_PT2_FUU(9),
     *                    KM_HE_PT2_STAT(9),  KM_HE_PT2_SYS(9),
     *                    KM_HE_PT2_FUUA(9),  KM_HE_PT2_DIS(9)

      COMMON /KM_HE_PT2/  KM_HE_PT2_VALUE,    KM_HE_PT2_NU,
     *                    KM_HE_PT2_Z,        KM_HE_PT2_Q2,
     *                    KM_HE_PT2_MULT,     KM_HE_PT2_FUU,
     *                    KM_HE_PT2_STAT,     KM_HE_PT2_SYS,
     *                    KM_HE_PT2_FUUA,     KM_HE_PT2_DIS

      INTEGER I_KM_HE_PT2

C-------NE
      REAL*8 ::           KM_NE_PT2_VALUE(9), KM_NE_PT2_NU(9),
     *                    KM_NE_PT2_Z(9),     KM_NE_PT2_Q2(9),
     *                    KM_NE_PT2_MULT(9),  KM_NE_PT2_FUU(9),
     *                    KM_NE_PT2_STAT(9),  KM_NE_PT2_SYS(9),
     *                    KM_NE_PT2_FUUA(9),  KM_NE_PT2_DIS(9)

      COMMON /KM_NE_PT2/  KM_NE_PT2_VALUE,    KM_NE_PT2_NU,
     *                    KM_NE_PT2_Z,        KM_NE_PT2_Q2,
     *                    KM_NE_PT2_MULT,     KM_NE_PT2_FUU,
     *                    KM_NE_PT2_STAT,     KM_NE_PT2_SYS,
     *                    KM_NE_PT2_FUUA,     KM_NE_PT2_DIS

      INTEGER I_KM_NE_PT2
C------- KR
      REAL*8 ::           KM_KR_PT2_VALUE(9), KM_KR_PT2_NU(9),
     *                    KM_KR_PT2_Z(9),     KM_KR_PT2_Q2(9),
     *                    KM_KR_PT2_MULT(9),  KM_KR_PT2_FUU(9),
     *                    KM_KR_PT2_STAT(9),  KM_KR_PT2_SYS(9),
     *                    KM_KR_PT2_FUUA(9),  KM_KR_PT2_DIS(9)

      COMMON /KM_KR_PT2/  KM_KR_PT2_VALUE,    KM_KR_PT2_NU,
     *                    KM_KR_PT2_Z,        KM_KR_PT2_Q2,
     *                    KM_KR_PT2_MULT,     KM_KR_PT2_FUU,
     *                    KM_KR_PT2_STAT,     KM_KR_PT2_SYS,
     *                    KM_KR_PT2_FUUA,     KM_KR_PT2_DIS

      INTEGER I_KM_KR_PT2
C------- XE
      REAL*8 ::           KM_XE_PT2_VALUE(9), KM_XE_PT2_NU(9),
     *                    KM_XE_PT2_Z(9),     KM_XE_PT2_Q2(9),
     *                    KM_XE_PT2_MULT(9),  KM_XE_PT2_FUU(9),
     *                    KM_XE_PT2_STAT(9),  KM_XE_PT2_SYS(9),
     *                    KM_XE_PT2_FUUA(9),  KM_XE_PT2_DIS(9)

      COMMON /KM_XE_PT2/  KM_XE_PT2_VALUE,    KM_XE_PT2_NU,
     *                    KM_XE_PT2_Z,        KM_XE_PT2_Q2,
     *                    KM_XE_PT2_MULT,     KM_XE_PT2_FUU,
     *                    KM_XE_PT2_STAT,     KM_XE_PT2_SYS,
     *                    KM_XE_PT2_FUUA,     KM_XE_PT2_DIS

      INTEGER I_KM_XE_PT2

C-----------PI0------------------------
C-----------------------------------------
C------PI0 NU DEPENDENT DATA------------------------
C------- HE
      REAL*8 ::          PI0_HE_NU_VALUE(9), PI0_HE_NU_Z(9),
     *                   PI0_HE_NU_Q2(9),    PI0_HE_NU_PT2(9),
     *                   PI0_HE_NU_MULT(9),  PI0_HE_NU_FUU(9),
     *                   PI0_HE_NU_STAT(9),  PI0_HE_NU_SYS(9),
     *                   PI0_HE_NU_FUUA(9),  PI0_HE_NU_DIS(9)

      COMMON /PI0_HE_NU/ PI0_HE_NU_VALUE,    PI0_HE_NU_Z,
     *                   PI0_HE_NU_Q2,       PI0_HE_NU_PT2,
     *                   PI0_HE_NU_MULT,     PI0_HE_NU_FUU,
     *                   PI0_HE_NU_STAT,     PI0_HE_NU_SYS,
     *                   PI0_HE_NU_FUUA,     PI0_HE_NU_DIS

      INTEGER I_PI0_HE_NU
C------- NE
      REAL*8 ::          PI0_NE_NU_VALUE(9), PI0_NE_NU_Z(9),
     *                   PI0_NE_NU_Q2(9),    PI0_NE_NU_PT2(9),
     *                   PI0_NE_NU_MULT(9),  PI0_NE_NU_FUU(9),
     *                   PI0_NE_NU_STAT(9),  PI0_NE_NU_SYS(9),
     *                   PI0_NE_NU_FUUA(9),  PI0_NE_NU_DIS(9)

      COMMON /PI0_NE_NU/ PI0_NE_NU_VALUE,    PI0_NE_NU_Z,
     *                   PI0_NE_NU_Q2,       PI0_NE_NU_PT2,
     *                   PI0_NE_NU_MULT,     PI0_NE_NU_FUU,
     *                   PI0_NE_NU_STAT,     PI0_NE_NU_SYS,
     *                   PI0_NE_NU_FUUA,     PI0_NE_NU_DIS


      INTEGER I_PI0_NE_NU
C------- KR
      REAL*8 ::          PI0_KR_NU_VALUE(9), PI0_KR_NU_Z(9),
     *                   PI0_KR_NU_Q2(9),    PI0_KR_NU_PT2(9),
     *                   PI0_KR_NU_MULT(9),  PI0_KR_NU_FUU(9),
     *                   PI0_KR_NU_STAT(9),  PI0_KR_NU_SYS(9),
     *                   PI0_KR_NU_FUUA(9),  PI0_KR_NU_DIS(9)

      COMMON /PI0_KR_NU/ PI0_KR_NU_VALUE,    PI0_KR_NU_Z,
     *                   PI0_KR_NU_Q2,       PI0_KR_NU_PT2,
     *                   PI0_KR_NU_MULT,     PI0_KR_NU_FUU,
     *                   PI0_KR_NU_STAT,     PI0_KR_NU_SYS,
     *                   PI0_KR_NU_FUUA,     PI0_KR_NU_DIS

      INTEGER I_PI0_KR_NU
C------- XE
      REAL*8 ::          PI0_XE_NU_VALUE(9), PI0_XE_NU_Z(9),
     *                   PI0_XE_NU_Q2(9),    PI0_XE_NU_PT2(9),
     *                   PI0_XE_NU_MULT(9),  PI0_XE_NU_FUU(9),
     *                   PI0_XE_NU_STAT(9),  PI0_XE_NU_SYS(9),
     *                   PI0_XE_NU_FUUA(9),  PI0_XE_NU_DIS(9)

      COMMON /PI0_XE_NU/ PI0_XE_NU_VALUE,    PI0_XE_NU_Z,
     *                   PI0_XE_NU_Q2,       PI0_XE_NU_PT2,
     *                   PI0_XE_NU_MULT,     PI0_XE_NU_FUU,
     *                   PI0_XE_NU_STAT,     PI0_XE_NU_SYS,
     *                   PI0_XE_NU_FUUA,     PI0_XE_NU_DIS

      INTEGER I_PI0_XE_NU

C------PI+ Z DEPENDENT DATA------------------------
C------- HE
      REAL*8 ::          PI0_HE_Z_VALUE(9), PI0_HE_Z_NU(9),
     *                   PI0_HE_Z_Q2(9),    PI0_HE_Z_PT2(9),
     *                   PI0_HE_Z_MULT(9),  PI0_HE_Z_FUU(9),
     *                   PI0_HE_Z_STAT(9),  PI0_HE_Z_SYS(9),
     *                   PI0_HE_Z_FUUA(9),  PI0_HE_Z_DIS(9)

      COMMON /PI0_HE_Z/  PI0_HE_Z_VALUE,    PI0_HE_Z_NU,
     *                   PI0_HE_Z_Q2,       PI0_HE_Z_PT2,
     *                   PI0_HE_Z_MULT,     PI0_HE_Z_FUU,
     *                   PI0_HE_Z_STAT,     PI0_HE_Z_SYS,
     *                   PI0_HE_Z_FUUA,     PI0_HE_Z_DIS


      INTEGER I_PI0_HE_Z
C------- NE
      REAL*8 ::          PI0_NE_Z_VALUE(9), PI0_NE_Z_NU(9),
     *                   PI0_NE_Z_Q2(9),    PI0_NE_Z_PT2(9),
     *                   PI0_NE_Z_MULT(9),  PI0_NE_Z_FUU(9),
     *                   PI0_NE_Z_STAT(9),  PI0_NE_Z_SYS(9),
     *                   PI0_NE_Z_FUUA(9),  PI0_NE_Z_DIS(9)

      COMMON /PI0_NE_Z/  PI0_NE_Z_VALUE,    PI0_NE_Z_NU,
     *                   PI0_NE_Z_Q2,       PI0_NE_Z_PT2,
     *                   PI0_NE_Z_MULT,     PI0_NE_Z_FUU,
     *                   PI0_NE_Z_STAT,     PI0_NE_Z_SYS,
     *                   PI0_NE_Z_FUUA,     PI0_NE_Z_DIS

      INTEGER I_PI0_NE_Z
C------- KR
      REAL*8 ::          PI0_KR_Z_VALUE(9), PI0_KR_Z_NU(9),
     *                   PI0_KR_Z_Q2(9),    PI0_KR_Z_PT2(9),
     *                   PI0_KR_Z_MULT(9),  PI0_KR_Z_FUU(9),
     *                   PI0_KR_Z_STAT(9),  PI0_KR_Z_SYS(9),
     *                   PI0_KR_Z_FUUA(9),  PI0_KR_Z_DIS(9)

      COMMON /PI0_KR_Z/  PI0_KR_Z_VALUE,    PI0_KR_Z_NU,
     *                   PI0_KR_Z_Q2,       PI0_KR_Z_PT2,
     *                   PI0_KR_Z_MULT,     PI0_KR_Z_FUU,
     *                   PI0_KR_Z_STAT,     PI0_KR_Z_SYS,
     *                   PI0_KR_Z_FUUA,     PI0_KR_Z_DIS

      INTEGER I_PI0_KR_Z
C------- XE
      REAL*8 ::          PI0_XE_Z_VALUE(9), PI0_XE_Z_NU(9),
     *                   PI0_XE_Z_Q2(9),    PI0_XE_Z_PT2(9),
     *                   PI0_XE_Z_MULT(9),  PI0_XE_Z_FUU(9),
     *                   PI0_XE_Z_STAT(9),  PI0_XE_Z_SYS(9),
     *                   PI0_XE_Z_FUUA(9),  PI0_XE_Z_DIS(9)

      COMMON /PI0_XE_NU/ PI0_XE_Z_VALUE,    PI0_XE_Z_NU,
     *                   PI0_XE_Z_Q2,       PI0_XE_Z_PT2,
     *                   PI0_XE_Z_MULT,     PI0_XE_Z_FUU,
     *                   PI0_XE_Z_STAT,     PI0_XE_Z_SYS,
     *                   PI0_XE_Z_FUUA,     PI0_XE_Z_DIS

      INTEGER I_PI0_XE_Z

C------PI+ Q2 DEPENDENT DATA------------------------
C------- HE
      REAL*8 ::          PI0_HE_Q2_VALUE(9), PI0_HE_Q2_NU(9),
     *                   PI0_HE_Q2_Z(9),     PI0_HE_Q2_PT2(9),
     *                   PI0_HE_Q2_MULT(9),  PI0_HE_Q2_FUU(9),
     *                   PI0_HE_Q2_STAT(9),  PI0_HE_Q2_SYS(9),
     *                   PI0_HE_Q2_FUUA(9),  PI0_HE_Q2_DIS(9)

      COMMON /PI0_HE_Q2/ PI0_HE_Q2_VALUE,    PI0_HE_Q2_NU,
     *                   PI0_HE_Q2_Z,        PI0_HE_Q2_PT2,
     *                   PI0_HE_Q2_MULT,     PI0_HE_Q2_FUU,
     *                   PI0_HE_Q2_STAT,     PI0_HE_Q2_SYS,
     *                   PI0_HE_Q2_FUUA,     PI0_HE_Q2_DIS

      INTEGER I_PI0_HE_Q2
C------- NE
      REAL*8 ::          PI0_NE_Q2_VALUE(9), PI0_NE_Q2_NU(9),
     *                   PI0_NE_Q2_Z(9),     PI0_NE_Q2_PT2(9),
     *                   PI0_NE_Q2_MULT(9),  PI0_NE_Q2_FUU(9),
     *                   PI0_NE_Q2_STAT(9),  PI0_NE_Q2_SYS(9),
     *                   PI0_NE_Q2_FUUA(9),  PI0_NE_Q2_DIS(9)

      COMMON /PI0_NE_Q2/ PI0_NE_Q2_VALUE,    PI0_NE_Q2_NU,
     *                   PI0_NE_Q2_Z,        PI0_NE_Q2_PT2,
     *                   PI0_NE_Q2_MULT,     PI0_NE_Q2_FUU,
     *                   PI0_NE_Q2_STAT,     PI0_NE_Q2_SYS,
     *                   PI0_NE_Q2_FUUA,     PI0_NE_Q2_DIS


      INTEGER I_PI0_NE_Q2
C------- KR
      REAL*8 ::          PI0_KR_Q2_VALUE(9), PI0_KR_Q2_NU(9),
     *                   PI0_KR_Q2_Z(9),     PI0_KR_Q2_PT2(9),
     *                   PI0_KR_Q2_MULT(9),  PI0_KR_Q2_FUU(9),
     *                   PI0_KR_Q2_STAT(9),  PI0_KR_Q2_SYS(9),
     *                   PI0_KR_Q2_FUUA(9),  PI0_KR_Q2_DIS(9)

      COMMON /PI0_KR_Q2/ PI0_KR_Q2_VALUE,    PI0_KR_Q2_NU,
     *                   PI0_KR_Q2_Z,        PI0_KR_Q2_PT2,
     *                   PI0_KR_Q2_MULT,     PI0_KR_Q2_FUU,
     *                   PI0_KR_Q2_STAT,     PI0_KR_Q2_SYS,
     *                   PI0_KR_Q2_FUUA,     PI0_KR_Q2_DIS

      INTEGER I_PI0_KR_Q2
C------- XE
      REAL*8 ::          PI0_XE_Q2_VALUE(9), PI0_XE_Q2_NU(9),
     *                   PI0_XE_Q2_Z(9),     PI0_XE_Q2_PT2(9),
     *                   PI0_XE_Q2_MULT(9),  PI0_XE_Q2_FUU(9),
     *                   PI0_XE_Q2_STAT(9),  PI0_XE_Q2_SYS(9),
     *                   PI0_XE_Q2_FUUA(9),  PI0_XE_Q2_DIS(9)

      COMMON /PI0_XE_Q2/ PI0_XE_Q2_VALUE,    PI0_XE_Q2_NU,
     *                   PI0_XE_Q2_Z,        PI0_XE_Q2_PT2,
     *                   PI0_XE_Q2_MULT,     PI0_XE_Q2_FUU,
     *                   PI0_XE_Q2_STAT,     PI0_XE_Q2_SYS,
     *                   PI0_XE_Q2_FUUA,     PI0_XE_Q2_DIS

      INTEGER I_PI0_XE_Q2

C------PI+ PT2 DATA------------------------
C------- HE
      REAL*8 ::           PI0_HE_PT2_VALUE(9), PI0_HE_PT2_NU(9),
     *                    PI0_HE_PT2_Z(9),     PI0_HE_PT2_Q2(9),
     *                    PI0_HE_PT2_MULT(9),  PI0_HE_PT2_FUU(9),
     *                    PI0_HE_PT2_STAT(9),  PI0_HE_PT2_SYS(9),
     *                    PI0_HE_PT2_FUUA(9),  PI0_HE_PT2_DIS(9)

      COMMON /PI0_HE_PT2/ PI0_HE_PT2_VALUE,    PI0_HE_PT2_NU,
     *                    PI0_HE_PT2_Z,        PI0_HE_PT2_Q2,
     *                    PI0_HE_PT2_MULT,     PI0_HE_PT2_FUU,
     *                    PI0_HE_PT2_STAT,     PI0_HE_PT2_SYS,
     *                    PI0_HE_PT2_FUUA,     PI0_HE_PT2_DIS

      INTEGER I_PI0_HE_PT2

C-------NE
      REAL*8 ::           PI0_NE_PT2_VALUE(9), PI0_NE_PT2_NU(9),
     *                    PI0_NE_PT2_Z(9),     PI0_NE_PT2_Q2(9),
     *                    PI0_NE_PT2_MULT(9),  PI0_NE_PT2_FUU(9),
     *                    PI0_NE_PT2_STAT(9),  PI0_NE_PT2_SYS(9),
     *                    PI0_NE_PT2_FUUA(9),  PI0_NE_PT2_DIS(9)

      COMMON /PI0_NE_PT2/ PI0_NE_PT2_VALUE,    PI0_NE_PT2_NU,
     *                    PI0_NE_PT2_Z,        PI0_NE_PT2_Q2,
     *                    PI0_NE_PT2_MULT,     PI0_NE_PT2_FUU,
     *                    PI0_NE_PT2_STAT,     PI0_NE_PT2_SYS,
     *                    PI0_NE_PT2_FUUA,     PI0_NE_PT2_DIS

      INTEGER I_PI0_NE_PT2
C------- KR
      REAL*8 ::           PI0_KR_PT2_VALUE(9), PI0_KR_PT2_NU(9),
     *                    PI0_KR_PT2_Z(9),     PI0_KR_PT2_Q2(9),
     *                    PI0_KR_PT2_MULT(9),  PI0_KR_PT2_FUU(9),
     *                    PI0_KR_PT2_STAT(9),  PI0_KR_PT2_SYS(9),
     *                    PI0_KR_PT2_FUUA(9),  PI0_KR_PT2_DIS(9)

      COMMON /PI0_KR_PT2/ PI0_KR_PT2_VALUE,    PI0_KR_PT2_NU,
     *                    PI0_KR_PT2_Z,        PI0_KR_PT2_Q2,
     *                    PI0_KR_PT2_MULT,     PI0_KR_PT2_FUU,
     *                    PI0_KR_PT2_STAT,     PI0_KR_PT2_SYS,
     *                    PI0_KR_PT2_FUUA,     PI0_KR_PT2_DIS

      INTEGER I_PI0_KR_PT2
C------- XE
      REAL*8 ::           PI0_XE_PT2_VALUE(9), PI0_XE_PT2_NU(9),
     *                    PI0_XE_PT2_Z(9),     PI0_XE_PT2_Q2(9),
     *                    PI0_XE_PT2_MULT(9),  PI0_XE_PT2_FUU(9),
     *                    PI0_XE_PT2_STAT(9),  PI0_XE_PT2_SYS(9),
     *                    PI0_XE_PT2_FUUA(9),  PI0_XE_PT2_DIS(9)

      COMMON /PI0_XE_PT2/ PI0_XE_PT2_VALUE,    PI0_XE_PT2_NU,
     *                    PI0_XE_PT2_Z,        PI0_XE_PT2_Q2,
     *                    PI0_XE_PT2_MULT,     PI0_XE_PT2_FUU,
     *                    PI0_XE_PT2_STAT,     PI0_XE_PT2_SYS,
     *                    PI0_XE_PT2_FUUA,     PI0_XE_PT2_DIS

      INTEGER I_PI0_XE_PT2

C--------JLAB DATA
      REAL*8 ::           CLAS1_PTLOW(22),  CLAS1_PTHIGH(22),
     *                    CLAS1_NULOW(22),  CLAS1_NUHIGH(22),
     *                    CLAS1_Q2LOW(22),  CLAS1_Q2HIGH(22),
     *                    CLAS1_RA_C(22),   CLAS1_RA_Cerr(22),
     *                    CLAS1_RA_FE(22),  CLAS1_RA_FEerr(22),
     *                    CLAS1_RA_PB(22),  CLAS1_RA_PBerr(22),
     *                    CLAS1_C_FUU(22),  CLAS1_C_FUUA(22),
     *                    CLAS1_FE_FUU(22), CLAS1_FE_FUUA(22),
     *                    CLAS1_PB_FUU(22), CLAS1_PB_FUUA(22),
     *                    CLAS1_C_DIS(22),  CLAS1_FE_DIS(22),
     *                    CLAS1_PB_DIS(22)

      COMMON /CLAS1/      CLAS1_PTLOW,  CLAS1_PTHIGH,
     *                    CLAS1_NULOW,  CLAS1_NUHIGH,
     *                    CLAS1_Q2LOW,  CLAS1_Q2HIGH,
     *                    CLAS1_RA_C,   CLAS1_RA_Cerr,
     *                    CLAS1_RA_FE,  CLAS1_RA_FEerr,
     *                    CLAS1_RA_PB,  CLAS1_RA_PBerr,
     *                    CLAS1_C_FUU,  CLAS1_C_FUUA,
     *                    CLAS1_FE_FUU, CLAS1_FE_FUUA,
     *                    CLAS1_PB_FUU, CLAS1_PB_FUUA,
     *                    CLAS1_C_DIS,  CLAS1_FE_DIS,
     *                    CLAS1_PB_DIS

      INTEGER I_CLAS1

      REAL*8 ::           CLAS2_PTLOW(22), CLAS2_PTHIGH(22),
     *                    CLAS2_NULOW(22), CLAS2_NUHIGH(22),
     *                    CLAS2_Q2LOW(22), CLAS2_Q2HIGH(22),
     *                    CLAS2_RA_C(22),  CLAS2_RA_Cerr(22),
     *                    CLAS2_RA_FE(22),  CLAS2_RA_FEerr(22),
     *                    CLAS2_RA_PB(22),  CLAS2_RA_PBerr(22),
     *                    CLAS2_C_FUU(22),  CLAS2_C_FUUA(22),
     *                    CLAS2_FE_FUU(22), CLAS2_FE_FUUA(22),
     *                    CLAS2_PB_FUU(22), CLAS2_PB_FUUA(22),
     *                    CLAS2_C_DIS(22),  CLAS2_FE_DIS(22),
     *                    CLAS2_PB_DIS(22)

      COMMON /CLAS2/      CLAS2_PTLOW, CLAS2_PTHIGH,
     *                    CLAS2_NULOW, CLAS2_NUHIGH,
     *                    CLAS2_Q2LOW, CLAS2_Q2HIGH,
     *                    CLAS2_RA_C,  CLAS2_RA_Cerr,
     *                    CLAS2_RA_FE,  CLAS2_RA_FEerr,
     *                    CLAS2_RA_PB,  CLAS2_RA_PBerr,
     *                    CLAS2_C_FUU,  CLAS2_C_FUUA,
     *                    CLAS2_FE_FUU, CLAS2_FE_FUUA,
     *                    CLAS2_PB_FUU, CLAS2_PB_FUUA,
     *                    CLAS2_C_DIS,  CLAS2_FE_DIS,
     *                    CLAS2_PB_DIS

      INTEGER I_CLAS2

      REAL*8 ::           CLAS3_PTLOW(22),  CLAS3_PTHIGH(22),
     *                    CLAS3_NULOW(22),  CLAS3_NUHIGH(22),
     *                    CLAS3_Q2LOW(22),  CLAS3_Q2HIGH(22),
     *                    CLAS3_RA_C(22),   CLAS3_RA_Cerr(22),
     *                    CLAS3_RA_FE(22),  CLAS3_RA_FEerr(22),
     *                    CLAS3_RA_PB(22),  CLAS3_RA_PBerr(22),
     *                    CLAS3_C_FUU(22),  CLAS3_C_FUUA(22),
     *                    CLAS3_FE_FUU(22), CLAS3_FE_FUUA(22),
     *                    CLAS3_PB_FUU(22), CLAS3_PB_FUUA(22),
     *                    CLAS3_C_DIS(22),  CLAS3_FE_DIS(22),
     *                    CLAS3_PB_DIS(22)

      COMMON /CLAS3/      CLAS3_PTLOW,  CLAS3_PTHIGH,
     *                    CLAS3_NULOW,  CLAS3_NUHIGH,
     *                    CLAS3_Q2LOW,  CLAS3_Q2HIGH,
     *                    CLAS3_RA_C,   CLAS3_RA_Cerr,
     *                    CLAS3_RA_FE,  CLAS3_RA_FEerr,
     *                    CLAS3_RA_PB,  CLAS3_RA_PBerr,
     *                    CLAS3_C_FUU,  CLAS3_C_FUUA,
     *                    CLAS3_FE_FUU, CLAS3_FE_FUUA,
     *                    CLAS3_PB_FUU, CLAS3_PB_FUUA,
     *                    CLAS3_C_DIS,  CLAS3_FE_DIS,
     *                    CLAS3_PB_DIS

      INTEGER I_CLAS3

      REAL*8 ::           CLAS4_PTLOW(22),  CLAS4_PTHIGH(22),
     *                    CLAS4_NULOW(22),  CLAS4_NUHIGH(22),
     *                    CLAS4_Q2LOW(22),  CLAS4_Q2HIGH(22),
     *                    CLAS4_RA_C(22),   CLAS4_RA_Cerr(22),
     *                    CLAS4_RA_FE(22),  CLAS4_RA_FEerr(22),
     *                    CLAS4_RA_PB(22),  CLAS4_RA_PBerr(22),
     *                    CLAS4_C_FUU(22),  CLAS4_C_FUUA(22),
     *                    CLAS4_FE_FUU(22), CLAS4_FE_FUUA(22),
     *                    CLAS4_PB_FUU(22), CLAS4_PB_FUUA(22),
     *                    CLAS4_C_DIS(22),  CLAS4_FE_DIS(22),
     *                    CLAS4_PB_DIS(22)

      COMMON /CLAS4/      CLAS4_PTLOW,  CLAS4_PTHIGH,
     *                    CLAS4_NULOW,  CLAS4_NUHIGH,
     *                    CLAS4_Q2LOW,  CLAS4_Q2HIGH,
     *                    CLAS4_RA_C,   CLAS4_RA_Cerr,
     *                    CLAS4_RA_FE,  CLAS4_RA_FEerr,
     *                    CLAS4_RA_PB,  CLAS4_RA_PBerr,
     *                    CLAS4_C_FUU,  CLAS4_C_FUUA,
     *                    CLAS4_FE_FUU, CLAS4_FE_FUUA,
     *                    CLAS4_PB_FUU, CLAS4_PB_FUUA,
     *                    CLAS4_C_DIS,  CLAS4_FE_DIS,
     *                    CLAS4_PB_DIS


      INTEGER I_CLAS4

      REAL*8 ::           CLAS5_PTLOW(22),  CLAS5_PTHIGH(22),
     *                    CLAS5_NULOW(22),  CLAS5_NUHIGH(22),
     *                    CLAS5_Q2LOW(22),  CLAS5_Q2HIGH(22),
     *                    CLAS5_RA_C(22),   CLAS5_RA_Cerr(22),
     *                    CLAS5_RA_FE(22),  CLAS5_RA_FEerr(22),
     *                    CLAS5_RA_PB(22),  CLAS5_RA_PBerr(22),
     *                    CLAS5_C_FUU(22),  CLAS5_C_FUUA(22),
     *                    CLAS5_FE_FUU(22), CLAS5_FE_FUUA(22),
     *                    CLAS5_PB_FUU(22), CLAS5_PB_FUUA(22),
     *                    CLAS5_C_DIS(22),  CLAS5_FE_DIS(22),
     *                    CLAS5_PB_DIS(22)

      COMMON /CLAS5/      CLAS5_PTLOW,  CLAS5_PTHIGH,
     *                    CLAS5_NULOW,  CLAS5_NUHIGH,
     *                    CLAS5_Q2LOW,  CLAS5_Q2HIGH,
     *                    CLAS5_RA_C,   CLAS5_RA_Cerr,
     *                    CLAS5_RA_FE,  CLAS5_RA_FEerr,
     *                    CLAS5_RA_PB,  CLAS5_RA_PBerr,
     *                    CLAS5_C_FUU,  CLAS5_C_FUUA,
     *                    CLAS5_FE_FUU, CLAS5_FE_FUUA,
     *                    CLAS5_PB_FUU, CLAS5_PB_FUUA,
     *                    CLAS5_C_DIS,  CLAS5_FE_DIS,
     *                    CLAS5_PB_DIS

      INTEGER I_CLAS5

      REAL*8 ::           CLAS6_PTLOW(22),  CLAS6_PTHIGH(22),
     *                    CLAS6_NULOW(22),  CLAS6_NUHIGH(22),
     *                    CLAS6_Q2LOW(22),  CLAS6_Q2HIGH(22),
     *                    CLAS6_RA_C(22),   CLAS6_RA_Cerr(22),
     *                    CLAS6_RA_FE(22),  CLAS6_RA_FEerr(22),
     *                    CLAS6_RA_PB(22),  CLAS6_RA_PBerr(22),
     *                    CLAS6_C_FUU(22),  CLAS6_C_FUUA(22),
     *                    CLAS6_FE_FUU(22), CLAS6_FE_FUUA(22),
     *                    CLAS6_PB_FUU(22), CLAS6_PB_FUUA(22),
     *                    CLAS6_C_DIS(22),  CLAS6_FE_DIS(22),
     *                    CLAS6_PB_DIS(22)

      COMMON /CLAS6/      CLAS6_PTLOW,  CLAS6_PTHIGH,
     *                    CLAS6_NULOW,  CLAS6_NUHIGH,
     *                    CLAS6_Q2LOW,  CLAS6_Q2HIGH,
     *                    CLAS6_RA_C,   CLAS6_RA_Cerr,
     *                    CLAS6_RA_FE,  CLAS6_RA_FEerr,
     *                    CLAS6_RA_PB,  CLAS6_RA_PBerr,
     *                    CLAS6_C_FUU,  CLAS6_C_FUUA,
     *                    CLAS6_FE_FUU, CLAS6_FE_FUUA,
     *                    CLAS6_PB_FUU, CLAS6_PB_FUUA,
     *                    CLAS6_C_DIS,  CLAS6_FE_DIS,
     *                    CLAS6_PB_DIS

      INTEGER I_CLAS6

      REAL*8 ::           CLAS7_PTLOW(22),  CLAS7_PTHIGH(22),
     *                    CLAS7_NULOW(22),  CLAS7_NUHIGH(22),
     *                    CLAS7_Q2LOW(22),  CLAS7_Q2HIGH(22),
     *                    CLAS7_RA_C(22),   CLAS7_RA_Cerr(22),
     *                    CLAS7_RA_FE(22),  CLAS7_RA_FEerr(22),
     *                    CLAS7_RA_PB(22),  CLAS7_RA_PBerr(22),
     *                    CLAS7_C_FUU(22),  CLAS7_C_FUUA(22),
     *                    CLAS7_FE_FUU(22), CLAS7_FE_FUUA(22),
     *                    CLAS7_PB_FUU(22), CLAS7_PB_FUUA(22),
     *                    CLAS7_C_DIS(22),  CLAS7_FE_DIS(22),
     *                    CLAS7_PB_DIS(22)


      COMMON /CLAS7/      CLAS7_PTLOW,  CLAS7_PTHIGH,
     *                    CLAS7_NULOW,  CLAS7_NUHIGH,
     *                    CLAS7_Q2LOW,  CLAS7_Q2HIGH,
     *                    CLAS7_RA_C,   CLAS7_RA_Cerr,
     *                    CLAS7_RA_FE,  CLAS7_RA_FEerr,
     *                    CLAS7_RA_PB,  CLAS7_RA_PBerr,
     *                    CLAS7_C_FUU,  CLAS7_C_FUUA,
     *                    CLAS7_FE_FUU, CLAS7_FE_FUUA,
     *                    CLAS7_PB_FUU, CLAS7_PB_FUUA,
     *                    CLAS7_C_DIS,  CLAS7_FE_DIS,
     *                    CLAS7_PB_DIS

      INTEGER I_CLAS7

      REAL*8 ::           CLAS8_PTLOW(22),  CLAS8_PTHIGH(22),
     *                    CLAS8_NULOW(22),  CLAS8_NUHIGH(22),
     *                    CLAS8_Q2LOW(22),  CLAS8_Q2HIGH(22),
     *                    CLAS8_RA_C(22),   CLAS8_RA_Cerr(22),
     *                    CLAS8_RA_FE(22),  CLAS8_RA_FEerr(22),
     *                    CLAS8_RA_PB(22),  CLAS8_RA_PBerr(22),
     *                    CLAS8_C_FUU(22),  CLAS8_C_FUUA(22),
     *                    CLAS8_FE_FUU(22), CLAS8_FE_FUUA(22),
     *                    CLAS8_PB_FUU(22), CLAS8_PB_FUUA(22),
     *                    CLAS8_C_DIS(22),  CLAS8_FE_DIS(22),
     *                    CLAS8_PB_DIS(22)


      COMMON /CLAS8/      CLAS8_PTLOW,  CLAS8_PTHIGH,
     *                    CLAS8_NULOW,  CLAS8_NUHIGH,
     *                    CLAS8_Q2LOW,  CLAS8_Q2HIGH,
     *                    CLAS8_RA_C,   CLAS8_RA_Cerr,
     *                    CLAS8_RA_FE,  CLAS8_RA_FEerr,
     *                    CLAS8_RA_PB,  CLAS8_RA_PBerr,
     *                    CLAS8_C_FUU,  CLAS8_C_FUUA,
     *                    CLAS8_FE_FUU, CLAS8_FE_FUUA,
     *                    CLAS8_PB_FUU, CLAS8_PB_FUUA,
     *                    CLAS8_C_DIS,  CLAS8_FE_DIS,
     *                    CLAS8_PB_DIS

      INTEGER I_CLAS8

      REAL*8 ::           CLAS9_PTLOW(22), CLAS9_PTHIGH(22),
     *                    CLAS9_NULOW(22), CLAS9_NUHIGH(22),
     *                    CLAS9_Q2LOW(22), CLAS9_Q2HIGH(22),
     *                    CLAS9_RA_C(22),  CLAS9_RA_Cerr(22),
     *                    CLAS9_RA_FE(22),  CLAS9_RA_FEerr(22),
     *                    CLAS9_RA_PB(22),  CLAS9_RA_PBerr(22),
     *                    CLAS9_C_FUU(22),  CLAS9_C_FUUA(22),
     *                    CLAS9_FE_FUU(22), CLAS9_FE_FUUA(22),
     *                    CLAS9_PB_FUU(22), CLAS9_PB_FUUA(22),
     *                    CLAS9_C_DIS(22),  CLAS9_FE_DIS(22),
     *                    CLAS9_PB_DIS(22)


      COMMON /CLAS9/      CLAS9_PTLOW,  CLAS9_PTHIGH,
     *                    CLAS9_NULOW,  CLAS9_NUHIGH,
     *                    CLAS9_Q2LOW,  CLAS9_Q2HIGH,
     *                    CLAS9_RA_C,   CLAS9_RA_Cerr,
     *                    CLAS9_RA_FE,  CLAS9_RA_FEerr,
     *                    CLAS9_RA_PB,  CLAS9_RA_PBerr,
     *                    CLAS9_C_FUU,  CLAS9_C_FUUA,
     *                    CLAS9_FE_FUU, CLAS9_FE_FUUA,
     *                    CLAS9_PB_FUU, CLAS9_PB_FUUA,
     *                    CLAS9_C_DIS,  CLAS9_FE_DIS,
     *                    CLAS9_PB_DIS

      INTEGER I_CLAS9

C--------JLAB2022 DATA
C--------pi+
      REAL*8 ::           pip2022_zlow(48),   pip2022_zup(48),
     *                    pip2022_bin(48),  
     *                    pip2022_ptlow(48),  pip2022_ptup(48),
     *                    pip2022_c(48), 
     *                    pip2022_cstat(48),  pip2022_csys(48),
     *                    pip2022_fe(48),
     *                    pip2022_festat(48), pip2022_fesys(48),
     *                    pip2022_pb(48),
     *                    pip2022_pbstat(48), pip2022_pbsys(48),
     *                    pip2022_qc(48),     pip2022_xc(48),
     *                    pip2022_qfe(48),    pip2022_xfe(48),
     *                    pip2022_qpb(48),    pip2022_xpb(48),
     *                    pip2022_cfuu(48),   pip2022_cfuua(48),
     *                    pip2022_cdis(48),
     *                    pip2022_fefuu(48),   pip2022_fefuua(48),
     *                    pip2022_fedis(48),
     *                    pip2022_pbfuu(48),   pip2022_pbfuua(48),
     *                    pip2022_pbdis(48)

      COMMON /pip2022/    pip2022_zlow,   pip2022_zup,
     *                    pip2022_bin,
     *                    pip2022_ptlow,  pip2022_ptup,
     *                    pip2022_c, 
     *                    pip2022_cstat,  pip2022_csys,
     *                    pip2022_fe,
     *                    pip2022_festat, pip2022_fesys,
     *                    pip2022_pb, 
     *                    pip2022_pbstat, pip2022_pbsys,
     *                    pip2022_qc,     pip2022_xc,
     *                    pip2022_qfe,    pip2022_xfe,
     *                    pip2022_qpb,    pip2022_xpb,
     *                    pip2022_cfuu,   pip2022_cfuua,
     *                    pip2022_cdis,
     *                    pip2022_fefuu,  pip2022_fefuua,
     *                    pip2022_fedis,
     *                    pip2022_pbfuu,  pip2022_pbfuua,
     *                    pip2022_pbdis

      INTEGER I_pip_2022

C--------pi-
      REAL*8 ::           pim2022_zlow(40),   pim2022_zup(40),
     *                    pim2022_bin(40),  
     *                    pim2022_ptlow(40),  pim2022_ptup(40),
     *                    pim2022_c(40), 
     *                    pim2022_cstat(40),  pim2022_csys(40),
     *                    pim2022_fe(40),
     *                    pim2022_festat(40), pim2022_fesys(40),
     *                    pim2022_pb(40),
     *                    pim2022_pbstat(40), pim2022_pbsys(40),
     *                    pim2022_qc(40),     pim2022_xc(40),
     *                    pim2022_qfe(40),    pim2022_xfe(40),
     *                    pim2022_qpb(40),    pim2022_xpb(40),
     *                    pim2022_cfuu(40),   pim2022_cfuua(40),
     *                    pim2022_cdis(40),
     *                    pim2022_fefuu(40),  pim2022_fefuua(40),
     *                    pim2022_fedis(40),
     *                    pim2022_pbfuu(40),  pim2022_pbfuua(40),
     *                    pim2022_pbdis(40)


      COMMON /pim2022/    pim2022_zlow,   pim2022_zup,
     *                    pim2022_bin,
     *                    pim2022_ptlow,  pim2022_ptup,
     *                    pim2022_c, 
     *                    pim2022_cstat,  pim2022_csys,
     *                    pim2022_fe,
     *                    pim2022_festat, pim2022_fesys,
     *                    pim2022_pb, 
     *                    pim2022_pbstat, pim2022_pbsys,
     *                    pim2022_qc,     pim2022_xc,
     *                    pim2022_qfe,    pim2022_xfe,
     *                    pim2022_qpb,    pim2022_xpb,
     *                    pim2022_cfuu,   pim2022_cfuua,
     *                    pim2022_cdis,
     *                    pim2022_fefuu,  pim2022_fefuua,
     *                    pim2022_fedis,
     *                    pim2022_pbfuu,  pim2022_pbfuua,
     *                    pim2022_pbdis

      INTEGER I_pim_2022 !2023
      INTEGER I_pip_2022_C 
      INTEGER I_pip_2022_Fe 
      INTEGER I_pip_2022_Pb
      INTEGER I_pim_2022_C  
      INTEGER I_pim_2022_Fe
      INTEGER I_pim_2022_Pb

      INTEGER I_EIC
    
C- - - - - - - - - - - - - - - - -
C----------------------------------
C------ ALL DATA
      COMMON /DATASET/ I_PIP_HE_NU,  I_PIP_NE_NU,  I_PIP_KR_NU,
     $ I_PIP_XE_NU,   I_PIP_HE_Z,    I_PIP_NE_Z,   I_PIP_KR_Z,
     $ I_PIP_XE_Z,    I_PIP_HE_Q2,   I_PIP_NE_Q2,  I_PIP_KR_Q2,
     $ I_PIP_XE_Q2,   I_KP_HE_NU,    I_KP_NE_NU,   I_KP_KR_NU,
     $ I_KP_XE_NU,    I_KP_HE_Z,     I_KP_NE_Z,    I_KP_KR_Z,
     $ I_KP_XE_Z,     I_KP_HE_Q2,    I_KP_NE_Q2,   I_KP_KR_Q2,
     $ I_KP_XE_Q2,    I_PIM_HE_NU,   I_PIM_NE_NU,  I_PIM_KR_NU,
     $ I_PIM_XE_NU,   I_PIM_HE_Z,    I_PIM_NE_Z,   I_PIM_KR_Z,
     $ I_PIM_XE_Z,    I_PIM_HE_Q2,   I_PIM_NE_Q2,  I_PIM_KR_Q2,
     $ I_PIM_XE_Q2,   I_KM_HE_NU,    I_KM_NE_NU,   I_KM_KR_NU,
     $ I_KM_XE_NU,    I_KM_HE_Z,     I_KM_NE_Z,    I_KM_KR_Z,
     $ I_KM_XE_Z,     I_KM_HE_Q2,    I_KM_NE_Q2,   I_KM_KR_Q2,
     $ I_KM_XE_Q2,    I_PIP_HE_PT2,  I_PIP_NE_PT2, I_PIP_KR_PT2,
     $ I_PIP_XE_PT2,  I_KP_HE_PT2,   I_KP_NE_PT2,  I_KP_KR_PT2,
     $ I_KP_XE_PT2,   I_PIM_HE_PT2,  I_PIM_NE_PT2, I_PIM_KR_PT2,
     $ I_PIM_XE_PT2,  I_KM_HE_PT2,   I_KM_NE_PT2,  I_KM_KR_PT2,
     $ I_KM_XE_PT2,   I_PI0_HE_NU,   I_PI0_NE_NU,  I_PI0_KR_NU,
     $ I_PI0_XE_NU,   I_PI0_HE_Z,    I_PI0_NE_Z,   I_PI0_KR_Z,
     $ I_PI0_XE_Z,    I_PI0_HE_Q2,   I_PI0_NE_Q2,  I_PI0_KR_Q2,
     $ I_PI0_XE_Q2,   I_PI0_HE_PT2,  I_PI0_NE_PT2, I_PI0_KR_PT2,
     $ I_PI0_XE_PT2,   I_E288_200,
     $                 I_E288_300,
     $                 I_E288_400,
     $                 I_E866_800,
     $                 I_E866_800q,
     $                 I_E772_800,
     $                 I_E605_800,
     $ I_CLAS1, I_CLAS2, I_CLAS3, I_CLAS4, I_CLAS5, I_CLAS6, I_CLAS7,
     $ I_CLAS8, I_CLAS9,
     $ I_RHIC_pp, I_RHIC_pAu1, I_RHIC_pAu2,
     $ I_RHIC_Ratio_pAu1, I_RHIC_Ratio_pAu2,
     $ I_ATLAS5_Y1,I_ATLAS5_Y2,I_ATLAS5_Y3,I_ATLAS5_Y_RAT,
     $ I_CMS5,
     $ I_ATLAS3,
     $ I_CMS8_DY_no,I_CMS8_DY_fid,I_CMS8_ZZ_no,I_CMS8_ZZ_fid,
     $ I_JLAB_PIP_12, I_JLAB_PIM_12,
     $ I_EIC_PRED, I_EICC_PRED, I_JLAB_PRED, 
     $ I_x_PRED, I_z_PRED, I_pht_PRED, 
     $ I_pip_2022, I_pim_2022, !2023
     $ I_pip_2022_C,
     $ I_pip_2022_Fe,
     $ I_pip_2022_Pb,
     $ I_pim_2022_C,
     $ I_pim_2022_Fe,
     $ I_pim_2022_Pb,
     $ I_EIC
