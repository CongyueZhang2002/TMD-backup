      PROGRAM MASTER
      IMPLICIT NONE

C-----VARIABLE DECLARATIONS
      REAL*8 VAL, MULT, STAT, SYS
      INTEGER I, J
      INCLUDE './tools/data-inc-DIS.f'
C------ VARIABLES FOR APFEL SUBROUTINES
      double precision Q0,Q
      double precision Nu, Z, PT2, xb, y, Q02, Q2
      double precision M, alpha_EM, S_EP
      double precision F2light, F2light_d, F2light_a
      double precision FLlight, FLlight_d, FLlight_a
      real*8 disd, disa, fuu, fuua, sidis_d, sidis_a
      real*8 nucl_target
      real*8 err
      real*8 mult_a, mult_d, R_a,dis
      real*8 xfmin, xfmax, ymin, ymax, Qbar, const, ybar
      real*8 pt2min, pt2max, numin, numax, Q2min, Q2max, zmin, zmax
      real*8 ptmin, ptmax
      real*8 R_C, Err_C, R_Fe, Err_Fe, R_Pb, Err_Pb
      real*8 fuu_D, fuu_C, fuu_Fe, fuu_Pb
      real*8 DIS_C, DIS_Fe, DIS_Pb
      double precision eps
      Character*72 HE, NE, XE, KR
      integer nfit,nloops,hop,nll
      integer N,fft,chi2_pof
      integer nloop,order,hopscheme
      integer usenc
      INTEGER IINPDF
      real*8 Qs ,zh ,pht !2023
      COMMON /nckern/ usenc
      COMMON /params/ nfit
      COMMON /scheme/ nloops,hop,nll
      COMMON /FT/ fft
      COMMON / INTINPDF / IINPDF
      parameter(eps=1d-10)

C---------CHOOSE FITTING SCHEME
      usenc = 0
      hop = 0
      nloops = 1
      nll=2
      fft = 1
      IINPDF = 0

C----------To not show the calls, although they are called
      CALL SetLHAPARM('SILENT')

C------ SET INITIAL SCALE Q0 --> Based on PDF Set (nCTEQ15)
      Q0 = 1.3d0 ! GeV
C----- SET PARAMETERS
      M = 0.938272046d0 ! GeV (Mass of Proton)
      alpha_EM = (1/137)
      S_EP = 52.7 ! GeV^2

C------NUCLEI
      HE = 'He'
      NE = 'Ne'
      XE = 'Xe'
      KR = 'Kr'

C-----CHOOSE EXP DATA SETS
      I_PIP_HE_Z   = 0
      I_PIP_NE_Z   = 0
      I_PIP_KR_Z   = 0
      I_PIP_XE_Z   = 0
      I_PIM_HE_Z   = 0
      I_PIM_NE_Z   = 0
      I_PIM_KR_Z   = 0
      I_PIM_XE_Z   = 0
      I_PI0_HE_Z   = 0
      I_PI0_NE_Z   = 0
      I_PI0_KR_Z   = 0
      I_PI0_XE_Z   = 0

      I_PI0_HE_PT2   = 0
      I_PI0_NE_PT2   = 0
      I_PI0_KR_PT2   = 0
      I_PI0_XE_PT2   = 0
      I_PIP_NE_PT2   = 0
      I_PIP_KR_PT2   = 0
      I_PIP_XE_PT2   = 0
      I_PIM_HE_PT2   = 0
      I_PIM_NE_PT2   = 0
      I_PIM_KR_PT2   = 0
      I_PIM_XE_PT2   = 0

      I_JLAB_PIP_12     = 0
      I_JLAB_PIM_12     = 0

      I_JLAB_PRED       = 1
      I_EIC_PRED        = 0
      I_EICC_PRED        = 0

      I_pip_2022 = 0
      I_pim_2022 = 0
   
      I_x_PRED = 0
      I_z_PRED = 0
      I_pht_PRED = 0

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

C------ READ THE DATA
      CALL READDATA

C-----triple
      IF(I_x_PRED.EQ.1) THEN
      OPEN(UNIT = 1,
     *STATUS='UNKNOWN', FILE ='expdata/PRED/DIS_x.dat')
      WRITE(1,*) 'Qs  ', 'xb ' , 'zh ', 'pht ', 'DIS '
      DO I=1,270 !2023
          Qs = x_PRED_Qs(I)  
          xb = x_PRED_xb(I) 
          zh = x_PRED_zh(I)  
          pht = x_PRED_pht(I)
          DIS = x_PRED_dis(I) 
          WRITE(1,*) Qs ,xb ,zh ,pht , DIS 
      ENDDO
      CLOSE(1)
      ENDIF

      IF(I_z_PRED.EQ.1) THEN
      OPEN(UNIT = 2,
     *STATUS='UNKNOWN',FILE ='expdata/PRED/DIS_z.dat')
      WRITE(2,*) 'Qs  ', 'xb ' , 'zh ', 'pht ', 'DIS '
      DO I=1,270
          Qs = z_PRED_Qs(I)  
          xb = z_PRED_xb(I) 
          zh = z_PRED_zh(I)  
          pht = z_PRED_pht(I)
          DIS = z_PRED_dis(I) 
          WRITE(2,*) Qs ,xb ,zh ,pht , DIS 
      ENDDO
      CLOSE(2)
      ENDIF

      IF(I_pht_PRED.EQ.1) THEN
      OPEN(UNIT = 3,
     *STATUS='UNKNOWN',FILE ='expdata/PRED/DIS_pht.dat')
      WRITE(3,*) 'Qs  ', 'xb ' , 'zh ', 'pht ', 'DIS '
      DO I=1,270
          Qs = pht_PRED_Qs(I)  
          xb = pht_PRED_xb(I) 
          zh = pht_PRED_zh(I)  
          pht = pht_PRED_pht(I)
          DIS = pht_PRED_dis(I) 
          WRITE(3,*) Qs ,xb ,zh ,pht , DIS 
      ENDDO
      CLOSE(3)
      ENDIF

C-----SIDIS (JLAB 12 preliminary data from Miguel)
      IF(I_JLAB_PIP_12.EQ.1) THEN
      OPEN(UNIT = 5,FILE ='expdata/JLAB_piplus.dat')
      WRITE(5,*) 'target  ', 'Z ' , 'PTMIN ', 'PTMAX ', 'MULT-RATIO ',
     *       'ERR ', 'fuu ', 'fuua ',  'DIS'
      DO I=1,180
          nucl_target = JLAB_PIP_TARGET(I)
          z           = JLAB_PIP_Z(I)/10d0
          ptmin       = JLAB_PIP_PTLOW(I)
          ptmax       = JLAB_PIP_PTHIGH(I)
          MULT        = JLAB_PIP_MULT_RATIO(I)
          ERR         = JLAB_PIP_ERR(I)
          fuu         = JLAB_PIP_FUU(I)
          fuua        = JLAB_PIP_FUUA(I)
          DIS         = JLAB_PIP_DIS(I)
          WRITE(5,*) nucl_target, z, ptmin, ptmax, MULT, ERR,
     *               fuu, fuua, DIS
      ENDDO
      CLOSE(5)
      ENDIF

C-----SIDIS (JLAB 12 preliminary data from Miguel)
      IF(I_JLAB_PIM_12.EQ.1) THEN
      OPEN(UNIT = 5,FILE ='expdata/JLAB_piminus.dat')
      WRITE(5,*) 'target  ', 'Z ' , 'PTMIN ', 'PTMAX ', 'MULT-RATIO ',
     *       'ERR ', 'fuu ', 'fuua ',  'DIS'
      DO I=1,180
          nucl_target = JLAB_PIM_TARGET(I)
          z           = JLAB_PIM_Z(I)/10d0
          ptmin       = JLAB_PIM_PTLOW(I)
          ptmax       = JLAB_PIM_PTHIGH(I)
          MULT        = JLAB_PIM_MULT_RATIO(I)
          ERR         = JLAB_PIM_ERR(I)
          fuu         = JLAB_PIM_FUU(I)
          fuua        = JLAB_PIM_FUUA(I)
          DIS         = JLAB_PIM_DIS(I)
          WRITE(5,*) nucl_target, z, ptmin, ptmax, MULT, ERR,
     *               fuu, fuua, DIS
      ENDDO
      CLOSE(5)
      ENDIF


C-----SIDIS (HERMES)
      IF(I_PIP_HE_Z.EQ.1) THEN
      OPEN(UNIT = 5,FILE ='expdata/SIDIS/HERMES_DIS/pip_he_z.dat')
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
          fuu = PIP_HE_Z_FUU(I)
          fuua = PIP_HE_Z_FUUA(I)
          DIS = PIP_HE_Z_DIS(I)
          R_a = fuua/DIS
          !CALL Ra(He, Q0, Q, xb, Nu, fuu, fuua, mult_d, mult_a, R_a)
          WRITE(5,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
      ENDDO
      CLOSE(5)
      ENDIF

      IF(I_PIP_NE_Z.EQ.1) THEN
      OPEN(UNIT = 66,FILE ='expdata/SIDIS/HERMES_DIS/pip_ne_z.dat')
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
          fuu = PIP_NE_Z_FUU(I)
          fuua = PIP_NE_Z_FUUA(I)
          DIS = PIP_NE_Z_DIS(I)
          R_a = fuua/DIS
          !CALL Ra(Ne, Q0, Q, xb, Nu, fuu, fuua, mult_d, mult_a, R_a)
          WRITE(66,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
      ENDDO
      CLOSE(66)
      ENDIF

      IF(I_PIP_KR_Z.EQ.1) THEN
      OPEN(UNIT = 7,FILE ='expdata/SIDIS/HERMES_DIS/pip_kr_z.dat')
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
          fuu = PIP_KR_Z_FUU(I)
          fuua = PIP_KR_Z_FUUA(I)
          DIS = PIP_KR_Z_DIS(I)
          R_a = fuua/DIS
          !CALL Ra(Kr, Q0, Q, xb, Nu, fuu, fuua, mult_d, mult_a, R_a)
          WRITE(7,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
      ENDDO
      CLOSE(7)
      ENDIF

      IF(I_PIP_XE_Z.EQ.1) THEN
      OPEN(UNIT = 8,FILE ='expdata/SIDIS/HERMES_DIS/pip_xe_z.dat')
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
          fuu  = PIP_XE_Z_FUU(I)
          fuua = PIP_XE_Z_FUUA(I)
          DIS  = PIP_XE_Z_DIS(I)
          R_a  = fuua/DIS
          !CALL Ra(Xe, Q0, Q, xb, Nu, fuu, fuua, mult_d, mult_a, R_a)
          WRITE(8,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
      ENDDO
      CLOSE(8)
      ENDIF

C------ PI -

      IF(I_PIM_HE_Z.EQ.1) THEN
      OPEN(UNIT = 41,FILE ='expdata/SIDIS/HERMES_DIS/pim_he_z.dat')
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
          fuu = PIM_HE_Z_FUU(I)
          fuua = PIM_HE_Z_FUUA(I)
          DIS = PIM_HE_Z_DIS(I)
          R_a = fuua/DIS
          !CALL Ra(He, Q0, Q, xb, Nu, fuu, fuua, mult_d, mult_a, R_a)
          WRITE(41,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
      ENDDO
      CLOSE(41)
      ENDIF

      IF(I_PIM_NE_Z.EQ.1) THEN
      OPEN(UNIT = 42,FILE ='expdata/SIDIS/HERMES_DIS/pim_ne_z.dat')
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
          fuu = PIM_NE_Z_FUU(I)
          fuua = PIM_NE_Z_FUUA(I)
          DIS = PIM_NE_Z_DIS(I)
          R_a = fuua/DIS
          !CALL Ra(Ne, Q0, Q, xb, Nu, fuu, fuua, mult_d, mult_a, R_a)
          WRITE(42,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
      ENDDO
      CLOSE(42)
      ENDIF

      IF(I_PIM_KR_Z.EQ.1) THEN
      OPEN(UNIT = 43,FILE ='expdata/SIDIS/HERMES_DIS/pim_kr_z.dat')
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
          fuu = PIM_KR_Z_FUU(I)
          fuua = PIM_KR_Z_FUUA(I)
          DIS = PIM_KR_Z_DIS(I)
          R_a = fuua/DIS
          !CALL Ra(Kr, Q0, Q, xb, Nu, fuu, fuua, mult_d, mult_a, R_a)
          WRITE(43,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
      ENDDO
      CLOSE(43)
      ENDIF

      IF(I_PIM_XE_Z.EQ.1) THEN
      OPEN(UNIT = 44,FILE ='expdata/SIDIS/HERMES_DIS/pim_xe_z.dat')
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
          fuu  = PIM_XE_Z_FUU(I)
          fuua = PIM_XE_Z_FUUA(I)
          DIS  = PIM_XE_Z_DIS(I)
          R_a  = fuua/DIS
          !CALL Ra(Xe, Q0, Q, xb, Nu, fuu, fuua, mult_d, mult_a, R_a)
          WRITE(44,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
      ENDDO
      CLOSE(44)
      ENDIF

      IF(I_PI0_HE_Z.EQ.1) THEN
      OPEN(UNIT = 105,FILE ='expdata/SIDIS/HERMES_DIS/pi0_he_z.dat')
      WRITE(105,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
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
          fuu = PI0_HE_Z_FUU(I)
          fuua = PI0_HE_Z_FUUA(I)
          DIS = PI0_HE_Z_DIS(I)
          R_a = fuua/DIS
          !CALL Ra(He, Q0, Q, xb, Nu, fuu, fuua, mult_d, mult_a, R_a)
          WRITE(105,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
      ENDDO
      CLOSE(105)
      ENDIF

      IF(I_PI0_NE_Z.EQ.1) THEN
      OPEN(UNIT = 106,FILE ='expdata/SIDIS/HERMES_DIS/pi0_ne_z.dat')
      WRITE(106,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
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
          fuu = PI0_NE_Z_FUU(I)
          fuua = PI0_NE_Z_FUUA(I)
          DIS = PI0_NE_Z_DIS(I)
          R_a = fuua/DIS
          !CALL Ra(Ne, Q0, Q, xb, Nu, fuu, fuua, mult_d, mult_a, R_a)
          WRITE(106,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
      ENDDO
      CLOSE(106)
      ENDIF

      IF(I_PI0_KR_Z.EQ.1) THEN
      OPEN(UNIT = 107,FILE ='expdata/SIDIS/HERMES_DIS/pi0_kr_z.dat')
      WRITE(107,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
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
          fuu  = PI0_KR_Z_FUU(I)
          fuua = PI0_KR_Z_FUUA(I)
          DIS  = PI0_KR_Z_DIS(I)
          R_a  = fuua/DIS
          !CALL Ra(Kr, Q0, Q, xb, Nu, fuu, fuua, mult_d, mult_a, R_a)
          WRITE(107,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
      ENDDO
      CLOSE(107)
      ENDIF

      IF(I_PI0_XE_Z.EQ.1) THEN
      OPEN(UNIT = 108,FILE ='expdata/SIDIS/HERMES_DIS/pi0_xe_z.dat')
      WRITE(108,*) 'Nu  ', 'Z ' , 'Q2 ', 'pt2 ', 'MULT-RATIO ',
     *       'STAT ', 'SYS ', 'fuu ', 'fuua ', 'DIS'
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
          fuu  = PI0_XE_Z_FUU(I)
          fuua = PI0_XE_Z_FUUA(I)
          DIS  = PI0_XE_Z_DIS(I)
          R_a  = fuua/DIS
          !CALL Ra(Xe, Q0, Q, xb, Nu, fuu, fuua, mult_d, mult_a, R_a)
          WRITE(108,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
      ENDDO
      CLOSE(108)
      ENDIF

C------- PI0 DATA
      IF(I_PI0_HE_PT2.EQ.1) THEN
      OPEN(UNIT = 203,FILE ='expdata/SIDIS/HERMES_DIS/pi0_he_pt2.dat')
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
          fuu = PI0_HE_PT2_FUU(I)
          fuua = PI0_HE_PT2_FUUA(I)
          DIS = PI0_HE_PT2_DIS(I)
          R_a = fuua/DIS
          !CALL Ra(He, Q0, Q, xb, Nu, fuu, fuua, mult_d, mult_a, R_a)
          WRITE(203,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
      ENDDO
      CLOSE(203)
      ENDIF

      IF(I_PI0_NE_PT2.EQ.1) THEN
      OPEN(UNIT = 204,FILE ='expdata/SIDIS/HERMES_DIS/pi0_ne_pt2.dat')
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
          fuu = PI0_NE_PT2_FUU(I)
          fuua = PI0_NE_PT2_FUUA(I)
          DIS = PI0_NE_PT2_DIS(I)
          R_a = fuua/DIS
          !CALL Ra(Ne, Q0, Q, xb, Nu, fuu, fuua, mult_d, mult_a, R_a)
          WRITE(204,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
      ENDDO
      CLOSE(204)
      ENDIF

      IF(I_PI0_KR_PT2.EQ.1) THEN
      OPEN(UNIT = 205,FILE ='expdata/SIDIS/HERMES_DIS/pi0_kr_pt2.dat')
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
          fuu = PI0_KR_PT2_FUU(I)
          fuua = PI0_KR_PT2_FUUA(I)
          DIS = PI0_KR_PT2_DIS(I)
          R_a = fuua/DIS
          !CALL Ra(Kr, Q0, Q, xb, Nu, fuu, fuua, mult_d, mult_a, R_a)
          WRITE(205,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
      ENDDO
      CLOSE(205)
      ENDIF

      IF(I_PI0_XE_PT2.EQ.1) THEN
      OPEN(UNIT = 206,FILE ='expdata/SIDIS/HERMES_DIS/pi0_xe_pt2.dat')
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
          fuu = PI0_XE_PT2_FUU(I)
          fuua = PI0_XE_PT2_FUUA(I)
          DIS = PI0_XE_PT2_DIS(I)
          R_a = fuua/DIS
          !CALL Ra(Xe, Q0, Q, xb, Nu, fuu, fuua, mult_d, mult_a, R_a)
          WRITE(206,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
      ENDDO
      CLOSE(206)
      ENDIF

C------------------
C------ PT2 DEPENDET DATA
C------- PI + DATA
      IF(I_PIP_HE_PT2.EQ.1) THEN
      OPEN(UNIT = 73,FILE ='expdata/SIDIS/HERMES_DIS/pip_he_pt2.dat')
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
          fuu = PIP_HE_PT2_FUU(I)
          fuua = PIP_HE_PT2_FUUA(I)
          DIS = PIP_HE_PT2_DIS(I)
          R_a = fuua/DIS
          !CALL Ra(He, Q0, Q, xb, Nu, fuu, fuua, mult_d, mult_a, R_a)
          WRITE(73,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
      ENDDO
      CLOSE(73)
      ENDIF

      IF(I_PIP_NE_PT2.EQ.1) THEN
      OPEN(UNIT = 74,FILE ='expdata/SIDIS/HERMES_DIS/pip_ne_pt2.dat')
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
          fuu = PIP_NE_PT2_FUU(I)
          fuua = PIP_NE_PT2_FUUA(I)
          DIS = PIP_NE_PT2_DIS(I)
          R_a = fuua/DIS
          !CALL Ra(Ne, Q0, Q, xb, Nu, fuu, fuua, mult_d, mult_a, R_a)
          WRITE(74,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
      ENDDO
      CLOSE(74)
      ENDIF

      IF(I_PIP_KR_PT2.EQ.1) THEN
      OPEN(UNIT = 75,FILE ='expdata/SIDIS/HERMES_DIS/pip_kr_pt2.dat')
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
          fuu = PIP_KR_PT2_FUU(I)
          fuua = PIP_KR_PT2_FUUA(I)
          DIS = PIP_KR_PT2_DIS(I)
          R_a = fuua/DIS
          !CALL Ra(Kr, Q0, Q, xb, Nu, fuu, fuua, mult_d, mult_a, R_a)
          WRITE(75,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
      ENDDO
      CLOSE(75)
      ENDIF

      IF(I_PIP_XE_PT2.EQ.1) THEN
      OPEN(UNIT = 76,FILE ='expdata/SIDIS/HERMES_DIS/pip_xe_pt2.dat')
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
          fuu = PIP_XE_PT2_FUU(I)
          fuua = PIP_XE_PT2_FUUA(I)
          DIS = PIP_XE_PT2_DIS(I)
          R_a = fuua/DIS
          !CALL Ra(Xe, Q0, Q, xb, Nu, fuu, fuua, mult_d, mult_a, R_a)
          WRITE(76,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
      ENDDO
      CLOSE(76)
      ENDIF

C------ PI - DATA
      IF(I_PIM_HE_PT2.EQ.1) THEN
      OPEN(UNIT = 85,FILE ='expdata/SIDIS/HERMES_DIS/pim_he_pt2.dat')
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
          fuu = PIM_HE_PT2_FUU(I)
          fuua = PIM_HE_PT2_FUUA(I)
          DIS = PIM_HE_PT2_DIS(I)
          R_a = fuua/DIS
          !CALL Ra(He, Q0, Q, xb, Nu, fuu, fuua, mult_d, mult_a, R_a)
          WRITE(85,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
      ENDDO
      CLOSE(85)
      ENDIF

      IF(I_PIM_NE_PT2.EQ.1) THEN
      OPEN(UNIT = 86,FILE ='expdata/SIDIS/HERMES_DIS/pim_ne_pt2.dat')
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
          fuu = PIM_NE_PT2_FUU(I)
          fuua = PIM_NE_PT2_FUUA(I)
          DIS = PIM_NE_PT2_DIS(I)
          R_a = fuua/DIS
          !CALL Ra(Ne, Q0, Q, xb, Nu, fuu, fuua, mult_d, mult_a, R_a)
          WRITE(86,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
      ENDDO
      CLOSE(86)
      ENDIF

      IF(I_PIM_KR_PT2.EQ.1) THEN
      OPEN(UNIT = 87,FILE ='expdata/SIDIS/HERMES_DIS/pim_kr_pt2.dat')
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
          fuu = PIM_KR_PT2_FUU(I)
          fuua = PIM_KR_PT2_FUUA(I)
          DIS = PIM_KR_PT2_DIS(I)
          R_a = fuua/DIS
          !CALL Ra(Kr, Q0, Q, xb, Nu, fuu, fuua, mult_d, mult_a, R_a)
          WRITE(87,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
      ENDDO
      CLOSE(87)
      ENDIF

      IF(I_PIM_XE_PT2.EQ.1) THEN
      OPEN(UNIT = 88,FILE ='expdata/SIDIS/HERMES_DIS/pim_xe_pt2.dat')
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
          fuu = PIM_XE_PT2_FUU(I)
          fuua = PIM_XE_PT2_FUUA(I)
          DIS = PIM_XE_PT2_DIS(I)
          R_a = fuua/DIS
          !CALL Ra(Xe, Q0, Q, xb, Nu, fuu, fuua, mult_d, mult_a, R_a)
          WRITE(88,*) Nu, Z, Q2, pt2, MULT,STAT,SYS,fuu, fuua, DIS
      ENDDO
      CLOSE(88)
      ENDIF

C-----JLAB2022
C------pi+
C-------C
      IF(I_pip_2022.EQ.1) THEN
      OPEN(UNIT = 89,FILE ='expdata/JLAB2022/pi+C.dat')
      WRITE(89,*) 'Z ' , 'PTMIN ', 'PTMAX ', 'MULT-RATIO ',
     *       'ERR ', 'fuu ', 'fuua ',  'DIS'
      DO I=1,48
          z           = 0.5d0*(pip2022_zlow(I)+pip2022_zup(I))
          ptmin       = pip2022_ptlow(I)
          ptmax       = pip2022_ptup(I)
          MULT        = pip2022_c(I)
          ERR         = sqrt(pip2022_cstat(I)**2+pip2022_csys(I)**2)
          fuu         = pip2022_cfuu(I)
          fuua        = pip2022_cfuua(I)
          DIS         = pip2022_cdis(I)
          WRITE(89,*) z, ptmin, ptmax, MULT, ERR,
     *               fuu, fuua, DIS
      print *, "pip c finished"
      ENDDO

      CLOSE(89)
      ENDIF
C-------Fe
      IF(I_pip_2022.EQ.1) THEN
      OPEN(UNIT = 90,FILE ='expdata/JLAB2022/pi+Fe.dat')
      WRITE(90,*) 'Z ' , 'PTMIN ', 'PTMAX ', 'MULT-RATIO ',
     *       'ERR ', 'fuu ', 'fuua ',  'DIS'
      DO I=1,48
          z           = 0.5d0*(pip2022_zlow(I)+pip2022_zup(I))
          ptmin       = pip2022_ptlow(I)
          ptmax       = pip2022_ptup(I)
          MULT        = pip2022_c(I)
          ERR         = sqrt(pip2022_festat(I)**2+pip2022_fesys(I)**2)
          fuu         = pip2022_fefuu(I)
          fuua        = pip2022_fefuua(I)
          DIS         = pip2022_fedis(I)
          WRITE(90,*) z, ptmin, ptmax, MULT, ERR,
     *               fuu, fuua, DIS
      print *, "pip fe finished"
      ENDDO
      CLOSE(90)
      ENDIF
C-------Pb
      IF(I_pip_2022.EQ.1) THEN
      OPEN(UNIT = 91,FILE ='expdata/JLAB2022/pi+Pb.dat')
      WRITE(91,*) 'Z ' , 'PTMIN ', 'PTMAX ', 'MULT-RATIO ',
     *       'ERR ', 'fuu ', 'fuua ',  'DIS'
      DO I=1,48
          z           = 0.5d0*(pip2022_zlow(I)+pip2022_zup(I))
          ptmin       = pip2022_ptlow(I)
          ptmax       = pip2022_ptup(I)
          MULT        = pip2022_pb(I)
          ERR         = sqrt(pip2022_pbstat(I)**2+pip2022_pbsys(I)**2)
          fuu         = pip2022_pbfuu(I)
          fuua        = pip2022_pbfuua(I)
          DIS         = pip2022_pbdis(I)
          WRITE(91,*) z, ptmin, ptmax, MULT, ERR,
     *               fuu, fuua, DIS
      print *, "pip pb finished"
      ENDDO
      CLOSE(91)
      ENDIF

C------pi-
C-------C
      IF(I_pim_2022.EQ.1) THEN
      OPEN(UNIT = 92,FILE ='expdata/JLAB2022/pi-C.dat')
      WRITE(92,*) 'Z ' , 'PTMIN ', 'PTMAX ', 'MULT-RATIO ',
     *       'ERR ', 'fuu ', 'fuua ',  'DIS'
      DO I=1,40
          z           = 0.5d0*(pim2022_zlow(I)+pim2022_zup(I))
          ptmin       = pim2022_ptlow(I)
          ptmax       = pim2022_ptup(I)
          MULT        = pim2022_c(I)
          ERR         = sqrt(pim2022_cstat(I)**2+pim2022_csys(I)**2)
          fuu         = pim2022_cfuu(I)
          fuua        = pim2022_cfuua(I)
          DIS         = pim2022_cdis(I)
          WRITE(92,*) z, ptmin, ptmax, MULT, ERR,
     *               fuu, fuua, DIS
      print *, "pim c finished"
      ENDDO
      CLOSE(92)
      ENDIF
C-------Fe
      IF(I_pim_2022.EQ.1) THEN
      OPEN(UNIT = 93,FILE ='expdata/JLAB2022/pi-Fe.dat')
      WRITE(93,*) 'Z ' , 'PTMIN ', 'PTMAX ', 'MULT-RATIO ',
     *       'ERR ', 'fuu ', 'fuua ',  'DIS'
      DO I=1,40
          z           = 0.5d0*(pim2022_zlow(I)+pim2022_zup(I))
          ptmin       = pim2022_ptlow(I)
          ptmax       = pim2022_ptup(I)
          MULT        = pim2022_c(I)
          ERR         = sqrt(pim2022_festat(I)**2+pim2022_fesys(I)**2)
          fuu         = pim2022_fefuu(I)
          fuua        = pim2022_fefuua(I)
          DIS         = pim2022_fedis(I)
          WRITE(93,*) z, ptmin, ptmax, MULT, ERR,
     *               fuu, fuua, DIS
      print *, "pim fe finished"
      ENDDO
      CLOSE(93)
      ENDIF
C-------Pb
      IF(I_pim_2022.EQ.1) THEN
      OPEN(UNIT = 94,FILE ='expdata/JLAB2022/pi-Pb.dat')
      WRITE(94,*) 'Z ' , 'PTMIN ', 'PTMAX ', 'MULT-RATIO ',
     *       'ERR ', 'fuu ', 'fuua ',  'DIS'
      DO I=1,40
          z           = 0.5d0*(pim2022_zlow(I)+pim2022_zup(I))
          ptmin       = pim2022_ptlow(I)
          ptmax       = pim2022_ptup(I)
          MULT        = pim2022_pb(I)
          ERR         = sqrt(pim2022_pbstat(I)**2+pim2022_pbsys(I)**2)
          fuu         = pim2022_pbfuu(I)
          fuua        = pim2022_pbfuua(I)
          DIS         = pim2022_pbdis(I)
          WRITE(94,*) z, ptmin, ptmax, MULT, ERR,
     *               fuu, fuua, DIS
      ENDDO
      CLOSE(94)
      ENDIF


      END PROGRAM MASTER
