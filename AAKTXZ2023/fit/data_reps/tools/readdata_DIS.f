      SUBROUTINE READDATA
      IMPLICIT NONE
C-----------STORED VARIABLE DECLARATIONS
C-----------DECLARE VARIABLES.
C-----------NLINES = NUMBER OF LINES SKIPPED.
C-----------NTOT = TOTAL NUMBER OF ROWS.
      REAL*8 X1(1500),X2(1500),X3(1500),X4(1500),X5(1500),X6(1500)
      REAL*8 X7(200),X8(200),X9(200),X10(200),X11(200),X12(200),X13(200)
      REAL*8 X14(200),X15(200),X16(200),X17(200),X18(200),X19(200)
      REAL*8 X20(200)
      REAL*8 Sep_hermes
      REAL*8 s, Qs, Nu, xb, zh, pht, fuu, fuua
      REAL*8 Sep_Jlab
      real*8 Sep_EIC, Sep_EICC
      real*8 Sep_JLAB2022
      double precision M
      CHARACTER*72 FILEEXP
      INTEGER I,J,K,L
      INTEGER NLINES ! NUMBER OF LINES SKIPPED
      INTEGER NTOT, NTOT2, NTOT3 ! NUMBER OF ROWS
      INTEGER NTOT_JLAB
      INTEGER IT,IH,IC
      CHARACTER*8 RT, RT2
      character*2 proc
      real*8 F2light_d,FLlight_d
      real*8 F2light_A,FLlight_A
      real*8 disA,disD
      real*8 F2light,FLlight
      real*8 Q0,Q2,Q
      real*8 nulow, nuhigh
      real*8 Q2low, Q2high
      real*8 ptlow, pthigh
      INTEGER LO, NLO
      INCLUDE './data-inc-DIS.f'
      RT = 'expdata/'

C-----APFEL Perturbative Order Scheme
      LO  = 0
      NLO = 1

      proc = "EM"
      call SetMassScheme("ZM-VFNS")
      call SetProcessDIS(proc)
      call SetProjectileDIS("electron")
      call SetTargetDIS("isoscalar")
      call SetPerturbativeOrder(NLO)

      CALL SETLHAPARM('SILENT')
      CALL EnableWelcomeMessage(.False.)

C-----------IMPORTATION OF ALL HERMES (HERA-DESY) SIDIS DATA
      NTOT = 9 ! 9 ROWS
      NTOT2 = 8! 8 ROWS
      NTOT3 = 22
      NLINES = 1 !SKIP THE FIRST LINE

C------------PARAMETERS
      Q0 = 1.3d0
      M = 0.938272046d0 ! GeV (Mass of Proton)
      Sep_hermes = 2*27.6*M + M**2 ! GeV^2
      Sep_JLAB   = 11.03 ! GeV^2
      Sep_EICC   = 112 ! GeV^2
      Sep_EIC    = 329 !2d0*2d0*10d0*110d0 ! GeV^2
      Sep_JLAB2022 = 2*5.014*M + M**2 ! 2023

      IF(I_EIC_PRED.EQ.1) THEN
      !OPEN(UNIT = 1001,FILE ='expdata/EIC.dat')
      !DO J=1,3 !vertical
        DO K=1,3 !horizontal
          DO L=1,10 !variable

            s = Sep_EIC
            Qs = 6.58
            Q = dsqrt(Qs)

            xb = 0.1+0.4*(L-1)/9!0.1 !2023
            ! Compute DIS CX for proton
            call SetPDFSet("CT18ANLO")
            call SetTargetDIS("isoscalar")
            call InitializeAPFEL_DIS
            CALL ComputeStructureFunctionsAPFEL(Q0,Q)
            F2light_D = F2light(xb)
            FLlight_D = FLlight(xb)
            call dis_cx(s,xb,Qs,F2light_D,FLlight_D,disD)
            ! Compute DIS CX for Au
            call SetPDFSet("EPPS21nlo_CT18Anlo_Au197")
            call SetTargetDIS("proton")
            call InitializeAPFEL_DIS
            CALL ComputeStructureFunctionsAPFEL(Q0,Q)
            F2light_A = F2light(xb)
            FLlight_A = FLlight(xb)
            call dis_cx(s,xb,Qs,F2light_A,FLlight_A,disA)
            !write(1001,*) "Qs DIS EIC"
            !write(1001,*)  Qs, disA/disD

            EIC_PRED_Qs((K-1)*10+L) = Qs
            EIC_PRED_xb((K-1)*10+L) = xb
            EIC_PRED_dis((K-1)*10+L) = disA/disD

          ENDDO
        ENDDO
      !ENDDO

      !CLOSE(1001)
      ENDIF

      IF(I_EICC_PRED.EQ.1) THEN
      !OPEN(UNIT = 1003,FILE ='expdata/EICC.dat')
      !DO J=1,3 !vertical
        DO K=1,3 !horizontal
          DO L=1,10 !variable

            s = Sep_EICC
            Qs = 5.61
            Q = dsqrt(Qs)
            xb = 0.1+0.4*(L-1)/9!0.1 !2023

            ! Compute DIS CX for proton
            call SetPDFSet("CT18ANLO")
            call SetTargetDIS("isoscalar")
            call InitializeAPFEL_DIS
            CALL ComputeStructureFunctionsAPFEL(Q0,Q)
            F2light_D = F2light(xb)
            FLlight_D = FLlight(xb)
            call dis_cx(s,xb,Qs,F2light_D,FLlight_D,disD)
            ! Compute DIS CX for Au
            call SetPDFSet("EPPS21nlo_CT18Anlo_Au197")
            call SetTargetDIS("proton")
            call InitializeAPFEL_DIS
            CALL ComputeStructureFunctionsAPFEL(Q0,Q)
            F2light_A = F2light(xb)
            FLlight_A = FLlight(xb)
            call dis_cx(s,xb,Qs,F2light_A,FLlight_A,disA)
            write(1003,*) "Qs DIS EICC"
            write(1003,*)  Qs, disA/disD

            EICC_PRED_Qs((K-1)*10+L) = Qs
            EICC_PRED_xb((K-1)*10+L) = xb
            EICC_PRED_dis((K-1)*10+L) = disA/disD

          ENDDO
        ENDDO
      !ENDDO

      !CLOSE(1003)
      ENDIF

      IF(I_JLAB_PRED.EQ.1) THEN
      !OPEN(UNIT = 1002,FILE ='expdata/JLAB.dat')
      !DO J=1,3 !vertical
        DO K=1,3 !horizontal
          DO L=1,10 !variable

            s = 11.03
            Qs = 1.65
            xb = 0.1+0.4*(L-1)/9!0.1 !2023

            Q = dsqrt(Qs)
            ! Compute DIS CX for proton
            call SetPDFSet("CT18ANLO")
            call SetTargetDIS("isoscalar")
            call InitializeAPFEL_DIS
            CALL ComputeStructureFunctionsAPFEL(Q0,Q)
            F2light_D = F2light(xb)
            FLlight_D = FLlight(xb)
            call dis_cx(s,xb,Qs,F2light_D,FLlight_D,disD)
            ! Compute DIS CX for Au
            call SetPDFSet("EPPS21nlo_CT18Anlo_Pb208")
            call SetTargetDIS("proton")
            call InitializeAPFEL_DIS
            CALL ComputeStructureFunctionsAPFEL(Q0,Q)
            F2light_A = F2light(xb)
            FLlight_A = FLlight(xb)
            call dis_cx(s,xb,Qs,F2light_A,FLlight_A,disA)
            write(1002,*) "Qs DIS JLAB"
            write(1002,*)  Qs, disA/disD

            JLAB_PRED_Qs((K-1)*10+L) = Qs
            JLAB_PRED_xb((K-1)*10+L) = xb
            JLAB_PRED_dis((K-1)*10+L) = disA/disD

          ENDDO
        ENDDO
      !ENDDO

      !CLOSE(1002)
      ENDIF

C-----triple plot

      s=802d0

      IF(I_x_PRED.EQ.1) THEN
            NLINES=1
            NTOT=270 !2023
            FILEEXP = RT//'PRED/x.dat'
            CALL READDATA6(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,NTOT)
            DO I=1,NTOT
                  x_PRED_Qs(I)  = X1(I)
                  x_PRED_xb(I)  = X2(I)
                  x_PRED_zh(I)  = X3(I)
                  x_PRED_pht(I) = X4(I)
                  Qs  =X1(I)
                  xb  =X2(I)
                  zh  =X3(I)
                  pht =X4(I)
                  IH = 1
                  IC = 1
                  IT = 197
                  call SetPDFSet("CT18ANLO")
                  CALL InitializeAPFEL_DIS
                  call SetTargetDIS("isoscalar")
                  CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
                  F2light_d = F2light(xb)
                  FLlight_d = FLlight(xb)
                  call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
                  call SetPDFSet("EPPS21nlo_CT18Anlo_Au197")
                  CALL InitializeAPFEL_DIS
                  call SetTargetDIS("proton")
                  CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
                  F2light_A = F2light(xb)
                  FLlight_A = FLlight(xb)
                  call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
                  x_PRED_dis(I)= disA/disD
            END DO
            print *, "x done"
      END IF

      IF(I_z_PRED.EQ.1) THEN
            NLINES=1
            NTOT=270
            FILEEXP = RT//'PRED/z.dat'
            CALL READDATA6(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,NTOT)
            DO I=1,NTOT
                  z_PRED_Qs(I)  = X1(I)
                  z_PRED_xb(I)  = X2(I)
                  z_PRED_zh(I)  = X3(I)
                  z_PRED_pht(I) = X4(I)
                  Qs  =X1(I)
                  xb  =X2(I)
                  zh  =X3(I)
                  pht =X4(I)
                  IH = 1
                  IC = 1
                  IT = 197
                  call SetPDFSet("CT18ANLO")
                  CALL InitializeAPFEL_DIS
                  call SetTargetDIS("isoscalar")
                  CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
                  F2light_d = F2light(xb)
                  FLlight_d = FLlight(xb)
                  call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
                  call SetPDFSet("EPPS21nlo_CT18Anlo_Au197")
                  CALL InitializeAPFEL_DIS
                  call SetTargetDIS("proton")
                  CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
                  F2light_A = F2light(xb)
                  FLlight_A = FLlight(xb)
                  call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
                  z_PRED_dis(I)= disA/disD
            END DO
            print *, "z done"
      END IF

      IF(I_pht_PRED.EQ.1) THEN
            NLINES=1
            NTOT=270
            FILEEXP = RT//'PRED/pht.dat'
            CALL READDATA6(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,NTOT)
            DO I=1,NTOT
                  pht_PRED_Qs(I)  = X1(I)
                  pht_PRED_xb(I)  = X2(I)
                  pht_PRED_zh(I)  = X3(I)
                  pht_PRED_pht(I) = X4(I)
                  Qs  =X1(I)
                  xb  =X2(I)
                  zh  =X3(I)
                  pht =X4(I)
                  IH = 1
                  IC = 1
                  IT = 197
                  call SetPDFSet("CT18ANLO")
                  CALL InitializeAPFEL_DIS
                  call SetTargetDIS("isoscalar")
                  CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
                  F2light_d = F2light(xb)
                  FLlight_d = FLlight(xb)
                  call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
                  call SetPDFSet("EPPS21nlo_CT18Anlo_Au197")
                  CALL InitializeAPFEL_DIS
                  call SetTargetDIS("proton")
                  CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
                  F2light_A = F2light(xb)
                  FLlight_A = FLlight(xb)
                  call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
                  pht_PRED_dis(I)= disA/disD
            END DO
            print *, "pht done"
      END IF

C-----JLAB-pre
      NTOT_JLAB = 180
      IF(I_JLAB_PIP_12.EQ.1) THEN
            FILEEXP = RT//'JLAB12/piplusCronin_singleHadron.csv'
            CALL READDATA6(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,NTOT_JLAB)
            DO I=1,NTOT_JLAB
               JLAB_PIP_TARGET(I)      =X1(I)
               JLAB_PIP_Z(I)           =X2(I)
               JLAB_PIP_PTLOW(I)       =X3(I)
               JLAB_PIP_PTHIGH(I)      =X4(I)
               JLAB_PIP_MULT_RATIO(I)  =X5(I)
               JLAB_PIP_ERR(I)         =X6(I)
               Q2LOW  = 1d0
               Q2HIGH = 4.1d0
               NULOW  = 2.2d0
               NUHIGH = 4.2d0
               PTLOW = X3(I)
               PTHIGH = X4(I)
               Qs = 0.5d0*(Q2LOW+Q2HIGH)
               Nu = 0.5d0*(NULOW+NUHIGH)
               xb  = Qs/(2*M*Nu)
               zh  = X2(I)/10d0
               pht = 0.5d0*(PTLOW+PTHIGH)
               IH = 1
               IC = 1
               IF(X1(I).eq.0) THEN
               IT = 12
               call SetPDFSet("EPPS21nlo_CT18ANLO_C12")
               ELSEIF(X1(I).eq.1) THEN
               IT = 56
               call SetPDFSet("EPPS21nlo_CT18ANLO_Fe56")
               ELSEIF(X1(I).eq.2) THEN
               IT = 208
               call SetPDFSet("EPPS21nlo_CT18ANLO_Pb208")
               ENDIF
               ! Computing DIS CX for e+A
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
               !Computing DIS CX for e+d
               call SetPDFSet("CT18ANLO")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
               ! Computing SIDIS CX for e+A
               CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
               ! Computing SIDIS CX for d+A
               CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
               JLAB_PIP_FUU(I) = fuu
               JLAB_PIP_FUUA(I) = fuua
               JLAB_PIP_DIS(I)= disA/disD
            END DO
      END IF


      NTOT_JLAB = 180
      IF(I_JLAB_PIM_12.EQ.1) THEN
            FILEEXP = RT//'JLAB12/piminusCronin_singleHadron.txt'
            CALL READDATA6(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,NTOT_JLAB)
            DO I=1,NTOT_JLAB
               JLAB_PIM_TARGET(I)      =X1(I)
               JLAB_PIM_Z(I)           =X2(I)
               JLAB_PIM_PTLOW(I)       =X3(I)
               JLAB_PIM_PTHIGH(I)      =X4(I)
               JLAB_PIM_MULT_RATIO(I)  =X5(I)
               JLAB_PIM_ERR(I)         =X6(I)
               Q2LOW  = 1d0
               Q2HIGH = 4.1d0
               NULOW  = 2.2d0
               NUHIGH = 4.2d0
               PTLOW = X3(I)
               PTHIGH = X4(I)
               Qs = 0.5d0*(Q2LOW+Q2HIGH)
               Nu = 0.5d0*(NULOW+NUHIGH)
               xb  = Qs/(2*M*Nu)
               zh  = X2(I)/10d0
               pht = 0.5d0*(PTLOW+PTHIGH)
               IH = 1
               IC = -1
               IF(X1(I).eq.0) THEN
               IT = 12
               call SetPDFSet("EPPS21nlo_CT18ANLO_C12")
               ELSEIF(X1(I).eq.1) THEN
               IT = 56
               call SetPDFSet("EPPS21nlo_CT18ANLO_Fe56")
               ELSEIF(X1(I).eq.2) THEN
               IT = 208
               call SetPDFSet("EPPS21nlo_CT18ANLO_Pb208")
               ENDIF
               ! Computing DIS CX for e+A
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
               !Computing DIS CX for e+d
               call SetPDFSet("CT18ANLO")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
               ! Computing SIDIS CX for e+A
               CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
               ! Computing SIDIS CX for d+A
               CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
               JLAB_PIM_FUU(I) = fuu
               JLAB_PIM_FUUA(I) = fuua
               JLAB_PIM_DIS(I)= disA/disD
            END DO
      END IF



      s = SEP_HERMES

C------PI+ Z DEPENDENT DATA------------------------
C------- HE
      IF(I_PIP_HE_Z.EQ.1) THEN
            FILEEXP = RT//'HERMES/pip_he_z.csv'
            CALL READDATA7(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7,NTOT)
            DO I=1,NTOT
               PIP_HE_Z_VALUE(I)   =X1(I)
               PIP_HE_Z_NU(I)      =X2(I)
               PIP_HE_Z_Q2(I)      =X3(I)
               PIP_HE_Z_PT2(I)     =X4(I)
               PIP_HE_Z_MULT(I)    =X5(I)
               PIP_HE_Z_STAT(I)    =X6(I)
               PIP_HE_Z_SYS(I)     =X7(I)
               Qs  = X3(I)
               Nu  = X2(I)
               xb  = Qs/(2*M*Nu)
               zh  = X1(I)
               pht = sqrt(X4(I))
               IT = 4
               IH = 1
               IC = 1
               CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIP_HE_Z_FUU(I) = fuu
               PIP_HE_Z_FUUA(I) = fuua
               call SetPDFSet("CT18ANLO")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS21nlo_CT18Anlo_He4")
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
               PIP_HE_Z_DIS(I)= disA/disD
            END DO
      END IF
C------- NE
      IF(I_PIP_NE_Z.EQ.1) THEN
            FILEEXP = RT//'HERMES/pip_ne_z.csv'
            CALL READDATA7(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7,NTOT)
            DO I=1,NTOT
               PIP_NE_Z_VALUE(I)   =X1(I)
               PIP_NE_Z_NU(I)       =X2(I)
               PIP_NE_Z_Q2(I)      =X3(I)
               PIP_NE_Z_PT2(I)     =X4(I)
               PIP_NE_Z_MULT(I)    =X5(I)
               PIP_NE_Z_STAT(I)    =X6(I)
               PIP_NE_Z_SYS(I)     =X7(I)
               Qs  = X3(I)
               Nu  = X2(I)
               xb  = Qs/(2*M*Nu)
               zh  = X1(I)
               pht = sqrt(X4(I))
               IT = 20
               IH = 1
               IC = 1
               CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIP_NE_Z_FUU(I) = fuu
               PIP_NE_Z_FUUA(I) = fuua
               call SetPDFSet("CT18ANLO")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS21nlo_CT18Anlo_Ne20")
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
               PIP_NE_Z_DIS(I)= disA/disD
            END DO
      END IF
C------- KR
      IF(I_PIP_KR_Z.EQ.1) THEN
            FILEEXP = RT//'HERMES/pip_kr_z.csv'
            CALL READDATA7(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7,NTOT)
            DO I=1,NTOT
               PIP_KR_Z_VALUE(I)   =X1(I)
               PIP_KR_Z_NU(I)      =X2(I)
               PIP_KR_Z_Q2(I)      =X3(I)
               PIP_KR_Z_PT2(I)     =X4(I)
               PIP_KR_Z_MULT(I)    =X5(I)
               PIP_KR_Z_STAT(I)    =X6(I)
               PIP_KR_Z_SYS(I)     =X7(I)
               Qs  = X3(I)
               Nu  = X2(I)
               xb  = Qs/(2*M*Nu)
               zh  = X1(I)
               pht = sqrt(X4(I))
               IT = 84
               IH = 1
               IC = 1
               CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIP_KR_Z_FUU(I) = fuu
               PIP_KR_Z_FUUA(I) = fuua
               call SetPDFSet("CT18ANLO")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS21nlo_CT18Anlo_Kr84")
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
               PIP_KR_Z_DIS(I)= disA/disD
            END DO
      END IF
C------- XE
      IF(I_PIP_XE_Z.EQ.1) THEN
            FILEEXP = RT//'HERMES/pip_xe_z.csv'
            CALL READDATA7(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7,NTOT)
            DO I=1,NTOT
               PIP_XE_Z_VALUE(I)   =X1(I)
               PIP_XE_Z_NU(I)      =X2(I)
               PIP_XE_Z_Q2(I)      =X3(I)
               PIP_XE_Z_PT2(I)     =X4(I)
               PIP_XE_Z_MULT(I)    =X5(I)
               PIP_XE_Z_STAT(I)    =X6(I)
               PIP_XE_Z_SYS(I)     =X7(I)
               Qs  = X3(I)
               Nu  = X2(I)
               xb  = Qs/(2*M*Nu)
               zh  = X1(I)
               pht = sqrt(X4(I))
               IT = 131
               IH = 1
               IC = 1
               CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIP_XE_Z_FUU(I) = fuu
               PIP_XE_Z_FUUA(I) = fuua
               call SetPDFSet("CT18ANLO")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS21nlo_CT18Anlo_Xe131")
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
               PIP_XE_Z_DIS(I)= disA/disD
            END DO
      END IF
C------PI- Z DEPENDENT DATA------------------------
C------- HE
      IF(I_PIM_HE_Z.EQ.1) THEN
            FILEEXP = RT//'HERMES/pim_he_z.csv'
            CALL READDATA7(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7,NTOT)
            DO I=1,NTOT
               PIM_HE_Z_VALUE(I)   =X1(I)
               PIM_HE_Z_NU(I)       =X2(I)
               PIM_HE_Z_Q2(I)      =X3(I)
               PIM_HE_Z_PT2(I)     =X4(I)
               PIM_HE_Z_MULT(I)    =X5(I)
               PIM_HE_Z_STAT(I)    =X6(I)
               PIM_HE_Z_SYS(I)     =X7(I)
               Qs  = X3(I)
               Nu  = X2(I)
               xb  = Qs/(2*M*Nu)
               zh  = X1(I)
               pht = sqrt(X4(I))
               IT = 4
               IH = 1
               IC = -1
               CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIM_HE_Z_FUU(I) = fuu
               PIM_HE_Z_FUUA(I) = fuua
               call SetPDFSet("CT18ANLO")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS21nlo_CT18Anlo_He4")
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
               PIM_HE_Z_DIS(I)= disA/disD
            END DO
      END IF
C------- NE
      IF(I_PIM_NE_Z.EQ.1) THEN
            FILEEXP = RT//'HERMES/pim_ne_z.csv'
            CALL READDATA7(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7,NTOT)
            DO I=1,NTOT
               PIM_NE_Z_VALUE(I)   =X1(I)
               PIM_NE_Z_NU(I)       =X2(I)
               PIM_NE_Z_Q2(I)      =X3(I)
               PIM_NE_Z_PT2(I)     =X4(I)
               PIM_NE_Z_MULT(I)    =X5(I)
               PIM_NE_Z_STAT(I)    =X6(I)
               PIM_NE_Z_SYS(I)     =X7(I)
               Qs  = X3(I)
               Nu  = X2(I)
               xb  = Qs/(2*M*Nu)
               zh  = X1(I)
               pht = sqrt(X4(I))
               IT = 20
               IH = 1
               IC = -1
               CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIM_NE_Z_FUU(I) = fuu
               PIM_NE_Z_FUUA(I) = fuua
               call SetPDFSet("CT18ANLO")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS21nlo_CT18Anlo_Ne20")
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
               PIM_NE_Z_DIS(I)= disA/disD
            END DO
      END IF
C------- KR
      IF(I_PIM_KR_Z.EQ.1) THEN
            FILEEXP = RT//'HERMES/pim_kr_z.csv'
            CALL READDATA7(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7,NTOT)
            DO I=1,NTOT
               PIM_KR_Z_VALUE(I)   =X1(I)
               PIM_KR_Z_NU(I)       =X2(I)
               PIM_KR_Z_Q2(I)      =X3(I)
               PIM_KR_Z_PT2(I)     =X4(I)
               PIM_KR_Z_MULT(I)    =X5(I)
               PIM_KR_Z_STAT(I)    =X6(I)
               PIM_KR_Z_SYS(I)     =X7(I)
               Qs  = X3(I)
               Nu  = X2(I)
               xb  = Qs/(2*M*Nu)
               zh  = X1(I)
               pht = sqrt(X4(I))
               IT = 84
               IH = 1
               IC = -1
               CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIM_KR_Z_FUU(I) = fuu
               PIM_KR_Z_FUUA(I) = fuua
               call SetPDFSet("CT18ANLO")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS21nlo_CT18Anlo_Kr84")
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
               PIM_KR_Z_DIS(I)= disA/disD
            END DO
      END IF
C------- XE
      IF(I_PIM_XE_Z.EQ.1) THEN
            FILEEXP = RT//'HERMES/pim_xe_z.csv'
            CALL READDATA7(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7,NTOT)
            DO I=1,NTOT
               PIM_XE_Z_VALUE(I)   =X1(I)
               PIM_XE_Z_NU(I)      =X2(I)
               PIM_XE_Z_Q2(I)      =X3(I)
               PIM_XE_Z_PT2(I)     =X4(I)
               PIM_XE_Z_MULT(I)    =X5(I)
               PIM_XE_Z_STAT(I)    =X6(I)
               PIM_XE_Z_SYS(I)     =X7(I)
               Qs  = X3(I)
               Nu  = X2(I)
               xb  = Qs/(2*M*Nu)
               zh  = X1(I)
               pht = sqrt(X4(I))
               IT = 131
               IH = 1
               IC = -1
               CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIM_XE_Z_FUU(I) = fuu
               PIM_XE_Z_FUUA(I) = fuua
               call SetPDFSet("CT18ANLO")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS21nlo_CT18Anlo_Xe131")
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
               PIM_XE_Z_DIS(I)= disA/disD
            END DO
      END IF
C------PI0 Z DEPENDENT DATA------------------------
C------- HE
      IF(I_PI0_HE_Z.EQ.1) THEN
            FILEEXP = RT//'HERMES/pi0_he_z.csv'
            CALL READDATA7(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7,NTOT)
            DO I=1,NTOT
               PI0_HE_Z_VALUE(I)   =X1(I)
               PI0_HE_Z_NU(I)      =X2(I)
               PI0_HE_Z_Q2(I)      =X3(I)
               PI0_HE_Z_PT2(I)     =X4(I)
               PI0_HE_Z_MULT(I)    =X5(I)
               PI0_HE_Z_STAT(I)    =X6(I)
               PI0_HE_Z_SYS(I)     =X7(I)
               Qs  = X3(I)
               Nu  = X2(I)
               xb  = Qs/(2*M*Nu)
               zh  = X1(I)
               pht = sqrt(X4(I))
               IT = 4
               IH = 1
               IC = 0
               CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PI0_HE_Z_FUU(I) = fuu
               PI0_HE_Z_FUUA(I) = fuua
               call SetPDFSet("CT18ANLO")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS21nlo_CT18Anlo_He4")
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
               PI0_HE_Z_DIS(I)= disA/disD
            END DO
      END IF
C------- NE
      IF(I_PI0_NE_Z.EQ.1) THEN
            FILEEXP = RT//'HERMES/pi0_ne_z.csv'
            CALL READDATA7(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7,NTOT)
            DO I=1,NTOT
               PI0_NE_Z_VALUE(I)   =X1(I)
               PI0_NE_Z_NU(I)       =X2(I)
               PI0_NE_Z_Q2(I)      =X3(I)
               PI0_NE_Z_PT2(I)     =X4(I)
               PI0_NE_Z_MULT(I)    =X5(I)
               PI0_NE_Z_STAT(I)    =X6(I)
               PI0_NE_Z_SYS(I)     =X7(I)
               Qs  = X3(I)
               Nu  = X2(I)
               xb  = Qs/(2*M*Nu)
               zh  = X1(I)
               pht = sqrt(X4(I))
               IT = 20
               IH = 1
               IC = 0
               CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PI0_NE_Z_FUU(I) = fuu
               PI0_NE_Z_FUUA(I) = fuua
               call SetPDFSet("CT18ANLO")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS21nlo_CT18Anlo_Ne20")
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
               PI0_NE_Z_DIS(I)= disA/disD
            END DO
      END IF
C------- KR
      IF(I_PI0_KR_Z.EQ.1) THEN
            FILEEXP = RT//'HERMES/pi0_kr_z.csv'
            CALL READDATA7(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7,NTOT)
            DO I=1,NTOT
               PI0_KR_Z_VALUE(I)   =X1(I)
               PI0_KR_Z_NU(I)       =X2(I)
               PI0_KR_Z_Q2(I)      =X3(I)
               PI0_KR_Z_PT2(I)     =X4(I)
               PI0_KR_Z_MULT(I)    =X5(I)
               PI0_KR_Z_STAT(I)    =X6(I)
               PI0_KR_Z_SYS(I)     =X7(I)
               Qs  = X3(I)
               Nu  = X2(I)
               xb  = Qs/(2*M*Nu)
               zh  = X1(I)
               pht = sqrt(X4(I))
               IT = 84
               IH = 1
               IC = 0
               CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PI0_KR_Z_FUU(I) = fuu
               PI0_KR_Z_FUUA(I) = fuua
               call SetPDFSet("CT18ANLO")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS21nlo_CT18Anlo_Kr84")
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
               PI0_KR_Z_DIS(I)= disA/disD
            END DO
      END IF
C------- XE
      IF(I_PI0_XE_Z.EQ.1) THEN
            FILEEXP = RT//'HERMES/pi0_xe_z.csv'
            CALL READDATA7(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7,NTOT)
            DO I=1,NTOT
               PI0_XE_Z_VALUE(I)   =X1(I)
               PI0_XE_Z_NU(I)      =X2(I)
               PI0_XE_Z_Q2(I)      =X3(I)
               PI0_XE_Z_PT2(I)     =X4(I)
               PI0_XE_Z_MULT(I)    =X5(I)
               PI0_XE_Z_STAT(I)    =X6(I)
               PI0_XE_Z_SYS(I)     =X7(I)
               Qs  = X3(I)
               Nu  = X2(I)
               xb  = Qs/(2*M*Nu)
               zh  = X1(I)
               pht = sqrt(X4(I))
               IT = 131
               IH = 1
               IC = 0
               CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PI0_XE_Z_FUU(I) = fuu
               PI0_XE_Z_FUUA(I) = fuua
               call SetPDFSet("CT18ANLO")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS21nlo_CT18Anlo_Xe131")
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
               PI0_XE_Z_DIS(I)= disA/disD
            END DO
      END IF

C-----------------------------------------
C------PI0 PT2 DATA-----------------------
C----------- HE
      IF(I_PI0_HE_PT2.EQ.1) THEN
            FILEEXP = RT//'HERMES/pi0_he_pt2.csv'
            CALL READDATA7(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7,NTOT)
            DO I=1,NTOT
               PI0_HE_PT2_VALUE(I)   =X1(I)
               PI0_HE_PT2_NU(I)      =X2(I)
               PI0_HE_PT2_Z(I)       =X3(I)
               PI0_HE_PT2_Q2(I)      =X4(I)
               PI0_HE_PT2_MULT(I)    =X5(I)
               PI0_HE_PT2_STAT(I)    =X6(I)
               PI0_HE_PT2_SYS(I)     =X7(I)
               Qs  = X4(I)
               Nu  = X2(I)
               xb  = Qs/(2*M*Nu)
               zh  = X3(I)
               pht = sqrt(X1(I))
               IT = 4
               IH = 1
               IC = 0
               CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PI0_HE_PT2_FUU(I) = fuu
               PI0_HE_PT2_FUUA(I) = fuua
               call SetPDFSet("CT18ANLO")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("isoscalar")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS21nlo_CT18Anlo_He4")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("proton")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
               PI0_HE_PT2_DIS(I)= disA/disD
            END DO
      END IF
C----------- NE
      IF(I_PI0_NE_PT2.EQ.1) THEN
            FILEEXP = RT//'HERMES/pi0_ne_pt2.csv'
            CALL READDATA7(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7,NTOT)
            DO I=1,NTOT
               PI0_NE_PT2_VALUE(I)   =X1(I)
               PI0_NE_PT2_NU(I)      =X2(I)
               PI0_NE_PT2_Z(I)       =X3(I)
               PI0_NE_PT2_Q2(I)      =X4(I)
               PI0_NE_PT2_MULT(I)    =X5(I)
               PI0_NE_PT2_STAT(I)    =X6(I)
               PI0_NE_PT2_SYS(I)     =X7(I)
               Qs  = X4(I)
               Nu  = X2(I)
               xb  = Qs/(2*M*Nu)
               zh  = X3(I)
               pht = sqrt(X1(I))
               IT = 20
               IH = 1
               IC = 0
               CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PI0_NE_PT2_FUU(I) = fuu
               PI0_NE_PT2_FUUA(I) = fuua
               call SetPDFSet("CT18ANLO")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("isoscalar")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS21nlo_CT18Anlo_Ne20")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("proton")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
               PI0_NE_PT2_DIS(I)= disA/disD
            END DO
      END IF
C----------- KR
      IF(I_PI0_KR_PT2.EQ.1) THEN
            FILEEXP = RT//'HERMES/pi0_kr_pt2.csv'
            CALL READDATA7(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7,NTOT)
            DO I=1,NTOT
               PI0_KR_PT2_VALUE(I)   =X1(I)
               PI0_KR_PT2_NU(I)      =X2(I)
               PI0_KR_PT2_Z(I)       =X3(I)
               PI0_KR_PT2_Q2(I)      =X4(I)
               PI0_KR_PT2_MULT(I)    =X5(I)
               PI0_KR_PT2_STAT(I)    =X6(I)
               PI0_KR_PT2_SYS(I)     =X7(I)
               Qs  = X4(I)
               Nu  = X2(I)
               xb  = Qs/(2*M*Nu)
               zh  = X3(I)
               pht = sqrt(X1(I))
               IT = 84
               IH = 1
               IC = 0
               CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PI0_KR_PT2_FUU(I) = fuu
               PI0_KR_PT2_FUUA(I) = fuua
               call SetPDFSet("CT18ANLO")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("isoscalar")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS21nlo_CT18Anlo_Kr84")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("proton")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
               PI0_KR_PT2_DIS(I)= disA/disD
            END DO
      END IF
C---------- XE
      IF(I_PI0_XE_PT2.EQ.1) THEN
            FILEEXP = RT//'HERMES/pi0_xe_pt2.csv'
            CALL READDATA7(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7,NTOT)
            DO I=1,NTOT
               PI0_XE_PT2_VALUE(I)   =X1(I)
               PI0_XE_PT2_NU(I)      =X2(I)
               PI0_XE_PT2_Z(I)       =X3(I)
               PI0_XE_PT2_Q2(I)      =X4(I)
               PI0_XE_PT2_MULT(I)    =X5(I)
               PI0_XE_PT2_STAT(I)    =X6(I)
               PI0_XE_PT2_SYS(I)     =X7(I)
               Qs  = X4(I)
               Nu  = X2(I)
               xb  = Qs/(2*M*Nu)
               zh  = X3(I)
               pht = sqrt(X1(I))
               IT = 131
               IH = 1
               IC = 0
               CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PI0_XE_PT2_FUU(I) = fuu
               PI0_XE_PT2_FUUA(I) = fuua
               call SetPDFSet("CT18ANLO")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("isoscalar")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS21nlo_CT18Anlo_Xe131")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("proton")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
               PI0_XE_PT2_DIS(I)= disA/disD
            END DO
      END IF

C-----------------------------------------
C------PI+ PT2 DATA-----------------------
C--------- HE
      IF(I_PIP_HE_PT2.EQ.1) THEN
            FILEEXP = RT//'HERMES/pip_he_pt2.csv'
            CALL READDATA7(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7,NTOT)
            DO I=1,NTOT
               PIP_HE_PT2_VALUE(I)   =X1(I)
               PIP_HE_PT2_NU(I)      =X2(I)
               PIP_HE_PT2_Z(I)       =X3(I)
               PIP_HE_PT2_Q2(I)      =X4(I)
               PIP_HE_PT2_MULT(I)    =X5(I)
               PIP_HE_PT2_STAT(I)    =X6(I)
               PIP_HE_PT2_SYS(I)     =X7(I)
               Qs  = X4(I)
               Nu  = X2(I)
               xb  = Qs/(2*M*Nu)
               zh  = X3(I)
               pht = sqrt(X1(I))
               IT = 4
               IH = 1
               IC = 1
               CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIP_HE_PT2_FUU(I) = fuu
               PIP_HE_PT2_FUUA(I) = fuua
               call SetPDFSet("CT18ANLO")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("isoscalar")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS21nlo_CT18Anlo_He4")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("proton")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
               PIP_HE_PT2_DIS(I)= disA/disD
            END DO
      END IF
C--------- NE
      IF(I_PIP_NE_PT2.EQ.1) THEN
            FILEEXP = RT//'HERMES/pip_ne_pt2.csv'
            CALL READDATA7(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7,NTOT)
            DO I=1,NTOT
               PIP_NE_PT2_VALUE(I)   =X1(I)
               PIP_NE_PT2_NU(I)      =X2(I)
               PIP_NE_PT2_Z(I)       =X3(I)
               PIP_NE_PT2_Q2(I)      =X4(I)
               PIP_NE_PT2_MULT(I)    =X5(I)
               PIP_NE_PT2_STAT(I)    =X6(I)
               PIP_NE_PT2_SYS(I)     =X7(I)
               Qs  = X4(I)
               Nu  = X2(I)
               xb  = Qs/(2*M*Nu)
               zh  = X3(I)
               pht = sqrt(X1(I))
               IT = 20
               IH = 1
               IC = 1
               CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIP_NE_PT2_FUU(I) = fuu
               PIP_NE_PT2_FUUA(I) = fuua
               call SetPDFSet("CT18ANLO")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("isoscalar")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS21nlo_CT18Anlo_Ne20")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("proton")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
               PIP_NE_PT2_DIS(I)= disA/disD
            END DO
      END IF
C----------- KR
      IF(I_PIP_KR_PT2.EQ.1) THEN
            FILEEXP = RT//'HERMES/pip_kr_pt2.csv'
            CALL READDATA7(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7,NTOT)
            DO I=1,NTOT
               PIP_KR_PT2_VALUE(I)   =X1(I)
               PIP_KR_PT2_NU(I)      =X2(I)
               PIP_KR_PT2_Z(I)       =X3(I)
               PIP_KR_PT2_Q2(I)      =X4(I)
               PIP_KR_PT2_MULT(I)    =X5(I)
               PIP_KR_PT2_STAT(I)    =X6(I)
               PIP_KR_PT2_SYS(I)     =X7(I)
               Qs  = X4(I)
               Nu  = X2(I)
               xb  = Qs/(2*M*Nu)
               zh  = X3(I)
               pht = sqrt(X1(I))
               IT = 84
               IH = 1
               IC = 1
               CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIP_KR_PT2_FUU(I) = fuu
               PIP_KR_PT2_FUUA(I) = fuua
               call SetPDFSet("CT18ANLO")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("isoscalar")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS21nlo_CT18Anlo_Kr84")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("proton")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
               PIP_KR_PT2_DIS(I)= disA/disD
            END DO
      END IF
C----------- XE
      IF(I_PIP_XE_PT2.EQ.1) THEN
            FILEEXP = RT//'HERMES/pip_xe_pt2.csv'
            CALL READDATA7(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7,NTOT)
            DO I=1,NTOT
               PIP_XE_PT2_VALUE(I)   =X1(I)
               PIP_XE_PT2_NU(I)      =X2(I)
               PIP_XE_PT2_Z(I)       =X3(I)
               PIP_XE_PT2_Q2(I)      =X4(I)
               PIP_XE_PT2_MULT(I)    =X5(I)
               PIP_XE_PT2_STAT(I)    =X6(I)
               PIP_XE_PT2_SYS(I)     =X7(I)
               Qs  = X4(I)
               Nu  = X2(I)
               xb  = Qs/(2*M*Nu)
               zh  = X3(I)
               pht = sqrt(X1(I))
               IT = 131
               IH = 1
               IC = 1
               CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIP_XE_PT2_FUU(I) = fuu
               PIP_XE_PT2_FUUA(I) = fuua
               call SetPDFSet("CT18ANLO")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("isoscalar")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS21nlo_CT18Anlo_Xe131")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("proton")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
               PIP_XE_PT2_DIS(I)= disA/disD
            END DO
      END IF
C
C------PI- PT2 DATA-----------------------
C------------ HE
      IF(I_PIM_HE_PT2.EQ.1) THEN
            FILEEXP = RT//'HERMES/pim_he_pt2.csv'
            CALL READDATA7(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7,NTOT)
            DO I=1,NTOT
               PIM_HE_PT2_VALUE(I)   =X1(I)
               PIM_HE_PT2_NU(I)      =X2(I)
               PIM_HE_PT2_Z(I)       =X3(I)
               PIM_HE_PT2_Q2(I)      =X4(I)
               PIM_HE_PT2_MULT(I)    =X5(I)
               PIM_HE_PT2_STAT(I)    =X6(I)
               PIM_HE_PT2_SYS(I)     =X7(I)
               Qs  = X4(I)
               Nu  = X2(I)
               xb  = Qs/(2*M*Nu)
               zh  = X3(I)
               pht = sqrt(X1(I))
               IT = 4
               IH = 1
               IC = -1
               CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIM_HE_PT2_FUU(I) = fuu
               PIM_HE_PT2_FUUA(I) = fuua
               call SetPDFSet("CT18ANLO")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("isoscalar")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS21nlo_CT18Anlo_He4")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("proton")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
               PIM_HE_PT2_DIS(I)= disA/disD
            END DO
      END IF
C------------- NE
      IF(I_PIM_NE_PT2.EQ.1) THEN
            FILEEXP = RT//'HERMES/pim_ne_pt2.csv'
            CALL READDATA7(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7,NTOT)
            DO I=1,NTOT
               PIM_NE_PT2_VALUE(I)   =X1(I)
               PIM_NE_PT2_NU(I)      =X2(I)
               PIM_NE_PT2_Z(I)       =X3(I)
               PIM_NE_PT2_Q2(I)      =X4(I)
               PIM_NE_PT2_MULT(I)    =X5(I)
               PIM_NE_PT2_STAT(I)    =X6(I)
               PIM_NE_PT2_SYS(I)     =X7(I)
               Qs  = X4(I)
               Nu  = X2(I)
               xb  = Qs/(2*M*Nu)
               zh  = X3(I)
               pht = sqrt(X1(I))
               IT = 20
               IH = 1
               IC = -1
               CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIM_NE_PT2_FUU(I) = fuu
               PIM_NE_PT2_FUUA(I) = fuua
               call SetPDFSet("CT18ANLO")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("isoscalar")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS21nlo_CT18Anlo_Ne20")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("proton")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
               PIM_NE_PT2_DIS(I)= disA/disD
            END DO
      END IF
C------------- KR
      IF(I_PIM_KR_PT2.EQ.1) THEN
            FILEEXP = RT//'HERMES/pim_kr_pt2.csv'
            CALL READDATA7(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7,NTOT)
            DO I=1,NTOT
               PIM_KR_PT2_VALUE(I)   =X1(I)
               PIM_KR_PT2_NU(I)      =X2(I)
               PIM_KR_PT2_Z(I)       =X3(I)
               PIM_KR_PT2_Q2(I)      =X4(I)
               PIM_KR_PT2_MULT(I)    =X5(I)
               PIM_KR_PT2_STAT(I)    =X6(I)
               PIM_KR_PT2_SYS(I)     =X7(I)
               Qs  = X4(I)
               Nu  = X2(I)
               xb  = Qs/(2*M*Nu)
               zh  = X3(I)
               pht = sqrt(X1(I))
               IT = 84
               IH = 1
               IC = -1
               CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIM_KR_PT2_FUU(I) = fuu
               PIM_KR_PT2_FUUA(I) = fuua
               call SetPDFSet("CT18ANLO")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("isoscalar")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS21nlo_CT18Anlo_Kr84")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("proton")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
               PIM_KR_PT2_DIS(I)= disA/disD
            END DO
      END IF
C------------- XE
      IF(I_PIM_XE_PT2.EQ.1) THEN
            FILEEXP = RT//'HERMES/pim_xe_pt2.csv'
            CALL READDATA7(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7,NTOT)
            DO I=1,NTOT
               PIM_XE_PT2_VALUE(I)   =X1(I)
               PIM_XE_PT2_NU(I)      =X2(I)
               PIM_XE_PT2_Z(I)       =X3(I)
               PIM_XE_PT2_Q2(I)      =X4(I)
               PIM_XE_PT2_MULT(I)    =X5(I)
               PIM_XE_PT2_STAT(I)    =X6(I)
               PIM_XE_PT2_SYS(I)     =X7(I)
               Qs  = X4(I)
               Nu  = X2(I)
               xb  = Qs/(2*M*Nu)
               zh  = X3(I)
               pht = sqrt(X1(I))
               IT = 131
               IH = 1
               IC = -1
               CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIM_XE_PT2_FUU(I) = fuu
               PIM_XE_PT2_FUUA(I) = fuua
               call SetPDFSet("CT18ANLO")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("isoscalar")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS21nlo_CT18Anlo_Xe131")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("proton")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
               PIM_XE_PT2_DIS(I)= disA/disD
            END DO
      END IF

C--------JLAB2022 DATA

      s = Sep_JLAB2022

C---------pi+
C----------C
      IF(I_pip_2022.EQ.1) THEN
            NLINES=1
            NTOT=48
            FILEEXP = RT//'JLAB2022/pi+.csv'
            CALL READDATA20(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, X11, X12, X13,
     *                     X14,X15,X16,X17,X18,X19,X20,NTOT)
            DO I=1,NTOT
                  pip2022_zlow(I)   =X1(I)
                  pip2022_zup(I)    =X2(I)
                  pip2022_ptlow(I)  =X4(I)
                  pip2022_ptup(I)   =X5(I)
                  pip2022_c(I)      =X6(I)
                  pip2022_cstat(I)  =X7(I)
                  pip2022_csys(I)   =X8(I)
                  pip2022_qc(I)     =X15(I)
                  pip2022_xc(I)     =x18(I) 
                  Qs = pip2022_qc(I)
                  xb  = pip2022_xc(I)
                  zh  = 0.5d0*(pip2022_zlow(I)+pip2022_zup(I))
                  pht = sqrt(0.5d0*(pip2022_ptlow(I)+pip2022_ptup(I)))
                  IH = 1
                  IC = 1
                  IT = 12
                  CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
                  CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
                  pip2022_cfuu(I) = fuu
                  pip2022_cfuua(I) = fuua
                  call SetPDFSet("CT18ANLO")
                  CALL InitializeAPFEL_DIS
                  call SetTargetDIS("isoscalar")
                  CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
                  F2light_d = F2light(xb)
                  FLlight_d = FLlight(xb)
                  call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
                  call SetPDFSet("EPPS21nlo_CT18Anlo_C12")
                  CALL InitializeAPFEL_DIS
                  call SetTargetDIS("proton")
                  CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
                  F2light_A = F2light(xb)
                  FLlight_A = FLlight(xb)
                  call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
                  pip2022_cdis(I)= disA/disD
            END DO
            print *, "C done"
      END IF
C----------Fe
      IF(I_pip_2022.EQ.1) THEN
            NLINES=1
            NTOT=48
            FILEEXP = RT//'JLAB2022/pi+.csv'
!            CALL READDATA20(FILEEXP,NLINES,X1,X2,X3,X4,X5,
!     *                      X6,X7, X8, X9, X10, X11, X12, X13,
!     *                     X14,X15,X16,X17,X18,X19,X20,NTOT)
            DO I=1,NTOT
		  print *,1
                  pip2022_zlow(I)   =X1(I)
                  pip2022_zup(I)    =X2(I)
                  pip2022_ptlow(I)  =X4(I)
                  pip2022_ptup(I)   =X5(I)
                  pip2022_fe(I)      =X9(I)
                  pip2022_festat(I)  =X10(I)
                  pip2022_fesys(I)   =X11(I)
                  pip2022_qfe(I)     =X16(I)
                  pip2022_xfe(I)     =X19(I)
                  print *,2
                  Qs = pip2022_qfe(I)
                  xb  = pip2022_xfe(I)
                  zh  = 0.5d0*(pip2022_zlow(I)+pip2022_zup(I))
                  pht = sqrt(0.5d0*(pip2022_ptlow(I)+pip2022_ptup(I)))
                  IH = 1
                  IC = 1
                  IT = 56
                  CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
                  CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
                  pip2022_fefuu(I) = fuu
                  pip2022_fefuua(I) = fuua
                  call SetPDFSet("CT18ANLO")
                  CALL InitializeAPFEL_DIS
                  call SetTargetDIS("isoscalar")
                  CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
                  F2light_d = F2light(xb)
                  FLlight_d = FLlight(xb)
                  call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
                  call SetPDFSet("EPPS21nlo_CT18Anlo_Fe56")
                  CALL InitializeAPFEL_DIS
                  call SetTargetDIS("proton")
                  CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
                  F2light_A = F2light(xb)
                  FLlight_A = FLlight(xb)
                  call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
                  pip2022_fedis(I)= disA/disD
            END DO
      END IF
C----------Pb
      IF(I_pip_2022.EQ.1) THEN
            NLINES=1
            NTOT=48
            FILEEXP = RT//'JLAB2022/pi+.csv'
!            CALL READDATA20(FILEEXP,NLINES,X1,X2,X3,X4,X5,
!     *                      X6,X7, X8, X9, X10, X11, X12, X13,
!     *                     X14,X15,X16,X17,X18,X19,X20,NTOT)
            DO I=1,NTOT
                  pip2022_zlow(I)   =X1(I)
                  pip2022_zup(I)    =X2(I)
                  pip2022_ptlow(I)  =X4(I)
                  pip2022_ptup(I)   =X5(I)
                  pip2022_pb(I)      =X12(I)
                  pip2022_pbstat(I)  =X13(I)
                  pip2022_pbsys(I)   =X14(I)
                  pip2022_qpb(I)     =X17(I)
                  pip2022_xpb(I)     =x20(I)
                  Qs = pip2022_qpb(I)
                  xb  = pip2022_xpb(I)
                  zh  = 0.5d0*(pip2022_zlow(I)+pip2022_zup(I))
                  pht = sqrt(0.5d0*(pip2022_ptlow(I)+pip2022_ptup(I)))
                  IH = 1
                  IC = 1
                  IT = 208
                  CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
                  CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
                  pip2022_pbfuu(I) = fuu
                  pip2022_pbfuua(I) = fuua
                  call SetPDFSet("CT18ANLO")
                  CALL InitializeAPFEL_DIS
                  call SetTargetDIS("isoscalar")
                  CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
                  F2light_d = F2light(xb)
                  FLlight_d = FLlight(xb)
                  call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
                  call SetPDFSet("EPPS21nlo_CT18Anlo_Pb208")
                  CALL InitializeAPFEL_DIS
                  call SetTargetDIS("proton")
                  CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
                  F2light_A = F2light(xb)
                  FLlight_A = FLlight(xb)
                  call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
                  pip2022_pbdis(I)= disA/disD
            END DO
      END IF

C---------pi-
C----------C
      IF(I_pim_2022.EQ.1) THEN
            NLINES=1
            NTOT=40
            FILEEXP = RT//'JLAB2022/pi-.csv'
            CALL READDATA20(FILEEXP,NLINES,X1,X2,X3,X4,X5,
     *                      X6,X7, X8, X9, X10, X11, X12, X13,
     *                     X14,X15,X16,X17,X18,X19,X20,NTOT)
            DO I=1,NTOT
                  pim2022_zlow(I)   =X1(I)
                  pim2022_zup(I)    =X2(I)
                  pim2022_ptlow(I)  =X4(I)
                  pim2022_ptup(I)   =X5(I)
                  pim2022_c(I)      =X6(I)
                  pim2022_cstat(I)  =X7(I)
                  pim2022_csys(I)   =X8(I)
                  pim2022_qc(I)     =X15(I)
                  pim2022_xc(I)     =x18(I)
                  Qs = pim2022_qc(I)
                  xb  = pim2022_xc(I)
                  zh  = 0.5d0*(pim2022_zlow(I)+pim2022_zup(I))
                  pht = sqrt(0.5d0*(pim2022_ptlow(I)+pim2022_ptup(I)))
                  IH = 1
                  IC = -1
                  IT = 12
                  CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
                  CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
                  pim2022_cfuu(I) = fuu
                  pim2022_cfuua(I) = fuua
                  call SetPDFSet("CT18ANLO")
                  CALL InitializeAPFEL_DIS
                  call SetTargetDIS("isoscalar")
                  CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
                  F2light_d = F2light(xb)
                  FLlight_d = FLlight(xb)
                  call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
                  call SetPDFSet("EPPS21nlo_CT18Anlo_C12")
                  CALL InitializeAPFEL_DIS
                  call SetTargetDIS("proton")
                  CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
                  F2light_A = F2light(xb)
                  FLlight_A = FLlight(xb)
                  call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
                  pim2022_cdis(I)= disA/disD
            END DO
      END IF
C----------Fe
      IF(I_pim_2022.EQ.1) THEN
	    print *, "Fe"
            NLINES=1
            NTOT=40
            FILEEXP = RT//'JLAB2022/pi-.csv'
!            CALL READDATA20(FILEEXP,NLINES,X1,X2,X3,X4,X5,
!     *                      X6,X7, X8, X9, X10, X11, X12, X13,
!     *                     X14,X15,X16,X17,X18,X19,X20,NTOT)
            DO I=1,NTOT
                  pim2022_zlow(I)   =X1(I)
                  pim2022_zup(I)    =X2(I)
                  pim2022_ptlow(I)  =X4(I)
                  pim2022_ptup(I)   =X5(I)
                  pim2022_fe(I)      =X9(I)
                  pim2022_festat(I)  =X10(I)
                  pim2022_fesys(I)   =X11(I)
                  pim2022_qfe(I)     =X16(I)
                  pim2022_xfe(I)     =x19(I)
                  Qs = pim2022_qfe(I)
                  xb  = pim2022_xfe(I)
                  zh  = 0.5d0*(pim2022_zlow(I)+pim2022_zup(I))
                  pht = sqrt(0.5d0*(pim2022_ptlow(I)+pim2022_ptup(I)))
                  IH = 1
                  IC = -1
                  IT = 56
                  CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
                  CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
                  pim2022_fefuu(I) = fuu
                  pim2022_fefuua(I) = fuua
                  call SetPDFSet("CT18ANLO")
                  CALL InitializeAPFEL_DIS
                  call SetTargetDIS("isoscalar")
                  CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
                  F2light_d = F2light(xb)
                  FLlight_d = FLlight(xb)
                  call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
                  call SetPDFSet("EPPS21nlo_CT18Anlo_Fe56")
                  CALL InitializeAPFEL_DIS
                  call SetTargetDIS("proton")
                  CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
                  F2light_A = F2light(xb)
                  FLlight_A = FLlight(xb)
                  call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
                  pim2022_fedis(I)= disA/disD
            END DO
      END IF
C----------Pb
      IF(I_pim_2022.EQ.1) THEN
            NLINES=1
            NTOT=40
            FILEEXP = RT//'JLAB2022/pi-.csv'
!            CALL READDATA20(FILEEXP,NLINES,X1,X2,X3,X4,X5,
!     *                      X6,X7, X8, X9, X10, X11, X12, X13,
!     *                     X14,X15,X16,X17,X18,X19,X20,NTOT)
            DO I=1,NTOT
                  pim2022_zlow(I)   =X1(I)
                  pim2022_zup(I)    =X2(I)
                  pim2022_ptlow(I)  =X4(I)
                  pim2022_ptup(I)   =X5(I)
                  pim2022_pb(I)      =X12(I)
                  pim2022_pbstat(I)  =X13(I)
                  pim2022_pbsys(I)   =X14(I)
                  pim2022_qpb(I)     =X17(I)
                  pim2022_xpb(I)     =x20(I)
                  Qs = pim2022_qpb(I)
                  xb  = pim2022_xpb(I)
                  zh  = 0.5d0*(pim2022_zlow(I)+pim2022_zup(I))
                  pht = sqrt(0.5d0*(pim2022_ptlow(I)+pim2022_ptup(I)))
                  IH = 1
                  IC = -1
                  IT = 208
                  CALL DISUU(s,Qs,xb,zh,pht,fuua,IT,IH,IC)
                  CALL DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
                  pim2022_pbfuu(I) = fuu
                  pim2022_pbfuua(I) = fuua
                  call SetPDFSet("CT18ANLO")
                  CALL InitializeAPFEL_DIS
                  call SetTargetDIS("isoscalar")
                  CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
                  F2light_d = F2light(xb)
                  FLlight_d = FLlight(xb)
                  call dis_cx(s, xb, Qs, F2light_d, FLlight_d, disD)
                  call SetPDFSet("EPPS21nlo_CT18Anlo_Pb208")
                  CALL InitializeAPFEL_DIS
                  call SetTargetDIS("proton")
                  CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
                  F2light_A = F2light(xb)
                  FLlight_A = FLlight(xb)
                  call dis_cx(s, xb, Qs, F2light_A, FLlight_A, disA)
                  pim2022_pbdis(I)= disA/disD
            END DO
      END IF

      END SUBROUTINE READDATA


      SUBROUTINE DISUUD(s,Qs,xb,zh,pht,fuu,IT,IH,IC)
      IMPLICIT NONE
              REAL*8 s, Qs, Nu, xb, zh, pht, fuu
              INTEGER IT, IH, IC, ID

              ID = 3
              CALL DISUU(s,Qs,xb,zh,pht,fuu,ID,IH,IC)
      END SUBROUTINE DISUUD

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
      DO
      NT = NT + 1
      IF (NT > NTOT) EXIT
      READ(91,*,END=99, ERR=100)
     +     X1(NT),X2(NT),X3(NT),X4(NT),X5(NT),X6(NT),
     +     X7(NT),X8(NT),X9(NT),X10(NT),X11(NT),X12(NT),X13(NT),
     +     X14(NT)
      END DO
        
99    CONTINUE
      CLOSE(91)
      RETURN
        
100   PRINT *, 'ERROR IN READING FILE    ', FILEEXP
      STOP
      END

c      SUBROUTINE READDATA14(FILEEXP,NLINES,X1,
c     $           X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12, X13,X14,NTOT)
c      CHARACTER*72 FILEEXP,LINE
c      REAL*8 X1(100),X2(100),X3(100),X4(100),X5(100),X6(100)
c      REAL*8 X7(100),X8(100),X9(100),X10(100),X11(100),X12(100),X13(100)
c      REAL*8 X14(100)
c      INTEGER NLINES,NTOT,I,NT

c      OPEN(unit = 91,FILE=FILEEXP,STATUS='unknown')

c      DO I = 1,NLINES,1
c       READ(91,'(A)') LINE    ! skip NLINES..
c      ENDDO

c      NT = 0
c   10   NT =  NT +1

c      READ(91,*,ERR = 99, end = 11)
c     +     X1(NT),X2(NT),X3(NT),X4(NT),X5(NT),X6(NT),
c     +     X7(NT),X8(NT),X9(NT),X10(NT),X11(NT),X12(NT),X13(NT),
c     +     x14(NT)

c      goto 10
c   11 if(NT.GT.NTOT) goto 101

c      CLOSE(91)

c      GOTO 101
c 99    PRINT *, 'ERROR IN READING FILE    ', FILEEXP
c 101   CONTINUE
c      RETURN
c      END

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
c     read data file which has 7 columns
c-----------------------------------------
      SUBROUTINE READDATA7(FILEEXP,NLINES,
     >                    X1,X2,X3,X4,X5,X6,X7,NTOT)
      CHARACTER*72 FILEEXP,LINE
      REAL*8 X1(100),X2(100),X3(100),X4(100),X5(100),X6(100)
      REAL*8 X7(100)
      INTEGER NLINES,NTOT,I,NT

      OPEN(unit = 91,FILE=FILEEXP,STATUS='unknown')

      DO I = 1,NLINES,1
         READ(91,'(A)') LINE    ! skip NLINES..
      ENDDO

      NT = 0
 10   NT =  NT +1

      READ(91,*,ERR = 99, end = 11)
     +     X1(NT),X2(NT),X3(NT),X4(NT),X5(NT),X6(NT),
     +     X7(NT)

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
     >                    X1,X2,X3,X4,X5,NTOT)
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
c     read data file which has 4 columns
c-----------------------------------------
      SUBROUTINE READDATA4(FILEEXP,NLINES,
     >                    X1,X2,X3,X4,NTOT)
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
