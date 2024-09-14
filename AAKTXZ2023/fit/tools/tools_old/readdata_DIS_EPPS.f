      SUBROUTINE READDATA
      IMPLICIT NONE
C-----------STORED VARIABLE DECLARATIONS
C-----------DECLARE VARIABLES.
C-----------NLINES = NUMBER OF LINES SKIPPED.
C-----------NTOT = TOTAL NUMBER OF ROWS.
      REAL*8 X1(200),X2(200),X3(200),X4(200),X5(200),X6(200)
      REAL*8 X7(200),X8(200),X9(200),X10(200),X11(200),X12(200),X13(200)
      REAL*8 Sep_hermes
      REAL*8 rts, Qs, Nu, xb, zh, pht, fuu, fuua
      REAL*8 rts_Jlab
      REAL*8 Sep_Jlab
      real*8 Sep_EIC
      double precision M
      CHARACTER*72 FILEEXP
      INTEGER I,J
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

C-----------IMPORTATION OF ALL HERMES (HERA-DESY) SIDIS DATA
      NTOT = 9 ! 9 ROWS
      NTOT2 = 8! 8 ROWS
      NTOT3 = 22
      NLINES = 1 !SKIP THE FIRST LINE

C------------PARAMETERS
      Q0 = 1.3d0
      M = 0.938272046d0 ! GeV (Mass of Proton)
      Sep_hermes = 2*27.6*M + M**2 ! GeV^2
      Sep_JLAB   = 2*12*M + M**2 ! GeV^2
      Sep_EIC    = 2d0*2d0*10d0*110d0 ! GeV^2


      IF(I_EIC_PRED.EQ.1) THEN
            DO I =1,3
            rts = Sep_EIC
            IF(I.EQ.1) THEN
              Qs = 4d0
            ELSEIF(I.EQ.2) THEN
              Qs = 25d0
            ELSEIF(I.EQ.3) THEN
              Qs = 100d0
            ENDIF
            Q = dsqrt(Qs)
            xb = 0.05
            ! Compute DIS CX for proton
            call SetPDFSet("CT14nlo")
            call SetTargetDIS("isoscalar")
            call InitializeAPFEL_DIS
            CALL ComputeStructureFunctionsAPFEL(Q0,Q)
            F2light_D = F2light(xb)
            FLlight_D = FLlight(xb)
            call dis_cx(rts,xb,Qs,F2light_D,FLlight_D,disD)
            ! Compute DIS CX for Au
            call SetPDFSet("EPPS16nlo_CT14nlo_Au197")
            call SetTargetDIS("proton")
            call InitializeAPFEL_DIS
            CALL ComputeStructureFunctionsAPFEL(Q0,Q)
            F2light_A = F2light(xb)
            FLlight_A = FLlight(xb)
            call dis_cx(rts,xb,Qs,F2light_A,FLlight_A,disA)
            write(1001,*) "Qs DIS"
            write(1001,*)  Qs, disA/disD
            ENDDO
      ENDIF

      IF(I_JLAB_PRED.EQ.1) THEN
            rts = Sep_JLAB
            Qs = 2.5
            xb = 0.4
            Q = dsqrt(Qs)
            ! Compute DIS CX for proton
            call SetPDFSet("CT14nlo")
            call SetTargetDIS("isoscalar")
            call InitializeAPFEL_DIS
            CALL ComputeStructureFunctionsAPFEL(Q0,Q)
            F2light_D = F2light(xb)
            FLlight_D = FLlight(xb)
            call dis_cx(rts,xb,Qs,F2light_D,FLlight_D,disD)
            ! Compute DIS CX for Au
            call SetPDFSet("EPPS16nlo_CT14nlo_Pb208")
            call SetTargetDIS("proton")
            call InitializeAPFEL_DIS
            CALL ComputeStructureFunctionsAPFEL(Q0,Q)
            F2light_A = F2light(xb)
            FLlight_A = FLlight(xb)
            call dis_cx(rts,xb,Qs,F2light_A,FLlight_A,disA)
            write(1002,*) "Qs DIS"
            write(1002,*)  Qs, disA/disD
      ENDIF



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
               call SetPDFSet("EPPS16nlo_CT14nlo_C12")
               ELSEIF(X1(I).eq.1) THEN
               IT = 56
               call SetPDFSet("EPPS16nlo_CT14nlo_Fe56")
               ELSEIF(X1(I).eq.2) THEN
               IT = 208
               call SetPDFSet("EPPS16nlo_CT14nlo_Pb208")
               ENDIF
               ! Computing DIS CX for e+A
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_A, FLlight_A, disA)
               !Computing DIS CX for e+d
               call SetPDFSet("CT14nlo")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_d, FLlight_d, disD)
               ! Computing SIDIS CX for e+A
               CALL DISUU(rts,Qs,xb,zh,pht,fuua,IT,IH,IC)
               ! Computing SIDIS CX for d+A
               CALL DISUUD(rts,Qs,xb,zh,pht,fuu,IT,IH,IC)
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
               call SetPDFSet("EPPS16nlo_CT14nlo_C12")
               ELSEIF(X1(I).eq.1) THEN
               IT = 56
               call SetPDFSet("EPPS16nlo_CT14nlo_Fe56")
               ELSEIF(X1(I).eq.2) THEN
               IT = 208
               call SetPDFSet("EPPS16nlo_CT14nlo_Pb208")
               ENDIF
               ! Computing DIS CX for e+A
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_A, FLlight_A, disA)
               !Computing DIS CX for e+d
               call SetPDFSet("CT14nlo")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_d, FLlight_d, disD)
               ! Computing SIDIS CX for e+A
               CALL DISUU(rts,Qs,xb,zh,pht,fuua,IT,IH,IC)
               ! Computing SIDIS CX for d+A
               CALL DISUUD(rts,Qs,xb,zh,pht,fuu,IT,IH,IC)
               JLAB_PIM_FUU(I) = fuu
               JLAB_PIM_FUUA(I) = fuua
               JLAB_PIM_DIS(I)= disA/disD
            END DO
      END IF



      rts = SEP_HERMES

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
               CALL DISUU(rts,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(rts,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIP_HE_Z_FUU(I) = fuu
               PIP_HE_Z_FUUA(I) = fuua
               call SetPDFSet("CT14nlo")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_He4")
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_A, FLlight_A, disA)
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
               CALL DISUU(rts,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(rts,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIP_NE_Z_FUU(I) = fuu
               PIP_NE_Z_FUUA(I) = fuua
               call SetPDFSet("CT14nlo")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_Ne20")
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_A, FLlight_A, disA)
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
               CALL DISUU(rts,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(rts,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIP_KR_Z_FUU(I) = fuu
               PIP_KR_Z_FUUA(I) = fuua
               call SetPDFSet("CT14nlo")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_Kr84")
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_A, FLlight_A, disA)
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
               CALL DISUU(rts,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(rts,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIP_XE_Z_FUU(I) = fuu
               PIP_XE_Z_FUUA(I) = fuua
               call SetPDFSet("CT14nlo")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_Xe131")
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_A, FLlight_A, disA)
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
               CALL DISUU(rts,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(rts,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIM_HE_Z_FUU(I) = fuu
               PIM_HE_Z_FUUA(I) = fuua
               call SetPDFSet("CT14nlo")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_He4")
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_A, FLlight_A, disA)
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
               CALL DISUU(rts,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(rts,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIM_NE_Z_FUU(I) = fuu
               PIM_NE_Z_FUUA(I) = fuua
               call SetPDFSet("CT14nlo")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_Ne20")
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_A, FLlight_A, disA)
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
               CALL DISUU(rts,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(rts,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIM_KR_Z_FUU(I) = fuu
               PIM_KR_Z_FUUA(I) = fuua
               call SetPDFSet("CT14nlo")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_Kr84")
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_A, FLlight_A, disA)
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
               CALL DISUU(rts,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(rts,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIM_XE_Z_FUU(I) = fuu
               PIM_XE_Z_FUUA(I) = fuua
               call SetPDFSet("CT14nlo")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_Xe131")
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_A, FLlight_A, disA)
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
               CALL DISUU(rts,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(rts,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PI0_HE_Z_FUU(I) = fuu
               PI0_HE_Z_FUUA(I) = fuua
               call SetPDFSet("CT14nlo")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_He4")
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_A, FLlight_A, disA)
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
               CALL DISUU(rts,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(rts,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PI0_NE_Z_FUU(I) = fuu
               PI0_NE_Z_FUUA(I) = fuua
               call SetPDFSet("CT14nlo")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_Ne20")
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_A, FLlight_A, disA)
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
               CALL DISUU(rts,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(rts,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PI0_KR_Z_FUU(I) = fuu
               PI0_KR_Z_FUUA(I) = fuua
               call SetPDFSet("CT14nlo")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_Kr84")
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_A, FLlight_A, disA)
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
               CALL DISUU(rts,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(rts,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PI0_XE_Z_FUU(I) = fuu
               PI0_XE_Z_FUUA(I) = fuua
               call SetPDFSet("CT14nlo")
               call SetTargetDIS("isoscalar")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_Xe131")
               call SetTargetDIS("proton")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_A, FLlight_A, disA)
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
               CALL DISUU(rts,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(rts,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PI0_HE_PT2_FUU(I) = fuu
               PI0_HE_PT2_FUUA(I) = fuua
               call SetPDFSet("CT14nlo")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("isoscalar")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_He4")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("proton")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_A, FLlight_A, disA)
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
               CALL DISUU(rts,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(rts,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PI0_NE_PT2_FUU(I) = fuu
               PI0_NE_PT2_FUUA(I) = fuua
               call SetPDFSet("CT14nlo")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("isoscalar")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_Ne20")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("proton")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_A, FLlight_A, disA)
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
               CALL DISUU(rts,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(rts,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PI0_KR_PT2_FUU(I) = fuu
               PI0_KR_PT2_FUUA(I) = fuua
               call SetPDFSet("CT14nlo")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("isoscalar")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_Kr84")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("proton")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_A, FLlight_A, disA)
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
               CALL DISUU(rts,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(rts,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PI0_XE_PT2_FUU(I) = fuu
               PI0_XE_PT2_FUUA(I) = fuua
               call SetPDFSet("CT14nlo")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("isoscalar")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_Xe131")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("proton")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_A, FLlight_A, disA)
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
               CALL DISUU(rts,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(rts,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIP_HE_PT2_FUU(I) = fuu
               PIP_HE_PT2_FUUA(I) = fuua
               call SetPDFSet("CT14nlo")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("isoscalar")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_He4")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("proton")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_A, FLlight_A, disA)
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
               CALL DISUU(rts,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(rts,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIP_NE_PT2_FUU(I) = fuu
               PIP_NE_PT2_FUUA(I) = fuua
               call SetPDFSet("CT14nlo")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("isoscalar")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_Ne20")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("proton")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_A, FLlight_A, disA)
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
               CALL DISUU(rts,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(rts,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIP_KR_PT2_FUU(I) = fuu
               PIP_KR_PT2_FUUA(I) = fuua
               call SetPDFSet("CT14nlo")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("isoscalar")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_Kr84")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("proton")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_A, FLlight_A, disA)
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
               CALL DISUU(rts,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(rts,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIP_XE_PT2_FUU(I) = fuu
               PIP_XE_PT2_FUUA(I) = fuua
               call SetPDFSet("CT14nlo")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("isoscalar")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_Xe131")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("proton")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_A, FLlight_A, disA)
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
               CALL DISUU(rts,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(rts,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIM_HE_PT2_FUU(I) = fuu
               PIM_HE_PT2_FUUA(I) = fuua
               call SetPDFSet("CT14nlo")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("isoscalar")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_He4")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("proton")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_A, FLlight_A, disA)
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
               CALL DISUU(rts,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(rts,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIM_NE_PT2_FUU(I) = fuu
               PIM_NE_PT2_FUUA(I) = fuua
               call SetPDFSet("CT14nlo")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("isoscalar")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_Ne20")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("proton")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_A, FLlight_A, disA)
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
               CALL DISUU(rts,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(rts,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIM_KR_PT2_FUU(I) = fuu
               PIM_KR_PT2_FUUA(I) = fuua
               call SetPDFSet("CT14nlo")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("isoscalar")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_Kr84")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("proton")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_A, FLlight_A, disA)
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
               CALL DISUU(rts,Qs,xb,zh,pht,fuua,IT,IH,IC)
               CALL DISUUD(rts,Qs,xb,zh,pht,fuu,IT,IH,IC)
               PIM_XE_PT2_FUU(I) = fuu
               PIM_XE_PT2_FUUA(I) = fuua
               call SetPDFSet("CT14nlo")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("isoscalar")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_Xe131")
               CALL InitializeAPFEL_DIS
               call SetTargetDIS("proton")
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, F2light_A, FLlight_A, disA)
               PIM_XE_PT2_DIS(I)= disA/disD
            END DO
      END IF



      END SUBROUTINE READDATA


      SUBROUTINE DISUUD(rts,Qs,xb,zh,pht,fuu,IT,IH,IC)
      IMPLICIT NONE
              REAL*8 rts, Qs, Nu, xb, zh, pht, fuu
              INTEGER IT, IH, IC, ID

              ID = 3
              CALL DISUU(rts,Qs,xb,zh,pht,fuu,ID,IH,IC)
      END SUBROUTINE DISUUD
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
