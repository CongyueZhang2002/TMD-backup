      SUBROUTINE READDATA
      IMPLICIT NONE
C-----------STORED VARIABLE DECLARATIONS
C-----------DECLARE VARIABLES.
C-----------NLINES = NUMBER OF LINES SKIPPED.
C-----------NTOT = TOTAL NUMBER OF ROWS.
      REAL*8 X1(100),X2(100),X3(100),X4(100),X5(100),X6(100)
      REAL*8 X7(100),X8(100),X9(100),X10(100),X11(100),X12(100),X13(100)
      REAL*8 Sep_hermes
      REAL*8 rts, Qs, Nu, xb, zh, pht, fuu, fuua
      REAL*8 rts_Jlab
      double precision M
      CHARACTER*72 FILEEXP
      INTEGER I,J
      INTEGER NLINES ! NUMBER OF LINES SKIPPED
      INTEGER NTOT, NTOT2, NTOT3 ! NUMBER OF ROWS
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
      Sep_hermes = 2*27.6*M + M**2 ! GeV

C-----------Initialize PDF Sets for HERMES
C----------- Deuteron (free nucleon)
      CALL InitPDFsetByNameM(1,"CT14nlo")
      CALL InitPDFM(1,0)

C-----------Helium (bound proton)
      CALL InitPDFsetByNameM(4,"EPPS16nlo_CT14nlo_He4")
      CALL InitPDFM(4,0)

C-----------Neon (bound proton)
      CALL InitPDFsetByNameM(20,"EPPS16nlo_CT14nlo_Ne20")
      CALL InitPDFM(20,0)


C-----------Krypton (bound proton)
      CALL InitPDFsetByNameM(84,"EPPS16nlo_CT14nlo_Kr84")
      CALL InitPDFM(84,0)

C-----------Xenon (bound proton)
      CALL InitPDFsetByNameM(131,"EPPS16nlo_CT14nlo_Xe131")
      CALL InitPDFM(131,0)

C-------- SETS for JLAB
C-----------Carbon (bound proton)
      CALL InitPDFsetByNameM(12,"EPPS16nlo_CT14nlo_C12")
      CALL InitPDFM(12,0)

C-----------Iron (bound proton)
      CALL InitPDFsetByNameM(56,"EPPS16nlo_CT14nlo_Fe56")
      CALL InitPDFM(56,0)

C-----------Lead (bound proton)
      CALL InitPDFsetByNameM(208,"EPPS16nlo_CT14nlo_Pb208")
      CALL InitPDFM(208,0)

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
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, Nu, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_He4")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, Nu, F2light_A, FLlight_A, disA)
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
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, Nu, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_Ne20")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, Nu, F2light_A, FLlight_A, disA)
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
               PIP_KR_Z_NU(I)       =X2(I)
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
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, Nu, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_Kr84")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, Nu, F2light_A, FLlight_A, disA)
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
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, Nu, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_Xe131")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, Nu, F2light_A, FLlight_A, disA)
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
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, Nu, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_He4")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, Nu, F2light_A, FLlight_A, disA)
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
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, Nu, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_Ne20")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, Nu, F2light_A, FLlight_A, disA)
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
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, Nu, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_Kr84")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, Nu, F2light_A, FLlight_A, disA)
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
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, Nu, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_Xe131")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, Nu, F2light_A, FLlight_A, disA)
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
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, Nu, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_He4")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, Nu, F2light_A, FLlight_A, disA)
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
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, Nu, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_Ne20")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, Nu, F2light_A, FLlight_A, disA)
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
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, Nu, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_Kr84")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, Nu, F2light_A, FLlight_A, disA)
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
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_d = F2light(xb)
               FLlight_d = FLlight(xb)
               call dis_cx(rts, xb, Qs, Nu, F2light_d, FLlight_d, disD)
               call SetPDFSet("EPPS16nlo_CT14nlo_Xe131")
               CALL InitializeAPFEL_DIS
               CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
               F2light_A = F2light(xb)
               FLlight_A = FLlight(xb)
               call dis_cx(rts, xb, Qs, Nu, F2light_A, FLlight_A, disA)
               PI0_XE_Z_DIS(I)= disA/disD
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
