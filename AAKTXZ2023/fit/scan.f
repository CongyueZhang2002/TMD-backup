      PROGRAM MASTER
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      integer i,j,k,l
      real*8 chi2, chisquare
      character*12 filenames(101)
      real*8 chi2s(101)
      integer nums(101)
      real*8 XX(6)

C-----HERMES
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

C-----RHIC (DY)
      real*8 CHI2_RHIC_Ratio_pAu1, CHI2_RHIC_Ratio_pAu2
      integer NUM_RHIC_Ratio_pAu1, NUM_RHIC_Ratio_pAu2

C-----LHC (DY)
      real*8 CHI2_ATLAS5_Y1,CHI2_ATLAS5_Y2,CHI2_ATLAS5_Y3
      real*8 CHI2_CMS5
      integer NUM_ATLAS5_Y1,NUM_ATLAS5_Y2,NUM_ATLAS5_Y3
      integer NUM_CMS5

C-----E866
      real*8  CHI2_E866_800q, CHI2_E772_800
      integer NUM_E866_800q , num_e772_800

      common /hermesch2/
     >  CHI2_HERMES,
     >  CHI2_PIP_HE_NU,
     >  CHI2_PIP_HE_Z,
     >  CHI2_PIP_HE_Q2,
     >  CHI2_PIP_NE_NU,
     >  CHI2_PIP_NE_Z,
     >  CHI2_PIP_NE_Q2,
     >  CHI2_PIP_KR_NU,
     >  CHI2_PIP_KR_Z,
     >  CHI2_PIP_KR_Q2,
     >  CHI2_PIP_XE_NU,
     >  CHI2_PIP_XE_Z,
     >  CHI2_PIP_XE_Q2,
     >  CHI2_PIM_HE_NU,
     >  CHI2_PIM_HE_Z,
     >  CHI2_PIM_HE_Q2,
     >  CHI2_PIM_NE_NU,
     >  CHI2_PIM_NE_Z,
     >  CHI2_PIM_NE_Q2,
     >  CHI2_PIM_KR_NU,
     >  CHI2_PIM_KR_Z,
     >  CHI2_PIM_KR_Q2,
     >  CHI2_PIM_XE_NU,
     >  CHI2_PIM_XE_Z,
     >  CHI2_PIM_XE_Q2,
     >  CHI2_PI0_HE_NU,
     >  CHI2_PI0_HE_Z,
     >  CHI2_PI0_HE_Q2,
     >  CHI2_PI0_NE_NU,
     >  CHI2_PI0_NE_Z,
     >  CHI2_PI0_NE_Q2,
     >  CHI2_PI0_KR_NU,
     >  CHI2_PI0_KR_Z,
     >  CHI2_PI0_KR_Q2,
     >  CHI2_PI0_XE_NU,
     >  CHI2_PI0_XE_Z,
     >  CHI2_PI0_XE_Q2,
     >  CHI2_PIP_HE_PT2,
     >  CHI2_PIM_HE_PT2,
     >  CHI2_PI0_HE_PT2,
     >  CHI2_PIP_NE_PT2,
     >  CHI2_PIM_NE_PT2,
     >  CHI2_PI0_NE_PT2,
     >  CHI2_PIP_KR_PT2,
     >  CHI2_PIM_KR_PT2,
     >  CHI2_PI0_KR_PT2,
     >  CHI2_PIP_XE_PT2,
     >  CHI2_PIM_XE_PT2,
     >  CHI2_PI0_XE_PT2,
     >  CHI2_KP_HE_NU,
     >  CHI2_KP_HE_Z,
     >  CHI2_KP_HE_Q2,
     >  CHI2_KP_NE_NU,
     >  CHI2_KP_NE_Z,
     >  CHI2_KP_NE_Q2,
     >  CHI2_KP_KR_NU,
     >  CHI2_KP_KR_Z,
     >  CHI2_KP_KR_Q2,
     >  CHI2_KP_XE_NU,
     >  CHI2_KP_XE_Z,
     >  CHI2_KP_XE_Q2,
     >  CHI2_KM_HE_NU,
     >  CHI2_KM_HE_Z,
     >  CHI2_KM_HE_Q2,
     >  CHI2_KM_NE_NU,
     >  CHI2_KM_NE_Z,
     >  CHI2_KM_NE_Q2,
     >  CHI2_KM_KR_NU,
     >  CHI2_KM_KR_Z,
     >  CHI2_KM_KR_Q2,
     >  CHI2_KM_XE_NU,
     >  CHI2_KM_XE_Z,
     >  CHI2_KM_XE_Q2,
     >  CHI2_KP_HE_PT2,
     >  CHI2_KM_HE_PT2,
     >  CHI2_KP_NE_PT2,
     >  CHI2_KM_NE_PT2,
     >  CHI2_KP_KR_PT2,
     >  CHI2_KM_KR_PT2,
     >  CHI2_KP_XE_PT2,
     >  CHI2_KM_XE_PT2,
     >  CHI2_PIP_NU,
     >  CHI2_PIP_Z,
     >  CHI2_PIP_Q2,
     >  CHI2_PIP_PT2,
     >  CHI2_PI0_NU,
     >  CHI2_PI0_Z,
     >  CHI2_PI0_Q2,
     >  CHI2_PI0_PT2,
     >  CHI2_PIM_NU,
     >  CHI2_PIM_Z,
     >  CHI2_PIM_Q2,
     >  CHI2_PIM_PT2,
     >  NUM_HERMES,
     >  NUM_PIP_HE_NU,
     >  NUM_PIP_HE_Z,
     >  NUM_PIP_HE_Q2,
     >  NUM_PIP_NE_NU,
     >  NUM_PIP_NE_Z,
     >  NUM_PIP_NE_Q2,
     >  NUM_PIP_KR_NU,
     >  NUM_PIP_KR_Z,
     >  NUM_PIP_KR_Q2,
     >  NUM_PIP_XE_NU,
     >  NUM_PIP_XE_Z,
     >  NUM_PIP_XE_Q2,
     >  NUM_PIM_HE_NU,
     >  NUM_PIM_HE_Z,
     >  NUM_PIM_HE_Q2,
     >  NUM_PIM_NE_NU,
     >  NUM_PIM_NE_Z,
     >  NUM_PIM_NE_Q2,
     >  NUM_PIM_KR_NU,
     >  NUM_PIM_KR_Z,
     >  NUM_PIM_KR_Q2,
     >  NUM_PIM_XE_NU,
     >  NUM_PIM_XE_Z,
     >  NUM_PIM_XE_Q2,
     >  NUM_PI0_HE_NU,
     >  NUM_PI0_HE_Z,
     >  NUM_PI0_HE_Q2,
     >  NUM_PI0_NE_NU,
     >  NUM_PI0_NE_Z,
     >  NUM_PI0_NE_Q2,
     >  NUM_PI0_KR_NU,
     >  NUM_PI0_KR_Z,
     >  NUM_PI0_KR_Q2,
     >  NUM_PI0_XE_NU,
     >  NUM_PI0_XE_Z,
     >  NUM_PI0_XE_Q2,
     >  NUM_PIP_HE_PT2,
     >  NUM_PIM_HE_PT2,
     >  NUM_PI0_HE_PT2,
     >  NUM_PIP_NE_PT2,
     >  NUM_PIM_NE_PT2,
     >  NUM_PI0_NE_PT2,
     >  NUM_PIP_KR_PT2,
     >  NUM_PIM_KR_PT2,
     >  NUM_PI0_KR_PT2,
     >  NUM_PIP_XE_PT2,
     >  NUM_PIM_XE_PT2,
     >  NUM_PI0_XE_PT2,
     >  NUM_KP_HE_NU,
     >  NUM_KP_HE_Z,
     >  NUM_KP_HE_Q2,
     >  NUM_KP_NE_NU,
     >  NUM_KP_NE_Z,
     >  NUM_KP_NE_Q2,
     >  NUM_KP_KR_NU,
     >  NUM_KP_KR_Z,
     >  NUM_KP_KR_Q2,
     >  NUM_KP_XE_NU,
     >  NUM_KP_XE_Z,
     >  NUM_KP_XE_Q2,
     >  NUM_KM_HE_NU,
     >  NUM_KM_HE_Z,
     >  NUM_KM_HE_Q2,
     >  NUM_KM_NE_NU,
     >  NUM_KM_NE_Z,
     >  NUM_KM_NE_Q2,
     >  NUM_KM_KR_NU,
     >  NUM_KM_KR_Z,
     >  NUM_KM_KR_Q2,
     >  NUM_KM_XE_NU,
     >  NUM_KM_XE_Z,
     >  NUM_KM_XE_Q2,
     >  NUM_KP_HE_PT2,
     >  NUM_KM_HE_PT2,
     >  NUM_KP_NE_PT2,
     >  NUM_KM_NE_PT2,
     >  NUM_KP_KR_PT2,
     >  NUM_KM_KR_PT2,
     >  NUM_KP_XE_PT2,
     >  NUM_KM_XE_PT2,
     >  NUM_PIP_NU,
     >  NUM_PIP_Z,
     >  NUM_PIP_Q2,
     >  NUM_PIP_PT2,
     >  NUM_PI0_NU,
     >  NUM_PI0_Z,
     >  NUM_PI0_Q2,
     >  NUM_PI0_PT2,
     >  NUM_PIM_NU,
     >  NUM_PIM_Z,
     >  NUM_PIM_Q2,
     >  NUM_PIM_PT2

      common /RHICch2/
     > CHI2_RHIC_Ratio_pAu1,
     > CHI2_RHIC_Ratio_pAu2,
     > NUM_RHIC_Ratio_pAu1,
     > NUM_RHIC_Ratio_pAu2

      common /LHCch2/
     > CHI2_ATLAS5_Y1,
     > CHI2_ATLAS5_Y2,
     > CHI2_ATLAS5_Y3,
     > CHI2_CMS5,
     > NUM_ATLAS5_Y1,
     > NUM_ATLAS5_Y2,
     > NUM_ATLAS5_Y3,
     > NUM_CMS5

      common /E866ch2/
     > CHI2_E866_800q,
     > CHI2_E772_800,
     > NUM_E866_800q,
     > num_e772_800

C-----Local Variable Declarations
      include "tools/data-inc.f"
      integer nloops,hop,nll,prepdf,preff
      common /scheme/ nloops,hop,nll,prepdf,preff


c-----FIT PARAMETERS
      nloops = 2 ! LO: 1 NLO: 2
      nll = 3
      prepdf = 1
      preff  = 1

!     CALL SETLHAPARM('SILENT') ! TO NOT SHOW THE CALLS, ALTHOUGH THEY ARE CALLED
      CALL SETLHAPARM('SILENT')

      filenames = (/
     > "scnn/001.dat",
     > "scnn/002.dat",
     > "scnn/003.dat",
     > "scnn/004.dat",
     > "scnn/005.dat",
     > "scnn/006.dat",
     > "scnn/007.dat",
     > "scnn/008.dat",
     > "scnn/009.dat",
     > "scnn/010.dat",
     > "scnn/011.dat",
     > "scnn/012.dat",
     > "scnn/013.dat",
     > "scnn/014.dat",
     > "scnn/015.dat",
     > "scnn/016.dat",
     > "scnn/017.dat",
     > "scnn/018.dat",
     > "scnn/019.dat",
     > "scnn/020.dat",
     > "scnn/021.dat",
     > "scnn/022.dat",
     > "scnn/023.dat",
     > "scnn/024.dat",
     > "scnn/025.dat",
     > "scnn/026.dat",
     > "scnn/027.dat",
     > "scnn/028.dat",
     > "scnn/029.dat",
     > "scnn/030.dat",
     > "scnn/031.dat",
     > "scnn/032.dat",
     > "scnn/033.dat",
     > "scnn/034.dat",
     > "scnn/035.dat",
     > "scnn/036.dat",
     > "scnn/037.dat",
     > "scnn/038.dat",
     > "scnn/039.dat",
     > "scnn/040.dat",
     > "scnn/041.dat",
     > "scnn/042.dat",
     > "scnn/043.dat",
     > "scnn/044.dat",
     > "scnn/045.dat",
     > "scnn/046.dat",
     > "scnn/047.dat",
     > "scnn/048.dat",
     > "scnn/049.dat",
     > "scnn/050.dat",
     > "scnn/051.dat",
     > "scnn/052.dat",
     > "scnn/053.dat",
     > "scnn/054.dat",
     > "scnn/055.dat",
     > "scnn/056.dat",
     > "scnn/057.dat",
     > "scnn/058.dat",
     > "scnn/059.dat",
     > "scnn/060.dat",
     > "scnn/061.dat",
     > "scnn/062.dat",
     > "scnn/063.dat",
     > "scnn/064.dat",
     > "scnn/065.dat",
     > "scnn/066.dat",
     > "scnn/067.dat",
     > "scnn/068.dat",
     > "scnn/069.dat",
     > "scnn/070.dat",
     > "scnn/071.dat",
     > "scnn/072.dat",
     > "scnn/073.dat",
     > "scnn/074.dat",
     > "scnn/075.dat",
     > "scnn/076.dat",
     > "scnn/077.dat",
     > "scnn/078.dat",
     > "scnn/079.dat",
     > "scnn/080.dat",
     > "scnn/081.dat",
     > "scnn/082.dat",
     > "scnn/083.dat",
     > "scnn/084.dat",
     > "scnn/085.dat",
     > "scnn/086.dat",
     > "scnn/087.dat",
     > "scnn/088.dat",
     > "scnn/089.dat",
     > "scnn/090.dat",
     > "scnn/091.dat",
     > "scnn/092.dat",
     > "scnn/093.dat",
     > "scnn/094.dat",
     > "scnn/095.dat",
     > "scnn/096.dat",
     > "scnn/097.dat",
     > "scnn/098.dat",
     > "scnn/099.dat",
     > "scnn/100.dat",
     > "scnn/101.dat"
     >            /)

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
!      CALL InitPDFsetByNameM(16,"LIKEnCC")
!      CALL InitPDFM(16,0)
      print *, "17"
!      CALL InitPDFsetByNameM(17,"LIKEnFE")
!      CALL InitPDFM(17,0)
      print *, "18"
!      CALL InitPDFsetByNameM(18,"LIKEnPB")
!      CALL InitPDFM(18,0)
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

      CALL READDATA

C-----CHOOSE EXP DATA SETS

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
      I_PIP_NE_PT2 = 1
      I_PIP_KR_PT2 = 1
      I_PIP_XE_PT2 = 1
      I_PI0_HE_PT2 = 0
      I_PI0_NE_PT2 = 1
      I_PI0_KR_PT2 = 1
      I_PI0_XE_PT2 = 1
      I_PIM_HE_PT2 = 0
      I_PIM_NE_PT2 = 1
      I_PIM_KR_PT2 = 1
      I_PIM_XE_PT2 = 1

      j = 10000
      do i = 1,101
      k = j+i
      open(unit= k,file=filenames(i))
      write(k ,*) 'GAMMA_FF ','chi2 ','num'
      enddo

      j = 10000
      do l = 1, 100
      print *, l
      XX(1)  = 1.5948389281365716E-002
      XX(2)  = 9.7265002153500910E-003
      if(l.lt.10) then
      XX(4) = 1d0/3d0 - 10**(-25d0/l)
      else
      XX(4) = 1d0/3d0 +10**(-25d0/(l-50d0))
      endif
      XX(3)  = 1/3d0
      XX(5)  = 0d0
      XX(6)  = 0d0
      chi2 = chisquare(XX)
      chi2s = (/
     >  CHI2_HERMES,
     >  CHI2_PIP_HE_NU,
     >  CHI2_PIP_HE_Z,
     >  CHI2_PIP_HE_Q2,
     >  CHI2_PIP_NE_NU,
     >  CHI2_PIP_NE_Z,
     >  CHI2_PIP_NE_Q2,
     >  CHI2_PIP_KR_NU,
     >  CHI2_PIP_KR_Z,
     >  CHI2_PIP_KR_Q2,
     >  CHI2_PIP_XE_NU,
     >  CHI2_PIP_XE_Z,
     >  CHI2_PIP_XE_Q2,
     >  CHI2_PIM_HE_NU,
     >  CHI2_PIM_HE_Z,
     >  CHI2_PIM_HE_Q2,
     >  CHI2_PIM_NE_NU,
     >  CHI2_PIM_NE_Z,
     >  CHI2_PIM_NE_Q2,
     >  CHI2_PIM_KR_NU,
     >  CHI2_PIM_KR_Z,
     >  CHI2_PIM_KR_Q2,
     >  CHI2_PIM_XE_NU,
     >  CHI2_PIM_XE_Z,
     >  CHI2_PIM_XE_Q2,
     >  CHI2_PI0_HE_NU,
     >  CHI2_PI0_HE_Z,
     >  CHI2_PI0_HE_Q2,
     >  CHI2_PI0_NE_NU,
     >  CHI2_PI0_NE_Z,
     >  CHI2_PI0_NE_Q2,
     >  CHI2_PI0_KR_NU,
     >  CHI2_PI0_KR_Z,
     >  CHI2_PI0_KR_Q2,
     >  CHI2_PI0_XE_NU,
     >  CHI2_PI0_XE_Z,
     >  CHI2_PI0_XE_Q2,
     >  CHI2_PIP_HE_PT2,
     >  CHI2_PIM_HE_PT2,
     >  CHI2_PI0_HE_PT2,
     >  CHI2_PIP_NE_PT2,
     >  CHI2_PIM_NE_PT2,
     >  CHI2_PI0_NE_PT2,
     >  CHI2_PIP_KR_PT2,
     >  CHI2_PIM_KR_PT2,
     >  CHI2_PI0_KR_PT2,
     >  CHI2_PIP_XE_PT2,
     >  CHI2_PIM_XE_PT2,
     >  CHI2_PI0_XE_PT2,
     >  CHI2_KP_HE_NU,
     >  CHI2_KP_HE_Z,
     >  CHI2_KP_HE_Q2,
     >  CHI2_KP_NE_NU,
     >  CHI2_KP_NE_Z,
     >  CHI2_KP_NE_Q2,
     >  CHI2_KP_KR_NU,
     >  CHI2_KP_KR_Z,
     >  CHI2_KP_KR_Q2,
     >  CHI2_KP_XE_NU,
     >  CHI2_KP_XE_Z,
     >  CHI2_KP_XE_Q2,
     >  CHI2_KM_HE_NU,
     >  CHI2_KM_HE_Z,
     >  CHI2_KM_HE_Q2,
     >  CHI2_KM_NE_NU,
     >  CHI2_KM_NE_Z,
     >  CHI2_KM_NE_Q2,
     >  CHI2_KM_KR_NU,
     >  CHI2_KM_KR_Z,
     >  CHI2_KM_KR_Q2,
     >  CHI2_KM_XE_NU,
     >  CHI2_KM_XE_Z,
     >  CHI2_KM_XE_Q2,
     >  CHI2_KP_HE_PT2,
     >  CHI2_KM_HE_PT2,
     >  CHI2_KP_NE_PT2,
     >  CHI2_KM_NE_PT2,
     >  CHI2_KP_KR_PT2,
     >  CHI2_KM_KR_PT2,
     >  CHI2_KP_XE_PT2,
     >  CHI2_KM_XE_PT2,
     >  CHI2_PIP_NU,
     >  CHI2_PIP_Z,
     >  CHI2_PIP_Q2,
     >  CHI2_PIP_PT2,
     >  CHI2_PI0_NU,
     >  CHI2_PI0_Z,
     >  CHI2_PI0_Q2,
     >  CHI2_PI0_PT2,
     >  CHI2_PIM_NU,
     >  CHI2_PIM_Z,
     >  CHI2_PIM_Q2,
     >  CHI2_PIM_PT2,
     >  CHI2_RHIC_Ratio_pAu1,
     >  CHI2_RHIC_Ratio_pAu2,
     >  CHI2_ATLAS5_Y1,
     >  CHI2_ATLAS5_Y2,
     >  CHI2_ATLAS5_Y3,
     >  CHI2_CMS5,
     >  CHI2_E866_800q,
     >  CHI2_E772_800 /)

      nums = (/
     >  NUM_HERMES,
     >  NUM_PIP_HE_NU,
     >  NUM_PIP_HE_Z,
     >  NUM_PIP_HE_Q2,
     >  NUM_PIP_NE_NU,
     >  NUM_PIP_NE_Z,
     >  NUM_PIP_NE_Q2,
     >  NUM_PIP_KR_NU,
     >  NUM_PIP_KR_Z,
     >  NUM_PIP_KR_Q2,
     >  NUM_PIP_XE_NU,
     >  NUM_PIP_XE_Z,
     >  NUM_PIP_XE_Q2,
     >  NUM_PIM_HE_NU,
     >  NUM_PIM_HE_Z,
     >  NUM_PIM_HE_Q2,
     >  NUM_PIM_NE_NU,
     >  NUM_PIM_NE_Z,
     >  NUM_PIM_NE_Q2,
     >  NUM_PIM_KR_NU,
     >  NUM_PIM_KR_Z,
     >  NUM_PIM_KR_Q2,
     >  NUM_PIM_XE_NU,
     >  NUM_PIM_XE_Z,
     >  NUM_PIM_XE_Q2,
     >  NUM_PI0_HE_NU,
     >  NUM_PI0_HE_Z,
     >  NUM_PI0_HE_Q2,
     >  NUM_PI0_NE_NU,
     >  NUM_PI0_NE_Z,
     >  NUM_PI0_NE_Q2,
     >  NUM_PI0_KR_NU,
     >  NUM_PI0_KR_Z,
     >  NUM_PI0_KR_Q2,
     >  NUM_PI0_XE_NU,
     >  NUM_PI0_XE_Z,
     >  NUM_PI0_XE_Q2,
     >  NUM_PIP_HE_PT2,
     >  NUM_PIM_HE_PT2,
     >  NUM_PI0_HE_PT2,
     >  NUM_PIP_NE_PT2,
     >  NUM_PIM_NE_PT2,
     >  NUM_PI0_NE_PT2,
     >  NUM_PIP_KR_PT2,
     >  NUM_PIM_KR_PT2,
     >  NUM_PI0_KR_PT2,
     >  NUM_PIP_XE_PT2,
     >  NUM_PIM_XE_PT2,
     >  NUM_PI0_XE_PT2,
     >  NUM_KP_HE_NU,
     >  NUM_KP_HE_Z,
     >  NUM_KP_HE_Q2,
     >  NUM_KP_NE_NU,
     >  NUM_KP_NE_Z,
     >  NUM_KP_NE_Q2,
     >  NUM_KP_KR_NU,
     >  NUM_KP_KR_Z,
     >  NUM_KP_KR_Q2,
     >  NUM_KP_XE_NU,
     >  NUM_KP_XE_Z,
     >  NUM_KP_XE_Q2,
     >  NUM_KM_HE_NU,
     >  NUM_KM_HE_Z,
     >  NUM_KM_HE_Q2,
     >  NUM_KM_NE_NU,
     >  NUM_KM_NE_Z,
     >  NUM_KM_NE_Q2,
     >  NUM_KM_KR_NU,
     >  NUM_KM_KR_Z,
     >  NUM_KM_KR_Q2,
     >  NUM_KM_XE_NU,
     >  NUM_KM_XE_Z,
     >  NUM_KM_XE_Q2,
     >  NUM_KP_HE_PT2,
     >  NUM_KM_HE_PT2,
     >  NUM_KP_NE_PT2,
     >  NUM_KM_NE_PT2,
     >  NUM_KP_KR_PT2,
     >  NUM_KM_KR_PT2,
     >  NUM_KP_XE_PT2,
     >  NUM_KM_XE_PT2,
     >  NUM_PIP_NU,
     >  NUM_PIP_Z,
     >  NUM_PIP_Q2,
     >  NUM_PIP_PT2,
     >  NUM_PI0_NU,
     >  NUM_PI0_Z,
     >  NUM_PI0_Q2,
     >  NUM_PI0_PT2,
     >  NUM_PIM_NU,
     >  NUM_PIM_Z,
     >  NUM_PIM_Q2,
     >  NUM_PIM_PT2,
     >  NUM_RHIC_Ratio_pAu1,
     >  NUM_RHIC_Ratio_pAu2,
     >  NUM_ATLAS5_Y1,
     >  NUM_ATLAS5_Y2,
     >  NUM_ATLAS5_Y3,
     >  NUM_CMS5,
     >  NUM_E866_800q,
     >  num_e772_800 /)
      do i = 1,101
      k = j+i
      write(k ,*) XX(3),chi2s(i),nums(i)
      enddo
      enddo

      j = 10000
      do i = 1,101
      k = j+i
      close(k)
      enddo

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

C-----HERMES
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

C-----RHIC (DY)
      real*8 CHI2_RHIC_Ratio_pAu1, CHI2_RHIC_Ratio_pAu2
      integer NUM_RHIC_Ratio_pAu1, NUM_RHIC_Ratio_pAu2

C-----LHC (DY)
      real*8 CHI2_ATLAS5_Y1,CHI2_ATLAS5_Y2,CHI2_ATLAS5_Y3
      real*8 CHI2_CMS5
      integer NUM_ATLAS5_Y1,NUM_ATLAS5_Y2,NUM_ATLAS5_Y3
      integer NUM_CMS5

C-----E866
      real*8  CHI2_E866_800q, CHI2_E772_800
      integer NUM_E866_800q , num_e772_800

      common /hermesch2/
     >  CHI2_HERMES,
     >  CHI2_PIP_HE_NU, CHI2_PIP_HE_Z, CHI2_PIP_HE_Q2,
     >  CHI2_PIP_NE_NU, CHI2_PIP_NE_Z, CHI2_PIP_NE_Q2,
     >  CHI2_PIP_KR_NU, CHI2_PIP_KR_Z, CHI2_PIP_KR_Q2,
     >  CHI2_PIP_XE_NU, CHI2_PIP_XE_Z, CHI2_PIP_XE_Q2,
     >  CHI2_PIM_HE_NU, CHI2_PIM_HE_Z, CHI2_PIM_HE_Q2,
     >  CHI2_PIM_NE_NU, CHI2_PIM_NE_Z, CHI2_PIM_NE_Q2,
     >  CHI2_PIM_KR_NU, CHI2_PIM_KR_Z, CHI2_PIM_KR_Q2,
     >  CHI2_PIM_XE_NU, CHI2_PIM_XE_Z, CHI2_PIM_XE_Q2,
     >  CHI2_PI0_HE_NU, CHI2_PI0_HE_Z, CHI2_PI0_HE_Q2,
     >  CHI2_PI0_NE_NU, CHI2_PI0_NE_Z, CHI2_PI0_NE_Q2,
     >  CHI2_PI0_KR_NU, CHI2_PI0_KR_Z, CHI2_PI0_KR_Q2,
     >  CHI2_PI0_XE_NU, CHI2_PI0_XE_Z, CHI2_PI0_XE_Q2,
     >  CHI2_PIP_HE_PT2, CHI2_PIM_HE_PT2, CHI2_PI0_HE_PT2,
     >  CHI2_PIP_NE_PT2, CHI2_PIM_NE_PT2, CHI2_PI0_NE_PT2,
     >  CHI2_PIP_KR_PT2, CHI2_PIM_KR_PT2, CHI2_PI0_KR_PT2,
     >  CHI2_PIP_XE_PT2, CHI2_PIM_XE_PT2, CHI2_PI0_XE_PT2,
     >  CHI2_KP_HE_NU, CHI2_KP_HE_Z, CHI2_KP_HE_Q2,
     >  CHI2_KP_NE_NU, CHI2_KP_NE_Z, CHI2_KP_NE_Q2,
     >  CHI2_KP_KR_NU, CHI2_KP_KR_Z, CHI2_KP_KR_Q2,
     >  CHI2_KP_XE_NU, CHI2_KP_XE_Z, CHI2_KP_XE_Q2,
     >  CHI2_KM_HE_NU, CHI2_KM_HE_Z, CHI2_KM_HE_Q2,
     >  CHI2_KM_NE_NU, CHI2_KM_NE_Z, CHI2_KM_NE_Q2,
     >  CHI2_KM_KR_NU, CHI2_KM_KR_Z, CHI2_KM_KR_Q2,
     >  CHI2_KM_XE_NU, CHI2_KM_XE_Z, CHI2_KM_XE_Q2,
     >  CHI2_KP_HE_PT2, CHI2_KM_HE_PT2,
     >  CHI2_KP_NE_PT2, CHI2_KM_NE_PT2,
     >  CHI2_KP_KR_PT2, CHI2_KM_KR_PT2,
     >  CHI2_KP_XE_PT2, CHI2_KM_XE_PT2,
     >  CHI2_PIP_NU,
     >  CHI2_PIP_Z,
     >  CHI2_PIP_Q2,
     >  CHI2_PIP_PT2,
     >  CHI2_PI0_NU,
     >  CHI2_PI0_Z,
     >  CHI2_PI0_Q2,
     >  CHI2_PI0_PT2,
     >  CHI2_PIM_NU,
     >  CHI2_PIM_Z,
     >  CHI2_PIM_Q2,
     >  CHI2_PIM_PT2,
     >  NUM_HERMES,
     >  NUM_PIP_HE_NU, NUM_PIP_HE_Z, NUM_PIP_HE_Q2,
     >  NUM_PIP_NE_NU, NUM_PIP_NE_Z, NUM_PIP_NE_Q2,
     >  NUM_PIP_KR_NU, NUM_PIP_KR_Z, NUM_PIP_KR_Q2,
     >  NUM_PIP_XE_NU, NUM_PIP_XE_Z, NUM_PIP_XE_Q2,
     >  NUM_PIM_HE_NU, NUM_PIM_HE_Z, NUM_PIM_HE_Q2,
     >  NUM_PIM_NE_NU, NUM_PIM_NE_Z, NUM_PIM_NE_Q2,
     >  NUM_PIM_KR_NU, NUM_PIM_KR_Z, NUM_PIM_KR_Q2,
     >  NUM_PIM_XE_NU, NUM_PIM_XE_Z, NUM_PIM_XE_Q2,
     >  NUM_PI0_HE_NU, NUM_PI0_HE_Z, NUM_PI0_HE_Q2,
     >  NUM_PI0_NE_NU, NUM_PI0_NE_Z, NUM_PI0_NE_Q2,
     >  NUM_PI0_KR_NU, NUM_PI0_KR_Z, NUM_PI0_KR_Q2,
     >  NUM_PI0_XE_NU, NUM_PI0_XE_Z, NUM_PI0_XE_Q2,
     >  NUM_PIP_HE_PT2, NUM_PIM_HE_PT2, NUM_PI0_HE_PT2,
     >  NUM_PIP_NE_PT2, NUM_PIM_NE_PT2, NUM_PI0_NE_PT2,
     >  NUM_PIP_KR_PT2, NUM_PIM_KR_PT2, NUM_PI0_KR_PT2,
     >  NUM_PIP_XE_PT2, NUM_PIM_XE_PT2, NUM_PI0_XE_PT2,
     >  NUM_KP_HE_NU, NUM_KP_HE_Z, NUM_KP_HE_Q2,
     >  NUM_KP_NE_NU, NUM_KP_NE_Z, NUM_KP_NE_Q2,
     >  NUM_KP_KR_NU, NUM_KP_KR_Z, NUM_KP_KR_Q2,
     >  NUM_KP_XE_NU, NUM_KP_XE_Z, NUM_KP_XE_Q2,
     >  NUM_KM_HE_NU, NUM_KM_HE_Z, NUM_KM_HE_Q2,
     >  NUM_KM_NE_NU, NUM_KM_NE_Z, NUM_KM_NE_Q2,
     >  NUM_KM_KR_NU, NUM_KM_KR_Z, NUM_KM_KR_Q2,
     >  NUM_KM_XE_NU, NUM_KM_XE_Z, NUM_KM_XE_Q2,
     >  NUM_KP_HE_PT2, NUM_KM_HE_PT2,
     >  NUM_KP_NE_PT2, NUM_KM_NE_PT2,
     >  NUM_KP_KR_PT2, NUM_KM_KR_PT2,
     >  NUM_KP_XE_PT2, NUM_KM_XE_PT2,
     >  NUM_PIP_NU,
     >  NUM_PIP_Z,
     >  NUM_PIP_Q2,
     >  NUM_PIP_PT2,
     >  NUM_PI0_NU,
     >  NUM_PI0_Z,
     >  NUM_PI0_Q2,
     >  NUM_PI0_PT2,
     >  NUM_PIM_NU,
     >  NUM_PIM_Z,
     >  NUM_PIM_Q2,
     >  NUM_PIM_PT2

      common /RHICch2/
     > CHI2_RHIC_Ratio_pAu1, CHI2_RHIC_Ratio_pAu2,
     > NUM_RHIC_Ratio_pAu1, NUM_RHIC_Ratio_pAu2

      common /LHCch2/
     > CHI2_ATLAS5_Y1,CHI2_ATLAS5_Y2,CHI2_ATLAS5_Y3,
     > CHI2_CMS5,
     > NUM_ATLAS5_Y1,NUM_ATLAS5_Y2,NUM_ATLAS5_Y3,
     > NUM_CMS5

      common /E866ch2/
     > CHI2_E866_800q, CHI2_E772_800,
     > NUM_E866_800q , num_e772_800

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
      data Sep_hermes/52.7d0/
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

      real*8 CHI2
      real*8 xb,zh,pht,tmp
      real*8 sigma
      real*8 dcorr,fN
      real*8 rts,y,Q,Qmin,Qmax,pt,fuu,Rds
      real*8 fuuA,fuuB, DIS, R_a
      real*8 R_dy
      real*8 xfmin, xfmax, ymin, ymax, Qbar, const, ybar
      real*8 xbmin, xbmax
      real*8 pi
      data pi/3.1415926535d0/
      real*8 M
      data M/0.938d0/
      real*8 CX_LHC,CX_RHIC
      integer i

      aN    = xx( 1)
      bN    = xx( 2)
      g2A   = xx( 3)
      a2    = xx( 4)
      g2B   = xx( 5)
      b2    = xx( 6)



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
      xbmin = 0.05d0
      xbmax = 0.3d0
      if (pt/Qbar.lt.qTdQcut) then
      call DY_overxbQ(4,rts,pt,xbmin,xbmax,Qmin,Qmax,fuuB,3)
      call DY_overxbQ(4,rts,pt,xbmin,xbmax,Qmin,Qmax,fuuA,12)
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
      xbmin = 0.05d0
      xbmax = 0.3d0
      if (pt/Qbar.lt.qTdQcut) then
      call DY_overxbQ(4,rts,pt,xbmin,xbmax,Qmin,Qmax,fuuB,3)
      call DY_overxbQ(4,rts,pt,xbmin,xbmax,Qmin,Qmax,fuuA,40)
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
      if (pt/Qbar.lt.qTdQcut) then
      call DY_overxbQ(4,rts,pt,xfmin,xfmax,Qmin,Qmax,fuuB,3)
      call DY_overxbQ(4,rts,pt,xfmin,xfmax,Qmin,Qmax,fuuA,56)
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
      if (pt/Qbar.lt.qTdQcut) then
      call DY_overxbQ(4,rts,pt,xfmin,xfmax,Qmin,Qmax,fuuB,3)
      call DY_overxbQ(4,rts,pt,xfmin,xfmax,Qmin,Qmax,fuuA,184)
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
      if (pt/Q.lt.qTdQcut) then
      call DY_overxF(2,rts,pt,xfmin,xfmax,Q,fuuB,9)
      call DY_overxF(2,rts,pt,xfmin,xfmax,Q,fuuA,56)
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
      if (pt/Q.lt.qTdQcut) then
      call DY_overxF(2,rts,pt,xfmin,xfmax,Q,fuuB,9)
      call DY_overxF(2,rts,pt,xfmin,xfmax,Q,fuuA,184)
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

      IF(I_PI0_HE_Z.EQ.1) THEN
      OPEN(UNIT = 105,FILE ='plot_data/HERMES/pi0_he_z.dat')
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
          rts = Sep_hermes
          IT = 20
          IH = 1
          IC = 0
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(rts,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(rts,Q2,xb,z,pht,fuu,IT,IH,IC)
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
          rts = Sep_hermes
          IT = 84
          IH = 1
          IC = 0
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(rts,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(rts,Q2,xb,z,pht,fuu,IT,IH,IC)
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
          rts = Sep_hermes
          IT = 131
          IH = 1
          IC = 0
          IF(  (PT2.LT.pt2cut).AND.( Z.LT.zcut )
     >    .AND.( (pht/(Q*z)).LT.pTdQdzcut) ) THEN
            CALL DISUU(rts,Q2,xb,z,pht,fuua,IT,IH,IC)
            IT = 3
            CALL DISUU(rts,Q2,xb,z,pht,fuu,IT,IH,IC)
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
         write(82,*)
         close(82)

      endif

      return
      end
