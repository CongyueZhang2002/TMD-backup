C***********************************************************************
C------This program generates the TMDPDF for a BOUND PROTON in a nucleus
C***********************************************************************
      program example
      implicit none
      real*8 x,Q, z
      integer Evolve3D_TMDPDF, Evolve3D_TMDFF
      integer plot3D_TMDPDF, plot3D_TMDFF
      real*8 up,ubp,dp,dbp,sp,sbp,glp
      real*8 u1,ub1,d1,db1,s1,sb1,gl1
      real*8 u2,ub2,d2,db2,s2,sb2,gl2
      real*8 kx,ky
      real*8 Q0
      integer i,j,k0, k1,k2,k3,k4,k5,k6,set
      integer usenc
      integer I_P, I_AU_Q0, I_AU_10, I_AU_MZ, I_P_10
      integer I_PB_Q0, I_PB_10
      integer I_P_XHI, I_AU_Q0_XHI, I_AU_10_XHI, I_AU_MZ_XHI, I_P_10_XHI
      integer I_PB_Q0_XHI, I_PB_10_XHI
      integer I_TEST
      Integer nucleon, Nucleus1, Nucleus2, Nucleus3
      common /nckern/ usenc
      integer nloops,hop,nll
      common /scheme/ nloops,hop,nll
      INTEGER FINI14
      COMMON / FRAGINI14 / FINI14
      real*8 kt

C-----perturbative orders:
      NLOOPS = 2
      HOP = 0
      NLL = 3

C-----Set Nuclei by MASS NUMBER:
      nucleon = 1
      Nucleus1 = 197 ! Au
      Nucleus2 = 208 ! Pb

C----- Set data to export
      I_P     = 0
      I_P_10  = 0
      I_AU_Q0 = 0
      I_AU_10 = 0
      I_AU_MZ = 0
      I_PB_Q0 = 1
      I_PB_10 = 1
C----- XHI
      I_P_XHI     = 1
      I_P_10_XHI  = 1
      I_AU_Q0_XHI = 1
      I_AU_10_XHI = 1
      I_AU_MZ_XHI = 0
      I_PB_Q0_XHI = 1
      I_PB_10_XHI = 1

      I_TEST  = 0


C-----INITIALIZE PDF SETS USING LHAPDF
      CALL InitPDFsetByNameM(1,"CT14nlo")
      CALL InitPDFM(1,0)
      CALL InitPDFsetByNameM(2,"EPPS16HE")
      CALL InitPDFM(2,0)
      CALL InitPDFsetByNameM(3,"EPPS16NE")
      CALL InitPDFM(3,0)
      CALL InitPDFsetByNameM(4,"EPPS16KR")
      CALL InitPDFM(4,0)
      CALL InitPDFsetByNameM(5,"EPPS16XE")
      CALL InitPDFM(5,0)
      CALL InitPDFsetByNameM(6,"LIKEnHE")
      CALL InitPDFM(6,0)
      CALL InitPDFsetByNameM(7,"LIKEnNE")
      CALL InitPDFM(7,0)
      CALL InitPDFsetByNameM(8,"LIKEnKR")
      CALL InitPDFM(8,0)
      CALL InitPDFsetByNameM(9,"LIKEnXE")
      CALL InitPDFM(9,0)
      CALL InitPDFsetByNameM(10,"EPPS16BE")
      CALL InitPDFM(10,0)
      CALL InitPDFsetByNameM(11,"EPPS16FE")
      CALL InitPDFM(11,0)
      CALL InitPDFsetByNameM(12,"EPPS16WW")
      CALL InitPDFM(12,0)
      CALL InitPDFsetByNameM(13,"EPPS16JCC")
      CALL InitPDFM(13,0)
      CALL InitPDFsetByNameM(14,"EPPS16JFE")
      CALL InitPDFM(14,0)
      CALL InitPDFsetByNameM(15,"EPPS16JPB")
      CALL InitPDFM(15,0)
      CALL InitPDFsetByNameM(20,"EPPS16AU")
      CALL InitPDFM(20,0)
      CALL InitPDFsetByNameM(21,"EPPS16PR")
      CALL InitPDFM(21,0)
      CALL InitPDFsetByNameM(22,"EPPS16CA")
      CALL InitPDFM(22,0)

      SET = 1
      call init_nTMDs(set)

      FINI14 = 0


C-----NUMBERING SYSTEM FOR PARTON TYPE FOR PDF_KT
C-----1: u   2:  d    3: s   4: ub   5: db   6: sb
C---------------------------------------------------

C---- Get TMPDFs vs. KT (Proton at Q0, x=0.1)
      IF(I_P.eq.1) then
      OPEN(UNIT = 4, FILE = 'PDF_p_Q0.csv')
      WRITE(4, *) 'kx ky u'
      x = 0.1
      Q = sqrt(2.4d0)
      DO i =1,60
      kx = -3 + 0.1*i
      DO J = 1,60
      ky = -3+ 0.1*j
      kt = sqrt(kx**2 + ky**2)
      CALL PDF_kt(kt,x,Q,1,nucleon,up)
      !CALL PDF_kt(kt,x,Q,2,nucleon,dp)
      WRITE(4,*) kx,ky, up
      ENDDO
      ENDDO
      CLOSE(4)
      ENDIF


C---- Get TMPDFs vs. KT (Proton at Q =10, x=0.1)
      IF(I_P_10.eq.1) then
      OPEN(UNIT = 4, FILE = 'PDF_p_10.csv')
      WRITE(4, *) 'kx ky u'
      x = 0.1
      Q = 10
      DO i =1,60
      kx = -3 + 0.1*i
      DO J = 1,60
      ky = -3 + 0.1*j
      kt = sqrt(kx**2 + ky**2)
      CALL PDF_kt(kt,x,Q,1,nucleon,up)
      !CALL PDF_kt(kt,x,Q,2,nucleon,dp)
      WRITE(4,*) kx,ky, up
      ENDDO
      ENDDO
      CLOSE(4)
      ENDIF



C------- AU
      IF(I_AU_Q0.eq.1) then
C---- Get TMPDFs vs. KT (Au at Q0, x = 0.1)
      OPEN(UNIT = 4, FILE = 'PDF_Au_Q0.csv')
      WRITE(4, *) 'kx ky u'
      x = 0.1
      Q = sqrt(2.4d0)
      DO i =1,60
      kx = -3 + 0.1*i
      DO J = 1,60
      ky = -3+ 0.1*j
      kt = sqrt(kx**2 + ky**2)
      CALL PDF_kt(kt,x,Q,1,Nucleus1,up)
      !CALL PDF_kt(kt,x,Q,2,Nucleus1,dp)
      WRITE(4,*) kx, ky, up
      ENDDO
      ENDDO
      CLOSE(4)
      ENDIF


      IF(I_AU_10.eq.1) then
C---- Get TMPDFs vs. KT (Au at 10 GeV, x = 0.1)
      OPEN(UNIT = 4, FILE = 'PDF_Au_10.csv')
      WRITE(4, *) 'kx ky u'
      x = 0.1
      Q = 10
      DO i =1,60
      kx = -3 + 0.1*i
      DO J = 1,60
      ky = -3+ 0.1*j
      kt = sqrt(kx**2 + ky**2)
      CALL PDF_kt(kt,x,Q,1,Nucleus1,up)
      !CALL PDF_kt(kt,x,Q,2,Nucleus1,dp)
      WRITE(4,*) kx,ky, up
      ENDDO
      ENDDO
      CLOSE(4)
      ENDIF


C------- AU
      IF(I_PB_Q0.eq.1) then
C---- Get TMPDFs vs. KT (Au at Q0, x = 0.1)
      OPEN(UNIT = 4, FILE = 'PDF_Pb_Q0.csv')
      WRITE(4, *) 'kx ky u'
      x = 0.1
      Q = sqrt(2.4d0)
      DO i =1,60
      kx = -3 + 0.1*i
      DO J = 1,60
      ky = -3+ 0.1*j
      kt = sqrt(kx**2 + ky**2)
      CALL PDF_kt(kt,x,Q,1,Nucleus2,up)
      !CALL PDF_kt(kt,x,Q,2,Nucleus1,dp)
      WRITE(4,*) kx, ky, up
      ENDDO
      ENDDO
      CLOSE(4)
      ENDIF


      IF(I_PB_10.eq.1) then
C---- Get TMPDFs vs. KT (Au at 10 GeV, x = 0.1)
      OPEN(UNIT = 4, FILE = 'PDF_Pb_10.csv')
      WRITE(4, *) 'kx ky u'
      x = 0.1
      Q = 10
      DO i =1,60
      kx = -3 + 0.1*i
      DO J = 1,60
      ky = -3+ 0.1*j
      kt = sqrt(kx**2 + ky**2)
      CALL PDF_kt(kt,x,Q,1,Nucleus2,up)
      !CALL PDF_kt(kt,x,Q,2,Nucleus1,dp)
      WRITE(4,*) kx,ky, up
      ENDDO
      ENDDO
      CLOSE(4)
      ENDIF


C---- Get TMPDFs vs. KT (Proton at Q0, x=0.1)
      IF(I_P_XHI.eq.1) then
      OPEN(UNIT = 4, FILE = 'XHI/PDF_p_Q0.csv')
      WRITE(4, *) 'kx ky u'
      x = 0.4
      Q = sqrt(2.4d0)
      DO i =1,20
      kx = -1 + 0.1*i
      DO J = 1,20
      ky = -1+ 0.1*j
      kt = sqrt(kx**2 + ky**2)
      CALL PDF_kt(kt,x,Q,1,nucleon,up)
      !CALL PDF_kt(kt,x,Q,2,nucleon,dp)
      WRITE(4,*) kx,ky, up
      ENDDO
      ENDDO
      CLOSE(4)
      ENDIF


C---- Get TMPDFs vs. KT (Proton at Q =10, x=0.1)
      IF(I_P_10_XHI.eq.1) then
      OPEN(UNIT = 4, FILE = 'XHI/PDF_p_10.csv')
      WRITE(4, *) 'kx ky u'
      x = 0.4
      Q = 10
      DO i =1, 20
      kx = -1 + 0.1*i
      DO J = 1, 20
      ky = -1 + 0.1*j
      kt = sqrt(kx**2 + ky**2)
      CALL PDF_kt(kt,x,Q,1,nucleon,up)
      !CALL PDF_kt(kt,x,Q,2,nucleon,dp)
      WRITE(4,*) kx,ky, up
      ENDDO
      ENDDO
      CLOSE(4)
      ENDIF

C----- XHI

C------- AU
      IF(I_AU_Q0_XHI.eq.1) then
C---- Get TMPDFs vs. KT (Au at Q0, x = 0.1)
      OPEN(UNIT = 4, FILE = 'XHI/PDF_Au_Q0.csv')
      WRITE(4, *) 'kx ky u'
      x = 0.4
      Q = sqrt(2.4d0)
      DO i =1, 20
      kx = -1 + 0.1*i
      DO J = 1,20
      ky = -1+ 0.1*j
      kt = sqrt(kx**2 + ky**2)
      CALL PDF_kt(kt,x,Q,1,Nucleus1,up)
      !CALL PDF_kt(kt,x,Q,2,Nucleus1,dp)
      WRITE(4,*) kx, ky, up
      ENDDO
      ENDDO
      CLOSE(4)
      ENDIF


      IF(I_AU_10_XHI.eq.1) then
C---- Get TMPDFs vs. KT (Au at 10 GeV, x = 0.1)
      OPEN(UNIT = 4, FILE = 'XHI/PDF_Au_10.csv')
      WRITE(4, *) 'kx ky u'
      x = 0.4
      Q = 10
      DO i =1,20
      kx = -1 + 0.1*i
      DO J = 1,20
      ky = -1+ 0.1*j
      kt = sqrt(kx**2 + ky**2)
      CALL PDF_kt(kt,x,Q,1,Nucleus1,up)
      !CALL PDF_kt(kt,x,Q,2,Nucleus1,dp)
      WRITE(4,*) kx,ky, up
      ENDDO
      ENDDO
      CLOSE(4)
      ENDIF


C------- AU
      IF(I_PB_Q0_XHI.eq.1) then
C---- Get TMPDFs vs. KT (Au at Q0, x = 0.1)
      OPEN(UNIT = 4, FILE = 'XHI/PDF_Pb_Q0.csv')
      WRITE(4, *) 'kx ky u'
      x = 0.4
      Q = sqrt(2.4d0)
      DO i =1,20
      kx = -1 + 0.1*i
      DO J = 1,20
      ky = -1+ 0.1*j
      kt = sqrt(kx**2 + ky**2)
      CALL PDF_kt(kt,x,Q,1,Nucleus2,up)
      !CALL PDF_kt(kt,x,Q,2,Nucleus1,dp)
      WRITE(4,*) kx, ky, up
      ENDDO
      ENDDO
      CLOSE(4)
      ENDIF


      IF(I_PB_10_XHI.eq.1) then
C---- Get TMPDFs vs. KT (Au at 10 GeV, x = 0.1)
      OPEN(UNIT = 4, FILE = 'XHI/PDF_Pb_10.csv')
      WRITE(4, *) 'kx ky u'
      x = 0.4
      Q = 10
      DO i =1,20
      kx = -1 + 0.1*i
      DO J = 1,20
      ky = -1+ 0.1*j
      kt = sqrt(kx**2 + ky**2)
      CALL PDF_kt(kt,x,Q,1,Nucleus2,up)
      !CALL PDF_kt(kt,x,Q,2,Nucleus1,dp)
      WRITE(4,*) kx,ky, up
      ENDDO
      ENDDO
      CLOSE(4)
      ENDIF


C---- Get TMPDFs vs. KT (Proton at Q0)
      IF(I_TEST.eq.1) then
      OPEN(UNIT = 4, FILE = 'x01_LOW.csv')
      WRITE(4, *) 'kt up u1'
      x = 0.1
      Q = sqrt(2.4)
      DO i =1,10
      kt = 0.1*i
      CALL PDF_kt(kt,x,Q,1,nucleon,up)
      CALL PDF_kt(kt,x,Q,1,nucleus1,u1)
      WRITE(4,*) kt, up, u1
      ENDDO
      CLOSE(4)
      ENDIF

      print*, "done"


C---- Get TMPDFs vs. KT (Proton at Q0)
      IF(I_TEST.eq.1) then
      OPEN(UNIT = 4, FILE = 'x01_HIGH.csv')
      WRITE(4, *) 'kt up u1'
      x = 0.1
      Q = 10
      DO i =1,10
      kt = 0.1*i
      CALL PDF_kt(kt,x,Q,1,nucleon,up)
      CALL PDF_kt(kt,x,Q,1,nucleus1,u1)
      WRITE(4,*) kt, up, u1
      ENDDO
      CLOSE(4)
      ENDIF

      print*, "done"























C---- Get TMPDFs vs. KT (Proton at Q0)
      IF(I_TEST.eq.1) then
      OPEN(UNIT = 4, FILE = 'x02_LOW.csv')
      WRITE(4, *) 'kt up u1'
      x = 0.2
      Q = sqrt(2.4)
      DO i =1,10
      kt = 0.1*i
      CALL PDF_kt(kt,x,Q,1,nucleon,up)
      CALL PDF_kt(kt,x,Q,1,nucleus1,u1)
      WRITE(4,*) kt, up, u1
      ENDDO
      CLOSE(4)
      ENDIF

      print*, "done"


C---- Get TMPDFs vs. KT (Proton at Q0)
      IF(I_TEST.eq.1) then
      OPEN(UNIT = 4, FILE = 'x02_HIGH.csv')
      WRITE(4, *) 'kt up u1'
      x = 0.2
      Q = 10
      DO i =1,10
      kt = 0.1*i
      CALL PDF_kt(kt,x,Q,1,nucleon,up)
      CALL PDF_kt(kt,x,Q,1,nucleus1,u1)
      WRITE(4,*) kt, up, u1
      ENDDO
      CLOSE(4)
      ENDIF

      print*, "done"


C---- Get TMPDFs vs. KT (Proton at Q0)
      IF(I_TEST.eq.1) then
      OPEN(UNIT = 4, FILE = 'x03_LOW.csv')
      WRITE(4, *) 'kt up u1'
      x = 0.3
      Q = sqrt(2.4)
      DO i =1,10
      kt = 0.1*i
      CALL PDF_kt(kt,x,Q,1,nucleon,up)
      CALL PDF_kt(kt,x,Q,1,nucleus1,u1)
      WRITE(4,*) kt, up, u1
      ENDDO
      CLOSE(4)
      ENDIF

      print*, "done"



C---- Get TMPDFs vs. KT (Proton at Q0)
      IF(I_TEST.eq.1) then
      OPEN(UNIT = 4, FILE = 'x03_HIGH.csv')
      WRITE(4, *) 'kt up u1'
      x = 0.3
      Q = 10
      DO i =1,10
      kt = 0.1*i
      CALL PDF_kt(kt,x,Q,1,nucleon,up)
      CALL PDF_kt(kt,x,Q,1,nucleus1,u1)
      WRITE(4,*) kt, up, u1
      ENDDO
      CLOSE(4)
      ENDIF
      print*, "done"


C---- Get TMPDFs vs. KT (Proton at Q0)
      IF(I_TEST.eq.1) then
      OPEN(UNIT = 4, FILE = 'x04_LOW.csv')
      WRITE(4, *) 'kt up u1'
      x = 0.4
      Q = sqrt(2.4)
      DO i =1,10
      kt = 0.1*i
      CALL PDF_kt(kt,x,Q,1,nucleon,up)
      CALL PDF_kt(kt,x,Q,1,nucleus1,u1)
      WRITE(4,*) kt, up, u1
      ENDDO
      CLOSE(4)
      ENDIF
      print*, "almost done"



C---- Get TMPDFs vs. KT (Proton at Q0)
      IF(I_TEST.eq.1) then
      OPEN(UNIT = 4, FILE = 'x04_HIGH.csv')
      WRITE(4, *) 'kt up u1'
      x = 0.4
      Q = 10
      DO i =1,10
      kt = 0.1*i
      CALL PDF_kt(kt,x,Q,1,nucleon,up)
      CALL PDF_kt(kt,x,Q,1,nucleus1,u1)
      WRITE(4,*) kt, up, u1
      ENDDO
      CLOSE(4)
      ENDIF

      print*, "done done done"





      END PROGRAM

      subroutine init_nTMDs(set)
      implicit none
      integer set, I
      real*8 Ans ,Bns, GDs, ADs, GPs, APs
      common /fitp/
     >       Ans ,Bns ,GDs, Ads, GPs, APs
      real*8 ANA(101), BNA(101), GDA(101), ADA(101), GPA(101), APA(101)
      COMMON /param_dat/ ANA, BNA, GDA, ADA, GPA, APA

      CALL READDATA

      ANs =       ANA(set)
      BNs =       BNA(set)
      GPs =       0d0
      APs =       0d0
      GDs =       0d0
      ADs =       0d0

      !DO I=1,101
         !PRINT*, I, ANA(I), BNA(I), GBA(I)
      !END DO

      end
