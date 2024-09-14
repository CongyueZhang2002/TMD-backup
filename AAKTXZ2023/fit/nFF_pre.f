      PROGRAM MASTER
      implicit none
      real*8 x,Q1,Q2
      real*8 U0, U1, ratio
      integer i, set, A1, A2, IH, IC
      integer nloops
      real*8 kt
      INTEGER FINI14
      COMMON / FRAGINI14 / FINI14

*     Input COMMON BLOCKS !2023
      double precision H_AA, A
      COMMON /HMASS/ H_AA  
      COMMON /meson/ IH,IC
      real*8 fu,fub,fd,fdb,fs,fsb,fg,fp,fpb1,fpb0,func
      real*8 VSTRT(12), mub, z

      REAL*8 gamma, g3f, g3D, Nq1, gq1, dq1,
     &       Nq2, gq2, dq2, p_10, p_11, p_12 !2023 
      COMMON /FITP/ gamma, g3f, g3D, Nq1, gq1, dq1,
     &              Nq2, gq2, dq2, p_10, p_11, p_12 !2023      
      REAL*8 gammaBEST,g3fBEST,g3DBEST,Nq1BEST,gq1BEST,dq1BEST,
     &       Nq2BEST,gq2BEST,dq2BEST,p_10BEST,p_11BEST,p_12BEST !2023

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
      REAL*8 r1, r2, r3, r4, r5, r6 !2023
      REAL*8 r7, r8, r9, r10, r11, r12
      integer seed,initseed
      real*8 r0,r01000

      IH = 1
      IC = 1
      A = 208
      H_AA = A

      INCLUDE "best-params-rand.f"

      gamma  = gammaBEST
      g3f  = g3fBEST
      g3D  = g3DBEST
      Nq1  = Nq1BEST
      gq1  = gq1BEST
      dq1  = dq1BEST
      Nq2  = Nq2BEST
      gq2  = gq2BEST
      dq2  = dq2BEST
      p_10 = p_10BEST
      p_11 = p_11BEST
      p_12 = p_12BEST

C-----Get FF

      OPEN(UNIT = 1000,
     *STATUS='UNKNOWN',FILE ='plot_data/TMD/nFF_pre.dat')
      WRITE(1000,*)  'x ', 'kt ' , 'FF(p) ', 'FF(Pb) ', 'ratio'

      DO i = 1, 100
        z = 0.2 + (i-1)*(0.4/99)  
        kt = 0.0
 
        ratio = 0
        U0 = 0
        U1 = func(2,z)

        WRITE(1000,*) z,kt,U0,U1,ratio
      ENDDO

      CLOSE(1000)
      End
