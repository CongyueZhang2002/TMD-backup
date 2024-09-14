      PROGRAM MASTER
      implicit none
      real*8 x,Q1,Q2
      real*8 U0, U1, ratio
      integer i, set, A1, A2, IH, IC, A
      integer nloops
      real*8 kt
      INTEGER FINI14
      COMMON / FRAGINI14 / FINI14
      integer replica_num
      COMMON /REPLICA_INFO/ REPLICA_NUM

*     Input COMMON BLOCKS !2023
      double precision H_AA
      COMMON /HMASS/ H_AA  
      COMMON /meson/ IH,IC
      real*8 fu,fub,fd,fdb,fs,fsb,fg,fp,fpb1,fpb0,func
      real*8 VSTRT(9), Q, z

      REAL*8 aN, bN, Ng1, Nq1, Ng2, Nq2, bg1, gg1, dg1 !2023
      COMMON /FITP/ aN, bN, Ng1, Nq1, Ng2, Nq2, bg1, gg1, dg1 !2023

      integer j
      integer, dimension(8) :: A_array = (/4,12,20,56,84,131,197,208/)

      integer IS
      real U1,UB1,D1,DB1,S1,SB1,C1,B1,GL1,U

      Q = 1.0
      nloops = 1
      IH = 1
      IC = 1

      call readdata
      print *, replica_num

      IF(REPLICA_NUM.le.120) THEN
      IS = 0
      ELSE 
      IS = (REPLICA_NUM - 120)

      ENDIF

C---------------------------------------------------------------
      OPEN(UNIT = 1000,
     *STATUS='UNKNOWN',FILE ='plot_data/TMD/nFF_col.dat')
      WRITE(1000,*)  'x ', 'kt ' , 'FF(p) ', 'FF(Pb) ', 'ratio'

      DO j = 1, 8
      DO i = 1, 100
        A = A_array(j)
        print *, A
        z = 0.01+(i-1)*(0.99/99)  
        kt = 0.0
      
        ratio = 0
        U0 = 0
        CALL LIKEnFF(z,Q,A,fu,fub,fd,fdb,fs,fsb,fg,nloops)
        U1 = z*fu

        WRITE(1000,*) z,kt,U0,U1,ratio
      ENDDO
      ENDDO

      CLOSE(1000)

C-------------------------------------------------------------
      OPEN(UNIT = 1001,
     *STATUS='UNKNOWN',FILE ='plot_data/TMD/DEHSS.dat')
      WRITE(1001,*)  'x ', 'kt ' , 'FF(p) ', 'FF(Pb) ', 'ratio'

      DO i = 1, 100
        z = 0.01+(i-1)*(0.99/99)  
        kt = 0.0
      
        ratio = 0
        U0 = 0
        CALL call fDSSH(IS,1,1,1,z,1.0,U1,UB1,D1,DB1,S1,SB1,C1,B1,GL1)
        U = z*U1

        WRITE(1001,*) z,kt,U0,U,ratio
      ENDDO

      END
