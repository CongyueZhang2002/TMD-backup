      PROGRAM MASTER
      implicit none
      real*8 x,Q1,Q2
      real*8 U0, ratio
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

      integer j,k
      integer, dimension(8) :: A_array = (/4,12,20,56,84,131,197,208/)

      integer IS
      real*8 U1,UB1,D1,DB1,S1,SB1,C1,B1,GL1,U
      
      FINI14 = 0
      Q = 1.0
      nloops = 1
      IH = 1
      IC = 1

      call readdata
      print *, replica_num

      IF (REPLICA_NUM .LE. 120) THEN
        IS = 0
      ELSE
        j = REPLICA_NUM - 120
        if (mod(j, 2) == 0) then
          k = j
          IS = k/2
        else
          k = j+1
          IS = -k/2
        endif
      ENDIF
      print *, 'IS=', IS
C-------------------------------------------------------------
      OPEN(UNIT = 1001,
     *STATUS='UNKNOWN',FILE ='plot_data/TMD/DEHSS.dat')
      WRITE(1001,*)  'x ', 'kt ' , 'FF(p) ', 'FF(Pb) ', 'ratio'

      DO i = 1, 100
        z = 0.06+(i-1)*(0.94/99)-0.0000001  
        kt = 0.0
      
        ratio = 0
        U0 = 0
        Q2 = 1.0000001
        CALL fDSSH14(IS,1,1,1,z,Q2,U1,UB1,D1,DB1,S1,SB1,C1,B1,GL1)
        U = U1

        WRITE(1001,*) z,kt,U0,U,ratio
      ENDDO

      END
