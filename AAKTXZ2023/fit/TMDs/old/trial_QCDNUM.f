        PROGRAM main
        implicit none
    
        INTEGER :: CL, IH, IC, SET, nloops, REPLICA_NUM, AA
        DOUBLE PRECISION :: Z, Q2, U, UB, D, DB, S, SB, C, CB, B, BB, GL
        DOUBLE PRECISION :: zz, mub, fu,fub,fd,fdb,fs,fsb,fg,H_AA 
        INTEGER :: INIT, i

        CHARACTER*20 file_name1, file_name2, file_name3, file_name4
        CHARACTER*20 file_name5, file_name6, file_name7, file_name8
        INTEGER unit_number


      DOUBLE PRECISION :: array(47), def(-6:6,12), qq(2),wt(2)
      DOUBLE PRECISION :: pdf(-6:6) 
      Double precision as0, r20, xmin, eps 
      integer :: iord, nfin, itype, iset, jset, iosp, nx
      integer :: nxin, nqin, lun, idbug, iqc, iqb, iq0, nq 
      double precision :: q2c, q2b, q0 
      double precision Qf 
      double precision Qf2  

C------------- X grid and mu2 grid paramaters
      data idbug/0/                                     !debug printpout
      data xmin/1.D-4/, nxin/100/                       !x grid, 
      data qq/1.D0,1.D5/, wt/1.D0,1.D0/, nqin/100/       !mu2 grid


C------ Z-array for output file 
      DATA array/0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
     4        0.095, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
     5        0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475,  0.5,
     6        0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7,
     7        0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9,
     8        0.925, 0.95, 0.975, 1.0/

C--------------------------------------------------------------   
      external func                       !input parton dists
      integer, external :: iqfrmq
      data def  /                         !flavor composition
C--   tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
C--   -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6  
     + 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.,         !d
     + 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.,         !u
     + 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0.,         !s
     + 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,         !dbar
     + 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,         !ubar
     + 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,         !sbar
     + 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,         !c
     + 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,         !cbar
     + 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,         !b
     + 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,         !bbar
     + 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,         !t
     + 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.  /       !tbar 
                                       !pdfout

                                       !pdfout

C--   Weight files
      character*26 fnam(3)

      data fnam /'weights/unpolarised.wgt',
     +           'weights/polarised.wgt  ',
     +           'weights/timelike.wgt   '/

        COMMON / INITIATE / INIT
        COMMON /HMASS/ H_AA 
        COMMON /REPLICA_INFO/ REPLICA_NUM
        COMMON /meson/ IH,IC
        ! parameter(eps=1d-10)            


        INIT = 0
    
        CL = 1  
        IH = 1  
        IC = 1  
        SET = 0    
        mub = 2.00D0
        Q2 = mub**2
        nloops = 1     
        !Q0 = dsqrt(1d0) - eps 
        Q0 = 1d0


        file_name1 = 'like_H.csv'
        file_name2 = 'like_Ne.csv'
        file_name3 = 'like_Kr.csv'
        file_name4 = 'like_Xe.csv'
        file_name5 = 'new_H.csv'
        file_name6 = 'new_Ne.csv'
        file_name7 = 'new_Kr.csv'
        file_name8 = 'new_Xe.csv'        

        OPEN (UNIT = 20, FILE = file_name1, STATUS = 'REPLACE')
        INIT = 0        
       write(20 ,*) 'Z ', 'U ', 'UB ', 'D ', 'DB ', 'S ', 'SB ', 'C ',
     9        'CB ', 'B ', 'BB ', 'GL '
        Do i=1,47
            Z = array(i)
            call LIKEn(CL,IH,IC,set,array(i),Q2,1,U,UB,D,DB,S,SB,C,CB,
     9      B,BB,GL)
            
        write(20, *), Z, U/Z, UB/Z, D/Z, DB/Z, S/Z, SB/Z, 
     9                      C/Z, CB/Z, B/Z, BB/Z, GL/Z 
            
        End do
        CLOSE (UNIT = 20)

        OPEN (UNIT = 20, FILE = file_name2, STATUS = 'REPLACE')
        INIT = 0
        write(20 ,*) 'Z ', 'U ', 'UB ', 'D ', 'DB ', 'S ', 'SB ', 'C ',
     9        'CB ', 'B ', 'BB ', 'GL '
        Do i=1,47
            Z = array(i)
            call LIKEn(CL,IH,IC,set,array(i),Q2,20,U,UB,D,DB,S,SB,C,CB,
     9      B,BB,GL)
            
        write(20, *), Z, U/Z, UB/Z, D/Z, DB/Z, S/Z, SB/Z, 
     9                      C/Z, CB/Z, B/Z, BB/Z, GL/Z 
            
        End do
        CLOSE (UNIT = 20)


        OPEN (UNIT = 20, FILE = file_name3, STATUS = 'REPLACE')
        INIT = 0
        write(20 ,*) 'Z ', 'U ', 'UB ', 'D ', 'DB ', 'S ', 'SB ', 'C ',
     9        'CB ', 'B ', 'BB ', 'GL '
        Do i=1,47
            Z = array(i)
            call LIKEn(CL,IH,IC,set,array(i),Q2,84,U,UB,D,DB,S,SB,C,CB,
     9      B,BB,GL)
            
        write(20, *), Z, U/Z, UB/Z, D/Z, DB/Z, S/Z, SB/Z, 
     9                      C/Z, CB/Z, B/Z, BB/Z, GL/Z 
            
        End do
        CLOSE (UNIT = 20)


        OPEN (UNIT = 20, FILE = file_name4, STATUS = 'REPLACE')
        INIT = 0
        write(20 ,*) 'Z ', 'U ', 'UB ', 'D ', 'DB ', 'S ', 'SB ', 'C ',
     9        'CB ', 'B ', 'BB ', 'GL '
        Do i=1,47
            Z = array(i)
            call LIKEn(CL,IH,IC,set,array(i),Q2,131,U,UB,D,DB,S,SB,C,CB,
     9      B,BB,GL)
            
        write(20, *), Z, U/Z, UB/Z, D/Z, DB/Z, S/Z, SB/Z, 
     9                      C/Z, CB/Z, B/Z, BB/Z, GL/Z 
            
        End do
        CLOSE (UNIT = 20)

C-------QCDNUM 
C----- Set-up -----------------------------------------------------------
        AA=1
        H_AA=1
        lun = 6                                   ! -6 supresses banner
        call qcinit(lun,' ')                       ! initialize
        iosp = 3 ! 2 for linear interpolarion, 3 for spline       
        call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
        call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
        itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
        call wtfile(itype,fnam(itype))   !calculate weights
        iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
        call setord(iord)                                 
        as0 = 0.364d0 
        r20 = 2.0d0 
        call setalf(as0,r20)                              !input alphas
        q0 = 1.0D0 
        q2c = (1.43d0)**2d0 
        q2b = (4.3d0)**2d0 
        iqc  = iqfrmq(q2c)    !Charm threshold
        iqb  = iqfrmq(q2b)    !Bottom threshold
        nfin = 1 
        call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
        iset = 1 
        jset = 10*iset+itype            !output pdf set and evolution type
        iq0  = iqfrmq(q0)                !start scale index 
        call setint('edbg',idbug)        !debug printout

        call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

        OPEN (UNIT = 10, FILE = file_name5, STATUS = 'Unknown')
         write(10 ,*) 'Z ', 'U ', 'UB ', 'D ', 'DB ', 'S ', 'SB ', 'GL'
            
        Do i=1,47
            zz = array(i)
            Call Nuclear_TMDFF(zz,mub,AA,fu,fub,fd,fdb,fs,fsb,fg,nloops)
                
            write(10, *), zz,fu,fub,fd,fdb,fs,fsb,fg
                
        End do

        CLOSE (UNIT = 10)

C----- Set-up -----------------------------------------------------------
        AA= 20
        H_AA=20
        lun = 6                                   ! -6 supresses banner
        call qcinit(lun,' ')                       ! initialize
        iosp = 3 ! 2 for linear interpolarion, 3 for spline       
        call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
        call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
        itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
        call wtfile(itype,fnam(itype))   !calculate weights
        iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
        call setord(iord)                                 
        as0 = 0.364d0 
        r20 = 2.0d0 
        call setalf(as0,r20)                              !input alphas
        q0 = 1.0D0 
        q2c = (1.43d0)**2d0 
        q2b = (4.3d0)**2d0 
        iqc  = iqfrmq(q2c)    !Charm threshold
        iqb  = iqfrmq(q2b)    !Bottom threshold
        nfin = 1 
        call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
        iset = 1 
        jset = 10*iset+itype            !output pdf set and evolution type
        iq0  = iqfrmq(q0)                !start scale index 
        call setint('edbg',idbug)        !debug printout

        call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

       OPEN (UNIT = 10, FILE = file_name6, STATUS = 'Unknown')
         write(10 ,*) 'Z ', 'U ', 'UB ', 'D ', 'DB ', 'S ', 'SB ', 'GL'
            
        Do i=1,47
            zz = array(i)
            Call Nuclear_TMDFF(zz,mub,AA,fu,fub,fd,fdb,fs,fsb,fg,nloops)
                
            write(10, *), zz,fu,fub,fd,fdb,fs,fsb,fg
                
        End do

        CLOSE (UNIT = 10)
C----- Set-up -----------------------------------------------------------
        AA=84
        H_AA=84
        lun = 6                                   ! -6 supresses banner
        call qcinit(lun,' ')                       ! initialize
        iosp = 3 ! 2 for linear interpolarion, 3 for spline       
        call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
        call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
        itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
        call wtfile(itype,fnam(itype))   !calculate weights
        iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
        call setord(iord)                                 
        as0 = 0.364d0 
        r20 = 2.0d0 
        call setalf(as0,r20)                              !input alphas
        q0 = 1.0D0 
        q2c = (1.43d0)**2d0 
        q2b = (4.3d0)**2d0 
        iqc  = iqfrmq(q2c)    !Charm threshold
        iqb  = iqfrmq(q2b)    !Bottom threshold
        nfin = 1 
        call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
        iset = 1 
        jset = 10*iset+itype            !output pdf set and evolution type
        iq0  = iqfrmq(q0)                !start scale index 
        call setint('edbg',idbug)        !debug printout

        call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
C------------------------------------------------------------------------

         OPEN (UNIT = 10, FILE = file_name7, STATUS = 'Unknown')
         write(10 ,*) 'Z ', 'U ', 'UB ', 'D ', 'DB ', 'S ', 'SB ', 'GL'
            
        Do i=1,47
            zz = array(i)
            Call Nuclear_TMDFF(zz,mub,AA,fu,fub,fd,fdb,fs,fsb,fg,nloops)
                
            write(10, *), zz,fu,fub,fd,fdb,fs,fsb,fg
                
        End do

        CLOSE (UNIT = 10)

C----- Set-up -----------------------------------------------------------
        AA=131
        H_AA=131
        lun = 6                                   ! -6 supresses banner
        call qcinit(lun,' ')                       ! initialize
        iosp = 3 ! 2 for linear interpolarion, 3 for spline       
        call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
        call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
        itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
        call wtfile(itype,fnam(itype))   !calculate weights
        iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
        call setord(iord)                                 
        as0 = 0.364d0 
        r20 = 2.0d0 
        call setalf(as0,r20)                              !input alphas
        q0 = 1.0D0 
        q2c = (1.43d0)**2d0 
        q2b = (4.3d0)**2d0 
        iqc  = iqfrmq(q2c)    !Charm threshold
        iqb  = iqfrmq(q2b)    !Bottom threshold
        nfin = 1 
        call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
C------Evolution --------------------------------------------------------
        iset = 1 
        jset = 10*iset+itype            !output pdf set and evolution type
        iq0  = iqfrmq(q0)                !start scale index 
        call setint('edbg',idbug)        !debug printout
C------------------------------------------------------------------------

        call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's

       OPEN (UNIT = 10, FILE = file_name8, STATUS = 'Unknown')
         write(10 ,*) 'Z ', 'U ', 'UB ', 'D ', 'DB ', 'S ', 'SB ', 'GL'
            
        Do i=1,47
            zz = array(i)
            Call Nuclear_TMDFF(zz,mub,AA,fu,fub,fd,fdb,fs,fsb,fg,nloops)
                
            write(10, *), zz,fu,fub,fd,fdb,fs,fsb,fg
                
        End do

        CLOSE (UNIT = 10)

        STOP
        END PROGRAM main

