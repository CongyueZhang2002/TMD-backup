C---------------------------------------------------------------------
      program exampleFF
C---------------------------------------------------------------------
C------- Basic QCDNUM example job
C------- All calculations are performed in the MS bar scheme
C     ----------------------------------------------------------------
      IMPLICIT NONE 
      DOUBLE PRECISION :: Z, Q2, U, UB, D, DB, S, SB, C, CB, B, BB, GL
      DOUBLE PRECISION :: array(47), def(-6:6,12), qq(2),wt(2)
      DOUBLE PRECISION :: pdf(-6:6) 
      INTEGER :: INIT, i
      COMMON / INITIATE / INIT
      Double precision as0, r20, xmin, eps 
      integer :: iord, nfin, itype, iset, jset, iosp, nx
      integer :: nxin, nqin, lun, idbug, iqc, iqb, iq0, nq 
      double precision :: q2c, q2b, q0 
      double precision Qf 
      double precision Qf2  

C------------- X grid and mu2 grid paramaters
      data idbug/0/                                     !debug printpout
      data xmin/1.D-4/, nxin/100/                       !x grid, 
      data qq/1.D0,1.D4/, wt/1.D0,1.D0/, nqin/60/       !mu2 grid


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
     + 52*0.    /
                                       !pdfout

C--   Weight files
      character*26 fnam(3)

      data fnam /'weights/unpolarised.wgt',
     +           'weights/polarised.wgt  ',
     +           'weights/timelike.wgt   '/


C----- Set-up -----------------------------------------------------------
      lun = 6                                   ! -6 supresses banner
      call qcinit(lun,' ')                       ! initialize
      iosp = 3 ! 2 for linear interpolarion, 3 for spline
      call gxmake(xmin,1,1,nxin, nx,  iosp)                         !x-grid
C----------- generates a logarithmic grid in x with exactly 100 points in the range :
C----------- xmin(1) = 1.D-4 to 1
      call gqmake(qq,wt,2,nqin,nq)                                !mu2-grid
C----------- generates a logarithmic grid in muF**2 on which the parton densities are evolved
C----------- qq(1) --> lower and upper end of the grid.
C----------- qq(1) should always be above 0.1 GeV^2 
c----------  qq(n) is the  upper end of the grid
C----------  If n >2 then additional points in qarr are put in the grid 
C----------- In this way you can incorporate a set of starting scales or thresholds

      itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
      call wtfile(itype,fnam(itype))   !calculate weights
C------------ Reads table of type 1 2 or 3 from a diskfile 
C------------ Upon failure, computes the tables from scratch and dump them onto the file 
C------------ Routine checks that the weights on an existing input file match with itype 

C----- Perturbative Order 
      iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
      call setord(iord)                                 

C------- Set Alphas as0 at renormalization scale r20 
      as0 = 0.364d0 
      r20 = 2.0d0 
      call setalf(as0,r20)                              !input alphas

C----- Initial scale and threshhold values 
      q0 = 1.0D0 
      q2c = (1.43d0)**2d0 
      q2b = (4.3d0)**2d0 
C------- iqc, iqb --> Calculate the index for the threshold. 
      iqc  = iqfrmq(q2c)    !Charm threshold
      iqb  = iqfrmq(q2b)    !Bottom threshold
C----------------------------------------------------------
C------  Set to VFNS scheme (variable flavour number scheme --> precisely what we want)
      nfin = 1 ! Set to 1 for nonzero but frozen distributions for c,b below threshhold
      ! nfin = 0 ! Set to 0 for zero distributions for c,b below threhshhold
      call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS


C------Evolution --------------------------------------------------------
      iset = 1 
      jset = 10*iset+itype            !output pdf set and evolution type
      iq0  = iqfrmq(q0)                !start scale index 
      call setint('edbg',idbug)        !debug printout


C--------- THIS IS THE CODE which performs the Evolution
      call evolfg(jset,func,def,iq0, eps)  !evolve all pdf's


C--   Get results   ----------------------------------------------------
C---- q2mine

      Qf = 1.001d0 
      Qf2 = Qf**2 


      OPEN (UNIT = 10, FILE = 'newFF_1.csv', STATUS = 'REPLACE')
        write(10 ,*) 'Z ', 'U ', 'UB ', 'D ', 'DB ', 'S ', 'SB ', 'C ',
     9        'CB ', 'B ', 'BB ', 'GL '
        INIT = 0        
        Do i=1,47
            Z = array(i)
            call allfxq(iset,Z,Qf2,pdf,0,1) 
            U   = pdf(2)/Z
            UB  = pdf(-2)/Z 
            D   = pdf(1)/Z 
            Db  = pdf(-1)/Z 
            S   = pdf(3)/Z 
            Sb  = pdf(-3)/Z 
            C   = pdf(4)/Z
            CB  = pdf(4)/Z
            B   = pdf(5)/Z
            BB  = pdf(-5)/Z
            GL  = pdf(0)/Z  
        write(10, *), Z, U, UB, D, DB, S, SB, C, CB, B, BB, GL 
            
        End do

        CLOSE (UNIT = 10)



      OPEN (UNIT = 20, FILE = 'LIKEn_1.csv', STATUS = 'REPLACE')
       write(20 ,*) 'Z ', 'U ', 'UB ', 'D ', 'DB ', 'S ', 'SB ', 'C ',
     9        'CB ', 'B ', 'BB ', 'GL '
       INIT = 0
        Do i=1,47
            Z = array(i)
            call LIKEn(1,1,1,0,array(i),Qf2,1,U,UB,D,DB,S,SB,C,CB,
     9      B,BB,GL)
            
       write(20, *), Z, U/Z, UB/Z, D/Z, DB/Z, S/Z, SB/Z, 
     9                      C/Z, CB/Z, B/Z, BB/Z, GL/Z 
            
        End do

        CLOSE(UNIT = 20)



C--   Get results   ----------------------------------------------------
C---- q2mine

      Qf = 2.00d0 
      Qf2 = Qf**2 


      OPEN (UNIT = 10, FILE = 'newFF_mc.csv', STATUS = 'REPLACE')
        write(10 ,*) 'Z ', 'U ', 'UB ', 'D ', 'DB ', 'S ', 'SB ', 'C ',
     9        'CB ', 'B ', 'BB ', 'GL '
        INIT = 0        
        Do i=1,47
            Z = array(i)
            call allfxq(iset,Z,Qf2,pdf,0,1) 
            U   = pdf(2)/Z
            UB  = pdf(-2)/Z 
            D   = pdf(1)/Z 
            Db  = pdf(-1)/Z 
            S   = pdf(3)/Z 
            Sb  = pdf(-3)/Z 
            C   = pdf(4)/Z
            CB  = pdf(4)/Z
            B   = pdf(5)/Z
            BB  = pdf(-5)/Z
            GL  = pdf(0)/Z  
        write(10, *), Z, U, UB, D, DB, S, SB, C, CB, B, BB, GL 
            
        End do

        CLOSE (UNIT = 10)



      OPEN (UNIT = 20, FILE = 'LIKEn_mc.csv', STATUS = 'REPLACE')
       write(20 ,*) 'Z ', 'U ', 'UB ', 'D ', 'DB ', 'S ', 'SB ', 'C ',
     9        'CB ', 'B ', 'BB ', 'GL '
       INIT = 0
        Do i=1,47
            Z = array(i)
            call LIKEn(1,1,1,0,array(i),Qf2,1,U,UB,D,DB,S,SB,C,CB,
     9      B,BB,GL)
            
       write(20, *), Z, U/Z, UB/Z, D/Z, DB/Z, S/Z, SB/Z, 
     9                      C/Z, CB/Z, B/Z, BB/Z, GL/Z 
            
        End do

        CLOSE(UNIT = 20)
       


      end
      
C     ----------------------------------------------------------------

C     ======================================
      double precision function func(ipdf,x)
C     ======================================

      implicit none
**
*     Input Variables
*
      double precision x, A, Q
    
**
*     Internal Variables
*
      integer ipdf, i
      double precision EBETA
      double precision utot,dtot,stot,ctot,btot,g
      double precision u,ub,d,db,s,sb,c,cb,b,bb
      double precision Ni(0:6),ai(0:6),bi(0:6),gi(0:6),di(0:6)
      double precision Nq1,Nq2,Ng1,Ng2 
      double precision aq1,bq1,gq1,dq1,ag1,bg1,gg1,dg1
      double precision aq2,bq2,gq2,dq2,ag2,bg2,gg2,dg2


*     Parameters of Baseline DEHSS
*
      A=1.0d0

*     utot
      Ni(1)  = 0.387d0 !/ 0.682937d0 ! Ni/(Beta[2 + ai, bi + 1] + gi Beta[2 + ai, bi + di + 1]) 
      ai(1)  = -0.388d0
      bi(1)  = 0.910d0
      gi(1)  = 7.15d0
      di(1)  = 3.96d0
C      Ni(1)  = Ni(1) / (EBETA(2.d0+ai(1), bi(1)+1.d0) 
C     1              + gi(1) * EBETA(2.d0+ai(1), bi(1)+di(1)+1.d0))


*     dtot
      Ni(2)  = 0.388d0 !/ 0.682937d0
      ai(2)  = ai(1)
      bi(2)  = bi(1)
      gi(2)  = gi(1)
      di(2)  = di(1)
C      Ni(2)  = Ni(2) / (EBETA(2.d0+ai(2), bi(2)+1.d0) 
C     1              + gi(2) * EBETA(2.d0+ai(2), bi(2)+di(2)+1.d0))

*     ubar = d
      Ni(3) = 0.105d0 !/ 0.0499286d0
      ai(3) = 1.649d0
      bi(3) = 3.286d0
      gi(3) = 49.95d0
      di(3) =  8.67d0
C      Ni(3)  = Ni(3) / (EBETA(2.d0+ai(3), bi(3)+1.d0) 
C     1              + gi(3) * EBETA(2.d0+ai(3), bi(3)+di(3)+1.d0))

*     stot
      Ni(4)  = 0.273d0 !/ 0.100144d0
      ai(4)  = 1.449d0
      bi(4)  = bi(3)
      gi(4)  = gi(3)
      di(4)  = di(3)
C     Ni(4)  = Ni(4) / (EBETA(2.d0+ai(4), bi(4)+1.d0) 
C     1              + gi(4) * EBETA(2.d0+ai(4), bi(4)+di(4)+1.d0))

*     ctot
      Ni(5)  = 0.306d0 !/ 0.0228321d0
      ai(5)  = 1.345d0
      bi(5)  = 5.519d0
      gi(5)  = 19.78d0
      di(5)  = 10.22d0
C      Ni(5)  = Ni(5) / (EBETA(2.d0+ai(5), bi(5)+1.d0) 
C     1              + gi(5) * EBETA(2.d0+ai(5), bi(5)+di(5)+1.d0))

*     btot
      Ni(6)  =  0.372d0 !/ 0.176549d0
      ai(6)  = -0.127d0 
      bi(6)  =  4.490d0 
      gi(6)  =  24.49d0 
      di(6)  =  12.80d0
C      Ni(6)  = Ni(6) / (EBETA(2.d0+ai(6), bi(6)+1.d0) 
C     1              + gi(6) * EBETA(2.d0+ai(6), bi(6)+di(6)+1.d0))

*     g
      Ni(0)   =  0.260d0 ! /0.00832656
      ai(0)   =  2.552d0
      bi(0)   =  6.194d0
      gi(0)   =  87.06d0
      di(0)   =  20.36d0  
C      Ni(0)  = Ni(0) / (EBETA(2.d0+ai(0), bi(0)+1.d0) 
C     1              + gi(0) * EBETA(2.d0+ai(0), bi(0)+di(0)+1.d0))
       
*
*     Pia Zurita Nuclear modification
*

       Ng1 = -0.0262d0  !param(1)
       Nq1 = 0.0322d0   !param(2) 
             
       Ng2 = 0.6654d0   !param(3)
       Nq2 = 0.4567d0   !param(4)

       ag1 = 0.d0
       ag2 = ag1

       aq1 =  ag1
       aq2 =  ag2

       bg1 = -0.0148d0 !param(5)
       bq1 =  bg1

       bg2 =  Ng2 
       bq2 =  Nq2

       gg1 = -0.1555d0 !param(6)
       gq1 = gg1

       gg2 = Ng2
       gq2 = Nq2

       dg1 = -0.0451d0 !param(7)
       dq1 = dg1

       dg2 = Ng2
       dq2 = Nq2


      A = 1.0d0  

C      print*, "nuclear mass: ", A 

      

      DO i=1,6

            Ni(i) = Ni(i) * (1.d0 +  Nq1*(1d0-A**Nq2))
            ai(i) = ai(i) + aq1 *(1d0-A**aq2)
            bi(i) = bi(i) + bq1 *(1d0-A**bq2)
            gi(i) = gi(i) + gq1 *(1d0-A**gq2)
            di(i) = di(i) + dq1 *(1d0-A**dq2)
            Ni(i)  = Ni(i) / (EBETA(2.d0+ai(i), bi(i)+1.d0) 
     1              + gi(i) * EBETA(2.d0+ai(i), bi(i)+di(i)+1.d0))
          
      ENDDO


            Ni(0) = Ni(0) * (1.d0 + Ng1*(1d0-A**Ng2))
            ai(0) = ai(0) + ag1 *(1d0-A**ag2)
            bi(0) = bi(0) + bg1 *(1d0-A**bg2)
            gi(0) = gi(0) + gg1 *(1d0-A**gg2)
            di(0) = di(0) + dg1 *(1d0-A**dg2)
            Ni(0) = Ni(0) / (EBETA(2.d0+ai(0), bi(0)+1.d0) 
     1              + gi(0) * EBETA(2.d0+ai(0), bi(0)+di(0)+1.d0))
          
*     User defined PDFs
*
      utot  = Ni(1)*x**ai(1)*(1-x)**bi(1)*(1.d0+gi(1)*(1-x)**di(1))
      dtot  = Ni(2)*x**ai(2)*(1-x)**bi(2)*(1.d0+gi(2)*(1-x)**di(2))
      ub    = Ni(3)*x**ai(3)*(1-x)**bi(3)*(1.d0+gi(3)*(1-x)**di(3))
      d     = ub
      stot  = Ni(4)*x**ai(4)*(1-x)**bi(4)*(1.d0+gi(4)*(1-x)**di(4))
      ctot  = Ni(5)*x**ai(5)*(1-x)**bi(5)*(1.d0+gi(5)*(1-x)**di(5))
      btot  = Ni(6)*x**ai(6)*(1-x)**bi(6)*(1.d0+gi(6)*(1-x)**di(6))
      g     = Ni(0)*x**ai(0)*(1-x)**bi(0)*(1.d0+gi(0)*(1-x)**di(0))


      u   = utot - ub
      db  = dtot - d

      s   = stot/2.d0
      sb  = s
      c   = ctot/2.d0
      cb  = c
      b   = btot/2.d0
      bb  = b

*
                     func = 0.D0
      if(ipdf.eq. 0) func = x*g
      if(ipdf.eq. 1) func = x*d
      if(ipdf.eq. 2) func = x*u
      if(ipdf.eq. 3) func = x*s 
      if(ipdf.eq. 4) func = x*db 
      if(ipdf.eq. 5) func = x*ub 
      if(ipdf.eq. 6) func = x*sb   
      if(ipdf.eq. 7) func = x*c
      if(ipdf.eq. 8) func = x*cb 
      if(ipdf.eq. 9) func = 0.d0
      if(ipdf.eq.10) func = 0.d0  
      if(ipdf.eq.11) func = 0.D0 
      if(ipdf.eq.12) func = 0.D0

      return
      end






      DOUBLE PRECISION FUNCTION EBETA(X,Z)

      implicit none
      double precision x, z , y

      Y = X + Z

      EBETA = DGAMMA(X) * DGAMMA(Z) / DGAMMA (Y)

      return 
      end



CS    REAL FUNCTION GAMMA(X)
      DOUBLE PRECISION FUNCTION DGAMMA(X)
C----------------------------------------------------------------------
C
C This routine calculates the GAMMA function for a real argument X.
C   Computation is based on an algorithm outlined in reference 1.
C   The program uses rational functions that approximate the GAMMA
C   function to at least 20 significant decimal digits.  Coefficients
C   for the approximation over the interval (1,2) are unpublished.
C   Those for the approximation for X .GE. 12 are from reference 2.
C   The accuracy achieved depends on the arithmetic system, the
C   compiler, the intrinsic functions, and proper selection of the
C   machine-dependent constants.
C
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants
C
C beta   - radix for the floating-point representation
C maxexp - the smallest positive power of beta that overflows
C XBIG   - the largest argument for which GAMMA(X) is representable
C          in the machine, i.e., the solution to the equation
C                  GAMMA(XBIG) = beta**maxexp
C XINF   - the largest machine representable floating-point number;
C          approximately beta**maxexp
C EPS    - the smallest positive floating-point number such that
C          1.0+EPS .GT. 1.0
C XMININ - the smallest positive floating-point number such that
C          1/XMININ is machine representable
C
C     Approximate values for some important machines are:
C
C                            beta       maxexp        XBIG
C
C CRAY-1         (S.P.)        2         8191        966.961
C Cyber 180/855
C   under NOS    (S.P.)        2         1070        177.803
C IEEE (IBM/XT,
C   SUN, etc.)   (S.P.)        2          128        35.040
C IEEE (IBM/XT,
C   SUN, etc.)   (D.P.)        2         1024        171.624
C IBM 3033       (D.P.)       16           63        57.574
C VAX D-Format   (D.P.)        2          127        34.844
C VAX G-Format   (D.P.)        2         1023        171.489
C
C                            XINF         EPS        XMININ
C
C CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
C Cyber 180/855
C   under NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
C IEEE (IBM/XT,
C   SUN, etc.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
C IEEE (IBM/XT,
C   SUN, etc.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
C IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
C VAX D-Format   (D.P.)   1.70D+38     1.39D-17    5.88D-39
C VAX G-Format   (D.P.)   8.98D+307    1.11D-16    1.12D-308
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  The program returns the value XINF for singularities or
C     when overflow would occur.  The computation is believed
C     to be free of underflow and overflow.
C
C
C  Intrinsic functions required are:
C
C     INT, DBLE, EXP, LOG, REAL, SIN
C
C
C References: "An Overview of Software Development for Special
C              Functions", W. J. Cody, Lecture Notes in Mathematics,
C              506, Numerical Analysis Dundee, 1975, G. A. Watson
C              (ed.), Springer Verlag, Berlin, 1976.
C
C              Computer Approximations, Hart, Et. Al., Wiley and
C              sons, New York, 1968.
C
C  Latest modification: October 12, 1989
C
C  Authors: W. J. Cody and L. Stoltz
C           Applied Mathematics Division
C           Argonne National Laboratory
C           Argonne, IL 60439
C
C----------------------------------------------------------------------
      INTEGER I,N
      LOGICAL PARITY
CS    REAL 
      DOUBLE PRECISION 
     1    C,CONV,EPS,FACT,HALF,ONE,P,PI,Q,RES,SQRTPI,SUM,TWELVE,
     2    TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      DIMENSION C(7),P(8),Q(8)
C----------------------------------------------------------------------
C  Mathematical constants
C----------------------------------------------------------------------
CS    DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/,
CS   1     SQRTPI/0.9189385332046727417803297E0/,
CS   2     PI/3.1415926535897932384626434E0/
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0D0,0.5D0,12.0D0,2.0D0,0.0D0/,
     1     SQRTPI/0.9189385332046727417803297D0/,
     2     PI/3.1415926535897932384626434D0/
C----------------------------------------------------------------------
C  Machine dependent parameters
C----------------------------------------------------------------------
CS    DATA XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/,
CS   1     XINF/3.4E38/
      DATA XBIG,XMININ,EPS/171.624D0,2.23D-308,2.22D-16/,
     1     XINF/1.79D308/
C----------------------------------------------------------------------
C  Numerator and denominator coefficients for rational minimax
C     approximation over (1,2).
C----------------------------------------------------------------------
CS    DATA P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,
CS   1       -3.79804256470945635097577E+2,6.29331155312818442661052E+2,
CS   2       8.66966202790413211295064E+2,-3.14512729688483675254357E+4,
CS   3       -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
CS    DATA Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,
CS   1      -1.01515636749021914166146E+3,-3.10777167157231109440444E+3,
CS   2        2.25381184209801510330112E+4,4.75584627752788110767815E+3,
CS   3      -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
      DATA P/-1.71618513886549492533811D+0,2.47656508055759199108314D+1,
     1       -3.79804256470945635097577D+2,6.29331155312818442661052D+2,
     2       8.66966202790413211295064D+2,-3.14512729688483675254357D+4,
     3       -3.61444134186911729807069D+4,6.64561438202405440627855D+4/
      DATA Q/-3.08402300119738975254353D+1,3.15350626979604161529144D+2,
     1      -1.01515636749021914166146D+3,-3.10777167157231109440444D+3,
     2        2.25381184209801510330112D+4,4.75584627752788110767815D+3,
     3      -1.34659959864969306392456D+5,-1.15132259675553483497211D+5/
C----------------------------------------------------------------------
C  Coefficients for minimax approximation over (12, INF).
C----------------------------------------------------------------------
CS    DATA C/-1.910444077728E-03,8.4171387781295E-04,
CS   1     -5.952379913043012E-04,7.93650793500350248E-04,
CS   2     -2.777777777777681622553E-03,8.333333333333333331554247E-02,
CS   3      5.7083835261E-03/
      DATA C/-1.910444077728D-03,8.4171387781295D-04,
     1     -5.952379913043012D-04,7.93650793500350248D-04,
     2     -2.777777777777681622553D-03,8.333333333333333331554247D-02,
     3      5.7083835261D-03/
C----------------------------------------------------------------------
C  Statement functions for conversion between integer and float
C----------------------------------------------------------------------
CS    CONV(I) = REAL(I)
      CONV(I) = DBLE(I)
      PARITY = .FALSE.
      FACT = ONE
      N = 0
      Y = X
      IF (Y .LE. ZERO) THEN
C----------------------------------------------------------------------
C  Argument is negative
C----------------------------------------------------------------------
            Y = -X
            Y1 = AINT(Y)
            RES = Y - Y1
            IF (RES .NE. ZERO) THEN
                  IF (Y1 .NE. AINT(Y1*HALF)*TWO) PARITY = .TRUE.
                  FACT = -PI / SIN(PI*RES)
                  Y = Y + ONE
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
      END IF
C----------------------------------------------------------------------
C  Argument is positive
C----------------------------------------------------------------------
      IF (Y .LT. EPS) THEN
C----------------------------------------------------------------------
C  Argument .LT. EPS
C----------------------------------------------------------------------
            IF (Y .GE. XMININ) THEN
                  RES = ONE / Y
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
         ELSE IF (Y .LT. TWELVE) THEN
            Y1 = Y
            IF (Y .LT. ONE) THEN
C----------------------------------------------------------------------
C  0.0 .LT. argument .LT. 1.0
C----------------------------------------------------------------------
                  Z = Y
                  Y = Y + ONE
               ELSE
C----------------------------------------------------------------------
C  1.0 .LT. argument .LT. 12.0, reduce argument if necessary
C----------------------------------------------------------------------
                  N = INT(Y) - 1
                  Y = Y - CONV(N)
                  Z = Y - ONE
            END IF
C----------------------------------------------------------------------
C  Evaluate approximation for 1.0 .LT. argument .LT. 2.0
C----------------------------------------------------------------------
            XNUM = ZERO
            XDEN = ONE
            DO 260 I = 1, 8
               XNUM = (XNUM + P(I)) * Z
               XDEN = XDEN * Z + Q(I)
  260       CONTINUE
            RES = XNUM / XDEN + ONE
            IF (Y1 .LT. Y) THEN
C----------------------------------------------------------------------
C  Adjust result for case  0.0 .LT. argument .LT. 1.0
C----------------------------------------------------------------------
                  RES = RES / Y1
               ELSE IF (Y1 .GT. Y) THEN
C----------------------------------------------------------------------
C  Adjust result for case  2.0 .LT. argument .LT. 12.0
C----------------------------------------------------------------------
                  DO 290 I = 1, N
                     RES = RES * Y
                     Y = Y + ONE
  290             CONTINUE
            END IF
         ELSE
C----------------------------------------------------------------------
C  Evaluate for argument .GE. 12.0,
C----------------------------------------------------------------------
            IF (Y .LE. XBIG) THEN
                  YSQ = Y * Y
                  SUM = C(7)
                  DO 350 I = 1, 6
                     SUM = SUM / YSQ + C(I)
  350             CONTINUE
                  SUM = SUM/Y - Y + SQRTPI
                  SUM = SUM + (Y-HALF)*LOG(Y)
                  RES = EXP(SUM)
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
      END IF
C----------------------------------------------------------------------
C  Final adjustments and return
C----------------------------------------------------------------------
      IF (PARITY) RES = -RES
      IF (FACT .NE. ONE) RES = FACT / RES
CS900 GAMMA = RES
  900 DGAMMA = RES
      RETURN
C ---------- Last line of GAMMA ----------
      END



