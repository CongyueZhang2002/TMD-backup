      PROGRAM nFF_test

            implicit none
*
      integer ilha
      double precision Q0,Q
      double precision Q02,Q2
      double precision AlphaQCD,AlphaQED
      double precision xPDF,xgamma
      double precision eps
      double precision xlha(18)
      double precision H_AA

      COMMON /HMASS/ H_AA

      parameter(eps=1d-10)
      data xlha / 5d-2,1d-1,1.5d-1,2d-1,2.5d-1,3d-1,3.5d-1,4d-1,4.5d-1,
     1            5d-1,5.5d-1,6d-1,6.5d-1,7d-1,7.5d-1,8d-1,8.5d-1,9d-1 /
*
*     Some examples ...
*
c      call SetFastEvolution(.true.)
c      call LockGrids(.true.)
       call SetTimeLikeEvolution(.true.)
c      call EnableEvolutionOperator(.true.)
c      call SetFFNS(4)
      call SetVFNS
      call SetNPVFNS(.true.)
c      call SetTheory("QavDS")
c      call SetTheory("QUniD")
       call SetPerturbativeOrder(1)
      call SetPDFEvolution("exactalpha")
c      call SetPDFSet("NNPDF23_nlo_as_0119_qed")
c      call SetPDFSet("MRST2004qed")
       call SetPDFSet("external")
c      call SetNumberOfGrids(1)
c      call SetGridParameters(1,30,3,1d-5)
c      call SetGridParameters(2,30,3,2d-1)
c      call SetGridParameters(3,30,3,8d-1)
c      call SetPDFSet("NNPDF30_nnlo_as_0118")
c      call SetAlphaQCDRef(0.118d0,91.2d0)
c      call SetAlphaEvolution("expanded")
c      call SetPDFEvolution("expandalpha")
c      call SetPoleMasses(1.43d0,4.3d0,173.03d0)
       call SetMSbarMasses(1.43d0,4.3d0,173.03d0)
c      call SetMaxFlavourPDFs(4)
c      call SetMaxFlavourAlpha(4)
*
*     Initializes integrals on the grids
*
      call InitializeAPFEL
*
*     Evolve PDFs on the grids
*
*      write(6,*) "Enter initial and final scale in GeV^2"
*      read(5,*) 

       Q02 = 1d0 

       Q2 =  1d0 
       Q2  =1.570000001d0**2
      Q2  =1000d0
*
      Q0 = dsqrt(Q02) !- eps
      Q  = dsqrt(Q2)

      H_AA = 1d0
      call EvolveAPFEL(Q0,Q)
*
*     Tabulate PDFs for the LHA x values
*
      write(6,*) "alpha_QCD(mu2F) =",AlphaQCD(Q)
      write(6,*) "alpha_QED(mu2F) =",AlphaQED(Q)
      write(6,*) "  "
*
      write(6,*) "Standard evolution:"
      write(6,'(a5,2a12,a12,a9,2a13,2a14)') "x",
     1         "u+ubar","d+dbar","ubar","c","b","gluon","photon"
      do ilha=1,18
         write(6,'(es7.1,7es12.4)') 
     1         xlha(ilha),
     2         xPDF(2,xlha(ilha)) + xPDF(-2,xlha(ilha)),
     3         xPDF(1,xlha(ilha)) + xPDF(-1,xlha(ilha)),
     4         xPDF(-2,xlha(ilha)),
     5         xPDF(4,xlha(ilha)) ,
     6         xPDF(5,xlha(ilha)) ,
     7         xPDF(0,xlha(ilha)),
     8         xgamma(xlha(ilha))
      enddo
      write(*,*) "  "
*
      end







************************************************************************
*
*     Collinear nuclear FF parametrization accrding to:
*
*     baseline according to DEHSS arXiv: 1410.6027
*
*     nuclear parametrization according to Pia Zurita arXiv: 2101.01088
*
************************************************************************

      subroutine ExternalSetAPFEL(x,Q,xpdf)
*
      implicit none
**
*     Input Variables
*
      double precision x, A, Q
      double precision H_AA

**
*     Input COMMON BLOCKS
*
      COMMON /HMASS/ H_AA 

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

**
*     Output Variables
*
      double precision xpdf(-6:7)
*
*     Parameters of Baseline DEHSS
*

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


      A = H_AA  

c      print*, "nuclear mass: ", A 

      

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
          
*
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
*     Initialize PDFs to zero
*
      do ipdf=-6,7
         xpdf(ipdf) = 0d0
      enddo
*
      if(x.gt.1d0) return
*
      xpdf(5)  = x*b
      xpdf(4)  = x*cb
      xpdf(3)  = x*s
      xpdf(2)  = x*u
      xpdf(1)  = x*d
      xpdf(0)  = x*g
      xpdf(-1) = x*db
      xpdf(-2) = x*ub
      xpdf(-3) = x*sb
      xpdf(-4) = x*cb
      xpdf(-5) = x*bb
*
C      print*,"x=",x
C      print*,"utot",utot
C      print*,"u=",u
C      print*,"ub=",ub

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