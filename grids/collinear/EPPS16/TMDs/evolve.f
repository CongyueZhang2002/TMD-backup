c--------------------------------------------------------------
c     alphas(q) -- NLO strong coupling constant (from CTEQ)
c--------------------------------------------------------------

      function alphas(q)
      implicit none
      real*8 q,q2,lambda,lambda2,alphas,b0,b1,tt,ktw,ptw
      real*8 c0,c02,pi,CF,CA,mb,alam4,alam5,bmax,kappa2,Qini
      integer nf
      common /cons/ c0,c02,pi,CF,CA,mb,alam4,alam5,
     >     bmax,ktw,ptw,Qini,kappa2

      if (q.le.mb) then
         nf=4
         lambda=alam4
      else
         nf=5
         lambda=alam5
      endif
      b0=11d0-2d0/3d0*nf
      b1=102d0-38d0/3d0*nf

      q2=q*q
      lambda2=lambda*lambda
      tt=dlog(q2/lambda2)
      alphas=4d0*pi/(b0*tt)*(1d0-b1/(b0*b0)*dlog(tt)/tt)
      return
      end

c--------------------------------------------------------------------------
c               / xf
c               |
c       value = | dx func(x)         n > 1, number of divisions for xf-xi
c               |
c              / xi
c--------------------------------------------------------------------------

      function qgauss(func,xi,xf,n)
      implicit none
      real*8 qgauss,func,xi,xf,xn,value,x1,x2,val
      integer i,n
      external func

      if(n.le.1) then           ! same as n=1
         x1=xi
         x2=xf
         call gauss(func,x1,x2,val)
         qgauss=val
         return
      endif

      xn=(xf-xi)/float(n)
      value=0d0
      x2=xi
      Do 100 i=1,n
         x1=x2
         x2=x1+xn
         call gauss(func,x1,x2,val)
         value=value+val
 100  continue

      qgauss=value
      return
      end

c-------------------------------------------------
      subroutine gauss(func,xi,xf,value)
      implicit none
      real*8 func,xi,xf,value,xm,xr,dx,x(8),w(8)
      real*8 eps
      integer j
      data eps /1.0d-25/
      data w
     1   / 0.02715 24594 11754 09485 17805 725D0,
     2     0.06225 35239 38647 89286 28438 370D0,
     3     0.09515 85116 82492 78480 99251 076D0,
     4     0.12462 89712 55533 87205 24762 822D0,
     5     0.14959 59888 16576 73208 15017 305D0,
     6     0.16915 65193 95002 53818 93120 790D0,
     7     0.18260 34150 44923 58886 67636 680D0,
     8     0.18945 06104 55068 49628 53967 232D0 /
      DATA X
     1   / 0.98940 09349 91649 93259 61541 735D0,
     2     0.94457 50230 73232 57607 79884 155D0,
     3     0.86563 12023 87831 74388 04678 977D0,
     4     0.75540 44083 55003 03389 51011 948D0,
     5     0.61787 62444 02643 74844 66717 640D0,
     6     0.45801 67776 57227 38634 24194 430D0,
     7     0.28160 35507 79258 91323 04605 015D0,
     8     0.09501 25098 37637 44018 53193 354D0 /

      xm=0.5d0*(xf+xi)
      xr=0.5d0*(xf-xi)
      if (abs(xr).lt.eps) print *,
     >     'WARNING: Too high accuracy required for QGAUSS!'

      value=0d0

      Do 100 j=1,8
         dx=xr*x(j)
         value=value+w(j)*(func(xm+dx)+func(xm-dx))
 100  continue

      value=xr*value
      return
      end

      BLOCK DATA
      IMPLICIT DOUBLE PRECISION(A-Z)
      real*8 bdiv
      integer icon,IOpdf,FINI,nbel,nbuu,nbut,nY,nQ,IH,IC,IO,lhapdf
      common /fitcon/ icon
      common /pdforder/ IOpdf
      common /fforder/ IO
      common /lhappdf/ lhapdf
      COMMON / FRAGINI / FINI
c-----icon=1(no evo),2(evo), IOpdf(PDF)/IO(FF)=0(LO),1(NLO)
      data icon/2/IOpdf/1/IO/1/FINI/0/
      common /ngbesl/ bdiv,nbel,nbuu,nbut,nQ,nY
      data nbel/4/nbuu/4/nbut/6/nY/2/nQ/2/bdiv/10d0/
      common /cons/ c0,c02,pi,CF,CA,mb,alam4,alam5,
     >     bmax,ktw,ptw,Qini,kappa2
      common /Wpar/ Vud,Vus,Mw,GF,sw,Vzu,Vzds,Azu,Azds,Mz,Nc,convert
      common /nucleus/ AA,ZZ
      data AA/64d0/ZZ/29d0/
      data Vud/0.97425d0/Vus/0.2252/Mw/80.385d0/
c-----Vzu = 1/2 - 4/3*sin^2(theta_w)
c     Vzds=-1/2 + 2/3*sin^2(theta_w)
c     Azu = 1/2, Azds = -1/2
      data Vzu/0.191787d0/Vzds/-0.345893d0/
      data Azu/0.5d0/Azds/-0.5d0/Mz/91.1876d0/
      data GF/1.1663787d-5/Nc/3d0/convert/0.3894d9/
      data c0/1.122919d0/c02/1.260947d0/ ! c02=c0^2, c0=2*exp(-\gamma_E)
      data pi/3.141592653589793d0/CF/1.333333d0/CA/3d0/mb/4.5d0/
      data alam4/0.326d0/alam5/0.226d0/
      data bmax/1.5d0/
      data Qini/1.54919d0/
      data ktw/0.424d0/ptw/0.168d0/kappa2/0.88d0/
      !data ktw/0.38d0/ptw/0.19d0/kappa2/0.08d0/

c----   NEW

      real*8 lambda3(3)    !lambda(order) for nf=3 and order=1,2,3
      real*8 lambda4(3)    !lambda(order) for nf=4 and order=1,2,3
      real*8 lambda5(3)    !lambda(order) for nf=5 and order=1,2,3
      common /lambda/lambda3,lambda4,lambda5
      data lambda3/0.266d0, 0.266d0, 0.266d0/
      data lambda4/0.246d0, 0.246d0, 0.246d0/
      data lambda5/0.224d0, 0.224d0, 0.224d0/

      real*8 Beta3(3)    !Beta(order) for nf=3 and order=1,2,3
      real*8 Beta4(3)    !Beta(order) for nf=4 and order=1,2,3
      real*8 Beta5(3)    !Beta(order) for nf=5 and order=1,2,3
      common /Beta/Beta3,Beta4,Beta5
      data Beta3/9.000000000000000d0, 64.000000000000000d0,
     >     643.8333333333334d0/
      data Beta4/8.333333333333333d0, 51.333333333333336d0,
     >     406.35185185185185d0/
      data Beta5/7.666666666666667d0, 38.666666666666664d0,
     >     180.90740740740742d0/

cccccc We use Padé approximation for cusp at 4 loops:
cccccc CuspQnf(4) = (CuspQnf(3))^2/CuspQnf(2)
      real*8 CuspQ3(4)      !Cusp quark coefficients for nf=3
      real*8 CuspQ4(4)      !Cusp quark coefficients for nf=4
      real*8 CuspQ5(4)      !Cusp quark coefficients for nf=5
      common /CuspQ/CuspQ3,CuspQ4,CuspQ5
      data CuspQ3/5.333333333333333d0, 48.69544319419009d0,
     >     618.2248693918798d0, 7848.82d0/
      data CuspQ4/5.333333333333333d0, 42.76951726826417d0,
     >     429.5065747522099d0, 4313.26d0/
      data CuspQ5/5.333333333333333d0, 36.84359134233824d0,
     >     239.20803319895987d0, 1553.06d0/

cccccc We use Padé approximation for cusp at 4 loops:
cccccc CuspGnf(4) = (CuspGnf(3))^2/CuspGnf(2)
      real*8 CuspG3(4)      !Cusp gluon coefficients for nf=3
      real*8 CuspG4(4)      !Cusp gluon coefficients for nf=4
      real*8 CuspG5(4)      !Cusp gluon coefficients for nf=5
      common /CuspG/CuspG3,CuspG4,CuspG5
      data CuspG3/12.00000000000000d0, 109.5647471869277d0,
     >     1391.0059561317298d0, 17659.9d0/
      data CuspG4/12.00000000000000d0, 96.23141385359438d0,
     >     966.3897931924722d0, 9704.83d0/
      data CuspG5/12.00000000000000d0, 82.89808052026105d0,
     >     538.2180746976597d0, 3494.4d0/

      real*8 NCuspQ3(3)     !Non-cusp quark coefficients for nf=3
      real*8 NCuspQ4(3)     !Non-cusp quark coefficients for nf=4
      real*8 NCuspQ5(3)     !Non-cusp quark coefficients for nf=5
      common /NCuspQ/NCuspQ3,NCuspQ4,NCuspQ5
      data NCuspQ3/-8.00000000000000d0, -29.243530284415503d0,
     >     -738.2562930508085d0/
      data NCuspQ4/-8.00000000000000d0, -14.050795508138547d0,
     >     -491.96573445169145d0/
      data NCuspQ5/-8.00000000000000d0, 1.1419392681384102d0,
     >     -249.38756710544408d0/

      real*8 NCuspG3(3)     !Non-cusp gluon coefficients for nf=3
      real*8 NCuspG4(3)     !Non-cusp gluon coefficients for nf=4
      real*8 NCuspG5(3)     !Non-cusp gluon coefficients for nf=5
      common /NCuspG/NCuspG3,NCuspG4,NCuspG5
      data NCuspG3/-18.000000000000000d0, -227.8995118764504d0,
     >     -3957.378545284555d0/
      data NCuspG4/-16.666666666666668d0, -200.7014703660655d0,
     >     -3206.1850399832547d0/
      data NCuspG5/-15.333333333333334d0, -173.50342885568062d0,
     >     -2485.5387880396356d0/

ccccc Coefficients of the D term for quarks, for nf=3,4,5 and for different powers
ccccc of the logarithem Lperp, as in 1604.07869
ccccc DtermQnm(nf) is for aspi^n and Lperp^m
      real*8 DtermQ11(3:5)
      real*8 DtermQ10(3:5)
      real*8 DtermQ22(3:5)
      real*8 DtermQ21(3:5)
      real*8 DtermQ20(3:5)
      real*8 DtermQ33(3:5)
      real*8 DtermQ32(3:5)
      real*8 DtermQ31(3:5)
      real*8 DtermQ30(3:5)
      common /DtermQ/DtermQ11,DtermQ10,DtermQ22,DtermQ21,DtermQ20,
     >               DtermQ33,DtermQ32,DtermQ31,DtermQ30
      data DtermQ11/2.666666666666667d0, 2.666666666666667d0,
     >     2.666666666666667d0/
      data DtermQ10/0d0,0d0,0d0/
      data DtermQ22/12.00000000000000d0, 11.11111111111111d0,
     >     10.22222222222222d0/
      data DtermQ21/24.34772159709504d0, 21.38475863413208d0,
     >     18.42179567116912d0/
      data DtermQ20/-15.75963102138172d0, -18.52506312014716d0,
     >     -21.29049521891259d0/
      data DtermQ33/72.00000000000000d0, 61.72839506172840d0,
     >     52.24691358024691d0/
      data DtermQ32/304.4628277071887d0, 246.6507663955451d0,
     >     192.7893223678521d0/
      data DtermQ31/25.43907631106866d0, -93.99776462634783d0,
     >     -206.8502434238466d0/
      data DtermQ30/-260.5342569653863d0, -305.8365169792609d0,
     >     -342.0455323786138d0/

ccccc Coefficients of the D term for gluons, for nf=3,4,5 and for different powers
ccccc of the logarithem Lperp, as in 1604.07869
ccccc DtermQnm(nf) is for aspi^n and Lperp^m
      real*8 DtermG11(3:5)
      real*8 DtermG10(3:5)
      real*8 DtermG22(3:5)
      real*8 DtermG21(3:5)
      real*8 DtermG20(3:5)
      real*8 DtermG33(3:5)
      real*8 DtermG32(3:5)
      real*8 DtermG31(3:5)
      real*8 DtermG30(3:5)
      common /DtermG/DtermG11,DtermG10,DtermG22,DtermG21,DtermG20,
     >               DtermG33,DtermG32,DtermG31,DtermG30
      data DtermG11/6.000000000000000d0, 6.000000000000000d0,
     >     6.000000000000000d0/
      data DtermG10/0d0,0d0,0d0/
      data DtermG22/27.00000000000000d0, 25.00000000000000d0,
     >     23.00000000000000d0/
      data DtermG21/54.78237359346385d0, 48.11570692679718d0,
     >     41.44904026013051d0/
      data DtermG20/-35.45916979810888d0, -41.68139202033110d0,
     >     -47.90361424255332d0/
      data DtermG33/162.0000000000000d0, 138.8888888888889d0,
     >     117.5555555555556d0/
      data DtermG32/685.0413623411746d0, 554.9642243899765d0,
     >     433.7759753276673d0/
      data DtermG31/57.23792169990449d0, -211.4949704092826d0,
     >     -465.4130477036549d0/
      data DtermG30/-586.2020781721192d0, -688.1321632033369d0,
     >     -769.6024478518811d0/

ccccc Legendre-Gauss quadrature integration
ccccc https://pomax.github.io/bezierinfo/legendre-gauss.html#n5
      real*8 Xi8(8),Wi8(8)     ! Points and weights
      common /gaussparam8/Xi8,Wi8
      data Xi8/
     >    -0.1834346424956498d0, 0.1834346424956498d0,
     >    -0.5255324099163290d0, 0.5255324099163290d0,
     >    -0.7966664774136267d0, 0.7966664774136267d0,
     >    -0.9602898564975363d0, 0.9602898564975363d0 /
      data Wi8/
     >     0.3626837833783620d0, 0.3626837833783620d0,
     >     0.3137066458778873d0, 0.3137066458778873d0,
     >     0.2223810344533745d0, 0.2223810344533745d0,
     >     0.1012285362903763d0, 0.1012285362903763d0 /

ccccc Legendre-Gauss quadrature integration
ccccc https://pomax.github.io/bezierinfo/legendre-gauss.html#n5
      real*8 Xi12(12),Wi12(12)     ! Points and weights
      common /gaussparam12/Xi12,Wi12
      data Xi12/
     >    -0.1252334085114689d0, 0.1252334085114689d0,
     >    -0.3678314989981802d0, 0.3678314989981802d0,
     >    -0.5873179542866175d0, 0.5873179542866175d0,
     >    -0.7699026741943047d0, 0.7699026741943047d0,
     >    -0.9041172563704749d0, 0.9041172563704749d0,
     >    -0.9815606342467192d0, 0.9815606342467192d0 /
      data Wi12/
     >     0.2491470458134028d0, 0.2491470458134028d0,
     >     0.2334925365383548d0, 0.2334925365383548d0,
     >     0.2031674267230659d0, 0.2031674267230659d0,
     >     0.1600783285433462d0, 0.1600783285433462d0,
     >     0.1069393259953184d0, 0.1069393259953184d0,
     >     0.0471753363865118d0, 0.0471753363865118d0 /

ccccc Legendre-Gauss quadrature integration
ccccc https://pomax.github.io/bezierinfo/legendre-gauss.html#n5
      real*8 Xi16(16),Wi16(16)     ! Points and weights
      common /gaussparam16/Xi16,Wi16
      data Xi16/
     &    -0.0950125098376374d0, 0.0950125098376374d0,
     &    -0.2816035507792589d0, 0.2816035507792589d0,
     &    -0.4580167776572274d0, 0.4580167776572274d0,
     &    -0.6178762444026438d0, 0.6178762444026438d0,
     &    -0.7554044083550030d0, 0.7554044083550030d0,
     &    -0.8656312023878318d0, 0.8656312023878318d0,
     &    -0.9445750230732326d0, 0.9445750230732326d0,
     &    -0.9894009349916499d0, 0.9894009349916499d0/
      data Wi16/
     &     0.1894506104550685d0, 0.1894506104550685d0,
     &     0.1826034150449236d0, 0.1826034150449236d0,
     &     0.1691565193950025d0, 0.1691565193950025d0,
     &     0.1495959888165767d0, 0.1495959888165767d0,
     &     0.1246289712555339d0, 0.1246289712555339d0,
     &     0.0951585116824928d0, 0.0951585116824928d0,
     &     0.0622535239386479d0, 0.0622535239386479d0,
     &     0.0271524594117541d0, 0.0271524594117541d0/

      end
