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

c----------------------------------------------------------------------
c     alpha_qed(q) -- QED coupling constant
c     eq23=e_u^2+e_d^2+e_s^2
c     eq24=eq23+e_c^2
c     eq25=eq24+e_b^2
c----------------------------------------------------------------------
      function alpha_qed(Q)
      implicit none
      real*8 alpha_qed,Q,mc,mtau,mb,mz,b0,eq23,eq24,eq25
      real*8 aem_mc,aem_mtau,aem_mb,aem_mz,aem,pi
      integer nl
      data mc,mtau,mb,mz /1.3d0,1.777d0,4.7d0,91.1876d0/
      data aem_mc   /0.007476296864d0/
      data aem_mtau /0.007496122052d0/
      data aem_mb   /0.007570302837d0/
      data aem_mz   /0.007874015748d0/
      eq23=6d0/9d0
      eq24=10d0/9d0
      eq25=11d0/9d0
      pi=4d0*datan(1d0)
      if(Q.le.mc) then
         nl=2
         b0=(nl+3d0*eq23)/(3d0*pi)
         aem=aem_mc/(1d0-b0*aem_mc*dlog(Q*Q/mc/mc))
      elseif(Q.ge.mc.and.Q.le.mtau) then
         nl=2
         b0=(nl+3d0*eq24)/(3d0*pi)
         aem=aem_mc/(1d0-b0*aem_mc*dlog(Q*Q/mc/mc))
      elseif(Q.ge.mtau.and.Q.le.mb) then
         nl=3
         b0=(nl+3d0*eq24)/(3d0*pi)
         aem=aem_mtau/(1d0-b0*aem_mtau*dlog(Q*Q/mtau/mtau))
      elseif(Q.ge.mb) then
         nl=3
         b0=(nl+3d0*eq25)/(3d0*pi)
         aem=aem_mz/(1d0-b0*aem_mb*dlog(Q*Q/mz/mz))
      endif
      alpha_qed=aem
      return
      end

c--------------------------------------------------------------
c     Perbative evolution kernel with bstar: Revo
c     iq = 1 / 0 for quark or gluon
c--------------------------------------------------------------
      subroutine evolve(iq,b,Q,Q0,Revo)
Cf2py intent(in)  iq
Cf2py intent(in)  b
Cf2py intent(in)  Q
Cf2py intent(in)  Q0
Cf2py intent(out) Revo
      implicit none
      real*8 Revo,b,Q0,Q,beta0,beta1,bmax,bstar,bdiv
      real*8 c0,pi,CF,CA,mb,alam,alam4,alam5,funcsp
      real*8 b0,b1b0,rnf,xq,xup,xlow,xuplog,xlowlog,sp1,sp2,sp
      integer nbel,nbuu,nz,nf,iq
      
      bdiv = 10d0
      nbel = 10 
      nbuu = 5
      nz   = 5

      c0 = 1.122919d0
      pi = 4d0*atan(1d0)
      CF = 4d0/3d0
      CA = 3d0
      mb = 4.5d0
      alam4 = 0.326d0
      alam5 = 0.226d0
      bmax = 0.5d0
      
      bstar = b/dsqrt(1d0+(b/bmax)**2d0)
      Q0 = c0/bstar
      if(Q.le.mb.and.Q0.le.mb) then
         nf = 4
         rnf = 4d0
         beta0 = 11d0-2d0*rnf/3d0
         beta1 = 102d0-38d0*rnf/3d0
         b0 = beta0
         b1b0  = beta1/(beta0**2d0)
         alam = alam4
         xq = 2d0*dlog(Q/alam)
         xup = xq
         xlow = 2d0 * dlog(Q0/alam)
         xuplog = dlog(xup)
         xlowlog = dlog(xlow)
         sp = funcsp(iq,xq,xup,xuplog,xlow,xlowlog,b0,b1b0,nf)
      elseif(Q.ge.mb.and.Q0.ge.mb) then
         nf = 5
         rnf = 5d0
         beta0 = 11d0-2d0*rnf/3d0
         beta1 = 102d0-38d0*rnf/3d0
         b0 = beta0
         b1b0  = beta1/(beta0**2d0)
         alam = alam5
         xq = 2d0*dlog(Q/alam)
         xup = xq
         xlow = 2d0 * dlog(Q0/alam)
         xuplog = dlog(xup)
         xlowlog = dlog(xlow)
         sp = funcsp(iq,xq,xup,xuplog,xlow,xlowlog,b0,b1b0,nf)
      elseif(Q.ge.mb.and.Q0.le.mb) then
         nf = 5
         rnf = 5d0
         beta0 = 11d0-2d0*rnf/3d0
         beta1 = 102d0-38d0*rnf/3d0
         b0 = beta0
         b1b0  = beta1/(beta0**2d0)
         alam = alam5
         xq = 2d0*dlog(Q/alam)
         xup = xq
         xlow =2d0*dlog(mb/alam)
         xuplog = dlog(xup)
         xlowlog = dlog(xlow)
         sp1 = funcsp(iq,xq,xup,xuplog,xlow,xlowlog,b0,b1b0,nf)

         nf = 4
         rnf = 4d0
         beta0 = 11d0-2d0*rnf/3d0
         beta1 = 102d0-38d0*rnf/3d0
         b0 = beta0
         b1b0  = beta1/(beta0**2d0)
         alam = alam4
         xq = 2d0*dlog(Q/alam)
         xup = 2d0*dlog(mb/alam)
         xlow = 2d0*dlog(Q0/alam)
         xuplog = dlog(xup)
         xlowlog = dlog(xlow)
         sp2 = funcsp(iq,xq,xup,xuplog,xlow,xlowlog,b0,b1b0,nf)
         sp = sp1 + sp2
      elseif(Q.le.mb.and.Q0.ge.mb) then
         nf = 4
         rnf = 4d0
         beta0 = 11d0-2d0*rnf/3d0
         beta1 = 102d0-38d0*rnf/3d0
         b0 = beta0
         b1b0  = beta1/(beta0**2d0)
         alam = alam4
         xq = 2d0*dlog(Q/alam)
         xup = xq
         xlow =2d0*dlog(mb/alam)
         xuplog = dlog(xup)
         xlowlog = dlog(xlow)
         sp1 = funcsp(iq,xq,xup,xuplog,xlow,xlowlog,b0,b1b0,nf)

         nf = 5
         rnf = 5d0
         beta0 = 11d0-2d0*rnf/3d0
         beta1 = 102d0-38d0*rnf/3d0
         b0 = beta0
         b1b0  = beta1/(beta0**2d0)
         alam = alam5
         xq = 2d0*dlog(Q/alam)
         xup = 2d0*dlog(mb/alam)
         xlow = 2d0*dlog(Q0/alam)
         xuplog = dlog(xup)
         xlowlog = dlog(xlow)
         sp2 = funcsp(iq,xq,xup,xuplog,xlow,xlowlog,b0,b1b0,nf)
         sp = sp1 + sp2
      endif

      Revo = dexp( -sp )
      
      return
      end

c--------------------------------------------------------------
c     Perturbative Sudakov exponent: integral of \gamma_F
c--------------------------------------------------------------

      function funcsp(iqg,xq,xup,xuplog,xlow,xlowlog,b0,b1b0,nf)
      implicit none
      real*8 funcsp,xq,xup,xuplog,xlow,xlowlog,b0,b1b0
      real*8 f00,f01,f11,f02,f12,f03,f13,f23,f04,f14,f24
      real*8 CF,CA,pi,G0,G1,G0V,sp0,sp0V,sp1,alam4,alam5,mb,c0,bmax
      integer nf,iqg

      c0 = 1.122919d0
      pi = 4d0*atan(1d0)
      CF = 4d0/3d0
      CA = 3d0
      mb = 4.5d0
      alam4 = 0.326d0
      alam5 = 0.226d0
      bmax = 0.5d0
      
      if(iqg.eq.1) then         ! for quark TMD
         G0  =  4d0 * CF
         G0V = -6d0 * CF
         G1  =  4d0*CF*( (67d0/9d0 - pi*pi/3d0)*CA - 10d0/9d0*nf )
      elseif(iqg.eq.0) then     ! for gluon TMD
         G0  =  4d0 * CA
         G0V = -2d0*(11d0/3d0*CA - 2d0/3d0*nf)
         G1  =  4d0*CA*( (67d0/9d0 - pi*pi/3d0)*CA - 10d0/9d0*nf )
      else
         print *, 'iqg should be 0 (1) for g (q): ', iqg
      endif
         
      f00 = xup - xlow
      f01 = xuplog - xlowlog
      f11 = ( xuplog**2d0 - xlowlog**2d0 ) / 2d0
      f02 = 1d0/xlow - 1d0/xup
      f12 = xlowlog/xlow - xuplog/xup + f02
      f03 = (1d0/xlow**2d0 - 1d0/xup**2d0) / 2d0
      f13 = (xlowlog/xlow**2d0 - xuplog/xup**2d0)/2d0 + f03/2d0
      f23 = ((xlowlog/xlow)**2d0 - (xuplog/xup)**2d0 )/2d0 + f13
      f04 = (1d0/xlow**3d0 - 1d0/xup**3d0)/3d0
      f14 = (xlowlog/xlow**3d0 - xuplog/xup**3d0)/3d0 + f04/3d0
      f24 = (xlowlog**2d0/xlow**3d0 - xuplog**2d0/xup**3d0)/3d0 
     >     + f14*2d0/3d0

      sp0 = G0 /(2d0*b0) * ( (xq * f01 - f00) - b1b0*(xq * f12 - f11) )
      sp0V = G0V /(2d0*b0) * ( f01 - b1b0 * f12 )
      sp1 = G1 /(2d0*b0*b0) * ( (xq*f02 - f01) - 2d0*b1b0*(xq*f13-f12)
     >     +b1b0**2d0*(xq*f24-f23) )

      funcsp = sp0 + sp0V + sp1

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
      data ktw/0.424d0/ptw/0.168d0/kappa2/0.84d0/
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
