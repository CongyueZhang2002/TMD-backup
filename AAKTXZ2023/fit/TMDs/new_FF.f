      subroutine Nuclear_TMDFF(zz,mub,AA,fu,fub,fd,fdb,fs,fsb,fg,nloops)
      IMPLICIT NONE
C--------- Inputs:  
C----- zz : The momentum fraction of pion from the parton it fragmented from
C----- mub: Natural scale for TMD evolution (defined in https://arxiv.org/pdf/2107.12401.pdf)
C----- AA: Nuclear size

C----- Outputs equation (2) in https://arxiv.org/pdf/2107.12401.pdf
C----- WITHOUT the exponential factor 

C------for each parton 
C----- i.e. "u->pion" is given by fu, etc.
      double precision Q0,Q
      double precision Q02,Q2
      double precision AlphaQCD,AlphaQED
      double precision xPDF,xgamma
      double precision eps

      REAL*8 zz, mub 
      REAL*8 u, ub, d, db, s, sb, c, b, cb, bb, g
      REAL*8 fu, fub, fd, fdb, fs, fsb, fg
      REAL*8 fun,fubn,fdn,fdbn,fsn,fsbn, fgn
      REAL*8 QQ, QQ2
      INTEGER AA,IH,IC,nloops,IO
      INTEGER INIT
      INTEGER CL, set
      INTEGER REPLICA_NUM
      INTEGER START_NUM
      double precision H_AA
      COMMON /HMASS/ H_AA 
      COMMON /REPLICA_INFO/ REPLICA_NUM
      COMMON /meson/ IH,IC

!      INTEGER FININUCLEAR
      REAL*8 fc,fb
      DOUBLE PRECISION Q2MIN, Q2MAX

      H_AA = AA 


C------  Switch notation from mub to QQ
      QQ = mub 
      QQ2 = QQ**2
      Q2MIN = 1.0D0
      Q2MAX = 1.D4
C      print*, QQ2

C----- Cutoff Q  so we dont get a headache
      IF(QQ2.LT.Q2MIN) THEN
      QQ2 = Q2MIN
      ENDIF

      IF(QQ2.GT.Q2MAX) THEN
      QQ2 = Q2MAX
      ENDIF

C------ Parameters for Evolution:
C------ Initial scale: Q02. Final scale, set to mub**2
C------ We set initial scale Q0 which the PDFs are defined at Q0 =1 GeV
C------ (Same as Zurita https://arxiv.org/pdf/2101.01088.pdf)
C----- Note: For charm quark, bottom quark Q0 must be different. 
C----- However: we do not consider charm or bottom here so we are good
C----- In setting a universal scale Q02 = 1 GeV 

      if (IC.eq.1) then
      fu  = xPDF(2,zz)/zz
      fub = xPDF(-2,zz)/zz
      fd  = xPDF(1,zz)/zz
      fdb = xPDF(-1,zz)/zz
      fs  = xPDF(3,zz)/zz
      fsb = xPDF(-3,zz)/zz
      fg  = xPDF(0,zz)/zz
      elseif (IC.eq.-1) then
      fub  = xPDF(2,zz)/zz
      fu   = xPDF(-2,zz)/zz
      fdb  = xPDF(1,zz)/zz
      fd   = xPDF(-1,zz)/zz
      fsb  = xPDF(3,zz)/zz
      fs   = xPDF(-3,zz)/zz
      fg   = xPDF(0,zz)/zz      
      elseif (IC.eq.0) then
      fu  =(xPDF(2,zz)/zz+ xPDF(-2,zz)/zz)/2d0
      fd  =(xPDF(1,zz)/zz+ xPDF(-1,zz)/zz)/2d0
      fs  =(xPDF(3,zz)/zz+ xPDF(-3,zz)/zz)/2d0
      fub =fu
      fdb =fd
      fsb =fs
      fg  = xPDF(0,zz)/zz
      else
      print *, 'charge not valid'
      endif

C----- Now perform the convolution ( only needed for NLO calculations )

      if (nloops.eq.2) then
          call nFF_NLO(zz,QQ,AA,fun,fubn,fdn,fdbn,fsn,fsbn,fgn)
          fu = fu +fun
          fub= fub+fubn
          fd = fd +fdn
          fdb= fdb+fdbn
          fs = fs +fsn
          fsb= fsb+fsbn
          fg= fg+fgn
      end if



      RETURN

      END


C---- Precalculated FF's (old --> Only use for Deuteron)
      subroutine Pre_FF_A(x,Q,A,U,D,UB,DB,SS,SB,GL)
      implicit none
      real*8 x,Q,U,D,UB,DB,SS,SB,GL,pdf(-6:6)
      real*8 Q2min,Q2max,Q2
      integer A,IC,IH
      common /meson/ IH,IC

      Q2    = Q**2d0
      Q2MIN = 1.01D0
      Q2MAX = 1.D5

      IF(Q2.LT.Q2MIN) THEN
      Q = dsqrt(Q2MIN)
      ELSEIF (Q2.GT.Q2MAX) THEN
      Q = dsqrt(Q2MAX)
      ENDIF

      u = 0d0
      d = 0d0
      ss= 0d0
      ub= 0d0
      db= 0d0
      sb= 0d0

      if (A.eq.3) then  ! Deuteron
      call evolvePDFM(23,x,Q,pdf)
      elseif (A.eq.4) then ! HE
      call evolvePDFM(6,x,Q,pdf)
      elseif (A.eq.20) then ! Ne
      call evolvePDFM(7,x,Q,pdf)
      elseif (A.eq.84) then  ! Kr
      call evolvePDFM(8,x,Q,pdf)
      elseif (A.eq.131) then  ! Xe
      call evolvePDFM(9,x,Q,pdf)
      elseif (A.eq.12) then ! C
      call evolvePDFM(16,x,Q,pdf)
      elseif (A.eq.56) then  ! Fe
      call evolvePDFM(17,x,Q,pdf)
      elseif (A.eq.197) then  ! Au
      call evolvePDFM(24,x,Q,pdf)
      elseif (A.eq.208) then  ! Pb
      call evolvePDFM(18,x,Q,pdf)
      endif

      if (IC.eq.1) then
      U=pdf(2)/x
      D=pdf(1)/x
      SS=pdf(3)/x
      UB=pdf(-2)/x
      DB=pdf(-1)/x
      SB=pdf(-3)/x
      GL=pdf(0)/x
      elseif (IC.eq.-1) then
      UB=pdf(2)/x
      DB=pdf(1)/x
      SB=pdf(3)/x
      U =pdf(-2)/x
      D =pdf(-1)/x
      SS=pdf(-3)/x
      GL=pdf(0)/x
      elseif (IC.eq.0) then
      U =(pdf( 2)/x + pdf(-2)/x)/2d0
      D =(pdf( 1)/x + pdf(-1)/x)/2d0
      SS=(pdf( 3)/x + pdf(-3)/x)/2d0
      UB=U
      DB=D
      SB=SS
      GL=pdf( 0)/x
      else
      print *, 'charge not valid'
      endif
      return
      end



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
      Ni(5)  = 0d0 !0.306d0 !/ 0.0228321d0
      ai(5)  = 1.345d0
      bi(5)  = 5.519d0
      gi(5)  = 19.78d0
      di(5)  = 10.22d0
C      Ni(5)  = Ni(5) / (EBETA(2.d0+ai(5), bi(5)+1.d0) 
C     1              + gi(5) * EBETA(2.d0+ai(5), bi(5)+di(5)+1.d0))

*     btot
      Ni(6)  =  0d0 ! 0.372d0 !/ 0.176549d0
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

C      print*, "nuclear mass: ", A 

      

      DO i=1,6
            Ni(i)  = Ni(i) / (EBETA(2.d0+ai(i), bi(i)+1.d0) 
     1              + gi(i) * EBETA(2.d0+ai(i), bi(i)+di(i)+1.d0))
            Ni(i) = Ni(i) * (1.d0 +  Nq1*(1d0-A**Nq2))
            ai(i) = ai(i) + aq1 *(1d0-A**aq2)
            bi(i) = bi(i) + bq1 *(1d0-A**bq2)
            gi(i) = gi(i) + gq1 *(1d0-A**gq2)
            di(i) = di(i) + dq1 *(1d0-A**dq2)

          
      ENDDO
            Ni(0) = Ni(0) / (EBETA(2.d0+ai(0), bi(0)+1.d0) 
     1              + gi(0) * EBETA(2.d0+ai(0), bi(0)+di(0)+1.d0))
            Ni(0) = Ni(0) * (1.d0 + Ng1*(1d0-A**Ng2))
            ai(0) = ai(0) + ag1 *(1d0-A**ag2)
            bi(0) = bi(0) + bg1 *(1d0-A**bg2)
            gi(0) = gi(0) + gg1 *(1d0-A**gg2)
            di(0) = di(0) + dg1 *(1d0-A**dg2)

          
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

c----------------------------------------------------------------
      subroutine FF(z,Q,fu,fub,fd,fdb,fs,fsb,fg,nloops)
      implicit none
      real*8 z,Q,Q2,fu,fub,fd,fdb,fs,fsb,fg
      real*8 pref,alphas,as,pi,pdf(-6:6)
      real*8 U ,UB ,D ,DB ,S ,SB ,C ,B ,GL
      real*8 U1,UB1,D1,DB1,S1,SB1,C1,B1,GL1
      real*8 U2,UB2,D2,DB2,S2,SB2,C2,B2,GL2
      integer IH,IC,IO,nloops
      common /meson/ IH,IC
      common /fforder/ IO
      data pi /3.14159265359d0/
      integer proc
      common /process/ proc
      real*8 hQ
      common /hardQ/ hQ

      fu =0.
      fub=0.
      fd =0.
      fdb=0.
      fs =0.
      fsb=0.
      fg =0.

      u =0.
      ub=0.
      d =0.
      db=0.
      s =0.
      sb=0.
      c =0.
      b =0.
      gl=0.

      if(Q.le.1d0) Q=1d0
      if(Q.gt.316d0) Q=316d0
      Q2=Q*Q
      if (IH.eq.1) then
      if (IC.eq.2) then
      call fDSSH14(0,0, 1,1,z,Q2,U1,UB1,D1,DB1,S1,SB1,C1,B1,GL1)
      call fDSSH14(0,0,-1,1,z,Q2,U2,UB2,D2,DB2,S2,SB2,C2,B2,GL2)
      U = U1 -U2
      D = D1 -D2
      S = S1 -S2
      UB= UB1-UB2
      DB= DB1-DB2
      SB= SB1-SB2
      else
      call fDSSH14(0,0,IC,1,z,Q2,U,UB,D,DB,S,SB,C,B,GL)
      endif
      elseif (IH.eq.2) then
      call fDSSH17(0,0,IC,1,z,Q2,U,UB,D,DB,S,SB,C,B,GL)
      else
      call fDSSH14(0,0,IC,1,z,Q2,U1,UB1,
     $             D1,DB1,S1,SB1,C1,B1,GL1)
      call fDSSH17(0,0,IC,1,z,Q2,U2,UB2,
     $             D2,DB2,S2,SB2,C2,B2,GL2)
      U = U1 +U2
      D = D1 +D2
      S = S1 +S2
      UB= UB1+UB2
      DB= DB1+DB2
      SB= SB1+SB2
      endif

      fu=U/z
      fub=UB/z
      fd=D/z
      fdb=DB/z
      fs=S/z
      fsb=SB/z
      fg=GL/z

      if (nloops.eq.2) then
          call FF_NLO(z,Q,u,ub,d,db,s,sb,gl)
          fu = fu +u
          fub= fub+ub
          fd = fd +d
          fdb= fdb+db
          fs = fs +s
          fsb= fsb+sb
          fg= fg+gl
      end if

      return
      end

c-----------------------------------------------------------------------
c FF at NLO
c-----------------------------------------------------------------------

      subroutine FF_NLO(z,Q,fun,fubn,fdn,fdbn,fssn,fsbn,fgn)
      implicit none
      real*8 z,Q
      real*8 fu ,fd ,fub ,fdb ,fss ,fsb ,fg
      real*8 fun,fdn,fubn,fdbn,fssn,fsbn,fgn
      real*8 cu,cd,cub,cdb,css,csb,cg

      call conv_qgauss_FF(10,z,Q,cu,cd,cub,cdb,css,csb,cg)
      fun =cu
      fdn =cd
      fubn=cub
      fdbn=cdb
      fssn=css
      fsbn=csb
      return

      end

c-----------------------------------------------------------------------
c NLO Coefficient at mu = mub and xi = mub**2
c-----------------------------------------------------------------------
      real*8 function qqff(z,Q)
      implicit none
      real*8 z,Nc,pi,euler,Q,alfas,Cf
      real*8, external :: alphas,Pqq

      pi=3.141592653589793d0

      Cf=4d0/3d0
      qqff=alphas(Q)/pi*Cf*((1d0-z)/2d0+(1d0+z*z)/(1d0-z)*dlog(z))/z
      end function qqff

      real*8 function qgff(z,Q)
      implicit none
      real*8 z,Nc,pi,euler,Q,alfas,Cf
      real*8, external :: alphas

      pi=3.141592653589793d0

      Cf=4d0/3d0
      qgff=alphas(Q)*Cf/pi*(z/2d0+(1d0+(1d0-z)**2d0)/z*dlog(z))/z

      end function qgff

c-----------------------------------------------------------------------
c     Int_x_i^x_f func(x) dx using gaussian integration (for convolving)
c-----------------------------------------------------------------------
      subroutine conv_qgauss_FF(n,x,Q,cu,cd,cub,cdb,cs,csb,cg)
      implicit none
      real*8 xi,xf,xn,x1,x2,x,Q,cu,cd,cub,cdb,cs,csb,cg
      real*8 cui,cdi,cubi,cdbi,csi,csbi,cgi
      integer i,n

      cu =0d0
      cub=0d0
      cd =0d0
      cdb=0d0
      cs =0d0
      csb=0d0
      cg=0d0

      if(n.le.1) then           ! same as n=1
         x1=x
         x2=1d0
         call conv_gauss_FF(x1,x2,x,Q,cui,cdi,cubi,cdbi,csi,csbi,cgi)
         cu =cui
         cd =cdi
         cub=cubi
         cdb=cdbi
         cs =csi
         csb=csbi
         cg=cgi
         return
      endif

      xi=x
      xf=1d0
      xn=(xf-xi)/float(n)
      x2=x
      Do 100 i=1,n
         x1=x2
         x2=x1+xn
         call conv_gauss_FF(x1,x2,x,Q,cui,cdi,cubi,cdbi,csi,csbi,cgi)
         cu =cu +cui
         cd =cd +cdi
         cub=cub+cubi
         cdb=cdb+cdbi
         cs =cs +csi
         csb=csb+csbi
         cg=cg+cgi
 100  continue
      return
      end

c-------------------------------------------------
      subroutine conv_gauss_FF(xi,xf,x,Q,cu,cd,cub,cdb,cs,csb,cg)
      implicit none
      real*8 cu,cd,cub,cdb,cs,csb,cg
      real*8 xi,xf,xm,xr,d,xd(8),w(8),eps,x,Q
      real*8 Up,Dp,UBp,DBp,SSp,SBp,GLp
      real*8 Um,Dm,UBm,DBm,SSm,SBm,GLm,gluon
      real*8 qq,gq,cp,cm
      real*8, external :: qqff,qgff
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
      DATA XD
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

      cu =0d0
      cd =0d0
      cub=0d0
      cdb=0d0
      cs =0d0
      csb=0d0

      Do 100 j=1,8
         d=xr*xd(j)
         call FF(x/(xm+d),Q,Up,UBp,Dp,DBp,SSp,SBp,GLp,1)
         call FF(x/(xm-d),Q,Um,UBm,Dm,DBm,SSm,SBm,GLm,1)
         gluon=qgff(xm+d,Q)*GLp+qgff(xm-d,Q)*GLm
         cp = qqff(xm+d,Q)
         cm = qqff(xm-d,Q)
         cu =cu +w(j)*(cp*Up +cm*Um +gluon)
         cd =cd +w(j)*(cp*Dp +cm*Dm +gluon)
         cub=cub+w(j)*(cp*UBp+cm*UBm+gluon)
         cdb=cdb+w(j)*(cp*DBp+cm*DBm+gluon)
         cs =cs +w(j)*(cp*SSp+cm*SSm+gluon)
         csb=csb+w(j)*(cp*SBp+cm*SBm+gluon)
 100  continue

      cu =xr*cu
      cub=xr*cub
      cd =xr*cd
      cdb=xr*cdb
      cs =xr*cs
      csb=xr*csb
      return
      end




c-----------------------------------------------------------------------
c FF at NLO
c-----------------------------------------------------------------------

      subroutine nFF_NLO(z,Q,A,fun,fubn,fdn,fdbn,fssn,fsbn,fgn)
      implicit none
      real*8 z,Q
      integer A
      real*8 fu ,fd ,fub ,fdb ,fss ,fsb ,fg
      real*8 fun,fdn,fubn,fdbn,fssn,fsbn,fgn
      real*8 cu,cd,cub,cdb,css,csb,cg

      call conv_qgauss_nFF(10,z,Q,A,cu,cd,cub,cdb,css,csb,cg)
      fun =cu
      fdn =cd
      fubn=cub
      fdbn=cdb
      fssn=css
      fsbn=csb
      return

      end

c-----------------------------------------------------------------------
c     Int_x_i^x_f func(x) dx using gaussian integration (for convolving)
c-----------------------------------------------------------------------
      subroutine conv_qgauss_nFF(n,x,Q,A,cu,cd,cub,cdb,cs,csb,cg)
      implicit none
      real*8 xi,xf,xn,x1,x2,x,Q,cu,cd,cub,cdb,cs,csb,cg
      real*8 cui,cdi,cubi,cdbi,csi,csbi,cgi
      integer i,n,A

      cu =0d0
      cub=0d0
      cd =0d0
      cdb=0d0
      cs =0d0
      csb=0d0
      cg=0d0

      if(n.le.1) then           ! same as n=1
         x1=x
         x2=1d0
         call conv_gauss_nFF(x1,x2,x,Q,A,cui,cdi,cubi,cdbi,csi,csbi,cgi)
         cu =cui
         cd =cdi
         cub=cubi
         cdb=cdbi
         cs =csi
         csb=csbi
         cg=cgi
         return
      endif

      xi=x
      xf=1d0
      xn=(xf-xi)/float(n)
      x2=x
      Do 100 i=1,n
         x1=x2
         x2=x1+xn
         call conv_gauss_nFF(x1,x2,x,Q,A,cui,cdi,cubi,cdbi,csi,csbi,cgi)
         cu =cu +cui
         cd =cd +cdi
         cub=cub+cubi
         cdb=cdb+cdbi
         cs =cs +csi
         csb=csb+csbi
         cg=cg+cgi
 100  continue
      return
      end

c-------------------------------------------------
      subroutine conv_gauss_nFF(xi,xf,x,Q,A,cu,cd,cub,cdb,cs,csb,cg)
      implicit none
      real*8 cu,cd,cub,cdb,cs,csb,cg
      real*8 xi,xf,xm,xr,d,xd(8),w(8),eps,x,Q
      real*8 Up,Dp,UBp,DBp,SSp,SBp,GLp
      real*8 Um,Dm,UBm,DBm,SSm,SBm,GLm,gluon
      real*8 qq,gq,cp,cm
      real*8, external :: qqff,qgff
      integer j,A
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
      DATA XD
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

      cu =0d0
      cd =0d0
      cub=0d0
      cdb=0d0
      cs =0d0
      csb=0d0

      Do 100 j=1,8
         d=xr*xd(j)
         call Nuclear_TMDFF(x/(xm+d),Q,A,Up,UBp,Dp,DBp,SSp,SBp,GLp,1)
         call Nuclear_TMDFF(x/(xm-d),Q,A,Um,UBm,Dm,DBm,SSm,SBm,GLm,1)
         gluon=qgff(xm+d,Q)*GLp+qgff(xm-d,Q)*GLm
         cp = qqff(xm+d,Q)
         cm = qqff(xm-d,Q)
         cu =cu +w(j)*(cp*Up +cm*Um +gluon)
         cd =cd +w(j)*(cp*Dp +cm*Dm +gluon)
         cub=cub+w(j)*(cp*UBp+cm*UBm+gluon)
         cdb=cdb+w(j)*(cp*DBp+cm*DBm+gluon)
         cs =cs +w(j)*(cp*SSp+cm*SSm+gluon)
         csb=csb+w(j)*(cp*SBp+cm*SBm+gluon)
 100  continue

      cu =xr*cu
      cub=xr*cub
      cd =xr*cd
      cdb=xr*cdb
      cs =xr*cs
      csb=xr*csb
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
