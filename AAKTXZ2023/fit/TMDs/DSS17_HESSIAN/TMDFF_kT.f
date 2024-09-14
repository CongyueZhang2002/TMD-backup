c-----f(z,kt,Q). the x in FF_kt is of-course z,
C-----Q: scale, kt: transverse momentum
C-----NUMBERING SYSTEM FOR PARTON TYPE FOR PDF_KT
C-----1: u   2:  d    3: s   4: ub   5: db   6: sb
      subroutine FF_kt(kt,x,Q,iq,IT,part)
      implicit none
      real*8 kt,x,Q
      integer iq
      integer IT, AA
      real*8 part
      real*8 u ,ub ,d ,db ,s ,sb
      real*8 u1,ub1,d1,db1,s1,sb1
      real*8 xx,QQ
      real*8 M
      integer iqq
      common /bint/ xx,QQ,iqq
      common /NuclearMass/ AA
      real*8 qFFb
      external qFFb

      AA = IT
      xx = x
      QQ = Q
      iqq= iq

      !write(1010,*) qFFb(1d0)
      call adogt(qFFb,kt/ x,Q,0,1d0,60d0,part)

      return
      end

      function qFFb(b)
      implicit none
      real*8 qFFb
      real*8 b,x,Q
      integer iq
      integer AA
      real*8 u,ub,d,db,s,sb
      real*8 xx,QQ
      integer iqq
      common /bint/ xx,QQ,iqq
      common /NuclearMass/ AA

      x = xx
      Q = QQ
      iq= iqq

      call FF_b(b,x,Q,iq,AA,qFFb)

      return
      end

! modification

      subroutine FF_b(b,x,Q,iq,AA,part)
      implicit none
      real*8 b,x,Q
      integer iq
      integer AA
      real*8 part
      real*8 parts
      real*8 fu,fub,fd,fdb,fs,fsb,fg !2023
      REAL*8 aN, bN, Ng1, Nq1, Ng2, Nq2, bg1, gg1, dg1, gamma, g3 !2023
      COMMON /FITP/ 
     &       aN, bN, Ng1, Nq1, Ng2, Nq2, bg1, gg1, dg1, gamma, g3 !2023
      real*8 c0,c02,pi,CF,CA,mb,alam4,alam5,
     >     bmax,ktw,ptw,Qini,kappa2
      common /cons/ c0,c02,pi,CF,CA,mb,alam4,alam5,
     >     bmax,ktw,ptw,Qini,kappa2
      integer nloops,hop,nll
      common /scheme/ nloops,hop,nll
      real*8 bstar,Q00,Revo,FNP,kappa1,pref,widpdf,widuu,widff,corff !2023

      bstar = b/dsqrt(1d0+(b/bmax)**2d0)
      Q00 = c0/bstar
      call kernel_q(bstar,Q00,Q,Q00,Q,nll,Revo)

      widff  = ptw+bN*(AA**gamma-1d0) !2023
      widuu = kappa2*dlog(Q/Qini)
      kappa1 = (widff/4/(x*x))
      corff = 0.5*g3*(AA**gamma-1d0)*(Qini/Q)**2 !2023
      FNP = dexp( -(kappa1)*b*b - dlog(b/bstar)*widuu/2 - corff) !2023

      pref = b*Revo*FNP/x/x

      if (AA.eq.1) then
      nloops = 2
            call FF_vacuum(x,Q,fu,fub,fd,fdb,fs,fsb,fg,nloops)
      else
      
            call Nuclear_TMDFF(x,Q,AA,fu,fub,fd,fdb,fs,fsb,fg,nloops) !2023
      endif

      if (iq.eq.1) then 
      parts = fu !2023
      elseif (iq.eq.2) then
      parts = fd !2023
      elseif (iq.eq.3) then
      parts = fs !2023
      elseif (iq.eq.-1) then
      parts = fub !2023
      elseif (iq.eq.-2) then
      parts = fdb !2023
      elseif (iq.eq.-3) then
      parts = fsb !2023
      endif

      part  = pref*parts

      !print *, part

      return
      end

      subroutine FF_vacuum(z,Q,fu,fub,fd,fdb,fs,fsb,fg,nloops)
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

      ! Set Pion to Pi +
      IH = 1
      IC = 1

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
      !call fDSS(IH,IC,1,z,Q2,U,UB,D,DB,S,SB,C,B,GL)
      call fDSSH17(0,0,IC,1,z,Q2,U,UB,D,DB,S,SB,C,B,GL)
      !if (IC.eq.1) then
      !call evolvePDFM(3,z,Q,pdf)
      !U = pdf( 2)
      !D = pdf( 1)
      !S = pdf( 3)
      !UB= pdf(-2)
      !DB= pdf(-1)
      !SB= pdf(-3)
      !else
      !call evolvePDFM(4,z,Q,pdf)
      !U	= pdf( 2)
      !D	= pdf( 1)
      !S	= pdf( 3)
      !UB= pdf(-2)
      !DB= pdf(-1)
      !SB= pdf(-3)
      !endif
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
          !as = alphas(hQ)
          call FF_NLO_vacuum(z,Q,u,ub,d,db,s,sb,gl)
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

      subroutine FF_NLO_vacuum(z,Q,fun,fubn,fdn,fdbn,fssn,fsbn,fgn)
      implicit none
      real*8 z,Q
      real*8 fu ,fd ,fub ,fdb ,fss ,fsb ,fg
      real*8 fun,fdn,fubn,fdbn,fssn,fsbn,fgn
      real*8 cu,cd,cub,cdb,css,csb,cg

      call conv_qgauss_FF_vacuum(10,z,Q,cu,cd,cub,cdb,css,csb,cg)
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
      subroutine conv_qgauss_FF_vacuum(n,x,Q,cu,cd,cub,cdb,cs,csb,cg)
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
         call conv_gauss_FF_vacuum(x1,x2,x,Q,cui,cdi,cubi,cdbi,csi,
     &                           csbi,cgi)
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
         call conv_gauss_FF_vacuum(x1,x2,x,Q,cui,cdi,cubi,cdbi,csi,
     &                           csbi,cgi)
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
      subroutine conv_gauss_FF_vacuum(xi,xf,x,Q,cu,cd,cub,cdb,cs,
     &                               csb,cg)
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
         call FF_vacuum(x/(xm+d),Q,Up,UBp,Dp,DBp,SSp,SBp,GLp,1)
         call FF_vacuum(x/(xm-d),Q,Um,UBm,Dm,DBm,SSm,SBm,GLm,1)
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
