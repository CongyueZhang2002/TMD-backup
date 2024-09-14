C----- nFF (nuclear fDSS)
      Subroutine LIKEnFF(zz,QQ,AA,fu,fub,fd,fdb,fs,fsb,fg,nloops)
      IMPLICIT NONE
      REAL*8 zz,QQ
      REAL*8 u, ub, d, db, s, sb, c, b, cb, bb, g
      REAL*8 fu, fub, fd, fdb, fs, fsb, fg
      REAL*8 fun,fubn,fdn,fdbn,fsn,fsbn, fgn
      REAL*8 QQ2
      INTEGER AA,IH,IC,nloops,IO
      INTEGER INIT
      integer REPLICA_NUM
      INTEGER CL, set
      COMMON /rep_num/ REPLICA_NUM
      common /meson/ IH,IC
!      INTEGER FININUCLEAR
      REAL*8 fc,fb
      DOUBLE PRECISION Q2MIN, Q2MAX
      COMMON / INITIATE / INIT

!      INIT = 0
      IH =1
      IC =1

      IO = nloops
      QQ2 = QQ**2
      Q2MIN = 1.01D0
      Q2MAX = 1.D5

C-----Confidence Levl
      CL = 1
      set = REPLICA_NUM

      IF(QQ2.LT.Q2MIN) THEN
      QQ2 = Q2MIN
      ENDIF

      IF(QQ2.GT.Q2MAX) THEN
      QQ2 = Q2MAX
      ENDIF

      CALL LIKEn(CL,IH,IC,set,ZZ,QQ2,AA,U,UB,D,DB,S,SB,C,CB,B,BB,G)

      fu=u/zz
      fub =UB/zz
      fd =D/zz
      fdb =DB/zz
      fs =S/zz
      fsb =SB/zz
      fg =g/zz

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
         call LIKEnFF(x/(xm+d),Q,A,Up,UBp,Dp,DBp,SSp,SBp,GLp,1)
         call LIKEnFF(x/(xm-d),Q,A,Um,UBm,Dm,DBm,SSm,SBm,GLm,1)
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
