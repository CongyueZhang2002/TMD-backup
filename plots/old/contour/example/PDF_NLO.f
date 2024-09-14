c----------------------------------------------------------------
c     BOUND PROTON PDF in nucleus ( EPPS16 )
c----------------------------------------------------------------
      subroutine EPPS16_proton(x,Q,AA,ZZ,U,D,UB,DB,SS,SB,GL,nloops)
      implicit none
      real*8 x,Q
      real*8 U,D,UB,DB,SS,SB,GL
      real*8 Un,Dn,UBn,DBn,SSn,SBn,GLn
      real*8 fU,fD,fUB,fDB,fSS,fSB,fGL
      real*8 ruv,rdv,ru,rd,rs,rc,rb,rg
      real*8 Qmin
      real*8 UB_p,DB_p,U_p,D_p
      integer order, pset, AA, ZZ, nloops

      pset = 1


      Qmin = 1.3d0

      IF(Q.lt.Qmin) THEN
      Q = Qmin
      ENDIF

      CALL EPPS16(order,pset,AA,x,Q,ruv,rdv,ru,rd,rs,rc,rb,rg)
      CALL PDF_p(x,Q,fU,fD,fUB,fDB,fSS,fSB,fGL,1)
      UB_p = ru*fUB
      DB_p = rd*fDB
      U_p = ruv*(fU-fUB)+UB_P
      D_p = rdv*(fD-fDB)+DB_p

      U  = U_p
      D  = D_p
      UB = UB_p
      DB = DB_p
      SS=rs*fSS
      SB=rs*fSB
      GL=rg*fGL

      if (nloops.eq.2) then
          call PDF_EPPS16_NLO(x,Q,AA,ZZ,Un,Dn,UBn,DBn,SSn,SBn,GLn)
          U =U +Un
          D =D +Dn
          SS=SS+SSn
          SB=SB+SBn
          UB=UB+UBn
          DB=DB+DBn
      end if

      return
      end
c-----------------------------------------------------------------------
c PDF (in a nucleus) at NLO using EPPS16
c-----------------------------------------------------------------------

      subroutine PDF_EPPS16_NLO(x,Q,AA,ZZ,UN,DN,UBN,DBN,SSN,SBN,GL)
      implicit none
      real*8 x,Q,cxf
      real*8 U,D,UB,DB,SS,SB,GL,UN,DN,UBN,DBN,SSN,SBN,GLN
      real*8 cu,cd,cub,cdb,css,csb,cg
      integer AA,ZZ


      call conv_qgauss_PDF_EPPS(10,x,Q,AA,ZZ,cu,cd,cub,cdb,css,csb,cg)

      UN =cu
      DN =cd
      UBN=cub
      DBN=cdb
      SSN=css
      SBN=csb
      return

      end

c-----------------------------------------------------------------------
c PDF (proton) at NLO
c-----------------------------------------------------------------------

      subroutine PDF_P_NLO(x,Q,UN,DN,UBN,DBN,SSN,SBN,GL)
      implicit none
      real*8 x,Q,cxf
      real*8 U,D,UB,DB,SS,SB,GL,UN,DN,UBN,DBN,SSN,SBN,GLN
      real*8 cu,cd,cub,cdb,css,csb,cg

      call conv_qgauss_PDF(10,x,Q,cu,cd,cub,cdb,css,csb,cg)

      UN =cu
      DN =cd
      UBN=cub
      DBN=cdb
      SSN=css
      SBN=csb
      return

      end

c-----------------------------------------------------------------------
c NLO Coefficient at mu = mub and xi = mub**2
c-----------------------------------------------------------------------
      real*8 function qq(x,Q)
      implicit none
      real*8 x,Nc,pi,euler,Q,alfas,Cf
      real*8, external :: alphas

      pi=3.141592653589793d0

      Cf=1.3333333d0
      qq=alphas(Q)*Cf/2d0/pi*(1d0-x)/x

      end function

      real*8 function gq(x,Q)
      implicit none
      real*8 x,Nc,pi,euler,Q,alfas,TR
      real*8, external :: alphas

      pi=3.141592653589793d0

      TR=0.5d0
      gq=alphas(Q)*TR/pi*x*(1d0-x)/x

      end function

c-----------------------------------------------------------------------
c     Int_x_i^x_f func(x) dx using gaussian integration (for convolving)
c-----------------------------------------------------------------------
      subroutine conv_qgauss_PDF(n,x,Q,cu,cd,cub,cdb,cs,csb,cg)
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
         call conv_gauss_PDF(x1,x2,x,Q,cui,cdi,cubi,cdbi,csi,csbi,cgi)
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
         call conv_gauss_PDF(x1,x2,x,Q,cui,cdi,cubi,cdbi,csi,csbi,cgi)
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

c-----------------------------------------------------------------------
c     Int_x_i^x_f func(x) dx using gaussian integration (for convolving)
c-----------------------------------------------------------------------
      subroutine conv_qgauss_PDF_EPPS(n,x,Q,A,Z,cu,cd,cub,
     >                               cdb,cs,csb,cg)
      implicit none
      real*8 xi,xf,xn,x1,x2,x,Q,cu,cd,cub,cdb,cs,csb,cg
      real*8 cui,cdi,cubi,cdbi,csi,csbi,cgi
      integer i,n
      integer A,Z

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
         call conv_gauss_PDF_EPPS(x1,x2,x,Q,A,Z,cui,
     >        cdi,cubi,cdbi,csi,csbi,cgi)
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
         call conv_gauss_PDF_EPPS(x1,x2,x,Q,A,Z,cui,
     >        cdi,cubi,cdbi,csi,csbi,cgi)
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
C-------------------------------------------------
      subroutine conv_gauss_PDF_EPPS(xi,xf,x,Q,A,Z,
     >           cu,cd,cub,cdb,cs,csb,cg)
      implicit none
      real*8 cu,cd,cub,cdb,cs,csb,cg
      real*8 xi,xf,xm,xr,d,xd(8),w(8),eps,x,Q
      real*8 Up,Dp,UBp,DBp,SSp,SBp,GLp
      real*8 Um,Dm,UBm,DBm,SSm,SBm,GLm,gluon
      real*8 qq,gq
      integer A,Z
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
         call EPPS16_proton(x/(xm+d),Q,A,Z,
     >                  Up,Dp,UBp,DBp,SSp,SBp,GLp,1)
         call EPPS16_proton(x/(xm-d),Q,A,Z,
     >                  Um,Dm,UBm,DBm,SSm,SBm,GLm,1)
         gluon=gq(xm+d,Q)*GLp+gq(xm-d,Q)*GLm
         cu =cu +w(j)*(qq(xm+d,Q)*Up +qq(xm-d,Q)*Um +gluon)
         cd =cd +w(j)*(qq(xm+d,Q)*Dp +qq(xm-d,Q)*Dm +gluon)
         cub=cub+w(j)*(qq(xm+d,Q)*UBp+qq(xm-d,Q)*UBm+gluon)
         cdb=cdb+w(j)*(qq(xm+d,Q)*DBp+qq(xm-d,Q)*DBm+gluon)
         cs =cs +w(j)*(qq(xm+d,Q)*SSp+qq(xm-d,Q)*SSm+gluon)
         csb=csb+w(j)*(qq(xm+d,Q)*SBp+qq(xm-d,Q)*SBm+gluon)
 100  continue

      cu =xr*cu
      cub=xr*cub
      cd =xr*cd
      cdb=xr*cdb
      cs =xr*cs
      csb=xr*csb
      return
      end

c-------------------------------------------------
      subroutine conv_gauss_PDF(xi,xf,x,Q,cu,cd,cub,cdb,cs,csb,cg)
      implicit none
      real*8 cu,cd,cub,cdb,cs,csb,cg
      real*8 xi,xf,xm,xr,d,xd(8),w(8),eps,x,Q
      real*8 Up,Dp,UBp,DBp,SSp,SBp,GLp
      real*8 Um,Dm,UBm,DBm,SSm,SBm,GLm,gluon
      real*8 qq,gq
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
         call PDF_P(x/(xm+d),Q,Up,Dp,UBp,DBp,SSp,SBp,GLp,1)
         call PDF_P(x/(xm-d),Q,Um,Dm,UBm,DBm,SSm,SBm,GLm,1)
         gluon=gq(xm+d,Q)*GLp+gq(xm-d,Q)*GLm
         cu =cu +w(j)*(qq(xm+d,Q)*Up +qq(xm-d,Q)*Um +gluon)
         cd =cd +w(j)*(qq(xm+d,Q)*Dp +qq(xm-d,Q)*Dm +gluon)
         cub=cub+w(j)*(qq(xm+d,Q)*UBp+qq(xm-d,Q)*UBm+gluon)
         cdb=cdb+w(j)*(qq(xm+d,Q)*DBp+qq(xm-d,Q)*DBm+gluon)
         cs =cs +w(j)*(qq(xm+d,Q)*SSp+qq(xm-d,Q)*SSm+gluon)
         csb=csb+w(j)*(qq(xm+d,Q)*SBp+qq(xm-d,Q)*SBm+gluon)
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
c PDF at NLO
c-----------------------------------------------------------------------

      subroutine PDF_N_NLO(x,Q,UN,DN,UBN,DBN,SSN,SBN,GL)
      implicit none
      real*8 UN,DN,UBN,DBN,SSN,SBN,GL,Q,x

      call PDF_P_NLO(x,Q,DN,UN,DBN,UBN,SSN,SBN,GL)

      end


c-----------------------------------------------------------------------
c PDF at NLO
c-----------------------------------------------------------------------

      subroutine PDF_D_NLO(x,Q,UN,DN,UBN,DBN,SSN,SBN,GL)
      implicit none
      real*8 U,D,UB,DB,SS,SB
      real*8 UN,DN,UBN,DBN,SSN,SBN,GL,Q,x

      call PDF_P_NLO(x,Q,U,D,UB,DB,SS,SB,GL)
      UN = (U+D)/2.0
      DN = (U+D)/2.0
      UBN = (UB+DB)/2.0
      DBN = (UB+DB)/2.0
      SSN = SS
      SBN = SB

      end
