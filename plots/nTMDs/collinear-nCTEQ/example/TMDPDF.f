c-----f(x,kt,Q). x: bjorken-x, Q: scale, kt: transverse momentum
C-----NUMBERING SYSTEM FOR PARTON TYPE FOR PDF_KT
C-----1: u   2:  d    3: s   4: ub   5: db   6: sb
      subroutine PDF_kt(kt,x,Q,iq,IT,part)
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
      real*8 qPDFb

      AA = IT
      xx = x
      QQ = Q
      iqq= iq

      write(1010,*) qPDFb(1d0)
      call adogt(qPDFb,kt,Q,0,1d0,199d0,part)

      return
      end

      function qPDFb(b)
      implicit none
      real*8 qPDFb
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

      call PDF_b(b,x,Q,iq,AA,qPDFb)

      return
      end

      subroutine PDF_b(b,x,Q,iq,AA,part)
      implicit none
      real*8 b,x,Q
      integer iq
      integer AA
      real*8 part
      real*8 parts
      real*8 u ,ub ,d ,db ,s ,sb, gl
      real*8 us,ubs,ds,dbs,ss,sbs
      real*8 Ans ,Bns, GDs, ADs, GPs, APs
      common /fitp/
     >       Ans ,Bns ,GDs, Ads, GPs, APs
      real*8 c0,c02,pi,CF,CA,mb,alam4,alam5,
     >     bmax,ktw,ptw,Qini,kappa2
      common /cons/ c0,c02,pi,CF,CA,mb,alam4,alam5,
     >     bmax,ktw,ptw,Qini,kappa2
      integer nloops,hop,nll
      integer ZZ
      common /scheme/ nloops,hop,nll
      real*8 bstar,Q0,Revo,FNP,kappa1,pref, widpdf, widuu

      bstar = b/dsqrt(1d0+(b/bmax)**2d0)
      Q0 = c0/bstar
      call kernel_q(bstar,Q0,Q,Q0,Q,nll,Revo)

      ! Quick Note: ktw, kappa2, ptw are in common.f!

      widpdf = ktw + Ans*(AA**(1d0/3d0)-1d0)
      widuu = kappa2*dlog(Q/Qini)
      kappa1 = widpdf/4
      FNP = dexp( -(kappa1)*b*b  - dlog(b/bstar)*widuu/2d0)

      pref = b*Revo*FNP

      if(AA.eq.1) then
      call Pre_PDF_A(x,Q,AA,U,D,UB,DB,SS,SB,GL)
      else
      call PDF_A_nCTEQ15(x,Q,AA,1,U,D,UB,DB,SS,SB,GL,nloops)
      endif

C      if(AA.eq.1)then
C      u = (u+d)/2
C      d = u
C      ub = (ub+db)/2
C      db = ub
C      endif

      if (iq.eq.1) then
      parts = u
      elseif (iq.eq.2) then
      parts = d
      elseif (iq.eq.3) then
      parts = s
      elseif (iq.eq.-1) then
      parts = ub
      elseif (iq.eq.-2) then
      parts = db
      elseif (iq.eq.-3) then
      parts = sb
      endif

      part  = pref*parts

      return
      end

      subroutine Pre_PDF_A(x,Q,AA,U,D,UB,DB,SS,SB,GL)
      implicit none
      real*8 x,Q,U,D,UB,DB,SS,SB,GL,pdf(-6:6)
      real*8 Qmin,Qmax
      integer AA

      Qmin = 1.3d0
      Qmax = dsqrt(100000000d0)

      IF(Q.lt.Qmin) THEN
      Q = Qmin
      ELSEIF (Q.gt.Qmax) then
      Q = Qmax
      ENDIF

      u = 0d0
      d = 0d0
      ss= 0d0
      ub= 0d0
      db= 0d0
      sb= 0d0

      if (AA.eq.4) then
      call evolvePDFM(2,x,Q,pdf)
      elseif (AA.eq.20) then
      call evolvePDFM(3,x,Q,pdf)
      elseif (AA.eq.84) then
      call evolvePDFM(4,x,Q,pdf)
      elseif (AA.eq.131) then
      call evolvePDFM(5,x,Q,pdf)
      elseif (AA.eq.9) then
      call evolvePDFM(10,x,Q,pdf)
      elseif (AA.eq.56) then
      call evolvePDFM(11,x,Q,pdf)
      elseif (AA.eq.184) then
      call evolvePDFM(12,x,Q,pdf)
      elseif (AA.eq.1) then
      call evolvePDFM(21,x,Q,pdf)
      endif

      U=pdf(2)/x
      D=pdf(1)/x
      SS=pdf(3)/x
      UB=pdf(-2)/x
      DB=pdf(-1)/x
      SB=pdf(-3)/x
      GL=pdf(0)/x

      return
      end

      subroutine PDF_A_nCTEQ15(x,Q,A,Z,U,D,UB,DB,SS,SB,GL,nloops)
      implicit none
      real*8 x,Q,U,D,UB,DB,SS,SB,GL,pdf(-6:6)
      real*8 Qmin,Qmax
      real*8 Un,Dn,UBn,DBn,SSn,SBn,GLn
      real*8 as,alphas,pref,pi
      data pi /3.14159265359d0/
      integer A,Z
      integer nloops
      integer proc
      common /process/ proc
      real*8 hQ
      common /hardQ/ hQ

      Qmin = 1.3d0
      Qmax = 1d4

      IF(Q.lt.Qmin) THEN
      Q = Qmin
      ELSEIF (Q.gt.Qmax) then
      Q = Qmax
      ENDIF

      u = 0d0
      d = 0d0
      ss= 0d0
      ub= 0d0
      db= 0d0
      sb= 0d0

      if (A.eq.1) then
      call evolvePDFM(21,x,Q,pdf) ! Proton
      elseif (A.eq.4) then ! He
      call evolvePDFM(2,x,Q,pdf)
      elseif (A.eq.20) then ! Ne
      call evolvePDFM(3,x,Q,pdf)
      elseif (A.eq.84) then ! Kr
      call evolvePDFM(4,x,Q,pdf)
      elseif (A.eq.131) then ! Xe
      call evolvePDFM(5,x,Q,pdf)
      elseif (A.eq.9) then ! Be
      call evolvePDFM(10,x,Q,pdf)
      elseif (A.eq.56) then ! Fe (E866)
      call evolvePDFM(11,x,Q,pdf)
      elseif (A.eq.184) then ! W
      call evolvePDFM(12,x,Q,pdf)
      elseif (A.eq.12) then ! C
      call evolvePDFM(13,x,Q,pdf)
      elseif (A.eq.57) then ! Fe (JLAB)
      call evolvePDFM(14,x,Q,pdf)
      elseif (A.eq.208) then ! Pb
      call evolvePDFM(15,x,Q,pdf)
      elseif (A.eq.209) then ! Pb (LHC)
      call evolvePDFM(19,x,Q,pdf)
      elseif (A.eq.197) then ! Au (RHIC)
      call evolvePDFM(20,x,Q,pdf)
      elseif (A.eq.40) then ! Ca
      call evolvePDFM(22,x,Q,pdf)
      endif

      U=pdf(2)/x
      D=pdf(1)/x
      SS=pdf(3)/x
      UB=pdf(-2)/x
      DB=pdf(-1)/x
      SB=pdf(-3)/x
      GL=pdf(0)/x


      if (nloops.eq.2) then
          call PDF_nCTEQ15_NLO(x,Q,A,Z,Un,Dn,UBn,DBn,SSn,SBn,GLn)
          U =U +Un
          D =D +Dn
          SS=SS+SSn
          SB=SB+SBn
          UB=UB+UBn
          DB=DB+DBn
      end if

      return
      end
C
c-----------------------------------------------------------------------
c PDF (in a nucleus) at NLO using EPPS16
c-----------------------------------------------------------------------

      subroutine PDF_nCTEQ15_NLO(x,Q,AA,ZZ,UN,DN,UBN,DBN,SSN,SBN,GL)
      implicit none
      real*8 x,Q,cxf
      real*8 U,D,UB,DB,SS,SB,GL,UN,DN,UBN,DBN,SSN,SBN,GLN
      real*8 cu,cd,cub,cdb,css,csb,cg
      integer AA,ZZ


      call conv_qgauss_PDF_nCTEQ15(10,x,Q,AA,ZZ,cu,cd,cub,cdb,css,csb,
     >    cg)

      UN =cu
      DN =cd
      UBN=cub
      DBN=cdb
      SSN=css
      SBN=csb
      return

      end


C
c-----------------------------------------------------------------------
c     Int_x_i^x_f func(x) dx using gaussian integration (for convolving)
c-----------------------------------------------------------------------
      subroutine conv_qgauss_PDF_nCTEQ15(n,x,Q,A,Z,cu,cd,cub,
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
         call conv_gauss_PDF_nCTEQ15(x1,x2,x,Q,A,Z,cui,
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
         call conv_gauss_PDF_nCTEQ15(x1,x2,x,Q,A,Z,cui,
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
      subroutine conv_gauss_PDF_nCTEQ15(xi,xf,x,Q,A,Z,
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
         call PDF_A_nCTEQ15(x/(xm+d),Q,A,Z,
     >                  Up,Dp,UBp,DBp,SSp,SBp,GLp,1)
         call PDF_A_nCTEQ15(x/(xm-d),Q,A,Z,
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
c---------------------------------------------------------------

      DOUBLE PRECISION FUNCTION DBFINT(NARG,ARG,NA,ENT,TABLE)
      implicit real*8 (a-h,o-z)
      INTEGER NA(NARG), INDEX(32)
      double precision
     +       ARG(NARG),ENT(32),TABLE(10),WEIGHT(32)
      DATA ZEROD/0.D0/ONED/1.D0/
C

           DBFINT =  ZEROD
           IF(NARG .LT. 1  .OR.  NARG .GT. 5)  RETURN
C
           LMAX      =  0
           ISTEP     =  1
           KNOTS     =  1
           INDEX(1)  =  1
           WEIGHT(1) =  ONED
           DO 100    N  =  1, NARG
              X     =  ARG(N)
              NDIM  =  NA(N)
              LOCA  =  LMAX
              LMIN  =  LMAX + 1
              LMAX  =  LMAX + NDIM
              IF(NDIM .GT. 2)  GOTO 10
              IF(NDIM .EQ. 1)  GOTO 100
              H  =  X - ENT(LMIN)
              IF(H .EQ. ZEROD)  GOTO 90
              ISHIFT  =  ISTEP
              IF(X-ENT(LMIN+1) .EQ. ZEROD)  GOTO 21
              ISHIFT  =  0
              ETA     =  H / (ENT(LMIN+1) - ENT(LMIN))
              GOTO 30
   10         LOCB  =  LMAX + 1
   11         LOCC  =  (LOCA+LOCB) / 2
              IF(X-ENT(LOCC))  12, 20, 13
   12         LOCB  =  LOCC
              GOTO 14
   13         LOCA  =  LOCC
   14         IF(LOCB-LOCA .GT. 1)  GOTO 11
              LOCA    =  MIN ( MAX (LOCA,LMIN), LMAX-1 )
              ISHIFT  =  (LOCA - LMIN) * ISTEP
              ETA     =  (X - ENT(LOCA)) / (ENT(LOCA+1) - ENT(LOCA))
              GOTO 30
   20         ISHIFT  =  (LOCC - LMIN) * ISTEP
   21         DO 22  K  =  1, KNOTS
                 INDEX(K)  =  INDEX(K) + ISHIFT
   22            CONTINUE
              GOTO 90
   30         DO 31  K  =  1, KNOTS
                 INDEX(K)         =  INDEX(K) + ISHIFT
                 INDEX(K+KNOTS)   =  INDEX(K) + ISTEP
                 WEIGHT(K+KNOTS)  =  WEIGHT(K) * ETA
                 WEIGHT(K)        =  WEIGHT(K) - WEIGHT(K+KNOTS)
   31            CONTINUE
              KNOTS  =  2*KNOTS
   90         ISTEP  =  ISTEP * NDIM
  100         CONTINUE
           DO 200    K  =  1, KNOTS
              I  =  INDEX(K)
              DBFINT =  DBFINT + WEIGHT(K) * TABLE(I)
  200         CONTINUE
           RETURN
           END
