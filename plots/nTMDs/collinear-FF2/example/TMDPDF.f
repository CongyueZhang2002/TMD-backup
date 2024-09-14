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

      ZZ = 1 ! this doesnt actually do anything.

      if(AA.eq.1) then
      call Pre_PDF_A(x,Q,AA,U,D,UB,DB,SS,SB,GL)
      else
      call PDF_A_EPPS16(x,Q,AA,ZZ,U,D,UB,DB,SS,SB,GL,nloops)
      endif

      if(AA.eq.1)then
      u = (u+d)/2
      d = u
      ub = (ub+db)/2
      db = ub
      endif

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
      call evolvePDFM(1,x,Q,pdf)
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
