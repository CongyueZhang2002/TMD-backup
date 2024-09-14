c-----f(x,kt,Q). x: bjorken-x, Q: scale, kt: transverse momentum
C-----NUMBERING SYSTEM FOR PARTON TYPE FOR PDF_KT
C-----1: u   2:  d    3: s   4: ub   5: db   6: sb
      subroutine PDF_kt(kt,x,Q,iq,IT,part)
      implicit none
      real*8 kt,x,Q
      integer iq
      integer IT, AA
      real*8 part
      real*8 xx,QQ
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

! modification

      subroutine PDF_b(b,x,Q,iq,AA,part)
      implicit none
      real*8 b,x,Q
      integer iq
      integer AA,ZZ
      real*8 part
      real*8 parts
      real*8 U, UB ,D ,DB ,S ,SB, GL !2023
      REAL*8 aN, bN, Ng1, Nq1, Ng2, Nq2, bg1, gg1, dg1, gamma, g3 !2023
      COMMON /FITP/ 
     &       aN, bN, Ng1, Nq1, Ng2, Nq2, bg1, gg1, dg1, gamma, g3 !2023
      real*8 c0,c02,pi,CF,CA,mb,alam4,alam5,
     >     bmax,ktw,ptw,Qini,kappa2
      common /cons/ c0,c02,pi,CF,CA,mb,alam4,alam5,
     >     bmax,ktw,ptw,Qini,kappa2
      integer nloops,hop,nll
      common /scheme/ nloops,hop,nll
      real*8 bstar,Q0,Revo,FNP,kappa1,pref,widpdf, widuu

      bstar = b/dsqrt(1d0+(b/bmax)**2d0)
      Q0 = c0/bstar
      call kernel_q(bstar,Q0,Q,Q0,Q,nll,Revo)

      widpdf = ktw+aN*(AA**(1d0/3d0)-1d0) !2023
      widuu = kappa2*dlog(Q/Qini)
      kappa1 = widpdf/4
      FNP = dexp( -(kappa1)*b*b  - dlog(b/bstar)*widuu/2d0)

      pref = b*Revo*FNP

      if (AA.eq.1) then !2023
        call PDF_P(x,Q,U,D,UB,DB,S,SB,GL,nloops)
      else
        ZZ = 82
        call EPPS16_proton(x,Q,AA,U,D,UB,DB,S,SB,GL,nloops)
      endif

      if (iq.eq.1) then
      parts = U
      elseif (iq.eq.2) then
      parts = D
      elseif (iq.eq.3) then
      parts = S
      elseif (iq.eq.-1) then
      parts = UB
      elseif (iq.eq.-2) then
      parts = DB
      elseif (iq.eq.-3) then
      parts = SB
      endif

      part  = pref*parts

      return
      end

c----------------------------------------------------------------
c     BOUND PROTON PDF in nucleus ( EPPS16 )
c----------------------------------------------------------------
      subroutine EPPS16_proton(x,Q,AA,U,D,UB,DB,SS,SB,GL,nloops)
      implicit none
      real*8 x,Q
      real*8 U,D,UB,DB,SS,SB,GL
      real*8 Un,Dn,UBn,DBn,SSn,SBn,GLn
      real*8 fU,fD,fUB,fDB,fSS,fSB,fGL
      real*8 ruv,rdv,ru,rd,rs,rc,rb,rg
      real*8 Qmin
      real*8 UB_p,DB_p,U_p,D_p
      integer order, pset, AA, ZZ, nloops
      integer REPLICA_NUM

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