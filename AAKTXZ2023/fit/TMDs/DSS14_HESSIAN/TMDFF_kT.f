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

      AA = IT
      xx = x
      QQ = Q
      iqq= iq

      write(1010,*) qFFb(1d0)
      call adogt(qFFb,kt/x,Q,0,1d0,199d0,part)

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
      real*8 bstar,Q0,Revo,FNP,kappa1,pref, widpdf, widuu, widff, cor !2023

      bstar = b/dsqrt(1d0+(b/bmax)**2d0)
      Q0 = c0/bstar
      call kernel_q(bstar,Q0,Q,Q0,Q,nll,Revo)

      widff  = ptw+bN*(AA**(1d0/3d0)-1d0) !2023
      widuu = kappa2*dlog(Q/Qini)
      kappa1 = (widff/4/(x*x))
      cor = g3*(AA**(1d0/3d0)-1d0)*(Qini/Q)**gamma !2023
      FNP = dexp( -(kappa1)*b*b  - dlog(b/bstar)*widuu/2 - cor) !2023

      pref = b*Revo*FNP/x/x

      if (AA.eq.1) then
      call FF(x,Q,fu,fub,fd,fdb,fs,fsb,fg,nloops)
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

      return
      end

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