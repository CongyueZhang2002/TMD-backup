c--------------------------------------------------------------
c     DIS X-section in {qt}: d[sigma]/d^2qt
c     integrated in Q2, xb, z
c--------------------------------------------------------------
      subroutine xzQintDISUU(rtss, QQ2min,QQ2max,xbbmin,xbbmax,
     >                     zzmin,zzmax, qtt,fuu,
     >                     ITT,IHH,ICC)
      implicit none
      real*8 rtss,QQ2,zhh,xbbmin,xbbmax, QQ2min, QQ2max
      real*8 zzmin, zzmax, zmin, zmax, Q2min, Q2max
      real*8 xbmin, xbmax
      real*8 qtt,fuu,bdiv,BdepDIS,qgauss
      real*8 rts,Q2,xb,zh,qt,fut,AN,BdepDIS_Siv
      real*8 BdepDIS2,BdepDIS_Siv2,pt
      real*8 fut_res,fuu_res
      integer ngauss,nbel,nbuu,nbut,nQ,nY,ITT,IT,IH,IC,IHH,ICC
      common /ngbesl/bdiv,nbel,nbuu,nbut,nQ,nY
      common /meson/ IH,IC

      common /xzQ2intsid/ rts, qt, xbmin, xbmax, zmin, zmax, IT
      integer fft
      common /FT/ fft
      external BdepDIS_Siv,BdepDIS,BdepDIS_Siv2,BdepDIS2,qgauss
      external BdepDIS_og, BdepDIS_Siv_og,TMD_FF2
      real*8 Q,Qi,b,bstar,Revo,Revo1,sig
      real*8 fun_DIS_xzQ
      external fun_DIS_xzQ
      integer i

      rts= rtss
      qt = qtt
      IT = ITT
      IH = IHH
      IC = ICC
      Q2min = QQ2min
      Q2max = QQ2max
      xbmin = xbbmin
      xbmax = xbbmax
      zmin = zzmin
      zmax = zzmax

      fuu = qgauss(fun_DIS_xzQ, Q2min, Q2max, 3)

      return
      end

c-----

      function fun_DIS_xzQ(Q2)
      implicit none
      real*8 Q2
      real*8 fun_DIS_xzQ
      real*8 rts, qt, xbmin, xbmax, zmin, zmax
      real*8 fuu
      integer IT
      integer IH,IC
      common /xzQ2intsid/ rts, qt, xbmin, xbmax, zmin, zmax, IT
      common /meson/ IH,IC

      call xzintDISUU(rts,Q2, xbmin, xbmax, zmin, zmax, qt, fuu,
     >                     IT,IH,IC)
      fun_DIS_xzQ = fuu

      return
      end

c--------------------------------------------------------------
c     DIS X-section in {Q2,qt}: d[sigma]/dQ^2 d^2qt
c     integrated in xb, z
c--------------------------------------------------------------
      subroutine xzintDISUU(rtss,QQ2,xbbmin,xbbmax,
     >                     zzmin,zzmax,qtt,fuu,
     >                     ITT,IHH,ICC)
      implicit none
      real*8 rtss,QQ2,xbbmin,xbbmax,zzmin,zzmax,qtt
      integer ITT,IHH,ICC
      real*8 fuu
      real*8 rts,Q2,zh,qt,zmin,zmax
      integer IT ,IH ,IC
      common /xzintsid/ rts,Q2,qt,zmin,zmax,IT,IH,IC
      real*8 qgauss,fun_DIS_xz
      external qgauss,fun_DIS_xz

      rts= rtss
      Q2 = QQ2
      qt = qtt
      zmin = zzmin
      zmax = zzmax
      IT = ITT
      IH = IHH
      IC = ICC

      fuu = qgauss(fun_DIS_xz,xbbmin,xbbmax,3)

      return
      end

c-----

      function fun_DIS_xz(x)
      implicit none
      real*8 x
      real*8 fun_DIS_xz
      real*8 rts,Q2,qt,zmin,zmax
      integer IT,IH,IC
      common /xzintsid/ rts,Q2,qt,zmin,zmax,IT,IH,IC

      call zintDISUU(rts,Q2,x,zmin,zmax,qt,fun_DIS_xz
     >              ,IT,IH,IC)

      return
      end

c--------------------------------------------------------------
c     DIS X-section in {zh,Q2,qt}: d[sigma]/dz_HdQ^2d^2qt
c     integrated in xb
c--------------------------------------------------------------
      subroutine xintDISUU(rtss,QQ2,xbbmin,xbbmax,
     >                     zhh,qtt,fuu,
     >                     ITT,IHH,ICC)
      implicit none
      real*8 rtss,QQ2,xbbmin,xbbmax,zhh,qtt
      integer ITT,IHH,ICC
      real*8 fuu
      real*8 rts,Q2,zh,qt
      integer IT,IH,IC
      real*8   qgauss,fun_sidis_xb
      external qgauss,fun_sidis_xb

      rts= rtss
      Q2 = QQ2
      zh = zhh
      qt = qtt
      IT = ITT
      IH = IHH
      IC = ICC

      fuu = qgauss(fun_sidis_xb,xbbmin,xbbmax,3)

      return
      end

c------

      function fun_sidis_xb(x)
      implicit none
      real*8 fun_sidis_xb
      real*8 x
      real*8 rts,Q2,zh,qt,fuu
      integer IT,IH,IC
      common /xbintsid/ rts,Q2,zh,qt,IT,IH,IC

      call DISUU(rts,Q2,x,zh,qt,fuu,IT,IH,IC)
      fun_sidis_xb=fuu

      return
      end

c--------------------------------------------------------------
c     DIS X-section in {xb,Q2,qt}: d[sigma]/dx_BdQ^2d^2qt
c     integrated in zh
c--------------------------------------------------------------
      subroutine zintDISUU(rtss,QQ2,xbb,zzmin,zzmax,qtt,fuu,
     >                     ITT,IHH,ICC)
      implicit none
      real*8 rtss,QQ2,xbb,zzmin,zzmax,qtt
      real*8 fuu
      integer ITT,IHH,ICC
      real*8 rts,Q2,xb,qt
      integer IT,IH,IC
      common /zhintsid/ rts,Q2,xb,qt,IT,IH,IC
      real*8 qgauss,fun_DIS_zh
      external qgauss,fun_DIS_zh

      rts= rtss
      Q2 = QQ2
      xb = xbb
      qt = qtt
      IT = ITT
      IH = IHH
      IC = ICC

      fuu = qgauss(fun_DIS_zh,zzmin,zzmax,3)

      return
      end

c------

      function fun_DIS_zh(z)
      implicit none
      real*8 z
      real*8 fun_DIS_zh
      real*8 rts,Q2,xb,qt
      integer IT,IH,IC
      common /zhintsid/ rts,Q2,xb,qt,IT,IH,IC

      call DISUU(rts,Q2,xb,z,qt,fun_DIS_zh,IT,IH,IC)

      return
      end

c--------------------------------------------------------------
c     DIS X-section in {xb,Q2,z,qt}: d[sigma]/dx_BdQ^2dzd^2qt
c--------------------------------------------------------------
      subroutine DISUU(rtss,QQ2,xbb,zz,qtt,fuu,ITT,IHH,ICC)
      implicit none
      real*8 rtss,QQ2,xbb,zz,qtt,fuu,bdiv,BdepDIS,qgauss
      real*8 rts,Q2,xb,z,qt,fut,AN,BdepDIS_Siv
      real*8 BdepDIS2,BdepDIS_Siv2,pt
      real*8 fut_res,fuu_res
      integer ngauss,nbel,nbuu,nbut,nQ,nY,ITT,IT,IH,IC,IHH,ICC
      common /bintdis/ rts,Q2,xb,z,qt,IT
      common /ngbesl/bdiv,nbel,nbuu,nbut,nQ,nY
      common /meson/ IH,IC
      integer fft
      common /FT/ fft
      external BdepDIS_Siv,BdepDIS,BdepDIS_Siv2,BdepDIS2,qgauss
      external BdepDIS_og, BdepDIS_Siv_og,TMD_FF2
      real*8 Q,Qi,b,bstar,Revo,Revo1,sig
      integer i

      rts= rtss
      Q2 = QQ2
      xb = xbb
      z  = zz
      qt = qtt
      IT = ITT
      IH = IHH
      IC = ICC

      call adogt(BdepDIS    ,qt,dsqrt(Q2),0,z,60d0,fuu)

      return
      end

c--------------------------------------------------------------
      function BdepDIS(b)
      implicit none
      real*8 BdepDIS,rts,Q2,xb,z,qt,b,sig
      integer IT,IH,IC
      common /bintdis/ rts,Q2,xb,z,qt,IT
      common /meson/ IH,IC

      call UU_DIS(rts,Q2,xb,z,qt,b,sig,IT,IH,IC)
      BdepDIS = sig

      return
      end

c--------------------------------------------------------------
c     DIS: F_{UU} in b-space (with evolution kernel R(b))
c--------------------------------------------------------------
      subroutine UU_DIS(rts,Q2,xb,z,qt,b,sig,IT,IH,IC)
      implicit none
      real*8 rts,Q2,Q,y,qt,sig,z,xb,b,Revo,Q0,bstar,eu,ed
      real*8 fu,fub,fd,fdb,fs,fsb,fg,bessel0,alfem,Qini
      real*8 U,D,UB,DB,SS,SB,GL,sigma0
      real*8 fu1,fub1,fd1,fdb1,fs1,fsb1,fg1
      real*8 fu2,fub2,fd2,fdb2,fs2,fsb2,fg2
      real*8 fun,fubn,fdn,fdbn,fsn,fsbn,fgn
      real*8 U1,D1,UB1,DB1,SS1,SB1,GL1
      real*8 U2,D2,UB2,DB2,SS2,SB2,GL2
      real*8 Un,Dn,UBn,DBn,SSn,SBn,GLn
      real*8 hlo,hnlo
      real*8 Vud,Vus,Mw,GF,sw,Vzu,Vzds,Azu,Azds,Mz,pi,Nc,convert,DIS,FNP
      real*8 c0,c02,CF,CA,mb,alam4,alam5,bmax,kappa1,kappa2,ktw,ptw
      integer IT,nloops,hop,nll,pre
      integer IH,IC
      integer A_FE
      real*8 Q02, Qaux
      REAL*8 widpdf, widff, widuu
      real*8 fourier
      real*8 fc, fb
      INTEGER REPLICA_NUM
      COMMON /REPLICA_INFO/ REPLICA_NUM
      common /Wpar/ Vud,Vus,Mw,GF,sw,Vzu,Vzds,Azu,Azds,Mz,Nc,convert
      common /cons/ c0,c02,pi,CF,CA,mb,alam4,alam5,
     >     bmax,ktw,ptw,Qini,kappa2
      common /scheme/ nloops,hop,nll,pre
      REAL*8 aN, bN, g2A, A2, g2B, B2
      COMMON /FITP/ aN, bN, g2A, A2, g2B, B2
      real*8 alphas,hQ
      common /hardQ/ hQ
      integer TIH,TIC
      common /ttarget/ TIH,TIC
      integer fft
      common /FT/ fft

*-----b* prescription
      bstar = b/dsqrt(1d0+(b/bmax)**2d0)
      Q0 = c0/bstar
      Q02 = Q0**2

*-----Hard scale
      Q = dsqrt(Q2)
      hQ = Q

*-----Perturbative evolution
      call kernel_q(bstar,Q0,hQ,Q0,hQ,nll,Revo) ! 1/2/3/4 for LL/NLL/NNLL/NNNLL

*-----Point-like scattering cross section
      alfem = 1d0/137d0
      y = Q2/(xb*rts)
      sigma0 = 2d0*pi*alfem**2d0/(Q2*Q2)*(1d0+(1d0-y)**2d0)/(z*z)

*-----Precalculated TMDFF
      !call PreFF(z,b,IH,IC,fu,fd,fs,fub,fdb,fsb)
      IF(REPLICA_NUM.le.36) THEN
      pre =1
      ELSE
      pre = 0
      ENDIF


*-----Calculated TMDFF
      if (IT.eq.3) then
        call Pre_FF_A(z,Q0,IT,fU,fD,fUB,fDB,fS,fSB,fG)
!        CALL FF(z,Q0,fu,fub,fd,fdb,fs,fsb,fg,nloops)
        !print*, 'FF: ', xb, Q0, U,D,UB,DB,SS,SB,GL
      else
         if (pre.eq.1) then
         call Pre_FF_A(z,Q0,IT,fU,fD,fUB,fDB,fS,fSB,fG)
         else
         ! note: code for NLO needs to be updated, use only at LO:
         CALL LIKEnFF(z,Q0,IT,fU,fUB,fD,fDB,fS,fSB,fG,nloops)
         !call Pre_FF_A(z,Q0,IT,fU,fD,fUB,fDB,fS,fSB,fG)
         endif
      endif

      pre = 0
      
*-----TMDPDF
      if(IT.eq.1) then          ! proton target
         TIH = 4
         TIC = 1
*--------Precalculated TMDPDF for proton
         !call getTMDPDFP(xb,b,u1,d1,ss1,ub1,db1,sb1)
         u = u1
         d = d1
         ss= ss1
         ub= ub1
         db= db1
         sb= sb1
*--------Calculated TMDPDF for proton
         call PDF_p(xb,Q0,U,D,UB,DB,SS,SB,GL,nloops)

      elseif(IT.eq.2) then      ! neutron target
         TIH = 4
         TIC = 1
         !call getTMDPDFP(xb,b,u1,d1,ss1,ub1,db1,sb1)
       	 u = d1
         d = u1
       	 ss= ss1
       	 ub= db1
       	 db= ub1
       	 sb= sb1
         call PDF_n(xb,Q0,U,D,UB,DB,SS,SB,GL,nloops)

      elseif(IT.eq.3) then      ! deuteron target
         TIH = 4
         TIC = 1
         !call getTMDPDFP(xb,b,u1,d1,ss1,ub1,db1,sb1)
       	 u = (u1+d1)/2d0
         d = u
       	 ss= ss1
       	 ub= (ub1+db1)/2d0
       	 db= ub
       	 sb= sb1
         call PDF_D(xb,Q0,U,D,UB,DB,SS,SB,GL,nloops)
         !print*, 'PDF: ', xb, Q0, U,D,UB,DB,SS,SB,GL
      elseif(IT.eq.4) then        ! helium target
          TIH = 4
          TIC = 1
          !call getTMDPDFP(xb,b,u1,d1,ss1,ub1,db1,sb1)
          u = (u1+d1)/2d0
          d = u
          ss= ss1
          ub= (ub1+db1)/2d0
          db= ub
          sb= sb1
          !if (xb.lt.0.355d0) then
            if (pre.eq.1) then
            call Pre_PDF_A(xb,Q0,IT,U,D,UB,DB,SS,SB,GL)
            else
            call PDF_A_LHAPDF(xb,Q0,IT,2,U,D,UB,DB,SS,SB,GL,nloops)
            !                       (A,Z)
            endif
      elseif(IT.eq.12) then        ! carbon target
        TIH = 4
        TIC = 1
        !call getTMDPDFP(xb,b,u1,d1,ss1,ub1,db1,sb1)
        u = (u1+d1)/2d0
        d = u
        ss= ss1
        ub= (ub1+db1)/2d0
        db= ub
        sb= sb1
        !if (xb.lt.0.355d0) then
          if (pre.eq.1) then
          call Pre_PDF_A(xb,Q0,IT,U,D,UB,DB,SS,SB,GL)
          else
          call PDF_A_LHAPDF(xb,Q0,IT,6,U,D,UB,DB,SS,SB,GL,nloops)
          !                      (A,Z)
          endif
      elseif(IT.eq.20) then        ! neon target
          TIH = 4
          TIC = 1
          !call getTMDPDFP(xb,b,u1,d1,ss1,ub1,db1,sb1)
         u = (u1+d1)/2d0
          d = u
          ss= ss1
          ub= (ub1+db1)/2d0
          db= ub
          sb= sb1
        !  if (xb.lt.0.355d0) then
            if (pre.eq.1) then
            call Pre_PDF_A(xb,Q0,IT,U,D,UB,DB,SS,SB,GL)
            else
            call PDF_A_LHAPDF(xb,Q0,IT,10,U,D,UB,DB,SS,SB,GL,nloops)
            !                      (A,Z)
            endif
      elseif(IT.eq.56) then        ! iron target
          TIH = 4
          TIC = 1
          !call getTMDPDFP(xb,b,u1,d1,ss1,ub1,db1,sb1)
          u = (u1+d1)/2d0
          d = u
          ss= ss1
          ub= (ub1+db1)/2d0
          db= ub
          sb= sb1
          A_FE = 57
          !if (xb.lt.0.355d0) then
          if (pre.eq.1) then
            call Pre_PDF_A(xb,Q0,A_FE,U,D,UB,DB,SS,SB,GL)
          else
            call PDF_A_LHAPDF(xb,Q0,IT,26,U,D,UB,DB,SS,SB,GL,nloops)
            !                      (A,Z)
          endif
      elseif(IT.eq.84) then        ! krypton target
          TIH = 4
          TIC = 1
          !call getTMDPDFP(xb,b,u1,d1,ss1,ub1,db1,sb1)
          u = (u1+d1)/2d0
          d = u
          ss= ss1
          ub= (ub1+db1)/2d0
          db= ub
          sb= sb1
          if (pre.eq.1) then
       	  call Pre_PDF_A(xb,Q0,IT,U,D,UB,DB,SS,SB,GL)
          else
          call PDF_A_LHAPDF(xb,Q0,IT,36,U,D,UB,DB,SS,SB,GL,nloops)
          endif
      elseif(IT.eq.131) then        ! xenon target
          TIH = 4
          TIC = 1
          !call getTMDPDFP(xb,b,u1,d1,ss1,ub1,db1,sb1)
          u = (u1+d1)/2d0
          d = u
          ss= ss1
          ub= (ub1+db1)/2d0
          db= ub
          sb= sb1
          if (pre.eq.1) then
       	  call Pre_PDF_A(xb,Q0,IT,U,D,UB,DB,SS,SB,GL)
          else
          call PDF_A_LHAPDF(xb,Q0,IT,54,U,D,UB,DB,SS,SB,GL,nloops)
          endif
        elseif(IT.eq.197) then        ! Gold target
            TIH = 4
            TIC = 1
            if (pre.eq.1) then
            call Pre_PDF_A(xb,Q0,IT,U,D,UB,DB,SS,SB,GL)
            else
            call PDF_A_LHAPDF(xb,Q0,IT,79,U,D,UB,DB,SS,SB,GL,nloops)
            endif
      elseif(IT.eq.208) then        ! lead target
          TIH = 4
          TIC = 1
          !call getTMDPDFP(xb,b,u1,d1,ss1,ub1,db1,sb1)
          u = (u1+d1)/2d0
          d = u
          ss= ss1
          ub= (ub1+db1)/2d0
          db= ub
          sb= sb1
          if (pre.eq.1) then
          call Pre_PDF_A(xb,Q0,IT,U,D,UB,DB,SS,SB,GL)
          else
          call PDF_A_LHAPDF(xb,Q0,IT,82,U,D,UB,DB,SS,SB,GL,nloops)
          endif
       endif

      pre = 1 !!! SET TO 1 FOR FIT

*-----EM charges for quarks
      eu = 4d0/9d0
      ed = 1d0/9d0
      if (Q.lt.1.59d0) then
      Qaux = 1.59d0
      else
      Qaux = Q
      endif

*-----Hard function from 1107.3973
      if (nloops.eq.2) then
      sigma0 = sigma0*(1d0+alphas(hQ)/2d0/pi*CF*(3d0*dlog(Q2/hQ/hQ)
     $                                    -dlog(Q2/hQ/hQ)**2d0-8d0))
      endif


*-----Sum over flavors
      DIS=sigma0*(eu*(U*fu+UB*fub)+ed*(D*fd+DB*fdb+SS*fs+SB*fsb))

*-----NP evolution from 1406.3073v2 for proton
*---- with Modifications to account for Nuclear Medium Effects:
      if (IT.eq.3) then
      widpdf = ktw+aN*(2.**(1d0/3d0)-1d0) !*(xb**g2a)*((1-xb)**A2)
      widff  = ptw+bN*(2.**(1d0/3d0)-1d0) !*(z**a2) !*((1-z)**B2)
      widuu = kappa2*dlog(Q/Qini)
C     >        + 0.5*gD*(2.**(1d0/3d0)-1d0)*((Qaux/Qini)**aD)
C     >        + 0.5*gP*(2.**(1d0/3d0)-1d0)*((Qaux/Qini)**aP)
      else
      widpdf = ktw+aN*(IT**(1d0/3d0)-1d0) !*(xb**g2a)*((1-xb)**A2)
      widff  = ptw+bN*(IT**(1d0/3d0)-1d0) !*(z**a2) !*((1-z)**B2)
      widuu = kappa2* dlog(Q/Qini)
C     >        + 0.5*gD*(IT**(1d0/3d0)-1d0)*((Qaux/Qini)**aD)
C     >        + 0.5*gP*(IT**(1d0/3d0)-1d0)*((Qaux/Qini)**aP)
      endif

      kappa1 = widpdf/4 + (widff/4/(z*z))

      FNP = dexp( -(kappa1)*b*b  - dlog(b/bstar)*widuu)


*-----Cross section
      sig = DIS*b*(Revo)**2d0*FNP

      return
      end
