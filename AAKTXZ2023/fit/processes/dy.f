c--------------------------------------------------------------
c     DY cross section in {qt,Q} only: integrate over {y}
c     d[sigma]/dqt unpolarized only
c--------------------------------------------------------------

      subroutine DY_overy(Ixxfy,rtss,qqt,yymin,yymax,
     >                     QQ, fuu,ITT)
      implicit none
      real*8 rtss,qqt,qt,fuu,fut,bdiv,qgauss,rts,AN
      real*8 fun_DY_y
      real*8 yymin,yymax,ymin,ymax
      real*8 QQ, Q
      integer nbel,nbuu,nbut,nQ,nY,ITT,IT,Ixxfy,Ixfy
      common /yint2/ rts,qt,Q,ymin,ymax,IT,Ixfy
      common /ngbesl/bdiv,nbel,nbuu,nbut,nQ,nY
      external fun_DY_yu,qgauss

      rts= rtss
      qt = qqt
      ymin = yymin
      ymax = yymax
      Q = QQ
      IT = ITT
      Ixfy = Ixxfy

      fuu = qgauss(fun_DY_yu,ymin,ymax,nY)

      return
      end

c------

      function fun_DY_yu(y)
      implicit none
      real*8 xf, ymin, ymax
      real*8 fun_DY_yu,fuu,rts,y,qt,Q
      integer IT,Ixfy
      common /yint2/ rts,qt,Q,ymin,ymax,IT,Ixfy

      call DY_avg(Ixfy,rts,xf,y,Q,qt,fuu,IT)
      fun_DY_yu=fuu

      return
      end


c--------------------------------------------------------------
c     DY cross section in {qt,Q} only: integrate over {xF}
c     d[sigma]/dqt unpolarized only
c--------------------------------------------------------------

      subroutine DY_overxF(Ixxfy,rtss,qqt,xxfmin,xxfmax,
     >                     QQ, fuu,ITT)
      implicit none
      real*8 rtss,qqt,xxfmin,xxfmax,QQ,fuu
      integer Ixxfy,ITT
      real*8 rts,qt,Q,xfmin,xfmax
      integer IT,Ixfy
      common /XFint/ rts,Q,qt,Ixfy,IT
      real*8 qgauss,fun_DY_xF
      external qgauss,fun_DY_xF
      integer nbel,nbuu,nbut,nQ,nY
      real*8 bdiv
      common /ngbesl/bdiv,nbel,nbuu,nbut,nQ,nY

      rts= rtss
      qt = qqt
      xfmin = xxfmin
      xfmax = xxfmax
      Q = QQ
      IT = ITT
      Ixfy = Ixxfy

      fuu = qgauss(fun_DY_xf,xfmin,xfmax,nY)

      return
      end

c------

      function fun_DY_XF(xf)
      implicit none
      real*8 xf
      real*8 fun_DY_XF
      real*8 rts,Q,qt
      integer Ixfy,IT
      real*8 y
      common /XFint/ rts,Q,qt,Ixfy,IT

      call DY_avg(Ixfy,rts,xf,y,Q,qt,fun_DY_XF,IT)

      return
      end

c--------------------------------------------------------------
c     DY cross section in {qt} only: integrate over {x2,Q}
C     note in the notation here xb=x2.
c     d[sigma]/dqt unpolarized only
c--------------------------------------------------------------

      subroutine DY_overxbQ(Ixxfy,rtss,qqt,xbmin,
     >     xbmax,QQmin,QQmax,fuu,ITT)
      implicit none
      real*8 rtss,qqt,xbmin,xbmax,QQmin,QQmax,fuu
      integer Ixxfy,ITT
      real*8 rts,qt,Qmin,Qmax
      integer Ixfy,IT
      common /Xbint/ rts,qt,Qmin,Qmax,Ixfy,IT
      integer nbel,nbuu,nbut,nQ,nY
      real*8 bdiv
      common /ngbesl/bdiv,nbel,nbuu,nbut,nQ,nY
      real*8 qgauss,fun_DY_xb
      external qgauss,fun_DY_xb

      rts= rtss
      qt = qqt
      Qmin = QQmin
      Qmax = QQmax
      Ixfy = Ixxfy
      IT   = ITT

      fuu = qgauss(fun_DY_xb,xbmin,xbmax,nY)

      return
      end

      function fun_DY_xb(xb)
      implicit none
      real*8 xb
      real*8 fun_DY_xb
      real*8 rts,qt,Qmin,Qmax
      integer IT,Ixfy
      common /Xbint/ rts,qt,Qmin,Qmax,Ixfy,IT
      real*8 y

      call DY_overQ(Ixfy,rts,xb,y,qt,Qmin,Qmax,
     >              fun_DY_xb,IT)

      return
      end

c--------------------------------------------------------------
c     DY cross section in {qt} only: integrate over {y,Q}
c     d[sigma]/dqt unpolarized only
c--------------------------------------------------------------

      subroutine DY_overyu(Ixxfy,rtss,qqt,yymin,
     >                     yymax,
     >                     QQmin,QQmax,fuu,ITT)
      implicit none
      real*8 rtss,qqt,yymin,yymax,QQmin,QQmax
      integer Ixxfy,ITT
      real*8 fuu
      real*8 rts ,qt,              Qmin ,Qmax
      integer IT,Ixfy
      common /yint/ rts,qt,Qmin,Qmax,IT,Ixfy
      real*8 fun_DY_y,qgauss
      external fun_DY_y,qgauss
      real*8 bdiv
      integer nbel,nbuu,nbut,nQ,nY
      common /ngbesl/bdiv,nbel,nbuu,nbut,nQ,nY

      rts= rtss
      qt = qqt
      Qmin = QQmin
      Qmax = QQmax
      IT = ITT
      Ixfy = Ixxfy

      fuu = qgauss(fun_DY_y,yymin,yymax,nY)

      return
      end

c------

      function fun_DY_y(y)
      implicit none
      real*8 y
      real*8 fun_DY_y
      real*8 rts,xf,qt,Qmin,Qmax
      integer Ixfy,IT
      common /yint/ rts,qt,Qmin,Qmax,IT,Ixfy

      call DY_overQ(Ixfy,rts,xf,y,qt,Qmin,Qmax,fun_DY_y,IT)

      return
      end

c--------------------------------------------------------------
c     DY Xsection in {y,qt} or {xf,qt} only: integrate over {Q}
c     d[sigma]/dyd^2qt: convert to pb*GeV^{-2}
c     d[sigma]/dxfd^2qt: convert to pb*GeV^{-2}
c--------------------------------------------------------------

      subroutine DY_overQ(Ixxfy,rtss,xxf,yy,qtt,QQmin,QQmax,fuu,ITT)
      implicit none
      real*8 rtss,xxf,yy,qtt,QQmin,QQmax
      integer Ixxfy,ITT
      real*8 fuu
      real*8 rts ,xf ,y ,qt ,Qmin ,Qmax
      integer Ixfy ,IT
      real*8 fun_DY_Q,qgauss
      external fun_DY_Q,qgauss
      common /Qint/ rts,xf,y,qt,Ixfy,IT
      real*8 bdiv
      integer nbel,nbuu,nbut,nQ,nY
      common /ngbesl/bdiv,nbel,nbuu,nbut,nQ,nY

      rts = rtss
      xf = xxf
      y = yy
      qt = qtt
      Ixfy = Ixxfy
      IT = ITT

      fuu = qgauss(fun_DY_Q,QQmin,QQmax,nQ)

      return
      end

c-----

      function fun_DY_Q(Q)
      implicit none
      real*8 Q
      real*8 fun_DY_Q
      real*8 rts,xf,y,qt
      integer IT,Ixfy
      real*8 fuu
      common /Qint/ rts,xf,y,qt,Ixfy,IT

      call DY_avg(Ixfy,rts,xf,y,Q,qt,fuu,IT)
      fun_DY_Q=fuu * 2d0 * Q

      return
      end

c--------------------------------------------------------------
c     DY boson cross section in {y, Q,qt}: d[sigma]/dQ^2dyd^2qt
c                         or in {xf,Q,qt}: d[sigma]/dQ^2dxfd^2qt
c     convert to pb*GeV^{-2}
c     Also the asymmetry A_N in {y,Q,qt} or {xf,Q,qt} space
c--------------------------------------------------------------
      subroutine DY_avg(Ixxfy,rtss,xxf,yy,QQ,qtt,fuu,ITT)
      implicit none
      real*8 rtss,xxf,yy,QQ,qtt
      real*8 fuu
      integer Ixxfy,ITT
      real*8 rts ,xf ,y ,Q ,qt
      integer Ixfy ,IT
      real*8 BdepDY,qgauss
      external BdepDY,qgauss
      common /bintdy/ rts,xf,y,Q,qt,Ixfy,IT

      rts = rtss
      xf = xxf
      y  = yy
      Q = QQ
      qt = qtt
      Ixfy = Ixxfy
      IT = ITT

      call adogt(BdepDY,qt,Q,0,1d0,20d0,fuu)

      return
      end

c--------------------------------------------------------------
      function BdepDY(b)
      implicit none
      real*8 b
      real*8 BdepDY
      real*8 rts,xf,y,Q,qt
      integer IT,Ixfy
      common /bintdy/ rts,xf,y,Q,qt,Ixfy,IT

      call UU_DY(Ixfy,rts,xf,y,Q,qt,b,BdepDY,IT)

      return
      end

c--------------------------------------------------------------
c     DY: F_{UU} in b-space (with evolution kernel R(b))
c--------------------------------------------------------------
      subroutine UU_DY(Ixfy,rts,xf,y,Q,qt,b,sig,IT)
      use leptoncutsdy
      implicit none
      real*8 rts,Q,xf,y,qt,sig,xa,xb,b,Revo,Q0,bstar,eu,ed, x2, Qbar
      real*8 U ,D ,UB ,DB ,SS ,SB ,GL
      real*8 U1,D1,UB1,DB1,SS1,SB1,GL1,bessel0,alfem,Qini,norm
      real*8 U2,D2,UB2,DB2,SS2,SB2,GL2,sigma0,Wplus,Wminus,kappa3
      real*8 Vud,Vus,Mw,GF,sw,Vzu,Vzds,Azu,Azds,Mz,pi,Nc,convert,DY,FNP
      real*8 c0,c02,CF,CA,mb,alam4,alam5,bmax,kappa1,kappa2,ktw,ptw
      real*8 g1pi,FNPp,FNPN
      real*8 Q2,wid, widuu, Qaux, cor
      integer IT,Ixfy,nloops,hop,nll,prepdf,preff,ZZ,pre
      common /Wpar/ Vud,Vus,Mw,GF,sw,Vzu,Vzds,Azu,Azds,Mz,Nc,convert
      common /cons/ c0,c02,pi,CF,CA,mb,alam4,alam5,
     >     bmax,ktw,ptw,Qini,kappa2

      REAL*8 gamma, g3f, g3D, Nq1, gq1, dq1,
     &       Nq2, gq2, dq2, p_10, p_11, p_12 !2023 
      COMMON /FITP/ gamma, g3f, g3D, Nq1, gq1, dq1,
     &              Nq2, gq2, dq2, p_10, p_11, p_12 !2023

      common /scheme/ nloops,hop,nll,pre
      real*8 hQ,alphas
      common /hardQ/ hQ
      integer TIH,TIC
      common /ttarget/ TIH,TIC
      integer fft
      INTEGER REPLICA_NUM
      COMMON /REPLICA_INFO/ REPLICA_NUM
      common /FT/ fft
      integer cutflag
      real*8 ptcut,etamin,etamax,cutfac
      common /cuts/ ptcut,etamin,etamax,cutflag
      real*8 cw,GZ
      real*8 zfugg,zfugz,zfuzz,zfdgg,zfdgz,zfdzz
      real*8 zlgg,zlgz,zlzz
      real*8 b0,aa,mm
      real*8 AlphaQED,alpha_qed
      real*8 Ggamma

      if(cutflag.eq.1) then
         call SetCutParameters(ptcut,etamin,etamax)
         cutfac = CutFactor3(qt,Q,y)
      else
         cutfac = 1.0d0
      endif

*-----Hard scale
      hQ = Q
      Q2 = Q*Q

*-----Get xa and xb from input variables
      if(Ixfy.eq.1) then        ! {Q,y,qt} are variables
         xa = Q/rts*dexp(y)     ! xa is the proton beam
         xb = Q/rts*dexp(-y)    ! xb is the nuclear target
         norm = 1d0
      elseif(Ixfy.eq.2) then    ! {Q,xf,qt} are variables
         xa = ( xf+dsqrt(xf*xf+4d0*Q*Q/rts/rts))/2d0
         xb = (-xf+dsqrt(xf*xf+4d0*Q*Q/rts/rts))/2d0
         norm = 1d0/(xa+xb)
      elseif(Ixfy.eq.3) then        ! {Q,xa,qt} are variables
         xa = xf                    ! xf is xa for proton beam
         xb = (Q*Q/(rts*rts))/xa    ! xb belongs to nuclear target
         norm = 1d0/xf
       elseif(Ixfy.eq.4) then       ! {Q,xb,qt} are variables
         xb = xf                    ! xf is xb for nuclear target
         xa = (Q*Q/(rts*rts))/xb    ! xa belongs to proton beam
         norm = 1d0/xf
      endif

*-----b* presciption
      bstar = b/dsqrt(1d0+(b/bmax)**2d0)
      Q0 = c0/bstar

*-----Perturbative evoltuion
      call kernel_q(bstar,Q0,hQ,Q0,hQ,nll,Revo)

*-----Point-like scattering cross section
      alfem = alpha_qed(Q)
      sigma0 = norm*4d0*pi*alfem**2d0/(3d0*Q*Q*Nc)/rts/rts

*-----TMD PDF1
      if(Q0.le. 1.3d0) then
        Q0 = 1.31d0
      endif

      if(pre.eq.1) then
      call Pre_PDF_A(xa,Q0,1,U1,D1,UB1,DB1,SS1,SB1,GL1)
      else
      call PDF_p(xa,Q0,U1,D1,UB1,DB1,SS1,SB1,GL1,nloops)
      endif



*-----TMD PDF2
      if (pre.eq.1) then
      if(IT.eq.1) then          ! proton target
         !call PDF_p(xb,Q0,U2,D2,UB2,DB2,SS2,SB2,GL2,nloops)
         call Pre_PDF_A(xb,Q0,IT,U2,D2,UB2,DB2,SS2,SB2,GL2)
      elseif(IT.eq.2) then      ! neutron target
         !call PDF_n(xb,Q0,U2,D2,UB2,DB2,SS2,SB2,GL2,nloops)
         call Pre_PDF_A(xb,Q0, 1,U ,D ,UB ,DB ,SS ,SB ,GL )
         U2 = D
         D2 = U
         UB2= DB
         DB2= UB
         SS2= SS
         SB2= SB
      elseif(IT.eq.3) then      ! deuteron target
         !call PDF_D(xb,Q0,U2,D2,UB2,DB2,SS2,SB2,GL2,nloops)
         call Pre_PDF_A(xb,Q0, 1,U ,D ,UB ,DB ,SS ,SB ,GL )
       	 U2 = (D +U )/2.
       	 D2 = (U +D )/2.
       	 UB2= (DB+UB)/2.
       	 DB2= (UB+DB)/2.
       	 SS2= (SS+SB)/2.
       	 SB2= (SS+SB)/2.
      elseif(IT.eq.4) then      ! He nucleus target
         call Pre_PDF_A(xb,Q0,IT,U2,D2,UB2,DB2,SS2,SB2,GL2)
         !                      (A,Z)
      elseif(IT.eq.9) then      ! Be nucleus target
         call Pre_PDF_A(xb,Q0,IT,U2,D2,UB2,DB2,SS2,SB2,GL2)
         !                      (A,Z)
      elseif(IT.eq.12) then     ! C nucleus target
         call Pre_PDF_A(xb,Q0,IT,U2,D2,UB2,DB2,SS2,SB2,GL2)
         !                       (A,Z)
      elseif(IT.eq.40) then     ! Ca nucleus target
         call Pre_PDF_A(xb,Q0,IT,U2,D2,UB2,DB2,SS2,SB2,GL2)
         !                       (A,Z)
      elseif(IT.eq.56) then     ! Fe nucleus target
         call Pre_PDF_A(xb,Q0,IT,U2,D2,UB2,DB2,SS2,SB2,GL2)
         !                       (A,Z)
      elseif(IT.eq.184) then    ! W nucleus target
         call Pre_PDF_A(xb,Q0,IT,U2,D2,UB2,DB2,SS2,SB2,GL2)
         !                     (A,Z)
      elseif(IT.eq.197) then    ! Au nucleus target
         call Pre_PDF_A(xb,Q0,IT,U2,D2,UB2,DB2,SS2,SB2,GL2)
         !                     (A,Z)
      elseif(IT.eq.208) then    ! Pb nucleus target
         call Pre_PDF_A(xb,Q0,IT,U2,D2,UB2,DB2,SS2,SB2,GL2)
         !                     (A,Z)
      endif
      else
        if(IT.eq.1) then          ! proton target
       call PDF_p(xb,Q0,U2,D2,UB2,DB2,SS2,SB2,GL2,nloops)
      elseif(IT.eq.2) then	! neutron target
       call PDF_n(xb,Q0,U2,D2,UB2,DB2,SS2,SB2,GL2,nloops)
      elseif(IT.eq.3) then	! deuteron target
       call PDF_D(xb,Q0,U2,D2,UB2,DB2,SS2,SB2,GL2,nloops)
      elseif(IT.eq.9) then  ! Be target
          ZZ = 4
      call PDF_A_EPPS21(xb,Q0,IT,ZZ,U2,D2,UB2,DB2,SS2,SB2,GL2,nloops)
        !                      (A,Z)
      elseif(IT.eq.12) then ! C target
          ZZ = 6
      call PDF_A_EPPS21(xb,Q0,IT,ZZ,U2,D2,UB2,DB2,SS2,SB2,GL2,nloops)
        !                      (A,Z)
      elseif(IT.eq.40) then ! Ca target
          ZZ = 20
      call PDF_A_EPPS21(xb,Q0,IT,ZZ,U2,D2,UB2,DB2,SS2,SB2,GL2,nloops)
        !                      (A,Z)
      elseif(IT.eq.56) then ! Fe target
          ZZ = 26
      call PDF_A_EPPS21(xb,Q0,IT,ZZ,U2,D2,UB2,DB2,SS2,SB2,GL2,nloops)
        !                      (A,Z)
      elseif(IT.eq.184) then ! W target
          ZZ = 74
      call PDF_A_EPPS21(xb,Q0,IT,ZZ,U2,D2,UB2,DB2,SS2,SB2,GL2,nloops)
        !                      (A,Z)
      elseif(IT.eq.197) then ! Au target
          ZZ = 79
      call PDF_A_EPPS21(xb,Q0,IT,ZZ,U2,D2,UB2,DB2,SS2,SB2,GL2,nloops)
        !                      (A,Z)
      elseif(IT.eq.208) then    ! Pb nucleus target
          ZZ = 82
      call PDF_A_EPPS21(xb,Q0,IT,ZZ,U2,D2,UB2,DB2,SS2,SB2,GL2,nloops)
      !                      (A,Z)
      endif
      endif

C      pre = 1 !!! SET TO 1 FOR FIT

*-----quark charge fractions
      eu = 2d0/3d0
      ed = 1d0/3d0

*-----Weinberg angle
      sw = dsqrt(0.2313d0)
      cw = dsqrt(1d0-sw*sw)

*-----Z boson width
      GZ = 2.5d0

*-----quark charge fractions
      zfugg = eu*eu
      zfdgg = ed*ed
      zfuzz = (1d0-4d0*eu*sw*sw+8d0*eu*eu*sw**4d0)/8d0/sw**2d0/cw**2d0
      zfdzz = (1d0-4d0*ed*sw*sw+8d0*ed*ed*sw**4d0)/8d0/sw**2d0/cw**2d0
      zfugz = eu*(1d0-4d0*eu*sw*sw)/4d0/sw/cw
      zfdgz = ed*(1d0-4d0*ed*sw*sw)/4d0/sw/cw

*-----quark charge fractions
      zlgg = 1d0
      zlzz = Q**4d0/((Q*Q-MZ*MZ)**2d0+GZ*GZ*MZ*MZ)
     >          *(1d0-4d0*sw*sw+8d0*sw**4d0)/8d0/sw**2d0/cw**2d0
      zlgz = 2d0*Q*Q*(Q*Q-MZ*MZ)/((Q*Q-MZ*MZ)**2d0+GZ*GZ*MZ*MZ)
     >          *(1d0-4d0*sw*sw)/4d0/sw/cw

*-----Hard scattering contribution from 1107.3973
      if (nloops.eq.2) then
      sigma0 = sigma0*(1d0+alphas(hQ)/2d0/pi*CF*(3d0*dlog(Q2/hQ/hQ)
     $                -dlog(Q2/hQ/hQ)**2d0+pi**2d0-8d0))
      endif



*-----Sum contributions from flavors
      DY=sigma0*(
     >   zfugg*zlgg*(U1*UB2+UB1*U2)+
     >   zfdgg*zlgg*(D1*DB2+DB1*D2+SS1*SB2+SB1*SS2)
     >  +zfuzz*zlzz*(U1*UB2+UB1*U2)+
     >   zfdzz*zlzz*(D1*DB2+DB1*D2+SS1*SB2+SB1*SS2)
     >  +zfugz*zlgz*(U1*UB2+UB1*U2)+
     >   zfdgz*zlgz*(D1*DB2+DB1*D2+SS1*SB2+SB1*SS2)
     >          )

      if (Q.lt.1.59d0) then
      Qaux = 1.59d0
      else
      Qaux = Q
      endif

*-----NP evolution from 1406.3073v2 for proton
*---- with Modifications to account for Nuclear Medium Effects:
      if (IT.eq.3) then
      wid = ktw 
      widuu = kappa2*dlog(Q/Qini)
      cor = 0.5*g3f*(2.**(1d0/3d0)-1d0)*b*b*(1/Q)**gamma
      else
      wid = ktw 
      widuu = kappa2*dlog(Q/Qini)
      cor = 0.5*g3f*(IT**(1d0/3d0)-1d0)*b*b*(1/Q)**gamma 
      endif

      FNPN=dexp(-wid/4d0*b*b-widuu/2d0*dlog(b/bstar) - cor)

*-----NP evolution from 1406.3073v2 for proton
      FNPp=dexp(-ktw/4d0*b*b-kappa2/2d0*dlog(b/bstar)*dlog(Q/Qini))

*-----Total NP evolution
      FNP = FNPp*FNPN

      sig = b*DY*(Revo)**(2d0)*FNP*convert*cutfac

      return
      end
