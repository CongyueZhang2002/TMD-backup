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
      call adogt(qFFb,kt/ x,Q,0,1d0,199d0,part)

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

*     Input COMMON BLOCKS !2023
      double precision H_AA
      COMMON /HMASS/ H_AA 

C-------QCDNUM initialization
      DOUBLE PRECISION :: array(47), def(-6:6,12), qq(2),wt(2)
      DOUBLE PRECISION :: pdf(-6:6) 
      Double precision as0, r20, xmin, eps 
      integer :: iord, nfin, itype, iset, jset, iosp, nx
      integer :: nxin, nqin, lun, idbug, iqc, iqb, iq0, nq 
      double precision :: q2c, q2b, q0 
      double precision Qf 
      double precision Qf2  

C------------- X grid and mu2 grid paramaters
      data idbug/0/                                     !debug printpout
      data xmin/1.D-4/, nxin/100/                       !x grid, 
      data qq/1.D0,1.D5/, wt/1.D0,1.D0/, nqin/60/       !mu2 grid


C------ Z-array for output file 
      DATA array/0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
     4        0.095, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
     5        0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475,  0.5,
     6        0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7,
     7        0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9,
     8        0.925, 0.95, 0.975, 1.0/

C--------------------------------------------------------------   
      external func                       !input parton dists
      integer, external :: iqfrmq
      data def  /                         !flavor composition
C--   tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
C--   -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6  
     + 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.,         !d
     + 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.,         !u
     + 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0.,         !s
     + 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,         !dbar
     + 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,         !ubar
     + 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,         !sbar
     + 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,         !c
     + 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,         !cbar
     + 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,         !b
     + 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,         !bbar
     + 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,         !t
     + 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.  /       !tbar 
                                       !pdfout

C--   Weight files
      character*26 fnam(3)

      data fnam /'weights/unpolarised.wgt',
     +           'weights/polarised.wgt  ',
     +           'weights/timelike.wgt   '/

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
            H_AA = 208
            lun = -6    ! -6 supresses banner
            call qcinit(lun,' ')                       ! initialize
            iosp = 2 ! 2 for linear interpolarion, 3 for spline       
            call gxmake(xmin,1,1,nxin,nx,iosp)                         !x-grid
            call gqmake(qq,wt,2,nqin,nq)                             !mu2-grid
            itype = 3  ! 1: Unpol PDF, 2: Pol PDF, 3: FF 
            call wtfile(itype,fnam(itype))   !calculate weights
            iord = 2  ! 1 : LO, 2 : NLO, 3 : NNLO
            call setord(iord)                                 
            as0 = 0.364d0 
            r20 = 2.0d0 
            call setalf(as0,r20)                              !input alphas
            q0 = 1.0D0 
            q2c = (1.43d0)**2d0 
            q2b = (4.3d0)**2d0 
            iqc  = iqfrmq(q2c)    !Charm threshold
            iqb  = iqfrmq(q2b)    !Bottom threshold
            nfin = 1 
            call setcbt(nfin,iqc,iqb,999)         ! Threshholds in the VFNS
            iset = 1 
            jset = 10*iset+itype            !output pdf set and evolution type
            iq0  = iqfrmq(q0)                !start scale index 
            call setint('edbg',idbug)        !debug printout
            call evolfg(jset,func,def,iq0,eps)  !evolve all pdf's
      
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