c----- gives -kt/M*f1T(x,kt)
      subroutine sivers_kt(kt,x,Q,iq,part)
      implicit none
      real*8 kt,x,Q
      integer iq
      real*8 part
      real*8 u ,ub ,d ,db ,s ,sb
      real*8 u1,ub1,d1,db1,s1,sb1
      real*8 xx,QQ
      real*8 M
      integer iqq
      common /bint/ xx,QQ,iqq
      real*8 qsivb

      xx = x
      QQ = Q
      iqq= iq

      write(1010,*) qsivb(1d0)
      call adogt(qsivb,kt,Q,1,1d0,199d0,part)

      return
      end

      function qsivb(b)
      implicit none
      real*8 qsivb
      real*8 b,x,Q
      integer iq
      real*8 u,ub,d,db,s,sb
      real*8 xx,QQ
      integer iqq
      common /bint/ xx,QQ,iqq

      x = xx
      Q = QQ
      iq= iqq

      call sivers_b(b,x,Q,iq,qsivb)

      return
      end

      subroutine sivers_b(b,x,Q,iq,part)
      implicit none
      real*8 b,x,Q
      integer iq
      real*8 part
      real*8 parts
      real*8 u ,ub ,d ,db ,s ,sb
      real*8 us,ubs,ds,dbs,ss,sbs
      real*8 Nus ,Aus ,Bus ,
     >       Nds ,Ads ,Bds ,
     >       Nss ,Ass ,Bss ,
     >       Nubs,Aubs,Bubs,
     >       Ndbs,Adbs,Bdbs,
     >       Nsbs,Asbs,Bsbs,
     >       ktss,gams
      common /fitp/
     >       Nus ,Aus ,Bus ,
     >       Nds ,Ads ,Bds ,
     >       Nss ,Ass ,Bss ,
     >       Nubs,Aubs,Bubs,
     >       Ndbs,Adbs,Bdbs,
     >       Nsbs,Asbs,Bsbs,
     >       ktss,gams
      real*8 c0,c02,pi,CF,CA,mb,alam4,alam5,
     >     bmax,ktw,ptw,Qini,kappa2
      common /cons/ c0,c02,pi,CF,CA,mb,alam4,alam5,
     >     bmax,ktw,ptw,Qini,kappa2
      integer nloops,hop,nll
      common /scheme/ nloops,hop,nll
      real*8 bstar,Q0,Revo,FNP,kappa1,pref

      bstar = b/dsqrt(1d0+(b/bmax)**2d0)
      Q0 = c0/bstar
      call kernel_q(bstar,Q0,Q,Q0,Q,nll,Revo)

      kappa1 = ktss/4d0
      FNP = dexp( -(kappa1)*b*b
     >        -kappa2/2d0*dlog(b/bstar)*dlog(Q/Qini))

      pref = -b*b*Revo*FNP/2d0

      if (iq.eq.1) then
      call sivers_mel_p_u(x,Q0,parts)
      elseif (iq.eq.2) then
      call sivers_mel_p_d(x,Q0,parts)
      elseif (iq.eq.3) then
      call sivers_mel_p_s(x,Q0,parts)
      elseif (iq.eq.-1) then
      call sivers_mel_p_ub(x,Q0,parts)
      elseif (iq.eq.-2) then
      call sivers_mel_p_db(x,Q0,parts)
      elseif (iq.eq.-3) then
      call sivers_mel_p_sb(x,Q0,parts)
      endif

      part  = pref*parts

      return
      end

      subroutine sivers_mel_p_u(x,Q,u)
      implicit none
      real*8 x,Q
      real*8 xx,QQ
      common /qsivint/ xx,QQ
      real*8 u,ub,d,db,s,sb
      real*8 siv_integrandu
      real*8 siv_integrandd
      real*8 siv_integrands
      real*8 siv_integrandub
      real*8 siv_integranddb
      real*8 siv_integrandsb
      external siv_integrandu
      external siv_integrandd
      external siv_integrands
      external siv_integrandub
      external siv_integranddb
      external siv_integrandsb
      real*8 ri,rf
      integer n
      double precision a,abserr,alist,area,area1,area12,area2,a1,a2,b,
     *  blist,b1,b2,dabs,defabs,defab1,defab2,dmax1,d1mach,elist,epmach,
     *  epsabs,epsrel,errbnd,errmax,error1,error2,erro12,errsum,f,
     *  resabs,result,rlist,uflow
      integer ier,iord,iroff1,iroff2,k,key,keyf,last,limit,maxerr,neval,
     *  nrmax
      dimension alist(10),blist(10),elist(10),iord(10),
     *  rlist(10)

      xx = x
      QQ = Q

      ri = 0.0
      rf = 10.
      epsabs = 0d0
      epsrel = 1d0
      key = 1
      limit = 1

      call dqage(siv_integrandu ,ri,rf,epsabs,epsrel,key,limit,u ,
     *   abserr,neval,ier,alist,blist,rlist,elist,iord,last)

      return
      end

      subroutine sivers_mel_p_d(x,Q,d)
      implicit none
      real*8 x,Q
      real*8 xx,QQ
      common /qsivint/ xx,QQ
      real*8 u,ub,d,db,s,sb
      real*8 siv_integrandu
      real*8 siv_integrandd
      real*8 siv_integrands
      real*8 siv_integrandub
      real*8 siv_integranddb
      real*8 siv_integrandsb
      external siv_integrandu
      external siv_integrandd
      external siv_integrands
      external siv_integrandub
      external siv_integranddb
      external siv_integrandsb
      real*8 ri,rf
      integer n
      double precision a,abserr,alist,area,area1,area12,area2,a1,a2,b,
     *  blist,b1,b2,dabs,defabs,defab1,defab2,dmax1,d1mach,elist,epmach,
     *  epsabs,epsrel,errbnd,errmax,error1,error2,erro12,errsum,f,
     *  resabs,result,rlist,uflow
      integer ier,iord,iroff1,iroff2,k,key,keyf,last,limit,maxerr,neval,
     *  nrmax
      dimension alist(10),blist(10),elist(10),iord(10),
     *  rlist(10)

      xx = x
      QQ = Q

      ri = 0.0
      rf = 10.
      epsabs = 0d0
      epsrel = 1d0
      key = 1
      limit = 1

      call dqage(siv_integrandd ,ri,rf,epsabs,epsrel,key,limit,d ,
     *   abserr,neval,ier,alist,blist,rlist,elist,iord,last)

      return
      end

      subroutine sivers_mel_p_s(x,Q,s)
      implicit none
      real*8 x,Q
      real*8 xx,QQ
      common /qsivint/ xx,QQ
      real*8 u,ub,d,db,s,sb
      real*8 siv_integrandu
      real*8 siv_integrandd
      real*8 siv_integrands
      real*8 siv_integrandub
      real*8 siv_integranddb
      real*8 siv_integrandsb
      external siv_integrandu
      external siv_integrandd
      external siv_integrands
      external siv_integrandub
      external siv_integranddb
      external siv_integrandsb
      real*8 ri,rf
      integer n
      double precision a,abserr,alist,area,area1,area12,area2,a1,a2,b,
     *  blist,b1,b2,dabs,defabs,defab1,defab2,dmax1,d1mach,elist,epmach,
     *  epsabs,epsrel,errbnd,errmax,error1,error2,erro12,errsum,f,
     *  resabs,result,rlist,uflow
      integer ier,iord,iroff1,iroff2,k,key,keyf,last,limit,maxerr,neval,
     *  nrmax
      dimension alist(10),blist(10),elist(10),iord(10),
     *  rlist(10)

      xx = x
      QQ = Q

      ri = 0.0
      rf = 10.
      epsabs = 0d0
      epsrel = 1d0
      key = 1
      limit = 1

      call dqage(siv_integrands ,ri,rf,epsabs,epsrel,key,limit,s ,
     *   abserr,neval,ier,alist,blist,rlist,elist,iord,last)

      return
      end

      subroutine sivers_mel_p_ub(x,Q,ub)
      implicit none
      real*8 x,Q
      real*8 xx,QQ
      common /qsivint/ xx,QQ
      real*8 u,ub,d,db,s,sb
      real*8 siv_integrandu
      real*8 siv_integrandd
      real*8 siv_integrands
      real*8 siv_integrandub
      real*8 siv_integranddb
      real*8 siv_integrandsb
      external siv_integrandu
      external siv_integrandd
      external siv_integrands
      external siv_integrandub
      external siv_integranddb
      external siv_integrandsb
      real*8 ri,rf
      integer n
      double precision a,abserr,alist,area,area1,area12,area2,a1,a2,b,
     *  blist,b1,b2,dabs,defabs,defab1,defab2,dmax1,d1mach,elist,epmach,
     *  epsabs,epsrel,errbnd,errmax,error1,error2,erro12,errsum,f,
     *  resabs,result,rlist,uflow
      integer ier,iord,iroff1,iroff2,k,key,keyf,last,limit,maxerr,neval,
     *  nrmax
      dimension alist(10),blist(10),elist(10),iord(10),
     *  rlist(10)

      xx = x
      QQ = Q

      ri = 0.0
      rf = 10.
      epsabs = 0d0
      epsrel = 1d0
      key = 1
      limit = 1

      call dqage(siv_integrandub,ri,rf,epsabs,epsrel,key,limit,ub,
     *   abserr,neval,ier,alist,blist,rlist,elist,iord,last)

      return
      end

      subroutine sivers_mel_p_db(x,Q,db)
      implicit none
      real*8 x,Q
      real*8 xx,QQ
      common /qsivint/ xx,QQ
      real*8 u,ub,d,db,s,sb
      real*8 siv_integrandu
      real*8 siv_integrandd
      real*8 siv_integrands
      real*8 siv_integrandub
      real*8 siv_integranddb
      real*8 siv_integrandsb
      external siv_integrandu
      external siv_integrandd
      external siv_integrands
      external siv_integrandub
      external siv_integranddb
      external siv_integrandsb
      real*8 ri,rf
      integer n
      double precision a,abserr,alist,area,area1,area12,area2,a1,a2,b,
     *  blist,b1,b2,dabs,defabs,defab1,defab2,dmax1,d1mach,elist,epmach,
     *  epsabs,epsrel,errbnd,errmax,error1,error2,erro12,errsum,f,
     *  resabs,result,rlist,uflow
      integer ier,iord,iroff1,iroff2,k,key,keyf,last,limit,maxerr,neval,
     *  nrmax
      dimension alist(10),blist(10),elist(10),iord(10),
     *  rlist(10)

      xx = x
      QQ = Q

      ri = 0.0
      rf = 10.
      epsabs = 0d0
      epsrel = 1d0
      key = 1
      limit = 1

      call dqage(siv_integranddb,ri,rf,epsabs,epsrel,key,limit,db,
     *   abserr,neval,ier,alist,blist,rlist,elist,iord,last)

      return
      end

      subroutine sivers_mel_p_sb(x,Q,sb)
      implicit none
      real*8 x,Q
      real*8 xx,QQ
      common /qsivint/ xx,QQ
      real*8 u,ub,d,db,s,sb
      real*8 siv_integrandu
      real*8 siv_integrandd
      real*8 siv_integrands
      real*8 siv_integrandub
      real*8 siv_integranddb
      real*8 siv_integrandsb
      external siv_integrandu
      external siv_integrandd
      external siv_integrands
      external siv_integrandub
      external siv_integranddb
      external siv_integrandsb
      real*8 ri,rf
      integer n
      double precision a,abserr,alist,area,area1,area12,area2,a1,a2,b,
     *  blist,b1,b2,dabs,defabs,defab1,defab2,dmax1,d1mach,elist,epmach,
     *  epsabs,epsrel,errbnd,errmax,error1,error2,erro12,errsum,f,
     *  resabs,result,rlist,uflow
      integer ier,iord,iroff1,iroff2,k,key,keyf,last,limit,maxerr,neval,
     *  nrmax
      dimension alist(10),blist(10),elist(10),iord(10),
     *  rlist(10)

      xx = x
      QQ = Q

      ri = 0.0
      rf = 10.
      epsabs = 0d0
      epsrel = 1d0
      key = 1
      limit = 1

      call dqage(siv_integrandsb,ri,rf,epsabs,epsrel,key,limit,sb,
     *   abserr,neval,ier,alist,blist,rlist,elist,iord,last)

      return
      end

      subroutine sivers_mel_p(x,Q,u,ub,d,db,s,sb)
      implicit none
      real*8 x,Q
      real*8 xx,QQ
      common /qsivint/ xx,QQ
      real*8 u,ub,d,db,s,sb
      real*8 siv_integrandu
      real*8 siv_integrandd
      real*8 siv_integrands
      real*8 siv_integrandub
      real*8 siv_integranddb
      real*8 siv_integrandsb
      external siv_integrandu
      external siv_integrandd
      external siv_integrands
      external siv_integrandub
      external siv_integranddb
      external siv_integrandsb
      real*8 ri,rf
      integer n
      double precision a,abserr,alist,area,area1,area12,area2,a1,a2,b,
     *  blist,b1,b2,dabs,defabs,defab1,defab2,dmax1,d1mach,elist,epmach,
     *  epsabs,epsrel,errbnd,errmax,error1,error2,erro12,errsum,f,
     *  resabs,result,rlist,uflow
      integer ier,iord,iroff1,iroff2,k,key,keyf,last,limit,maxerr,neval,
     *  nrmax
      dimension alist(10),blist(10),elist(10),iord(10),
     *  rlist(10)

      xx = x
      QQ = Q

      ri = 0.0
      rf = 10.
      epsabs = 0d0
      epsrel = 1d0
      key = 1
      limit = 1

      call dqage(siv_integrandub,ri,rf,epsabs,epsrel,key,limit,ub,
     *   abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      call dqage(siv_integranddb,ri,rf,epsabs,epsrel,key,limit,db,
     *   abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      call dqage(siv_integrandsb,ri,rf,epsabs,epsrel,key,limit,sb,
     *   abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      call dqage(siv_integrandu ,ri,rf,epsabs,epsrel,key,limit,u ,
     *   abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      call dqage(siv_integrandd ,ri,rf,epsabs,epsrel,key,limit,d ,
     *   abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      call dqage(siv_integrands ,ri,rf,epsabs,epsrel,key,limit,s ,
     *   abserr,neval,ier,alist,blist,rlist,elist,iord,last)

      return
      end

      subroutine sivers_mel_n(x,Q,u,ub,d,db,s,sb)
      implicit none
      real*8 x,Q
      real*8 xx,QQ
      common /qsivint/ xx,QQ
      real*8 u ,ub ,d ,db ,s ,sb
      real*8 up,ubp,dp,dbp,sp,sbp
      real*8 siv_integrandu
      real*8 siv_integrandd
      real*8 siv_integrands
      real*8 siv_integrandub
      real*8 siv_integranddb
      real*8 siv_integrandsb
      external siv_integrandu
      external siv_integrandd
      external siv_integrands
      external siv_integrandub
      external siv_integranddb
      external siv_integrandsb
      real*8 ri,rf
      integer n
      double precision a,abserr,alist,area,area1,area12,area2,a1,a2,b,
     *  blist,b1,b2,dabs,defabs,defab1,defab2,dmax1,d1mach,elist,epmach,
     *  epsabs,epsrel,errbnd,errmax,error1,error2,erro12,errsum,f,
     *  resabs,result,rlist,uflow
      integer ier,iord,iroff1,iroff2,k,key,keyf,last,limit,maxerr,neval,
     *  nrmax
      dimension alist(10),blist(10),elist(10),iord(10),
     *  rlist(10)

      xx = x
      QQ = Q

      ri = 0.0
      rf = 10.
      epsabs = 0d0
      epsrel = 1d0
      key = 1
      limit = 1

      call dqage(siv_integrandu ,ri,rf,epsabs,epsrel,key,limit,up ,
     *   abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      call dqage(siv_integrandd ,ri,rf,epsabs,epsrel,key,limit,dp ,
     *   abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      call dqage(siv_integrands ,ri,rf,epsabs,epsrel,key,limit,sp ,
     *   abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      call dqage(siv_integrandub,ri,rf,epsabs,epsrel,key,limit,ubp,
     *   abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      call dqage(siv_integranddb,ri,rf,epsabs,epsrel,key,limit,dbp,
     *   abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      call dqage(siv_integrandsb,ri,rf,epsabs,epsrel,key,limit,sbp,
     *   abserr,neval,ier,alist,blist,rlist,elist,iord,last)

      d  = up
      u  = dp
      s  = sp
      db = ubp
      ub = dbp
      sb = sbp

      return
      end

      subroutine sivers_mel_d(x,Q,u,ub,d,db,s,sb)
      implicit none
      real*8 x,Q
      real*8 xx,QQ
      common /qsivint/ xx,QQ
      real*8 u ,ub ,d ,db ,s ,sb
      real*8 up,ubp,dp,dbp,sp,sbp
      real*8 siv_integrandu
      real*8 siv_integrandd
      real*8 siv_integrands
      real*8 siv_integrandub
      real*8 siv_integranddb
      real*8 siv_integrandsb
      external siv_integrandu
      external siv_integrandd
      external siv_integrands
      external siv_integrandub
      external siv_integranddb
      external siv_integrandsb
      real*8 ri,rf
      integer n
      double precision a,abserr,alist,area,area1,area12,area2,a1,a2,b,
     *  blist,b1,b2,dabs,defabs,defab1,defab2,dmax1,d1mach,elist,epmach,
     *  epsabs,epsrel,errbnd,errmax,error1,error2,erro12,errsum,f,
     *  resabs,result,rlist,uflow
      integer ier,iord,iroff1,iroff2,k,key,keyf,last,limit,maxerr,neval,
     *  nrmax
      dimension alist(10),blist(10),elist(10),iord(10),
     *  rlist(10)

      xx = x
      QQ = Q

      ri = 0.0
      rf = 10.
      epsabs = 0d0
      epsrel = 1d0
      key = 1
      limit = 1

      call dqage(siv_integrandu ,ri,rf,epsabs,epsrel,key,limit,up ,
     *   abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      call dqage(siv_integrandd ,ri,rf,epsabs,epsrel,key,limit,dp ,
     *   abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      call dqage(siv_integrands ,ri,rf,epsabs,epsrel,key,limit,sp ,
     *   abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      call dqage(siv_integrandub,ri,rf,epsabs,epsrel,key,limit,ubp,
     *   abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      call dqage(siv_integranddb,ri,rf,epsabs,epsrel,key,limit,dbp,
     *   abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      call dqage(siv_integrandsb,ri,rf,epsabs,epsrel,key,limit,sbp,
     *   abserr,neval,ier,alist,blist,rlist,elist,iord,last)

      d  = up +dp
      u  = dp +up
      s  = sp +dbp
      db = ubp+dbp
      ub = dbp+ubp
      sb = sbp+sp

      return
      end

      function siv_integrandsb(r)
      implicit none
      real*8 r
      real*8 siv_integrandsb
      real*8 x,Q
      real*8 c0
      common /qsivint/ x,Q
      double complex z
      real*8 pi,phi
      double complex Revo
      double complex siv_initsb
      double complex sbi
      integer nloops,hop,nll
      common /scheme/ nloops,hop,nll
      real*8 alphas

      pi = 4.*atan(1.)
      phi = 3./4.*pi
      c0 = 3.

      z = c0 + r*zexp(complex(0d0,phi))

      sbi = siv_initsb(z)
      if (nloops.eq.2) then
      sbi = sbi*(1d0-alphas(Q)/4d0/pi/3./z/(z+1d0))
      endif

      call evo(Q,z,Revo)

      siv_integrandsb = aimag(1./pi*zexp(complex(0.,phi))
     *    *x**(-z)*
     *    sbi*Revo)

      return
      end

      function siv_integranddb(r)
      implicit none
      real*8 r
      real*8 siv_integranddb
      real*8 x,Q
      common /qsivint/ x,Q
      double complex z
      real*8 pi,phi,c0
      double complex Revo
      double complex siv_initdb
      double complex dbi
      integer nloops,hop,nll
      common /scheme/ nloops,hop,nll
      real*8 alphas

      pi = 4.*atan(1.)
      phi = 3./4.*pi
      c0 = 3.

      z = c0 + r*zexp(complex(0d0,phi))

      dbi = siv_initdb(z)
      if (nloops.eq.2) then
      dbi = dbi*(1d0-alphas(Q)/4d0/pi/3./z/(z+1d0))
      endif

      call evo(Q,z,Revo)

      siv_integranddb = aimag(1./pi*zexp(complex(0.,phi))
     *    *x**(-z)*
     *    dbi*Revo)

      return
      end

      function siv_integrandub(r)
      implicit none
      real*8 r
      real*8 siv_integrandub
      real*8 x,Q
      common /qsivint/ x,Q
      double complex z
      real*8 pi,phi,c0
      double complex Revo
      double complex siv_initub
      double complex ubi
      integer nloops,hop,nll
      common /scheme/ nloops,hop,nll
      real*8 alphas

      pi = 4.*atan(1.)
      phi = 3./4.*pi
      c0 = 3.

      z = c0 + r*zexp(complex(0d0,phi))

      ubi = siv_initub(z)
      if (nloops.eq.2) then
      ubi = ubi*(1d0-alphas(Q)/4d0/pi/3./z/(z+1d0))
      endif

      call evo(Q,z,Revo)

      siv_integrandub = aimag(1./pi*zexp(complex(0.,phi))
     *    *x**(-z)*
     *    ubi*Revo)

      return
      end

      function siv_integrands(r)
      implicit none
      real*8 r
      real*8 siv_integrands
      real*8 x,Q
      common /qsivint/ x,Q
      double complex z
      real*8 pi,phi,c0
      double complex Revo
      double complex siv_inits
      double complex si
      integer nloops,hop,nll
      common /scheme/ nloops,hop,nll
      real*8 alphas

      pi = 4.*atan(1.)
      phi = 3./4.*pi
      c0 = 3.

      z = c0 + r*zexp(complex(0d0,phi))

      si = siv_inits(z)
      if (nloops.eq.2) then
      si = si*(1d0-alphas(Q)/4d0/pi/3./z/(z+1d0))
      endif

      call evo(Q,z,Revo)

      siv_integrands = aimag(1./pi*zexp(complex(0.,phi))
     *    *x**(-z)*
     *    si *Revo)

      return
      end

      function siv_integrandd(r)
      implicit none
      real*8 r
      real*8 siv_integrandd
      real*8 x,Q
      common /qsivint/ x,Q
      double complex z
      real*8 pi,phi,c0
      double complex Revo
      double complex siv_initd
      double complex di
      integer nloops,hop,nll
      common /scheme/ nloops,hop,nll
      real*8 alphas

      pi = 4.*atan(1.)
      phi = 3./4.*pi
      c0 = 3.

      z = c0 + r*zexp(complex(0d0,phi))

      di = siv_initd(z)
      if (nloops.eq.2) then
      di = di*(1d0-alphas(Q)/4d0/pi/3./z/(z+1d0))
      endif

      call evo(Q,z,Revo)

      siv_integrandd = aimag(1./pi*zexp(complex(0.,phi))
     *    *x**(-z)*
     *    di *Revo)

      return
      end

      function siv_integrandu(r)
      implicit none
      real*8 r
      real*8 siv_integrandu
      real*8 x,Q
      common /qsivint/ x,Q
      double complex z
      real*8 pi,phi,c0
      double complex Revo
      double complex siv_initu
      double complex ui
      integer nloops,hop,nll
      common /scheme/ nloops,hop,nll
      real*8 alphas

      pi = 4.*atan(1.)
      phi = 3./4.*pi
      c0 = 3.

      z = c0 + r*zexp(complex(0d0,phi))

      ui = siv_initu(z)
      if (nloops.eq.2) then
      ui = ui*(1d0-alphas(Q)/4d0/pi/3./z/(z+1d0))
      endif

      call evo(Q,z,Revo)

      siv_integrandu = aimag(1./pi*zexp(complex(0.,phi))
     *    *x**(-z)*
     *    ui *Revo)

      return
      end

      function siv_inits(z)
      implicit none
      double complex siv_inits
      double complex z,cgamma
      double complex u,ub,d,db,s,sb
      real*8 Auv,Buv,Cuv,Duv,Euv
      real*8 Adv,Bdv,Cdv,Ddv,Edv
      real*8 Aub,Bub,Cub,Dub,Eub
      real*8 Adb,Bdb,Cdb,Ddb,Edb
      real*8 fs
      real*8 Nus ,Aus ,Bus ,
     >       Nds ,Ads ,Bds ,
     >       Nss ,Ass ,Bss ,
     >       Nubs,Aubs,Bubs,
     >       Ndbs,Adbs,Bdbs,
     >       Nsbs,Asbs,Bsbs,
     >       ktss,gams
      common /fitp/
     >       Nus ,Aus ,Bus ,
     >       Nds ,Ads ,Bds ,
     >       Nss ,Ass ,Bss ,
     >       Nubs,Aubs,Bubs,
     >       Ndbs,Adbs,Bdbs,
     >       Nsbs,Asbs,Bsbs,
     >       ktss,gams

      Adb =  0.176
      Bdb = -0.172
      Cdb =  4.880
      Ddb =  0.000
      Edb =  0.000

      fs  =  0.4

      siv_inits = (Ass+Bss)**(Ass+Bss)
     *        /(Ass**Ass)/(Bss**Bss)*Nss*fs*(
     *    Adb*cgamma(complex(1.+Bss+Cdb,0.))*cgamma(-1.+Ass+Bdb+z)/
     *    cgamma(Ass+Bdb+Bss+Cdb+z) +
     *    (
     *    Adb*(Edb*(Ass+Bdb+z)+Ddb*(1.+Ass+Bdb+Bss+Cdb+z))*
     *    cgamma(complex(1.+Bss+Cdb,0.))*cgamma(Ass+Bdb+z)
     *    )/
     *    cgamma(2.+Ass+Bdb+Bss+Cdb+z)
     *          )

      return
      end

      function siv_initsb(z)
      implicit none
      double complex siv_initsb
      double complex z,cgamma
      double complex u,ub,d,db,s,sb
      real*8 Auv,Buv,Cuv,Duv,Euv
      real*8 Adv,Bdv,Cdv,Ddv,Edv
      real*8 Aub,Bub,Cub,Dub,Eub
      real*8 Adb,Bdb,Cdb,Ddb,Edb
      real*8 fs
      real*8 Nus ,Aus ,Bus ,
     >       Nds ,Ads ,Bds ,
     >       Nss ,Ass ,Bss ,
     >       Nubs,Aubs,Bubs,
     >       Ndbs,Adbs,Bdbs,
     >       Nsbs,Asbs,Bsbs,
     >       ktss,gams
      common /fitp/
     >       Nus ,Aus ,Bus ,
     >       Nds ,Ads ,Bds ,
     >       Nss ,Ass ,Bss ,
     >       Nubs,Aubs,Bubs,
     >       Ndbs,Adbs,Bdbs,
     >       Nsbs,Asbs,Bsbs,
     >       ktss,gams

      Adb =  0.176
      Bdb = -0.172
      Cdb =  4.880
      Ddb =  0.000
      Edb =  0.000

      fs  =  0.4

      siv_initsb = (Asbs+Bsbs)**(Asbs+Bsbs)
     *        /(Asbs**Asbs)/(Bsbs**Bsbs)*Nsbs*fs*(
     *    Adb*cgamma(complex(1.+Bsbs+Cdb,0.))*cgamma(-1.+Asbs+Bdb+z)/
     *    cgamma(Asbs+Bdb+Bsbs+Cdb+z) +
     *    (
     *    Adb*(Edb*(Asbs+Bdb+z)+Ddb*(1.+Asbs+Bdb+Bsbs+Cdb+z))*
     *    cgamma(complex(1.+Bsbs+Cdb,0.))*cgamma(Asbs+Bdb+z)
     *    )/
     *    cgamma(2.+Asbs+Bdb+Bsbs+Cdb+z)
     *          )

      return
      end

      function siv_initdb(z)
      implicit none
      double complex siv_initdb
      double complex z,cgamma
      double complex u,ub,d,db,s,sb
      real*8 Auv,Buv,Cuv,Duv,Euv
      real*8 Adv,Bdv,Cdv,Ddv,Edv
      real*8 Aub,Bub,Cub,Dub,Eub
      real*8 Adb,Bdb,Cdb,Ddb,Edb
      real*8 fs
      real*8 Nus ,Aus ,Bus ,
     >       Nds ,Ads ,Bds ,
     >       Nss ,Ass ,Bss ,
     >       Nubs,Aubs,Bubs,
     >       Ndbs,Adbs,Bdbs,
     >       Nsbs,Asbs,Bsbs,
     >       ktss,gams
      common /fitp/
     >       Nus ,Aus ,Bus ,
     >       Nds ,Ads ,Bds ,
     >       Nss ,Ass ,Bss ,
     >       Nubs,Aubs,Bubs,
     >       Ndbs,Adbs,Bdbs,
     >       Nsbs,Asbs,Bsbs,
     >       ktss,gams

      Adb =  0.176
      Bdb = -0.172
      Cdb =  4.880
      Ddb =  0.000
      Edb =  0.000

      fs  =  0.4

      siv_initdb = (Adbs+Bdbs)**(Adbs+Bdbs)
     *        /(Adbs**Adbs)/(Bdbs**Bdbs)*Ndbs*(1.-fs)*(
     *    Adb*cgamma(complex(1.+Bdbs+Cdb,0.))*cgamma(-1.+Adbs+Bdb+z)/
     *    cgamma(Adbs+Bdb+Bdbs+Cdb+z) +
     *    (
     *    Adb*(Edb*(Adbs+Bdb+z)+Ddb*(1.+Adbs+Bdb+Bdbs+Cdb+z))*
     *    cgamma(complex(1.+Bdbs+Cdb,0.))*cgamma(Adbs+Bdb+z)
     *    )/
     *    cgamma(2.+Adbs+Bdb+Bdbs+Cdb+z)
     *          )

      return
      end

      function siv_initub(z)
      implicit none
      double complex siv_initub
      double complex z,cgamma
      double complex u,ub,d,db,s,sb
      real*8 Auv,Buv,Cuv,Duv,Euv
      real*8 Adv,Bdv,Cdv,Ddv,Edv
      real*8 Aub,Bub,Cub,Dub,Eub
      real*8 Adb,Bdb,Cdb,Ddb,Edb
      real*8 fs
      real*8 Nus ,Aus ,Bus ,
     >       Nds ,Ads ,Bds ,
     >       Nss ,Ass ,Bss ,
     >       Nubs,Aubs,Bubs,
     >       Ndbs,Adbs,Bdbs,
     >       Nsbs,Asbs,Bsbs,
     >       ktss,gams
      common /fitp/
     >       Nus ,Aus ,Bus ,
     >       Nds ,Ads ,Bds ,
     >       Nss ,Ass ,Bss ,
     >       Nubs,Aubs,Bubs,
     >       Ndbs,Adbs,Bdbs,
     >       Nsbs,Asbs,Bsbs,
     >       ktss,gams

      Aub =  0.105
      Bub = -0.172
      Cub =  8.060
      Dub =  11.90
      Eub =  0.000

      siv_initub = (Aubs+Bubs)**(Aubs+Bubs)
     *        /(Aubs**Aubs)/(Bubs**Bubs)*Nubs*(
     *    Aub*cgamma(complex(1.+Bubs+Cub,0.))*cgamma(-1.+Aubs+Bub+z)/
     *    cgamma(Aubs+Bub+Bubs+Cub+z) +
     *    (
     *    Aub*(Eub*(Aubs+Bub+z)+Dub*(1.+Aubs+Bub+Bubs+Cub+z))*
     *    cgamma(complex(1.+Bubs+Cub,0.))*cgamma(Aubs+Bub+z)
     *    )/
     *    cgamma(2.+Aubs+Bub+Bubs+Cub+z)
     *          )

      return
      end

      function siv_initd(z)
      implicit none
      double complex siv_initd
      double complex z,cgamma
      double complex u,ub,d,db,s,sb
      real*8 Auv,Buv,Cuv,Duv,Euv
      real*8 Adv,Bdv,Cdv,Ddv,Edv
      real*8 Aub,Bub,Cub,Dub,Eub
      real*8 Adb,Bdb,Cdb,Ddb,Edb
      real*8 fs
      real*8 x,y
      double complex arg1,arg2,arg3,arg4,arg5
      double complex arg6,arg7,arg8,arg9,arg10
      real*8 Nus ,Aus ,Bus ,
     >       Nds ,Ads ,Bds ,
     >       Nss ,Ass ,Bss ,
     >       Nubs,Aubs,Bubs,
     >       Ndbs,Adbs,Bdbs,
     >       Nsbs,Asbs,Bsbs,
     >       ktss,gams
      common /fitp/
     >       Nus ,Aus ,Bus ,
     >       Nds ,Ads ,Bds ,
     >       Nss ,Ass ,Bss ,
     >       Nubs,Aubs,Bubs,
     >       Ndbs,Adbs,Bdbs,
     >       Nsbs,Asbs,Bsbs,
     >       ktss,gams

      Adv =  3.150
      Bdv =  0.806
      Cdv =  4.080
      Ddv =  0.000
      Edv =  0.000

      Adb =  0.176
      Bdb = -0.172
      Cdb =  4.880
      Ddb =  0.000
      Edb =  0.000

      fs  =  0.4

      siv_initd = (Ads+Bds)**(Ads+Bds)
     *        /(Ads**Ads)/(Bds**Bds)*Nds*(
     *    Adb*(1.-fs)*cgamma(complex(1.+Bds+Cdb,0.))
     *    *cgamma(-1.+Ads+Bdb+z)/
     *    cgamma(Ads+Bdb+Bds+Cdb+z) +
     *    (
     *    Adb*(1.-fs)*(Edb*(Ads+Bdb+z)+Ddb*(1.+Ads+Bdb+Bds+Cdb+z))*
     *    cgamma(complex(1.+Bds+Cdb,0.))*cgamma(Ads+Bdb+z)
     *    )/
     *    cgamma(2.+Ads+Bdb+Bds+Cdb+z)+
     *    Adv*cgamma(complex(1.+Bds+Cdv,0.))*cgamma(-1.+Ads+Bdv+z)/
     *    cgamma(Ads+Bds+Bdv+Cdv+z)+
     *    (
     *    Adv*(Edv*(Ads+Bdv+z)+Ddv*(1.+Ads+Bds+Bdv+Cdv+z))*
     *    cgamma(complex(1.+Bds+Cdv,0.))*cgamma(Ads+Bdv+z)
     *    )/
     *    cgamma(2.+Ads+Bds+Bdv+Cdv+z)
     *        )

      return
      end

      function siv_initu(z)
      implicit none
      real*8 x, y
      double complex siv_initu
      double complex z,cgamma
      double complex u,ub,d,db,s,sb
      double complex arg1,arg2,arg3,arg4,arg5
      double complex arg6,arg7,arg8,arg9,arg10
      real*8 Auv,Buv,Cuv,Duv,Euv
      real*8 Adv,Bdv,Cdv,Ddv,Edv
      real*8 Aub,Bub,Cub,Dub,Eub
      real*8 Adb,Bdb,Cdb,Ddb,Edb
      real*8 fs
      real*8 Nus ,Aus ,Bus ,
     >       Nds ,Ads ,Bds ,
     >       Nss ,Ass ,Bss ,
     >       Nubs,Aubs,Bubs,
     >       Ndbs,Adbs,Bdbs,
     >       Nsbs,Asbs,Bsbs,
     >       ktss,gams
      common /fitp/
     >       Nus ,Aus ,Bus ,
     >       Nds ,Ads ,Bds ,
     >       Nss ,Ass ,Bss ,
     >       Nubs,Aubs,Bubs,
     >       Ndbs,Adbs,Bdbs,
     >       Nsbs,Asbs,Bsbs,
     >       ktss,gams

      Auv =  4.070
      Buv =  0.714
      Cuv =  4.840
      Duv =  0.000
      Euv =  13.40

      Aub =  0.105
      Bub = -0.172
      Cub =  8.060
      Dub =  11.90
      Eub =  0.000

      x = real(z)
      y = aimag(z)

      arg1 = complex(1.+Bus+Cub,0.)
      arg2 = complex(-1.+Aus+Bub+x,y)
      arg3 = complex(Aus+Bub+Bus+Cub+x,y)
      arg4 = complex(Aus+Bub+x,y)
      arg5 = complex(2.+Aus+Bub+Bus+Cub+x,y)
      arg6 = complex(1.+Bus+Cuv,0.)
      arg7 = complex(-1.+Aus+Buv+x,y)
      arg8 = complex(Aus+Bus+Buv+Cuv+x,y)
      arg9 = complex(Aus+Buv+x,y)
      arg10= complex(2.+Aus+Bus+Buv+Cuv+x,y)

      siv_initu = (Aus+Bus)**(Aus+Bus)
     *        /(Aus**Aus)/(Bus**Bus)*Nus*(
     *    Aub*cgamma(arg1)
     *       *cgamma(arg2)/
     *    cgamma(arg3) +
     *    (
     *    Aub*(Eub*(Aus+Bub+z)+Dub*(1.+Aus+Bub+Bus+Cub+z))*
     *    cgamma(arg1)*cgamma(arg4)
     *    )/
     *    cgamma(arg5)+
     *    Auv*cgamma(arg6)*cgamma(arg7)/
     *    cgamma(arg8)+
     *    (
     *    Auv*(Euv*(Aus+Buv+z)+Duv*(1.+Aus+Bus+Buv+Cuv+z))*
     *    cgamma(arg6)*cgamma(arg9)
     *    )/
     *    cgamma(arg10)
     *        )

      return
      end

      function cgamma(z)
      implicit none
      double complex z,cgamma
      external CGAMA
      real*8 x,y,gr,gi
      common /reim/ x,y

      x = real(z)
      y = aimag(z)

      write(1010,*) x

      call CGAMA(x,y,1,GR,GI)

      cgamma = zexp(complex(GR,GI))

      return
      end

      SUBROUTINE CGAMA(X,Y,KF,GR,GI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION A(10)

      PI=3.141592653589793D0

      DATA A/8.333333333333333D-02,-2.777777777777778D-03,
     &       7.936507936507937D-04,-5.952380952380952D-04,
     &       8.417508417508418D-04,-1.917526917526918D-03,
     &       6.410256410256410D-03,-2.955065359477124D-02,
     &       1.796443723688307D-01,-1.39243221690590D+00/


      IF (X.LT.0.0D0) THEN

        X1=X

        Y1=Y

        X=-X

        Y=-Y

      ENDIF

      X0=X

      IF (X.LE.7.0) THEN

        NA=INT(7-X)

        X0=X+NA

      ENDIF

      Z1=DSQRT(X0*X0+Y*Y)

      TH=DATAN(Y/X0)

      GR=(X0-.5D0)*DLOG(Z1)-TH*Y-X0+0.5D0*DLOG(2.0D0*PI)

      GI=TH*(X0-0.5D0)+Y*DLOG(Z1)-Y

      DO K=1,10

        T=Z1**(1-2*K)

        GR=GR+A(K)*T*DCOS((2.0D0*K-1.0D0)*TH)

        GI=GI-A(K)*T*DSIN((2.0D0*K-1.0D0)*TH)
      END DO

      IF (X.LE.7.0) THEN

        GR1=0.0D0

        GI1=0.0D0

        DO J=0,NA-1

          GR1=GR1+.5D0*DLOG((X+J)**2+Y*Y)

          GI1=GI1+DATAN(Y/(X+J))
        END DO

        GR=GR-GR1

        GI=GI-GI1

      ENDIF

      IF (X1.LT.0.0D0) THEN

        Z1=DSQRT(X*X+Y*Y)

        TH1=DATAN(Y/X)

        SR=-DSIN(PI*X)*DCOSH(PI*Y)

        SI=-DCOS(PI*X)*DSINH(PI*Y)

        Z2=DSQRT(SR*SR+SI*SI)

        TH2=DATAN(SI/SR)

        IF (SR.LT.0.0D0) TH2=PI+TH2

        GR=DLOG(PI/(Z1*Z2))-GR

        GI=-TH1-TH2-GI

        X=X1

        Y=Y1

      ENDIF

      RETURN

      END

c------------------------------------------------------------------
c     evolution factor
c------------------------------------------------------------------
      subroutine evo(Q,z,Revo)
      implicit none
      real*8 Q
      double complex z, Revo, pqqsiv
      real*8 as0,asi,asQ
      real*8   alphas,beta0
      external alphas,pqqsiv,beta0
      real*8 c0,c02,pi,CF,CA,mb,alam4,alam5,
     >     bmax,ktw,ptw,Qini,kappa2
      common /cons/ c0,c02,pi,CF,CA,mb,alam4,alam5,
     >     bmax,ktw,ptw,Qini,kappa2
      real*8 Q0,mc

      Q0 = dsqrt(1.9d0)

      mc = 1.275d0

      if (Q.lt.Q0) then
      Q = Q0
      endif

      if (Q.eq.Q0) then
      Revo = 1d0

      elseif (Q.gt.mb) then
      as0 = alphas(Q0)
      asi = alphas(mb)
      asQ = alphas(Q )

      Revo = (asi/as0)**(-pqqsiv(z)/beta0(Q0))*
     *       (asQ/asi)**(-pqqsiv(z)/beta0(Q ))
      else
      as0 = alphas(Q0)
      asQ = alphas(Q )

      Revo = (asQ/as0)**(-pqqsiv(z)/beta0(Q))

      endif

      end


c----------------------------------------------------------------------
      FUNCTION beta0(Q)
      implicit none
      real*8 beta0,Q
      real*8 numflav
      external numflav

      beta0 = 11./2.-numflav(Q)/3.

      end


c----------------------------------------------------------------------
      FUNCTION numflav(Q)
      implicit none
      real*8 numflav
      real*8 Q
      real*8 c0,c02,pi,CF,CA,mb,alam4,alam5,
     >     bmax,ktw,ptw,Qini,kappa2
      common /cons/ c0,c02,pi,CF,CA,mb,alam4,alam5,
     >     bmax,ktw,ptw,Qini,kappa2
      real*8 mc

      mc = 1.275d0

      if (Q.lt.mc) then
      numflav = 3.
      elseif (mc.lt.Q.and.Q.lt.mb) then
      numflav = 4.
      else
      numflav = 5.
      endif

      end

c----------------------------------------------------------------------
      function pqqsiv(z)
      implicit none
      double complex z,pqqsiv,psim
      integer usenc
      common /nckern/ usenc
      real*8 euler,nc
      real*8 c0,c02,pi,CF,CA,mb,alam4,alam5,
     >     bmax,ktw,ptw,Qini,kappa2
      common /cons/ c0,c02,pi,CF,CA,mb,alam4,alam5,
     >     bmax,ktw,ptw,Qini,kappa2
      real*8 Nu ,Au ,Bu ,
     >       Nd ,Ad ,Bd ,
     >       Ns ,As ,Bs ,
     >       Nub,Aub,Bub,
     >       Ndb,Adb,Bdb,
     >       Nsb,Asb,Bsb,
     >       kts,gam
      common /fitp/ Nu ,Au ,Bu ,
     >              Nd ,Ad ,Bd ,
     >              Ns ,As ,Bs ,
     >              Nub,Aub,Bub,
     >              Ndb,Adb,Bdb,
     >              Nsb,Asb,Bsb,
     >              kts,
     >              gam


      euler = 0.577216d0
      nc = ca

      pqqsiv = cf*(3./2.-2.*(euler+psim(z+1.))+1./(z*(z+1.)))
     >         -nc*usenc
      end


c----------------------------------------------------------------------
c      PSI - FUNCTION FOR COMPLEX ARGUMENT
c----------------------------------------------------------------------
       DOUBLE COMPLEX FUNCTION PSIM(Z)
       DOUBLE COMPLEX Z, ZZ, RZ, DZ, SUB
       SUB = DCMPLX (0.D0,0.D0)
       ZZ = Z
  1    CONTINUE
       IF (DREAL (ZZ) .LT. 10.) THEN
         SUB = SUB - 1./ ZZ
         ZZ = ZZ + 1.
         GOTO 1
       END IF
       RZ = 1./ ZZ
       DZ = RZ * RZ
       PSIM = SUB + LOG(ZZ) - RZ/2.- DZ/2520. * ( 210.+ DZ * (-21.+
     1         10.*DZ ))
       RETURN
       END


c--------------------------------------------------------------
c     alphas(q) -- NLO strong coupling constant (from CTEQ)
c--------------------------------------------------------------

      function alphas(q)
      implicit none
      real*8 q,q2,lambda,lambda2,alphas,b0,b1,tt,ktw,ptw
      real*8 c0,c02,pi,CF,CA,mb,alam4,alam5,bmax,kappa2,Qini
      integer nf
      common /cons/ c0,c02,pi,CF,CA,mb,alam4,alam5,
     >     bmax,ktw,ptw,Qini,kappa2

      if (q.le.mb) then
         nf=4
         lambda=alam4
      else
         nf=5
         lambda=alam5
      endif
      b0=11d0-2d0/3d0*nf
      b1=102d0-38d0/3d0*nf

      q2=q*q
      lambda2=lambda*lambda
      tt=dlog(q2/lambda2)
      alphas=4d0*pi/(b0*tt)*(1d0-b1/(b0*b0)*dlog(tt)/tt)
      return
      end
