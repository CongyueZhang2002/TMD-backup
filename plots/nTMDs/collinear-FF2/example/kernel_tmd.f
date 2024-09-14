c--------------------------------------------------------------
c     Test subroutine to check the integration routine qgauss2
c--------------------------------------------------------------
      subroutine integration(res)
      implicit none
      real*8 test,x,x1,x2,qgauss,qgauss2,res
      integer steps
      external test

      steps = 1
      x1    = 0d0
      x2    = 3.1415926535d0
      res   = qgauss2(test,x1,x2,steps)

      return
      end
c--------------------------------------------------------------

c--------------------------------------------------------------
      function test(x)
      implicit none
      real*8 test,x

      test = dsin(10d0*x)*dexp(-x*x/2d0)

      return
      end
c--------------------------------------------------------------






c--------------------------------------------------------------
c     alphas/(4pi) in terms of Lambda_QCD at NNLL
c--------------------------------------------------------------
      function aspi(Q,order)
      implicit none
      real*8 aspi,Q,x,xlog
      integer nf,order
      real*8 b0,b1,b2,b1b0,b2b0,alphas
      real*8 c0,c02,pi,CF,CA,mc,mb,lambda
      common /cons/c0,c02,pi,CF,CA,mc,mb
      integer ord,flag
      real*8 lambda3(3)    !lambda(order) for nf=3 and order=1,2,3
      real*8 lambda4(3)    !lambda(order) for nf=4 and order=1,2,3
      real*8 lambda5(3)    !lambda(order) for nf=5 and order=1,2,3
      common /lambda/lambda3,lambda4,lambda5
      real*8 Beta3(3)    !Beta(order) for nf=3 and order=1,2,3
      real*8 Beta4(3)    !Beta(order) for nf=4 and order=1,2,3
      real*8 Beta5(3)    !Beta(order) for nf=5 and order=1,2,3
      common /Beta/Beta3,Beta4,Beta5

      flag = 2  ! 1 for analytic and 2 for the one of the PDF set

      if(flag.eq.1) then

      ord = order
      ord = 3 !order I set the order always to the highest

      if(Q.lt.mc) then
        lambda=lambda3(ord)
        b0=Beta3(1)
        b1=Beta3(2)
        b2=Beta3(3)
      elseif(Q.ge.mc.and.Q.lt.mb) then
        lambda=lambda4(ord)
        b0=Beta4(1)
        b1=Beta4(2)
        b2=Beta4(3)
      elseif(Q.ge.mb) then
        lambda=lambda5(ord)
        b0=Beta5(1)
        b1=Beta5(2)
        b2=Beta5(3)
      endif

      b1b0 = b1/(b0**2d0)
      b2b0 = b2/(b0*b0*b0)

      x = 2d0*dlog(Q/lambda)
      xlog = dlog(x)

      if(ord.eq.1) then ! LO
        aspi = 1d0/b0/x
      elseif(ord.eq.2) then ! NLO
        aspi = 1d0/b0/x - b1b0/b0*(xlog/x/x)
      elseif(ord.eq.3) then ! NNLO
        aspi = 1d0/b0/x - b1b0/b0*(xlog/x/x)
     >       + b1b0*b1b0/b0*(xlog*xlog-xlog-1d0)/(x*x*x)
     >       + b2b0/b0/(x*x*x)
      elseif(ord.gt.3) then
        write(*,*) 'Order in aspi larger than 3',ord
        stop
      endif


      elseif(flag.eq.2) then
c     calling the subroutine initialized in the main file
c     alphas(mu_r) returns the coupling constant (see alphas.f)
      aspi=alphas(Q)/(4d0*pi)

      endif

      return
      end function aspi



c---------------------------------------------------------
c     aspi_int(mui,muf,power,order) =
c     = \int_{mui}^{muf}dln(mu) aspi(mu,order)^power
c---------------------------------------------------------
      function aspi_int(mui,muf,power,order)
      real*8 aspi_int,mui,muf,aspi_int_aux,qgauss2
      integer power, pow, order, ord, steps
      common /order_aspi/ ord,pow
      external aspi_int_aux

      ord = order
      pow = power
      steps = 4
      aspi_int = qgauss2(aspi_int_aux,mui,muf,steps)

      return
      end function aspi_int
c---------------------------------------------------------
c     Here I define the integrand to be input of qgauss
      function aspi_int_aux(mu)
      real*8 aspi_int_aux,aspi,mu
      integer ord,pow
      common /order_aspi/ ord,pow

      aspi_int_aux = (1d0/mu)
     &               *(aspi(mu,ord))**(float(pow))

      return
      end function aspi_int_aux
c---------------------------------------------------------


c---------------------------------------------------------
c     aspi_log_int(mui,muf,power,order) =
c     = \int_{mui}^{muf}dln(mu) aspi(mu,order)^power ln(muf^2/mu^2)
c---------------------------------------------------------
      function aspi_log_int(mui,muf,Q,power,order)
      real*8 aspi_log_int,mui,muf,aspi_log_int_aux,qgauss2,muff
      real*8 Q,QQ
      integer power, pow, order, ord, steps
      common /order_aspi2/ ord,pow,muff,QQ
      external aspi_log_int_aux

      ord = order
      pow = power
      muff=muf
      QQ=Q
      steps = 4
      aspi_log_int = qgauss2(aspi_log_int_aux,mui,muf,steps)

      return
      end function aspi_log_int
c---------------------------------------------------------
c     Here I define the integrand to be input of qgauss
      function aspi_log_int_aux(mu)
      real*8 aspi_log_int_aux,aspi,mu,muff,QQ
      integer ord,pow
      common /order_aspi2/ ord,pow,muff,QQ

      aspi_log_int_aux = (1d0/mu)
     &                   *(aspi(mu,ord))**(float(pow))
     &                   *2d0*dlog(QQ/mu)

      return
      end function aspi_log_int_aux
c---------------------------------------------------------




c--------------------------------------------------------------
c     alpha_e/(4pi) at LO
c--------------------------------------------------------------
      function aspie(Q)
      implicit none
      real*8 aspie,Q,aspie0,b0
      real*8 c0,c02,pi,CF,CA,mc,mb,mtau
      real*8 Beta0e(4)
      common /Beta0qed/Beta0e
      data aspie0/0.000596532770210d0/ ! aspie0=alpha_e0/(4pi) = (1/133.4)/(4pi)
      data mtau/1.777d0/
      common /cons/c0,c02,pi,CF,CA,mc,mb

      if(Q.lt.mtau) then
        b0=Beta0e(1)
      elseif(Q.ge.mtau.and.Q.lt.mc) then
        b0=Beta0e(2)
      elseif(Q.ge.mc.and.Q.lt.mb) then
        b0=Beta0e(3)
      elseif(Q.ge.mb) then
        b0=Beta0e(4)
      endif


      aspie = aspie0/(1d0+aspie0*b0*2d0*dlog(Q/mtau))

      return
      end function aspie


c---------------------------------------------------------
c     aspi_int_e(mui,muf,power,order) =
c     = \int_{mui}^{muf}dln(mu) aspie(mu,order)^power
c---------------------------------------------------------
      function aspi_int_e(mui,muf,power)
      real*8 aspi_int_e,mui,muf,aspi_int_aux_e,qgauss2
      integer power, pow, steps
      common /order_aspie/ pow
      external aspi_int_aux_e

      pow = power
      steps = 4
      aspi_int_e = qgauss2(aspi_int_aux_e,mui,muf,steps)

      return
      end function aspi_int_e
c---------------------------------------------------------
c     Here I define the integrand to be input of qgauss
      function aspi_int_aux_e(mu)
      real*8 aspi_int_aux_e,aspie,mu
      integer pow
      common /order_aspie/ pow

      aspi_int_aux_e = (1d0/mu)
     &               *(aspie(mu))**(float(pow))

      return
      end function aspi_int_aux_e
c---------------------------------------------------------


c---------------------------------------------------------
c     aspi_log_int_e(mui,muf,power,order) =
c     = \int_{mui}^{muf}dln(mu) aspie(mu,order)^power ln(muf^2/mu^2)
c---------------------------------------------------------
      function aspi_log_int_e(mui,muf,Q,power)
      real*8 aspi_log_int_e,mui,muf,aspi_log_int_aux_e,qgauss2,muff
      real*8 Q,QQ
      integer power, pow, steps
      common /order_aspie2/ muff,QQ,pow
      external aspi_log_int_aux_e

      pow = power
      muff=muf
      QQ=Q
      steps = 4
      aspi_log_int_e = qgauss2(aspi_log_int_aux_e,mui,muf,steps)

      return
      end function aspi_log_int_e
c---------------------------------------------------------
c     Here I define the integrand to be input of qgauss
      function aspi_log_int_aux_e(mu)
      real*8 aspi_log_int_aux_e,aspie,mu,muff,QQ
      integer pow
      common /order_aspie2/ muff,QQ,pow

      aspi_log_int_aux_e = (1d0/mu)
     &                   *(aspie(mu))**(float(pow))
     &                   *2d0*dlog(QQ/mu)

      return
      end function aspi_log_int_aux_e
c---------------------------------------------------------


c---------------------------------------------------------
c     FOR QUARK TMDs
c     \exp[ - \int_{mui}^{muf} dln(mu) ( A \ln(Q^2/mu^2) + B ) ]
c     Inputs: mui,muf
c     The integrand is taken consistently at LL,NLL,NNLL and NNNLL
c     depending on the value of the variable order: 1,2,3,4
c---------------------------------------------------------
      function ExpGamma_q(muii,muff,order)
      real*8 ExpGamma_q,mui,muf,muii,muff,value
      integer order
      real*8 aspi_int, aspi_log_int
      integer i
      real*8 c0,c02
      real*8 pi,CF,CA,mc,mb
      common /cons/c0,c02,pi,CF,CA,mc,mb
      external aspi_int,aspi_log_int
      real*8 CuspQ3(4)      !Cusp quark coefficients for nf=3
      real*8 CuspQ4(4)      !Cusp quark coefficients for nf=4
      real*8 CuspQ5(4)      !Cusp quark coefficients for nf=5
      common /CuspQ/CuspQ3,CuspQ4,CuspQ5
      real*8 NCuspQ3(3)     !Non-cusp quark coefficients for nf=3
      real*8 NCuspQ4(3)     !Non-cusp quark coefficients for nf=4
      real*8 NCuspQ5(3)     !Non-cusp quark coefficients for nf=5
      common /NCuspQ/NCuspQ3,NCuspQ4,NCuspQ5


      if(muff.gt.muii) then
        muf=muff
        mui=muii
      elseif(muff.lt.muii) then
        muf=muii
        mui=muff
      else
        ExpGamma_q=1.0d0
        return
      end if

      if (mui.lt.mc.and.muf.lt.mc) then
        if(order.gt.0) then
          value = CuspQ3(1)*aspi_log_int(mui,muf,muff,1,order)
        endif
        if(order.gt.1) then
          do i=2,order
          value=value
     &          +CuspQ3(i)*aspi_log_int(mui,muf,muff,i,order)
     &          +NCuspQ3(i-1)*aspi_int(mui,muf,i-1,order)
          end do
        endif
      elseif(mui.ge.mc .and. muf.ge.mc
     &       .and. mui.lt.mb .and. muf.lt.mb) then
        if(order.gt.0) then
          value = CuspQ4(1)*aspi_log_int(mui,muf,muff,1,order)
        endif
        if(order.gt.1) then
          do i=2,order
          value=value
     &          +CuspQ4(i)*aspi_log_int(mui,muf,muff,i,order)
     &          +NCuspQ4(i-1)*aspi_int(mui,muf,i-1,order)
          end do
        endif
      elseif(mui.ge.mb .and. muf.ge.mb) then
        if(order.gt.0) then
          value = CuspQ5(1)*aspi_log_int(mui,muf,muff,1,order)
        endif
        if(order.gt.1) then
          do i=2,order
          value=value
     &          +CuspQ5(i)*aspi_log_int(mui,muf,muff,i,order)
     &          +NCuspQ5(i-1)*aspi_int(mui,muf,i-1,order)
          end do
        endif
      elseif(mui.lt.mc .and. muf.ge.mc .and. muf.lt.mb) then
        if(order.gt.0) then
          value = CuspQ3(1)*aspi_log_int(mui,mc,muff,1,order)
     &            +CuspQ4(1)*aspi_log_int(mc,muf,muff,1,order)
        endif
        if(order.gt.1) then
          do i=2,order
          value=value
     &          +CuspQ3(i)*aspi_log_int(mui,mc,muff,i,order)
     &          +CuspQ4(i)*aspi_log_int(mc,muf,muff,i,order)
     &          +NCuspQ3(i-1)*aspi_int(mui,mc,i-1,order)
     &          +NCuspQ4(i-1)*aspi_int(mc,muf,i-1,order)
          end do
        endif
      elseif(mui.lt.mc .and. muf.ge.mb) then
        if(order.gt.0) then
          value = CuspQ3(1)*aspi_log_int(mui,mc,muff,1,order)
     &            +CuspQ4(1)*aspi_log_int(mc,mb,muff,1,order)
     &            +CuspQ5(1)*aspi_log_int(mb,muf,muff,1,order)
        endif
        if(order.gt.1) then
          do i=2,order
          value=value
     &          +CuspQ3(i)*aspi_log_int(mui,mc,muff,i,order)
     &          +CuspQ4(i)*aspi_log_int(mc,mb,muff,i,order)
     &          +CuspQ5(i)*aspi_log_int(mb,muf,muff,i,order)
     &          +NCuspQ3(i-1)*aspi_int(mui,mc,i-1,order)
     &          +NCuspQ4(i-1)*aspi_int(mc,mb,i-1,order)
     &          +NCuspQ5(i-1)*aspi_int(mb,muf,i-1,order)
          end do
        endif
      elseif(mui.ge.mc .and. mui.lt.mb .and. muf.ge.mb) then
        if(order.gt.0) then
          value = CuspQ4(1)*aspi_log_int(mui,mb,muff,1,order)
     &            +CuspQ5(1)*aspi_log_int(mb,muf,muff,1,order)
        endif
        if(order.gt.1) then
          do i=2,order
          value=value
     &          +CuspQ4(i)*aspi_log_int(mui,mb,muff,i,order)
     &          +CuspQ5(i)*aspi_log_int(mb,muf,muff,i,order)
     &          +NCuspQ4(i-1)*aspi_int(mui,mb,i-1,order)
     &          +NCuspQ5(i-1)*aspi_int(mb,muf,i-1,order)
          end do
        endif
      end if

      if(muff.gt.muii) then
        ExpGamma_q=dexp(-1d0*value) ! Because we calculate exp[-gamma] !!!
      elseif(muff.lt.muii) then
        ExpGamma_q=dexp(+1.0d0*value)   ! When limit are inverted
      else
        ExpGamma_q=1.0d0
        return
      end if

      return
      end function ExpGamma_q
c---------------------------------------------------------



c---------------------------------------------------------
c     FOR GLUON TMDs
c     \exp[ - \int_{mui}^{muf} dln(mu) ( A \ln(Q^2/mu^2) + B ) ]
c     Inputs: mui,muf
c     The integrand is taken consistently at LL,NLL,NNLL and NNNLL
c     depending on the value of the variable order: 1,2,3,4
c---------------------------------------------------------
      function ExpGamma_g(muii,muff,order)
      real*8 ExpGamma_g,mui,muf,muii,muff,value
      integer order
      real*8 aspi_int, aspi_log_int
      integer i
      real*8 c0,c02
      real*8 pi,CF,CA,mc,mb
      common /cons/c0,c02,pi,CF,CA,mc,mb
      external aspi_int,aspi_log_int
      real*8 CuspG3(4)      !Cusp quark coefficients for nf=3
      real*8 CuspG4(4)      !Cusp quark coefficients for nf=4
      real*8 CuspG5(4)      !Cusp quark coefficients for nf=5
      common /CuspG/CuspG3,CuspG4,CuspG5
      real*8 NCuspG3(3)     !Non-cusp quark coefficients for nf=3
      real*8 NCuspG4(3)     !Non-cusp quark coefficients for nf=4
      real*8 NCuspG5(3)     !Non-cusp quark coefficients for nf=5
      common /NCuspG/NCuspG3,NCuspG4,NCuspG5



      if(muff.gt.muii) then
        muf=muff
        mui=muii
      elseif(muff.lt.muii) then
        muf=muii
        mui=muff
      else
        ExpGamma_g=1.0d0
        return
      end if

      if (mui.lt.mc.and.muf.lt.mc) then
        if(order.gt.0) then
          value = CuspG3(1)*aspi_log_int(mui,muf,muff,1,order)
        endif
        if(order.gt.1) then
          do i=2,order
          value=value
     &          +CuspG3(i)*aspi_log_int(mui,muf,muff,i,order)
     &          +NCuspG3(i-1)*aspi_int(mui,muf,i-1,order)
          end do
        endif
      elseif(mui.ge.mc .and. muf.ge.mc
     &       .and. mui.lt.mb .and. muf.lt.mb) then
        if(order.gt.0) then
          value = CuspG4(1)*aspi_log_int(mui,muf,muff,1,order)
        endif
        if(order.gt.1) then
          do i=2,order
          value=value
     &          +CuspG4(i)*aspi_log_int(mui,muf,muff,i,order)
     &          +NCuspG4(i-1)*aspi_int(mui,muf,i-1,order)
          end do
        endif
      elseif(mui.ge.mb .and. muf.ge.mb) then
        if(order.gt.0) then
          value = CuspG5(1)*aspi_log_int(mui,muf,muff,1,order)
        endif
        if(order.gt.1) then
          do i=2,order
          value=value
     &          +CuspG5(i)*aspi_log_int(mui,muf,muff,i,order)
     &          +NCuspG5(i-1)*aspi_int(mui,muf,i-1,order)
          end do
        endif
      elseif(mui.lt.mc .and. muf.ge.mc .and. muf.lt.mb) then
        if(order.gt.0) then
          value = CuspG3(1)*aspi_log_int(mui,mc,muff,1,order)
     &            +CuspG4(1)*aspi_log_int(mc,muf,muff,1,order)
        endif
        if(order.gt.1) then
          do i=2,order
          value=value
     &          +CuspG3(i)*aspi_log_int(mui,mc,muff,i,order)
     &          +CuspG4(i)*aspi_log_int(mc,muf,muff,i,order)
     &          +NCuspG3(i-1)*aspi_int(mui,mc,i-1,order)
     &          +NCuspG4(i-1)*aspi_int(mc,muf,i-1,order)
          end do
        endif
      elseif(mui.lt.mc .and. muf.ge.mb) then
        if(order.gt.0) then
          value = CuspG3(1)*aspi_log_int(mui,mc,muff,1,order)
     &            +CuspG4(1)*aspi_log_int(mc,mb,muff,1,order)
     &            +CuspG5(1)*aspi_log_int(mb,muf,muff,1,order)
        endif
        if(order.gt.1) then
          do i=2,order
          value=value
     &          +CuspG3(i)*aspi_log_int(mui,mc,muff,i,order)
     &          +CuspG4(i)*aspi_log_int(mc,mb,muff,i,order)
     &          +CuspG5(i)*aspi_log_int(mb,muf,muff,i,order)
     &          +NCuspG3(i-1)*aspi_int(mui,mc,i-1,order)
     &          +NCuspG4(i-1)*aspi_int(mc,mb,i-1,order)
     &          +NCuspG5(i-1)*aspi_int(mb,muf,i-1,order)
          end do
        endif
      elseif(mui.ge.mc .and. mui.lt.mb .and. muf.ge.mb) then
        if(order.gt.0) then
          value = CuspG4(1)*aspi_log_int(mui,mb,muff,1,order)
     &            +CuspG5(1)*aspi_log_int(mb,muf,muff,1,order)
        endif
        if(order.gt.1) then
          do i=2,order
          value=value
     &          +CuspG4(i)*aspi_log_int(mui,mb,muff,i,order)
     &          +CuspG5(i)*aspi_log_int(mb,muf,muff,i,order)
     &          +NCuspG4(i-1)*aspi_int(mui,mb,i-1,order)
     &          +NCuspG5(i-1)*aspi_int(mb,muf,i-1,order)
          end do
        endif
      end if

      if(muff.gt.muii) then
        ExpGamma_g=dexp(-1d0*value) ! Because we calculate exp[-gamma] !!!
      elseif(muff.lt.muii) then
        ExpGamma_g=dexp(+1.0d0*value)   ! When limit are inverted
      else
        ExpGamma_g=1.0d0
        return
      end if

      end function ExpGamma_g
c---------------------------------------------------------



c---------------------------------------------------------
c     FOR QUARK TMDs in QED
c     \exp[ - \int_{mui}^{muf} dln(mu) ( A_qed \ln(Q^2/mu^2) + B_qed ) ]
c     Inputs: mui,muf
c---------------------------------------------------------
      function ExpGamma_q_e(muii,muff,f)
      real*8 ExpGamma_q_e,mui,muf,muii,muff,value
      real*8 aspi_int_e, aspi_log_int_e
      integer i,f
      real*8 mtau
      real*8 c0,c02,pi,CF,CA,mc,mb
      common /cons/c0,c02,pi,CF,CA,mc,mb
      external aspi_int_e,aspi_log_int_e
      real*8 CuspQ1e(4),CuspQ2e(4)
      common /CuspQqed/CuspQ1e,CuspQ2e
      real*8 NCuspQe
      common /NCuspQqed/NCuspQe
      data mtau/1.777d0/


      if(muff.gt.muii) then
        muf=muff
        mui=muii
      elseif(muff.lt.muii) then
        muf=muii
        mui=muff
      else
        ExpGamma_q_e=1.0d0
        return
      end if


      if(mui.lt.mtau.and.muf.lt.mtau) then
        value = CuspQ1e(1)*aspi_log_int_e(mui,muf,muff,1)
     &        + CuspQ2e(1)*aspi_log_int_e(mui,muf,muff,2)
     &        + NCuspQe*aspi_int_e(mui,muf,1)
      elseif(mui.ge.mtau .and. muf.ge.mtau
     &       .and. mui.lt.mc .and. muf.lt.mc) then
        value = CuspQ1e(2)*aspi_log_int_e(mui,muf,muff,1)
     &        + CuspQ2e(2)*aspi_log_int_e(mui,muf,muff,2)
     &        + NCuspQe*aspi_int_e(mui,muf,1)
      elseif(mui.ge.mc .and. muf.ge.mc
     &       .and. mui.lt.mb .and. muf.lt.mb) then
        value = CuspQ1e(3)*aspi_log_int_e(mui,muf,muff,1)
     &        + CuspQ2e(3)*aspi_log_int_e(mui,muf,muff,2)
     &        + NCuspQe*aspi_int_e(mui,muf,1)
      elseif(mui.ge.mb .and. muf.ge.mb) then
        value = CuspQ1e(4)*aspi_log_int_e(mui,muf,muff,1)
     &        + CuspQ2e(4)*aspi_log_int_e(mui,muf,muff,2)
     &        + NCuspQe*aspi_int_e(mui,muf,1)
      elseif(mui.lt.mtau .and. muf.ge.mtau .and. muf.lt.mc) then
        value = CuspQ1e(1)*aspi_log_int_e(mui,mtau,muff,1)
     &        + CuspQ1e(2)*aspi_log_int_e(mtau,muf,muff,1)
     &        + CuspQ2e(1)*aspi_log_int_e(mui,mtau,muff,2)
     &        + CuspQ2e(2)*aspi_log_int_e(mtau,muf,muff,2)
     &        + NCuspQe*aspi_int_e(mui,muf,1)
      elseif(mui.lt.mtau .and. muf.ge.mc .and. muf.lt.mb) then
        value = CuspQ1e(1)*aspi_log_int_e(mui,mtau,muff,1)
     &        + CuspQ1e(2)*aspi_log_int_e(mtau,mc,muff,1)
     &        + CuspQ1e(3)*aspi_log_int_e(mc,muf,muff,1)
     &        + CuspQ2e(1)*aspi_log_int_e(mui,mtau,muff,2)
     &        + CuspQ2e(2)*aspi_log_int_e(mtau,mc,muff,2)
     &        + CuspQ2e(3)*aspi_log_int_e(mc,muf,muff,2)
     &        + NCuspQe*aspi_int_e(mui,muf,1)
      elseif(mui.lt.mtau .and. muf.ge.mb) then
        value = CuspQ1e(1)*aspi_log_int_e(mui,mtau,muff,1)
     &        + CuspQ1e(2)*aspi_log_int_e(mtau,mc,muff,1)
     &        + CuspQ1e(3)*aspi_log_int_e(mc,mb,muff,1)
     &        + CuspQ1e(4)*aspi_log_int_e(mb,muf,muff,1)
     &        + CuspQ2e(1)*aspi_log_int_e(mui,mtau,muff,2)
     &        + CuspQ2e(2)*aspi_log_int_e(mtau,mc,muff,2)
     &        + CuspQ2e(3)*aspi_log_int_e(mc,mb,muff,2)
     &        + CuspQ2e(4)*aspi_log_int_e(mb,muf,muff,2)
     &        + NCuspQe*aspi_int_e(mui,muf,1)
      elseif(mui.ge.mtau .and. mui.lt.mc
     &       .and. muf.ge.mc .and. muf.lt.mb) then
        value = CuspQ1e(2)*aspi_log_int_e(mui,mc,muff,1)
     &        + CuspQ1e(3)*aspi_log_int_e(mc,muf,muff,1)
     &        + CuspQ2e(2)*aspi_log_int_e(mui,mc,muff,2)
     &        + CuspQ2e(3)*aspi_log_int_e(mc,muf,muff,2)
     &        + NCuspQe*aspi_int_e(mui,muf,1)
      elseif(mui.ge.mtau .and. mui.lt.mc
     &       .and. muf.ge.mb) then
        value = CuspQ1e(2)*aspi_log_int_e(mui,mc,muff,1)
     &        + CuspQ1e(3)*aspi_log_int_e(mc,mb,muff,1)
     &        + CuspQ1e(4)*aspi_log_int_e(mb,muf,muff,1)
     &        + CuspQ2e(2)*aspi_log_int_e(mui,mc,muff,2)
     &        + CuspQ2e(3)*aspi_log_int_e(mc,mb,muff,2)
     &        + CuspQ2e(4)*aspi_log_int_e(mb,muf,muff,2)
     &        + NCuspQe*aspi_int_e(mui,muf,1)
      elseif(mui.ge.mc .and. mui.lt.mb
     &       .and. muf.ge.mb) then
        value = CuspQ1e(3)*aspi_log_int_e(mui,mb,muff,1)
     &        + CuspQ1e(4)*aspi_log_int_e(mb,muf,muff,1)
     &        + CuspQ2e(3)*aspi_log_int_e(mui,mb,muff,2)
     &        + CuspQ2e(4)*aspi_log_int_e(mb,muf,muff,2)
     &        + NCuspQe*aspi_int_e(mui,muf,1)
      endif


ccccc I need to multiply the cusp and the noncusp by the charge^2 of the fermion
ccccc This is true for cuspQ1e, cuspQ2e and NcuspQe, at this particular order
      if(f.eq.2 .or. f.eq.4 .or. f.eq.6
     &   .or. f.eq.-2 .or. f.eq.-4 .or. f.eq.-6) then
        value = value*(2d0/3d0)**2d0
      elseif(f.eq.1 .or. f.eq.3 .or. f.eq.5
     &   .or. f.eq.-1 .or. f.eq.-3 .or. f.eq.-5) then
        value = value*(1d0/3d0)**2d0
      endif


      if(muff.gt.muii) then
        ExpGamma_q_e=dexp(-1d0*value) ! Because we calculate exp[-gamma] !!!
      elseif(muff.lt.muii) then
        ExpGamma_q_e=dexp(+1.0d0*value)   ! When limits are inverted
      else
        ExpGamma_q_e=1.0d0
        return
      end if

      return
      end function ExpGamma_q_e
c---------------------------------------------------------


c--------------------------------------------------------------
c     D function for quarks at fixed order,
c     with order=1,2,3,4 for LL,NLL,NNLL,NNNLL
c--------------------------------------------------------------
      subroutine Dterm_q(bt,mu,order,res)
      implicit none
      real*8 bt,mu,res
      integer order,nf
      real*8 lperp,alf,alf2,alf3,aspi
      external aspi
      real*8 DtermQ11(3:5)
      real*8 DtermQ10(3:5)
      real*8 DtermQ22(3:5)
      real*8 DtermQ21(3:5)
      real*8 DtermQ20(3:5)
      real*8 DtermQ33(3:5)
      real*8 DtermQ32(3:5)
      real*8 DtermQ31(3:5)
      real*8 DtermQ30(3:5)
      common /DtermQ/DtermQ11,DtermQ10,DtermQ22,DtermQ21,DtermQ20,
     &               DtermQ33,DtermQ32,DtermQ31,DtermQ30
      real*8 c0,c02,pi,CF,CA,mc,mb
      common /cons/c0,c02,pi,CF,CA,mc,mb

      if(mu.lt.mc) then
        nf=3
      elseif(mu.ge.mc.and.mu.lt.mb) then
        nf=4
      elseif(mu.ge.mb) then
        nf=5
      endif

      lperp = dlog(mu*mu*bt*bt/c02)
      alf  = aspi(mu,order)
      alf2 = alf*alf
      alf3 = alf2*alf

      res = 0d0
c        write(*,*) order,alf,res

      if(order.ge.2) then
        res = alf*(DtermQ11(nf)*lperp)
      end if

      if(order.ge.3) then
        res = res +
     &            alf2*(DtermQ22(nf)*lperp*lperp
     &                  + DtermQ21(nf)*lperp
     &                  + DtermQ20(nf))
      end if


      if(order.eq.4) then
        res = res +
     &            alf3*(DtermQ33(nf)*lperp*lperp*lperp
     &                  + DtermQ32(nf)*lperp*lperp
     &                  + DtermQ31(nf)*lperp + DtermQ30(nf))
      end if
c        write(*,*) 'D32=',nf,DtermQ32(nf)

      return
      end subroutine
c--------------------------------------------------------------


c--------------------------------------------------------------
c     D function for gluons at fixed order,
c     with order=1,2,3,4 for LL,NLL,NNLL,NNNLL
c--------------------------------------------------------------
      subroutine Dterm_g(bt,mu,order,res)
      implicit none
      real*8 bt,mu,res
      integer order,nf
      real*8 lperp,alf,alf2,alf3,aspi
      external aspi
      real*8 DtermG11(3:5)
      real*8 DtermG10(3:5)
      real*8 DtermG22(3:5)
      real*8 DtermG21(3:5)
      real*8 DtermG20(3:5)
      real*8 DtermG33(3:5)
      real*8 DtermG32(3:5)
      real*8 DtermG31(3:5)
      real*8 DtermG30(3:5)
      common /DtermG/DtermG11,DtermG10,DtermG22,DtermG21,DtermG20,
     &               DtermG33,DtermG32,DtermG31,DtermG30
      real*8 c0,c02
      real*8 pi,CF,CA,mc,mb
      common /cons/c0,c02,pi,CF,CA,mc,mb

      if(mu.lt.mc) then
        nf=3
      elseif(mu.ge.mc.and.mu.lt.mb) then
        nf=4
      elseif(mu.ge.mb) then
        nf=5
      endif

      lperp = dlog(mu*mu*bt*bt/c02)
      alf  = aspi(mu,order)
      alf2 = alf*alf
      alf3 = alf2*alf

      res = 0d0

      if(order.ge.2) then
        res = alf*(DtermG11(nf)*lperp)
      end if

      if(order.ge.3) then
      res = res +
     &      alf2*(DtermG22(nf)*lperp*lperp
     &            + DtermG21(nf)*lperp
     &            + DtermG20(nf))
      end if

      if(order.ge.4) then
      res = res +
     &      alf3*(DtermG33(nf)*lperp*lperp*lperp
     &            + DtermG32(nf)*lperp*lperp
     &            + DtermG31(nf)*lperp + DtermG30(nf))
      end if

      end subroutine
c--------------------------------------------------------------



c--------------------------------------------------------------
c     D function for quarks at fixed order in QED
c     with order=1,2,3,4 for LL,NLL,NNLL,NNNLL
c--------------------------------------------------------------
      subroutine Dterm_q_e(bt,mu,res,f)
      implicit none
      real*8 bt,mu,res
      integer order,nf,f
      real*8 lperp,alfe,aspie
      external aspie
      real*8 c0,c02,pi,CF,CA,mc,mb
      common /cons/c0,c02,pi,CF,CA,mc,mb
      real*8 mtau
      data mtau/1.777d0/

      lperp = dlog(mu*mu*bt*bt/c02)
      alfe  = aspie(mu)

        res = alfe*2d0*lperp

ccccc I need to multiply it by the charge^2 of the fermion, since I replace
ccccc CF by Qi^2
      if(f.eq.2 .or. f.eq.4 .or. f.eq.6
     &   .or. f.eq.-2 .or. f.eq.-4 .or. f.eq.-6) then
        res = res*(2d0/3d0)**2d0
      elseif(f.eq.1 .or. f.eq.3 .or. f.eq.5
     &   .or. f.eq.-1 .or. f.eq.-3 .or. f.eq.-5) then
        res = res*(1d0/3d0)**2d0
      endif



      return
      end subroutine
c--------------------------------------------------------------




c---------------------------------------------------------
c     KERNEL FOR QUARK TMDs
c     The kernel is taken consistently at LL,NLL,NNLL and NNNLL
c     depending on the value of the variable order: 1,2,3,4
c---------------------------------------------------------
      subroutine kernel_q(b,mui,muf,Qi,Qf,order,res)
      real*8 b,mui,muf,Qi,Qf,res
      integer order
      real*8 ExpGamma_q,res1,res2

      res1 = ExpGamma_q(mui,muf,order)
      call Dterm_q(b,mui,order,res2)
      !print *, res1,res2

      res = res1*dexp(-res2*2d0*dlog(Qf/Qi))

      return
      end subroutine
c---------------------------------------------------------




c---------------------------------------------------------
c     KERNEL FOR GLUON TMDs
c     The kernel is taken consistently at LL,NLL,NNLL and NNNLL
c     depending on the value of the variable order: 1,2,3,4
c---------------------------------------------------------
      subroutine kernel_g(b,mui,muf,Qi,Qf,order,res)
      real*8 b,mui,muf,Qi,Qf,res
      integer order
      real*8 ExpGamma_g,res1,res2

      res1 = ExpGamma_g(mui,muf,order)
      call Dterm_g(b,mui,order,res2)

      res = res1*dexp(-res2*2d0*dlog(Qf/Qi))

      return
      end subroutine
c---------------------------------------------------------







c--------------------------------------------------------------------------
c               / xf
c               |
c       value = | dx func(x)         n > 1, number of divisions for xf-xi
c               |
c              / xi
c--------------------------------------------------------------------------

      function qgauss2(func,xi,xf,n)
      implicit none
      real*8 qgauss2,func,xi,xf,xii,xff,xn,value,x1,x2,val
      integer i,n
      external func

      if(xf.gt.xi) then
        xff = xf
        xii = xi
      elseif(xf.lt.xi) then
        xff = xi
        xii = xf
      endif

      if(n.le.1) then           ! same as n=1
         x1=xii
         x2=xff
         call gauss2(func,x1,x2,val)
         qgauss2=val
         return
      endif

      xn=(xff-xii)/float(n)
      value=0.d0
      x2=xii
      Do i=1,n
         x1=x2
         x2=x1+xn
         call gauss2(func,x1,x2,val)
         value=value+val
      end do

      qgauss2=value

      if(xff.gt.xii) then
        qgauss2 = qgauss2
      elseif(xff.lt.xii) then
        qgauss2 = -1d0*qgauss2
      endif


      return
      end

c-------------------------------------------------
      subroutine gauss2(func,xi,xf,value)
      implicit none
      real*8 func,xi,xf,xii,xff,value,xm,xr,dx
      real*8 eps
      integer j,flag
      data eps /1.0d-25/
      real*8 Xi8(8),Wi8(8)     ! Points and weights for 8 points
      common /gaussparam8/Xi8,Wi8
      real*8 Xi12(12),Wi12(12)     ! Points and weights for 12 points
      common /gaussparam12/Xi12,Wi12
      real*8 Xi16(16),Wi16(16)     ! Points and weights for 16 points
      common /gaussparam16/Xi16,Wi16

      if(xf.gt.xi) then
        xff = xf
        xii = xi
      elseif(xf.lt.xi) then
        xff=xi
        xii=xf
      endif

      xm=0.5d0*(xff+xii)
      xr=0.5d0*(xff-xii)
      if (abs(xr).lt.eps) print *,
     $     'WARNING: Too high accuracy required for GAUSS2!'

      value=0.d0

      flag = 16  ! flag=8/12/16 for 8/12/16 integration points
      if(flag.eq.8) then ! Integration with 8 points
        do j=1,8
         dx=xr*Xi8(j)
         value=value + Wi8(j)*func(xm+dx)
        enddo
      elseif(flag.eq.12) then ! Integration with 12 points
        do j=1,12
         dx=xr*Xi12(j)
         value=value + Wi12(j)*func(xm+dx)
        enddo
      elseif(flag.eq.16) then ! Integration with 16 points
        do j=1,16
         dx=xr*Xi16(j)
         value=value + Wi16(j)*func(xm+dx)
        enddo
      endif

      value=xr*value

      if(xff.gt.xii) then
        value = value
      elseif(xff.lt.xii) then
        value = -1d0*value
      endif


      return
      end


