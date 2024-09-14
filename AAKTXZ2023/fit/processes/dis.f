c--------------------------------------------------------------
c     DIS X-section in integrated in xb, Q2
c--------------------------------------------------------------
      subroutine DIS_OVERNUQ2(rtss, QQ0,
     >               QQ2low, QQ2high, NNulow, NNuhigh,
     >               fuu)
      implicit none
      real*8 rtss,QQ0,QQ2low, QQ2high, NNulow, NNuhigh
      real*8 fuu
      real*8 rts , Q0,                  Nulow, Nuhigh
      common /Nuintsid/ rts, Q0, Nulow, Nuhigh
      real*8 qgauss, fun_DIS_Q2
      external qgauss, fun_DIS_Q2

      rts= rtss
      Q0 = QQ0
      Nulow = NNulow
      Nuhigh = NNuhigh

      fuu = qgauss(fun_DIS_Q2,QQ2low,QQ2high,3)

      return
      end

c-----

      function fun_DIS_Q2(Q2)
      implicit none
      real*8 Q2
      real*8 fun_DIS_Q2
      real*8 rts, Q0, Nulow, Nuhigh
      common /Nuintsid/ rts, Q0, Nulow, Nuhigh
      real*8 fuu

      call DIS_OVERNU(rts, Q0, Q2, Nulow, Nuhigh, fuu)
      fun_DIS_Q2 = fuu

      return
      end

c--------------------------------------------------------------
c     DIS X-section in {Q2}: d[sigma]/dQ^2
c     integrated in xb
c--------------------------------------------------------------
      SUBROUTINE DIS_OVERNU(rrts, QQ0, QQs, Nulow, Nuhigh, fuu)
      IMPLICIT NONE
      real*8 rrts, QQ0, QQs, Nulow, Nuhigh
      real*8 fuu
      real*8 Q0, Qs, rts
      common /DIS_KIN/ Q0, Qs, rts
      real*8 M_proton,xblow,xbhigh
      real*8   qgauss,fun_DIS_xb
      external qgauss,fun_DIS_xb

      Q0 = QQ0
      Qs = QQs
      rts = rrts

      M_proton = 0.938272046d0
      xblow  = 0.5d0*Qs/(M_proton*Nuhigh)
      xbhigh = 0.5d0*Qs/(M_proton*Nulow)
      fuu = qgauss(fun_DIS_xb,xblow,xbhigh,3)

      return
      END SUBROUTINE

C-----


      function fun_DIS_xb(xb)
      implicit none
      real*8 xb
      real*8 fun_DIS_xb
      real*8 Q0, Qs, rts
      real*8 F2light_A, FLlight_A, F2light, FLlight
      real*8 dis
      real*8 xmax

      common /DIS_KIN/ Q0, Qs, rts
      xmax = 0.98d0
      if(xb.gt.xmax) then
        xb = xmax
      endif

      CALL ComputeStructureFunctionsAPFEL(Q0,dsqrt(Qs))
      F2light_A = F2light(xb)
      FLlight_A = FLlight(xb)
      CALL dis_cx(rts, xb, Qs, F2light_A,FLlight_A, dis)
      fun_DIS_xb = dis


      return
      end


C---------------------------------------------------------
C     Gives dsigma/dxB/dQ^2
C---------------------------------------------------------
      SUBROUTINE dis_cx(rts, xb, Q2, F2, FL, dis)
      IMPLICIT NONE
      real*8 rts,xb, Q2, F2, FL
      real*8 dis
      double precision alpha_EM, PI
      double precision KQ, Yplus
      real*8 y

C-----------Parameters
      alpha_EM = 7.297351d-3
      PI = 4.d0*DATAN(1.D0)

C-----------Calculation of dsigmadxbdQ2
      y = Q2/(rts*xb)
      Yplus = 1d0+(1d0-y)**2d0
      KQ = 2*PI*(alpha_EM**2)/(xb*Q2**2)
      dis = KQ*(Yplus*F2-FL*(y**2))
      return
      END SUBROUTINE
