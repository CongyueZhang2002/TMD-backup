      program master
      implicit none
      real*8 x,Q
      real*8 U1,D,UB,DB,SS,SB,GL
      real*8 U2
      integer AA, nloops
      integer set,k0
      integer i

      nloops =1
      AA = 197


      call SetLHAPARM('SILENT') ! To not show the calls, although they are called


C-----INITIALIZE PDF SETS USING LHAPDF
      print *, "1"
      CALL InitPDFsetByNameM(1,"CT14nlo")
      CALL InitPDFM(1,0)
      print *, "2"
      CALL InitPDFsetByNameM(2,"nCTEQ15_197") ! Au
      CALL InitPDFM(2,0)
      print *, "3"
      CALL InitPDFsetByNameM(3,"EPPS16nlo_CT14nlo_Au197") ! Au
      CALL InitPDFM(3,0)


      print*, "x,U_CTEQ,U_EPPS,R = U_CTEQ/U_EPPS"

      Q = sqrt(2.4d0)
      x = 0.004d0
      ! nCTEQ
      set = 2
      call PDF_LHAPDF(x,Q,U1,D,UB,DB,SS,SB,GL,set)
      ! EPPS
      set = 3
      call EPPS16_proton(x,Q,AA,U2,D,UB,DB,SS,SB,GL,nloops)
      print*, x,U1,U2,U1/U2

      x = 0.04d0
      ! nCTEQ
      set = 2
      call PDF_LHAPDF(x,Q,U1,D,UB,DB,SS,SB,GL,set)
      ! EPPS
      set = 3
      call EPPS16_proton(x,Q,AA,U2,D,UB,DB,SS,SB,GL,nloops)
      print*, x,U1,U2,U1/U2

      x = 0.4d0
      ! nCTEQ
      set = 2
      call PDF_LHAPDF(x,Q,U1,D,UB,DB,SS,SB,GL,set)
      ! EPPS
      set = 3
      call EPPS16_proton(x,Q,AA,U2,D,UB,DB,SS,SB,GL,nloops)
      print*, x,U1,U2,U1/U2


      k0 = 20

      WRITE(k0, *) 'x Q U_CTEQ U_EPPS R'

      DO i = 1, 100
C-----PT RANGE:
      x = 0.004 ! from pt = 0.01 to pt = 1.0 GeV
C-----Set z-value for nTMDFFs:
      x = 0.004+ 0.004*(i-1)

      set = 2
      call PDF_LHAPDF(x,Q,U1,D,UB,DB,SS,SB,GL,set)
      set = 3
      call EPPS16_proton(x,Q,AA,U2,D,UB,DB,SS,SB,GL,nloops)

      WRITE(k0,*) x,Q,U1,U2,U1/U2
      ENDDO
      CLOSE(k0)


      return
      end program
C

C----------------------------------------------------------------------------------
c     PDF in proton: nloops = 1 is collinear part ONLY, nloops = 2 has convolutions
c----------------------------------------------------------------------------------

            subroutine PDF_LHAPDF(x,Q,U,D,UB,DB,SS,SB,GL,set)
            implicit none
            real*8 x,Q,U,D,UB,DB,SS,SB,GL, pdf(-6:6)
            real*8 Un,Dn,UBn,DBn,SSn,SBn,GLn
            real*8 as,alphas,pref,pi
            integer set
            data pi /3.14159265359d0/
            integer nloops
            integer proc
            common /process/ proc
            real*8 hQ
            common /hardQ/ hQ

            if (Q.le.1d0) then
                Q=1d0
            endif

            pdf(-6) =	0d0
            pdf(-5) = 0d0
            pdf(-4) = 0d0
            pdf(-3) =	0d0
            pdf(-2) = 0d0
            pdf(-1) =	0d0
            pdf( 0) = 0d0
            pdf( 1) =	0d0
            pdf( 2) = 0d0
            pdf( 3) = 0d0
            pdf( 4) = 0d0
            pdf( 5) = 0d0
            pdf( 6) = 0d0

            call evolvePDFM(set,x,Q,pdf)

            U=pdf(2)/x
            D=pdf(1)/x
            SS=pdf(3)/x
            UB=pdf(-2)/x
            DB=pdf(-1)/x
            SB=pdf(-3)/x
            GL=pdf(0)/x

            end subroutine

C
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

      pset =1
      order = 1



      Qmin = 1.3d0

      IF(Q.lt.Qmin) THEN
      Q = Qmin
      ENDIF

      CALL EPPS16(order,pset,AA,x,Q,ruv,rdv,ru,rd,rs,rc,rb,rg)
      call PDF_LHAPDF(x,Q,fU,fD,fUB,fDB,fSS,fSB,fGL,1)
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


      return
      end
