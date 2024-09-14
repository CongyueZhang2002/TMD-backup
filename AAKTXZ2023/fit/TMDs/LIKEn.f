*************************************************************************
*
*  LIKEn21: COLLINEAR NUCLEAR FRAGMENTATION FUNCTIONS
*
*  P. Zurita
*
*  arXiv:2101.01088 [hep-ph]
*
*     CALL LIKEn (CL,IH,IC,SET,Z, Q2, A, U, UB, D, DB, S, SB, C, CB, B, BB, GL)
*
*  INPUT:
*  CL = CONFIDENCE LEVEL    1: 68%, DeltaChi^2=11
*                           2: 90%, DeltaChi^2=32
*
*  IH = hadron type    1: PION      ! only one for now
*                      2: KAON
*                      3: PROTON
*
*  IC = Hadron Charge  0: 0 (as average of + and -)
*                      1: +
*                     -1: -
*
*  SET = central (=0) or uncertainty set (1 to 14, 1-7 +, 8-14 -)
*
*
*  Z                    (between  0.01   and  0.9999)
*  Q2 = scale in GeV**2 (between  1.0    and  1.D5)
*             (for values outside the allowed range the program
*              writes a warning and extrapolates to the z and
*              Q2 values requested)
*
*  A = nucleus required. Right now only A = 1, 4, 12, 20, 56, 84, 131, 197, 208
*
*
*   OUTPUT: U, UB, D, DB, S, SB,   C, CB,    B,   BB,      GL
*           U Ubar D Dbar S Sbar Charm=Cbar Bottom=Bbar Gluon
*           Always Z times the distribution is returned
*
*   COMMON:  The main program or the calling routine has to have
*            a common block  COMMON / INITIATE / INIT , and  INIT
*            has always to be zero when LIKEn is called for the
*            first time or when the SET (central/eigenvector) has been changed.
*
*   Doubts, questions, comments to Pia Zurita:
*
*   pia@df.uba.ar
*   mzurita@bnl.gov
*   maria.zurita@ur.de
*   pia.zurita@usc.es
*
********************************************************************

      SUBROUTINE LIKEn(CL,IH,IC,set,Z,Q2,A,U,UB,D,DB,S,SB,C,CB,B,BB,GL)

       Implicit none
       INTEGER INIT,NX,I,IS,A,NQ,NARG,NPART,set,neig,CL,IH,IC,iiread,iq
      PARAMETER (NPART=11, NX=47, NQ=24, NARG=2, neig=7)
      double precision XUTOTF(NX,NQ), XDTOTF(NX,NQ), XSTOTF(NX,NQ)
      double precision XCTOTF(NX,NQ), XBTOTF(NX,NQ), XGF(NX,NQ)
      double precision XUVALF(NX,NQ), XDVALF(NX,NQ), XSVALF(NX,NQ)
      double precision XCVALF(NX,NQ), XBVALF(NX,NQ)

      double precision PARTON (NPART,NQ,NX-1)!, PROTON (NPART,NQ,NX)
      double precision QS(NQ), ZZ(NX), ZT(NARG), ARRF(NX+NQ)
      double precision Z0, Z1, Z, Q2,dum
      double precision U,UB,D,DB,S,SB,C,CB,B,BB,GL!,RU,RUB,RD,RDB,RS,RSB,RC,RB
      double precision BBP,BP,BTOT,CBP,CP,CTOT,DBP,DP,DTOT,DVAL,SVAL
      double precision UBP,UP,UTOT,SBP,SP,UVAL,CVAL,BVAL,STOT
       integer IX,N,M, NA(NARG),Anuc(9)

      COMMON / INITIATE / INIT
      SAVE XUTOTF,XDTOTF,XSTOTF,XCTOTF,XBTOTF,XGF,XUVALF,XDVALF,XSVALF,
     1     XCVALF,XBVALF
      SAVE ZZ, NA, ARRF
      double precision dummy
      character*128 fname,path
      character*7 charCL,charA,chariset
      character*20 charset,charspecies
      double precision FINTLIKE
      common/gridpath/path

*...Q**2 VALUES OF THE GRID.
       DATA QS / 1.d0, 1.25D0, 1.5D0, 2.5D0,
     1           4.0D0, 6.4D0, 1.0D1, 1.5D1, 2.5D1, 4.0D1, 6.4D1,
     2           1.0D2, 1.8D2, 3.2D2, 5.8D2, 1.0D3, 1.8D3,
     3           3.2D3, 5.8D3, 1.0D4, 1.8D4, 3.2D4, 5.8D4, 1.0D5/
*...Z VALUES OF THE GRID.

       DATA ZZ/0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
     4        0.095, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
     5        0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475,  0.5,
     6        0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7,
     7        0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9,
     8        0.925, 0.95, 0.975, 1.0/

*...A values available

       DATA Anuc/1, 4, 12, 20, 56, 84, 131, 197, 208/


C-----specify path here:
      path = '../fit/TMDs/grids/7param/'


*...CHECK OF Z, Q2 and A values:
C      IF ( (Z.LT.0.01D0) .OR. (Z.GT.1.0D0) ) THEN
C          WRITE(6,91)
C 91       FORMAT (2X,'PARTON INTERPOLATION: Z OUT OF RANGE. STOP')
C         STOP
C      ENDIF
C      IF ( (Q2.LT.1.D0) .OR. (Q2.GT.1.D5) ) THEN
C          WRITE(6,92)
C 92       FORMAT (2X,'PARTON INTERPOLATION: Q2 OUT OF RANGE. STOP')
C          STOP
C      ENDIF

       IF ((A.ne.1).and.(A.NE.4).and.(A.NE.12).and.(A.NE.20).and.
     1 (A.NE.56).and.(A.NE.84).and.(A.NE.131).and.(A.ne.197).and.
     2 (A.ne.208)) THEN
           WRITE(6,93)
  93   FORMAT (2X,
     1 'Nucleus not available. Ask P. Zurita to run it for you. STOP')
           STOP
       ENDIF

*...INITIALIZATION :

       IF (INIT.NE.0) GOTO 16

       if (A.eq.1) then
       charA='vacuum_'
       else
       write(charA, '(i0)') A
       charA=trim(charA)//'_'
       endif

        if (IH.eq.1) then
       charspecies='pion/'
       elseif (IH.eq.2) then
       charspecies='kaon/'
       else
           WRITE(6,99)
  99       FORMAT (2X,'Hadronic species not available. Stop.')
           STOP
       endif

       if (set.eq.0) then
       charset='central'
       elseif ((set.gt.0).and.(set.le.neig)) then
       write(chariset, '(i0)') set
       charset='p0'//trim(chariset)
       elseif (set.gt.neig) then
       write(chariset, '(i0)') set-neig
       charset='m0'//trim(chariset)
       endif

       if (set.ne.0) then   ! if an eigenvector is requested
       if (CL.eq.1) then
       charcl='_68CL'
       elseif (CL.eq.2) then
       charcl='_90CL'
       else
           WRITE(6,97)
  97       FORMAT (2X,'Confidence level not available. Stop.')
           STOP
       endif
       endif

       fname=trim(path)//trim(charspecies)//trim(charA)//trim(charset)

       if (set.eq.0) then
       fname=trim(fname)//'.txt'
       else
       fname=trim(fname)//trim(charcl)//'.txt'
       endif

c now read the files

       IIREAD=12
       OPEN(IIREAD,FILE=fname)

       DO N = 1, NQ
       DO M = 1, NX-1
       READ(IIREAD,90)
     1 PARTON(1,N,M), PARTON(2,N,M), PARTON(3,N,M), PARTON(4,N,M),
     1 PARTON(5,N,M), PARTON(6,N,M), PARTON(7,N,M), PARTON(8,N,M),
     2 PARTON(9,N,M), PARTON(10,N,M), PARTON(11,N,M)
  90   FORMAT (11(2x,G12.6))
       enddo
       enddo
       CLOSE(IIREAD)
C
      INIT = 1
*....ARRAYS FOR THE INTERPOLATION SUBROUTINE :
      DO IQ = 1, NQ
      DO IX = 1, NX-1
        Z0 = ZZ(IX)
        Z1 = 1.D0-ZZ(IX)
        XGF(IX,IQ)    = PARTON(1,IQ,IX) / (Z1**4 * Z0**0.3)
        XUTOTF(IX,IQ) = PARTON(2,IQ,IX) / (Z1**4 * Z0**0.5)
        XDTOTF(IX,IQ) = PARTON(3,IQ,IX) / (Z1**4 * Z0**0.5)
        XSTOTF(IX,IQ) = PARTON(4,IQ,IX) / (Z1**4 * Z0**0.5)
        XCTOTF(IX,IQ) = PARTON(5,IQ,IX) / (Z1**7 * Z0**0.3)
        XBTOTF(IX,IQ) = PARTON(6,IQ,IX) / (Z1**7 * Z0**0.3)
        XUVALF(IX,IQ) = PARTON(7,IQ,IX) / (Z1**4 * Z0**0.5)   ! this is u-ubar
        XDVALF(IX,IQ) = PARTON(8,IQ,IX) / (Z1**4 * Z0**0.5)   ! this is d-dbar
        XSVALF(IX,IQ) = PARTON(9,IQ,IX) / (Z1**4 * Z0**0.5)   ! this is s-sbar
        XCVALF(IX,IQ) = PARTON(10,IQ,IX) / (Z1**7 * Z0**0.3)   ! this is c-cbar
        XBVALF(IX,IQ) = PARTON(11,IQ,IX) / (Z1**7 * Z0**0.3)   ! this is b-bbar
      enddo
        XGF(NX,IQ)    = 0.D0
        XUTOTF(NX,IQ) = 0.D0
        XDTOTF(NX,IQ) = 0.D0
        XSTOTF(NX,IQ) = 0.D0
        XCTOTF(NX,IQ) = 0.D0
        XBTOTF(NX,IQ) = 0.D0
        XUVALF(NX,IQ) = 0.D0
        XDVALF(NX,IQ) = 0.D0
        XSVALF(NX,IQ) = 0.D0
        XCVALF(NX,IQ) = 0.D0
        XBVALF(NX,IQ) = 0.D0
       enddo


      NA(1) = NX
      NA(2) = NQ
      DO IX = 1, NX
        ARRF(IX) = DLOG(ZZ(IX))
      enddo
      DO IQ = 1, NQ
        ARRF(NX+IQ) = DLOG(QS(IQ))
      enddo
  16  CONTINUE
*...INTERPOLATION :
      ZT(1) = DLOG(Z)
      ZT(2) = DLOG(Q2)

      GL   = FINTLIKE(NARG,ZT,NA,ARRF,XGF)    * (1.D0-Z)**4 * Z**0.3
      UTOT = FINTLIKE(NARG,ZT,NA,ARRF,XUTOTF) * (1.D0-Z)**4 * Z**0.5
      DTOT = FINTLIKE(NARG,ZT,NA,ARRF,XDTOTF) * (1.D0-Z)**4 * Z**0.5
      STOT = FINTLIKE(NARG,ZT,NA,ARRF,XSTOTF) * (1.D0-Z)**4 * Z**0.5
      CTOT = FINTLIKE(NARG,ZT,NA,ARRF,XCTOTF) * (1.D0-Z)**7 * Z**0.3
      BTOT = FINTLIKE(NARG,ZT,NA,ARRF,XBTOTF) * (1.D0-Z)**7 * Z**0.3
      UVAL = FINTLIKE(NARG,ZT,NA,ARRF,XUVALF) * (1.D0-Z)**4 * Z**0.5
      DVAL = FINTLIKE(NARG,ZT,NA,ARRF,XDVALF) * (1.D0-Z)**4 * Z**0.5
      SVAL = FINTLIKE(NARG,ZT,NA,ARRF,XSVALF) * (1.D0-Z)**4 * Z**0.5
      CVAL = FINTLIKE(NARG,ZT,NA,ARRF,XCVALF) * (1.D0-Z)**7 * Z**0.3
      BVAL = FINTLIKE(NARG,ZT,NA,ARRF,XBVALF) * (1.D0-Z)**7 * Z**0.3

       Up  = (UTOT+UVAL)/2.
       UBp = (UTOT-UVAL)/2.
       Dp  = (DTOT+DVAL)/2.
       DBp = (DTOT-DVAL)/2.
       Sp  = (STOT+SVAL)/2.
       SBp = (STOT-SVAL)/2.
       Cp  = (CTOT+CVAL)/2.
       CBp = (CTOT-CVAL)/2.
       Bp  = (BTOT+BVAL)/2.
       BBp = (BTOT-BVAL)/2.

       IF (IC.EQ.1) THEN
       U  = Up
       UB = UBp
       D  = Dp
       DB = DBp
       S  = Sp
       SB = SBp
       C  = Cp
       CB = CBp
       B  = Bp
       BB = BBp
       GL=GL
       ELSEIF (IC.EQ.-1) THEN
       U  = UBp
       UB = Up
       D  = DBp
       DB = Dp
       S  = SBp
       SB = Sp
       C  = CBp
       CB = Cp
       B  = Bp
       BB = BBp
       GL=GL
       ELSEIF (IC.EQ.0) THEN
       U  = (UBp+Up)/2.
       UB =  U
       D  = (DBp+Dp)/2.
       DB =  D
       S  = (SBp+Sp)/2.
       SB =  S
       C  = (CBp+Cp)/2.
       CB = C
       B  = (BBp+Bp)/2.
       BB = B
       GL=GL
       ELSE
         WRITE(6,94)
 94      FORMAT (2X,' WRONG CHARGE')
         STOP
       END IF


 60   RETURN
       END
*
*...CERN LIBRARY ROUTINE E104 (INTERPOLATION) :
*
      FUNCTION FINTLIKE(NARG,ARG,NENT,ENT,TABLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CMS      DIMENSION ARG(5),NENT(5),ENT(63),TABLE(882)
      DIMENSION ARG(2),NENT(2),ENT(71),TABLE(1128)
      DIMENSION D(5),NCOMB(5),IENT(5)
      KD=1
      M=1
      JA=1
      DO I=1,NARG
      NCOMB(I)=1
      JB=JA-1+NENT(I)
         DO 2 J=JA,JB
      IF (ARG(I).LE.ENT(J)) GO TO 3
    2 CONTINUE
      J=JB
    3 IF (J.NE.JA) GO TO 4
      J=J+1
    4 JR=J-1
      D(I)=(ENT(J)-ARG(I))/(ENT(J)-ENT(JR))
      IENT(I)=J-JA
      KD=KD+IENT(I)*M
      M=M*NENT(I)
      JA=JB+1
      enddo
      FINTLIKE=0.D0
   10 FAC=1.D0
      IADR=KD
      IFADR=1
         DO I=1,NARG
      IF (NCOMB(I).EQ.0) GO TO 12
      FAC=FAC*(1.D0-D(I))
      GO TO 15
   12 FAC=FAC*D(I)
      IADR=IADR-IFADR
   15 IFADR=IFADR*NENT(I)
      enddo
      FINTLIKE=FINTLIKE+FAC*TABLE(IADR)
      IL=NARG
   40 IF (NCOMB(IL).EQ.0) GO TO 80
      NCOMB(IL)=0
      IF (IL.EQ.NARG) GO TO 10
      IL=IL+1
         DO K=IL,NARG
      NCOMB(K)=1
      enddo
      GO TO 10
   80 IL=IL-1
      IF(IL.NE.0) GO TO 40
      RETURN
      END
