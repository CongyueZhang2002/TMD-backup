      program master
      implicit none 
      integer NX, NB
      parameter (NX = 104,NB=51)
      real*8 Q, BB(NB), XB(NX), zh
      integer K1,K3
      integer A,Z
      real*8 b
      real*8 u,d,s,ub,db,sb,g
      integer nloops, TIC, TIH
      COMMON / ttarget / TIH,TIC
      integer IIREAD

      nloops = 2

      call SetLHAPARM('SILENT') ! To not show the calls, although they are called
      call InitPDFsetByName("CT14nlo")
      call InitPDF(0)

       DATA BB /1.00d-4,1.26d-4,1.56d-4,2.02d-4,2.56d-4,3.24d-4,4.09d-4,
     $          5.18d-4,6.55d-4,8.29d-4,1.05d-3,1.33d-3,1.68d-3,2.12d-3,
     $          2.68d-3,3.39d-3,4.29d-3,5.43d-3,6.87d-3,8.69d-3,1.10d-2,
     $          1.39d-2,1.76d-2,2.22d-2,2.81d-2,3.56d-2,4.50d-2,5.69d-2,
     $          7.20d-2,9.10d-2,1.15d-1,1.46d-1,1.84d-1,2.33d-1,2.95d-1,
     $          3.73d-1,4.71d-1,5.96d-1,7.54d-1,9.54d-1,1.21d0 ,1.53d0 ,
     $   1.73d0,1.93d0 ,2.44d0 ,3.09d0 ,3.91d0 ,4.94d0 ,6.25d0 ,7.91d0 ,
     $          1.00e1/

      A = 4
      Z = 2

      OPEN(UNIT = 1, FILE = 'xb_values/he_xb.csv',STATUS = 'unknown')
      DO K1 = 1, NX
        READ(1,*) XB(K1)
      ENDDO
      CLOSE(1)

      OPEN(UNIT = 2, FILE = 'grids/he.grid')

c-----generate Q2 TMD grid
      DO 1 K1 = 1, NX-1
            zh = XB(K1)
            DO 3 K3 = 1, NB
               B = BB(K3)
               Q = 1.12292d0*dsqrt(1d0+4d0/9d0*b*b)/b
               call PDF_A_EPPS16(zh,Q,A,Z,U,D,UB,DB,S,SB,G,nloops)
!               call PDF_P(zh,Q,U,D,UB,DB,S,SB,G,nloops)
               WRITE(2,90) zh,dlog(b),
     *         b*u,b*d,b*s,b*g,b*ub,b*db,b*sb
*     
 3          CONTINUE
 1    CONTINUE
      close(2)


*
  90  FORMAT (9(E14.5))
*
      
      return
      end


