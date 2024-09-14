c-----ad, bN, NG1, NQ1, NG2, NQ2, BG1, GG1, DG1
      gammaMIN = 0.0 !BL( 1)
      gammaMAX = 0.0 !BUP(1)
c-----
      g3fMIN = 0.004 !BL( 2)
      g3fMAX = 0.012 !BUP(2)
c-----
      g3DMIN = 0.003 !BL( 3)
      g3DMAX = 0.01 !BUP(3)
c-----
      Nq1MIN = 0.0
      Nq1MAX = 0.5
c-----
      gq1MIN = 0.0
      gq1MAX = 0.2
c-----
      dq1MIN = 0.0
      dq1MAX = 0.003
c-----
      Nq2MIN = 0.0
      Nq2MAX = 0.5
c-----
      gq2MIN = 0.0
      gq2MAX = 0.8
c-----
      dq2MIN = 0.0
      dq2MAX = 2.0
c-----
      p_10MIN = 0.0 !BL(10)
      p_10MAX = 0.0 !BUP(10)
c-----
      p_11MIN = 0.0 !BL(11)
      p_11MAX = 0.0 !BUP(11)
c-----
      p_12MIN = 0.0 !BL(12)
      p_12MAX = 0.0 !BUP(12)

      CALL ran(r0 ,0) 
      r01000 = 1000d0*r0
      initseed = NINT(r01000)

      seed = initseed

      CALL ran(r1 ,seed)
      seed = seed+1
      CALL ran(r2 ,seed)
      seed = seed+1
      CALL ran(r3 ,seed)
      seed = seed+1
      CALL ran(r4 ,seed)
      seed = seed+1
      CALL ran(r5 ,seed)
      seed = seed+1
      CALL ran(r6 ,seed)
      seed = seed+1
      CALL ran(r7 ,seed)
      seed = seed+1
      CALL ran(r8 ,seed)
      seed = seed+1
      CALL ran(r9 ,seed)
      seed = seed+1
      CALL ran(r10 ,seed)
      seed = seed+1
      CALL ran(r11  ,seed)
      seed = seed+1
      CALL ran(r12  ,seed)
      seed = seed+1

C-----0
      !gammaBEST = 0.0 !gammaMIN  +(gammaMAX  -gammaMIN  )*r1
      !g3fBEST   = g3fMIN  +(g3fMAX  -g3fMIN  )*r2
      !g3DBEST   = g3DMIN  +(g3DMAX  -g3DMIN  )*r3
      !Nq1BEST   = Nq1MIN  +(Nq1MAX  -Nq1MIN  )*r4
      !gq1BEST   = gq1MIN  +(gq1MAX  -gq1MIN  )*r5
      !dq1BEST   = dq1MIN  +(dq1MAX  -dq1MIN  )*r6
      !Nq2BEST   = Nq2MIN  +(Nq2MAX  -Nq2MIN  )*r7
      !gq2BEST   = gq2MIN  +(gq2MAX  -gq2MIN  )*r8
      !dq2BEST   = dq2MIN  +(dq2MAX  -dq2MIN  )*r9
      !p_10BEST  = 0.0 !p_10MIN  +(p_10MAX  -p_10MIN  )*r10
      !p_11BEST  = 0.0 !p_11MIN  +(p_11MAX  -p_11MIN  )*r11
      !p_12BEST  = 0.0 !p_12MIN  +(p_12MAX  -p_12MIN  )*r12

C-----0 fitb
      gammaBEST = 2.1953732304972053     
      g3fBEST   = 0.21999864082553472     
      g3DBEST   = 1.9005551041957672E-002
      Nq1BEST   = 0.25573866909683662     
      gq1BEST   = 6.3511479101485691E-003
      dq1BEST   = 0.18404030632272672
      Nq2BEST   = 0.15627533798111459     
      gq2BEST   = 1.1501650731427362
      dq2BEST   = 0.47402585933913333 
      p_10BEST  = 0.0 
      p_11BEST  = 0.0 
      p_12BEST  = 0.0 

C-----2 fita
      !gammaBEST = 0.0
      !g3fBEST   = 8.1048849994903895E-003
      !g3DBEST   = 6.4951483140125487E-003
      !Nq1BEST   = 0.23778925560331921
      !gq1BEST   = 8.1915234408331372E-002
      !dq1BEST   = 1.6894577731978444E-003
      !Nq2BEST   = 0.25898511416172854
      !gq2BEST   = 0.38838810011739255
      !dq2BEST   = 1.0329952753649483
      !p_10BEST  = 0.0
      !p_11BEST  = 0.0
      !p_12BEST  = 0.0
  
C-----1(old)
      !gammaBEST = 0.0
      !g3fBEST   = 0.008
      !g3DBEST   = 0.0048
      !Nq1BEST   = 0.0322
      !gq1BEST   = -0.1555
      !dq1BEST   = -0.0451
      !Nq2BEST   = 0.4567
      !gq2BEST   = 0.4567
      !dq2BEST   = 0.4567
      !p_10BEST  = -0.0148
      !p_11BEST  = 0.4567
      !p_12BEST  = 0.0

      VSTRT(1)  =  gammaBEST 
      VSTRT(2)  =  g3fBEST
      VSTRT(3)  =  g3DBEST
      VSTRT(4)  =  Nq1BEST
      VSTRT(5)  =  gq1BEST
      VSTRT(6)  =  dq1BEST
      VSTRT(7)  =  Nq2BEST
      VSTRT(8)  =  gq2BEST
      VSTRT(9)  =  dq2BEST
      VSTRT(10)  =  p_10BEST
      VSTRT(11)  =  p_11BEST
      VSTRT(12)  =  p_12BEST
