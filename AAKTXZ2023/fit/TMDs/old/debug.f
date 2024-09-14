        PROGRAM main
        IMPLICIT NONE
    
        INTEGER :: CL, IH, IC, SET, nloops, REPLICA_NUM, AA
        DOUBLE PRECISION :: Q2
        DOUBLE PRECISION :: Q0, zz, mub, fu,fub,fd,fdb,fs,fsb,fg,H_AA 
        INTEGER :: INIT, i

        COMMON / INITIATE / INIT
        COMMON /HMASS/ H_AA 
        COMMON /REPLICA_INFO/ REPLICA_NUM
        COMMON /meson/ IH,IC          
    
        INIT = 0
    
        CL = 1  
        IH = 1  
        IC = 1  
        SET = 0    
        Q2 = 1.0D0
        mub = Q2
        nloops = 1     
        Q0 = 1d0

        AA=1
        H_AA=1

        call SetTimeLikeEvolution(.true.)
        call SetVFNS
        call SetPerturbativeOrder(1)
        call SetPDFEvolution("exactalpha")      
        call SetPDFSet("external")  
        call SetMaxFlavourPDFs(5)
        call SetMaxFlavourAlpha(5)
        call InitializeAPFEL
        
        call EvolveAPFEL(Q0,mub)
            
        Call Nuclear_TMDFF(0.1,mub,AA,fu,fub,fd,fdb,fs,fsb,fg,nloops)
        
        STOP
        END PROGRAM main
