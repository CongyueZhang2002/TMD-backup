        PROGRAM main
        IMPLICIT NONE
    
        INTEGER :: CL, IH, IC, SET, nloops, REPLICA_NUM, AA
        DOUBLE PRECISION :: Z, Q2, U, UB, D, DB, S, SB, C, CB, B, BB, GL
        DOUBLE PRECISION :: Q0, zz, mub, fu,fub,fd,fdb,fs,fsb,fg,H_AA 
        INTEGER :: INIT, i
        DOUBLE PRECISION :: array(47)
        double precision eps

        CHARACTER*20 file_name1, file_name2, file_name3, file_name4
        CHARACTER*20 file_name5, file_name6, file_name7, file_name8
        INTEGER unit_number

        COMMON / INITIATE / INIT
        COMMON /HMASS/ H_AA 
        COMMON /REPLICA_INFO/ REPLICA_NUM
        COMMON /meson/ IH,IC
        parameter(eps=1d-10)            

    
        INIT = 0
    
        CL = 1  
        IH = 1  
        IC = 1  
        SET = 0    
        mub = 1.001D0
        Q2 = mub**2
        nloops = 1     
        !Q0 = dsqrt(1d0) - eps 
        Q0 = 1d0


        DATA array/0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
     4        0.095, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
     5        0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475,  0.5,
     6        0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7,
     7        0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9,
     8        0.925, 0.95, 0.975, 1.0/

        unit_number = 10

        file_name1 = 'like_H.csv'
        file_name2 = 'like_He.csv'
        file_name3 = 'like_C.csv'
        file_name4 = 'like_Ne.csv'
        file_name5 = 'new_H.csv'
        file_name6 = 'new_He.csv'
        file_name7 = 'new_C.csv'
        file_name8 = 'new_Ne.csv'        

        OPEN (UNIT = unit_number, FILE = file_name1, STATUS = 'REPLACE')
        INIT = 0        
        Do i=1,47
            Z = array(i)
            call LIKEn(CL,IH,IC,set,array(i),Q2,1,U,UB,D,DB,S,SB,C,CB,
     9      B,BB,GL)
            
        write(10,'(*(G0.6,:,","))') U/Z,UB/Z,D/Z,DB/Z,S/Z,SB/Z,GL/Z
            
        End do

        CLOSE (UNIT = unit_number)

        OPEN (UNIT = unit_number, FILE = file_name2, STATUS = 'REPLACE')
        INIT = 0
        Do i=1,47
            Z = array(i)
            call LIKEn(CL,IH,IC,set,array(i),Q2,4,U,UB,D,DB,S,SB,C,CB,
     9      B,BB,GL)
            
        write(10,'(*(G0.6,:,","))') U/Z,UB/Z,D/Z,DB/Z,S/Z,SB/Z,GL/Z
            
        End do

        CLOSE (UNIT = unit_number)

        OPEN (UNIT = unit_number, FILE = file_name3, STATUS = 'REPLACE')
        INIT = 0
        Do i=1,47
            Z = array(i)
            call LIKEn(CL,IH,IC,set,array(i),Q2,12,U,UB,D,DB,S,SB,C,CB,
     9      B,BB,GL)
            
        write(10,'(*(G0.6,:,","))') U/Z,UB/Z,D/Z,DB/Z,S/Z,SB/Z,GL/Z
            
        End do

        CLOSE (UNIT = unit_number)

        OPEN (UNIT = unit_number, FILE = file_name4, STATUS = 'REPLACE')
        INIT = 0
        Do i=1,47
            Z = array(i)
            call LIKEn(CL,IH,IC,set,array(i),Q2,20,U,UB,D,DB,S,SB,C,CB,
     9      B,BB,GL)
            
        write(10,'(*(G0.6,:,","))') U/Z,UB/Z,D/Z,DB/Z,S/Z,SB/Z,GL/Z
            
        End do

        CLOSE (UNIT = unit_number)

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

        OPEN (UNIT = unit_number, FILE = file_name5, STATUS = 'Unknown')
            
        Do i=1,47
            zz = array(i)
            Call Nuclear_TMDFF(zz,mub,AA,fu,fub,fd,fdb,fs,fsb,fg,nloops)
                
            write(10,'(*(G0.6,:,","))') fu,fub,fd,fdb,fs,fsb,fg
                
        End do

        CLOSE (UNIT = unit_number)

        AA=4
        H_AA=4

        call SetTimeLikeEvolution(.true.)
        call SetVFNS
        call SetPerturbativeOrder(1)
        call SetPDFEvolution("exactalpha")      
        call SetPDFSet("external")  
        call SetMaxFlavourPDFs(5)
        call SetMaxFlavourAlpha(5)
        call InitializeAPFEL
        
        call EvolveAPFEL(Q0,mub)

        OPEN (UNIT = unit_number, FILE = file_name6, STATUS = 'Unknown')
            
        Do i=1,47
            zz = array(i)
            Call Nuclear_TMDFF(zz,mub,AA,fu,fub,fd,fdb,fs,fsb,fg,nloops)
                
            write(10,'(*(G0.6,:,","))') fu,fub,fd,fdb,fs,fsb,fg
                
        End do

        CLOSE (UNIT = unit_number)

        AA=12
        H_AA=12

        call SetTimeLikeEvolution(.true.)
        call SetVFNS
        call SetPerturbativeOrder(1)
        call SetPDFEvolution("exactalpha")      
        call SetPDFSet("external")  
        call SetMaxFlavourPDFs(5)
        call SetMaxFlavourAlpha(5)
        call InitializeAPFEL
        
        call EvolveAPFEL(Q0,mub)

        OPEN (UNIT = unit_number, FILE = file_name7, STATUS = 'Unknown')
            
        Do i=1,47
            zz = array(i)
            Call Nuclear_TMDFF(zz,mub,AA,fu,fub,fd,fdb,fs,fsb,fg,nloops)
                
            write(10,'(*(G0.6,:,","))') fu,fub,fd,fdb,fs,fsb,fg
                
        End do

        CLOSE (UNIT = unit_number)

        AA=20
        H_AA=20

        call SetTimeLikeEvolution(.true.)
        call SetVFNS
        call SetPerturbativeOrder(1)
        call SetPDFEvolution("exactalpha")      
        call SetPDFSet("external")  
        call SetMaxFlavourPDFs(5)
        call SetMaxFlavourAlpha(5)
        call InitializeAPFEL
        
        call EvolveAPFEL(Q0,mub)

        OPEN (UNIT = unit_number, FILE = file_name8, STATUS = 'Unknown')
            
        Do i=1,47
            zz = array(i)
            Call Nuclear_TMDFF(zz,mub,AA,fu,fub,fd,fdb,fs,fsb,fg,nloops)
                
            write(10,'(*(G0.6,:,","))') fu,fub,fd,fdb,fs,fsb,fg
                
        End do

        CLOSE (UNIT = unit_number)

        STOP
        END PROGRAM main
