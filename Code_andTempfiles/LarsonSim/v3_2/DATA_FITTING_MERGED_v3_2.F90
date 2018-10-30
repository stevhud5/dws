!==================================================================
!                    SHARED PARAMETER
!==================================================================

MODULE SHARED

  SAVE

  REAL(8) :: CLIST(10500,400),CONST,CREP,CFL,T(5000),G(5000),GF(5000)
  REAL(8) :: tMIN,tMAX,TAGET_MIU(40),TAGET_TAO(40),PC,PM,Ls,d,MERR(2)
  REAL(8) :: FQE(2000),GGE(2000),GGGE(2000),GGF(2000),GGGF(2000),FRLP
  REAL(8) :: FQC1,GGC1,GGGC1,FQM,GGM,GGGM,FQC2,FQFC1,GGFC1,GGGFC1,FQFC2
  REAL(8) :: FQFM,GGFM,GGGFM,LpMI,LpMA,APHMA,APHMI,RATIOMI,RATIOMA,ZEMI
  REAL(8) :: ZEMA,G0MA,G0MI,Kb,Vs,TEM,PHI,PI,RELAX,GERR(3),t0,VERR,AERR
  REAL(8) :: V0,RANG,INDEX,BDIND,DT,TRY_ERRG,TRY_MIU(40),TRY_TAO(40)
  real(8) :: g_rot(2000),gg_rot(2000),g_rouse_u(2000),gg_rouse_u(2000)
  real(8) :: g_bend_u(2000),gg_bend_u(2000),g_rouse(2000),gg_rouse(2000)
  real(8) :: g_bend(2000),gg_bend(2000)

  INTEGER :: LMAP(10500),LBREATH,SLP,NUM,NUMG,NUMIN,FLAGPR,FLAGD
  INTEGER :: NN,SAMPLE,GENERATION,NSEED,SERR(2),FLAGMR,II,FLAGIT
  INTEGER :: NIND,TRY_NN,FLAGN,FLAGFC2,FLAGMGT,Nt0,NUM1

END MODULE SHARED

!==================================================================
!				         MAIN PROGRAM START
!==================================================================

PROGRAM MAIN

  USE SHARED
  IMPLICIT none

  REAL(8) :: RATIOF,LpF,ZEF,APH,RATIO,Lp,ZE,FMI,FMA,PERR(6),FERR(5),G0F,WB
  REAL(8) :: MLIST1(400),MLIST2(400),MLIST3(400),GGG(5000),GG(5000),DW
  REAL(8) :: REMAIN,BR,RND,G1,G2,ERRG,GGS,GGGS,GGC,GGGC,D0,W,La,MARK,EAPH
  REAL(8) :: tbr,trep,te,Le,G0,ZETA,FIC,FIC2,TL,tp,LBAR,ST,OS,APHU,APHL
  REAL(8) :: FIC3,tR,FIC4,CLENTH,G0R,Dr,VOF,LENGTH(5050),Li,APH0,RATIO0
  REAL(8) :: G21,G22,C1,C2,DGMAX,APH1,GRATIO,FRATIO,DTL,SUM
  real(8) :: X,Y,XY,X2

  INTEGER :: I,J,JJ,JJJ,K,KK,KKK,TT,N,M,LSEGL,LENTH,SUM1,SUM2,INV,REM
  INTEGER :: FLAGC2,FLAGH,FLAGF,FLAGIM,FLAGG,IDC1,IDM,IDH,LDR,LSEG,FLAGSS
  INTEGER :: LROUSE,LBEND,NUMOU,OUTF,OUTW,OUTU,IDATE(8),NLFU,NLFL,JMAX
  INTEGER :: LBR,NEWL1,NEWL2,NEWL3,LPI(3),XCH,NM,Z,NLF,ITER,ISIZE,FLAG
  integer :: OUT_sep,OUT_fit,mode
  INTEGER,ALLOCATABLE :: ISEED(:)

  CHARACTER :: ANS1,ANS2,ANS3
  character(len=15) :: title

  !==================================================================
  !				         INPUT PARAMETER
  !==================================================================

  Kb=1.381E-23          ! BOLTZMANN CONSTANT
  PI=3.141592653   	! CONSTANT PI
  Vs=0.000891        	! SOLVENT VISCOSITY
  PHI=0.0662            ! VOLUME FRACTION	
  TEM=298.15		! SYSTEM TEMPERATURE

  CALL DATE_AND_TIME(VALUES=IDATE)
  CALL RANDOM_SEED(SIZE=ISIZE)
  ALLOCATE(ISEED(ISIZE))
  ISEED=IDATE(2)+IDATE(3)+IDATE(4)
  ISEED=ISEED*(IDATE(6)+IDATE(5)+IDATE(7))
  CALL RANDOM_SEED(PUT=ISEED)

  APHMA=25.		       ! THE UPPPER BOUNDARY OF Le/Lp
  APHMI=0.1		       ! THE LOWER BOUNDARY OF Le/Lp
  RATIOMA=2000.		       ! THE UPPPER BOUNDARY OF BREAKAGE RATE (ZETA)
  RATIOMI=1E-4                 ! THE LOWER BOUNDARY OF BREAKAGE RATE (ZETA)
  LpMA=2.E-7      	       ! THE UPPPER BOUNDARY OF Lp
  LpMI=1.5E-8     	       ! THE LOWER BOUNDARY OF Lp
  ZEMA=500.		       ! THE UPPPER BOUNDARY OF Z
  ZEMI=0.5 		       ! THE LOWER BOUNDARY OF Z

  DTL=1E-6    		       ! TIME STEP CONSTANT
  TT=10**8		       ! THE LARGEST NUMBER OF TIME STEP
  NLFU=200		       ! MAXIMUM INCREASE OF MICELLE NUMBER
  NLFL=-200		       ! MAXIMUM DECREASE OF MICELLE NUMBER
  RANG=5.		       ! THE RANGE OF MICELLE LENGTH DISTRIBUTION
  Ls=0.2		       ! THE RATIO OF SEGMENTAL LENGTH TO Lp

  LBREATH=1		       ! THE SWITCHER FOR CLF
  LDR=1			       ! THE SWITCHER FOR DOUBLE REPTATION
  LROUSE=1		       ! THE SWITCHER FOR ROUSE MODES
  LBEND=1		       ! THE SWITCHER FOR BENDING MODES
  SLP=0			       ! THE SWITCHER FOR MODIFTING lp DURING ITERATIONS

  PC=0.75		       ! PROBABILITY OF PERMUTATION
  PM=0.2  	 	       ! PROBABILITY OF MUTATION
  GENERATION=20000 	       ! GENERATION NUMBER
  NN=20			       ! APPROXIMATE PARAMETER NUMBER
  SAMPLE=3*NN		       ! SAMPLE SIZE NUMBER 
  NSEED=9 		       ! SEEDS NNUMBER

  JJJ=0			       ! TESTING NUMBER
  NUMOU=500		       ! NUMBER OF OUTPUT DATA
  INDEX=0.01 	               ! DATA RECORDING INTERVAL
  INV=5*10**5		       ! PROCESS RECORDING INTERVAL

  FLAGF=0		       ! FLAG FOR INCLUDING FREQUENCY INFORMATION				  
  FLAGMR=0		       ! FLAG FOR PURE MECHANICAL DATA
  FLAGMGT=0		       ! FLAG FOR USING MAGNITUDE

  OUTF=0	  	       ! THE SWITCHER FOR GENERATING FILE OF GT+GF
  OUTW=0		       ! THE SWITCHER FOR GENERATING FILE OF GW
  OUTU=0		       ! THE SWITCHER FOR GENERATING FILE OF GW+UNENTANGLED
  OUT_sep=0                    ! Switch to generate separate files for bending/Rouse modes
  OUT_fit=0                    ! Switch to run fitting procedure or not
  mode=0                       ! Switch to allow single exponential relaxation

  NUM=5000		       ! ENSEMBLE SIZE
  ITER=25		       ! ITERATIONS FOR PARAMETER FITTING
  FLAGIT=1		       ! FLAG FOR CHANGING PLATEAU MODULUS/PERSISTENCE LENGTH
  RELAX=0.75
  REM=0

  FLAGD=0		       ! FLAG FOR ENDING THE ITERATIONS AUTOMATICALLY
  FLAGG=0		       ! FLAG FOR FIXING TAO_I FOR GA
  FLAGN=0		       ! FLAG FOR INCREASING ENSEMBLE SIZE
  FLAGSS=0		       ! FLAG FOR STARTING SIMULATION FROM SPECIFIED PARAMETER VALUE

  OPEN(20,FILE='SIMULATION OUTPUT.DAT')	 ! PARAMETER AND FITTING ERROR EVOLUTION OUTPUT FILE
  OPEN(23,FILE='TIME_FREQUENCY TRANSFORMATION.DAT')	! MIU_I&TAO_I OUTPUT FILE 
  OPEN(24,FILE='SIMULATION MONITOR.DAT')  ! PROCESSION OF SIMULATION OUTPUT FILE
  OPEN(10,FILE='INTRADATA.DAT')	! SCRATCH FILE

  !==================================================================
  !				         READING INPUT FILE
  !==================================================================

  OPEN(16,FILE='INPUT_v3_2.DAT')	! EXPERIMENTAL DATA INPUT FILE
  READ(16,*) TITLE
  READ(16,*) TEM,PHI,Vs

  READ (16,*) ANS1,ANS2
  IF (ANS1.EQ.'Y') THEN
     FLAGF=1						 			          
  END IF
  IF (ANS2.EQ.'Y') THEN
     FLAGMR=1					 
  END IF

  READ(16,*) ANS1,ANS2,ANS3
  IF (ANS1.EQ.'Y') THEN
     OUTF=1
  END IF
  IF (ANS2.EQ.'Y') THEN
     OUTW=1
  END IF
  IF (ANS3.EQ.'Y') THEN
     OUTU=1
  END IF

  read(16,*) ANS1,ANS2,ANS3
  if (ans1.eq.'Y') then
     out_sep=1
  end if
  if (ans2.eq.'Y') then
     out_fit=1
  else
     out_fit=0
     iter=10
  end if
  if (ans3.eq.'Y') then
     mode=1
  end if

  READ(16,*) ANS1,Lp
  IF (ANS1.EQ.'Y') THEN
     SLP=1
  END IF

  READ(16,*) V0,Li,d
  Li=Li*10**(-6.)		     ! AVERAGE MICELLE LENGTH
  Lp=Lp*10**(-9.)		     ! PERSISTENT LENGTH 
  d=d*10**(-9.)			     ! MICELLE DIAMETER

  READ(16,*) NUM1,ITER
  READ(16,*) ANS1,RATIO0,APH0
  IF (ANS1.EQ.'Y') THEN
     FLAGSS=1
  END IF
  NLFU=NUM1/10			     ! MAXIMUM INCREASE OF MICELLE NUMBER
  NLFL=-NUM1/10			     ! MAXIMUM DECREASE OF MICELLE NUMBER
  
  NUMIN=0
  IF (FLAGF.EQ.1) THEN
     DO I=1,1000
	READ (16,*,END=10) FQE(I),GGE(I),GGGE(I)
	NUMIN=NUMIN+1
     END DO
  END IF

10 CLOSE(16)
  WRITE(24,*) 'READ DATA FINISH'
  WRITE(*,*) 'READ DATA FINISH'

  IF (FLAGF.EQ.0)	THEN
     WRITE(24,*) 'ERROR: FREQUENCY INFORMATION NOT EXSIST'
     WRITE(*,*) 'ERROR: FREQUENCY INFORMATION NOT EXSIST'
     STOP
  END IF
  IF (FLAGMR.EQ.1) THEN
     WRITE(24,*) 'WARNING: MECHANICAL RHEOMETRIC HIGH FREQUENCY DATA NEGLECTED'
     WRITE(*,*) 'WARNING: MECHANICAL RHEOMETRIC HIGH FREQUENCY DATA NEGLECTED'
  END IF

  WRITE(24,*) 'THE NUMBER OF DATA POINTS IS ',NUMIN 
  WRITE(*,*) 'THE NUMBER OF DATA POINTS IS ',NUMIN	  

  !==================================================================
  !	      FINDING EXPERIMENTAL LOCAL RHEOLOGICAL FEATURES 
  !==================================================================

  FLAGIM=0			  ! FLAG FOR FINDING EXPERIMENTAL INTERMEDIATE POINT
  FLAGH=0			  ! FLAG FOR HIGH FREQUENCY EXPERIMENTAL DATA AVAILABLE
  FLAGC2=0			  ! FLAG FOR FINDING EXPERIMENTAL 2nd CROSSOVER POINT

  FQC1=FQE(1)	            	  ! EXPERIMENTAL 1st CROSSOVER FREQUENCY 
  FQM=FQE(NUMIN)		  ! EXPERIMENTAL INTERMEDIATE FREQUENCY 
  GGC1=GGE(1)			  ! EXPERIMENTAL VALUE OF G' FOR 1st CROSSOVER  
  GGGC1=GGGE(1)			  ! EXPERIMENTAL VALUE OF G" FOR 1st CROSSOVER
  GGM=GGE(NUMIN)		  ! EXPERIMENTAL VALUE OF G' FOR INTERMEDIATE FREQUENCY 
  GGGM=GGGE(NUMIN)	    	  ! EXPERIMENTAL VALUE OF G" FOR INTERMEDIATE FREQUENCY

  IDC1=1			  ! INDEX OF DATA POINT FOR 1st CROSSOVER FREQUENCY
  IDM=NUMIN			  ! INDEX OF DATA POINT FOR INTERMEDIATE FREQUENCY
  IDH=NUMIN			  ! INDEX OF DATA POINT FOR REACHING HIGH FREQUENCY REGION

  DGMAX=(GGE(1)-GGGE(1))/GGGE(1)
  JMAX=1
  DO J=1,NUMIN
     IF ((GGE(J)-GGGE(J))/GGGE(J).GT.DGMAX) THEN
	DGMAX=(GGE(J)-GGGE(J))/GGGE(J)
	JMAX=J
     END IF
  END DO

  IF (DGMAX.LT.1.0) THEN
     FLAGMGT=1
  END IF

  GRATIO=1.0

  DO J=2,NUMIN-1
     IF ((GRATIO*GGE(J-1).LE.GGGE(J-1)).AND.(GRATIO*GGE(J+1).GE.GGGE(J+1))) THEN	
	FQC1=FQE(J)
	IDC1=J
	GGC1=GGE(J)
	GGGC1=GGGE(J)
	EXIT
     END IF
  END DO

  FRATIO=FQE(JMAX)/FQC1

  IF (FLAGMR.EQ.0) THEN
     DO J=IDC1+2,NUMIN-1
	IF ((GGGE(J-1).LE.GGE(J-1)).AND.(GGGE(J+1).GE.GGE(J+1))) THEN
           FQC2=FQE(J)
           FLAGC2=1
           EXIT
	END IF
     END DO
  END IF

  IF (FLAGMGT.EQ.0) THEN
     FQM=FRATIO*FQC1
     IDM=JMAX
     GGM=GGE(JMAX)
     GGGM=GGGE(JMAX)
     FLAGIM=1
  ELSE
     FQM=SQRT(FQC1*FQC2)
     DO J=2,NUMIN-1
        IF ((FQE(J-1).LE.FQM).AND.(FQE(J+1).GE.FQM)) THEN
           IDM=J
           GGM=GGE(J)
           GGGM=GGGE(J)
           FLAGIM=1
           EXIT
	END IF
     END DO
  END IF

  DO J=2,NUMIN-1
     IF ((FQE(J-1).LE.10*FQM).AND.(FQE(J+1).GE.10*FQM)) THEN
	IDH=J
	FLAGH=1
	EXIT
     END IF
  END DO

  IF (FLAGF.EQ.1)	THEN
     !     FMI=FLOOR(LOG10(FQE(1)))            ! THE MINIMUM FREQUENCY OF EXPERIMENTAL DATA 
     !     FMA=CEILING(LOG10(FQE(NUMIN)))      ! THE MAXIMUM FREQUENCY OF EXPERIMENTAL DATA
     FMI=-2.
     FMA=5.
     DW=(FMA-FMI)/NUMOU			        ! INTERVAL OF FREQUENCY 
  END IF

  IF (FLAGIM.EQ.0) THEN
     WRITE(24,*) 'ERROR: BROADER FREQUENCY RANGE REQUIRED'
     WRITE(*,*) 'ERROR: BROADER FREQUENCY RANGE REQUIRED'
     STOP
  END IF
  IF (FLAGH.EQ.0)	THEN
     WRITE(24,*) 'WARNING: EXPERIMENT NOT REACH HIGH FREQUENCY ZONE'
     WRITE(24,*) 'INITIAL GUESS USED FOR PERSISTENCE LENGTH'
     WRITE(*,*) 'WARNING: EXPERIMENT NOT REACH HIGH FREQUENCY ZONE'
     WRITE(*,*) 'INITIAL GUESS USED FOR PERSISTENCE LENGTH'
     SLP=0
  END IF
  IF (FLAGC2.EQ.0) THEN
     IF (FLAGMGT.EQ.0) THEN
	WRITE(24,*) 'WARNING: 2nd CROSSOVER FREQUENCY NOT DETECTED'
	WRITE(24,*) 'ACCURACY OF PERSISTENCE LENGTH NOT GUARANTEED'
	WRITE(*,*) 'WARNING: 2nd CROSSOVER FREQUENCY NOT DETECTED'
	WRITE(*,*) 'ACCURACY OF PERSISTENCE LENGTH NOT GUARANTEED'
     ELSE
 	WRITE(24,*) 'ERROR: CROSSOVER FREQUENCY NOT DETECTED'
	WRITE(24,*) 'PLEASE INCLUDE HIGH FREQEUNCY DATA OR CONSIDER OTHER CHARACTERIZATION METHODS'
	WRITE(*,*) 'ERROR: CROSSOVER FREQUENCY NOT DETECTED'
	WRITE(*,*) 'PLEASE INCLUDE HIGH FREQEUNCY DATA OR CONSIDER OTHER CHARACTERIZATION METHODS'
     END IF
  END IF

  !==================================================================
  !				           INITIAL GUESSES
  !==================================================================

  G0MA=GGM*10. 		! THE UPPPER BOUNDARY OF PLATEAU MODULUS
  G0MI=GGC1*1.	        ! THE LOWER BOUNDARY OF PLATEAU MODULUS

  C1=9.75*Kb*TEM/Lp**3./G0MA
  C2=3.*28./5./PI*Kb*TEM*PHI/d**2./G0MA/Lp
  APHL=1.
  APH1=APHL
  DO I=1,50
     IF ((C1*APHL**2.2+C2-3.*APHL).LT.(APHMI**4.)) THEN
	APHL=APHMI
     ELSE
	APHL=(C1*APHL**2.2+C2-3.*APHL)**0.25
     END IF
     APHL=(APH1+APHL)/2.
     APH1=APHL
  END DO
  IF (APHL.GT.APHMI) THEN
     APHMI=APHL
  END IF

  C1=9.75*Kb*TEM/Lp**3./G0MI
  C2=3.*28./5./PI*Kb*TEM*PHI/d**2./G0MI/Lp
  APHU=2.
  APH1=APHU
  DO I=1,50
     IF ((C1*APHU**2.2+C2-3.*APHU).GT.(APHMA**4.)) THEN
	APHU=APHMA
     ELSE
	APHU=(C1*APHU**2.2+C2-3.*APHU)**0.25
     END IF
     APHU=(APH1+APHU)/2.
     APH1=APHU
  END DO
  IF (APHU.LT.APHMA) THEN
     APHMA=APHU
  END IF

  APH=3.			    ! RATIO OF Le TO Lp
  EAPH=10.
  K=0
  FLAG=1
  Z01:DO WHILE (FLAG.EQ.1)
     EAPH=APH
     Le=APH*Lp       	            ! ENTANGLEMENT LENGTH 
     ZETA=Le**0.6*Lp**0.4           ! MESH SIZE
     FIC=2.*PI*Vs/LOG(ZETA/d)	    ! PARALLEL FRICTION COEFFICIENT 
     D0=Kb*TEM/FIC		    ! MOBILITY OF MICELLE
     TL=SQRT(APH/2.)	            ! RATIO OF MICELLE LENGTH TO TUBE LENGTH 
     IF (APH.LT.2) THEN
	TL=1.
     END IF

     ZE=Li/APH/Lp		    ! NUMBER OF ENTANGLEMENT
     La=Li/TL	    	            ! AVERAGE TUBE LENGTH
     trep=La**3./D0/PI**2.*TL	    ! REPTATION TIME
     IF (FLAGMGT.EQ.0) THEN
	RATIO=(APH**3./trep/FQC1)**1.5      ! RATIO OF BREAKAGE TO REPTATION TIME
     ELSE
	RATIO=10.
     END IF

     IF (FLAGMGT.EQ.0) THEN
 	G0=GGC1/(0.265-0.0657*LOG10(RATIO)) ! PLATEAU MODULUS
     ELSE
 	C1=9.75*Kb*TEM/Lp**3./APH**1.8
	C2=28./5./PI*Kb*TEM*PHI/d**2./APH/Lp
	G0=(C1*APH**3.+C2*3.)/(3.+APH**3.)	
     END IF
     IF (G0.GT.G0MA) THEN
	G0=G0MA
     END IF
     IF (G0.LT.G0MI) THEN
	G0=G0MI
     END IF

     C1=9.75*Kb*TEM/Lp**3./G0
     C2=3.*28./5./PI*Kb*TEM*PHI/d**2./G0/Lp
     APH1=APH						 
     DO I=1,50					  
	IF ((C1*APH**2.2+C2-3.*APH).GT.(APHU**4.)) THEN
           APH=APHU
	ELSE
           IF ((C1*APH**2.2+C2-3.*APH).LT.(APHL**4.)) THEN
              APH=APHL
           ELSE
              APH=(C1*APH**2.2+C2-3.*APH)**0.25
           END IF
	END IF
	APH=(APH1+APH)/2.
	APH1=APH
     END DO

     EAPH=ABS(EAPH-APH)/APH	   ! CORRECTION OF ALPHA
     K=K+1
     IF (FLAGMGT.EQ.0) THEN
	IF ((EAPH.LT.0.01).OR.(K.GT.10)) THEN
           FLAG=0
	END IF
     ELSE
	FLAG=0
     END IF

  END DO Z01

  !==================================================================
  !				      ITERATION PART
  !==================================================================

  JJ=1
  II=1
  SERR(1)=0.
  SERR(2)=0.
  MERR(1)=10.
  MERR(2)=10.
  AERR=0.

  Z00:DO WHILE (JJ.LE.ITER)	! ITERATION STARTS				 			
     FLAGPR=0				! FLAG FOR OVER-WRITING PARAMETER & OUTPUT FILE			        

     C1=9.75*Kb*TEM/Lp**3./G0
     C2=3.*28./5./PI*Kb*TEM*PHI/d**2./G0/Lp
     APH=2.						 ! INTIAL GUESS FOR RATIO OF Le TO Lp
     APH1=APH
     DO I=1,50					  
	IF ((C1*APH**2.2+C2-3.*APH).GT.(APHMA**4.)) THEN
           APH=APHMA
	ELSE
           IF ((C1*APH**2.2+C2-3.*APH).LT.(APHMI**4.)) THEN
              APH=APHMI
           ELSE
              APH=(C1*APH**2.2+C2-3.*APH)**0.25
           END IF
	END IF
	APH=(APH1+APH)/2.
	APH1=APH
     END DO

     IF ((FLAGSS.EQ.1).AND.(JJ.EQ.1).AND.(II.EQ.1)) THEN
	APH=APH0
	RATIO=RATIO0
	ZE=Li/APH/Lp
	C1=9.75*Kb*TEM/Lp**3./APH**1.8
	C2=28./5./PI*Kb*TEM*PHI/d**2./APH/Lp
	G0=(C1*APH**3.+C2*3.)/(3.+APH**3.)
     END IF

     TL=SQRT(APH/2.)	        ! RATIO OF MICELLE LENGTH TO TUBE LENGTH 
     IF (APH.LT.2) THEN
	TL=1.
     END IF
     La=ZE*Lp*APH/TL		! AVERAGE TUBE LENGTH     
     Le=APH*Lp    	        ! ENTANGLEMENT LENGTH 

     ZETA=Le**0.6*Lp**0.4       ! MESH SIZE
     LBAR=La/Lp			! RATIO OF TUBE LENGTH TO Lp
     FIC=2.*PI*Vs/LOG(ZETA/d)	! PARALLEL FRICTION COEFFICIENT 
     D0=Kb*TEM/FIC	        ! MOBILITY OF MICELLE
     trep=La**3./D0/PI**2.*TL	! REPTATION TIME	
     tbr=RATIO*trep	        ! BREAKAGE TIME
     te=2.*Le**2.*Lp/3/PI**2./D0	! EQUILIBRATION TIME
     FIC2=4.*PI*Vs/LOG(0.6*ZETA/d)      ! PERPENDICULAR FRICTION COEFFICIENT
     FIC3=3.*Kb*TEM/PI/Vs		! FRICTION OF MICELLAR RODS 

     IF (JJ.EQ.ITER) THEN
	GOTO  42
     END IF

     IF (JJ.GT.10) THEN
        if (flagn.eq.0) then
           serr(1)=0.
           serr(2)=0.
           merr(1)=10.
           merr(2)=10.
        end if
	FLAGN=1
     END IF

     IF (FLAGN.EQ.1) THEN
	NUM=NUM1
     ELSE 
	NUM=NUM1/2
     END IF

     WRITE(20,*) 'ITERATION: ',JJ,II
     WRITE(20,*) 'PLATEAU:',G0
     WRITE(20,*) 'RATIO:',RATIO
     WRITE(20,*) 'Trep:',trep
     WRITE(20,*) 'Ze:',ZE   
     WRITE(20,*) 'MICELLE LENGTH:',La*TL
     WRITE(20,*) 'ALPHA:',APH
     WRITE(20,*) 'PERSISTENCE:',Lp
     WRITE(20,*) 'DIAMETER:',d

     WRITE(*,*) 'ITERATION: ',JJ,II
     WRITE(*,*) 'PLATEAU:',G0
     WRITE(*,*) 'RATIO:',RATIO
     WRITE(*,*)	'Trep:',trep
     WRITE(*,*)	'Ze:',ZE
     WRITE(*,*) 'MICELLE LENGTH:',La*TL
     WRITE(*,*) 'ALPHA:',APH
     WRITE(*,*) 'PERSISTENCE:',Lp
     WRITE(*,*) 'DIAMETER:',d

     !==================================================================
     !			  DISCRITIZE MICELLE LENGTH DISTRIBUTION
     !==================================================================

     DT=tbr/NUM/2.
     IF (DT.GT.DTL) THEN
	DT=DTL
     END IF
     IF (La.LT.1E-6) THEN
        DT=DT/2.
     END IF
     IF (La.LT.5E-7) THEN
        DT=DT/2.
     END IF
     MARK=tbr/DT/NUM/2.	

     N=1			     ! INTERVAL OF DISCRETIZED MICELLE TUBE LENGTH
     Ls=EXP(RANG)/NUM*LBAR
     M=CEILING(RANG*LBAR/Ls)	     ! MAXIMUM LENGTH IN THE ENSEMBLE
     DO J=1,10000
  	SUM1=0
        DO I=1,M		     ! GENERATE MICELLE FOLLOWED BY ITS LENGTH DISTRIBUTION
           SUM1=SUM1+CEILING(NUM*N*Ls*EXP(-1.*I*N*Ls/LBAR)/LBAR)
	END DO
	REMAIN=SUM1-NUM	
        IF (ABS(REMAIN).LE.1) THEN	 		 
           EXIT
	ELSE
           Ls=Ls/NUM*SUM1
           M=CEILING(RANG*LBAR/Ls)
        END IF
     END DO

21   WRITE(24,*) 'LINEAR LENGTH DISCRETIZATION FINISH'
     WRITE(*,*) 'LINEAR LENGTH DISCRETIZATION FINISH'

     WRITE(24,*) 'THE MINIMUM TUBE LENGTH IS',Lp*N*Ls
     WRITE(24,*) 'THE MAXIMUM TUBE LENGTH IS',M*N*Lp*Ls
     WRITE(*,*) 'THE MINIMUM TUBE LENGTH IS',Lp*N*Ls
     WRITE(*,*) 'THE MAXIMUM TUBE LENGTH IS',M*N*Lp*Ls

     SUM2=0.
     LSEG=0
     LSEGL=0
     A00:DO I=1,M	
	LENGTH(I)=CEILING(NUM*N*Ls*EXP(-1.*I*N*Ls/LBAR)/LBAR)
	LSEG=LSEG+I*N*LENGTH(I)
	IF ((I.EQ.M).AND.(SUM2+LENGTH(I).NE.NUM)) THEN
           LENGTH(I)=NUM-SUM2
	END IF

	SUM1=1+SUM2
        SUM2=SUM2+LENGTH(I)
	IF (I*N*Ls*TL.GE.APH) THEN
           LSEGL=LSEGL+I*N*LENGTH(I)
	END IF

	DO J=SUM1,SUM2
           LMAP(J)=I*N
           IF (I*N*Ls*TL.LT.APH) THEN
              CLIST(J,1)=0.
           ELSE			
              CLIST(J,1)=1.
              CLIST(J,2)=0.
              CLIST(J,3)=LMAP(J)
           END IF
	END DO
     END DO A00
     WRITE(24,*) 'CHAIN LIST INITIATION FINISH'
     WRITE(*,*) 'CHAIN LIST INITIATION FINISH'

     !==================================================================
     !	       LOW FREQUENCY SIMULATION FOR ENTANGLED MICELLES
     !==================================================================

     KKK=1
     NLF=0
     G2=1.   ! *
     Nt0=TT  ! *

     CREP=2.*D0*DT/(Lp*Ls)**3./TL
     CONST=2.*(CREP*TL/3./PI/Ls)**0.25
     CFL=4.*CREP*TL/3./PI/Ls

     A01:DO K=1,TT		! EVOLUTION OF ENSEMBLE DUE TO RELAXATION AND BREAKAGE 
        OPEN(21,FILE='TEMP.DAT')

	IF (K.EQ.1) THEN
           KK=0
	END IF

	DO I=1,NUM+NLF
           LENTH=2*INT(CLIST(I,1))+1
           DO J=1,LENTH				 
              MLIST1(J)=CLIST(I,J)
           END DO
           NEWL1=LMAP(I)
           CALL MPDR(MLIST1,NEWL1,MLIST2)	
           LENTH=2*INT(MLIST2(1))+1
           DO J=1,LENTH
              CLIST(I,J)=MLIST2(J)
           END DO
	END DO

	IF (MOD(K,INV).EQ.0) THEN
           WRITE(24,*) K,'REPTATION FINISH'
           WRITE(*,*) K,'REPTATION FINISH'
	END IF

        B01:IF (K.GT.MARK*KK) THEN		! BREAKAGE AND REJOINING OCCURS
           KK=KK+1  	 
           CALL RANDOM_NUMBER(BR)		! BREAKAGE 

           C01:IF ((BR.LT.0.5D0).OR.(NLF.LE.NLFL)) THEN
31            CALL RANDOM_NUMBER(RND)
              LBR=INT(RND*LSEG)
              IF (LBR.LT.1) THEN
                 GOTO 31
              END IF
              DO I=1,NUM+NLF
                 LBR=LBR-LMAP(I)
                 IF (LBR.EQ.0) THEN
                    GOTO 31
                 END IF
                 IF (LBR.LT.0) THEN
                    LBR=LBR+LMAP(I)
                    LPI(1)=I
                    EXIT
                 END IF
              END DO

              NEWL1=LMAP(LPI(1))
              LENTH=2*INT(CLIST(LPI(1),1))+1
              DO I=1,LENTH				 
                 MLIST1(I)=CLIST(LPI(1),I)
              END DO
              CALL BRKG(LBR,MLIST1,NEWL1,MLIST2,NEWL2,MLIST3,NEWL3)

              LMAP(LPI(1))=NEWL2
              LENTH=2*INT(MLIST2(1))+1
              DO I=1,LENTH
                 CLIST(LPI(1),I)=MLIST2(I)
              END DO
              LMAP(NUM+NLF+1)=NEWL3
              LENTH=2*INT(MLIST3(1))+1
              DO I=1,LENTH
                 CLIST(NUM+NLF+1,I)=MLIST3(I)
              END DO

              IF (MOD(K,INV).EQ.0) THEN
                 WRITE(24,*) 'CHAIN BREAKAGE OCCUR:',K
                 WRITE(*,*) 'CHAIN BREAKAGE OCCUR:',K
              END IF
              NLF=NLF+1

           END IF C01

           C02:IF ((BR.GT.0.5D0).OR.(NLF.GE.NLFU)) THEN	! REJOINING
32            CALL RANDOM_NUMBER(RND)
              LPI(1)=INT(RND*(NUM+NLF))
              IF (LPI(1).LT.1) THEN
                 GOTO 32
              END IF
33            CALL RANDOM_NUMBER(RND)
              LPI(2)=INT(RND*(NUM+NLF))
              IF (LPI(2).LT.1) THEN
                 GOTO 33
              END IF

              IF (LPI(2).EQ.LPI(1)) THEN
                 GOTO 32
              END IF
              IF ((LMAP(LPI(1))+LMAP(LPI(2))).GT.M*N) THEN
                 GOTO 32
              END IF
              if ((clist(lpi(1),1)+clist(lpi(2),1)).gt.199) then
                 goto 32
              end if

              NEWL1=LMAP(LPI(1))
              LENTH=2*INT(CLIST(LPI(1),1))+1
              DO I=1,LENTH
                 MLIST1(I)=CLIST(LPI(1),I)
              END DO
              NEWL2=LMAP(LPI(2))
              LENTH=2*INT(CLIST(LPI(2),1))+1
              DO I=1,LENTH
                 MLIST2(I)=CLIST(LPI(2),I)
              END DO
              CALL REJN(MLIST1,NEWL1,MLIST2,NEWL2,MLIST3,NEWL3)

              LMAP(LPI(1))=NEWL3
              LENTH=2*INT(MLIST3(1))+1
              DO I=1,LENTH
                 CLIST(LPI(1),I)=MLIST3(I)
              END DO

              DO I=LPI(2),NUM+NLF-1
                 LMAP(I)=LMAP(I+1)
                 LENTH=2*INT(CLIST(I+1,1))+1
                 DO J=1,LENTH
                    CLIST(I,J)=CLIST(I+1,J)
                 END DO
              END DO

              IF (MOD(K,INV).EQ.0) THEN
                 WRITE(24,*) 'CHAIN RECOMBATION OCCUR:',K
                 WRITE(*,*) 'CHAIN RECOMBATION OCCUR:',K
              END IF
              NLF=NLF-1

           END IF C02

	END IF B01

	G21=G2			! CALCULATE STRESS RELAXATION FUNCTION
	G1=0.
	DO I=1,NUM+NLF
           IF (INT(CLIST(I,1)).GT.0) THEN
              DO J=1,INT(CLIST(I,1))
                 G1=G1+(CLIST(I,2*J+1)-CLIST(I,2*J))/LSEGL
              END DO
           END IF
	END DO
	IF (LDR.EQ.1) THEN
           G2=G1**2.
	ELSE
           G2=G1
	END IF
	G22=G2

	IF ((G21.GE.0.9).AND.(G22.LE.0.9)) THEN	 ! *
           Nt0=K	  ! *
           NUMG=1	  ! *
	END IF	  ! *

	IF (K.GE.Nt0) THEN
           CALL GGT(K,G21,G22)	! PICK APPROPRIATE NUMBER OF THE DATA FROM TEMPORARY FILE
	END IF

	IF (MOD(K,INV).EQ.0) THEN
           WRITE(24,*) K,NLF,G2
           WRITE(*,*) K,NLF,G2
	END IF

	IF (G2.LT.0.001) THEN
           REM=K
           WRITE(21,*)	-1,K*DT,G2
           WRITE(*,*) 'REM',K
           WRITE(24,*)	'REM',K
           CLOSE(21)
           EXIT
	END IF

     END DO A01

     NUMG=NUMG-2
     WRITE(*,*) 'INITIATION OF GA FINISHED'
     WRITE(24,*) 'INITIATION OF GA FINISHED'

     !==================================================================
     !				TIME TO FREQUENCY TRANSFORMATION
     !==================================================================

     tMAX=LOG(T(NUMG))
     tMIN=LOG(T(1))
     ERRG=10.

     X=0.
     Y=0.
     XY=0.
     X2=0.

     if (mode.eq.1) then
        do I=1,NUMG
           X=X+T(I)
           Y=Y+log(G(I))
           XY=XY+T(I)*log(G(I))
           X2=X2+T(I)*T(I)
        end do
        NN=1
        TAGET_TAO(1)=(X*Y-real(NUMG)*XY)/(real(NUMG)*X2-X*X)
        TAGET_MIU(1)=(Y+TAGET_TAO(1)*X)/real(NUMG)
        write(*,*) 'tao:',TAGET_TAO(1)
        write(*,*) 'weight:',exp(TAGET_MIU(1))
        errg=0.
        do I=1,NUMG
           errg=errg+abs(G(I)-exp(TAGET_MIU(1))*exp(-TAGET_TAO(1)*T(I)))/G(I)
           GF(I)=exp(TAGET_MIU(1))*exp(-TAGET_TAO(1)*T(I))
        end do
        errg=errg/real(NUMG)
        write(*,*) 'error:',errg
        write(*,*) NUMG,X,Y,XY,X2
     end if
     
     if ((ERRG.gt.0.05).or.(mode.eq.0)) then
        DO I=1,3
           CALL GEAM(FLAGG,TRY_NN,TRY_MIU,TRY_TAO,TRY_ERRG)
           IF (TRY_ERRG.LT.ERRG) THEN
              NN=TRY_NN
              DO J=1,NN
                 TAGET_MIU(J)=TRY_MIU(J)
                 TAGET_TAO(J)=TRY_TAO(J)
              END DO
              ERRG=TRY_ERRG
           END IF
           WRITE(*,*) 'TRY',I,' ERR:',TRY_ERRG
           WRITE(23,*) 'TRY',I,' ERR:',TRY_ERRG
           IF (ERRG.LT.1E-2) THEN
              EXIT
           END IF
        END DO
     end if

     WRITE(23,*) 'ITERATION: ',JJ,II
     DO I=1,NN
	WRITE(23,*) EXP(TAGET_MIU(I)),TAGET_TAO(I)
	WRITE(*,*) EXP(TAGET_MIU(I)),TAGET_TAO(I)
     END DO

     WRITE(23,*) 'GA ERROR IS',ERRG
     WRITE(24,*) 'GA ERROR IS',ERRG
     WRITE(*,*) 'GA ERROR IS',ERRG

     WRITE(24,*) 'RELAXATION TIME CALCULATION FINISH'
     WRITE(*,*) 'RELAXATION TIME CALCULATION FINISH'

     IF (OUTF.EQ.1) THEN
 	OPEN(7,FILE='GF_t.DAT')	  ! WRITE OUTPUT FILE FOR FITTED G(t) GENERATE BY GA
        WRITE(7,*) 'ITERATION: ',JJ,II	 
 	WRITE(7,*) 'PARAMETERS PART'
 	WRITE(7,*) 'RATIO',RATIO
 	WRITE(7,*) 'TIME INTERVAL',DT	
 	WRITE(7,*) 'MICELLE LENGTH',LBAR*TL 	
        WRITE(7,*) 'PLATAEU',G0
 	WRITE(7,*) 'Le',Le 
        WRITE(7,*) 'Lp',Lp    
 	WRITE(7,*) 'DATA PART'
        DO I=1,NUMG
           WRITE(7,*) T(I),G(I),GF(I)
        END DO
	CLOSE(7)
     END IF

     VOF=0.			! VOLUME FRACTION FOR ENTANGLED MICELLES 
     DO I=1,M
	CLENTH=I*N*Ls*TL
	IF (CLENTH.GE.APH) THEN
           VOF=VOF+LENGTH(I)/NUM*CLENTH/TL/LBAR
	END IF
     END DO

     DO I=1,NUMOU		! GENERATE SIMULATED G'&G" FOR LOW FREQUENCY
	W=10**(FMI+DW*I)
	GG(I)=0.
	GGG(I)=0.
	DO J=1,NN
           GG(I)=GG(I)+G0*EXP(TAGET_MIU(J))*(W*TAGET_TAO(J))**2./(1.+(W*TAGET_TAO(J))**2.)
           GGG(I)=GGG(I)+G0*EXP(TAGET_MIU(J))*(W*TAGET_TAO(J))/(1.+(W*TAGET_TAO(J))**2.)
	END DO
     END DO

     DO I=1,NUMIN
	GGF(I)=0.
	GGGF(I)=0.
	DO J=1,NN
           GGF(I)=GGF(I)+G0*EXP(TAGET_MIU(J))*(FQE(I)*TAGET_TAO(J))**2./(1.+(FQE(I)*TAGET_TAO(J))**2.)
           GGGF(I)=GGGF(I)+G0*EXP(TAGET_MIU(J))*(FQE(I)*TAGET_TAO(J))/(1.+(FQE(I)*TAGET_TAO(J))**2.)
	END DO
     END DO

     IF (OUTW.EQ.1) THEN
	OPEN(12,FILE='GW.DAT')	   ! WRITE OUTPUT FILE FOR NORMALIZED LOW FREQUENCY G'&G"
        WRITE(12,*) 'ITERATION: ',JJ,II   
	WRITE(12,*) 'PARAMETERS PART'
	WRITE(12,*) 'RATIO',RATIO
	WRITE(12,*) 'TIME INTERVAL',DT	
	WRITE(12,*) 'MICELLE LENGTH',LBAR*TL 	
        WRITE(12,*) 'PLATAEU',G0
	WRITE(12,*) 'Le',Le 
        WRITE(12,*) 'Lp',Lp    
	WRITE(12,*) 'DATA PART'
        DO I=1,NUMOU
           WRITE(12,*) 10**(FMI+DW*I),GG(I),GGG(I)
        END DO
	CLOSE(12)
     END IF

     !==================================================================
     !			 CONTRIBUTION FROM UNENTANGLED MICELLES
     !==================================================================

     DO I=1,NUMIN	 ! ROTARY RELAXATION FOR UNENTANGLED MICELLES (L<Le)
	DO J=1,M
           CLENTH=J*N*Ls*TL
           IF ((CLENTH.GE.3.*d/Lp).AND.(CLENTH.LT.APH)) THEN
              Dr=FIC3*(LOG(CLENTH*Lp/d)-0.8)/(CLENTH*Lp)**3.
              VOF=LENGTH(J)/NUM*CLENTH/TL/LBAR
              G0R=4.*PHI/d**2./PI*Kb*TEM/(CLENTH*Lp)*VOF
              GGF(I)=GGF(I)+(FQE(I)/Dr)**2./(1.+(FQE(I)/Dr)**2.)*G0R
              GGGF(I)=GGGF(I)+(FQE(I)/Dr)/(1.+(FQE(I)/Dr)**2.)*G0R
           END IF
        END DO
     END DO

     g_rot=0.
     gg_rot=0.
     DO I=1,NUMOU		  
	W=10**(FMI+DW*I)	  
	DO J=1,M
           CLENTH=J*N*Ls*TL
           IF ((CLENTH.GE.3.*d/Lp).AND.(CLENTH.LT.APH)) THEN
              Dr=FIC3*(LOG(CLENTH*Lp/d)-0.8)/(CLENTH*Lp)**3.
              VOF=LENGTH(J)/NUM*CLENTH/TL/LBAR
              G0R=4.*PHI/d**2./PI*Kb*TEM/(CLENTH*Lp)*VOF
              GG(I)=GG(I)+(W/Dr)**2./(1.+(W/Dr)**2.)*G0R
              GGG(I)=GGG(I)+(W/Dr)/(1.+(W/Dr)**2.)*G0R
              g_rot(I)=g_rot(I)+(W/Dr)**2./(1.+(W/Dr)**2.)*G0R
              gg_rot(I)=gg_rot(I)+(W/Dr)/(1.+(W/Dr)**2.)*G0R
           END IF
        END DO
     END DO

     IF (OUT_sep.EQ.1) THEN
	OPEN(28,FILE='rotary.dat')  ! WRITE OUTPUT FILE FOR NORMALIZED G'&G" WITH CONTRIBUTIONS FROM UNETANGLED MICELLES
 write(28,*) 'clenth',M,N*Ls*TL,3.*d/Lp,APH
	WRITE(28,*) 'ITERATION: ',JJ,II
	WRITE(28,*) 'PARAMETERS PART'
	WRITE(28,*) 'RATIO',RATIO
	WRITE(28,*) 'TIME INTERVAL',DT	
	WRITE(28,*) 'MICELLE LENGTH',LBAR*TL 	
        WRITE(28,*) 'PLATAEU',G0 
        WRITE(28,*) 'Le',Le 
        WRITE(28,*) 'Lp',Lp 
 	WRITE(28,*) 'DATA PART'
        DO I=1,NUMOU
           WRITE(28,*) 10**(FMI+DW*I),g_rot(I),gg_rot(I)
        END DO
 	CLOSE(28)
     END IF

     
     A04:DO I=1,NUMIN		 ! ROUSE MODES FOR UNENTANGLED MICELLES (bK<L<Le) 
	GGS=0.
        GGGS=0.

	DO J=1,M
           CLENTH=J*N*Ls*TL
           NM=NINT(CLENTH/2.)
           FIC4=2.*PI*Vs/LOG(CLENTH*Lp/d)
           IF ((CLENTH.GE.2.).AND.(CLENTH.LT.APH)) THEN
              tR=2.*(CLENTH*Lp)**2.*Lp/3/PI**2./Kb/TEM*FIC4
              VOF=LENGTH(J)/NUM*CLENTH/TL/LBAR
              G0R=4.*PHI/d**2./PI*Kb*TEM/(CLENTH*Lp)*VOF
              GGC=0.
              GGGC=0.
              DO K=1,NM
                 GGC=GGC+(FQE(I)*tR/2/K**2.)**2./(1+(FQE(I)*tR/2/K**2.)**2.)
                 GGGC=GGGC+FQE(I)*tR/2/K**2./(1+(FQE(I)*tR/2/K**2.)**2.)
              END DO
              GGS=GGS+GGC*G0R
              GGGS=GGGS+GGGC*G0R
           END IF
	END DO

	GGF(I)=GGF(I)+GGS
	GGGF(I)=GGGF(I)+GGGS
     END DO A04

     g_rouse_u=0.
     gg_rouse_u=0.
     
     A05:DO I=1,NUMOU		   
	W=10**(FMI+DW*I)
	GGS=0.
        GGGS=0	

	DO J=1,M
           CLENTH=J*N*Ls*TL
           NM=NINT(CLENTH/2.)
           FIC4=2.*PI*Vs/LOG(CLENTH*Lp/d)
           IF ((CLENTH.GE.2.).AND.(CLENTH.LT.APH)) THEN
              tR=2.*(CLENTH*Lp)**2.*Lp/3/PI**2./Kb/TEM*FIC4
              VOF=LENGTH(J)/NUM*CLENTH/TL/LBAR
              G0R=4.*PHI/d**2./PI*Kb*TEM/(CLENTH*Lp)*VOF
              GGC=0.
              GGGC=0.
              DO K=1,NM
                 GGC=GGC+(W*tR/2/K**2.)**2./(1+(W*tR/2/K**2.)**2.)
                 GGGC=GGGC+W*tR/2/K**2./(1+(W*tR/2/K**2.)**2.)
              END DO
              GGS=GGS+GGC*G0R
              GGGS=GGGS+GGGC*G0R
           END IF
	END DO

 	GG(I)=GG(I)+GGS
        GGG(I)=GGG(I)+GGGS
        g_rouse_u(I)=g_rouse_u(I)+GGS
        gg_rouse_u(I)=gg_rouse_u(I)+GGGS
     END DO A05

     IF (OUT_sep.EQ.1) THEN
	OPEN(29,FILE='rouse_u.dat')  ! WRITE OUTPUT FILE FOR NORMALIZED G'&G" WITH CONTRIBUTIONS FROM UNETANGLED MICELLES
	WRITE(29,*) 'ITERATION: ',JJ,II
	WRITE(29,*) 'PARAMETERS PART'
	WRITE(29,*) 'RATIO',RATIO
	WRITE(29,*) 'TIME INTERVAL',DT	
	WRITE(29,*) 'MICELLE LENGTH',LBAR*TL 	
        WRITE(29,*) 'PLATAEU',G0 
        WRITE(29,*) 'Le',Le 
        WRITE(29,*) 'Lp',Lp 
 	WRITE(29,*) 'DATA PART'
        DO I=1,NUMOU
           WRITE(29,*) 10**(FMI+DW*I),g_rouse_u(I),gg_rouse_u(I)
        END DO
 	CLOSE(29)
     END IF
     
     DO I=1,NUMIN		  ! BENDING MODES FOR UNENTANGLED MICELLES (L<Le)
	DO J=1,M
           CLENTH=J*N*Ls*TL
           VOF=LENGTH(J)/NUM*CLENTH/TL/LBAR
           IF ((CLENTH.GE.3.*d/Lp).AND.(CLENTH.LT.APH)) THEN
              Dr=FIC3*(LOG(CLENTH*Lp/d)-0.8)/(CLENTH*Lp)**3.
              WB=(FQE(I)/Dr)**2./(1.+(FQE(I)/Dr)**2.)	 ! *
              IF (CLENTH.LT.1.) THEN	! BENDING TIME
                 tp=2*Lp**(5./3.)*CLENTH**4.*FIC2/Kb/TEM
              ELSE
                 tp=2*Lp**(5./3.)*FIC2/Kb/TEM 
              END IF
              ST=0.3827*Kb*TEM*tp**0.75/15./PI*4.*d**(-2.)*PHI
              OS=0.9239/0.3827*ST  		        
              GGF(I)=GGF(I)+(VOF*ST*FQE(I)**0.75)*WB
              GGGF(I)=GGGF(I)+(VOF*(OS*FQE(I)**0.75+FQE(I)*Vs))*WB
           END IF
	END DO
     END DO

     g_bend_u=0.
     gg_bend_u=0.
     
     DO I=1,NUMOU
	W=10**(FMI+DW*I)					  
	DO J=1,M
           CLENTH=J*N*Ls*TL
           VOF=LENGTH(J)/NUM*CLENTH/TL/LBAR
           IF ((CLENTH.GE.3.*d/Lp).AND.(CLENTH.LT.APH)) THEN
              Dr=FIC3*(LOG(CLENTH*Lp/d)-0.8)/(CLENTH*Lp)**3.
              WB=(W/Dr)**2./(1.+(W/Dr)**2.)	
              IF (CLENTH.LT.1.) THEN	! BENDING TIME
                 tp=2*Lp**(5./3.)*CLENTH**4.*FIC2/Kb/TEM
              ELSE
                 tp=2*Lp**(5./3.)*FIC2/Kb/TEM
              END IF
              ST=0.3827*Kb*TEM*tp**0.75/15./PI*4.*d**(-2.)*PHI 
              OS=0.9239/0.3827*ST          
              GG(I)=GG(I)+(VOF*ST*W**0.75)*WB
              GGG(I)=GGG(I)+(VOF*(OS*W**0.75+W*Vs))*WB
              g_bend_u(I)=g_bend_u(I)+(VOF*ST*W**0.75)*WB
              gg_bend_u(I)=gg_bend_u(I)+(VOF*(OS*W**0.75+W*Vs))*WB
           END IF
	END DO
     END DO

     IF (OUT_sep.EQ.1) THEN
	OPEN(30,FILE='bend_u.dat')  ! WRITE OUTPUT FILE FOR NORMALIZED G'&G" WITH CONTRIBUTIONS FROM UNETANGLED MICELLES
	WRITE(30,*) 'ITERATION: ',JJ,II
	WRITE(30,*) 'PARAMETERS PART'
	WRITE(30,*) 'RATIO',RATIO
	WRITE(30,*) 'TIME INTERVAL',DT	
	WRITE(30,*) 'MICELLE LENGTH',LBAR*TL 	
        WRITE(30,*) 'PLATAEU',G0 
        WRITE(30,*) 'Le',Le 
        WRITE(30,*) 'Lp',Lp 
 	WRITE(30,*) 'DATA PART'
        DO I=1,NUMOU
           WRITE(30,*) 10**(FMI+DW*I),g_bend_u(I),gg_bend_u(I)
        END DO
 	CLOSE(30)
     END IF
     
     IF (OUTU.EQ.1) THEN
	OPEN(9,FILE='UNENTANGLE.DAT')  ! WRITE OUTPUT FILE FOR NORMALIZED G'&G" WITH CONTRIBUTIONS FROM UNETANGLED MICELLES
	WRITE(9,*) 'ITERATION: ',JJ,II
	WRITE(9,*) 'PARAMETERS PART'
	WRITE(9,*) 'RATIO',RATIO
	WRITE(9,*) 'TIME INTERVAL',DT	
	WRITE(9,*) 'MICELLE LENGTH',LBAR*TL 	
        WRITE(9,*)  'PLATAEU',G0 
        WRITE(9,*)  'Le',Le 
        WRITE(9,*)  'Lp',Lp 
 	WRITE(9,*) 'DATA PART'
        DO I=1,NUMOU
           WRITE(9,*) 10**(FMI+DW*I),GG(I),GGG(I)
        END DO
 	CLOSE(9)
     END IF

     !==================================================================
     !				    HIGH FREQUENCY SIMULATION
     !==================================================================

     A02:IF (LROUSE.EQ.1) THEN  

        g_rouse=0.
        gg_rouse=0.
        
        B02:DO I=1,NUMOU		   ! ROUSE MODES FOR ENTANGLED MICELLES
           W=10**(FMI+DW*I)
           GGS=0.
           GGGS=0.

           DO J=1,M
              CLENTH=J*N*Ls*TL
              Z=NINT(CLENTH/APH)
              NM=NINT(CLENTH/2.)
              IF (CLENTH.GE.APH) THEN
                 tR=2.*(CLENTH*Lp)**2.*Lp/3/PI**2./Kb/TEM*FIC
                 VOF=LENGTH(J)/NUM*CLENTH/TL/LBAR
                 GGC=0.
                 GGGC=0.
                 DO K=Z,NM
                    GGC=GGC+(W*tR/2/K**2.)**2./(1+(W*tR/2/K**2.)**2.)
                    GGGC=GGGC+W*tR/2/K**2./(1+(W*tR/2/K**2.)**2.)
                 END DO
                 GGS=GGS+GGC*5./4.*G0*VOF/Z
                 GGGS=GGGS+GGGC*5./4.*G0*VOF/Z
              END IF
           END DO

           GG(I)=GG(I)+GGS
           GGG(I)=GGG(I)+GGGS
           g_rouse(I)=g_rouse(I)+GGS
           gg_rouse(I)=gg_rouse(I)+GGGS
	END DO B02


        B03:DO I=1,NUMIN
           GGS=0.
           GGGS=0.

           DO J=1,M
              CLENTH=J*N*Ls*TL
              Z=NINT(CLENTH/APH)
              NM=NINT(CLENTH/2.)
              IF (CLENTH.GE.APH) THEN
                 tR=2.*(CLENTH*Lp)**2.*Lp/3/PI**2./Kb/TEM*FIC
                 VOF=LENGTH(J)/NUM*CLENTH/TL/LBAR
                 GGC=0.
                 GGGC=0.
                 DO K=Z,NM
                    GGC=GGC+(FQE(I)*tR/2/K**2.)**2./(1+(FQE(I)*tR/2/K**2.)**2.)
                    GGGC=GGGC+FQE(I)*tR/2/K**2./(1+(FQE(I)*tR/2/K**2.)**2.)
                 END DO
                 GGS=GGS+GGC*5./4.*G0*VOF/Z
                 GGGS=GGGS+GGGC*5./4.*G0*VOF/Z
              END IF
           END DO

           GGF(I)=GGF(I)+GGS
           GGGF(I)=GGGF(I)+GGGS

	END DO B03

	WRITE(24,*) 'ROUSE MODES INCLUDE'
        WRITE(*,*) 'ROUSE MODES INCLUDE'

        if (out_sep.eq.1) then
           OPEN(31,FILE='rouse.dat')  ! output file for Rouse modes
           WRITE(31,*) 'ITERATION: ',JJ,II
           WRITE(31,*) 'PARAMETERS PART'
           WRITE(31,*) 'RATIO',RATIO
           WRITE(31,*) 'TIME INTERVAL',DT	
           WRITE(31,*) 'MICELLE LENGTH',LBAR*TL 	
           WRITE(31,*)  'PLATAEU',G0 
           WRITE(31,*)  'Le',Le 
           WRITE(31,*)  'Lp',Lp 
           WRITE(31,*) 'DATA PART'
           DO I=1,NUMOU
              WRITE(31,*) 10**(FMI+DW*I),g_rouse(I),gg_rouse(I)
           END DO
           CLOSE(31)
        end if

     END IF A02

     g_bend=0.
     gg_bend=0.
     
     IF (LBEND.EQ.1) THEN    	! BENDING MODES
	DO I=1,NUMOU
           W=10**(FMI+DW*I)
           DO J=1,M
              CLENTH=J*N*Ls*TL
              VOF=LENGTH(J)/NUM*CLENTH/TL/LBAR
              IF (CLENTH.GE.APH) THEN	    ! BENDING TIME
                 tR=te*(CLENTH/APH)**2.	 
                 WB=(W*tR)**2./(1.+(W*tR)**2.)      
                 tp=2*Lp**(5./3.)*FIC2/Kb/TEM
                 ST=0.3827*Kb*TEM*tp**0.75/15./PI*4.*d**(-2.)*PHI 
                 OS=0.9239/0.3827*ST          
                 GG(I)=GG(I)+(VOF*ST*W**0.75)*WB
                 GGG(I)=GGG(I)+(VOF*(OS*W**0.75+W*Vs))*WB
                 g_bend(I)=g_bend(I)+(VOF*ST*W**0.75)*WB
                 gg_bend(I)=gg_bend(I)+(VOF*(OS*W**0.75+W*Vs))*WB
              END IF
           END DO
	END DO

	DO I=1,NUMIN
           DO J=1,M
              CLENTH=J*N*Ls*TL
              VOF=LENGTH(J)/NUM*CLENTH/TL/LBAR
              IF (CLENTH.GE.APH) THEN	
                 tR=te*(CLENTH/APH)**2.	
                 WB=(FQE(I)*tR)**2./(1.+(FQE(I)*tR)**2.) 
                 tp=2*Lp**(5./3.)*FIC2/Kb/TEM
                 ST=0.3827*Kb*TEM*tp**0.75/15./PI*4.*d**(-2.)*PHI 
                 OS=0.9239/0.3827*ST          
                 GGF(I)=GGF(I)+(VOF*ST*FQE(I)**0.75)*WB
                 GGGF(I)=GGGF(I)+(VOF*(OS*FQE(I)**0.75+FQE(I)*Vs))*WB
              END IF
           END DO
	END DO
	WRITE(24,*) 'BENDING MODES INCLUDE'
        WRITE(*,*) 'BENDING MODES INCLUDE'

        if (out_sep.eq.1) then
           OPEN(32,FILE='bending.dat')  ! WRITE OUTPUT FILE FOR NORMALIZED G'&G" WITH CONTRIBUTIONS FROM UNETANGLED MICELLES
           WRITE(32,*) 'ITERATION: ',JJ,II
           WRITE(32,*) 'PARAMETERS PART'
           WRITE(32,*) 'RATIO',RATIO
           WRITE(32,*) 'TIME INTERVAL',DT	
           WRITE(32,*) 'MICELLE LENGTH',LBAR*TL 	
           WRITE(32,*)  'PLATAEU',G0 
           WRITE(32,*)  'Le',Le 
           WRITE(32,*)  'Lp',Lp 
           WRITE(32,*)	'DATA PART'
           DO I=1,NUMOU
              WRITE(32,*) 10**(FMI+DW*I),g_bend(I),gg_bend(I)
           END DO
           CLOSE(32)
        end if

     END IF

     WRITE(24,*) JJ,II,'LOOP FINISH'
     WRITE(24,*) 'REM',REM
     WRITE(*,*) JJ,II,'LOOP FINISH'
     WRITE(*,*) 'REM',REM

     !==================================================================
     !			 FINDING SIMULATED LOCAL RHEOLOGICAL FEATURE
     !==================================================================

     VERR=0.
     DO I=1,3
	VERR=VERR+SQRT(GGGF(I)**2.+GGF(I)**2.)/FQE(I)/3.
     END DO
     VERR=LOG(VERR/V0)

     FLAGFC2=0			 ! FLAG FOR FINDING SIMULATED 2nd CROSSOVER POINT

     FQFC1=10**(FMI+DW)	         ! SIMULATED 1st CROSSOVER FREQUECNY 
     FQFM=10**(FMI+DW*NUMOU)	 ! SIMULATED INTERMEDIATE FREQUECNY 
     GGFC1=GG(1)		 ! SIMULATED VALUE OF G' FOR 1st CROSSOVER  
     GGGFC1=GGG(1)		 ! SIMULATED VALUE OF G" FOR 1st CROSSOVER
     GGFM=GG(NUMOU)		 ! SIMULATED VALUE OF G' FOR INTERMEDIATE FREQUECNY 
     GGGFM=GGG(NUMOU)	    	 ! SIMULATED VALUE OF G" FOR INTERMEDIATE FREQUECNY

     DO J=2,NUMOU-1
	IF ((GRATIO*GG(J-1).LE.GGG(J-1)).AND.(GRATIO*GG(J+1).GE.GGG(J+1))) THEN	
           FQFC1=10**(FMI+DW*J)
           GGFC1=GG(J)
           GGGFC1=GGG(J)
           XCH=J
           EXIT
	END IF
     END DO

     IF (FLAGMR.EQ.0) THEN
	DO J=XCH+1,NUMOU-1
           IF ((GGG(J-1).LE.GG(J-1)).AND.(GGG(J+1).GE.GG(J+1))) THEN
              FQFC2=10**(FMI+DW*J)
              FLAGFC2=1
              EXIT
           END IF
	END DO
     END IF

     IF (FLAGMGT.EQ.0) THEN
	FQFM=FRATIO*FQFC1
     ELSE
	FQFM=SQRT(FQFC1*FQFC2)
     END IF
     DO J=2,NUMOU-1
 	SUM=10**(FMI+DW*J)
        IF ((SUM/10**DW.LE.FQFM).AND.(SUM*10**DW.GE.FQFM)) THEN
           GGFM=GG(J)
           GGGFM=GGG(J)
           EXIT
	END IF
     END DO

     !==================================================================
     !			CALCULATING FITTING ERRORS AND DIFFERENCE
     !==================================================================				

     CALL CFER(IDC1,IDM,IDH,FERR,PERR)
     IF (FLAGMR.EQ.1) THEN
	FRLP=0.5*FERR(3)
     END IF

     WRITE(20,*) 'FERROR1:',FERR(1)
     WRITE(20,*) 'FERROR2:',FERR(2)
     WRITE(20,*) 'FERROR3:',FERR(3)
     WRITE(20,*) 'FERROR4:',FERR(4)

     IF (FLAGN.EQ.0) THEN
	WRITE(20,*) 'VERROR:',VERR
     ELSE
	WRITE(20,*) '*VERROR:',VERR	
     END IF

     WRITE(20,*) 'PERROR1:',PERR(1)
     WRITE(20,*) 'PERROR2:',PERR(2)
     WRITE(20,*) 'PERROR3:',PERR(3)
     WRITE(20,*) 'PERROR4:',PERR(4)
     WRITE(20,*) 'PERROR5:',PERR(5)

     WRITE(*,*) 'FERROR1:',FERR(1)
     WRITE(*,*) 'FERROR2:',FERR(2)
     WRITE(*,*) 'FERROR3:',FERR(3)
     WRITE(*,*) 'FERROR4:',FERR(4)

     IF (FLAGN.EQ.0) THEN
	WRITE(*,*) 'VERROR:',VERR
     ELSE
	WRITE(*,*) '*VERROR:',VERR	
     END IF

     WRITE(*,*) 'PERROR1:',PERR(1)
     WRITE(*,*) 'PERROR2:',PERR(2)
     WRITE(*,*) 'PERROR3:',PERR(3)
     WRITE(*,*) 'PERROR4:',PERR(4)
     WRITE(*,*) 'PERROR5:',PERR(5)

     IF (FLAGPR.EQ.1) THEN

	OPEN(27,FILE='NEW_INPUT.DAT')	! OUTPUT FILE WITH MINIMUM FITTING ERROR
	WRITE(27,*) TITLE  
        WRITE(27,*) TEM,PHI,Vs

	ANS1='N'
	ANS2='N'
 	IF (FLAGF.EQ.1) THEN	
           ANS1='Y'					 			          
	END IF
	IF (FLAGMR.EQ.1) THEN
           ANS2='Y'						 
	END IF
        WRITE(27,*) ANS1,'	',ANS2

	ANS1='N'
	ANS2='N'
	ANS3='N'
	IF (OUTF.EQ.1) THEN
           ANS1='Y'
	END IF
 	IF (OUTW.EQ.1) THEN
           ANS2='Y'
	END IF
	IF (OUTU.EQ.1) THEN
           ANS3='Y'
	END IF
	WRITE(27,*) ANS1,'	',ANS2,'	',ANS3

	ANS1='N'
	IF (SLP.EQ.1) THEN
           ANS1='Y'
	END IF
	WRITE(27,*) ANS1,'	',Lp*10.**9.

	WRITE(27,*) V0,LBAR*TL*Lp*10.**6.,d*10.**9.
	WRITE(27,*) NUM,ITER
	WRITE(27,*) 'Y','	',RATIO,APH
	CLOSE(27)

     END IF

     WRITE(10,*) 'ITERATION: ',JJ,II
     WRITE(10,*) 'PARAMETERS PART'
     WRITE(10,*) 'RATIO',RATIO
     WRITE(10,*) 'TIME INTERVAL',DT
     WRITE(10,*) 'REPTATION TIME',trep
     WRITE(10,*) 'BREAKAGE TIME',trep*RATIO
     WRITE(10,*) 'MICELLE LENGTH',LBAR*TL*Lp
     WRITE(10,*) 'PLATAEU',G0	
     WRITE(10,*) 'Le',Le 
     WRITE(10,*) 'Lp',Lp  
     WRITE(10,*) 'd',d 	  
     WRITE(10,*) 'DATA PART'
     DO I=1,NUMOU
	WRITE(10,*) 10**(FMI+DW*I),GG(I),GGG(I)
     END DO

     if (flagpr.eq.1) then
        open(11,file='result_fit.dat')
        WRITE(11,*) 'ITERATION: ',JJ,II
        WRITE(11,*) 'PARAMETERS PART'
        WRITE(11,*) 'RATIO',RATIO
        WRITE(11,*) 'TIME INTERVAL',DT
        WRITE(11,*) 'REPTATION TIME',trep
        WRITE(11,*) 'BREAKAGE TIME',trep*RATIO
        WRITE(11,*) 'MICELLE LENGTH',LBAR*TL*Lp
        WRITE(11,*) 'PLATAEU',G0	
        WRITE(11,*) 'Le',Le 
        WRITE(11,*) 'Lp',Lp  
        WRITE(11,*) 'd',d 	  
        WRITE(11,*) 'DATA PART'
        DO I=1,NUMOU
           WRITE(11,*) 10**(FMI+DW*I),GG(I),GGG(I)
        END DO
        close(11)
     end if
     
     IF (FLAGD.EQ.1) THEN
	WRITE(*,*) JJ,II,' ITERATION SUCCEED'
	WRITE(24,*) JJ,II,' ITERATION SUCCEED'
	WRITE(20,*) 'THE RESULT CONVERGED HERE'
	GOTO 42
     END IF

     !==================================================================
     !		    	OPTIMIZATION OF MICELLE PARAMETERS
     !==================================================================

     RATIOF=RATIO
     ZEF=ZE
     LpF=Lp
     G0F=G0

     A03:IF (JJ.LT.ITER) THEN	
  	CALL PRDF(G0F,RATIOF,ZEF,LpF,G0,RATIO,ZE,Lp,FLAGIT)	

	IF (FLAGIT.EQ.0) THEN
           II=1
           JJ=JJ+1
           FLAGIT=1
	END IF

	WRITE(20,*) 'G0:',G0
	WRITE(*,*) 'G0:',G0
	WRITE(20,*) 'RATIO:',RATIO
	WRITE(*,*) 'RATIO:',RATIO
	WRITE(20,*) 'Ze:',ZE
	WRITE(*,*) 'Ze:',ZE
	WRITE(20,*) 'Lp:',Lp
	WRITE(*,*) 'Lp:',Lp
     END IF A03

     if (out_fit.eq.0) then
        goto 42
     else
        WRITE(*,*) 'GOTO NEXT LOOP'
        WRITE(24,*) 'GOTO NEXT LOOP'
     end if

  END DO Z00

  !==================================================================
  !		    	FINAL RESULT OF FITTING PARAMETERS
  !==================================================================

42 IF ((FLAGD.EQ.1).OR.(JJ.EQ.ITER).or.(out_fit.eq.0)) THEN
     WRITE(*,*) 'FINAL RESULT:'
     WRITE(*,*) 'ITERATION: ',JJ,II
     WRITE(*,*) 'PLATAEU:',G0
     WRITE(*,*) 'RATIO:',RATIO	
     WRITE(*,*) 'MICELLE LENGTH:',LBAR*TL*Lp	
     WRITE(*,*) 'Le:',Le 
     WRITE(*,*) 'Lp:',Lp 
     WRITE(*,*) 'd:',d

     OPEN(25,FILE='RESULT.DAT')	! OUTPUT FILE WITH MINIMUM FITTING ERROR
     WRITE(25,*) 'FINAL RESULT:'  
     WRITE(25,*) 'ITERATION: ',JJ,II
     WRITE(25,*) 'PARAMETERS PART'
     WRITE(25,*) 'RATIO',RATIO
     WRITE(25,*) 'TIME INTERVAL',DT
     WRITE(25,*) 'REPTATION TIME',trep
     WRITE(25,*) 'BREAKAGE TIME',trep*RATIO
     WRITE(25,*) 'MICELLE LENGTH',LBAR*TL*Lp
     WRITE(25,*) 'PLATAEU',G0	
     WRITE(25,*) 'Le',Le 
     WRITE(25,*) 'Lp',Lp  
     WRITE(25,*) 'd',d 	  
     WRITE(25,*) 'DATA PART'
     DO I=1,NUMOU
	WRITE(25,*) 10**(FMI+DW*I),GG(I),GGG(I)
     END DO
     CLOSE(25)

     WRITE(20,*) 'FINAL RESULT:'
     WRITE(20,*) 'ITERATION: ',JJ,II
     WRITE(20,*) 'PLATAEU:',G0
     WRITE(20,*) 'RATIO:',RATIO	
     WRITE(20,*) 'MICELLE LENGTH:',LBAR*TL*Lp	
     WRITE(20,*) 'Le:',Le 
     WRITE(20,*) 'Lp:',Lp 
     WRITE(20,*) 'd:',d

  ELSE
     WRITE(*,*) 'SIMULATION UNCONVERGED'
     WRITE(*,*) 'NEW INPUT GENERATED'
     WRITE(20,*) 'SIMULATION UNCONVERGED'
     WRITE(20,*) 'NEW INPUT GENERATED'		
  END IF

  CLOSE(20)
  CLOSE(23)
  CLOSE(24)
  CLOSE(10)

  STOP

END PROGRAM MAIN

!==================================================================
!                          END HERE
!==================================================================






























!==================================================================
!      SUBROUTINE 1: MOVEMENT OF POINTERS DUE TO RELAXATION
!==================================================================

SUBROUTINE MPDR(CLISTi1,LMAPi1,CLISTo1)

  USE SHARED

  REAL(8) :: DL1,DL2(2),FL(2),RNDs1
  REAL(8), INTENT(IN) :: CLISTi1(400)
  REAL(8), INTENT(OUT) :: CLISTo1(400)
  INTEGER :: Is1,Js1,ITERs1,LENTHs1,NUMPs1
  INTEGER, INTENT(IN) :: LMAPi1

  ITERs1=20
  DO Is1=1,2*INT(CLISTi1(1))+1 ! THE POSTION OF POINTERS 
     CLISTo1(Is1)=CLISTi1(Is1)
  END DO

  S01:IF (CLISTo1(1).GT.0) THEN	 
     LENTHs1=2*INT(CLISTo1(1))+1

     IF (LBREATH.EQ.1) THEN		 ! CLFS
	FL(1)=CLISTo1(2)
        FL(2)=LMAPi1-CLISTo1(LENTHs1)
	DO Is1=1,2
           IF (FL(Is1).LE.0.2) THEN
              DL2(Is1)=CONST
           ELSE
              DL2(Is1)=FL(Is1)
              DO Js1=1,ITERs1
                 DL2(Is1)=CFL/(FL(Is1)+DL2(Is1)/2.)**3.
              END DO
           END IF
	END DO
	CLISTo1(2)=CLISTo1(2)+DL2(1)
	CLISTo1(LENTHs1)=CLISTo1(LENTHs1)-DL2(2)
     END IF

     IF (CLISTo1(2).GT.CLISTo1(3)) THEN
	DO Is1=2,INT(CLISTo1(1))
           CLISTo1(2*Is1-2)=CLISTo1(2*Is1)
           CLISTo1(2*Is1-1)=CLISTo1(2*Is1+1)
	END DO
	CLISTo1(1)=CLISTo1(1)-1
     END IF

     IF (CLISTo1(LENTHs1-1).GT.CLISTo1(LENTHs1)) THEN
	CLISTo1(1)=CLISTo1(1)-1
     END IF
  END IF S01

  DL1=SQRT(CREP/LMAPi1)
  CALL RANDOM_NUMBER(RNDs1)
  S02:IF (CLISTo1(1).GT.0) THEN
     IF (RNDs1.LT.0.5) THEN		 ! REPTATION

	NUMPs1=0
	DO Is1=2,2*INT(CLISTo1(1))+1 
           IF (CLISTo1(Is1)+DL1.LE.LMAPi1) THEN
              CLISTo1(Is1)=CLISTo1(Is1)+DL1
              NUMPs1=NUMPs1+1
           END IF
	END DO
	IF (MOD(NUMPs1,2).NE.0) THEN
           LENTHs1=NUMPs1+2
           CLISTo1(LENTHs1)=LMAPi1
           CLISTo1(1)=(NUMPs1+1)/2
	ELSE
           CLISTo1(1)=NUMPs1/2
	END IF

     ELSE
	NUMPs1=0
	DO Is1=2,2*INT(CLISTo1(1))+1 
           IF (CLISTo1(Is1)-DL1.GE.0.) THEN
              CLISTo1(Is1)=CLISTo1(Is1)-DL1
              NUMPs1=NUMPs1+1
           END IF
	END DO
	IF (MOD(NUMPs1,2).NE.0) THEN
           CLISTo1(2)=0.
           LENTHs1=2*INT(CLISTo1(1))-NUMPs1
           DO Is1=LENTHs1+2,2*INT(CLISTo1(1))+1
              CLISTo1(Is1-LENTHs1+1)=CLISTo1(Is1)
           END DO
           CLISTo1(1)=(NUMPs1+1)/2
	ELSE
           LENTHs1=2*INT(CLISTo1(1))-NUMPs1
           DO Is1=LENTHs1+2,2*INT(CLISTo1(1))+1
              CLISTo1(Is1-LENTHs1)=CLISTo1(Is1)
           END DO
           CLISTo1(1)=NUMPs1/2
	END IF

     END IF
  END IF S02

  S03:DO Is1=1,INT(CLISTo1(1))		  ! ANNILATION OF POINTERS
     IF (CLISTo1(1).GT.0) THEN

	IF (CLISTo1(2*Is1).GT.CLISTo1(2*Is1+1)) THEN
           DO Js1=Is1,INT(CLISTo1(1)-1)
              CLISTo1(2*Js1)=CLISTo1(2*(Js1+1))
              CLISTo1(2*Js1+1)=CLISTo1(2*(Js1+1)+1)
           END DO
           CLISTo1(1)=CLISTo1(1)-1
	END IF

	IF ((CLISTo1(2*Is1+1).GT.CLISTo1(2*Is1+2)).AND.(Is1.LT.CLISTo1(1))) THEN
           DO Js1=Is1,INT(CLISTo1(1)-2)
              CLISTo1(2*Js1+1)=CLISTo1(2*(Js1+1)+1)
              CLISTo1(2*Js1+2)=CLISTo1(2*(Js1+1)+2)
           END DO
           LENTHs1=2*INT(CLISTo1(1))+1
           CLISTo1(LENTHs1-2)=CLISTo1(LENTHs1)
           CLISTo1(1)=CLISTo1(1)-1
	END IF

     END IF
  END DO S03

  RETURN

END SUBROUTINE MPDR

!==================================================================
!      SUBROUTINE 2: BREAKAGE
!==================================================================

SUBROUTINE BRKG(NUMi2,CLISTi2,LMAPi2,CLIST1o2,LMAP1o2,CLIST2o2&
     &,LMAP2o2)

  USE SHARED

  REAL(8), INTENT(IN) :: CLISTi2(400)
  REAL(8), INTENT(OUT) :: CLIST1o2(400),CLIST2o2(400)
  INTEGER :: Is2,LENTHs2,NUMP1s2,NUMP2s2,FLAGs2
  INTEGER, INTENT(IN) :: NUMi2,LMAPi2
  INTEGER, INTENT(OUT) :: LMAP1o2,LMAP2o2

  LMAP1o2=NUMi2
  LMAP2o2=LMAPi2-NUMi2

  IF (CLISTi2(1).EQ.0) THEN
     CLIST1o2(1)=0
     CLIST2o2(1)=0
  END IF

  S04:IF (CLISTi2(1).GT.0) THEN

     NUMP1s2=0
     LENTHs2=2*INT(CLISTi2(1))+1
     DO Is2=2,LENTHs2
	IF (CLISTi2(Is2).LE.NUMi2) THEN
           NUMP1s2=NUMP1s2+1
	END IF
     END DO

     DO Is2=2,NUMP1s2+1
	CLIST1o2(Is2)=CLISTi2(Is2)
     END DO

     IF (MOD(NUMP1s2,2).NE.0) THEN
	CLIST1o2(NUMP1s2+2)=LMAP1o2
	CLIST1o2(1)=(NUMP1s2+1)/2
     ELSE
	CLIST1o2(1)=NUMP1s2/2
     END IF

     NUMP2s2=2*INT(CLISTi2(1))-NUMP1s2
     IF (MOD(NUMP2s2,2).NE.0) THEN
	CLIST2o2(2)=0.
	CLIST2o2(1)=(NUMP2s2+1)/2
	FLAGs2=1
     ELSE
	CLIST2o2(1)=NUMP2s2/2
	FLAGs2=0
     END IF

     DO Is2=NUMP1s2+2,LENTHs2
 	CLIST2o2(Is2-NUMP1s2+FLAGs2)=CLISTi2(Is2)-NUMi2
     END DO

  END IF S04

  RETURN

END SUBROUTINE BRKG

!==================================================================
!      SUBROUTINE 3: REJOINING
!==================================================================

SUBROUTINE REJN(CLIST1i3,LMAP1i3,CLIST2i3,LMAP2i3,CLISTo3,LMAPo3)

  USE SHARED

  REAL(8), INTENT(IN) :: CLIST1i3(400),CLIST2i3(400)
  REAL(8), INTENT(OUT) :: CLISTo3(400)
  INTEGER :: Is3,LENTH1s3,LENTH2s3
  INTEGER, INTENT(IN) :: LMAP1i3,LMAP2i3
  INTEGER, INTENT(OUT) :: LMAPo3

  LMAPo3=LMAP1i3+LMAP2i3
  CLISTo3(1)=CLIST1i3(1)+CLIST2i3(1)

  LENTH1s3=2*INT(CLIST1i3(1))+1
  IF (CLIST1i3(1).GT.0) THEN	
     DO Is3=2,LENTH1s3
	CLISTo3(Is3)=CLIST1i3(Is3)
     END DO
  ELSE
  END IF

  LENTH2s3=2*INT(CLIST2i3(1))+1
  IF (CLIST2i3(1).GT.0) THEN	
     DO Is3=2,LENTH2s3
	CLISTo3(LENTH1s3+Is3-1)=CLIST2i3(Is3)+LMAP1i3
     END DO
  END IF

  RETURN

END SUBROUTINE REJN

!==================================================================
!      SUBROUTINE 4: GENERATE G(t) 
!==================================================================

SUBROUTINE GGT(Ni4,G1i4,G2i4)

  USE SHARED

  REAL(8), INTENT(IN) :: G1i4,G2i4
  INTEGER :: FLAGs4
  INTEGER, INTENT(IN) :: Ni4

  S05:IF (G2i4.GT.SQRT(0.1)) THEN
     FLAGs4=0
     IF (Ni4.EQ.Nt0) THEN
	BDIND=Nt0*DT
	NIND=0
     ELSE
	IF (FLOOR(LOG10(Ni4*DT)).EQ.CEILING(LOG10((Ni4-1)*DT)))	THEN
           BDIND=Ni4*DT
           NIND=0
	END IF
     END IF

  ELSE
     FLAGs4=1
     IF (G1i4.GT.SQRT(0.1)) THEN 	
	BDIND=G2i4
	NIND=0
     ELSE
 	IF (FLOOR(LOG10(G2i4)).EQ.CEILING(LOG10(G1i4))) THEN
           BDIND=G2i4
           NIND=0
	END IF
     END IF
  END IF S05

  IF (FLAGs4.EQ.0) THEN				 ! WRITTING TIME DOMAIN DATA IN THE TEMPORARY FILE
     IF (Ni4*DT.GE.BDIND*10.**(INDEX*NIND)) THEN	
	WRITE(21,*)	Ni4,Ni4*DT,G2i4
	NIND=NIND+1
	T(NUMG)=Ni4*DT
	G(NUMG)=G2i4
	NUMG=NUMG+1
     END IF
  ELSE
     IF (G2i4.LE.BDIND/10.**(INDEX*NIND)) THEN	
	WRITE(21,*)	Ni4,Ni4*DT,G2i4
	NIND=NIND+1
	T(NUMG)=Ni4*DT
	G(NUMG)=G2i4
	NUMG=NUMG+1
     END IF
  END IF

  RETURN 

END	SUBROUTINE GGT

!==================================================================
!      SUBROUTINE 5: GENETIC ALGORITHM 
!==================================================================

SUBROUTINE GEAM(FLAGi5,No5,Mo5,To5,ERRMo5)

  USE SHARED

  REAL(8) :: RNDs5,MIU(505,40),IND1s5,IND2s5,TAO(505,40),OPT,Ds5,ERRSD(20)
  REAL(8) :: RSUMs5,ERR(505),MINs5,R1s5,R2s5,CRs5,RA1s5,RA2s5,BDRs5(2,40)
  REAL(8) :: XCHMs5(2,2,40),SEEDs5(20,2,40),OPT0,CNSTRs5
  REAL(8), INTENT(OUT) :: Mo5(40),To5(40),ERRMo5
  INTEGER :: Is5,Js5,Ks5,IIs5,JJs5,RESCALE1,RESCALE2,CHECK,POINT,N0s5
  INTEGER :: XCHs5,Ms5,PT1(2),PT2(2),NMOD(505),LOCs5(2),FLAGs5,NMSD(20)
  INTEGER :: COUNT
  INTEGER, INTENT(IN) :: FLAGi5
  INTEGER, INTENT(OUT) :: No5

  CRs5=3.
  RA1s5=LOG(G(NUMG)/10.) ! *
  RA2s5=LOG(1.)	 ! *
  N0s5=FLOOR(0.5*(tMAX-tMIN))
  CNSTRs5=0.01

  S06:DO Is5=1,SAMPLE			  ! GENERATE RANDOM ENSEMBLE FOR GA FITTING PARAMETERS
     CALL RANDOM_NUMBER(RNDs5)
     NMOD(Is5)=CEILING(RNDs5*(N0s5-2)+1.25)
     DO Js5=1,NMOD(Is5)     
	CALL RANDOM_NUMBER(RNDs5)
	MIU(Is5,Js5)=RNDs5*RA2s5+(1-RNDs5)*RA1s5
	IF (FLAGi5.EQ.0) THEN

02         CALL RANDOM_NUMBER(RNDs5)
           RSUMs5=RNDs5*tMAX+(1-RNDs5)*tMIN
           XCHMs5=0.
           DO IIs5=1,Js5-1
              XCHMs5(1,1,IIs5)=TAO(Is5,IIs5)
           END DO

           IIs5=Js5-1
           DO WHILE (IIs5.GT.0)
              IF (RSUMs5.EQ.XCHMs5(1,1,IIs5)) THEN
                 GOTO 02
              END IF
              IF (RSUMs5.LT.XCHMs5(1,1,IIs5)) THEN
                 XCHMs5(1,1,IIs5+1)=XCHMs5(1,1,IIs5)
                 IF (IIs5.EQ.1) THEN
                    XCHMs5(1,1,1)=RSUMs5
                 END IF
              END IF
              IF (RSUMs5.GT.XCHMs5(1,1,IIs5)) THEN
                 XCHMs5(1,1,IIs5+1)=RSUMs5
                 EXIT
              END IF
              IIs5=IIs5-1
           END DO

           IF (Js5.EQ.1) THEN
              TAO(Is5,Js5)=RSUMs5
           ELSE
              DO IIs5=1,Js5
                 TAO(Is5,IIs5)=XCHMs5(1,1,IIs5)
              END DO
           END IF

	ELSE
           TAO(Is5,Js5)=EXP((tMAX-tMIN)/N0s5*(Js5-0.5)+tMIN) 
	END IF
     END DO
     DO Js5=1,NMOD(Is5)
	TAO(Is5,Js5)=EXP(TAO(Is5,Js5))
     END DO
  END DO S06


  ERRSD=1000.
  S07:DO Ks5=1,GENERATION 	! EVOLUATION OF THE FITTING ENSEMBLE 


     S08:DO Is5=1,SAMPLE			! NORMALIZE THE FITTING ENSEMBLE 
	RESCALE1=1
	RESCALE2=1
	CHECK=1

	JJs5=0
        S12:DO WHILE ((CHECK+RESCALE1+RESCALE2).NE.0)
           JJs5=JJs5+1

           DO Js5=1,NMOD(Is5)
              IF (Js5.EQ.1) THEN
                 BDRs5(1,Js5)=tMIN
              ELSE
                 BDRs5(1,Js5)=BDRs5(2,Js5-1)
              END IF
              IF (Js5.EQ.NMOD(Is5)) THEN
                 BDRs5(2,Js5)=tMAX
              ELSE
                 BDRs5(2,Js5)=(LOG(TAO(Is5,Js5))+LOG(TAO(Is5,Js5+1)))/2.
              END IF
           END DO

           CALL RANDOM_NUMBER(RNDs5)
           POINT=CEILING(RNDs5*(NMOD(Is5)-1)+0.25)

           MIU(Is5,1)=MAX(MIU(Is5,1),RA1s5)
           MIU(Is5,1)=MIN(MIU(Is5,1),RA2s5)
           MIU(Is5,NMOD(Is5))=MAX(MIU(Is5,NMOD(Is5)),RA1s5)
           MIU(Is5,NMOD(Is5))=MIN(MIU(Is5,NMOD(Is5)),RA2s5)

           DO Js5=1,POINT-1
              Ds5=ABS(MIU(Is5,Js5+1)-MIU(Is5,Js5))
              IF (Ds5.GT.CRs5) THEN
                 RESCALE1=1
                 CALL RANDOM_NUMBER(RNDs5)
                 Ds5=CRs5*(RNDs5*2.-1.)
                 MIU(Is5,Js5+1)=Ds5+MIU(Is5,Js5)
                 MIU(Is5,Js5+1)=MAX(MIU(Is5,Js5+1),RA1s5)
                 MIU(Is5,Js5+1)=MIN(MIU(Is5,Js5+1),RA2s5)
              END IF
           END DO

           DO Js5=1,NMOD(Is5)-POINT
              Ds5=ABS(MIU(Is5,NMOD(Is5)+1-Js5)-MIU(Is5,NMOD(Is5)-Js5))
              IF (Ds5.GT.CRs5) THEN
                 RESCALE1=1
                 CALL RANDOM_NUMBER(RNDs5)
                 Ds5=CRs5*(RNDs5*2.-1.)
                 MIU(Is5,NMOD(Is5)-Js5)=Ds5+MIU(Is5,NMOD(Is5)-Js5)
                 MIU(Is5,NMOD(Is5)-Js5)=MAX(MIU(Is5,NMOD(Is5)-Js5),RA1s5)
                 MIU(Is5,NMOD(Is5)-Js5)=MIN(MIU(Is5,NMOD(Is5)-Js5),RA2s5)
              END IF
           END DO
           CHECK=0		

           RSUMs5=0.
           DO Js5=1,NMOD(Is5)
              !	IF (MIU(Is5,Js5)-T(1)/TAO(Is5,Js5).GT.-20.) THEN   ! *
              RSUMs5=RSUMs5+EXP(MIU(Is5,Js5))  ! *
              !	END IF							 ! *
           END DO

           IF (ABS(RSUMs5-1.).LT.CNSTRs5) THEN	   ! *
              RESCALE1=0
           ELSE
              RESCALE2=1
              R1s5=LOG(1./RSUMs5)	   ! *

              DO Js5=1,NMOD(Is5)
                 !	R2s5=T(1)/TAO(Is5,Js5) ! *
                 !	CALL RANDOM_NUMBER(RNDs5)  ! *
                 !	TAO(Is5,Js5)=EXP((1.-RNDs5)*BDRs5(1,Js5)+RNDs5*BDRs5(2,Js5)) ! *

                 !	R2s5=R2s5-T(1)/TAO(Is5,Js5)	! *
                 MIU(Is5,Js5)=MIU(Is5,Js5)+R1s5	 ! *
                 MIU(Is5,Js5)=MIN(MIU(Is5,Js5),RA2s5)
                 MIU(Is5,Js5)=MAX(MIU(Is5,Js5),RA1s5)
              END DO
           END IF

           RSUMs5=0.
           DO Js5=1,NMOD(Is5)
              IF (MIU(Is5,Js5)-T(NUMG)/TAO(Is5,Js5).GT.-20.) THEN
                 RSUMs5=RSUMs5+EXP(MIU(Is5,Js5)-T(NUMG)/TAO(Is5,Js5))
              END IF
           END DO

           DO WHILE (RSUMs5.EQ.0.)
              TAO(Is5,NMOD(Is5))=EXP((LOG(TAO(Is5,NMOD(Is5)))+tMAX)/2.)
              RSUMs5=0.
              DO Js5=1,NMOD(Is5)
                 IF (MIU(Is5,Js5)-T(NUMG)/TAO(Is5,Js5).GT.-20.) THEN
                    RSUMs5=RSUMs5+EXP(MIU(Is5,Js5)-T(NUMG)/TAO(Is5,Js5))
                 END IF
              END DO
           END DO

           IF (ABS(RSUMs5-G(NUMG))/G(NUMG).LT.CNSTRs5) THEN
              RESCALE2=0	

           ELSE
              R1s5=LOG(G(NUMG)/RSUMs5)	

              DO Js5=1,NMOD(Is5)
                 IF (MIU(Is5,Js5)-T(NUMG)/TAO(Is5,Js5).GT.-20.) THEN
                    IND1s5=MAX(RA1s5-MIU(Is5,Js5),MIN(R1s5,0.))  ! *	
                    IND2s5=MIN(RA2s5-MIU(Is5,Js5),MAX(R1s5,0.))	! *
                    CALL RANDOM_NUMBER(RNDs5)	
                    R2s5=(1.-RNDs5)*IND1s5+RNDs5*IND2s5
                    MIU(Is5,Js5)=MIU(Is5,Js5)+R2s5	
                    TAO(Is5,Js5)=1./(1./TAO(Is5,Js5)-(R1s5-R2s5)/T(NUMG)) 
                    TAO(Is5,Js5)=MIN(TAO(Is5,Js5),EXP(BDRs5(2,Js5)))
                    TAO(Is5,Js5)=MAX(TAO(Is5,Js5),EXP(BDRs5(1,Js5)))
                    CHECK=1	
                 END IF
              END DO
           END IF

           !	WRITE(*,*) Ks5,Is5,JJs5,'  2: ',RSUMs5
           !	IF (JJs5.GT.100) THEN
           !	WRITE(31,*)	R1s5
           !	DO Js5=1,NMOD(Is5)
           !	WRITE(31,*) EXP(MIU(Is5,Js5)),TAO(Is5,Js5),EXP(BDRs5(1,Js5)),EXP(BDRs5(2,Js5))
           !	END DO
           !	CLOSE(31)
           !	END IF
           !	PAUSE

	END DO S12

     END DO S08

     DO Is5=1,SAMPLE		! CALCULATE THE FITNESS 
	ERR(Is5)=0.
	DO Js5=1,NUMG
           GF(Js5)=0.
           DO IIs5=1,NMOD(Is5)
              R1s5=MIU(Is5,IIs5)-T(Js5)/TAO(Is5,IIs5)
              IF (R1s5.GT.-20.) THEN
                 GF(Js5)=GF(Js5)+EXP(MIU(Is5,IIs5)-T(Js5)/TAO(Is5,IIs5))
              END IF
           END DO
           IF (Ks5.EQ.1) THEN
              ERR(Is5)=ERR(Is5)+ABS(G(Js5)-GF(Js5))/G(Js5)
           ELSE
              ERR(Is5)=ERR(Is5)+(G(Js5)-GF(Js5))**2./G(Js5)**2./OPT*NUMG
           END IF
	END DO
     END DO

     DO Is5=1,SAMPLE-1		  ! REORDER THE ENSEMBLE BASED ON THE FITNESS 
	MINs5=ERR(Is5)
	Ms5=Is5
	DO Js5=Is5,SAMPLE
           IF (ERR(Js5).LT.MINs5) THEN
              MINs5=ERR(Js5)
              ERR(Js5)=ERR(Is5)
              ERR(Is5)=MINs5
              Ms5=Js5
           END IF
	END DO
	XCHs5=NMOD(Is5)
	DO Js5=1,XCHs5
           XCHMs5(1,1,Js5)=MIU(Is5,Js5)
           XCHMs5(1,2,Js5)=TAO(Is5,Js5)
	END DO
 NMOD(Is5)=NMOD(Ms5)
	DO Js5=1,NMOD(Is5)
           MIU(Is5,Js5)=MIU(Ms5,Js5)
           TAO(Is5,Js5)=TAO(Ms5,Js5)
	END DO
 NMOD(Ms5)=XCHs5
	DO Js5=1,NMOD(Ms5)
           MIU(Ms5,Js5)=XCHMs5(1,1,Js5)
           TAO(Ms5,Js5)=XCHMs5(1,2,Js5)
	END DO
     END DO

     FLAGs5=1
     IIs5=1
     DO WHILE (FLAGs5.EQ.1)	     ! KEEP THE SEEDS
	FLAGs5=0
	DO Is5=1,NSEED
           IF (ERR(IIs5).LE.ERRSD(Is5)) THEN

              NMSD(Is5)=NMOD(IIs5)
              ERRSD(Is5)=ERR(IIs5)
              DO Js5=1,NMSD(Is5)
                 SEEDs5(Is5,1,Js5)=MIU(IIs5,Js5)
                 SEEDs5(Is5,2,Js5)=TAO(IIs5,Js5)
              END DO
              FLAGs5=1
              EXIT

           END IF
	END DO
	IIs5=IIs5+1
     END DO

     OPT=ERRSD(1)
     No5=NMSD(1)
     DO Is5=1,No5
	Mo5(Is5)=SEEDs5(1,1,Is5)
	To5(Is5)=SEEDs5(1,2,Is5)
     END DO

     S15:DO Is5=(SAMPLE/2)+1,SAMPLE		! SELECTION
 	CALL RANDOM_NUMBER(RNDs5)
        NMOD(Is5)=CEILING(RNDs5*(N0s5-2)+1.25)
	DO Js5=1,NMOD(Is5)     
           CALL RANDOM_NUMBER(RNDs5)
           MIU(Is5,Js5)=RNDs5*RA2s5+(1-RNDs5)*RA1s5
           IF (FLAGi5.EQ.0) THEN

04            CALL RANDOM_NUMBER(RNDs5)
              RSUMs5=RNDs5*tMAX+(1-RNDs5)*tMIN
              XCHMs5=0.
              DO IIs5=1,Js5-1
                 XCHMs5(1,1,IIs5)=TAO(Is5,IIs5)
              END DO

              IIs5=Js5-1
              DO WHILE (IIs5.GT.0)
                 IF (RSUMs5.EQ.XCHMs5(1,1,IIs5)) THEN
                    GOTO 04
                 END IF
                 IF (RSUMs5.LT.XCHMs5(1,1,IIs5)) THEN
                    XCHMs5(1,1,IIs5+1)=XCHMs5(1,1,IIs5)
                    IF (IIs5.EQ.1) THEN
                       XCHMs5(1,1,1)=RSUMs5
                    END IF
                 END IF
                 IF (RSUMs5.GT.XCHMs5(1,1,IIs5)) THEN
                    XCHMs5(1,1,IIs5+1)=RSUMs5
                    EXIT
                 END IF
                 IIs5=IIs5-1
              END DO

              IF (Js5.EQ.1) THEN
                 TAO(Is5,Js5)=RSUMs5
              ELSE
                 DO IIs5=1,Js5
                    TAO(Is5,IIs5)=XCHMs5(1,1,IIs5)
                 END DO
              END IF

           ELSE
              TAO(Is5,Js5)=EXP((tMAX-tMIN)/N0s5*(Js5-0.5)+tMIN) 
           END IF
	END DO
	DO Js5=1,NMOD(Is5)
           TAO(Is5,Js5)=EXP(TAO(Is5,Js5))
	END DO
     END DO S15

     IIs5=NSEED
     DO WHILE (IIs5.GT.0)
  	CALL RANDOM_NUMBER(RNDs5)
        Ms5=CEILING(RNDs5*(SAMPLE/2.-1.)+0.25)+SAMPLE/2
        NMOD(Ms5)=NMSD(IIs5)
	DO Js5=1,NMOD(Ms5)
           MIU(Ms5,Js5)=SEEDs5(IIs5,1,Js5)
           TAO(Ms5,Js5)=SEEDs5(IIs5,2,Js5)
	END DO
	IIs5=IIs5-1
     END DO

     S09:DO Is5=1,SAMPLE/2 	    ! CROSS-OVER 
	CALL RANDOM_NUMBER(RNDs5)
	IF (RNDs5.LE.PC) THEN

           DO IIs5=1,2
              LOCs5(IIs5)=2*Is5-2+IIs5
           END DO
           IF (NMOD(LOCs5(1)).LT.NMOD(LOCs5(2))) THEN
              Ms5=1
           ELSE
              Ms5=2
           END IF

           CALL RANDOM_NUMBER(RNDs5)
           PT1(Ms5)=CEILING(RNDs5*(NMOD(LOCs5(Ms5))-1.)+0.25)
           CALL RANDOM_NUMBER(RNDs5)
           PT2(Ms5)=CEILING(RNDs5*(NMOD(LOCs5(Ms5))-1.)+0.25)
           IF (PT1(Ms5).GT.PT2(Ms5)) THEN
              XCHs5=PT1(Ms5)
              PT1(Ms5)=PT2(Ms5)
              PT2(Ms5)=XCHs5
           END IF

           IF (PT1(Ms5).EQ.1) THEN
              IND1s5=tMIN
           ELSE
              IND1s5=(LOG(TAO(LOCs5(Ms5),PT1(Ms5)-1))+LOG(TAO(LOCs5(Ms5),PT1(Ms5))))/2.
           END IF
           IF (PT2(Ms5).EQ.NMOD(LOCs5(Ms5))) THEN
              IND2s5=tMAX
           ELSE
              IND2s5=(LOG(TAO(LOCs5(Ms5),PT2(Ms5)))+LOG(TAO(LOCs5(Ms5),PT2(Ms5)+1)))/2.
           END IF

           FLAGs5=0
           DO Js5=1,NMOD(LOCs5(3-Ms5))
              IF ((TAO(LOCs5(3-Ms5),Js5).GE.EXP(IND1s5)).AND.(FLAGs5.EQ.0)) THEN
                 PT1(3-Ms5)=Js5
                 PT2(3-Ms5)=NMOD(LOCs5(3-Ms5))
                 FLAGs5=1
              END IF
              IF ((TAO(LOCs5(3-Ms5),Js5).GT.EXP(IND2s5)).AND.(FLAGs5.EQ.1)) THEN
                 PT2(3-Ms5)=Js5-1
                 EXIT
              END IF
           END DO
           IF ((PT1(3-Ms5).GT.PT2(3-Ms5)).OR.(FLAGs5.EQ.0)) THEN
              GOTO 03
           END IF

           DO IIs5=1,2
              XCHs5=PT2(3-IIs5)-PT1(3-IIs5)-PT2(IIs5)+PT1(IIs5)
              IF (NMOD(LOCs5(IIs5))+XCHs5.EQ.1) THEN
                 GOTO 03
              END IF
           END DO

           DO IIs5=1,2
              DO Js5=1,NMOD(LOCs5(IIs5))
                 XCHMs5(IIs5,1,Js5)=MIU(LOCs5(IIs5),Js5)
                 XCHMs5(IIs5,2,Js5)=TAO(LOCs5(IIs5),Js5)
              END DO
           END DO

           S13:DO IIs5=1,2
              JJs5=0	
              DO Js5=1,PT1(IIs5)-1
                 JJs5=JJs5+1
                 MIU(LOCs5(IIs5),JJs5)=XCHMs5(IIs5,1,Js5)
                 TAO(LOCs5(IIs5),JJs5)=XCHMs5(IIs5,2,Js5)
              END DO

              DO Js5=PT1(3-IIs5),PT2(3-IIs5)
                 JJs5=JJs5+1
                 MIU(LOCs5(IIs5),JJs5)=XCHMs5(3-IIs5,1,Js5)
                 TAO(LOCs5(IIs5),JJs5)=XCHMs5(3-IIs5,2,Js5)
              END DO

              DO Js5=PT2(IIs5)+1,NMOD(LOCs5(IIs5))
                 JJs5=JJs5+1
                 MIU(LOCs5(IIs5),JJs5)=XCHMs5(IIs5,1,Js5)
                 TAO(LOCs5(IIs5),JJs5)=XCHMs5(IIs5,2,Js5)
              END DO
              NMOD(LOCs5(IIs5))=JJs5
           END DO S13

03	END IF

     END DO S09

     S10:DO Is5=1,SAMPLE		  ! MUTATION

	DO Js5=1,NMOD(Is5)
           IF (Js5.EQ.1) THEN
              BDRs5(1,Js5)=tMIN
           ELSE
              BDRs5(1,Js5)=BDRs5(2,Js5-1)
           END IF
           IF (Js5.EQ.NMOD(Is5)) THEN
              BDRs5(2,Js5)=tMAX
           ELSE
              BDRs5(2,Js5)=(LOG(TAO(Is5,Js5))+LOG(TAO(Is5,Js5+1)))/2.
           END IF
	END DO

	Js5=1
        S14:DO WHILE (Js5.LE.NMOD(Is5))
           CALL RANDOM_NUMBER(RNDs5)
           IF (RNDs5.LE.PM) THEN

              FLAGs5=0
              IF (NMOD(Is5).LE.2) THEN
                 FLAGs5=1
              END IF
              IF (NMOD(Is5).EQ.N0s5) THEN
                 FLAGs5=-1
              END IF
              CALL RANDOM_NUMBER(RNDs5)
              IF ((RNDs5.LE.0.33).AND.(NMOD(Is5).GT.2).AND.(Js5.LT.NMOD(Is5))) THEN
                 FLAGs5=-1
              END IF
              IF ((RNDs5.GT.0.66).AND.(NMOD(Is5).LT.N0s5)) THEN
                 FLAGs5=1
              END IF

              IF (FLAGs5.EQ.0) THEN
                 CALL RANDOM_NUMBER(RNDs5)
                 MIU(Is5,Js5)=RNDs5*RA2s5+(1-RNDs5)*RA1s5
                 IF (FLAGi5.EQ.0) THEN
                    CALL RANDOM_NUMBER(RNDs5)
                    TAO(Is5,Js5)=EXP((BDRs5(2,Js5)-BDRs5(1,Js5))*RNDs5+BDRs5(1,Js5))
                 END IF
              END IF

              IF (FLAGs5.EQ.-1) THEN
                 MIU(Is5,Js5)=MIU(Is5,Js5)+MIU(Is5,Js5+1)
                 TAO(Is5,Js5)=EXP(BDRs5(2,Js5))
                 DO JJs5=Js5+1,NMOD(Is5)-1
                    MIU(Is5,JJs5)=MIU(Is5,JJs5+1)
                    TAO(Is5,JJs5)=TAO(Is5,JJs5+1)
                 END DO
                 DO JJs5=Js5,NMOD(Is5)-1
                    BDRs5(2,JJs5)=BDRs5(2,JJs5+1)
                 END DO
                 DO JJs5=Js5+1,NMOD(Is5)-1
                    BDRs5(1,JJs5)=BDRs5(1,JJs5+1)
                 END DO
                 NMOD(Is5)=NMOD(Is5)-1
              END IF

              IF (FLAGs5.EQ.1) THEN
                 R1s5=MIU(Is5,Js5)
                 R2s5=TAO(Is5,Js5)
                 JJs5=NMOD(Is5)
                 DO WHILE (JJs5.GE.Js5+1)
                    MIU(Is5,JJs5+1)=MIU(Is5,JJs5)
                    TAO(Is5,JJs5+1)=TAO(Is5,JJs5)
                    BDRs5(1,JJs5+1)=BDRs5(1,JJs5)
                    BDRs5(2,JJs5+1)=BDRs5(2,JJs5)
                    JJs5=JJs5-1
                 END DO
                 MIU(Is5,Js5)=R1s5/2.
                 MIU(Is5,Js5+1)=R1s5/2.
                 TAO(Is5,Js5)=SQRT(R2s5*EXP(BDRs5(1,Js5)))
                 TAO(Is5,Js5+1)=SQRT(R2s5*EXP(BDRs5(2,Js5)))
                 BDRs5(1,Js5+1)=LOG(R2s5)
                 BDRs5(2,Js5+1)=BDRs5(2,Js5)
                 BDRs5(2,Js5)=LOG(R2s5)
                 NMOD(Is5)=NMOD(Is5)+1
                 Js5=Js5+1
              END IF

           END IF
           Js5=Js5+1
	END DO	S14
     END DO	S10

     OPT=0.							! *
     DO Is5=1,NUMG					! *
	GF(Is5)=0.						! *
	DO Js5=1,No5					! *
           R1s5=Mo5(Js5)-T(Is5)/To5(Js5)	! *
           IF (R1s5.GT.-20.) THEN			! *
              GF(Is5)=GF(Is5)+EXP(Mo5(Js5)-T(Is5)/To5(Js5)) ! *
           END IF							! *
	END DO							! *
	OPT=OPT+ABS(G(Is5)-GF(Is5))/G(Is5)	 ! *
     END DO							! *

     IF (MOD(Ks5,GENERATION).EQ.0) THEN
        !	IF (MOD(Ks5,100).EQ.0) THEN
	WRITE(*,*) Ks5,' GENERATION FINISH',OPT/NUMG
        write(*,*) 'nmod = ',nmod(1)
     END IF

     IF (Ks5.GT.1) THEN
	IF ((OPT0-OPT)/OPT0.GT.0.01)	THEN
           COUNT=0
           PM=0.3
           PC=0.5
	ELSE
           COUNT=COUNT+1
	END IF
	IF (COUNT.GT.500) THEN
           PM=0.1
           PC=1.0
	END IF
     END IF
     OPT0=OPT

  END DO S07

  ERRMo5=0.
  DO Is5=1,NUMG					! CALCULATE THE FITTING
     GF(Is5)=0.
     DO Js5=1,No5
	R1s5=Mo5(Js5)-T(Is5)/To5(Js5)
	IF (R1s5.GT.-20.) THEN
           GF(Is5)=GF(Is5)+EXP(Mo5(Js5)-T(Is5)/To5(Js5))
	END IF
     END DO
     ERRMo5=ERRMo5+ABS(G(Is5)-GF(Is5))/G(Is5)
  END DO

  ERRMo5=ERRMo5/NUMG

  RETURN 

END SUBROUTINE GEAM

!==================================================================
!          SUBROUTINE 6: CALCULATE FITTING ERROR 
!==================================================================

SUBROUTINE CFER(WLi6,WMi6,WHi6,FERo6,PERo6)

  USE SHARED

  REAL(8) :: MERs6(6)
  REAL(8), INTENT(OUT) :: FERo6(4),PERo6(5)
  INTEGER :: Is6,Js6,SERs6(6),FLAGs6
  INTEGER, INTENT(IN) :: WLi6,WMi6,WHi6	

  SERs6=0					! THE NUMBER OF REGIONS WHOSE FITTING ERROR IS LESS THAN 0.1
  FERo6=0.				! THE AVERGAE ERROR OF DIFFERENT REGIONS
  PERo6=0.				! THE ERROR OF SPECIFIC POINTS
  MERs6=0.				! THE OVERALL ERROR

  DO Is6=1,WLi6			! THE AVERAGE FITTING ERROR FOR LOW FREQUENCY REGION
     FERo6(1)=FERo6(1)+0.25*(LOG(GGF(Is6)/GGE(Is6))+0.75*LOG(GGGF(Is6)/GGGE(Is6)))
  END DO
  FERo6(1)=FERo6(1)/WLi6
  IF (ABS(FERo6(1)).LT.0.1) THEN
     SERs6(1)=SERs6(1)+1
  ELSE
     MERs6(1)=MERs6(1)+ABS(FERo6(1))
  END IF

  DO Is6=WLi6+1,WMi6		! THE AVERAGE FITTING ERROR FOR INTERMEDIATE LOW REGION
     FERo6(2)=FERo6(2)+(LOG(GGF(Is6)/GGE(Is6))+LOG(GGGF(Is6)/GGGE(Is6)))
  END DO
  FERo6(2)=FERo6(2)/2./(WMi6-WLi6)
  IF (ABS(FERo6(2)).LT.0.1) THEN
     SERs6(1)=SERs6(1)+1
  ELSE
     MERs6(1)=MERs6(1)+ABS(FERo6(2))
  END IF

  AERR=FERo6(2)

  DO Is6=WMi6+1,WHi6		! THE AVERAGE FITTING ERROR FOR INTERMEDIATE HIGH REGION
     FERo6(3)=FERo6(3)+(LOG(GGF(Is6)/GGE(Is6))+LOG(GGGF(Is6)/GGGE(Is6)))
  END DO
  FERo6(3)=FERo6(3)/2./(WHi6-WMi6)
  IF (FLAGMR.EQ.0) THEN
     IF (ABS(FERo6(3)).LT.0.1) THEN
	SERs6(1)=SERs6(1)+1
     ELSE
	MERs6(1)=MERs6(1)+ABS(FERo6(3))
     END IF
  END IF

  IF (FLAGMR.EQ.0) THEN				
     DO Is6=WHi6+1,NUMIN		! THE AVERAGE FITTING ERROR FOR HIGH FREQUENCY REGION
	FERo6(4)=FERo6(4)+(LOG(GGF(Is6)/GGE(Is6))+LOG(GGGF(Is6)/GGGE(Is6)))
     END DO
     FERo6(4)=FERo6(4)/2./(NUMIN-WHi6)
  END IF

  PERo6(1)=LOG(GGFC1/GGC1)+LOG(GGGFC1/GGGC1)    
  PERo6(1)=PERo6(1)/2.	! MAGITUDE DIFFERNECE OF 1st CROSSOVER POINT 
  IF (ABS(PERo6(1)).LT.0.1) THEN
     SERs6(2)=SERs6(2)+1
  ELSE
     MERs6(2)=MERs6(2)+ABS(PERo6(1))
  END IF

  PERo6(2)=LOG(FQFC1/FQC1) ! FREQUENCY DIFFERNECE OF 1st CROSSOVER POINT

  PERo6(3)=LOG(GGGFM/GGGM*GGM/GGFM)	! MAGITUDE DIFFERNECE OF INTERMEDIATE POINT
  IF (ABS(PERo6(3)).LT.0.1)	THEN
     SERs6(2)=SERs6(2)+1
  ELSE
     MERs6(2)=MERs6(2)+ABS(PERo6(3))
  END IF

  PERo6(4)=LOG(FQFM/FQM)	! FREQUENCY DIFFERNECE OF INTERMEDIATE POINT
  IF (ABS(PERo6(4)).LT.0.1)	THEN
     SERs6(4)=SERs6(4)+1
  ELSE
     MERs6(4)=MERs6(4)+ABS(PERo6(4))
  END IF

  PERo6(5)=LOG(FQFC2/FQC2)	! FREQUENCY DIFFERENCE OF 2nd CROSSOVER POINT 

  IF (SERs6(1).GT.SERR(1)) THEN
     FLAGPR=1
     SERR(1)=SERs6(1)
     SERR(2)=SERs6(2)
     MERR(1)=MERs6(1)
     MERR(2)=MERS6(2)
  END IF
  IF ((SERs6(1).EQ.SERR(1)).AND.(SERs6(2).GT.SERR(2))) THEN
     FLAGPR=1
     SERR(2)=SERs6(2)
     MERR(1)=MERs6(1)
     MERR(2)=MERs6(2)
  END IF
  IF ((SERs6(1).EQ.SERR(1)).AND.(SERs6(2).EQ.SERR(2))) THEN
     IF (MERs6(1).LT.MERR(1)) THEN
	FLAGPR=1
        MERR(1)=MERs6(1)
        MERR(2)=MERs6(2)
     END IF
     IF ((MERs6(1).EQ.MERR(1)).AND.(MERs6(2).LT.MERR(2))) THEN
	FLAGPR=1
	MERR(2)=MERs6(2)
     END IF
  END IF

  IF (FLAGMR.EQ.1) THEN
     Js6=1
  ELSE
     Js6=0
  END IF

  FLAGs6=1
  S14:IF (FLAGMGT.EQ.0) THEN
     DO Is6=1,2
	IF (ABS(FERo6(Is6)).GE.0.1) THEN
           FLAGs6=0
           GOTO 04
	END IF
     END DO
     DO Is6=1,5-Js6
	IF (ABS(PERo6(Is6)).GE.0.1) THEN
           FLAGs6=0
           GOTO 04
	END IF
     END DO
     IF (ABS(VERR).GE.0.1) THEN
	FLAGs6=0
	GOTO 04
     END IF

  ELSE
     DO Is6=1,2
	IF (ABS(FERo6(Is6)).GE.0.1) THEN
           FLAGs6=0
           GOTO 04
	END IF
     END DO
     DO Is6=1,5
	IF ((Is6.NE.3).AND.(ABS(PERo6(Is6)).GE.0.1)) THEN
           FLAGs6=0
           GOTO 04
	END IF
     END DO
     IF (ABS(VERR).GE.0.1) THEN
	FLAGs6=0
	GOTO 04
     END IF
  END IF S14

04 IF (FLAGs6.EQ.1) THEN
      IF (FLAGN.EQ.1) THEN	
         FLAGD=1
      ELSE
         FLAGN=1
         SERR(1)=0.
         SERR(2)=0.
         MERR(1)=10.
         MERR(2)=10.
         FLAGIT=0
         WRITE(*,*) 'INCREASE ENSEMBLE SIZE'
         WRITE(24,*) 'INCREASE ENSEMBLE SIZE'
      END IF
   END IF

  RETURN 

END SUBROUTINE CFER

!==================================================================
!      SUBROUTINE 7: CALCULATE FITTING ERROR 
!==================================================================

SUBROUTINE PRDF(G0i7,RATIOi7,ZEi7,Lpi7,G0o7,RATIOo7,ZEo7,Lpo7,FLAGo7)

  USE SHARED

  REAL(8) :: C1s7,C2s7,APHo7,APHi7,APH1s7
  REAL(8), INTENT(IN) :: G0i7,RATIOi7,ZEi7,Lpi7
  REAL(8), INTENT(OUT) :: G0o7,RATIOo7,ZEo7,Lpo7	  
  INTEGER :: Is7
  INTEGER, INTENT(OUT) :: FLAGo7

  GERR=0.	          ! THE FITTING DEVIATION FOR Z AND RATIO	          
  FLAGo7=1

  IF (FLAGMGT.EQ.0) THEN
     GERR(1)=(GGGFM/GGFM-GGGM/GGM)/(GGGM/GGM)
     GERR(2)=(FQFC1-FQC1)/FQC1
  ELSE
     GERR(1)=(FQFM/FQFC1-FQM/FQC1)/(FQM/FQC1)
     GERR(2)=(GGGFC1-GGGC1)/GGGC1
  END IF

  IF ((ABS(GERR(1)).LT.0.1).AND.(ABS(GERR(2)).LT.0.1)) THEN
     FLAGo7=0
  END IF
  IF (II.GT.6) THEN
     FLAGo7=0
  END IF

  C1s7=9.75*Kb*TEM/Lpi7**3./G0i7
  C2s7=3.*28./5./PI*Kb*TEM*PHI/d**2./G0i7/Lpi7
  APHi7=2.	
  APH1s7=APHi7					
  DO Is7=1,50				  
     IF ((C1s7*APHi7**2.2+C2s7-3.*APHi7).GT.(APHMA**4.)) THEN
	APHi7=APHMA
     ELSE
 	IF ((C1s7*APHi7**2.2+C2s7-3.*APHi7).LT.(APHMI**4.)) THEN
           APHi7=APHMI
	ELSE
           APHi7=(C1s7*APHi7**2.2+C2s7-3.*APHi7)**0.25
	END IF
     END IF
     APHi7=(APHi7+APH1s7)/2.
     APH1s7=APHi7
  END DO

  IF (FLAGMGT.EQ.0) THEN
     ZEo7=ZEi7*(GGM/GGFM*GGGFM/GGGM)**RELAX
  ELSE
     ZEo7=(ZEi7-2.)*(FQM/FQFM*FQFC1/FQC1)**(0.8*RELAX)+2.
  END IF
  IF (ZEo7.GT.ZEMA) THEN
     ZEo7=ZEMA
  END IF
  IF (ZEo7.LT.ZEMI) THEN
     ZEo7=ZEMI
  END IF
  ZEo7=ZEo7**RELAX*ZEi7**(1.-RELAX)

  S11:IF (SLP.EQ.1) THEN 
     IF (FLAGo7.EQ.0) THEN

	IF (FLAGMR.EQ.0) THEN
           Lpo7=Lpi7*(FQFC2/FQC2)**(0.33*RELAX)
	ELSE
           Lpo7=Lpi7*EXP(-0.8*FRLP*RELAX)
	END IF
        IF (Lpo7.GT.LpMA) THEN
           Lpo7=LpMA
	END IF
	IF (Lpo7.LT.LpMI) THEN
           Lpo7=LpMI
	END IF
	Lpo7=Lpo7**RELAX*Lpi7**(1.-RELAX)

     ELSE
	Lpo7=Lpi7
     END IF
  END IF S11

  IF ((FLAGo7.EQ.0).OR.(FLAGMGT.EQ.1)) THEN
     G0o7=G0i7*(GGGC1/GGGFC1)**RELAX
     IF ((ABS(LOG(GGGC1/GGGFC1)).LT.0.1).AND.(FLAGMGT.EQ.0)) THEN
	G0o7=G0i7/EXP(VERR*(1.-RELAX))
     END IF
  ELSE
     G0o7=G0i7
  END IF
  IF (G0o7.GT.G0MA) THEN
     G0o7=G0MA
  END IF
  IF (G0o7.LT.G0MI) THEN
     G0o7=G0MI
  END IF
  G0o7=G0o7**RELAX*G0i7**(1.-RELAX)

  C1s7=9.75*Kb*TEM/Lpo7**3./G0o7
  C2s7=3.*28./5./PI*Kb*TEM*PHI/d**2./G0o7/Lpo7
  APHo7=APHi7					
  APH1s7=APHo7				
  DO Is7=1,50				  
     IF ((C1s7*APHo7**2.2+C2s7-3.*APHo7).GT.(APHMA**4.)) THEN
	APHo7=APHMA
     ELSE
 	IF ((C1s7*APHo7**2.2+C2s7-3.*APHo7).LT.(APHMI**4.)) THEN
           APHo7=APHMI
	ELSE
           APHo7=(C1s7*APHo7**2.2+C2s7-3.*APHo7)**0.25
	END IF
     END IF
     APHo7=(APHo7+APH1s7)/2.
     APH1s7=APHo7
  END DO

  IF ((FLAGo7.EQ.0).OR.(FLAGMGT.EQ.0)) THEN
     RATIOo7=RATIOi7*(FQFC1/FQC1*(Lpi7*ZEi7/Lpo7/ZEo7)**3.)**(1.5*RELAX)
     IF ((ABS(LOG(FQFC1/FQC1)).LT.0.1).AND.(FLAGMGT.EQ.1)) THEN
	RATIOo7=RATIOi7/EXP(3.5*VERR*(1.-RELAX))
     END IF
  ELSE
     RATIOo7=RATIOi7
  END IF
  IF (RATIOo7.GT.RATIOMA) THEN
     RATIOo7=RATIOMA
  END IF
  IF (RATIOo7.LT.RATIOMI) THEN
     RATIOo7=RATIOMI
  END IF
  RATIOo7=RATIOo7**RELAX*RATIOi7**(1.-RELAX)

  II=II+1

  RETURN 

END SUBROUTINE PRDF
