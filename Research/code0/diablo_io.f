C This file contains subroutines for inputting and outputting data in
C Diablo as well as all subroutines called directly from diablo.f
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INITIALIZE
      INCLUDE 'header'

      REAL    VERSION, CURRENT_VERSION
      logical RESET_TIME
      INTEGER I, J, K, N

      OPEN (11,file='input.dat',form='formatted',status='old')      

C Read input file.
C   (Note - if you change the following section of code, update the
C    CURRENT_VERSION number to make obsolete previous input files!)


      CURRENT_VERSION=1.0
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*) FLAVOR, FLAVORS, FLAVORI, VERSION
      IF (VERSION .NE. CURRENT_VERSION) STOP 'Wrong input data format.'
      READ(11,*)
      READ(11,*) USE_MPI,    LES
      READ(11,*)
      READ(11,*) NU, LX, LY, LZ
      READ(11,*)
      READ(11,*) NUM_PER_DIR, CREATE_NEW_FLOW
      READ(11,*)
      READ(11,*) N_TIME_STEPS, DELTA_T, RESET_TIME, VARIABLE_DT, CFL
     &            , UPDATE_DT
      READ(11,*)
      READ(11,*) VERBOSITY, SAVE_FLOW_INT, SAVE_STATS_INT, MOVIE
      READ(11,*)
! Read in the parameters for the N_TH scalars
      DO N=1,N_TH
        READ(11,*)
        READ(11,*) CREATE_NEW_TH(N)
        READ(11,*)
        READ(11,*) FILTER_TH(N), FILTER_INT(N)
        READ(11,*)
        READ(11,*) RI_TAU(N), PR(N), REACTION(N)
      END DO

C      IF ((FLAVOR.NE.'Turb').OR.(FLAVOR.NE.'Bick').OR.(FLAVOR.NE.'Gyre')
C     &.OR.(FLAVOR.NE.'NFlo').OR.(FLAVOR.NE.'Beta').OR.(FLAVOR.NE.'Gyr2')
C     &) THEN
C         WRITE(*,*)'WRONG FLAVOR', FLAVOR
C         STOP
C      END IF

C If we are using MPI, then Initialize the MPI Variables
      IF (USE_MPI) THEN
        CALL INIT_MPI
      END IF

C Initialize case-specific packages.
      IF (NUM_PER_DIR.EQ.3) THEN
        CALL INPUT_PER
        CALL CREATE_GRID_PER
        CALL INIT_PER
      ELSEIF (NUM_PER_DIR.EQ.2) THEN
        CALL INPUT_CHAN
        CALL CREATE_GRID_CHAN
        IF (USE_MPI) THEN
          CALL INIT_CHAN_MPI
        ELSE 
          CALL INIT_CHAN
        END IF
      ELSEIF (NUM_PER_DIR.EQ.1) THEN
        CALL INPUT_DUCT
        CALL CREATE_GRID_DUCT
        CALL INIT_DUCT
      ELSEIF (NUM_PER_DIR.EQ.0) THEN
        CALL INPUT_CAV 
        CALL CREATE_GRID_CAV
        CALL INIT_CAV
      END IF

C Initialize grid
      IF (FLAVOR.NE.'Ensemble') THEN
      write(6,'(A17,A8,A18,A6)') ' Flavor of flow: ',FLAVOR, 
     *            'Flavor of scalar: ',FLAVORS
      WRITE(6,'(A17,I3,A7,I3,A7,I3,A1)') ' Grid size: NX = ',NX, 
     *            ', NY = ',NY,', NZ = ',NZ,'.'
      WRITE(6,*) 'Viscosity = ', NU
      DO N=1,N_TH
        WRITE(6,*) 'Scalar number: ',N
        WRITE(6,*) 'Ri=',RI_TAU(N),'Pr=',PR(N)
        IF (FLAVOR .EQ. 'CHEMOTAXIS')THEN
        WRITE(6,*) 'Da/Pe=',C0(N),'Vm=',V_S(N)
        ENDIF
      END DO
      END IF
      NXM=NX-1
      NYM=NY-1
      NZM=NZ-1

C Initialize storage arrays.
      DO K=0,NZ+1
        DO I=0,NX+1 
          DO J=0,NY+1
            U1(I,K,J)=0.
            U3(I,K,J)=0.
            U2(I,K,J)=0.
            P(I,K,J)=0.
            R1(I,K,J)=0.
            R2(I,K,J)=0.
            R3(I,K,J)=0.
            F1(I,K,J)=0.
            F2(I,K,J)=0.
            F3(I,K,J)=0.
          END DO
        END DO
      END DO

C Initialize FFT package (includes defining the wavenumber vectors).
      CALL INIT_FFT
	
C Initialize RKW3 parameters.
      H_BAR(1)=DELTA_T*(8.0/15.0)
      H_BAR(2)=DELTA_T*(2.0/15.0)
      H_BAR(3)=DELTA_T*(5.0/15.0)
      BETA_BAR(1)=1.0
      BETA_BAR(2)=25.0/8.0
      BETA_BAR(3)=9.0/4.0
      ZETA_BAR(1)=0.0
      ZETA_BAR(2)=-17.0/8.0
      ZETA_BAR(3)=-5.0/4.0

C Initialize values for reading of scalars
      NUM_READ_TH=0
      DO N=1,N_TH
        IF (CREATE_NEW_TH(N)) THEN
          NUM_READ_TH=NUM_READ_TH 
        ELSE
          NUM_READ_TH=NUM_READ_TH + 1
          READ_TH_INDEX(NUM_READ_TH)=N
        END IF
      END DO
      IF (NUM_PER_DIR.EQ.2) THEN
        CALL CREATE_TH_CHAN 
      ELSE IF (NUM_PER_DIR.EQ.3) THEN
        CALL CREATE_TH_PER
	  IF (TRUCK) CALL INIT_TRUCK
      END IF 
	
C Hook case-specific initializations
	IF (FLAVOR.EQ.'Ensemble') THEN
	  CALL INIT_ENSEM
	ELSEIF (FLAVOR.EQ.'Adjoint') THEN
	  CALL INIT_ADJOINT
	END IF
C Create flow.
      IF (CREATE_NEW_FLOW) THEN
        IF (NUM_PER_DIR.EQ.3) THEN
           IF (FLAVOR .EQ. 'Turb')THEN
!              CALL CREATE_FLOW_PER
           ENDIF
        ELSEIF (NUM_PER_DIR.EQ.2) THEN
          CALL CREATE_FLOW_CHAN
        ELSEIF (NUM_PER_DIR.EQ.1) THEN
          CALL CREATE_FLOW_DUCT
        ELSEIF (NUM_PER_DIR.EQ.0) THEN
          CALL CREATE_FLOW_CAV
       END IF
        IF (FLAVOR.NE.'Ensemble') THEN
	    write(*,*) 'A new flowfield has been created'
	  END IF
        IF (FLAVOR.EQ.'Basic') THEN
	    CALL SAVE_STATS(.FALSE.)
	    CALL SAVE_FLOW(.FALSE.)
	 END IF
      ELSE
        write(*,*) 'Reading flow...'
        CALL READ_FLOW
        write(*,*) 'Done reading flow'

	 
! TEMPORARY, ZERO MEAN FLOW
      write(*,*),'****Removing mean flow****'
      CU1(0,0,0)=0.d0
      CU2(0,0,0)=0.d0
      CU3(0,0,0)=0.d0
 
C Initialize flow.
      IF (RESET_TIME .OR. CREATE_NEW_FLOW) THEN
        PREVIOUS_TIME_STEP=0
        TIME_STEP=0
        TIME=0
      END IF

        CALL SAVE_STATS(.FALSE.)
        IF (NUM_PER_DIR.EQ.3) THEN
          CALL POISSON_P_PER
        ELSEIF (NUM_PER_DIR.EQ.2) THEN
          CALL POISSON_P_CHAN
        ELSEIF (NUM_PER_DIR.EQ.1) THEN
          CALL POISSON_P_DUCT
        ELSEIF (NUM_PER_DIR.EQ.0) THEN
          CALL POISSON_P_CAV
        END IF

      END IF

      RETURN
      END



C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_STATS(FINAL)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      LOGICAL FINAL

      IF (NUM_PER_DIR.EQ.3) THEN
        CALL SAVE_STATS_PER(FINAL)          
      ELSEIF (NUM_PER_DIR.EQ.2) THEN
        CALL SAVE_STATS_CHAN(FINAL)          
      ELSEIF (NUM_PER_DIR.EQ.1) THEN
        CALL SAVE_STATS_DUCT(FINAL)          
      ELSEIF (NUM_PER_DIR.EQ.0) THEN
        CALL SAVE_STATS_CAV(FINAL)          
      END IF
	
	if (FLAVOR.eq.'Basic') then
        write(*,*) 'done save_stats diablo'
	end if

      IF (FINAL) THEN
        IF (NUM_PER_DIR.EQ.3) THEN
          CALL VIS_FLOW_PER         
        ELSEIF (NUM_PER_DIR.EQ.2) THEN
          CALL VIS_FLOW_CHAN         
        ELSEIF (NUM_PER_DIR.EQ.1) THEN
          CALL VIS_FLOW_DUCT          
        ELSEIF (NUM_PER_DIR.EQ.0) THEN
          CALL VIS_FLOW_CAV         
        END IF
      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE READ_FLOW
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      CHARACTER*35 FNAME
      CHARACTER*35 FNAME_TH(N_TH)
      INTEGER I, J, K, N, NUM_PER_DIR_T

      FNAME='diablo.start'
      DO N=1,N_TH
        FNAME_TH(N)='diablo_th'
     &        //CHAR(MOD(N,100)/10+48)
     &        //CHAR(MOD(N,10)+48) // '.start'
      END DO
      WRITE(6,*)   'Reading flow from ',FNAME

      OPEN(UNIT=10,FILE=FNAME,STATUS="OLD",FORM="UNFORMATTED")
      READ (10) NX_T,NY_T,NZ_T,NUM_PER_DIR_T,TIME,TIME_STEP

!     Between data from 32 bit and 64 bit systems
!      NX_T=NX_T1
!      NY_T=NY_T1
!      NZ_T=NZ_T1
!      NUM_PER_DIR_T=NUM_PER_DIR_T1
!      TIME_STEP=TIME_STEP1

      WRITE(6,'(A16,I3,A7,I3,A7,I3,A1)')'Grid size: NX = ',NX,
     *            ', NY = ',NY,', NZ = ',NZ,'.'
      WRITE(6,*) 'Time Step = ', TIME_STEP

      IF ((NX .NE. NX_T) .OR. (NY .NE. NY_T) .OR. (NZ .NE. NZ_T))
     *     STOP 'Error: old flowfield wrong dimensions. '
      IF (NUM_PER_DIR .NE. NUM_PER_DIR_T)
     *     STOP 'Error: old flowfield wrong NUM_PER_DIR. '

      write(*,*) 'READING FLOW'
      IF (NUM_PER_DIR.EQ.3) THEN
        READ (10) (((CU1(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY),
     *            (((COm(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY),
     *            (((CU3(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY)
        DO N=1,NUM_READ_TH
! Specify in input.dat which scalars are to be read
          OPEN(UNIT=11,FILE=FNAME_TH(READ_TH_INDEX(N)),STATUS="OLD"
     &           ,FORM="UNFORMATTED")
          READ (11) NX_T, NY_T, NZ_T, NUM_PER_DIR_T, TIME, TIME_STEP
          READ (11) (((CTH(I,K,J,READ_TH_INDEX(N))
     &           ,I=0,NKX),K=0,TNKZ),J=0,TNKY)
         CLOSE(11)
        END DO
      ELSEIF (NUM_PER_DIR.EQ.2) THEN
        READ (10) (((CU1(I,K,J),I=0,NKX),K=0,TNKZ),J=1,NY),
     *            (((CU2(I,K,J),I=0,NKX),K=0,TNKZ),J=2,NY),
     *            (((CU3(I,K,J),I=0,NKX),K=0,TNKZ),J=1,NY)

        DO N=1,NUM_READ_TH
! Specify in input.dat which scalars are to be read
          OPEN(UNIT=11,FILE=FNAME_TH(READ_TH_INDEX(N)),STATUS="OLD"
     &           ,FORM="UNFORMATTED")
          READ (11) NX_T, NY_T, NZ_T, NUM_PER_DIR_T, TIME, TIME_STEP
          READ (11) (((CTH(I,K,J,READ_TH_INDEX(N))
     &           ,I=0,NKX),K=0,TNKZ),J=1,NY)
         CLOSE(11)
        END DO
      ELSEIF (NUM_PER_DIR.EQ.1) THEN
        READ (10) (((CU1(I,K,J),I=0,NKX),K=1,NZ),J=1,NY),
     *            (((CU2(I,K,J),I=0,NKX),K=1,NZ),J=2,NY),
     *            (((CU3(I,K,J),I=0,NKX),K=2,NZ),J=1,NY)
      ELSEIF (NUM_PER_DIR.EQ.0) THEN
        READ (10) (((U1(I,K,J),I=2,NX),K=1,NZ),J=1,NY),
     *            (((U2(I,K,J),I=1,NX),K=1,NZ),J=2,NY),
     *            (((U3(I,K,J),I=1,NX),K=2,NZ),J=1,NY)
      END IF
      CLOSE(10)
      CLOSE(11)

! get the initial energy in low wavenumbers
!      EK0=0.d0
!      DO J=0,TNKY
!        DO K=0,TNKZ
!          DO I=0,NKX
!              IF ((SQRT(KX2(I)+KY2(J)+KZ2(K)).le.2.5d0)) THEN
!                EK0=EK0+CU1(I,K,J)*CONJG(CU1(I,K,J))
!     &               +CU2(I,K,J)*CONJG(CU2(I,K,J))
!     &               +CU3(I,K,J)*CONJG(CU3(I,K,J))
!              ENDIF
!          END DO
!        END DO
!      END DO
!      write(*,*) 'EK0: ',EK0

C Apply initial boundary conditions, set ghost cells
      IF (USE_MPI) THEN
        call APPLY_BC_VEL_MPI
      ELSE
        call APPLY_BC_VEL_LOWER
        call APPLY_BC_VEL_UPPER
      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_FLOW(FINAL)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      CHARACTER*65 FNAME
      CHARACTER*65 FNAME_TH(N_TH)
      INTEGER      I, J, K, N
      LOGICAL      FINAL
 
      IF (FINAL) THEN
         FNAME='diablo.res'
         DO N=1,N_TH
           FNAME_TH(N)='diablo_th'
     &        //CHAR(MOD(N,100)/10+48)
     &        //CHAR(MOD(N,10)+48) // '.res'
         END DO
      ELSE
         FNAME='diablo.saved.'
     &        //CHAR(MOD(TIME_STEP,10000000)/1000000+48)
     &        //CHAR(MOD(TIME_STEP,1000000)/100000+48)
     &        //CHAR(MOD(TIME_STEP,100000)/10000+48)
     &        //CHAR(MOD(TIME_STEP,10000)/1000+48)
     &        //CHAR(MOD(TIME_STEP,1000)/100+48)
     &        //CHAR(MOD(TIME_STEP,100)/10+48)
     &        //CHAR(MOD(TIME_STEP,10)+48)
         DO N=1,N_TH
           FNAME_TH(N)='diablo_th'
     &        //CHAR(MOD(N,100)/10+48)
     &        //CHAR(MOD(N,10)+48) // '.saved.'
     &        //CHAR(MOD(TIME_STEP,10000000)/1000000+48)
     &        //CHAR(MOD(TIME_STEP,1000000)/100000+48)
     &        //CHAR(MOD(TIME_STEP,100000)/10000+48)
     &        //CHAR(MOD(TIME_STEP,10000)/1000+48)
     &        //CHAR(MOD(TIME_STEP,1000)/100+48)
     &        //CHAR(MOD(TIME_STEP,100)/10+48)
     &        //CHAR(MOD(TIME_STEP,10)/1+48)
         END DO
      END IF
      WRITE(6,'(A16, A22)') 'Saving to file ',FNAME

      OPEN(UNIT=10,FILE=FNAME,STATUS="UNKNOWN",FORM="UNFORMATTED")
      WRITE(10) NX, NY, NZ, NUM_PER_DIR, TIME, TIME_STEP


      IF (NUM_PER_DIR.EQ.3) THEN
        WRITE(10) (((CU1(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY),
     *            (((COm(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY),
     *            (((CU3(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY)
 
        DO N=1,N_TH
          OPEN(UNIT=11,FILE=FNAME_TH(N),STATUS="UNKNOWN"
     &       ,FORM="UNFORMATTED")
          WRITE(11) NX, NY, NZ, NUM_PER_DIR, TIME, TIME_STEP
          WRITE(11) (((CTH(I,K,J,N),I=0,NKX),K=0,TNKZ),J=0,TNKY)
          CLOSE(11)
        END DO
      ELSEIF (NUM_PER_DIR.EQ.2) THEN
        WRITE(10) (((CU1(I,K,J),I=0,NKX),K=0,TNKZ),J=1,NY),
     *            (((CU2(I,K,J),I=0,NKX),K=0,TNKZ),J=2,NY),
     *            (((CU3(I,K,J),I=0,NKX),K=0,TNKZ),J=1,NY)
        DO N=1,N_TH
          OPEN(UNIT=11,FILE=FNAME_TH(N),STATUS="UNKNOWN"
     &       ,FORM="UNFORMATTED")
          WRITE(11) NX, NY, NZ, NUM_PER_DIR, TIME, TIME_STEP
          WRITE(11) (((CTH(I,K,J,N),I=0,NKX),K=0,TNKZ),J=0,NY)
          CLOSE(11)
        END DO
      ELSEIF (NUM_PER_DIR.EQ.1) THEN
        WRITE(10) (((CU1(I,K,J),I=0,NKX),K=1,NZ),J=1,NY),
     *            (((CU2(I,K,J),I=0,NKX),K=1,NZ),J=2,NY),
     *            (((CU3(I,K,J),I=0,NKX),K=2,NZ),J=1,NY)
      ELSEIF (NUM_PER_DIR.EQ.0) THEN
        WRITE(10) (((U1(I,K,J),I=2,NX),K=1,NZ),J=1,NY),
     *            (((U2(I,K,J),I=1,NX),K=1,NZ),J=2,NY),
     *            (((U3(I,K,J),I=1,NX),K=2,NZ),J=1,NY)
      END IF
      CLOSE(10)
      CLOSE(11)

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FILTER(n)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      integer n

      IF (NUM_PER_DIR.EQ.3) THEN
        CALL FILTER_PER
      ELSEIF (NUM_PER_DIR.EQ.2) THEN
        CALL FILTER_CHAN(n)
      END IF

      RETURN
      END

