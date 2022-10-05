      program res_to_dat
! This program converts a diablo restart file to a format readable by vapor
! A specified number of wavelet transforms are performed to produce
! compressed versions of the data 
! Note, as of version 1.1.2, Vapor cannot handle stretched grids

! grid_def should contain the grid size and number of scalars
      include 'header'

      character*155 FNAME
      character*155 FNAME_TH(1:N_TH)

! Define Variables
      integer i,j,k,n
      real*8 buffer2(0:NX-1,0:NZ-1,1:NY)
      real*8 buffer3(0:NX-1,0:NZ-1,0:NY-1)
      real*8 buffercgrad(0:NX-1,0:NZ-1,0:NY-1,1:N_TH)
      real*8 ux(Nx,NZ),uy(Nx,Nz),uz(Nx,Nz)
      real*8 dumt(1000),minth,maxth
! The sizes of these strings may need to be changed
      character*11 size_str
      character*23 len_str
      character*3 N_TIMES
      character*4 cases
      character*120 savename,stt

! Variables for movie:
      integer TIME1,TIME_START,TIME_END,TIME_INDEX,DUMMY_TIME_STEP
      integer THYES
      real*8 DUMMY_TIME
      TIME_START=100000
      TIME_END=101000
      TIME_STEP=10

! **** User Input *****
! Number of periodic directions used in the simulation
      NUM_PER_DIR=3
! This string should contain the size of the buffer array
      LX=2
      LZ=2

      NXM=NX-1
      NYM=NY-1
      NZM=NZ-1
      
      call init_fft
 
      THYES=1

      TIME_INDEX=TIME_START/100-1
      do TIME1=TIME_START,TIME_END,TIME_STEP 
      write(*,*) 'TIME: ',TIME1
      TIME_INDEX=TIME_INDEX+1

      cases='Lint'

! Name of restart file(s)
         FNAME='../Cases/'//cases//'/diablo.saved.'
     &        //CHAR(MOD(TIME1,10000000)/1000000+48)
     &        //CHAR(MOD(TIME1,1000000)/100000+48)
     &        //CHAR(MOD(TIME1,100000)/10000+48)
     &        //CHAR(MOD(TIME1,10000)/1000+48)
     &        //CHAR(MOD(TIME1,1000)/100+48)
     &        //CHAR(MOD(TIME1,100)/10+48)
     &        //CHAR(MOD(TIME1,10)+48)
         stt='../Cases/'//cases//'/pos.'
     &        //CHAR(MOD(TIME1,10000000)/1000000+48)
     &        //CHAR(MOD(TIME1,1000000)/100000+48)
     &        //CHAR(MOD(TIME1,100000)/10000+48)
     &        //CHAR(MOD(TIME1,10000)/1000+48)
     &        //CHAR(MOD(TIME1,1000)/100+48)
     &        //CHAR(MOD(TIME1,100)/10+48)
     &        //CHAR(MOD(TIME1,10)+48)

      IF (THYES) THEN
         DO N=1,N_TH
      FNAME_TH(N)='../Cases/'//cases//'/diablo_th'
     &        //CHAR(MOD(N,100)/10+48)
     &        //CHAR(MOD(N,10)+48) // '.saved.'
     &        //CHAR(MOD(TIME1,10000000)/1000000+48)
     &        //CHAR(MOD(TIME1,1000000)/100000+48)
     &        //CHAR(MOD(TIME1,100000)/10000+48)
     &        //CHAR(MOD(TIME1,10000)/1000+48)
     &        //CHAR(MOD(TIME1,1000)/100+48)
     &        //CHAR(MOD(TIME1,100)/10+48)
     &        //CHAR(MOD(TIME1,10)+48)
         END DO
      ENDIF
      WRITE(6,*)   'Reading flow from ', FNAME_TH(1)

      OPEN(UNIT=10,FILE=FNAME
     &   ,STATUS="OLD",FORM="UNFORMATTED")
      READ (10) NX_T, NY_T, NZ_T, NUM_PER_DIR_T, DUMMY_TIME
     &        , DUMMY_TIME_STEP

      write(*,*) NX_T, NY_T, NZ_T, NUM_PER_DIR_T, DUMMY_TIME
     &        , DUMMY_TIME_STEP

      IF ((NX .NE. NX_T) .OR. (NY .NE. NY_T) .OR. (NZ .NE. NZ_T)) then
         write(*,*) 'NX,NY,NZ,NX_T,NY_T,NZ_T:',NX,NY,NZ,NX_T,NY_T,NZ_T
         STOP 'Error: old flowfield wrong dimensions. '
      END IF
      IF (NUM_PER_DIR .NE. NUM_PER_DIR_T)
     *     STOP 'Error: old flowfield wrong NUM_PER_DIR. '

!      stop
           
      write(*,*) 'READING FLOW'
      IF (NUM_PER_DIR.EQ.3) THEN
        READ (10) (((CU1(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY),
     *            (((COm(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY),
     *            (((CU3(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY)
! TEMPORARY: ZERO MEAN MODES:
         write(*,*) 'Zeroing mean flow...'
         CU1(0,0,0)=0.d0
         COm(0,0,0)=0.d0
         CU3(0,0,0)=0.d0

        DO J=0,TNKY
        DO K=0,TNKZ
        DO I=0,NKX
           CPsi(I,K,J)=-COm(I,K,J)/(KX2(I)+KY2(J)+KZ2(K))
        END DO
        END DO
        END DO
        CPsi(0,0,0)=0.d0

         DO J=0,TNKY
         DO K=0,TNKZ
         DO I=0,NKX
            CU1(I,K,J)=-CIKZ(K)*CPsi(I,K,J)
            CU3(I,K,J)=CIKX(I)*CPsi(I,K,J)
         END DO
         END DO
         END DO 



       IF (THYES) THEN
        DO N=1,N_TH
          OPEN(UNIT=11,FILE=FNAME_TH(N)
     &           ,STATUS="OLD"
     &           ,FORM="UNFORMATTED")
          WRITE(6,*) 'Reading density from ', FNAME_TH(N)
          READ (11) NX_T, NY_T, NZ_T, NUM_PER_DIR_T, DUMMY_TIME
     &          , DUMMY_TIME_STEP
          READ (11) (((CTH(I,K,J,N)
     &           ,I=0,NKX),K=0,TNKZ),J=0,TNKY)
         CLOSE(11)
        END DO
       ENDIF
      ENDIF
      CLOSE(10)

!      open(unit=12,file=stt,form="unformatted")
!      read(12) xpos,zpos
!      close(12)


      write(*,*) 'Done reading restart file'
      

       IF (NUM_PER_DIR.eq.2) THEN
!       do j=1,NY
!         write(22,*) J,GYF(J),CU1(0,0,J)
!         U1_BAR(j)=CU1(0,0,J)
!         DO N=1,N_TH
!           TH_BAR(J,N)=CTH(0,0,J,N)
!         END DO
!       end do

        call fft_xz_to_physical(CU1,U1,0,NY+1)
        call fft_xz_to_physical(COm,Om,0,NY+1)
        call fft_xz_to_physical(CU3,U3,0,NY+1)
        DO N=1,N_TH
          CALL FFT_XZ_TO_PHYSICAL(CTH(0,0,0,N),TH(0,0,0,N),0,NY+1)
        END DO
        ELSE IF (NUM_PER_DIR.eq.3) THEN
          CALL fft_xzy_to_physical(CU1,U1)
          CALL fft_xzy_to_physical(COm,Om)
          CALL fft_xzy_to_physical(CU3,U3)
          IF (THYES) THEN
          DO N=1,N_TH

            CALL fft_xzy_to_physical(CTH(0,0,0,N),TH(0,0,0,N))
            minth=TH(0,0,0,N)
            maxth=TH(0,0,0,N)
      	    DO J=0,NYM
            DO K=0,NZM
            DO I=0,NXM
               minth=min(TH(I,K,J,N),minth)
               maxth=max(TH(I,K,J,N),maxth)
            ENDDO
            ENDDO
            ENDDO
            write(*,*) 'extrema ', minth, maxth
          END DO
          ENDIF
        END IF

       savename='../Cases/'//cases//'/bin'
     &               //CHAR(MOD(TIME_INDEX,100000)/10000+48)
     &               //CHAR(MOD(TIME_INDEX,10000)/1000+48)
     &               //CHAR(MOD(TIME_INDEX,1000)/100+48)
     &               //CHAR(MOD(TIME_INDEX,100)/10+48)
     &               //CHAR(MOD(TIME_INDEX,10)+48)
       write(*,*) 'Saving to file: ', savename
       do i=1,NX
         do j=1,1
           do k=1,NZ
             ux(i,k)=U1(i-1,k-1,j-1)
             uy(i,k)=Om(i-1,k-1,j-1)
             uz(i,k)=U3(i-1,k-1,j-1)
           enddo
         enddo
       enddo       
!       DO I=1,NX
!       write(*,*) U1(i,256,0),U2(i,256,0),U3(i,256,0)
!       enddo


       open(22, file=savename, form='unformatted')
       write(22) ux, uz, uy
       close(22)
       savename='../Cases/'//cases//'/pos'
     &               //CHAR(MOD(TIME_INDEX,100000)/10000+48)
     &               //CHAR(MOD(TIME_INDEX,10000)/1000+48)
     &               //CHAR(MOD(TIME_INDEX,1000)/100+48)
     &               //CHAR(MOD(TIME_INDEX,100)/10+48)
     &               //CHAR(MOD(TIME_INDEX,10)+48)
!       open(22, file=savename, form='unformatted')
!       write(22) xpos,zpos
!       close(22)
!       write(*,*) NX,NY,NZ
       IF (THYES) THEN
       do i=1,NX
         do j=1,1
           do k=1,NZ
             ux(i,k)=TH(i-1,k-1,j-1,1)
!             uy(i,k)=TH(i-1,k-1,j-1,2)
!             uz(i,k)=TH(i-1,k-1,j-1,3)
!             write(*,*) i,j,k
           enddo
         enddo
       enddo       



      savename='../Cases/'//cases//'/Th'
     &               //CHAR(MOD(TIME_INDEX,100000)/10000+48)
     &               //CHAR(MOD(TIME_INDEX,10000)/1000+48)
     &               //CHAR(MOD(TIME_INDEX,1000)/100+48)
     &               //CHAR(MOD(TIME_INDEX,100)/10+48)
     &               //CHAR(MOD(TIME_INDEX,10)+48)
       write(*,*) 'Saving to file: ', savename
       open(22, file=savename, form='unformatted')
       write(22) ux!,uy,uz

!       do i=1,NX
!         do j=1,NY
!           do k=1,NZ
!             ux(i,k,j)=buffercgrad(i-1,j-1,k-1,1)
!             uy(i,k,j)=buffercgrad(i-1,j-1,k-1,2)
!           enddo
!         enddo
!       enddo

!       write(22) ux,uy

       close(22)
       ENDIF
      end do

!     We are done
      stop
      end

