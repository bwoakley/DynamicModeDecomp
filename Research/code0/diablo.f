C******************************************************************************|
C diablo.f -> DNS In A Box, Laptop Optimized                       VERSION 0.9
C
C This Fortran 77 code computes incompressible flow in a box.
C
C Primative variables (u,v,w,p) are used, and continuity is enforced with a
C fractional step algorithm.
C
C SPATIAL DERIVATIVES:
C   0, 1, 2, or 3 directions are taken to be periodic and handled spectrally
C   (these cases are referred to as the "periodic", "channel", "duct", and
C    "cavity" cases respectively).
C   The remaining directions are taken to be bounded by walls and handled with
C   momentum- and energy-conserving second-order central finite differences.
C
C TIME ADVANCEMENT
C   Two main approaches are implemented:
C     1. RKW3 on nonlinear terms and CN on viscous terms over each RK substep.
C     2. RKW3 on y-derivative terms and CN on other terms over each RK substep.
C
C The emphasis in this introductory code is on code simplicity:
C   -> All variables are in core.
C   -> The code is not explicitly designed for use with either MPI or SMP.
C   -> Overindexing is not used.
C A few simple high-performance programming constructs are used:
C   -> The inner 2 loops are broken out in such a way as to enable out-of-order
C      execution of these loops as much as possible, thereby leveraging
C      vector and superscalar CPU architectures.
C   -> The outer loops are fairly long (including as many operations as
C      possible inside on a single J plane of data) in order to make effective
C      use of cache.
C Multiple time advancement algorithms are implemented for the periodic,
C channel, duct, and cavity cases in order to compare their efficiency for
C various flows on different computational architectures.  In a few of the
C algorithms, the code is broken out fully into a PASS1/PASS2 architecture
C to maximize the efficient use of cache.
C
C This code was developed as a joint project in MAE 223 (CFD), taught by
C Thomas Bewley, at UC San Diego (spring 2001, spring 2005).
C Primary contributions follow:
C Thomas Bewley was the chief software architect
C John R. Taylor wrote the channel flow solvers
C******************************************************************************|
C
C This code is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the
C Free Software Foundation; either version 2 of the License, or (at your
C option) any later version. This code is distributed in the hope that it
C will be useful, but WITHOUT ANY WARRANTY; without even the implied
C warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
C GNU General Public License for more details. You should have received a 
C copy of the GNU General Public License along with this code; if not,
C write to the Free Software Foundation, Inc., 59 Temple Place - Suite
C 330, Boston, MA 02111-1307, USA.
C
C******************************************************************************|

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      PROGRAM DIABLO
      INCLUDE 'header'
      INTEGER N
      real*8 T1,T2,GT !Double Gyre variables
      real*8 beta, Delta, c1, c2, k1, k2, e1, e2, OmB !Bickley variables
      real*8 epsilonLinear !Linear variables
      real*8 rrr,RNUM1,RNUM2
      complex compi,rrth,cphase
      real*8 ugyre, utide, uvert, epsi, o1, al, c2b, o2, zhy !Gyre2 variables
      real*8 ff


      WRITE(6,*) 
      WRITE(6,*) 
      WRITE(6,*) 
      WRITE(6,*) '             ****** WELCOME TO DIABLO ******'
      WRITE(6,*)
      WRITE(6,*) 'Note that this code is distributed under the ',
     *           'GNU General Public License.'
      WRITE(6,*) 'No warranty is expressed or implied.'
      WRITE(6,*)
      WRITE(6,*) 
      WRITE(6,*)

      PI=4.*atan(1.0)

      CALL INITIALIZE
! Initialize START_TIME for run timing
      CALL DATE_AND_TIME (VALUES=TIME_ARRAY)
      START_TIME=TIME_ARRAY(5)*3600+TIME_ARRAY(6)*60
     &   +TIME_ARRAY(7)+TIME_ARRAY(8)*0.001


! Save flow at initial time
      IF (CREATE_NEW_TH(1))THEN
      CALL SAVE_FLOW(.FALSE.)
      ENDIF

      

C A flag to determine if we are considering the first time-step
      FIRST_TIME=.TRUE.
      TS0=TIME_STEP
      DO TIME_STEP = TIME_STEP+1, TIME_STEP+N_TIME_STEPS
!        WRITE(6,*) 'Now beginning TIME_STEP = ',TIME_STEP
        T1=TIME+0.
        T2=T1
        DO RK_STEP=1,3
           T2=T2+H_BAR(RK_STEP)
          IF (NUM_PER_DIR.EQ.3) THEN

!     Flow is 2D periodic for all cases, handled by 3D periodic in Diablo with one y layer
!     FLAVOR gives flavor of the background flow type
	IF (FLAVOR .EQ. 'Bick') THEN !Bickley jet, -2.5pi ---> 2.5pi, -3 ---> 3
               beta=.6144
               Delta=sqrt(1-beta*1.5)
               c1=(1.+Delta)/3.
               c2=(1.-Delta)/3.
               k1=(6*c1)**.5
               k2=(6*c2)**.5
               e1=0.1
               e2=.3
               OmB=k1*(c1-c2)
               DO J=0,NYM
               DO K=0,NZM
               DO I=0,NXM
            U1(I,K,J)=1/cosh(GZ(K))**2*(1+2*e2*cos(k2*GX(I))*tanh(GZ(K))
     &                +2*e1*cos(k1*GX(I)-OmB*T2)*tanh(GZ(K)))-c2
            U3(I,K,J)=1/cosh(GZ(K))**2
     &                *(-k2*e2*sin(k2*GX(I))-k1*e1*sin(k1*GX(I)-OmB*T2))
               END DO
               END DO
               END DO
!     At every step, turn velocity into Fourier
               CALL FFT_XZY_TO_FOURIER(U1,CU1)
               CALL FFT_XZY_TO_FOURIER(U3,CU3)
            ELSEIF (FLAVOR .EQ. 'Linear') THEN !Linear flow u = -x, v = y with dampening (mult by bump) so that u,v=0 at boundary, -1 ---> 1, -1 ---> 1
               epsilonLinear = .1D0
               DO J=0,NYM
               DO K=0,NZM
               DO I=0,NXM
               IF (K == 0 .or. I == 0 .or. K == NZM .or. I == NXM) THEN 
                        U1(I,K,J)=0.D0
               		U3(I,K,J)=0.D0
               ELSE 
                    	U1(I,K,J)=-1.*(GX(I)-1.)*EXP(2.*epsilonLinear)
     &                  *EXP(-1.*epsilonLinear/(1.-(GX(I)-1.)**2.))
     &			*EXP(-1.*epsilonLinear/(1.-(GZ(K)-1.)**2.))
                   	U3(I,K,J)=(GZ(K)-1.)*EXP(2.*epsilonLinear)
     &			*EXP(-1.*epsilonLinear/(1.-(GX(I)-1.)**2.))
     &			*EXP(-1.*epsilonLinear/(1.-(GZ(K)-1.)**2.))
               ENDIF
               END DO
               END DO
               END DO
!     At every step, turn velocity into Fourier
               CALL FFT_XZY_TO_FOURIER(U1,CU1)
               CALL FFT_XZY_TO_FOURIER(U3,CU3)
	    ELSEIF (FLAVOR .EQ. 'Lineart') THEN !Linear flow with time dependence u = -lambda x, v = lambda y with dampening (mult by bump) so that u,v=0 at boundary, -1 ---> 1, -1 ---> 1
               epsilonLinear = .05D0
!               WRITE(*,*) 'grid is', NXM	
               DO J=0,NYM
               DO K=0,NZM
               DO I=0,NXM
               IF (K == 0 .or. I == 0 .or. K == NZM .or. I == NXM) THEN 
                        U1(I,K,J)=0.D0
               		U3(I,K,J)=0.D0
               ELSE 
                    	U1(I,K,J)=-1.*(GX(I)-1.)*EXP(2.*epsilonLinear)
     &                  *EXP(-1.*epsilonLinear/(1.-(GX(I)-1.)**2.))
     &			*EXP(-1.*epsilonLinear/(1.-(GZ(K)-1.)**2.))
                   	U3(I,K,J)=(GZ(K)-1.)*EXP(2.*epsilonLinear)
     &			*EXP(-1.*epsilonLinear/(1.-(GX(I)-1.)**2.))
     &			*EXP(-1.*epsilonLinear/(1.-(GZ(K)-1.)**2.))
               ENDIF
               END DO
               END DO
               END DO
!     At every step, turn velocity into Fourier
               CALL FFT_XZY_TO_FOURIER(U1,CU1)
               CALL FFT_XZY_TO_FOURIER(U3,CU3)
            ELSEIF (FLAVOR .EQ. 'Gyre')THEN !Double Gyre, 2x2, vary eps1, 2
               GT=0.3*sin(4*PI*T2)+0.1*sin(2*T2)
               DO J=0,NYM
               DO K=0,NZM
               DO I=0,NXM
               U1(I,K,J)=-sin(PI*(GX(I)-GT))*cos(PI*GZ(K))*PI 
               U3(I,K,J)=PI*cos(PI*(GX(I)-GT))*sin(PI*GZ(K)) 
               END DO
               END DO
               END DO
               CALL FFT_XZY_TO_FOURIER(U1,CU1)
               CALL FFT_XZY_TO_FOURIER(U3,CU3)
            ELSEIF (FLAVOR .EQ. 'Gyr2')THEN !Double Gyre in Rom Kedar paper
               ugyre=1.
               utide=1.
               uvert=0.
               epsi=.3
               o1=36.*PI/19.
               al=10.
               c2b=1.
               o2=PI
               zhy=.5
               c1=sin(o1*T2)
               c2=.5*(c2b+4./PI*sin(o2*T2))
               DO J=0, NYM
               DO K=0, NZM
               DO I=0, NXM
       U1(I,K,J)=ugyre*(.5*(2./LZ*cos(2.*PI*GZ(K)/LZ)+1./LZ
     &                  *cos(PI*GZ(K)/LZ))*sin(PI*GX(I)/(LX/2.)))
     &      +epsi*utide*c1/LZ*cos(2.*PI*GZ(K)/LZ)*sin(PI*GX(I)/(LX/2.))
       U3(I,K,J)=-ugyre*.5*(sin(2.*PI*GZ(K)/LZ)+sin(PI*GZ(K)/LZ))
     &         /(LX/2.)*cos(PI*GX(I)/(LX/2.))-epsi*utide*c1*.5/(LX/2.)
     &         *sin(2.*PI*GZ(K)/LZ)*cos(PI*GX(I)/(LX/2))+epsi*uvert*c2
     &         *.5*tanh(GZ(K)*al/LZ)*tanh((zhy-GZ(K))*al
     &                  /LZ)*tanh((GZ(K)-LZ)*al/LZ)
               ENDDO
               ENDDO
               ENDDO
               CALL FFT_XZY_TO_FOURIER(U1,CU1)
               CALL FFT_XZY_TO_FOURIER(U3,CU3)
            ELSEIF (FLAVOR .EQ. 'NFlo')THEN !No flow, pure diffusion/diffusion-reaction
               DO J=0,NYM
               DO K=0,NZM
               DO I=0,NXM
               U1(I,K,J)=0.D0
               U3(I,K,J)=0.D0
               END DO
               END DO
               END DO
!     At every step, turn velocity into Fourier
               CALL FFT_XZY_TO_FOURIER(U1,CU1)
               CALL FFT_XZY_TO_FOURIER(U3,CU3)

            ELSEIF (FLAVOR .EQ. 'Turb')THEN !2D turbulence simulation on f and beta planes
               IF((RK_STEP.eq.1).and.(Time_STEP .eq.TS0+1)) THEN !Read in Om, compute Psi at first step
!     psi-Omega formulation, initialize from fully random noise
                  ff=101.3d0
                  if (0) then !Control spin up
                   compi=(0.,1.)
                   CALL RANDOM_SEED
                   DO J=0,TNKY
                   DO K=0,TNKZ
                   DO I=0,NKX
                   CALL RANDOM_NUMBER(RNUM1) !magnitude
                   CALL RANDOM_NUMBER(RNUM2) !phase
                   rrr=(RNUM1*2.-1.)
!                   if ((K.eq.4) .and. (I.eq.4))then
!                      rrr=10.d0
!                   endif
                   cphase=RNUM2*8.*ATAN(1.0)*compi
                   rrth=cexp(cphase)
                   COm(I,K,J)=rrr*rrth
                   END DO
                   END DO
                   END DO
                  endif

                  DO J=0,TNKY
                  DO K=0,TNKZ
                  DO I=0,NKX
!                     if (I.eq.0) .and.0 K.eq.0 .and. J.eq.0
                     CPsi(I,K,J)=-COm(I,K,J)/(KX2(I)+KY2(J)+KZ2(K))
                  END DO
                  END DO
                  END DO
!                  write(*,*)CPsi(0,0,0)
!                  stop
                  CPsi(0,0,0)=0.
!     Now and at every time step Om and Psi are both Fourier
                  do I=0,NX-1
                  do K=0,NZ-1
                     xpos(I,K)=GX(I)
                     zpos(I,K)=GZ(K)
                  enddo
                  enddo
               ENDIF
            ENDIF

            IF (TIME_AD_METH.EQ.1) CALL RK_PER_1
            IF (TIME_AD_METH.EQ.2) CALL RK_PER_2            
          ELSEIF (NUM_PER_DIR.EQ.2) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_CHAN_1
            IF (TIME_AD_METH.EQ.2) CALL RK_CHAN_2            
          ELSEIF (NUM_PER_DIR.EQ.1) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_DUCT_1
            IF (TIME_AD_METH.EQ.2) CALL RK_DUCT_2            
          ELSEIF (NUM_PER_DIR.EQ.0) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_CAV_1
            IF (TIME_AD_METH.EQ.2) CALL RK_CAV_2            
          END IF
        END DO
        TIME=TIME+DELTA_T
        FIRST_TIME=.FALSE.
! Save statistics to an output file
!             write(*,*) 'flowmax', U1(32,32,0)
        IF (MOD(TIME_STEP,SAVE_STATS_INT).EQ.0) THEN
            write(*,*) 'save_stats called'
            CALL SAVE_STATS(.FALSE.)
        END IF
! Save the flow to a restart file
        IF (MOD(TIME_STEP,SAVE_FLOW_INT).EQ.0) THEN
          CALL SAVE_FLOW(.FALSE.)
        END IF
! Save the full three-dimensional fields in NETCDF format to vis.nc
! (optional)
        IF (MOD(TIME_STEP,SAVE_FLOW_INT).EQ.0) THEN
          CALL NETCDF_OPEN_VIS
          CALL NETCDF_WRITE_VIS
          CALL NETCDF_CLOSE_VIS
        END IF

! Filter the scalar field
        DO N=1,N_TH
          IF (FILTER_TH(N)
     &       .AND.(MOD(TIME_STEP,FILTER_INT(N)).EQ.0)) THEN
          write(*,*) 'Filtering...'
          CALL FILTER(N)
          END IF 
        END DO
        
      END DO

! Calculate and display the runtime for the simulation
      CALL DATE_AND_TIME (VALUES=TIME_ARRAY)
      END_TIME=TIME_ARRAY(5)*3600+TIME_ARRAY(6)*60
     &   +TIME_ARRAY(7)+TIME_ARRAY(8)*0.001
      WRITE(*,*) 'Elapsed Time (sec): ',end_time-start_time
      WRITE(*,*) 'Seconds per Iteration: '
     &     ,(end_time-start_time)/N_TIME_STEPS

      TIME_STEP=TIME_STEP-1
      CALL SAVE_FLOW(.TRUE.)
      CALL SAVE_STATS(.TRUE.)
      WRITE(6,*)
      WRITE(6,*) '        ****** Hello world!  Have a nice day! ******'
      WRITE(6,*)
      END


