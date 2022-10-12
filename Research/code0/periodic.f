C chemotaxis2 
C These solvers were written by Tom Bewley and John Taylor.
C 2D omega-psi and reaction were written by Wenbo Tang.
C******************************************************************************|

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INIT_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_PER_1
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Main time-stepping algorithm for the fully periodic case
C INPUTS  (in Fourier space):  CUi, P, and (if k>1) CFi at (k-1)  (for i=1,2,3)
C OUTPUTS (in Fourier space):  CUi, P, and (if k<3) CFi at (k)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      REAL*8 TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, TEMP6, maxu
      INTEGER I,J,K,N,ii,jj
      REAL*8 ETA1,tt,thetaa(1001)
      INTEGER RandI
      real*8 rrr,RNUM1,RNUM2
      complex compi,rrth,cphase
      real*8 x,L,T,NPZa,NPZb,NPZc,NPZd,NPZe,NPZk,NPZr,NPZs,N0,NPZal !NPZ variables
      real*8 NPZbe,NPZgam,NPZlam,NPZmu,bb,fl,ff,innerp,alph
      real*8 NPZD0,NPZD1,NPZD2
      real*8 S3(0:NXM,0:NZM,0:NYM),STH
      real*8 velox(0:NX-1,0:NZ-1),veloz(0:NX-1,0:NZ-1)

      real*8 Thi(0:N_TIME_STEPS/100,N_TH)!Synthetic field, 0:time end
      character*20 stt




      if ((TIME_STEP.eq.TS0+1).and.(RK_STEP.eq.1))then
        DO N=1,N_TH
           CALL FFT_XZY_TO_PHYSICAL(CTH(0,0,0,N),TH(0,0,0,N))
           STH=0
           Do ii=0,255
             DO jj = 0,255
                STH = STH + TH(ii,jj,0,N)
             END DO 
           END DO
           Thi(0,N)=STH/256.D0/256.D0
           CALL FFT_XZY_TO_Fourier(TH(0,0,0,N),CTH(0,0,0,N))
        END DO
      endif



      compi=(0.0,1.0)
      IF (FLAVOR.EQ.'Beta') THEN
         bb=1.612d0!nondim beta
         ff=101.03d0
         fl=1.d0!flux rate
      ENDIF

C Compute the RHS in Fourier Space, CRi.
C If using moving source, propogate the source
      IF (TRUCK) CALL RK_SOURCE

C First, define the constants used for time-stepping
      TEMP1=NU*H_BAR(RK_STEP)/2.0
      TEMP2=BETA_BAR(RK_STEP)*H_BAR(RK_STEP)
      TEMP3=ZETA_BAR(RK_STEP)*H_BAR(RK_STEP)
      TEMP4=H_BAR(RK_STEP)
      TEMP5=0.0

      IF ((FLAVOR.eq.'Batch') .AND. (delta_t.lt.0.0)) THEN
         CALL backwards_constants(TEMP1,TEMP2,TEMP3,TEMP4,TEMP5)
      ENDIF
      DO J=0,TNKY
      DO K=0,TNKZ
      DO I=0,NKX
C Start with the explicit part of the Crank-Nicolson viscous term and
C  the pressure gradient treated with Explicit Euler:
         TEMP5=1.-TEMP1*(KX2(I)+KY2(J)+KZ2(K))**1.
C The following hook should be used for regularizing the NS equations.
C If marches are backwards-in-time (ie delta_t<0) then it will be used.
         IF (DELTA_T.LT.0) CALL quasi_rev_per(TEMP5,i,j,k)
         IF ((FLAVOR.EQ.'Turb').or.(FLAVOR.EQ.'Beta'))THEN !Turbulence on f or beta
            CR1(I,K,J)=TEMP5*COm(I,K,J) !CR1 used to evolve omega
         ENDIF
C For each scalar, start with the explict part of the Crank-Nicolson
C diffusive term for each scalar
         DO N=1,N_TH
            TEMP6=1-(TEMP1/PR(N))*(KX2(I)+KY2(J)+KZ2(K))
            CRTH(I,K,J,N)=TEMP6*CTH(I,K,J,N)
         END DO

         IF (RK_STEP .GT. 1) THEN
C Add the term: ZETA_BAR(RK_STEP)*R(U(RK_STEP-1))
            IF ((FLAVOR.EQ.'Turb').or.(FLAVOR.EQ.'Beta'))THEN
               CR1(I,K,J)=CR1(I,K,J)+TEMP3*CF1(I,K,J)
            ENDIF
C Do the same for each scalar:
            DO N=1,N_TH
               CRTH(I,K,J,N)=CRTH(I,K,J,N)+TEMP3*CFTH(I,K,J,N)
            END DO
         ENDIF
      END DO
      END DO
      END DO
C If we are considering a linear background scalar gradient then add
C the term owing to advection of the background state.
C This allows us to consider a mean scalar gradient (ie stratification)
C even though the vertical boundary conditions are periodic.
C (In this case the passive scalar is a perturbation from a linear
C gradient. This gradient and the vertical domain size are used to
C make the passive scalar nondimensional, so here the nondimensional
C gradient is equal to one
      DO N=1,N_TH
      IF (BACKGROUND_GRAD(N)) THEN
C If there is a background scalar gradient add advection term:
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX 
            CFTH(I,K,J,N)=-CU2(I,K,J)
          END DO
        END DO
      END DO 
      ELSE
C Otherwise don_t
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            CFTH(I,K,J,N)=0.D0
          END DO
        END DO
      END DO 
      END IF
      END DO


C Inverse transform to physical space to compute the nonlinear terms

      IF ((FLAVOR.EQ.'Turb').or.(FLAVOR.EQ.'Beta'))THEN
!     Compute U1, U3 based on psi, U2=0.
         DO J=0,TNKY
         DO K=0,TNKZ
         DO I=0,NKX
            CU1(I,K,J)=-CIKZ(K)*CPsi(I,K,J)
            CU3(I,K,J)=CIKX(I)*CPsi(I,K,J)
         END DO
         END DO
         END DO 
      ENDIF
!     Turn into physical for nonlinear advection
      CALL FFT_XZY_TO_PHYSICAL(CU1,U1)
      CALL FFT_XZY_TO_PHYSICAL(CU3,U3)

      call interpolation(velox,U1)
      call interpolation(veloz,U3)
      DO I=0,NX-1
      DO K=0,NZ-1
         xpos(I,K)=xpos(I,K)+velox(I,K)*TEMP4
         zpos(I,K)=zpos(I,K)+veloz(I,K)*TEMP4
      ENDDO
      ENDDO

      if ((RK_STEP.eq.3).and.(Mod(TIME_STEP,100).eq.0))then
!      write(*,*)xpos(5,5),zpos(5,5)
         stt='pos.'
     &        //CHAR(MOD(TIME_STEP,10000000)/1000000+48)
     &        //CHAR(MOD(TIME_STEP,1000000)/100000+48)
     &        //CHAR(MOD(TIME_STEP,100000)/10000+48)
     &        //CHAR(MOD(TIME_STEP,10000)/1000+48)
     &        //CHAR(MOD(TIME_STEP,1000)/100+48)
     &        //CHAR(MOD(TIME_STEP,100)/10+48)
     &        //CHAR(MOD(TIME_STEP,10)+48)
         open(unit=11,file=stt,form="unformatted")
         write(11) xpos,zpos
         close(11)
      endif
      

C Here we compute various terms in reaction, nonlinear, physical space
C Note that it is assumed that the nutrient concentration is in scalar #1
      DO N=1,N_TH 
        CALL FFT_XZY_TO_PHYSICAL(CTH(0,0,0,N),TH(0,0,0,N))
      END DO

!      write(*,*)TH(316,198,0,1),TH(316,198,0,2),TH(316,198,0,3)


      IF (FLAVORS .eq. 'FKPP')THEN!FKPP only 1 species
         ETA1=10.0d0
         DO J=0,NYM
         DO K=0,NZM
         DO I=0,NXM
            S1(I,K,J)=ETA1*TH(I,K,J,1)*(1.d0-TH(I,K,J,1))
         END DO
         END DO
         END DO
         CALL FFT_XZY_TO_FOURIER(S1,CS1)
         DO J=0,TNKY
         DO K=0,TNKZ
         DO I=0,NKX 
            CFTH(I,K,J,1)=CFTH(I,K,J,1) + CS1(I,K,J)
         END DO
         END DO
         END DO
      ENDIF
      
      IF (FLAVORS .eq. 'BSTA')THEN
         NPZal = 0.0 
         NPZk = 10.d0
         DO J=0,NYM
         DO K=0,NZM
         DO I=0,NXM
       S1(I,K,J)=NPZk*(TH(I,K,J,1)+1.d-2)*(NPZal-TH(I,K,J,1))
     &           *(TH(I,K,J,1)-1.D0)
         END DO
         END DO
         END DO
         CALL FFT_XZY_TO_FOURIER(S1,CS1)
         DO J=0,TNKY
         DO K=0,TNKZ
         DO I=0,NKX
            CFTH(I,K,J,1)=CFTH(I,K,J,1) + CS1(I,K,J)
         END DO
         END DO
         END DO
      ENDIF

      IF (FLAVORS .eq. 'NPZ')THEN!NPZ Limit cycle based on Edwards & Brindley
         x=86400.d0
         L=10000.d0
         T=100000.d0
         NPZa=0.2d0*L*T/x
         NPZb=0.2d0*L
         NPZc=0.4d0*L
         NPZd=1.5d0*T/x
         NPZe=.03d0
         NPZk=.05d0*T/x
         NPZr=.15d0*T/x
         NPZs=.04d0*T/x
         N0=.6d0
         NPZal=.25d0
         NPZbe=.33d0
         NPZgam=.5d0
         NPZlam=.6d0*T/x
         NPZmu=0.035d0
!     Nutrient
         DO J=0,NYM
         DO K=0,NZM
         DO I=0,NXM
             S1(I,K,J)=-NPZa*TH(I,K,J,1)*TH(I,K,J,2)/(NPZe+TH(I,K,J,1))
     & /(NPZb+NPZc*TH(I,K,J,2))    +NPZr*TH(I,K,J,2)   +NPZbe*NPZlam
     & *TH(I,K,J,3)*TH(I,K,J,2)**2./(NPZmu**2.+TH(I,K,J,2)**2.)
     & +NPZgam*NPZd*TH(I,K,J,3)**2.+NPZk*(N0-TH(I,K,J,1))
         END DO
         END DO
         END DO      
         IF (FLAVOR.eq.'Beta')THEN!beta plane vertical flux
            S1(I,K,J)=S1(I,K,J)-fl*bb/ff*U3(I,K,J)*(N0-TH(I,K,J,1))
         ENDIF
         CALL FFT_XZY_TO_FOURIER(S1,CS1)
         DO J=0,TNKY
         DO K=0,TNKZ
         DO I=0,NKX 
            CFTH(I,K,J,1)=CFTH(I,K,J,1) + CS1(I,K,J)
         END DO
         END DO
         END DO
!     Phytoplankton
         DO J=0,NYM
         DO K=0,NZM
         DO I=0,NXM
            S1(I,K,J)=NPZa*TH(I,K,J,1)*TH(I,K,J,2)/(NPZe+TH(I,K,J,1))
     & /(NPZb+NPZc*TH(I,K,J,2))    -NPZr*TH(I,K,J,2)   -NPZlam
     & *TH(I,K,J,3)*TH(I,K,J,2)**2./(NPZmu**2.+TH(I,K,J,2)**2.)
     & -(NPZs+NPZk)*TH(I,K,J,2)
         END DO
         END DO
         END DO      
         CALL FFT_XZY_TO_FOURIER(S1,CS1)
         DO J=0,TNKY
         DO K=0,TNKZ
         DO I=0,NKX 
            CFTH(I,K,J,2)=CFTH(I,K,J,2) + CS1(I,K,J)
         END DO
         END DO
         END DO
!     Zooplankton
         DO J=0,NYM
         DO K=0,NZM
         DO I=0,NXM
            S1(I,K,J)= NPZal*NPZlam
     & *TH(I,K,J,3)*TH(I,K,J,2)**2./(NPZmu**2.+TH(I,K,J,2)**2.)
     & -NPZd*TH(I,K,J,3)**2.
         END DO
         END DO
         END DO      
         CALL FFT_XZY_TO_FOURIER(S1,CS1)
         DO J=0,TNKY
         DO K=0,TNKZ
         DO I=0,NKX 
            CFTH(I,K,J,3)=CFTH(I,K,J,3) + CS1(I,K,J)
         END DO
         END DO
         END DO
      ENDIF

      IF (FLAVORS.eq.'NPZ2') THEN!Bistable NPZ (LC+FP) by Zhang & Wang
         NPZa = 1.d0
         NPZb = 1.d0
         NPZc = 5.d0
         NPZd = 0.1d0
         NPZD0= 0.1d0
         NPZD1= 0.2d0
         NPZD2= 2.1d0
         N0   = 9.96d0
         NPZal= 1.d0
         NPZbe= 0.5d0
         DO J=0,NYM
         DO K=0,NZM
         DO I=0,NXM
            IF ( TH(I,K,J,1) .LE.  NPZa) THEN
               S3(I,K,J) = NPZb/NPZa*TH(I,K,J,1)
            ELSE
               S3(I,K,J) = NPZb
            END IF
         ENDDO
         ENDDO
         ENDDO
C     Nutrient
         DO J=0,NYM
         DO K=0,NZM
         DO I=0,NXM
            S1(I,K,J)= NPZD0*(N0 - TH(I,K,J,1)) - S3(I,K,J)*TH(I,K,J,2)
         END DO
         END DO
	 END DO
         IF (FLAVOR.eq.'Beta')THEN!beta plane vertical flux
            S1(I,K,J)=S1(I,K,J)-fl*bb/ff*U3(I,K,J)*(N0-TH(I,K,J,1))
         ENDIF
         CALL FFT_XZY_TO_FOURIER(S1,CS1)
         DO J=0,TNKY
         DO K=0,TNKZ
         DO I=0,NKX
            CFTH(I,K,J,1)=CFTH(I,K,J,1) + CS1(I,K,J)
         END DO
         END DO
         END DO
C     Phytoplankton
         DO J=0,NYM
         DO K=0,NZM
         DO I=0,NXM
            S1(I,K,J)= NPZal*S3(I,K,J)*TH(I,K,J,2) -
     &           ((NPZc*TH(I,K,J,2))/(1.D0 + NPZd*TH(I,K,J,2)))*
     &           TH(I,K,J,3) - NPZD1*Th(I,K,J,2)
         END DO
         END DO
         END DO
         CALL FFT_XZY_TO_FOURIER(S1,CS1)
         DO J=0,TNKY
         DO K=0,TNKZ
         DO I=0,NKX
            CFTH(I,K,J,2)=CFTH(I,K,J,2) + CS1(I,K,J)
         END DO
         END DO
         END DO
C     Zooplankton
         DO J=0,NYM
         DO K=0,NZM
         DO I=0,NXM
          S1(I,K,J)= NPZbe*((NPZc*TH(I,K,J,2))/(1.D0+NPZd*TH(I,K,J,2)))*
     &           TH(I,K,J,3) - NPZD2*TH(I,K,J,3)
         END DO
         END DO
         END DO
         CALL FFT_XZY_TO_FOURIER(S1,CS1)
         DO J=0,TNKY
         DO K=0,TNKZ
         DO I=0,NKX
            CFTH(I,K,J,3)=CFTH(I,K,J,3) + CS1(I,K,J)
         END DO
         END DO
         END DO                     
      END IF

      IF (FLAVORS .eq. 'chmo') THEN!chemotaxis reaction

C Change Nutrient back to spectral
      CALL FFT_XZY_TO_FOURIER(TH(0,0,0,1),CTH(0,0,0,1))

C Here, calculate each of the chemotaxis terms
      DO N = 2,N_TH
C First, CHI*d/dx(B*dC/dx):
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            CS1(I,K,J)=CIKX(I)*CTH(I,K,J,1)
          END DO
        END DO
      END DO
      CALL FFT_XZY_TO_PHYSICAL(CS1,S1)
      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=S1(I,K,J)*TH(I,K,J,N) !N=2
          END DO
        END DO
      END DO
      CALL FFT_XZY_TO_FOURIER(S1,CS1)
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            CFTH(I,K,J,N)=CFTH(I,K,J,N)-CIKX(I)*CS1(I,K,J)*V_S(2)
          END DO
        END DO
      END DO

C Finally, calculate the CHI*d/dz(B*dC/dz) term:
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            CS1(I,K,J)=CIKZ(K)*CTH(I,K,J,1)
          END DO 
        END DO
      END DO
      CALL FFT_XZY_TO_PHYSICAL(CS1,S1)
      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=S1(I,K,J)*TH(I,K,J,N)!N=2
          END DO
        END DO
      END DO
      CALL FFT_XZY_TO_FOURIER(S1,CS1)
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            CFTH(I,K,J,N)=CFTH(I,K,J,N)-CIKZ(K)*CS1(I,K,J)*V_S(2)
          END DO
        END DO
      END DO
      END DO
C End loop over bacterial species

C End if CHEMOTAXIS


C Now that we are done with the Chemotaxis terms, we need to transform the
C nutrient concentration to physical space for calculating the nonlinear
C advection terms

      CALL FFT_XZY_TO_PHYSICAL(CTH(0,0,0,1),TH(0,0,0,1))
!      CALL FFT_XZY_TO_PHYSICAL(CS1,S1)
C Take care of uptake in the C equation


          DO J=0,NYM
            DO K=0,NZM
              DO I=0,NXM 
                S1(I,K,J)=0.d0
                IF ((TH(I,K,J,2).gt.0.d0)
     &             .and.(TH(I,K,J,1).gt.0.d0)) THEN
                  S1(I,K,J)=C0(2)*TH(I,K,J,1)*TH(I,K,J,2)
                END IF
              END DO
            END DO
          END DO
! Add an extra uptake term for the uniform non-motile bacteria  
          DO j=0,NYM
            DO k=0,NZM
              DO i=0,NXM
                IF (TH(I,K,J,1).gt.0.d0)THEN
                S1(I,K,J)=S1(I,K,J)+C0(2)*TH(I,K,J,1)*1.d0
                ENDIF
              END DO
            END DO
          END DO
          CALL FFT_XZY_TO_FOURIER(S1,CS1)
          DO J=0,TNKY
            DO K=0,TNKZ
              DO I=0,NKX
                CFTH(I,K,J,1)=CFTH(I,K,J,1)-CS1(I,K,J)
              END DO
            END DO
          END DO
      END IF


! Pure diffusion, do nothing
      IF (FLAVORS.eq.'diff')then
      ENDIF




C Compute the nonlinear terms for the passive scalar equation
C Do this before the nonlinear momentum terms to use Fi as a working
C array before using it for the momentum equation.
      IF (TRUCK) CALL INIT_FORCING
      DO N=1,N_TH
        DO J=0,NYM
          DO K=0,NZM
            DO I=0,NXM
              F1(I,K,J)=U1(I,K,J)*TH(I,K,J,N)
              F3(I,K,J)=U3(I,K,J)*TH(I,K,J,N)
            END DO
          END DO
        END DO
        CALL FFT_XZY_TO_FOURIER(F1,CF1)
        CALL FFT_XZY_TO_FOURIER(F3,CF3)
        DO J=0,TNKY
          DO K=0,TNKZ
            DO I=0,NKX
              CFTH(I,K,J,N)=CFTH(I,K,J,N)-CIKX(I)*CF1(I,K,J)
     &                      -CIKZ(K)*CF3(I,K,J)
C Add the forcing due to the moving source
              IF (TRUCK) CALL FORCING(I,K,J)
C Add R-K terms for the TH equation to the RHS
              CRTH(I,K,J,N)=CRTH(I,K,J,N)+TEMP2*CFTH(I,K,J,N)
            END DO
          END DO
        END DO
      END DO
C The RHS vector for the TH equation is now ready
     
C Compute the nonlinear terms for the momentum equations

      IF ((FLAVOR.eq.'Turb').or.(FLAVOR.eq.'Beta'))THEN
         CALL FFT_XZY_TO_PHYSICAL(COm,Om)
         DO J=0,NYM
         DO K=0,NZM
         DO I=0,NXM
            F1(I,K,J)=U1(I,K,J)*Om(I,K,J)
            F3(I,K,J)=U3(I,K,J)*Om(I,K,J)
         END DO
         END DO
         END DO
         CALL FFT_XZY_TO_FOURIER(Om,COm)
         CALL FFT_XZY_TO_FOURIER(F1,CF1)
         CALL FFT_XZY_TO_FOURIER(F3,CF3)
C Here we start constructing the R-K terms in CFi
C Note, that the order of the following operations are important
         DO J=0,TNKY
         DO K=0,TNKZ
         DO I=0,NKX
            IF (GUSTS) CALL GUST(I,K,J)
            CF1(I,K,J)=-CIKX(I)*CF1(I,K,J)
     &           -CIKZ(K)*CF3(I,K,J)
         END DO
         END DO
         END DO
         IF (FLAVOR.EQ.'Beta')THEN !-beta*v
            CALL FFT_XZY_TO_FOURIER(U3,CU3)
            DO J=0,TNKY
            DO K=0,TNKZ
            DO I=0,NKX
            CF1(I,K,J)=CF1(I,K,J)-bb*CU3(I,K,J)
            END DO
            END DO
            END DO
            CALL FFT_XZY_TO_PHYSICAL(CU3,U3)         
         ENDIF
      ENDIF
C Done with the computation of nonlinear terms


C If the scalar is active (RI_TAU NE 0), add the bouyancy forcing term
C as explicit R-K
       DO N=1,N_TH
C First, convert back to Fourier space
       CALL FFT_XZY_TO_FOURIER(TH(0,0,0,N),CTH(0,0,0,N))
       DO J=0,TNKY
         DO K=0,TNKZ
           DO I=0,NKX
             CF2(I,K,J)=CF2(I,K,J)+RI_TAU(N)*CTH(I,K,J,N)
           END DO
         END DO
       END DO
       END DO

C Now, add the R-K terms to the RHS
C Note, this has already been done for the scalar field TH
       DO J=0,TNKY
         DO K=0,TNKZ
           DO I=0,NKX
             CR1(I,K,J)=CR1(I,K,J)+TEMP2*CF1(I,K,J)!Just Om
           END DO
         END DO
       END DO

C Computation of CRi complete.
C If Variable timestepping and done with one full R-K step, update
C DELTA_T based on the specified CFL number
! Note, this change will not take effect until the next timestep
! since the TEMP variables have already been defined
      IF ((VARIABLE_DT).and.(RK_STEP.eq.3)
     &        .and.(MOD(TIME_STEP,UPDATE_DT).EQ.0)) THEN
        CALL COURANT
      END IF
C Now solve the implicit system for the intermediate field.
C (In the fully-periodic case, this is easy!)
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
             kkk=SQRT(KX2(I)+KY2(J)+KZ2(K))
!           Loss is largest at the largest scale
            TEMP5=1+TEMP1*(KX2(I)+KY2(J)+KZ2(K))**1.
     &              +1./(kkk**2)*TEMP4 !Last term energy loss at large scale)
C The following hook should be used for regularizing the NS equations.
C If marches are backwards-in-time (ie delta_t<0) then it will be used.		
            IF (DELTA_T.LT.0) CALL quasi_rev_per(TEMP5,i,j,k)
            IF ((FLAVOR.EQ.'Turb').or.(FLAVOR.eq.'Beta'))THEN
            COm(I,K,J)=CR1(I,K,J)/TEMP5
            CPsi(I,K,J)=-COm(I,K,J)/(KX2(I)+KY2(J)+KZ2(K))
            ENDIF
            DO N=1,N_TH
              TEMP6=1+(TEMP1/PR(N))*(KX2(I)+KY2(J)+KZ2(K)) 
              CTH(I,K,J,N)=CRTH(I,K,J,N)/TEMP6
            END DO
          END DO
        END DO
      END DO
      IF ((FLAVOR.EQ.'Turb').or.(FLAVOR.eq.'Beta')) CPsi(0,0,0)=0.      


      if ((mod(TIME_STEP,10000000).eq.0).and.(RK_STEP.eq.3))then
      CALL FFT_XZY_TO_PHYSICAL(COm,Om)
      do N=1,N_TH
      CALL FFT_XZY_TO_PHYSICAL(CTH(0,0,0,N),TH(0,0,0,N))
      enddo
      open(1001,name='test',form='formatted')
      do ii=0,NXM
         do jj=0,NZM
            write(1001,'(6D16.8)') U1(ii,jj,0),U3(ii,jj,0),
     &   TH(ii,jj,0,1),TH(ii,jj,0,2),TH(ii,jj,0,3),
     &   Om(ii,jj,0)
         enddo
      enddo
      close(1001)
      CALL FFT_XZY_TO_FOURIER(Om,COm)
      do N=1,N_TH
      CALL FFT_XZY_TO_FOURIER(TH(0,0,0,N),CTH(0,0,0,N))
      enddo
      
      endif


      IF ((FLAVOR.eq.'Turb').AND.(RK_STEP.eq.1)!Initial energy at intermediate length
     &     .AND. (TIME_STEP.eq.TS0+1)) THEN
        EK0=0.d0
        DO J=0,TNKY
          DO K=0,TNKZ
            DO I=0,NKX
             IF ((SQRT(KX2(I)+KY2(J)+KZ2(K)).le.10d0).and.
     &           (SQRT(KX2(I)+KY2(J)+KZ2(K)).gt.6d0)) THEN
                EK0=EK0+COm(I,K,J)*CONJG(COm(I,K,J))
             END IF
            END DO
          END DO
        END DO
      ENDIF
      EK0=2.5d0

      IF ((FLAVOR.eq.'Turb').AND.(RK_STEP.eq.3)) THEN


! Add forcing to maintain a constant energy at intermediate wavenumbers
! First, calculate the energy at intermediate wavenumbers

        EK=0.d0
        DO J=0,TNKY
          DO K=2,2!0,TNKZ
            DO I=2,2!0,NKX
!               IF ((SQRT(KX2(I)+KY2(J)+KZ2(K)).le.5.d0).and.
!     &             (SQRT(KX2(I)+KY2(J)+KZ2(K)).gt.3.d0)) THEN
                EK=EK+COm(I,K,J)*CONJG(COm(I,K,J))
!               END IF
            END DO
          END DO
        END DO
! Scale it to the original value at start of simulation
        DO J=0,TNKY
          DO K=2,2!0,TNKZ
            DO I=2,2!0,NKX
!              IF ((SQRT(KX2(I)+KY2(J)+KZ2(K)).le.5.d0).and.
!     &            (SQRT(KX2(I)+KY2(J)+KZ2(K)).gt.3.d0)) THEN
!                   CALL RANDOM_NUMBER(RNUM1) !magnitude
                   CALL RANDOM_NUMBER(RNUM2) !phase
                   cphase=(RNUM2-.5)*8.*ATAN(1.0)/50.*compi
                   rrth=cexp(cphase)
          innerp=real(COm(I,K,J))*real(rrth)+imag(COm(I,K,J))*imag(rrth)
          alph=sqrt(innerp**2.-COm(I,K,J)*CONJG(COm(I,K,J))+EK0)-innerp

!         write(*,*)COm(I,K,J),rrth*alph,COm(I,K,J)+rrth*alph
!          stop

        COm(I,K,J)=COm(I,K,J)*sqrt(EK0/EK)
!        COm(I,K,J)=COm(I,K,J)+rrth*alph
!              END IF
            END DO
          END DO
        END DO
        EK1=0.d0
        DO J=0,TNKY
          DO K=2,2!0,TNKZ
            DO I=2,2!0,NKX
!               IF ((SQRT(KX2(I)+KY2(J)+KZ2(K)).le.20.d0).and.
!     &             (SQRT(KX2(I)+KY2(J)+KZ2(K)).gt.15.d0)) THEN
                EK1=EK1+COm(I,K,J)*CONJG(COm(I,K,J))
!              END IF
            END DO
          END DO
        END DO
!        write(*,*) EK0, EK, EK1
!        stop

      END IF



      
      CALL FFT_XZY_TO_FOURIER(U1,CU1)
      CALL FFT_XZY_TO_FOURIER(U3,CU3)

      if ((mod(TIME_STEP,100).eq.0).and.(RK_STEP.eq.3))then
        DO N=1,N_TH
           CALL FFT_XZY_TO_PHYSICAL(CTH(0,0,0,N),TH(0,0,0,N))
           STH=0
           Do ii=0,255
              DOjj = 0,255
              STH = STH + TH(ii,jj,0,N)
                   END DO
           END DO
           Thi((Time_Step-TS0)/100,N)=STH/256.D0/256.D0
           CALL FFT_XZY_TO_Fourier(TH(0,0,0,N),CTH(0,0,0,N))
        END DO
        if (TIME_STEP.eq.TS0+N_TIME_STEPS)then
!            write(*,*)Thi
            open(22,file='Thi',form='unformatted')
            write(22)Thi
            close(22)
         endif

      endif




      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_PER_2
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Alternative time-stepping algorithm for the fully periodic case.
C INPUTS  (in Fourier space):  CUi, P, and (if k>1) CFi at (k-1)  (for i=1,2,3)
C OUTPUTS (in Fourier space):  CUi, P, and (if k<3) CFi at (k)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE REM_DIV_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER I,J,K
      REAL*8  TEMP5

C Compute phi, store in the variable CR1.
C Note the coefficient H_BAR is absorbed into phi.
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            TEMP5=-(KX2(I)+KY2(J)+KZ2(K)+EPS)
            CR1(I,K,J)=(CIKX(I)*CU1(I,K,J)+CIKY(J)*CU2(I,K,J)+
     *                  CIKZ(K)*CU3(I,K,J))/TEMP5
          END DO
        END DO
C Then update the CUi to make velocity field divergence-free.
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,J)=CU1(I,K,J)-CIKX(I)*CR1(I,K,J)
            CU2(I,K,J)=CU2(I,K,J)-CIKY(J)*CR1(I,K,J)
            CU3(I,K,J)=CU3(I,K,J)-CIKZ(K)*CR1(I,K,J)
          END DO
        END DO
      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE POISSON_P_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER I,J,K

      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            CR1(I,K,J)=CIKX(I)*CU1(I,K,J)
            CR2(I,K,J)=CIKZ(K)*CU3(I,K,J)
          END DO
        END DO
      END DO

      CALL FFT_XZY_TO_PHYSICAL(CR1,R1)
      CALL FFT_XZY_TO_PHYSICAL(CR2,R2)

      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
            P(I,K,J)=R1(I,K,J)*R1(I,K,J)+R1(I,K,J)*R2(I,K,J)+
     *               R2(I,K,J)*R2(I,K,J)
          END DO
        END DO
      END DO
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            CF1(I,K,J)=CIKY(J)*CU1(I,K,J)
            CF2(I,K,J)=CIKX(I)*CU2(I,K,J)
            CF3(I,K,J)=CIKZ(K)*CU1(I,K,J)
            CR1(I,K,J)=CIKX(I)*CU3(I,K,J)
            CR2(I,K,J)=CIKZ(K)*CU2(I,K,J)
            CR3(I,K,J)=CIKY(J)*CU3(I,K,J)
          END DO
        END DO
      END DO
      CALL FFT_XZY_TO_PHYSICAL(CF1,F1)
      CALL FFT_XZY_TO_PHYSICAL(CF2,F2)
      CALL FFT_XZY_TO_PHYSICAL(CF3,F3)
      CALL FFT_XZY_TO_PHYSICAL(CR1,R1)
      CALL FFT_XZY_TO_PHYSICAL(CR2,R2)
      CALL FFT_XZY_TO_PHYSICAL(CR3,R3)
      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
            P(I,K,J)=2*(P(I,K,J)+F1(I,K,J)*F2(I,K,J)
     *                          +F3(I,K,J)*R1(I,K,J)
     *                          +R2(I,K,J)*R3(I,K,J))
          END DO
        END DO
      END DO

      CALL FFT_XZY_TO_FOURIER(P,CP)
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            CP(I,K,J)=CP(I,K,J)/(KX2(I)+KY2(J)+KZ2(K)+EPS)
          END DO
        END DO
      END DO
      CP(0,0,0)=0.0

      RETURN
      END
	
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine forcing(i,k,j)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	include 'header'
	integer i,j,k,n
	
	if (num_per_dir.eq.3) then
	  
C Translate periodic forcing to source location
      do n=1,N_TH
	  CFTH(i,k,j,n) = CFTH(i,k,j,n) - CS2(i,k,j) + CF(i,k,j) * 
     &		   exp(-ci*(KX(i)*pos(1)+KY(j)*pos(2)+KZ(k)*pos(3)))
      end do
	else
	  write(6,*) 'Non-periodic cases not supported'
	end if
	
	
	return
	end
	
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine init_forcing
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	include 'header'
	integer n,i,k,j
	
	if (num_per_dir.eq.3) then
        do n=1,n_th
	    do j=0,nym
	      do k=0,nzm
		  do i=0,nxm
		    S2(i,k,j) = D(i,k,j)*TH(i,k,j,n)
	        end do
	      end do
	    end do
	  end do
	  call fft_xzy_to_fourier(S2,CS2)
	else
	  write(6,*) 'Non-periodic cases not supported'
	end if
	
	
	return
	end
	
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine gust(i,k,j)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	include 'header'
	integer i,j,k
	real*8 mag,zbqlnor,mu,var
	
	if (num_per_dir.eq.3) then
	    mu = 0.
          var = 70.
	    mag = zbqlnor(mu,sqrt(var))
	    CF1(i,k,j) = CF1(i,k,j) + mag*KZ(k)*exp(-abs(KZ(k)))  
	else
	  write(6,*) 'Non-periodic cases not supported'
	end if
	
	
	return
	end
	
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine rk_source
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	include 'header'
	real*8 temp2, temp3, b,a,gammam,gammae
	real*8 zbqlnor,mu,mag, tpos(3)
	integer i
	
	gammam = 1
	gammae = 10.
	
	temp2 = beta_bar(rk_step)*h_bar(rk_step)
	temp3 = zeta_bar(rk_step)*h_bar(rk_step)
	mu = 0.
	
	  if (rk_step .gt. 1) then
	    do i=1,3
	      pos(i) = pos(i) + temp3*fpos(i)
	    end do
	    theta = theta + temp3*ftheta
	  end if
	  
	  tpos(1) = pos(1)-LX/2.
	  tpos(2) = pos(2)
	  tpos(3) = pos(3)-LZ/2.
	  if (tpos(1).ne.0) then
	    b = atan(tpos(3)/tpos(1))
	    if (tpos(1).lt.0) b = b+pi
	  else
	    b = tpos(3)/abs(tpos(3))*pi
	  end if
	  
	  
	  a = gammam*((max(abs(tpos(1)),abs(tpos(3)))
     &		   /(LX*width*0.95/2.))**gammae)

	  mag = a * (pi-mod(theta-b,2.*pi))
	  ftheta = zbqlnor(mu,sqrt(theta_var)) + mag
	  mag = zbqlnor(v0,sqrt(mag_var))
	  
	  fpos(1) = mag*cos(theta)
	  fpos(2) = 0.
	  fpos(3) = mag*sin(theta)
	  
	  do i=1,3
	    pos(i) = pos(i) + temp2*fpos(i)
	  end do
        theta = theta + temp2*ftheta	  

	
	
	return
	end
	
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INIT_TRUCK
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      INCLUDE 'header'
	real version, current_version
	real*8 a,b,px,py,pz,peak,slope,rad
	integer i,j,k,norm


	open(12,file='input_truck.dat',form='formatted',status='old')  
	
	current_version = 1.0
	read(12,*)
	read(12,*)
	read(12,*)
	read(12,*)
	read(12,*) version
	if (version.ne.current_version) stop 'Wrong input data format.'
	read(12,*)
	read(12,*) width
	read(12,*)
	read(12,*) a, b
	read(12,*)
	read(12,*) mag_var, theta_var
	read(12,*)
	read(12,*) px, py, pz, v0, theta
	
	close(12)
      
	
C Initialize forcing
	pos(1) = px+.0000001
	pos(2) = py
	pos(3) = pz+.0000001
	do j=0,NY
	  do k=0,NZ
	    do i=0,NX
	      F(i,k,j) = b*exp(-a*((gx(i)-LX/2.)**2+(gy(j)-LY/2.)**2
     &	                    +(gz(k)-LZ/2.)**2))
	    end do
	  end do
	end do
	
	call fft_xzy_to_fourier(F,CF)
	
	do j=0,tnky
	  do k=0,tnkz
	    do i=0,nkx
	      CF(i,k,j) = CF(i,k,j) * 
     &		   exp(ci*(KX(i)*LX/2+KY(j)*LY/2.+KZ(k)*LZ/2))
	    end do
	  end do
	end do
C Initialize perimeter damping
      peak = 10
	slope = 30
	norm = 10
      do j=0,NY
	  do k=0,NZ
	    do i=0,NX
	     rad = ((GX(i)-LX/2.)**norm + (GZ(k)-LZ/2.)**norm)**(1./norm)
	      D(i,k,j) = peak*(1./pi*atan(slope*(rad-LX*0.45))+0.5)
	    end do
	  end do
	end do


      RETURN
      END



C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_FLOW_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER I, J, K
      REAL*8 RNUM1,RNUM2,RNUM3

C For an initial vortex, define the location of the centerline
      REAL*8 XC(0:NY+2),ZC(0:NY+2)

      IF ((FLAVOR .EQ. 'Basic').OR.(FLAVOR.eq.'Turb')) THEN
	
      WRITE(6,*) 'Creating new flow from scratch.'

C Initialize random number generator
      CALL RANDOM_SEED


      IF (IC_TYPE.eq.0) THEN
C Initizlize the flow using a Taylor-Green vortex
C Nondimensionalize with U0 and 1/kappa
      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
! Add a random phase
!            U1(I,K,J)=cos(2*pi*(GY(J))/LY)
!     &               *cos(2*pi*(GX(I))/LX)
!     &               *SIN(2*pi*(GZ(K))/LZ)
!            U2(I,K,J)=0.d0
!            U3(I,K,J)=-cos(2*pi*(GY(J))/LY)
!     &               *sin(2*pi*(GX(I))/LX)
!     &               *COS(2*pi*(GZ(K))/LZ)
            U1(I,K,J)=sin(2*pi*(GZ(K))/LZ)
     &               *cos(2*pi*(GX(I))/LX)
!     &               *SIN(2*pi*(GZ(K))/LZ)
            U2(I,K,J)=0.d0
            U3(I,K,J)=-cos(2*pi*(GZ(K))/LZ)
     &               *sin(2*pi*(GX(I))/LX)
!     &               *COS(2*pi*(GZ(K))/LZ)
          END DO
        END DO
      END DO
      ELSE IF (IC_TYPE.eq.1) THEN
C Start with an ideal vortex centered in the domain
      DO J=0,NYM
        XC(J)=LX/2.+(LX/10.)*sin(2*PI*GY(J)/LY)
!        XC(J)=LX/2.
        ZC(J)=LZ/2.
        DO K=0,NZM
          DO I=0,NXM
            IF ((GX(I)-XC(j))**2.+(GZ(K)-ZC(j))**2..gt.0.1) then
! If we aren't too close to the vortex center
              U1(I,K,J)=-1.d0*(GZ(K)-ZC(j))
     &                /((GX(I)-XC(j))**2.+(GZ(K)-ZC(j))**2.)
              U3(I,K,J)=1.d0*(GX(I)-XC(j))
     &                /((GX(I)-XC(j))**2.+(GZ(K)-ZC(j))**2.)
              U2(I,K,J)=0.d0
            ELSE
! Otherwise:
              U1(I,K,J)=-1.d0*(GZ(K)-ZC(j))
     &                /0.1
              U3(I,K,J)=1.d0*(GX(I)-XC(j))
     &                /0.1
              U2(I,K,J)=0.d0
            END IF 
            CALL RANDOM_NUMBER(RNUM1)
            CALL RANDOM_NUMBER(RNUM2)
            CALL RANDOM_NUMBER(RNUM3)
            U1(I,K,J)=U1(I,K,J)+(RNUM1-0.5)*KICK
            U2(I,K,J)=U2(I,K,J)+(RNUM1-0.5)*KICK
            U3(I,K,J)=U3(I,K,J)+(RNUM1-0.5)*KICK
          END DO
        END DO
      END DO
      ELSE
        WRITE(*,*) 'Warning, Undefined Initial conditions in periodic.f'
      END IF


      CALL FFT_XZY_TO_FOURIER(U1,CU1)
      CALL FFT_XZY_TO_FOURIER(U2,CU2)
      CALL FFT_XZY_TO_FOURIER(U3,CU3)
      DO J=0,2
        DO K=0,2
          DO I=0,2
            CALL RANDOM_NUMBER(RNUM1)
            CALL RANDOM_NUMBER(RNUM2)
            CALL RANDOM_NUMBER(RNUM3)
            CU1(I,K,J)=CU1(I,K,J)+(RNUM1-0.5)*KICK
            CU2(I,K,J)=CU2(I,K,J)+(RNUM1-0.5)*KICK
            CU3(I,K,J)=CU3(I,K,J)+(RNUM1-0.5)*KICK
          end do
        end do
      end do
	
	ELSE IF(FLAVOR .EQ. 'Ensemble') THEN
	
	CALL CREATE_FLOW_ENSEM
	
	ELSE
	write(6,*) 'Unknown flavor, flow-field not created'
	
	endif

! get the initial energy in low wavenumbers
      CALL FFT_XZY_TO_PHYSICAL(CU1,U1)
      CALL FFT_XZY_TO_PHYSICAL(CU2,U2)
      CALL FFT_XZY_TO_PHYSICAL(CU3,U3)
      EK0=0.d0
      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
              EK0=EK0+U1(I,K,J)**2.d0+U2(I,K,J)**2.d0+U3(I,K,J)**2.d0
          END DO
        END DO
      END DO
      write(*,*) 'EK0: ',EK0
      EPSILON_TARGET=((3.d0/DX(1))**4.d0)*(NU**3.d0)*(PR(1))**(-2.d0)
      write(*,*) 'EPSILON_TARGET: ',EPSILON_TARGET
      CALL FFT_XZY_TO_FOURIER(U1,CU1)
      CALL FFT_XZY_TO_FOURIER(U2,CU2)
      CALL FFT_XZY_TO_FOURIER(U3,CU3)


      CALL REM_DIV_PER
      CALL POISSON_P_PER

      CALL SAVE_STATS_PER(.FALSE.)

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_TH_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER I,J,K,N

      INTEGER si,sj,sk,indG
      real*8 GXs(0:NXM), GYs(0:NYM), GZs(0:NZM),rr
      real*8 sigmaPR

C Note, Since stratification is not permitted in the periodic flow field
C Any background stratification must be added to the governing equations


! HYP: 32, 128
! EDD: 72, 48

! Eddy 65, 135
! Hyp  55,195

!For turb: eddy: 160, 180, hyp: 100, 88

      DO N=1,N_TH
      IF (CREATE_NEW_TH(N)) THEN
      si=100
         sj=1   
      sk=88
         WRITE(*,*) 'Blob center coordinate',si,sk
         DO I=0,NXM
            indG=mod(I+(NX/2-si),NX)
            if (indG.lt.0)then
               indG=indG+NX
            endif
      !      GXs(I)=GX(indG)
           GXs(I)=GX(I)     !The blob should now be centered at the origin.
         ENDDO
         DO I=0,NZM
            indG=mod(I+(NZ/2-sk),NZ)
            if (indG.lt.0)then
               indG=indG+NZ
            endif
      !      GZs(I)=GZ(indG)
             GZs(I)=GZ(I)     !The blob should now be centered at the origin.
         ENDDO

         DO J=0,NYM
         DO K=0,NZM
            IF((FLAVORS.eq.'NPZ2').and.(FLAVORI.eq.'HD')) THEN!large release
               DO I=481,NXM     !NX/2,NXM                                               
                  IF (N.eq.1) TH(I,K,J,N)=.9793d0               
                  IF (N.eq.2) TH(I,K,J,N)=.917d0
                  IF (N.EQ.3) TH(I,K,J,N)=.1702d0
               END DO
               DO I =0,480      !NX/2 -1                                                 
                  IF (N.eq.1) TH(I,K,J,N)=1.8303d0
                  IF (N.eq.2) TH(I,K,J,N)=.81303d0
                  IF (N.EQ.3) TH(I,K,J,N)=.032441d0
               END DO
            ENDIF
            IF((FLAVORS.eq.'BSTA').and.(FLAVORI.eq.'HD')) THEN!large release
               DO I=NX/2,NXM                                               
                  IF (N.eq.1) TH(I,K,J,N)=1.d0               
               END DO
               DO I =0,NX/2 -1                                                 
                  IF (N.eq.1) TH(I,K,J,N)=0.d0
               END DO
            ENDIF
            IF((FLAVORS.eq.'FKPP').and.(FLAVORI.eq.'HD')) THEN!large release
               DO I=NX/2,NXM                                               
                  IF (N.eq.1) TH(I,K,J,N)=1.d0               
               END DO
               DO I =0,NX/2 -1                                                 
                  IF (N.eq.1) TH(I,K,J,N)=0.d0
               END DO
            ENDIF
            IF (FLAVORI.eq.'PR')THEN!Patch release
               DO I=0,NXM
               rr=((GXs(I)-LX/2)**2.+(GZs(K)-LZ/2)**2.)**.5
               IF (FLAVOR.eq.'Bick')THEN!Bickley jet different coord
                  rr=((GX(I)-GX(si))**2.+(GZ(K)-GZ(sk))**2.)**.5
               ENDIF
               IF ((FLAVORS.eq.'diff').or.(FLAVORS.eq.'FKPP'))THEN
                  sigmaPR = .05
                  TH(I,K,J,N)=EXP(-rr**2.d0/(2*(sigmaPR)**2.d0)) !Gaussian, 1 scalar
               ENDIF
               IF (FLAVORS.eq.'chmo')THEN
                  IF (N.eq.1) THEN
                     TH(I,K,J,N)=EXP(-rr**2.d0*(2*PI)**2.d0)
                  ELSE
                     TH(I,K,J,N)=1.d0
                  ENDIF
               ENDIF
               IF (FLAVORS.eq.'NPZ')THEN
                  IF (N.eq.1) THEN
                  TH(I,K,J,N)=0.093745798194359
     &                 +0.124295112039992*EXP(-rr**8.d0*(6000)**2.d0)
                  ENDIF
                  IF(N.eq.2)THEN
                  TH(I,K,J,N)=0.089948486576180
     &                 -0.045933519156255*EXP(-rr**8.d0*(6000)**2.d0)
                  ENDIF
                  IF(N.eq.3)THEN
                  TH(I,K,J,N)=0.081344095155286
     &                 -0.021915144461863*EXP(-rr**8.d0*(6000)**2.d0)
                  ENDIF
               ENDIF
! Currently: Eddy_LC
               IF ((FLAVORS.eq.'NPZ2'))THEN
                  IF (N.eq.1) THEN
                  TH(I,K,J,N)= 1.83027
     &                    - 0.851*EXP(-rr**8.d0*(6000)**2.d0)
                  ENDIF
                  IF(N.eq.2)THEN
                  TH(I,K,J,N)= .813027
     &                    + 0.104*EXP(-rr**8.d0*(6000)**2.d0)
                  ENDIF
                  IF(N.eq.3)THEN
                  TH(I,K,J,N)= .032441
     &                    + 0.13776*EXP(-rr**8.d0*(6000)**2.d0)
                  ENDIF
               ENDIF
               IF ((FLAVORS.eq.'BSTA').or.(FLAVORS.eq.'FKPP'))THEN
                  IF (N.eq.1) THEN
                     TH(I,K,J,N)=EXP(-rr**2.d0*(20.d0)**2.d0)
                  ENDIF
               ENDIF
               END DO
            ENDIF
         END DO
         END DO

!--------------------------------------------------------------
         IF (0) THEN
         si=70
         sj=1
         sk=300
         WRITE(*,*) 'Blob center coordinate',si,sk
         DO I=0,NXM
            indG=mod(I+(NX/2-si),NX)
            if (indG.lt.0)then
               indG=indG+NX
            endif
            GXs(I)=GX(indG)
         ENDDO
         DO I=0,NZM
            indG=mod(I+(NZ/2-sk),NZ)
            if (indG.lt.0)then
               indG=indG+NZ
            endif
            GZs(I)=GZ(indG)
         ENDDO


         DO J=0,NYM
         DO K=0,NZM

            IF (FLAVORI.eq.'PR')THEN
               DO I=0,NXM
                  IF ((FLAVORS.eq.'NPZ2'))THEN
                     rr=((GXs(I)-LX/2)**2.+(GZs(K)-LZ/2)**2.)**.5
                     IF (N.eq.1) THEN
                        TH(I,K,J,N)= TH(I,K,J,N)
     &                       - 0.851*EXP(-rr**8.d0*(6000)**2.d0)
                     ENDIF
                     IF(N.eq.2)THEN
                        TH(I,K,J,N)= TH(I,K,J,N)
     &                       + 0.104*EXP(-rr**8.d0*(6000)**2.d0)
                     ENDIF
                     IF(N.eq.3)THEN
                        TH(I,K,J,N)= TH(I,K,J,N)
     &                       + 0.13776*EXP(-rr**8.d0*(6000)**2.d0)
                     ENDIF
                  ENDIF
                  IF ((FLAVRS.eq.'BSTA'))THEN
                     rr=((GXs(I)-LX/2)**2.+(GZs(K)-LZ/2)**2.)**.5
                     IF (N.eq.1) THEN
                     TH(I,K,J,N)=TH(I,K,J,N)+EXP(-rr**8.d0*(6000)**2.d0)
                     ENDIF
                  ENDIF
               ENDDO
            END IF

         END DO
         END DO

!--------------------------------------------------------------

         si=498
         sj=1
         sk=88
         WRITE(*,*) 'Blob center coordinate',si,sk
         DO I=0,NXM
            indG=mod(I+(NX/2-si),NX)
            if (indG.lt.0)then
               indG=indG+NX
            endif
            GXs(I)=GX(indG)
         ENDDO
         DO I=0,NZM
            indG=mod(I+(NZ/2-sk),NZ)
            if (indG.lt.0)then
               indG=indG+NZ
            endif
            GZs(I)=GZ(indG)
         ENDDO


         DO J=0,NYM
         DO K=0,NZM

            IF (FLAVORI.eq.'PR')THEN
               DO I=0,NXM
                  IF ((FLAVORS.eq.'NPZ2'))THEN
                     rr=((GXs(I)-LX/2)**2.+(GZs(K)-LZ/2)**2.)**.5
                     IF (N.eq.1) THEN
                        TH(I,K,J,N)= TH(I,K,J,N)
     &                       - 0.851*EXP(-rr**8.d0*(6000)**2.d0)
                     ENDIF
                     IF(N.eq.2)THEN
                        TH(I,K,J,N)= TH(I,K,J,N)
     &                       + 0.104*EXP(-rr**8.d0*(6000)**2.d0)
                     ENDIF
                     IF(N.eq.3)THEN
                        TH(I,K,J,N)= TH(I,K,J,N)
     &                       + 0.13776*EXP(-rr**8.d0*(6000)**2.d0)
                     ENDIF
                  ENDIF
                  IF ((FLAVORS.eq.'BSTA'))THEN
                     rr=((GXs(I)-LX/2)**2.+(GZs(K)-LZ/2)**2.)**.5
                     IF (N.eq.1) THEN
                     TH(I,K,J,N)=TH(I,K,J,N)+EXP(-rr**8.d0*(6000)**2.d0)
                     ENDIF
                  ENDIF
               ENDDO
            END IF

         END DO
         END DO

!--------------------------------------------------------------

         si=180
         sj=1
         sk=506
         WRITE(*,*) 'Blob center coordinate',si,sk
         DO I=0,NXM
            indG=mod(I+(NX/2-si),NX)
            if (indG.lt.0)then
               indG=indG+NX
            endif
            GXs(I)=GX(indG)
         ENDDO
         DO I=0,NZM
            indG=mod(I+(NZ/2-sk),NZ)
            if (indG.lt.0)then
               indG=indG+NZ
            endif
            GZs(I)=GZ(indG)
         ENDDO


         DO J=0,NYM
         DO K=0,NZM

            IF (FLAVORI.eq.'PR')THEN
               DO I=0,NXM
                  IF ((FLAVORS.eq.'NPZ2'))THEN
                     rr=((GXs(I)-LX/2)**2.+(GZs(K)-LZ/2)**2.)**.5
                     IF (N.eq.1) THEN
                        TH(I,K,J,N)= TH(I,K,J,N)
     &                       - 0.851*EXP(-rr**8.d0*(6000)**2.d0)
                     ENDIF
                     IF(N.eq.2)THEN
                        TH(I,K,J,N)= TH(I,K,J,N)
     &                       + 0.104*EXP(-rr**8.d0*(6000)**2.d0)
                     ENDIF
                     IF(N.eq.3)THEN
                        TH(I,K,J,N)= TH(I,K,J,N)
     &                       + 0.13776*EXP(-rr**8.d0*(6000)**2.d0)
                     ENDIF
                  ENDIF
                  IF ((FLAVORS.eq.'BSTA'))THEN
                     rr=((GXs(I)-LX/2)**2.+(GZs(K)-LZ/2)**2.)**.5
                     IF (N.eq.1) THEN
                     TH(I,K,J,N)=TH(I,K,J,N)+EXP(-rr**8.d0*(6000)**2.d0)
                     ENDIF
                  ENDIF
               ENDDO
            END IF

         END DO
         END DO

!--------------------------------------------------------------

         si=380
         sj=1
         sk=405
         WRITE(*,*) 'Blob center coordinate',si,sk
         DO I=0,NXM
            indG=mod(I+(NX/2-si),NX)
            if (indG.lt.0)then
               indG=indG+NX
            endif
            GXs(I)=GX(indG)
         ENDDO
         DO I=0,NZM
            indG=mod(I+(NZ/2-sk),NZ)
            if (indG.lt.0)then
               indG=indG+NZ
            endif
            GZs(I)=GZ(indG)
         ENDDO


         DO J=0,NYM
         DO K=0,NZM

            IF (FLAVORI.eq.'PR')THEN
               DO I=0,NXM
                  IF ((FLAVORS.eq.'NPZ2'))THEN
                     rr=((GXs(I)-LX/2)**2.+(GZs(K)-LZ/2)**2.)**.5
                     IF (N.eq.1) THEN
                        TH(I,K,J,N)= TH(I,K,J,N)
     &                       - 0.851*EXP(-rr**8.d0*(6000)**2.d0)
                     ENDIF
                     IF(N.eq.2)THEN
                        TH(I,K,J,N)= TH(I,K,J,N)
     &                       + 0.104*EXP(-rr**8.d0*(6000)**2.d0)
                     ENDIF
                     IF(N.eq.3)THEN
                        TH(I,K,J,N)= TH(I,K,J,N)
     &                       + 0.13776*EXP(-rr**8.d0*(6000)**2.d0)
                     ENDIF
                  ENDIF
                  IF ((FLAVORS.eq.'BSTA'))THEN
                     rr=((GXs(I)-LX/2)**2.+(GZs(K)-LZ/2)**2.)**.5
                     IF (N.eq.1) THEN
                     TH(I,K,J,N)=TH(I,K,J,N)+EXP(-rr**8.d0*(6000)**2.d0)
                     ENDIF
                  ENDIF
               ENDDO
            END IF

         END DO
         END DO
      END IF !Single patch

!      write(*,*)TH(316,198,0,1)

!     FFT to Fourier newly created scalar field
         CALL FFT_XZY_TO_FOURIER(TH(0,0,0,N),CTH(0,0,0,N)) 
!     Endif create new scalar field
       END IF
!     Loop scalar species
       END DO

       RETURN
       END
      

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_GRID_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER I,J,K

         IF (FLAVOR.NE.'Ensemble') WRITE (6,*) 'Fourier in X, Z, Y'
         DO I=0,NX
           GX(I)=(I*LX)/NX
           DX(I)=LX/NX
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GX(',I,') = ',GX(I)
         END DO
!         IF (FLAVOR.NE.'Ensemble') WRITE (6,*) 'Fourier in Z'
         DO K=0,NZ
           GZ(K)=(K*LZ)/NZ
           DZ(K)=LZ/NZ
           IF (FLAVOR.EQ.'Bick') GZ(K)=GZ(K)-3.
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GZ(',K,') = ',GZ(K)
         END DO
!         IF (FLAVOR.NE.'Ensemble') WRITE (6,*) 'Fourier in Y'
         DO J=0,NY
           GY(J)=(J*LY)/NY
           DY(J)=LY/NY
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GY(',J,') = ',GY(J)
         END DO

         RETURN
         END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INPUT_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      REAL    VERSION, CURRENT_VERSION

! Read in input parameters specific for channel flow case
      OPEN (11,file='input_per.dat',form='formatted',
     &      status='old')
C Read input file.

      CURRENT_VERSION=1.0
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*) VERSION
      IF (VERSION .NE. CURRENT_VERSION)
     &         STOP 'Wrong input data format input_chan'
      READ(11,*)
      READ(11,*) TIME_AD_METH
      READ(11,*)
      READ(11,*) LES_MODEL_TYPE
      READ(11,*)
      READ(11,*) IC_TYPE, KICK
      READ(11,*)
      READ(11,*) TRUCK, GUSTS
      DO N=1,N_TH
        READ(11,*)
        READ(11,*) BACKGROUND_GRAD(N)
        READ(11,*)
        READ(11,*) CHI(N), C0(N), V_S(N)
      END DO

      RETURN
      END



C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE VIS_FLOW_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_STATS_PER(FINAL)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      CHARACTER*35 FNAME
      LOGICAL FINAL
          real*8 mean_u1(0:NYM)
          real*8 mean_u2(0:NYM)
          real*8 mean_u3(0:NYM)
          real*8 mean_th(0:NYM,1:N_TH), thbarall(1:N_TH)
          real*8 mean_p(0:NYM)


      integer i,j,k,n
      real*8 uc, ubulk, integL(1:N_TH), EC
	
	IF ((FLAVOR.EQ.'Basic').or.(FLAVOR.eq.'CHEMOTAXIS')) THEN

! Note that this routine uses CFi and CFTH for storage, so it should
! only be used between full R-K timesteps

      WRITE(6,*) 'Saving flow statistics.'

      if (FINAL) then
! We are done with the simulation
! Write out statistics to a file
        open(20,file='stats.txt',form='formatted',status='unknown')
        do j=0,NYM
          write(20,201) j,GYF(j),UBAR(j),VBAR(j),WBAR(j)
        end do
201     format(I3,',',F16.9,',',F16.9,',',F16.9,',',F16.9)
        do n=1,N_TH
        do j=0,NYM
          write(20,202) j,GYF(j),THBAR(j,n)
        end do
        end do
202     format(I3,',',F16.9,',',F16.9)

        else
! We are in the middle of a run, compile statistics

! Store the velocity in Fourier space in CRi
      do j=0,TNKY
        do k=0,TNKZ
          do i=0,NKX
            CR1(i,k,j)=CU1(i,k,j)
            CR2(i,k,j)=CU2(i,k,j)
            CR3(i,k,j)=CU3(i,k,j)
          end do
        end do
      end do

! Compute the vertical gradients in fourier space, store in CFi
      do j=0,TNKY
        do k=0,TNKZ
          do i=0,NKX
            CF1(i,k,j)=CIKY(j)*CU1(i,k,j)
            CF2(i,k,j)=CIKY(j)*CU2(i,k,j)
            CF3(i,k,j)=CIKY(j)*CU3(i,k,j)
          end do
        end do
      end do
! Save the  scalar in fourier space in CFTH
      do n=1,N_TH
        do j=0,TNKY
          do k=0,TNKZ
            do i=0,NKX
              CFTH(i,k,j,n)=CTH(i,k,j,n)
            end do
          end do
        end do
      end do







! Now, convert the velocity and vertical gradients to physical space
      CALL FFT_XZY_TO_PHYSICAL(CU1,U1)
      CALL FFT_XZY_TO_PHYSICAL(CU2,U2)
      CALL FFT_XZY_TO_PHYSICAL(CU3,U3)
      do n=1,N_TH
        CALL FFT_XZY_TO_PHYSICAL(CTH(0,0,0,n),TH(0,0,0,n))
      end do
      CALL FFT_XZY_TO_PHYSICAL(CF1,F1)
      CALL FFT_XZY_TO_PHYSICAL(CF2,F2)
      CALL FFT_XZY_TO_PHYSICAL(CF3,F3)

!          write(*,*) 'th1', TH(128,128,0,1)
!          write(*,*) 'th2', TH(128,128,0,2)
!          write(*,*) 'th3', TH(128,128,0,3)
! First get the number of samples taken so far
      NSAMPLES=NSAMPLES+1
! Get the mean velocity

! First, get the mean in physical space
      do j=0,NYM
        mean_u1(j)=0.d0
        mean_u2(j)=0.d0
        mean_u3(j)=0.d0
        do n=1,N_TH
          mean_th(j,n)=0.d0
        end do
        mean_p(j)=0.d0
        do i=0,NXM
          do k=0,NZM
            mean_u1(j)=mean_u1(j)+U1(i,k,j)
            mean_u2(j)=mean_u2(j)+U2(i,k,j)
            mean_u3(j)=mean_u3(j)+U3(i,k,j)
            do n=1,N_TH
              mean_th(j,n)=mean_th(j,n)+TH(i,k,j,n)
            end do
            mean_p(j)=mean_p(j)+P(i,k,j)
          end do
        end do
        mean_u1(j)=mean_u1(j)/dble(NX*NZ)
        mean_u2(j)=mean_u2(j)/dble(NX*NZ)
        mean_u3(j)=mean_u3(j)/dble(NX*NZ)
         do n=1,N_TH
          mean_th(j,n)=mean_th(j,n)/dble(NX*NZ)
        end do
        mean_p(j)=mean_p(j)/dble(NX*NZ)
      end do
     
      do j=0,NYM
        UBAR(j)=(1./float(NSAMPLES))*mean_u1(j)
     &      +((float(NSAMPLES)-1.)/float(NSAMPLES))*UBAR(j)
        VBAR(j)=(1./float(NSAMPLES))*mean_u2(j)
     &      +((float(NSAMPLES)-1.)/float(NSAMPLES))*VBAR(j)
        WBAR(j)=(1./float(NSAMPLES))*mean_u3(J)
     &      +((float(NSAMPLES)-1.)/float(NSAMPLES))*WBAR(j)
        do n=1,N_TH
          THBAR(j,n)=(1./float(NSAMPLES))*mean_th(j,n)
     &      +((float(NSAMPLES)-1.)/float(NSAMPLES))*THBAR(j,n)
        end do
      end do

      do n=1,N_TH
         thbarall(n)=0.d0
         do j=0,NYM
            thbarall(n)=thbarall(n)+mean_th(j,n)/dble(NY)
         end do
      end do

! Get the turbulent kinetic energy at each level 
      do j=0,NYM
        urms(j)=0.
        vrms(j)=0.
        wrms(j)=0.
      do k=0,NZM
      do i=0,NXM 
        urms(j)=urms(j)+(abs(U1(i,k,j)-mean_u1(j)))**2.
        vrms(j)=vrms(j)+(abs(U2(i,k,j)-mean_u2(j)))**2.
        wrms(j)=wrms(j)+(abs(U3(i,k,j)-mean_u3(j)))**2.
      end do
      end do
        urms(j)=sqrt(urms(j)/(float(NZ)*float(NX)))
        vrms(j)=sqrt(vrms(j)/(float(NZ)*float(NX)))
        wrms(j)=sqrt(wrms(j)/(float(NZ)*float(NX)))
      end do 
! Get the bulk rms value
      urms_b=0.
      do j=1,NYM
        urms_b=urms_b+0.5*(urms(j)+urms(j-1))*(GY(j)-GY(j-1))
      end do
      urms_b=urms_b/LY
! If we are in 2d:
      if (NY.eq.1) urms_b=urms(0)
       
! Compute the Reynolds stress and mean velocity gradient
      do j=0,NYM
        uv(j)=0. 
        uw(j)=0.
        wv(j)=0.
      do k=0,NZM
      do i=0,NXM
        uv(j)=uv(j)+(U1(i,k,j)-mean_u1(j))
     +    *(U2(i,k,j)-mean_u2(j))
        wv(j)=wv(j)+(U3(i,k,j)-mean_u3(j))
     +    *(U2(i,k,j)-mean_u2(j))
        uw(j)=uw(j)+(U1(i,k,j)-mean_u1(j))
     +    *(U3(i,k,j)-mean_u3(j))
      end do
      end do
        uv(j)=uv(j)/(float(NZ)*float(NX))
        uw(j)=uw(j)/(float(NZ)*float(NX))
        wv(j)=wv(j)/(float(NZ)*float(NX))
      end do
              
! Get the y-derivative of the mean velocity at GYF points
      do j=1,NY-2
        dudy(j)=(mean_u1(j+1)-mean_u1(j-1))/(GY(j+1)-GY(j-1))
        dwdy(j)=(mean_u3(j+1)-mean_u3(j-1))/(GY(j+1)-GY(j-1))
      end do
      j=0
        dudy(j)=(mean_u1(j+1)-mean_u1(NYM))/(2.d0*(GY(j+1)-GY(j)))
        dwdy(j)=(mean_u3(j+1)-mean_u3(NYM))/(2.d0*(GY(j+1)-GY(j)))
      j=NYM
        dudy(j)=(mean_u1(0)-mean_u1(j-1))/(2.d0*(GY(j)-GY(j-1)))
        dwdy(j)=(mean_u3(0)-mean_u3(j-1))/(2.d0*(GY(j)-GY(j-1)))

! Get the mean square shear
      do j=0,NYM
        shear(j)=0.d0
        do k=0,NZM
          do i=0,NXM
            shear(j)=shear(j)+F1(i,k,j)**2.d0+F3(i,k,j)**2.d0
          end do
        end do
        shear(j)=shear(j)/dble(NX*NZ)
      end do

! Write out the bulk rms velocity
!      write(*,*) '<U_rms>: ',urms_b
!      write(*,*) 'VERBOSITY: ',VERBOSITY

! Write out the mean statistics at each time
      open(40,file='mean.txt',form='formatted',status='unknown')
      write(40,*) TIME_STEP,TIME,DELTA_T,UBULK
      do j=0,NYM
        write(40,401) j,GY(J),mean_u1(j)
     +      ,mean_u2(j)
     +      ,mean_u3(j),urms(j),vrms(j),wrms(j)
     +      ,uv(j),uw(j),wv(j),dudy(j),dwdy(j),mean_p(j),shear(j)
      end do

401   format(I3,' ',14(F20.9,' '))



      if (MOVIE) then
! Output a 2d slice through the velocity field for animation in matlab
        open(79,file='movie_vel_xz.txt',status='unknown'
     &      ,form='formatted')
        do i=0,NXM
        do k=0,NZM
          write(79,*) U1(I,K,NY/2)**2.+U3(I,K,NY/2)**2.+U2(I,K,NY/2)**2.
        end do
        end do

        open(80,file='movie_vel_xy.txt',status='unknown'
     &        ,form='formatted')
        do i=0,NXM
        do j=0,NYM
          write(80,*) U1(I,NZ/2,J)**2.+U3(I,NZ/2,J)**2.+U2(I,NZ/2,J)**2.
        end do
        end do

! This file will contain a single plane and is used in conjunction with
! the matlab script 'realtime_movie' to visualize data during
! simulation
        open (76,file='temp.txt',status='unknown',form='formatted')
        do K=0,NZM
          write(76,*) gz(k)
        end do
        do I=0,NXM
        do K=0,NZM
          write(76,*) U1(I,K,NY/4)**2.+U3(I,K,NY/4)**2.+U2(I,K,NY/4)**2.
        end do
        end do
        close (76)
        CALL SYSTEM('mv temp.txt ../post_process/matlab/latest_slice.txt
     &')

      end if 

! Do over the number of passive scalars
      do n=1,N_TH
!        write(*,*) 'ths = ', TH(128,128,0,n)
       do j=0,NYM
         thrms(j,n)=0.
       do k=0,NZM
       do i=0,NXM
        thrms(j,n)=thrms(j,n)+(abs(TH(i,k,j,n)-mean_th(j,n)))**2.
       end do
       end do
        thrms(j,n)=sqrt(thrms(j,n)/(float(NZ)*float(NX)))
       end do
! Compute the Reynolds stress and mean velocity gradient
       do j=0,NYM
         thv(j,n)=0.
       do k=0,NZM
       do i=0,NXM
        thv(j,n)=thv(j,n)+(TH(i,k,j,n)-mean_th(j,n))
     &           *(U2(i,k,n)-mean_u2(j))
       end do
       end do
        thv(j,n)=thv(j,n)/(float(NZ)*float(NX))
       end do

! Get the y-derivative of the mean velocity at GYF points
      do j=2,NY
        dthdy(j,n)=(CRTH(0,0,j+1,n)-CRTH(0,0,j-1,n))/(2.*DYF(j))
      end do
      do j=1,NY-2
        dthdy(j,n)=(mean_th(j+1,n)-mean_th(j-1,n))/(GY(j+1)-GY(j-1))
      end do
      j=0
       dthdy(j,n)=(mean_th(j+1,n)-mean_th(NYM,n))/(2.d0*(GY(j+1)-GY(j)))
      j=NYM
       dthdy(j,n)=(mean_th(0,n)-mean_th(j-1,n))/(2.d0*(GY(j)-GY(j-1)))

! Get the bacterial/nutrient correlation then compute integral length
!       write(*,*) 'th', n, TH(128,128,0,n)
!       write(*,*) 'thbarall', n,  thbarall(n)
      thth(n)=0.d0
      thvar(n)=0.d0     
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        thvar(n)=thvar(n)+(TH(i,k,j,n)-thbarall(n))
     &            *(TH(i,k,j,n)-thbarall(n))*DX(I)*DY(J)*DZ(K)
      end do
      end do
      end do
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        thth(n)=thth(n)+(TH(i,k,j,n)-thbarall(n))
     &            *(TH(i,k,j,1)-thbarall(1))*DX(I)*DY(J)*DZ(K)
     &            /sqrt(thvar(n)*thvar(1))
 
      end do
      end do
      end do
      if (thvar(n).EQ.0.d0) then
         thth(n)=0.d0
      end if
      
      EC=0.d0
      do j=0,TNKY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=CIKX(i)*CFTH(i,k,j,n)
            EC=EC+CS1(I,K,J)*CONJG(CS1(I,K,J))*DX(I)*DY(J)*DZ(K)
            CS1(i,k,j)=CIKY(j)*CFTH(i,k,j,n)
            EC=EC+CS1(I,K,J)*CONJG(CS1(I,K,J))*DX(I)*DY(J)*DZ(K)
            CS1(i,k,j)=CIKZ(k)*CFTH(i,k,j,n)
            EC=EC+CS1(I,K,J)*CONJG(CS1(I,K,J))*DX(I)*DY(J)*DZ(K)
          end do
        end do
      end do
!      integL(n)=(thvar(n)/EC)**0.5d0


      write(*,*) 'n,thth(n),thvar(n),integL(n): ',n,thth(n),thvar(n),
     &           integL(n)

      if (MOVIE) then
! Output a 2d slice through the scalar field for animation in matlab
        if (n.eq.1) then
! Chose which scalar is to be outputted
        open(75,file='movie1_xy.txt',status='unknown',form='formatted')
        do i=0,NXM
        do j=0,NYM
          write(75,*) TH(I,NZ/2,J,n)
        end do
        end do
        open(77,file='movie1_xz.txt',status='unknown',form='formatted')
        do i=0,NXM
        do k=0,NZM
          write(77,*) TH(I,K,NY/2,n)
        end do
        end do
        end if
        if (n.eq.2) then
! Chose which scalar is to be outputted
        open(85,file='movie2_xy.txt',status='unknown',form='formatted')
         
        do i=0,NXM
        do j=0,NYM
          write(85,*) TH(I,NZ/2,J,n)
        end do
        end do
        open(87,file='movie2_xz.txt',status='unknown',form='formatted')
        do i=0,NXM
        do k=0,NZM
          write(87,*) TH(I,K,NY/2,n)
        end do
        end do
        end if
        if (n.eq.3) then
! Chose which scalar is to be outputted
        open(95,file='movie3_xy.txt',status='unknown',form='formatted')
        do i=0,NXM
        do j=0,NYM
          write(95,*) TH(I,NZ/2,J,n)
        end do
        end do
        open(97,file='movie3_xz.txt',status='unknown',form='formatted')
        do i=0,NXM
        do k=0,NZM
          write(97,*) TH(I,K,NY/2,n)
        end do
        end do
        end if
      end if 


! End do over number of passive scalars, n
      end do

! Convert back to Fourier space
      do n=1,N_TH
        call fft_xzy_to_fourier(TH(0,0,0,n),CTH(0,0,0,n))
      end do

! Write out the mean statistics at each time
      open(41,file='mean_th.txt',form='formatted',status='unknown')
      write(41,*) TIME_STEP,TIME,DELTA_T,UBULK
      do n=1,N_TH 
      do j=0,NYM
        write(41,402) j,GYF(J),mean_th(j,n)
     &      ,dthdy(j,n),thrms(j,n),thv(j,n),pe_diss(j,n)
     &      ,thth(n),thvar(n)
      end do
      end do

402   format(I3,' ',8(F30.25,' '))

!      write(*,*) 'VERBOSITY: ',VERBOSITY
      if (VERBOSITY.gt.4) then 
      write(*,*) 'Outputting info for gnuplot...'
      open (unit=10, file="solution")
      do i=2,NXM
        do j=2,NYM
          write (10,*) i, j, U1(i,0,j)
        end do
        write (10,*) ""
      end do
      close (10)
      call system ('gnuplot <gnuplot.in') 
      end if

C Convert velocity back to Fourier space
      call fft_xzy_to_fourier(U1,CU1)
      call fft_xzy_to_fourier(U2,CU2)
      call fft_xzy_to_fourier(U3,CU3)

      end if

 
      call tkebudget_per
	
	ELSE IF(FLAVOR.EQ.'Ensemble') THEN
	
	call save_stats_ensem(final)
	
	END IF


      RETURN
      END



C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FILTER_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This subroutine applies a filter to the highest wavenumbers
C It should be applied to the scalars in Fourier space
C The filter used is a sharpened raised cosine filter

      include 'header'

      integer I,J,K,N

! Variables for horizontal filtering
      real*8 sigma(0:NKX,0:TNKZ,0:TNKY),sigma0

C Set the filtering constants for the all directions
      DO i=0,NKX
       DO k=0,TNKZ
        DO j=0,TNKY
          sigma0=0.5d0*(1.d0+
     &       cos(sqrt((KX(i)*LX*1.d0/float(NX))**2.d0
     &            +(KZ(k)*LZ*1.d0/float(NZ))**2.d0
     &            +(KY(j)*LY*1.d0/float(NY))**2.d0)))
! Apply a sharpened raised cosine filter
         sigma(i,k,j)=sigma0**4.d0*(35.d0-84.d0*sigma0
     &        +70.d0*sigma0**2.d0-20.d0*sigma0**3.d0)
        END DO
       END DO
      END DO

C Apply the spectral filtering 
 
      DO N=1,N_TH
        DO J=0,TNKY
          DO K=0,TNKZ
            DO I=0,NKX
              CTH(I,K,J,N)=CTH(I,K,J,N)*sigma(i,k,j)
            END DO
          END DO
        END DO
      END DO

       return
       end





C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE TKEBUDGET_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! NOte, it is important to only run this routine after complete R-K
!  time advancement since F1 is overwritten which is needed between R-K steps

      include 'header'

      integer i,j,k

! Compute the turbulent dissipation rate, epsilon=nu*<du_i/dx_j du_i/dx_j>
      do j=0,NYM
        epsilon(j)=0.
      end do
! Store du/dx in CS1
      do j=0,TNKY
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKX(i)*CU1(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xzy_to_physical(CS1,S1)
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Store dv/dx in CS1
      do j=0,TNKY
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKX(i)*CU2(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xzy_to_physical(CS1,S1)
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Compute du/dy at GYF gridpoints, note remove mean
      do j=0,TNKY
      do k=0,TNKZ
      do i=0,NKX
        CF1(i,k,j)=CIKY(j)*CU1(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xzy_to_physical(CF1,F1)
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(F1(i,k,j)**2.0)
! Cross term dvdx*dudy
        epsilon(j)=epsilon(j)+(S1(i,k,j)*F1(i,k,j))
      end do
      end do
      end do
! Store dw/dx in CS1
      do j=0,TNKY
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKX(i)*CU3(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xzy_to_physical(CS1,S1)
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Compute du/dz at GYF gridpoints, note remove mean
! Store du/dz in CS1
      do j=0,TNKY
      do k=0,TNKZ
      do i=0,NKX
        CF1(i,k,j)=CIKZ(k)*CU1(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xzy_to_physical(CF1,F1)
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(F1(i,k,j)**2.0)
! Cross term dudz*dwdx
        epsilon(j)=epsilon(j)+S1(i,k,j)*F1(i,k,j)
      end do
      end do
      end do
! Compute dv/dy at GYF gridpoints, note remove mean
      do j=0,TNKY
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKY(j)*CU2(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xzy_to_physical(CS1,S1)
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Compute dw/dy at GYF gridpoints, note remove mean
      do j=0,TNKY
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKY(j)*CU3(i,k,j)
      end do
      end do
      end do
! Convert to physical space 
      call fft_xzy_to_physical(CS1,S1)
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Store dv/dz in CF1
      do j=0,TNKY
      do k=0,TNKZ
      do i=0,NKX
        CF1(i,k,j)=CIKZ(k)*CU2(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xzy_to_physical(CF1,F1)
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(F1(i,k,j)**2.0)
! Cross term dvdz*dwdy
        epsilon(j)=epsilon(j)+S1(i,k,j)*F1(i,k,j)
      end do
      end do
      end do
! Store dw/dz in CS1
      do j=0,TNKY
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKZ(k)*CU3(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xzy_to_physical(CS1,S1)
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
      do j=0,NYM
        epsilon(j)=epsilon(j)/float(NX*NZ)
      end do


! Write out the bulk rms velocity
!      write(*,*) '<U_rms>: ',urms_b


! Write out the mean statistics at each time
      open(45,file='tke.txt',form='formatted',status='unknown')
      write(45,*) TIME_STEP,TIME,DELTA_T
      do j=0,NYM
        write(45,401) j,GY(J),epsilon(j)
      end do
401   format(I3,' ',2(F20.9,' '))


! Get Kolmogorov wavelength
      epsilon_mean=NU*SUM(epsilon(0:NYM))/dble(NY) 

      k_eta=2.d0*PI*(NU**3.d0/epsilon_mean)**(-0.25d0)

      write(*,*) 'Kolmogorov scale: ',2.d0*PI/k_eta

      return
      end



      SUBROUTINE Interpolation(velo, V)
      include 'header'
      INTEGER    I, K, L, n1, n3
      REAL*8     resx, resz, LenX,LenZ,ddx,ddz
      REAL*8     f000, f100, f010,f001
      INTEGER    indx1, indx2, indz1, indz2
      REAL*8     velo(0:Nx-1,0:Nz-1)
      REAL*8     xpos1(0:Nx-1,0:Nz-1),zpos1(0:Nx-1,0:Nz-1)
      REAL*8     V(0:Nx+1,0:Nz+1,0:NY+1)
      REAL*8     xs1(0:Nx+1),zs1(0:Nz+1)

         xpos1=xpos
         zpos1=zpos
         xs1=Gx
         zs1=Gz
         ddx=Gx(2)-Gx(1)
         ddz=Gz(2)-Gz(1)
         n1=Nx
         n3=Nz
         LenX=2*3.1415926
         LenZ=LenX
         DO i = 0, n1-1
            DO K = 0, n3-1
               IF(xpos(I,K) .GT. LenX .OR. xpos(I,K) .LT. 0.0) THEN
                  xpos1(I,K)=MOD(xpos(I,K)+20.0*LenX,LenX)
               ENDIF
               IF(zpos(I,K) .GT. LenZ .OR. zpos(I,K) .LT. 0.0) THEN
                  zpos1(I,K)=MOD(zpos(I,K)+20.0*LenZ,LenZ)
               ENDIF
               L = 1
               DO WHILE (xpos1(I,K)-xs1(L) .GE. 0 .AND. L .LE. Nx)
                  L = L + 1
               ENDDO
               indx1 = L - 1
               indx2 = L
               IF (indx2 .EQ. Nx+1) THEN
                  indx2 = 1
               ENDIF
               resx = (xpos1(I,K)-xs1(indx1))/ddx
               L = 1
               DO WHILE (zpos1(I,K)-zs1(L) .GE. 0 .AND. L .LE. Nz)
                  L = L + 1
               ENDDO
               indz1 = L - 1
               indz2 = L
               IF (indz2 .EQ. Nz+1) THEN
                  indz2 = 1
               ENDIF
               resz = (zpos1(I,K)-zs1(indz1))/ddz
               f000 = V(indx1,indz1,0)
               f100 = V(indx2,indz1,0)
               f010 = V(indx1,indz2,0)
               f001 = V(indx1,indz1,0)
               velo(I,K) = f000*(1-resx)*(1-resz) 
     &              + f100*resx*(1-resz) 
     &              + f010*(1-resx)*resz 
     &              + f001*resx*resz 
            ENDDO
         ENDDO
         RETURN
         end

