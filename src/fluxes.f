        SUBROUTINE  FLUXES(CMU, CWT, FBEAM, GC, KK, LAYRU, LL, LYRCUT,  &
     &                    MAXULV, MXULV, NCUT, NN, NSTR, NTAU,          &
     &                    PI, PRNT, SSALB, TAUCPR, UMU0, UTAU, UTAUPR,  &
     &                    XR0, XR1, ZZ, ZPLK0, ZPLK1, DFDT, FLUP,       &
     &                    FLDN, FLDIR, RFLDIR, RFLDN, UAVG, U0C)

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!       CALCULATES THE RADIATIVE FLUXES, MEAN INTENSITY, AND FLUX
!       DERIVATIVE WITH RESPECT TO OPTICAL DEPTH FROM THE M=0 INTENSITY
!       COMPONENTS (THE AZIMUTHALLY-AVERAGED INTENSITY)

!                   I N P U T     V A R I A B L E S:

!       CMU      :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE
!       CWT      :  WEIGHTS FOR GAUSS QUADRATURE OVER ANGLE COSINE
!       GC       :  EIGENVECTORS AT POLAR QUADRATURE ANGLES, SC(1)
!       KK       :  EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)
!       LAYRU    :  LAYER NUMBER OF USER LEVEL -UTAU-
!       LL       :  CONSTANTS OF INTEGRATION IN EQ. SC(1), OBTAINED
!                   BY SOLVING SCALED VERSION OF EQ. SC(5);
!                   EXPONENTIAL TERM OF EQ. SC(12) NOT INCLUDED
!       LYRCUT   :  LOGICAL FLAG FOR TRUNCATION OF COMPUT. LAYER
!       NN       :  ORDER OF DOUBLE-GAUSS QUADRATURE (NSTR/2)
!       NCUT     :  NUMBER OF COMPUTATIONAL LAYER WHERE ABSORPTION
!                     OPTICAL DEPTH EXCEEDS -ABSCUT-
!       TAUCPR   :  CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED)
!       UTAUPR   :  OPTICAL DEPTHS OF USER OUTPUT LEVELS IN DELTA-M
!                     COORDINATES;  EQUAL TO  -UTAU- IF NO DELTA-M
!       XR0      :  EXPANSION OF THERMAL SOURCE FUNCTION IN EQ. SS(14)
!       XR1      :  EXPANSION OF THERMAL SOURCE FUNCTION EQS. SS(16)
!       ZZ       :  BEAM SOURCE VECTORS IN EQ. SS(19)
!       ZPLK0    :  THERMAL SOURCE VECTORS -Z0-, BY SOLVING EQ. SS(16)
!       ZPLK1    :  THERMAL SOURCE VECTORS -Z1-, BY SOLVING EQ. SS(16)
!       (REMAINDER ARE 'DISORT' INPUT VARIABLES)

!                   O U T P U T     V A R I A B L E S:

!       U0C      :  AZIMUTHALLY AVERAGED INTENSITIES
!                   (AT POLAR QUADRATURE ANGLES)
!       (RFLDIR, RFLDN, FLUP, DFDT, UAVG ARE 'DISORT' OUTPUT VARIABLES)

!                   I N T E R N A L       V A R I A B L E S:

!       DIRINT   :  DIRECT INTENSITY ATTENUATED
!       FDNTOT   :  TOTAL DOWNWARD FLUX (DIRECT + DIFFUSE)
!       FLDIR    :  DIRECT-BEAM FLUX (DELTA-M SCALED)
!       FLDN     :  DIFFUSE DOWN-FLUX (DELTA-M SCALED)
!       FNET     :  NET FLUX (TOTAL-DOWN - DIFFUSE-UP)
!       FACT     :  EXP(- UTAUPR / UMU0)
!       PLSORC   :  PLANCK SOURCE FUNCTION (THERMAL)
!       ZINT     :  INTENSITY OF m=0 CASE, IN EQ. SC(1)
!+---------------------------------------------------------------------+

        LOGICAL LYRCUT, PRNT(*)
        REAL    DFDT(*), FLUP(*), FLDIR(*), FLDN(*),RFLDIR(*),RFLDN(*), &
     &        U0C(MXCMU,MXULV), UAVG(*)
        INTEGER LAYRU(*)
        REAL    CMU(*), CWT(*), GC(MXCMU,MXCMU,*), KK(MXCMU,*),         &
     &        LL(MXCMU,*), SSALB(*), TAUCPR(0:*),                       &
     &        UTAU(*), UTAUPR(*), XR0(*), XR1(*), ZZ(MXCMU,*),          &
     &        ZPLK0(MXCMU,*), ZPLK1(MXCMU,*)
!MOD
!MOD  BEGIN MODTRAN INSERT:
!MOD
!MOD  COMMONS:
      INCLUDE 'IFIL.h'

!     /LSTART/
!       TOLABS   OUTPUT WARNING THRESHOLD FOR EXCESS ABSORPTION.
!       LLORMN   OUTPUT WARNING FLAG FOR SMALL LORENTZ HALF WIDTH.
!       LLORMX   OUTPUT WARNING FLAG FOR LARGE LORENTZ HALF WIDTH.
!       LDOPMN   OUTPUT WARNING FLAG FOR SMALL DOPPLER HALF WIDTH.
!       LDOPMX   OUTPUT WARNING FLAG FOR LARGE DOPPLER HALF WIDTH.
!       LLINMN   OUTPUT WARNING FLAG FOR SMALL EFFECTIVE LINE NUMBER.
!       LSSA     OUTPUT WARNING FLAG FOR SCATTERING ALBEDO NEAR 1.
!       LAEREX   OUTPUT WARNING FLAG FOR ROUTINE AEREXT.
!       LBMCRK   OUTPUT WARNING FLAG FOR ROUTINE BMCRKS.
!       LFLUXS   OUTPUT WARNING FLAG FOR ROUTINE FLUXES.
!       LSSRAD   OUTPUT WARNING FLAG FOR ROUTINE SSRAD.
!       LTRLAY   OUTPUT WARNING FLAG FOR ROUTINE TRLAY.
!       LNVERS   OUTPUT WARNING FLAG FOR ROUTINE DENFAC.
!       LGEOM    OUTPUT WARNING FLAG FOR INDEX OF REFRACTION GRADIENT.
!       LO3TRN   OUTPUT WARNING FLAG FOR OZONE CURTIS-GODSON PROBLEM.
      REAL TOLABS
      LOGICAL LLORMN,LLORMX,LDOPMN,LDOPMX,LLINMN,LSSA,                  &
     &  LAEREX,LBMCRK,LFLUXS,LSSRAD,LTRLAY,LNVERS,LGEOM,LO3TRN
      COMMON/LSTART/TOLABS,LLORMN,LLORMX,LDOPMN,LDOPMX,LLINMN,LSSA,     &
     &  LAEREX,LBMCRK,LFLUXS,LSSRAD,LTRLAY,LNVERS,LGEOM,LO3TRN
      SAVE /LSTART/
!MOD
!MOD    END MODTRAN INSERT:

      IF(.NOT.LJMASS .AND. PRNT(2))WRITE(*,'(//24X,A,/(A))')            &
     &   '<----------------------- FLUXES ----------------------->',    &
     &   '      OPTICAL  COMPU    DOWNWARD    DOWNWARD    DOWNWARD     '&
     &   //' UPWARD                    MEAN      PLANCK   D(NET FLUX)', &
     &   '        DEPTH  LAYER      DIRECT     DIFFUSE       TOTAL     '&
     &   //'DIFFUSE         NET   INTENSITY      SOURCE   / D(OP DEP)'
!                                          ** ZERO DISORT OUTPUT ARRAYS
        CALL  ZEROIT(U0C(1,1), MXULV*MXCMU)
        CALL  ZEROIT(RFLDIR(1), MAXULV)
        CALL  ZEROIT(FLDIR(1),  MAXULV)
        CALL  ZEROIT(RFLDN(1),  MAXULV)
        CALL  ZEROIT(FLDN(1),   MAXULV)
        CALL  ZEROIT(FLUP(1),   MAXULV)
        CALL  ZEROIT(UAVG(1),   MAXULV)
        CALL  ZEROIT(DFDT(1),   MAXULV)
!                                        ** LOOP OVER USER LEVELS
        DO LU=1, NTAU

           LYU=LAYRU(LU)

           IF (LYRCUT .AND. LYU.GT.NCUT) THEN
!                                                ** NO RADIATION REACHES
!                                                ** THIS LEVEL
              FDNTOT=0.0
              FNET  =0.0
              PLSORC=0.0
              GO TO 90
           END IF

           IF (FBEAM.GT.0.0)  THEN
              FACT =EXP(- UTAUPR(LU) / UMU0)
              DIRINT=FBEAM * FACT
              FLDIR( LU)=UMU0 * DIRINT
              RFLDIR(LU)=UMU0 * FBEAM * EXP(- UTAU(LU) / UMU0)
           ELSE
              DIRINT=0.0
              FLDIR( LU)=0.0
              RFLDIR(LU)=0.0
           END IF

           DO IQ=1, NN

              ZINT=0.0
              DO JQ=1, NN
                 IF(ABS(GC(IQ,JQ,LYU)).GT.1E-30)                        &
     &             ZINT=ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *             &
     &             EXP(- KK(JQ,LYU) * (UTAUPR(LU) - TAUCPR(LYU)))
              ENDDO
              DO JQ=NN+1, NSTR
                 IF(ABS(GC(IQ,JQ,LYU)).GT.1E-30)                        &
     &             ZINT=ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *             &
     &             EXP(- KK(JQ,LYU) * (UTAUPR(LU) - TAUCPR(LYU-1)))
              ENDDO

              U0C(IQ,LU)=ZINT
              IF (FBEAM.GT.0.0)  U0C(IQ,LU)=ZINT + ZZ(IQ,LYU) * FACT
              U0C(IQ,LU)=U0C(IQ,LU) + ZPLK0(IQ,LYU)                     &
     &                     + ZPLK1(IQ,LYU) * UTAUPR(LU)
              UAVG(LU)=UAVG(LU) + CWT(NN+1-IQ) * U0C(IQ,LU)
              FLDN(LU)=FLDN(LU) + CWT(NN+1-IQ)*CMU(NN+1-IQ) * U0C(IQ,LU)
           ENDDO
           DO IQ=NN+1, NSTR

              ZINT=0.0
              DO JQ=1, NN
                 IF(ABS(GC(IQ,JQ,LYU)).GT.1E-30)                        &
     &             ZINT=ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *             &
     &             EXP(- KK(JQ,LYU) * (UTAUPR(LU) - TAUCPR(LYU)))
              ENDDO
              DO JQ=NN+1, NSTR
                 IF(ABS(GC(IQ,JQ,LYU)).GT.1E-30)                        &
     &             ZINT=ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *             &
     &             EXP(- KK(JQ,LYU) * (UTAUPR(LU) - TAUCPR(LYU-1)))
              ENDDO

              U0C(IQ,LU)=ZINT
              IF (FBEAM.GT.0.0)  U0C(IQ,LU)=ZINT + ZZ(IQ,LYU) * FACT
              U0C(IQ,LU)=U0C(IQ,LU) + ZPLK0(IQ,LYU)                     &
     &                     + ZPLK1(IQ,LYU) * UTAUPR(LU)
              UAVG(LU)=UAVG(LU) + CWT(IQ-NN) * U0C(IQ,LU)
              FLUP(LU)=FLUP(LU) + CWT(IQ-NN) * CMU(IQ-NN) * U0C(IQ,LU)
           ENDDO

           FLUP(LU) =2.0 * PI * FLUP(LU)
           IF(FLUP(LU).LT.0.)FLUP(LU)=0.
           FLDN(LU) =2.0 * PI * FLDN(LU)
           FDNTOT=FLDN(LU) + FLDIR(LU)
           FNET  =FDNTOT - FLUP(LU)
           RFLDN(LU)=FDNTOT - RFLDIR(LU)
           IF(RFLDN(LU).LT.0.)THEN
               IF(NPR.EQ.-1 .AND. LFLUXS .AND. .NOT.LJMASS)THEN
                   LFLUXS=.FALSE.
                   WRITE(IPR,'(/2A,/22X,2A)')                           &
     &               ' WARNING from FLUXES:  A DISORT diffuse',         &
     &               ' downward flux value was negative and',           &
     &               ' reset to zero.  This warning will not',          &
     &               ' be repeated.'
               ENDIF
               RFLDN(LU)=0.
           ENDIF
           UAVG(LU)=(2*PI*UAVG(LU)+DIRINT)/(4*PI)
           PLSORC= XR0(LYU) + XR1(LYU) * UTAUPR(LU)
           DFDT(LU)=(1.0-SSALB(LYU)) * 4.*PI* (UAVG(LU) - PLSORC)
   90     CONTINUE
          IF(.NOT.LJMASS .AND. PRNT(2))WRITE(*,                         &
     &      '(1P,G13.6,I7,7E12.3,E14.3)')UTAU(LU),LYU,RFLDIR(LU),       &
     &      RFLDN(LU),FDNTOT,FLUP(LU),FNET,UAVG(LU),PLSORC,DFDT(LU)     &
      ENDDO

      IF(.NOT.PRNT(3) .OR. LJMASS)RETURN
      WRITE(*,'(//A)')' ******** AZIMUTHALLY AVERAGED INTENSITIES'      &
     &  //' (AT POLAR QUADRATURE ANGLES) *******'
      DO LU=1, NTAU
          WRITE(*,'(/A,1P,G13.6,//2(5X,A))')'OPTICAL DEPTH =',UTAU(LU), &
     &      'ANGLE (DEG)   COS(ANGLE)     INTENSITY  D(INTENSITY)',     &
     &      'ANGLE (DEG)   COS(ANGLE)     INTENSITY  D(INTENSITY)'
          DO IQ=1, NN
              ANG1=180/PI * ACOS(CMU(2*NN-IQ+1))
              ANG2=180/PI * ACOS(CMU(IQ))
              IF(LU.EQ.1)THEN
                  WRITE(*,'(0P,F16.4,F13.5,1P,E14.3,14X,                &
     &                      0P,F16.4,F13.5,1P,2E14.3)')                 &
     &              ANG1,CMU(2*NN-IQ+1),U0C(IQ,LU),ANG2,CMU(IQ),        &
     &              U0C(IQ+NN,LU),U0C(IQ+NN,LU)-U0C(IQ+NN,LU+1)         &
     &                *EXP((UTAU(LU)-UTAU(LU+1))/CMU(IQ))
              ELSEIF(LU.EQ.NTAU)THEN
                  WRITE(*,'(2(0P,F16.4,F13.5,1P,2E14.3))')              &
     &              ANG1,CMU(2*NN-IQ+1),U0C(IQ,LU),U0C(IQ,LU)           &
     &              -U0C(IQ,LU-1)*EXP((UTAU(LU)-UTAU(LU-1))             &
     &              /CMU(2*NN-IQ+1)),ANG2,CMU(IQ),U0C(IQ+NN,LU)
              ELSE
                  WRITE(*,'(2(0P,F16.4,F13.5,1P,2E14.3))')              &
     &              ANG1,CMU(2*NN-IQ+1),U0C(IQ,LU),U0C(IQ,LU)           &
     &              -U0C(IQ,LU-1)*EXP((UTAU(LU)-UTAU(LU-1))             &
     &              /CMU(2*NN-IQ+1)),ANG2,CMU(IQ),U0C(IQ+NN,LU),        &
     &              U0C(IQ+NN,LU)-U0C(IQ+NN,LU+1)                       &
     &                *EXP((UTAU(LU)-UTAU(LU+1))/CMU(IQ))
              ENDIF
          ENDDO
      ENDDO

      RETURN
      END
