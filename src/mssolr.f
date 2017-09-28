      SUBROUTINE MSSOLR(CMU,CWT,FBEAM,GC,GU,KK,LAYRU,LL,LYRCUT,MXCMU,   &
     &                  MLOS,MXUMU,NCUT,NN,NSTR,NTAU,NUMU,PI,TAUCPR,    &
     &                  UMU0,UTAU,UTAUPR,BEAMMS,ZZ,FDNSRT,S0CMS)

      REAL FACT,FBEAM,FDNSRT,FLDIR,PI,RFLDIR,UMU0,ZINT
      INTEGER IQ,IU,JQ,LU,LYU,MXCMU,MXUMU,NCUT,NN,NSTR,NTAU,NUMU
!    I N P U T     V A R I A B L E S:

!       CMU      :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE
!       CWT      :  WEIGHTS FOR GAUSS QUADRATURE OVER ANGLE COSINE
!       GC       :  EIGENVECTORS AT POLAR QUADRATURE ANGLES, SC(1)
!       GU       :  EIGENVECTORS INTERPOLATED TO USER POLAR ANGLES
!                   (I.E., -G- IN EQ. SC(1))
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
!       UTAU     :  OPTICAL DEPTHS OF USER OUTPUT LEVELS
!       UTAUPR   :  OPTICAL DEPTHS OF USER OUTPUT LEVELS IN DELTA-M
!                     COORDINATES;  EQUAL TO  -UTAU- IF NO DELTA-M
!       ZZ       :  BEAM SOURCE VECTORS IN EQ. SS(19)

!   I N T E R N A L       V A R I A B L E S:

!       FLDIR    :  DIRECT-BEAM FLUX (DELTA-M SCALED)
!       FACT     :  EXP( - UTAUPR / UMU0 )
!       ZINT     :  INTENSITY OF M = 0 CASE, IN EQ. SC(1)
!      RFLDIR    :  DIRECT-BEAM FLUX (NOT DELTA-M SCALED)

!   O U T P U T    V A R I A B L E S:

!       S0CMS    :  MULTIPLE SCATTERING SOLAR SOURCE FUNCTION
!      FDNSRT    :  DOWNWARD DIFFUSE SOLAR FLUX AT SURFACE

      LOGICAL LYRCUT
      REAL S0CMS(MLOS,*)
      INTEGER LAYRU(*)
      REAL CMU(*),CWT(*),GC(MXCMU,MXCMU,*),GU(0:MXUMU,MXCMU,*),         &
     &       KK(MXCMU,*),LL(MXCMU,*),TAUCPR(0:*),UTAU(*),UTAUPR(*),     &
     &       BEAMMS(MXUMU,*),ZZ(MXCMU,*)

!                                               ** LOOP OVER USER LEVELS

      DO 140 LU = 1,NTAU

         LYU = LAYRU(LU)
         IF(LYRCUT .AND. LYU.GT.NCUT)THEN

!                                     ** NO RADIATION REACHES THIS LEVEL
            DO IU = 1,NUMU
               S0CMS(IU,LU)=0.
            ENDDO
         ELSE

!                                        ** RADIATION REACHES THIS LEVEL
            FACT = EXP(-UTAUPR(LU)/UMU0)

!                                               ** LOOP OVER USER ANGLES
            DO 130 IU = 1,NUMU
               ZINT = 0.0
               DO 110 JQ = 1,NN
                  ZINT = ZINT+GU(IU,JQ,LYU)*LL(JQ,LYU)                  &
     &                   *EXP(-KK(JQ,LYU)*(UTAUPR(LU)-TAUCPR(LYU)))
  110          CONTINUE
               DO 120 JQ = NN+1,NSTR
                  ZINT = ZINT+GU(IU,JQ,LYU)*LL(JQ,LYU)                  &
     &                   *EXP(-KK(JQ,LYU)*(UTAUPR(LU)-TAUCPR(LYU-1)))
  120          CONTINUE

!                   **  MS SOURCE FUNCTIONS CALCULATED STW(30) M.S. TERM

               S0CMS(IU,LU) = ZINT+BEAMMS(IU,LYU)*FACT
  130       CONTINUE
         ENDIF
  140 CONTINUE

!*************VINCENT ROSS REMOVED FOR SOURCE FUNCTION IMPROVEMENT*****************
#ifndef DISORT_BOUND_SRC
!                                  **  LAYER AVERAGE OF S0CMS AS MODTRAN
      DO 160 IU = 1,NUMU
         DO 150 LU = 1,NTAU-1
            S0CMS(IU,LU) = (S0CMS(IU,LU)+S0CMS(IU,LU+1))/2.
  150    CONTINUE
  160 CONTINUE
#endif  
!*************END VINCENT ROSS ADDITON*********************************************

      LU = NTAU
      LYU = LAYRU(LU)
      FDNSRT = 0.

      IF(.NOT.(LYRCUT.AND.LYU.GT.NCUT))THEN
!                                           ** RADIATION REACHES SURFACE
         FACT = EXP(-UTAUPR(LU)/UMU0)
         FLDIR = UMU0*(FBEAM*FACT)
         RFLDIR = UMU0*FBEAM*EXP(-UTAU(LU)/UMU0)

         DO 190 IQ = 1,NN
            ZINT = 0.0
            DO 170 JQ = 1,NN
               ZINT = ZINT+GC(IQ,JQ,LYU)*LL(JQ,LYU)                     &
     &                *EXP(-KK(JQ,LYU)*(UTAUPR(LU)-TAUCPR(LYU)))
  170       CONTINUE
            DO 180 JQ = NN+1,NSTR
               ZINT = ZINT+GC(IQ,JQ,LYU)*LL(JQ,LYU)                     &
     &                *EXP(-KK(JQ,LYU)*(UTAUPR(LU)-TAUCPR(LYU-1)))
  180       CONTINUE
            FDNSRT = FDNSRT+CWT(NN+1-IQ)*CMU(NN+1-IQ)                   &
     &               *(ZINT+ZZ(IQ,LYU)*FACT)
  190    CONTINUE

         FDNSRT = 2.0*PI*FDNSRT+FLDIR-RFLDIR
         IF(FDNSRT.LT.0.)FDNSRT = 0.
      ENDIF
      RETURN
      END
