      SUBROUTINE  CMPINT(FBEAM, GC, KK, LAYRU, LL, LYRCUT, MAZIM,       &
     &                    MXCMU, MXUMU, NCUT, NN, NSTR, PLANK, NTAU,    &
     &                    TAUCPR, UMU0, UTAUPR, ZZ, ZPLK0, ZPLK1,       &
     &                    UUM)

!          CALCULATES THE FOURIER INTENSITY COMPONENTS AT THE QUADRATURE
!              ANGLES FOR AZIMUTHAL EXPANSION TERMS (MAZIM) IN EQ. SD(2)

!                  I N P U T    V A R I A B L E S:

!       KK      :  EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)
!       GC      :  EIGENVECTORS AT POLAR QUADRATURE ANGLES IN EQ. SC(1)
!       LL      :  CONSTANTS OF INTEGRATION IN EQ. SC(1), OBTAINED
!                  BY SOLVING SCALED VERSION OF EQ. SC(5);
!                  EXPONENTIAL TERM OF EQ. SC(12) NOT INCLUDED
!       LYRCUT  :  LOGICAL FLAG FOR TRUNCATION OF COMPUTATIONAL LAYER
!       MAZIM   :  ORDER OF AZIMUTHAL COMPONENT
!       NCUT    :  NUMBER OF COMPUTATIONAL LAYER WHERE ABSORPTION
!                  OPTICAL DEPTH EXCEEDS -ABSCUT-
!       NN      :  ORDER OF DOUBLE-GAUSS QUADRATURE (NSTR/2)
!       TAUCPR  :  CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED)
!       UTAUPR  :  OPTICAL DEPTHS OF USER OUTPUT LEVELS IN DELTA-M
!                  COORDINATES;  EQUAL TO -UTAU- IF NO DELTA-M
!       ZZ      :  BEAM SOURCE VECTORS IN EQ. SS(19)
!       ZPLK0   :  THERMAL SOURCE VECTORS -Z0-, BY SOLVING EQ. SS(16)
!       ZPLK1   :  THERMAL SOURCE VECTORS -Z1-, BY SOLVING EQ. SS(16)
!       (REMAINDER ARE 'DisORT' INPUT VARIABLES)

!                  O U T P U T   V A R I A B L E S:

!       UUM     :  FOURIER COMPONENTS OF THE INTENSITY IN EQ. SD(12)
!                  (AT POLAR QUADRATURE ANGLES)

!                  I N T E R N A L   V A R I A B L E S:

!       FACT    :  EXP(- UTAUPR / UMU0)
!       ZINT    :  INTENSITY OF M=0 CASE, IN EQ. SC(1)
!+----------------------------------------------------------------------

         LOGICAL  LYRCUT, PLANK
         INTEGER  LAYRU(*)
         REAL     UUM(MXUMU,*)
         REAL     GC(MXCMU,MXCMU,*), KK(MXCMU,*), LL(MXCMU,*),          &
     &          TAUCPR(0:*), UTAUPR(*), ZZ(MXCMU,*),                    &
     &          ZPLK0(MXCMU,*), ZPLK1(MXCMU,*)

!                                                  LOOP OVER USER LEVELS
         DO 40 LU=1, NTAU

            LYU=LAYRU(LU)
            IF (LYRCUT .AND. LYU.GT.NCUT)  GO TO 40

            DO 30 IQ=1, NSTR
               ZINT=0.0
               DO 10 JQ=1, NN
                 ZINT=ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *               &
     &                    EXP(- KK(JQ,LYU)*(UTAUPR(LU) - TAUCPR(LYU)))
10           CONTINUE
               DO 20 JQ=NN+1, NSTR
                  ZINT=ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *              &
     &                  EXP(- KK(JQ,LYU)*(UTAUPR(LU) - TAUCPR(LYU-1)))
20           CONTINUE

               UUM(IQ,LU)=ZINT
               IF (FBEAM.GT.0.0)                                        &
     &            UUM(IQ,LU)=ZINT +                                     &
     &                         ZZ(IQ,LYU) * EXP(- UTAUPR(LU) / UMU0)
               IF (PLANK .AND. MAZIM.EQ.0)                              &
     &            UUM(IQ,LU)=UUM(IQ,LU) + ZPLK0(IQ,LYU) +               &
     &                         ZPLK1(IQ,LYU) * UTAUPR(LU)
30        CONTINUE

40    CONTINUE

        RETURN
        END
