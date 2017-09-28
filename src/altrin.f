      SUBROUTINE  ALTRIN(GU, KK, LL, MXCMU, MXUMU, NLYR,                &
     &                    NN, NSTR, NUMU, TAUCPR, UMU, U0U, WK)

!       COMPUTES AZIMUTHALLY-AVERAGED INTENSITY AT TOP AND BOTTOM
!       OF MEDIUM (RELATED TO ALBEDO AND TRANSMISSION OF MEDIUM BY
!       RECIPROCITY PRINCIPLES;  SEE REF S2).  USER POLAR ANGLES ARE
!       USED AS INCIDENT BEAM ANGLES. (THIS IS A VERY SPECIALIZED
!       VERSION OF 'USRINT')

!       ** NOTE **  USER INPUT VALUES OF -UMU- (ASSUMED POSITIVE) ARE
!                   TEMPORARILY IN UPPER LOCATIONS OF  -UMU-  AND
!                   CORRESPONDING NEGATIVES ARE IN LOWER LOCATIONS
!                   (THIS MAKES -GU- COME OUT RIGHT).  I.E. THE CONTENTS
!                   OF THE TEMPORARY -UMU- ARRAY ARE:

!                     -UMU(NUMU),..., -UMU(1), UMU(1),..., UMU(NUMU)

!   I N P U T    V A R I A B L E S:

!       GU     :  EIGENVECTORS INTERPOLATED TO USER POLAR ANGLES
!                 (I.E., G IN EQ. SC(1))
!       KK     :  EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)
!       LL     :  CONSTANTS OF INTEGRATION IN EQ. SC(1), OBTAINED
!                 BY SOLVING SCALED VERSION OF EQ. SC(5);
!                 EXPONENTIAL TERM OF EQ. SC(12) NOT INCLUDED
!       NN     :  ORDER OF DOUBLE-GAUSS QUADRATURE (NSTR/2)
!       TAUCPR :  CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED)
!       (REMAINDER ARE 'DISORT' INPUT VARIABLES)

!   O U T P U T    V A R I A B L E:

!       U0U  :    DIFFUSE AZIMUTHALLY-AVERAGED INTENSITY AT TOP AND
!                 BOTTOM OF MEDIUM (DIRECTLY TRANSMITTED COMPONENT,
!                 CORRESPONDING TO -BNDINT- IN 'USRINT', IS OMITTED).

!   I N T E R N A L    V A R I A B L E S:

!       DTAU   :  OPTICAL DEPTH OF A COMPUTATIONAL LAYER
!       PALINT :  NON-BOUNDARY-FORCED INTENSITY COMPONENT
!       UTAUPR :  OPTICAL DEPTHS OF USER OUTPUT LEVELS (DELTA-M SCALED)
!       WK     :  SCRATCH VECTOR FOR SAVING 'EXP' EVALUATIONS
!       ALL THE EXPONENTIAL FACTORS (I.E., EXP1, EXPN,... ETC.)
!       COME FROM THE SUBSTITUTION OF CONSTANTS OF INTEGRATION IN
!       EQ. SC(12) INTO EQS. S1(8-9).  ALL HAVE NEGATIVE ARGUMENTS.
!+---------------------------------------------------------------------+

      REAL     UTAUPR(2)
      REAL GU(0:MXUMU,MXCMU,*),KK(MXCMU,*),LL(MXCMU,*),MU,              &
     &         TAUCPR(0:*), UMU(*), U0U(MXUMU,*), WK(*)

      UTAUPR(1)=0.0
      UTAUPR(2)=TAUCPR(NLYR)
      DO 100  LU=1, 2
         IF (LU.EQ.1)  THEN
            IUMIN=NUMU / 2 + 1
            IUMAX=NUMU
            SGN=1.0
         ELSE
            IUMIN=1
            IUMAX=NUMU / 2
            SGN=- 1.0
         END IF
!                                   ** LOOP OVER POLAR ANGLES AT WHICH
!                                   ** ALBEDOS/TRANSMISSIVITIES DESIRED
!                                   ** (UPWARD ANGLES AT TOP BOUNDARY,
!                                   ** DOWNWARD ANGLES AT BOTTOM)
         DO 50  IU=IUMIN, IUMAX
            MU=UMU(IU)
!                                     ** INTEGRATE FROM TOP TO BOTTOM
!                                     ** COMPUTATIONAL LAYER
            PALINT=0.0
            DO 30  LC=1, NLYR

               DTAU=TAUCPR(LC) - TAUCPR(LC-1)
               EXP1= EXP((UTAUPR(LU) - TAUCPR(LC-1)) / MU)
               EXP2= EXP((UTAUPR(LU) - TAUCPR(LC)) / MU)

!                                      ** -KK- IS NEGATIVE
               DO 20  IQ=1, NN
                  WK(IQ)=EXP(KK(IQ,LC) * DTAU)
                  DENOM=1.0 + MU * KK(IQ,LC)
                  IF (ABS(DENOM).LT.0.0001) THEN
!                                                   ** L'HOSPITAL LIMIT
                     EXPN=DTAU / MU * EXP2
                  ELSE
                     EXPN=(EXP1 * WK(IQ) - EXP2) * SGN / DENOM
                  END IF
                  PALINT=PALINT + GU(IU,IQ,LC) * LL(IQ,LC) * EXPN
20             CONTINUE
!                                      ** -KK- IS POSITIVE
               DO 21  IQ=NN+1, NSTR
                  DENOM=1.0 + MU * KK(IQ,LC)
                  IF (ABS(DENOM).LT.0.0001) THEN
                     EXPN=- DTAU / MU * EXP1
                  ELSE
                     EXPN=(EXP1 - EXP2 * WK(NSTR+1-IQ)) *SGN / DENOM
                  END IF
                  PALINT=PALINT + GU(IU,IQ,LC) * LL(IQ,LC) * EXPN
21             CONTINUE

30          CONTINUE

            U0U(IU, LU)=PALINT

 50      CONTINUE
100   CONTINUE

      RETURN
      END
