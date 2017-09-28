      SUBROUTINE  SOLEIG(AMB, APB, ARRAY, CMU, CWT, GL, MAZIM,          &
     &                   NN, NSTR, WK, YLMC, CC, EVECC, EVAL, KK, GC)
!         SOLVES EIGENVALUE/VECTOR PROBLEM NECESSARY TO CONSTRUCT
!         HOMOGENEOUS PART OF DISCRETE ORDINATE SOLUTION; STWJ(8B)
!         ** NOTE ** EIGENVALUE PROBLEM IS DEGENERATE WHEN SINGLE
!                    SCATTERING ALBEDO=1;  PRESENT WAY OF DOING IT
!                    SEEMS NUMERICALLY MORE STABLE THAN ALTERNATIVE
!                    METHODS THAT WE TRIED

!     ROUTINES CALLED:  ASYMTX

!   I N P U T     V A R I A B L E S:

!       GL     :  DELTA-M SCALED LEGENDRE COEFFICIENTS OF PHASE FUNCTION
!                    (INCLUDING FACTORS 2L+1 AND SINGLE-SCATTER ALBEDO)
!       CMU    :  COMPUTATIONAL POLAR ANGLE COSINES
!       CWT    :  WEIGHTS FOR QUADRATURE OVER POLAR ANGLE COSINE
!       MAZIM  :  ORDER OF AZIMUTHAL COMPONENT
!       NN     :  HALF THE TOTAL NUMBER OF STREAMS
!       YLMC   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL
!                 AT THE QUADRATURE ANGLES -CMU-
!       (REMAINDER ARE 'DISORT' INPUT VARIABLES)

!   O U T P U T    V A R I A B L E S:

!       CC     :  CAPITAL-C-SUB-IJ IN EQ. SS(5); NEEDED IN SS(15&18)
!       EVAL   :  -NN- EIGENVALUES OF EQ. SS(12) ON RETURN FROM 'ASYMTX'
!                    BUT THEN SQUARE ROOTS TAKEN
!       EVECC  :  -NN- EIGENVECTORS  (G+) - (G-)  ON RETURN
!                    FROM 'ASYMTX' (COLUMN J CORRESPONDS TO -EVAL(J)-)
!                    BUT THEN  (G+) + (G-)  IS CALCULATED FROM SS(10),
!                    G+  AND  G-  ARE SEPARATED, AND  G+  IS STACKED ON
!                    TOP OF  G-  TO FORM -NSTR- EIGENVECTORS OF SS(7)
!       GC     :  PERMANENT STORAGE FOR ALL -NSTR- EIGENVECTORS, BUT
!                    IN AN ORDER CORRESPONDING TO -KK-
!       KK     :  PERMANENT STORAGE FOR ALL -NSTR- EIGENVALUES OF SS(7),
!                    BUT RE-ORDERED WITH NEGATIVE VALUES FIRST (SQUARE
!                    ROOTS OF -EVAL- TAKEN AND NEGATIVES ADDED)

!   I N T E R N A L   V A R I A B L E S:

!       AMB,APB :  MATRICES (ALPHA-BETA), (ALPHA+BETA) IN REDUCED
!                    EIGENVALUE PROBLEM
!       ARRAY   :  COMPLETE COEFFICIENT MATRIX OF REDUCED EIGENVALUE
!                    PROBLEM: (ALFA+BETA)*(ALFA-BETA)
!       GPPLGM  :  (G+) + (G-) (CF. EQS. SS(10-11))
!       GPMIGM  :  (G+) - (G-) (CF. EQS. SS(10-11))
!       WK      :  SCRATCH ARRAY REQUIRED BY 'ASYMTX'
!+---------------------------------------------------------------------+
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     ARGUMENTS:
      INTEGER MAZIM,NN,NSTR
      REAL    AMB(MI,*), APB(MI,*), ARRAY(MI,*), CC(MXCMU,*),           &
     &        CMU(*), CWT(*), EVAL(*), EVECC(MXCMU,*), GC(MXCMU,*),     &
     &        GL(0:*), KK(*), WK(*), YLMC(0:MXCMU,*)

!     LOCAL VARIABLES:
      INTEGER IQ,JQ,L,KQ,IER
      REAL SUMVAL,ALPHA,BETA,GPPLGM,GPMIGM

!                             ** CALCULATE QUANTITIES IN EQS. SS(5-6)
      DO IQ =1, NN

         DO JQ=1, NSTR
            SUMVAL=0.0
            DO L=MAZIM, NSTR-1
               SUMVAL=SUMVAL + GL(L) * YLMC(L,IQ) * YLMC(L,JQ)
            ENDDO
            CC(IQ,JQ)=0.5 * SUMVAL * CWT(JQ)
         ENDDO

         DO JQ=1, NN
!                             ** FILL REMAINDER OF ARRAY USING SYMMETRY
!                             ** RELATIONS  C(-MUI,MUJ)=C(MUI,-MUJ)
!                             ** AND        C(-MUI,-MUJ)=C(MUI,MUJ)

            CC(IQ+NN,JQ)=CC(IQ,JQ+NN)
            CC(IQ+NN,JQ+NN)=CC(IQ,JQ)
!                                      ** GET FACTORS OF COEFF. MATRIX
!                                      ** OF REDUCED EIGENVALUE PROBLEM
            ALPHA=  CC(IQ,JQ) / CMU(IQ)
            BETA=CC(IQ,JQ+NN) / CMU(IQ)
            AMB(IQ,JQ)=ALPHA - BETA
            APB(IQ,JQ)=ALPHA + BETA
         ENDDO
         AMB(IQ,IQ)=AMB(IQ,IQ) - 1.0 / CMU(IQ)
         APB(IQ,IQ)=APB(IQ,IQ) - 1.0 / CMU(IQ)

      ENDDO
!                      ** FINISH CALCULATION OF COEFFICIENT MATRIX OF
!                      ** REDUCED EIGENVALUE PROBLEM:  GET MATRIX
!                      ** PRODUCT (ALFA+BETA)*(ALFA-BETA); SS(12)
      DO IQ=1, NN
         DO JQ=1, NN
            SUMVAL=0.
            DO KQ=1, NN
               SUMVAL=SUMVAL + APB(IQ,KQ) * AMB(KQ,JQ)
            ENDDO
            ARRAY(IQ,JQ)=SUMVAL
         ENDDO
      ENDDO
!                      ** FIND (REAL) EIGENVALUES AND EIGENVECTORS

      CALL  ASYMTX(ARRAY, EVECC, EVAL, NN, MI, MXCMU, IER, WK)
      IF (IER.GT.0)  THEN
         IF(.NOT.LJMASS)WRITE(*, '(//,A,I4,A)')                         &
     &     ' ASYMTX--EIGENVALUE NO. ', IER,                             &
     &     '  DIDNT CONVERGE.  LOWER-NUMBERED EIGENVALUES WRONG.'
         CALL  ERRMSG('ASYMTX--CONVERGENCE PROBLEMS', .TRUE.)
      END IF

!DIR$ IVDEP
      DO IQ=1, NN
         EVAL(IQ)=SQRT(ABS(EVAL(IQ)))
         KK(IQ+NN)=EVAL(IQ)
!                                             ** ADD NEGATIVE EIGENVALUE
         KK(NN+1-IQ)=- EVAL(IQ)
      ENDDO
!                          ** FIND EIGENVECTORS (G+) + (G-) FROM SS(10)
!                          ** AND STORE TEMPORARILY IN -APB- ARRAY
      DO JQ=1, NN
         DO IQ=1, NN
            SUMVAL=0.
            DO KQ=1,NN
               SUMVAL=SUMVAL + AMB(IQ,KQ) * EVECC(KQ,JQ)
            ENDDO
            APB(IQ,JQ)=SUMVAL / EVAL(JQ)
         ENDDO
      ENDDO

      DO JQ=1, NN
!DIR$ IVDEP
         DO IQ=1, NN
            GPPLGM=APB(IQ,JQ)
            GPMIGM=EVECC(IQ,JQ)
!                                ** RECOVER EIGENVECTORS G+,G- FROM
!                                   THEIR SUM AND DIFFERENCE; STACK THEM
!                                   TO GET EIGENVECTORS OF FULL SYSTEM
!                                   SS(7) (JQ=EIGENVECTOR NUMBER)

            EVECC(IQ,      JQ)=0.5 * (GPPLGM + GPMIGM)
            EVECC(IQ+NN,   JQ)=0.5 * (GPPLGM - GPMIGM)

!                                ** EIGENVECTORS CORRESPONDING TO
!                                ** NEGATIVE EIGENVALUES (CORRESP. TO
!                                ** REVERSING SIGN OF 'K' IN SS(10))
            GPPLGM=- GPPLGM
            EVECC(IQ,   JQ+NN)=0.5 * (GPPLGM + GPMIGM)
            EVECC(IQ+NN,JQ+NN)=0.5 * (GPPLGM - GPMIGM)
            GC(IQ+NN,   JQ+NN)  =EVECC(IQ,    JQ)
            GC(NN+1-IQ, JQ+NN)  =EVECC(IQ+NN, JQ)
            GC(IQ+NN,   NN+1-JQ)=EVECC(IQ,    JQ+NN)
            GC(NN+1-IQ, NN+1-JQ)=EVECC(IQ+NN, JQ+NN)
         ENDDO
      ENDDO

      RETURN
      END
