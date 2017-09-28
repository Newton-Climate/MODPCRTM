      SUBROUTINE  SGBCO( ABD, LDA, N, ML, MU, IPVT, RCOND, Z )

!         FACTORS A REAL BAND MATRIX BY GAUSSIAN ELIMINATION
!         AND ESTIMATES THE CONDITION OF THE MATRIX.

!         REVISION DATE:  8/1/82
!         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

!     IF  RCOND  IS NOT NEEDED, SGBFA IS SLIGHTLY FASTER.
!     TO SOLVE  A*X = B , FOLLOW SBGCO BY SGBSL.

!     INPUT:

!        ABD     REAL(LDA, N)
!                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS
!                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND
!                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
!                ML+1 THROUGH 2*ML+MU+1 OF  ABD .
!                SEE THE COMMENTS BELOW FOR DETAILS.

!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  ABD .
!                LDA MUST BE .GE. 2*ML + MU + 1 .

!        N       INTEGER
!                THE ORDER OF THE ORIGINAL MATRIX.

!        ML      INTEGER
!                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
!                0 .LE. ML .LT. N .

!        MU      INTEGER
!                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
!                0 .LE. MU .LT. N .
!                MORE EFFICIENT IF  ML .LE. MU .

!     ON RETURN

!        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
!                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
!                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
!                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
!                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.

!        IPVT    INTEGER(N)
!                AN INTEGER VECTOR OF PIVOT INDICES.

!        RCOND   REAL
!                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
!                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
!                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
!                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
!                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
!                           1.0 + RCOND .EQ. 1.0
!                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
!                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
!                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
!                UNDERFLOWS.

!        Z       REAL(N)
!                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
!                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
!                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .

!     BAND STORAGE

!           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT
!           WILL SET UP THE INPUT.

!                   ML = (BAND WIDTH BELOW THE DIAGONAL)
!                   MU = (BAND WIDTH ABOVE THE DIAGONAL)
!                   M = ML + MU + 1
!                   DO 20 J = 1, N
!                      I1 = MAX0(1, J-MU)
!                      I2 = MIN0(N, J+ML)
!                      DO 10 I = I1, I2
!                         K = I - J + M
!                         ABD(K,J) = A(I,J)
!                10    CONTINUE
!                20 CONTINUE

!           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD .
!           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR
!           ELEMENTS GENERATED DURING THE TRIANGULARIZATION.
!           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 .
!           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE
!           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED.

!     EXAMPLE:  IF THE ORIGINAL MATRIX IS

!           11 12 13  0  0  0
!           21 22 23 24  0  0
!            0 32 33 34 35  0
!            0  0 43 44 45 46
!            0  0  0 54 55 56
!            0  0  0  0 65 66

!      THEN  N = 6, ML = 1, MU = 2, LDA .GE. 5  AND ABD SHOULD CONTAIN

!            *  *  *  +  +  +  , * = NOT USED
!            *  * 13 24 35 46  , + = USED FOR PIVOTING
!            * 12 23 34 45 56
!           11 22 33 44 55 66
!           21 32 43 54 65  *

!     ROUTINES CALLED:  FROM LINPACK: SGBFA
!                       FROM BLAS:    SAXPY, SDOT, SSCAL, SASUM
!                       FROM FORTRAN: ABS, AMAX1, MAX0, MIN0, SIGN

      INTEGER  LDA, N, ML, MU, IPVT(*)
      REAL     ABD(LDA,*), Z(*)
      REAL     RCOND

      REAL     SDOT, EK, T, WK, WKM
      REAL     ANORM, S, SASUM, SM, YNORM
      INTEGER  IS, INFO, J, JU, K, KB, KP1, L, LA, LM, LZ, M, MM

!                       ** COMPUTE 1-NORM OF A
      ANORM = 0.0E0
      L = ML + 1
      IS = L + MU
      DO 10 J = 1, N
         ANORM = AMAX1(ANORM, SASUM(L,ABD(IS,J), 1))
         IF (IS .GT. ML + 1) IS = IS - 1
         IF (J .LE. MU) L = L + 1
         IF (J .GE. N - ML) L = L - 1
   10 CONTINUE
!                                               ** FACTOR
      CALL SGBFA(ABD, LDA, N, ML, MU, IPVT, INFO)

!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
!     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
!     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
!     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
!     OVERFLOW.

!                     ** SOLVE TRANS(U)*W = E
      EK = 1.0E0
      DO 20 J = 1, N
         Z(J) = 0.0E0
   20 CONTINUE

      M = ML + MU + 1
      JU = 0
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0E0) EK = SIGN(EK, -Z(K))
         IF (ABS(EK-Z(K)) .GT. ABS(ABD(M,K))) THEN
            S = ABS(ABD(M,K))/ABS(EK-Z(K))
            CALL SSCAL(N, S, Z, 1)
            EK = S*EK
         ENDIF
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = ABS(WK)
         SM = ABS(WKM)
         IF (ABD(M,K) .NE. 0.0E0) THEN
            WK  = WK /ABD(M,K)
            WKM = WKM/ABD(M,K)
         ELSE
            WK  = 1.0E0
            WKM = 1.0E0
         ENDIF
         KP1 = K + 1
         JU = MIN0(MAX0(JU, MU+IPVT(K)), N)
         MM = M
         IF (KP1 .LE. JU) THEN
            DO 60 J = KP1, JU
               MM = MM - 1
               SM = SM + ABS(Z(J)+WKM*ABD(MM,J))
               Z(J) = Z(J) + WK*ABD(MM,J)
               S = S + ABS(Z(J))
   60       CONTINUE
            IF (S .LT. SM) THEN
               T = WKM - WK
               WK = WKM
               MM = M
               DO 70 J = KP1, JU
                  MM = MM - 1
                  Z(J) = Z(J) + T*ABD(MM,J)
   70          CONTINUE
            ENDIF
         ENDIF
         Z(K) = WK
  100 CONTINUE

      S = 1.0E0 / SASUM(N, Z, 1)
      CALL SSCAL(N, S, Z, 1)

!                         ** SOLVE TRANS(L)*Y = W
      DO 120 KB = 1, N
         K = N + 1 - KB
         LM = MIN0(ML, N-K)
         IF (K .LT. N) Z(K) = Z(K) + SDOT(LM, ABD(M+1,K), 1, Z(K+1), 1)
         IF (ABS(Z(K)) .GT. 1.0E0) THEN
            S = 1.0E0 / ABS(Z(K))
            CALL SSCAL(N, S, Z, 1)
         ENDIF
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE

      S = 1.0E0 / SASUM(N, Z, 1)
      CALL SSCAL(N, S, Z, 1)

      YNORM = 1.0E0
!                         ** SOLVE L*V = Y
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         LM = MIN0(ML, N-K)
         IF(K.LT.N)CALL SAXPY(LM,T,ABD(M+1,K),Z(K+1))
         IF (ABS(Z(K)) .GT. 1.0E0) THEN
            S = 1.0E0 / ABS(Z(K))
            CALL SSCAL(N, S, Z, 1)
            YNORM = S*YNORM
         ENDIF
  140 CONTINUE

      S = 1.0E0/SASUM(N, Z, 1)
      CALL SSCAL(N, S, Z, 1)
      YNORM = S*YNORM
!                           ** SOLVE  U*Z = W
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (ABS(Z(K)) .GT. ABS(ABD(M,K))) THEN
            S = ABS(ABD(M,K)) / ABS(Z(K))
            CALL SSCAL(N, S, Z, 1)
            YNORM = S*YNORM
         ENDIF
         IF (ABD(M,K) .NE. 0.0E0) Z(K) = Z(K)/ABD(M,K)
         IF (ABD(M,K) .EQ. 0.0E0) Z(K) = 1.0E0
         LM = MIN0(K, M) - 1
         LA = M - LM
         LZ = K - LM
         T = -Z(K)
         CALL SAXPY(LM,T,ABD(LA,K),Z(LZ))
  160 CONTINUE
!                              ** MAKE ZNORM = 1.0
      S = 1.0E0 / SASUM(N, Z, 1)
      CALL SSCAL(N, S, Z, 1)
      YNORM = S*YNORM

      IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
      RETURN
      END
