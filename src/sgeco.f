      SUBROUTINE  SGECO( A, LDA, N,IPVT, RCOND, Z )

!         FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION
!         AND ESTIMATES THE CONDITION OF THE MATRIX.

!         REVISION DATE:  8/1/82
!         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

!         IF  RCOND  IS NOT NEEDED, SGEFA IS SLIGHTLY FASTER.
!         TO SOLVE  A*X = B , FOLLOW SGECO BY SGESL.

!     ON ENTRY

!        A       REAL(LDA, N)
!                THE MATRIX TO BE FACTORED.

!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .

!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .

!     ON RETURN

!        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
!                WHICH WERE USED TO OBTAIN IT.
!                THE FACTORIZATION CAN BE WRITTEN  A = L*U , WHERE
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

!     ROUTINES CALLED:  FROM LINPACK: SGEFA
!                       FROM BLAS:    SAXPY, SDOT, SSCAL, SASUM
!                       FROM FORTRAN: ABS, AMAX1, SIGN

      INTEGER  LDA, N, IPVT(*)
      REAL     A(LDA,*), Z(*)
      REAL     RCOND

      REAL     SDOT,EK,T,WK,WKM
      REAL     ANORM,S,SASUM,SM,YNORM
      INTEGER  INFO,J,K,KB,KP1,L

!                        ** COMPUTE 1-NORM OF A
      ANORM = 0.0E0
      DO 10 J = 1, N
         ANORM = AMAX1( ANORM, SASUM(N,A(1,J),1) )
   10 CONTINUE
!                                      ** FACTOR
      CALL SGEFA(A,LDA,N,IPVT,INFO)

!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
!     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
!     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
!     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
!     OVERFLOW.

!                        ** SOLVE TRANS(U)*W = E
      EK = 1.0E0
      DO 20 J = 1, N
         Z(J) = 0.0E0
   20 CONTINUE

      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0E0) EK = SIGN(EK, -Z(K))
         IF (ABS(EK-Z(K)) .GT. ABS(A(K,K))) THEN
            S = ABS(A(K,K)) / ABS(EK-Z(K))
            CALL SSCAL(N, S, Z, 1)
            EK = S*EK
         ENDIF
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = ABS(WK)
         SM = ABS(WKM)
         IF (A(K,K) .NE. 0.0E0) THEN
            WK  = WK  / A(K,K)
            WKM = WKM / A(K,K)
         ELSE
            WK  = 1.0E0
            WKM = 1.0E0
         ENDIF
         KP1 = K + 1
         IF (KP1 .LE. N) THEN
            DO 60 J = KP1, N
               SM = SM + ABS(Z(J)+WKM*A(K,J))
               Z(J) = Z(J) + WK*A(K,J)
               S = S + ABS(Z(J))
   60       CONTINUE
            IF (S .LT. SM) THEN
               T = WKM - WK
               WK = WKM
               DO 70 J = KP1, N
                  Z(J) = Z(J) + T*A(K,J)
   70          CONTINUE
            ENDIF
         ENDIF
         Z(K) = WK
  100 CONTINUE

      S = 1.0E0 / SASUM(N, Z, 1)
      CALL SSCAL(N, S, Z, 1)
!                                ** SOLVE TRANS(L)*Y = W
      DO 120 KB = 1, N
         K = N + 1 - KB
         IF (K .LT. N) Z(K) = Z(K) + SDOT(N-K, A(K+1,K), 1, Z(K+1), 1)
         IF (ABS(Z(K)) .GT. 1.0E0) THEN
            S = 1.0E0/ABS(Z(K))
            CALL SSCAL(N, S, Z, 1)
         ENDIF
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE

      S = 1.0E0 / SASUM(N, Z, 1)
      CALL SSCAL(N, S, Z, 1)
!                                 ** SOLVE L*V = Y
      YNORM = 1.0E0
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         IF(K.LT.N)CALL SAXPY(N-K,T,A(K+1,K),Z(K+1))
         IF (ABS(Z(K)) .GT. 1.0E0) THEN
            S = 1.0E0/ABS(Z(K))
            CALL SSCAL(N, S, Z, 1)
            YNORM = S*YNORM
         ENDIF
  140 CONTINUE

      S = 1.0E0 / SASUM(N, Z, 1)
      CALL SSCAL(N, S, Z, 1)
!                                  ** SOLVE  U*Z = V
      YNORM = S*YNORM
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (ABS(Z(K)) .GT. ABS(A(K,K))) THEN
            S = ABS(A(K,K))/ABS(Z(K))
            CALL SSCAL(N, S, Z, 1)
            YNORM = S*YNORM
         ENDIF
         IF (A(K,K) .NE. 0.0E0) Z(K) = Z(K)/A(K,K)
         IF (A(K,K) .EQ. 0.0E0) Z(K) = 1.0E0
         T = -Z(K)
         CALL SAXPY(K-1,T,A(1,K),Z(1))
  160 CONTINUE
!                                   ** MAKE ZNORM = 1.0
      S = 1.0E0 / SASUM(N, Z, 1)
      CALL SSCAL(N, S, Z, 1)
      YNORM = S*YNORM

      IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
      RETURN
      END
