      SUBROUTINE  SGESL( A, LDA, N,IPVT, B, JOB )

!         SOLVES THE REAL SYSTEM
!            A * X = B  OR  TRANS(A) * X = B
!         USING THE FACTORS COMPUTED BY SGECO OR SGEFA.

!         REVISION DATE:  8/1/82
!         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

!     ON ENTRY

!        A       REAL(LDA, N)
!                THE OUTPUT FROM SGECO OR SGEFA.

!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .

!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .

!        IPVT    INTEGER(N)
!                THE PIVOT VECTOR FROM SGECO OR SGEFA.

!        B       REAL(N)
!                THE RIGHT HAND SIDE VECTOR.

!        JOB     INTEGER
!                = 0         TO SOLVE  A*X = B ,
!                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE
!                            TRANS(A)  IS THE TRANSPOSE.

!     ON RETURN

!        B       THE SOLUTION VECTOR  X .

!     ERROR CONDITION

!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
!        ZERO ON THE DIAGONAL.  TECHNICALLY, THIS INDICATES SINGULARITY,
!        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
!        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
!        CALLED CORRECTLY AND IF SGECO HAS SET RCOND .GT. 0.0
!        OR SGEFA HAS SET INFO .EQ. 0 .

!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
!     WITH  P  COLUMNS
!           CALL SGECO(A,LDA,N,IPVT,RCOND,Z)
!           IF (RCOND IS TOO SMALL) GO TO ...
!           DO 10 J = 1, P
!              CALL SGESL(A,LDA,N,IPVT,C(1,J),0)
!        10 CONTINUE

!     ROUTINES CALLED:  FROM BLAS:    SAXPY, SDOT

      INTEGER  LDA, N, IPVT(*), JOB
      REAL     A(LDA,*), B(*)

      REAL     SDOT,T
      INTEGER  K,KB,L,NM1

      NM1 = N - 1
      IF (JOB .EQ. 0) THEN
!                                 ** JOB = 0 , SOLVE  A * X = B
!                                     ** FIRST SOLVE  L*Y = B
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .NE. K) THEN
               B(L) = B(K)
               B(K) = T
            ENDIF
            CALL SAXPY(N-K,T,A(K+1,K),B(K+1))
   20    CONTINUE
!                                    ** NOW SOLVE  U*X = Y
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K) / A(K,K)
            T = -B(K)
            CALL SAXPY(K-1,T,A(1,K),B(1))
   40    CONTINUE

      ELSE
!                         ** JOB = NONZERO, SOLVE  TRANS(A) * X = B
!                                    ** FIRST SOLVE  TRANS(U)*Y = B
         DO 60 K = 1, N
            T = SDOT( K-1, A(1,K), 1, B(1), 1 )
            B(K) = (B(K) - T) / A(K,K)
   60    CONTINUE
!                                    ** NOW SOLVE  TRANS(L)*X = Y
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + SDOT( N-K, A(K+1,K), 1, B(K+1), 1 )
            L = IPVT(K)
            IF (L .NE. K) THEN
               T = B(L)
               B(L) = B(K)
               B(K) = T
            ENDIF
   80    CONTINUE

      ENDIF

      RETURN
      END
