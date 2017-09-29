      SUBROUTINE  SGEFA( A, LDA, N, IPVT, INFO )

!         FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION.

!         REVISION DATE:  8/1/82
!         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

!     SGEFA IS USUALLY CALLED BY SGECO, BUT IT CAN BE CALLED
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
!     (TIME FOR SGECO) = (1 + 9/N)*(TIME FOR SGEFA) .

!     INPUT:  SAME AS 'SGECO'

!     ON RETURN:

!        A,IPVT  SAME AS 'SGECO'

!        INFO    INTEGER
!                = 0  NORMAL VALUE.
!                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
!                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
!                     INDICATE THAT SGESL OR SGEDI WILL DIVIDE BY ZERO
!                     IF CALLED.  USE  RCOND  IN SGECO FOR A RELIABLE
!                     INDICATION OF SINGULARITY.

!     ROUTINES CALLED:  FROM BLAS:    SAXPY, SSCAL, ISAMAX

      INTEGER  LDA, N, IPVT(*), INFO
      REAL     A(LDA,*)

      REAL     T
      INTEGER  ISAMAX,J,K,KP1,L,NM1

!                      ** GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
      INFO = 0
      NM1 = N - 1
      DO 60 K = 1, NM1
         KP1 = K + 1
!                                            ** FIND L = PIVOT INDEX
         L = ISAMAX( N-K+1, A(K,K), 1) + K-1
         IPVT(K) = L

         IF (A(L,K) .EQ. 0.0E0) THEN

!                                     ** ZERO PIVOT IMPLIES THIS COLUMN
!                                     ** ALREADY TRIANGULARIZED
            INFO = K
         ELSE
!                                     ** INTERCHANGE IF NECESSARY
            IF (L .NE. K) THEN
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
            ENDIF
!                                     ** COMPUTE MULTIPLIERS
            T = -1.0E0 / A(K,K)
            CALL SSCAL( N-K, T, A(K+1,K), 1 )

!                              ** ROW ELIMINATION WITH COLUMN INDEXING
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .NE. K) THEN
                  A(L,J) = A(K,J)
                  A(K,J) = T
               ENDIF
               CALL SAXPY(N-K,T,A(K+1,K),A(K+1,J))
   30       CONTINUE

         ENDIF

   60 CONTINUE

      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0E0) INFO = N
      RETURN
      END
