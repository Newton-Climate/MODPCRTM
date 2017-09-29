      SUBROUTINE  SGBFA( ABD, LDA, N, ML, MU, IPVT, INFO )

!         FACTORS A REAL BAND MATRIX BY ELIMINATION.

!         REVISION DATE:  8/1/82
!         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

!     SGBFA IS USUALLY CALLED BY SBGCO, BUT IT CAN BE CALLED
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.

!     INPUT:  SAME AS 'SGBCO'

!     ON RETURN:

!        ABD,IPVT    SAME AS 'SGBCO'

!        INFO    INTEGER
!                = 0  NORMAL VALUE.
!                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
!                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
!                     INDICATE THAT SGBSL WILL DIVIDE BY ZERO IF
!                     CALLED.  USE  RCOND  IN SBGCO FOR A RELIABLE
!                     INDICATION OF SINGULARITY.

!     (SEE 'SGBCO' FOR DESCRIPTION OF BAND STORAGE MODE)

!     ROUTINES CALLED:  FROM BLAS:    SAXPY, SSCAL, ISAMAX
!                       FROM FORTRAN: MAX0, MIN0

      INTEGER  LDA, N, ML, MU, IPVT(*), INFO
      REAL     ABD(LDA,*)

      REAL     T
      INTEGER  I,ISAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1

      M = ML + MU + 1
      INFO = 0
!                        ** ZERO INITIAL FILL-IN COLUMNS
      J0 = MU + 2
      J1 = MIN0(N, M) - 1
      DO 20 JZ = J0, J1
         I0 = M + 1 - JZ
         DO 10 I = I0, ML
            ABD(I,JZ) = 0.0E0
   10    CONTINUE
   20 CONTINUE
      JZ = J1
      JU = 0

!                       ** GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
      NM1 = N - 1
      DO 120 K = 1, NM1
         KP1 = K + 1
!                                  ** ZERO NEXT FILL-IN COLUMN
         JZ = JZ + 1
         IF (JZ .LE. N) THEN
            DO 40 I = 1, ML
               ABD(I,JZ) = 0.0E0
   40       CONTINUE
         ENDIF
!                                  ** FIND L = PIVOT INDEX
         LM = MIN0(ML, N-K)
         L = ISAMAX(LM+1, ABD(M,K), 1) + M - 1
         IPVT(K) = L + K - M

         IF (ABD(L,K) .EQ. 0.0E0) THEN
!                                      ** ZERO PIVOT IMPLIES THIS COLUMN
!                                      ** ALREADY TRIANGULARIZED
            INFO = K
         ELSE
!                                ** INTERCHANGE IF NECESSARY
            IF (L .NE. M) THEN
               T = ABD(L,K)
               ABD(L,K) = ABD(M,K)
               ABD(M,K) = T
            ENDIF
!                                   ** COMPUTE MULTIPLIERS
            T = -1.0E0 / ABD(M,K)
            CALL SSCAL(LM, T, ABD(M+1,K), 1)

!                               ** ROW ELIMINATION WITH COLUMN INDEXING

            JU = MIN0(MAX0(JU, MU+IPVT(K)), N)
            MM = M
            DO 80 J = KP1, JU
               L = L - 1
               MM = MM - 1
               T = ABD(L,J)
               IF (L .NE. MM) THEN
                  ABD(L,J) = ABD(MM,J)
                  ABD(MM,J) = T
               ENDIF
               CALL SAXPY(LM,T,ABD(M+1,K),ABD(MM+1,J))
   80       CONTINUE

         ENDIF

  120 CONTINUE

      IPVT(N) = N
      IF (ABD(M,N) .EQ. 0.0E0) INFO = N
      RETURN
      END
