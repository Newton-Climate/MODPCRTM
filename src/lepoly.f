      SUBROUTINE  LEPOLY(NMU, M, MAXMU, TWONM1, MU, YLM)

!       COMPUTES THE NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL,
!       DEFINED IN TERMS OF THE ASSOCIATED LEGENDRE POLYNOMIAL
!       PLM=P-SUB-L-SUPER-M AS

!             YLM(MU)=SQRT((L-M)!/(L+M)!) * PLM(MU)

!       FOR FIXED ORDER -M- AND ALL DEGREES FROM L=M TO TWONM1.
!       WHEN M.GT.0, ASSUMES THAT Y-SUB(M-1)-SUPER(M-1) IS AVAILABLE
!       FROM A PRIOR CALL TO THE ROUTINE.

!       REFERENCE: Dave, J.V. and B.H. Armstrong, Computations of
!                  High-Order Associated Legendre Polynomials,
!                  J. Quant. Spectrosc. Radiat. Transfer 10,
!                  557-562, 1970.  (hereafter D/A)

!       METHOD: Varying degree recurrence relationship.

!       NOTE 1: The D/A formulas are transformed by
!               setting  M=n-1; L=k-1.
!       NOTE 2: Assumes that routine is called first with  M=0,
!               then with  M=1, etc. up to  M=TWONM1.
!       NOTE 3: Loops are written in such a way as to vectorize.

!  I N P U T     V A R I A B L E S:

!       NMU    :  NUMBER OF ARGUMENTS OF -YLM-
!       M      :  ORDER OF -YLM-
!       MAXMU  :  FIRST DIMENSION OF -YLM-
!       TWONM1 :  MAX DEGREE OF -YLM-
!       MU(I)  :  I=1 TO NMU, ARGUMENTS OF -YLM-
!       IF M.GT.0, YLM(M-1,I) FOR I=1 TO NMU IS REQUIRED

!  O U T P U T     V A R I A B L E:

!       YLM(L,I) :  L=M TO TWONM1, NORMALIZED ASSOCIATED LEGENDRE
!                   POLYNOMIALS EVALUATED AT ARGUMENT -MU(I)-
!+---------------------------------------------------------------------+
      REAL     MU(*), YLM(0:MAXMU,*)
      INTEGER  M, NMU, TWONM1
      PARAMETER  (MAXSQT=1000)
      REAL     SQT(MAXSQT)
      LOGICAL  PASS1
      SAVE  SQT, PASS1
      DATA  PASS1 / .TRUE. /

      IF (PASS1)  THEN
         PASS1=.FALSE.
         DO NS=1, MAXSQT
            SQT(NS)=SQRT(FLOAT(NS))
         ENDDO
      ENDIF

      IF (2*TWONM1 .GT. MAXSQT)                                         &
     &   CALL ERRMSG('LEPOLY--NEED TO INCREASE PARAM MAXSQT', .TRUE.)

      IF (M .EQ. 0)  THEN
!                             ** UPWARD RECURRENCE FOR ORDINARY
!                             ** LEGENDRE POLYNOMIALS
         DO  I=1, NMU
            YLM(0,I)=1.
            YLM(1,I)=MU(I)
         ENDDO

         DO  L=2, TWONM1
            DO  I=1, NMU
               YLM(L,I)=((2*L-1) * MU(I) * YLM(L-1,I)                   &
     &                      - (L-1) * YLM(L-2,I)) / L
            ENDDO
         ENDDO

      ELSE

         DO  I=1, NMU
!                               ** Y-SUB-M-SUPER-M; DERIVED FROM
!                               ** D/A EQS. (11,12)

            YLM(M,I)=- SQT(2*M-1) / SQT(2*M)                            &
     &                  * SQRT(1. - MU(I)**2) * YLM(M-1,I)

!                              ** Y-SUB-(M+1)-SUPER-M; DERIVED FROM
!                              ** D/A EQS. (13,14) USING EQS. (11,12)

            YLM(M+1,I)=SQT(2*M+1) * MU(I) * YLM(M,I)
         ENDDO
!                                   ** UPWARD RECURRENCE; D/A EQ. (10)
         DO  L=M+2, TWONM1
            TMP1=SQT(L-M) * SQT(L+M)
            TMP2=SQT(L-M-1) * SQT(L+M-1)
            DO  I=1, NMU
               YLM(L,I)=((2*L-1) * MU(I) * YLM(L-1,I)                   &
     &                        - TMP2 * YLM(L-2,I)) / TMP1
            ENDDO
         ENDDO

      END IF

      RETURN
      END
