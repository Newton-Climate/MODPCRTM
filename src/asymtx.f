      SUBROUTINE  ASYMTX(AA, EVEC, EVAL, M, IA, IEVEC, IER, WK)

!    =======  S I N G L E    P R E C I S I O N    V E R S I O N  ======

!       SOLVES EIGENFUNCTION PROBLEM FOR REAL ASYMMETRIC MATRIX
!       FOR WHICH IT IS KNOWN A PRIORI THAT THE EIGENVALUES ARE REAL.

!       THIS IS AN ADAPTATION OF A SUBROUTINE EIGRF IN THE IMSL
!       LIBRARY TO USE REAL INSTEAD OF COMPLEX ARITHMETIC, ACCOUNTING
!       FOR THE KNOWN FACT THAT THE EIGENVALUES AND EIGENVECTORS IN
!       THE DISCRETE ORDINATE SOLUTION ARE REAL.  OTHER CHANGES INCLUDE
!       PUTTING ALL THE CALLED SUBROUTINES IN-LINE, DELETING THE
!       PERFORMANCE INDEX CALCULATION, UPDATING MANY DO-LOOPS
!       TO FORTRAN77, AND IN CALCULATING THE MACHINE PRECISION
!       DRIGHT INSTEAD OF SPECIFYING IT IN A DATA STATEMENT.

!       EIGRF IS BASED PRIMARILY ON EISPACK ROUTINES.  THE MATRIX IS
!       FIRST BALANCED USING THE PARLETT-REINSCH ALGORITHM.  THEN
!       THE MARTIN-WILKINSON ALGORITHM IS APPLIED.

!       REFERENCES:
!          DONGARRA, J. AND C. MOLER, EISPACK -- A PACKAGE FOR SOLVING
!             MATRIX EIGENVALUE PROBLEMS, IN COWELL, ED., 1984:
!             SOURCES AND DEVELOPMENT OF MATHEMATICAL SOFTWARE,
!             PRENTICE-HALL, ENGLEWOOD CLIFFS, NJ
!         PARLETT AND REINSCH, 1969: BALANCING A MATRIX FOR CALCULATION
!             OF EIGENVALUES AND EIGENVECTORS, NUM. MATH. 13, 293-304
!         WILKINSON, J., 1965: THE ALGEBRAIC EIGENVALUE PROBLEM,
!             CLARENDON PRESS, OXFORD

!   I N P U T    V A R I A B L E S:

!        A    :  INPUT ASYMMETRIC MATRIX, DESTROYED AFTER SOLVED
!        M    :  ORDER OF  A
!       IA    :  FIRST DIMENSION OF  A
!    IEVEC    :  FIRST DIMENSION OF  EVEC

!   O U T P U T    V A R I A B L E S:

!       EVEC  :  (UNNORMALIZED) EIGENVECTORS OF  A
!                   ( COLUMN J CORRESPONDS TO EVAL(J) )

!       EVAL  :  (UNORDERED) EIGENVALUES OF  A ( DIMENSION AT LEAST M )

!       IER   :  IF .NE. 0, SIGNALS THAT EVAL(IER) FAILED TO CONVERGE;
!                   IN THAT CASE EIGENVALUES IER+1,IER+2,...,M  ARE
!                   CORRECT BUT EIGENVALUES 1,...,IER ARE SET TO ZERO.

!   S C R A T C H   V A R I A B L E S:

!       AAD   :  DOUBLE PRECISION STAND-IN FOR -A-
!       EVECD :  DOUBLE PRECISION STAND-IN FOR -EVEC-
!       EVALD :  DOUBLE PRECISION STAND-IN FOR -EVAL-
!       WKD   :  DOUBLE PRECISION STAND-IN FOR -WK-
!+---------------------------------------------------------------------+

      REAL              AA(IA,*),  WK(*),  EVAL(*), EVEC(IEVEC,*)
      LOGICAL           NOCONV, NOTLAS

!     /RCNSTN/
!       PI       THE CONSTANT PI.
!       DEG      NUMBER OF DEGREES IN ONE RADIAN.
!       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.
!       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.
!       RRIGHT   SMALLEST SINGLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
      REAL PI,DEG,BIGNUM,BIGEXP,RRIGHT
      COMMON/RCNSTN/PI,DEG,BIGNUM,BIGEXP,RRIGHT
      DATA     C1/ 0.4375 /, C2/ 0.5 /, C3/ 0.75 /, C4/ 0.95 /,         &
     &         C5/ 16.0 /, C6/ 256.0 /, ZERO / 0.0 /, ONE / 1.0 /

      IER=0
      IF (M.LT.1 .OR. IA.LT.M .OR. IEVEC.LT.M)                          &
     &     CALL ERRMSG('ASYMTX--BAD INPUT VARIABLE(S)', .TRUE.)

!                           ** HANDLE 1X1 AND 2X2 SPECIAL CASES
      IF (M.EQ.1)  THEN
         EVAL(1)=AA(1,1)
         EVEC(1,1)=ONE
         RETURN

      ELSE IF (M.EQ.2)  THEN
         DISCRI=(AA(1,1) - AA(2,2))**2 + 4. * AA(1,2) * AA(2,1)
         IF (DISCRI.LT.ZERO)                                            &
     &        CALL ERRMSG('ASYMTX--COMPLEX EVALS IN 2X2 CASE', .TRUE.)
         SGN=ONE
         IF (AA(1,1).LT.AA(2,2))  SGN=- ONE
         EVAL(1)=0.5 * (AA(1,1) + AA(2,2) + SGN*SQRT(DISCRI))
         EVAL(2)=0.5 * (AA(1,1) + AA(2,2) - SGN*SQRT(DISCRI))
         EVEC(1,1)=ONE
         EVEC(2,2)=ONE
         IF (AA(1,1).EQ.AA(2,2) .AND.                                   &
     &        (AA(2,1).EQ.ZERO.OR.AA(1,2).EQ.ZERO))  THEN
            RNORM=  ABS(AA(1,1)) + ABS(AA(1,2)) + ABS(AA(2,1))          &
     &              + ABS(AA(2,2))
            W=RRIGHT * RNORM
            EVEC(2,1)=  AA(2,1) / W
            EVEC(1,2)=- AA(1,2) / W
         ELSE
            EVEC(2,1)=AA(2,1) / (EVAL(1) - AA(2,2))
            EVEC(1,2)=AA(1,2) / (EVAL(2) - AA(1,1))
         ENDIF

         RETURN

      END IF
!                                        ** INITIALIZE OUTPUT VARIABLES
      IER=0
      DO I=1, M
         EVAL(I)=ZERO
         DO J=1, M
            EVEC(I,J)=ZERO
         ENDDO
         EVEC(I,I)=ONE
      ENDDO
!                  ** BALANCE THE INPUT MATRIX AND REDUCE ITS NORM BY
!                  ** DIAGONAL SIMILARITY TRANSFORMATION STORED IN WK;
!                  ** THEN SEARCH FOR ROWS ISOLATING AN EIGENVALUE
!                  ** AND PUSH THEM DOWN
      RNORM=ZERO
      L =1
      K =M

30    KKK=K
         DO J=KKK, 1, -1
            ROW=ZERO
            DO I=1, K
               IF (I.NE.J) ROW=ROW + ABS(AA(J,I))
            ENDDO
            IF (ROW.EQ.ZERO) THEN
               WK(K)=J
               IF (J.NE.K) THEN
                  DO I=1, K
                     REPL  =AA(I,J)
                     AA(I,J)=AA(I,K)
                     AA(I,K)=REPL
                  ENDDO
                  DO I=L, M
                     REPL  =AA(J,I)
                     AA(J,I)=AA(K,I)
                     AA(K,I)=REPL
                  ENDDO
               END IF
               K=K - 1
               GO TO 30
            END IF
         ENDDO
!                                     ** SEARCH FOR COLUMNS ISOLATING AN
!                                       ** EIGENVALUE AND PUSH THEM LEFT
80    LLL=L
         DO J=LLL, K
            COL=ZERO
            DO I=L, K
               IF (I.NE.J) COL=COL + ABS(AA(I,J))
            ENDDO
            IF (COL.EQ.ZERO) THEN
               WK(L)=J
               IF (J.NE.L) THEN
                  DO I=1, K
                     REPL  =AA(I,J)
                     AA(I,J)=AA(I,L)
                     AA(I,L)=REPL
                  ENDDO
                  DO I=L, M
                     REPL  =AA(J,I)
                     AA(J,I)=AA(L,I)
                     AA(L,I)=REPL
                  ENDDO
               END IF
               L=L + 1
               GO TO 80
            END IF
         ENDDO
!                           ** BALANCE THE SUBMATRIX IN ROWS L THROUGH K
      DO I=L, K
         WK(I)=ONE
      ENDDO

140   NOCONV=.FALSE.
         DO I=L, K
            COL=ZERO
            ROW=ZERO
            DO J=L, K
               IF (J.NE.I) THEN
                  COL=COL + ABS(AA(J,I))
                  ROW=ROW + ABS(AA(I,J))
               END IF
            ENDDO
            F=ONE
            G=ROW / C5
            H=COL + ROW
160         IF (COL.LT.G) THEN
               F  =F * C5
               COL=COL * C6
               GO TO 160
            END IF
            G=ROW * C5
!lex170     IF (COL.GE.G) THEN
170         IF (COL.GT.G) THEN
               F  =F / C5
               COL=COL / C6
               GO TO 170
            END IF
!                                                         ** NOW BALANCE
            IF ((COL+ROW) / F .LT. C4 * H) THEN
               WK(I) =WK(I) * F
               NOCONV=.TRUE.
               DO J=L, M
                  AA(I,J)=AA(I,J) / F
               ENDDO
               DO J=1, K
                  AA(J,I)=AA(J,I) * F
               ENDDO
            END IF
         ENDDO

      IF (NOCONV) GO TO 140
!                                  ** IS -A- ALREADY IN HESSENBERG FORM?
      IF (K-1.LT.L+1) GO TO 350
!                                   ** TRANSFER -A- TO A HESSENBERG FORM
      DO N=L+1, K-1
         H      =ZERO
         WK(N+M)=ZERO
         SCAL  =ZERO
!                                                        ** SCALE COLUMN
         DO I=N, K
            SCAL=SCAL + ABS(AA(I,N-1))
         ENDDO
         IF (SCAL.NE.ZERO) THEN
            DO I=K, N, -1
               WK(I+M)=AA(I,N-1) / SCAL
               H=H + WK(I+M) * WK(I+M)
            ENDDO
            G=- SIGN(SQRT(H),WK(N+M))
            H=H - WK(N+M) * G
            WK(N+M)=WK(N+M) - G
!                                                 ** FORM (I-(U*UT)/H)*A
            DO J=N, M
               F=ZERO
               DO I=K, N, -1
                  F=F + WK(I+M) * AA(I,J)
               ENDDO
               DO I=N, K
                  AA(I,J)=AA(I,J) - WK(I+M) * F / H
               ENDDO
            ENDDO
!                                    ** FORM (I-(U*UT)/H)*A*(I-(U*UT)/H)
            DO I=1, K
               F=ZERO
               DO J=K, N, -1
                  F=F + WK(J+M) * AA(I,J)
               ENDDO
               DO J=N, K
                  AA(I,J)=AA(I,J) - WK(J+M) * F / H
               ENDDO
            ENDDO
            WK(N+M) =SCAL * WK(N+M)
            AA(N,N-1)=SCAL * G
         END IF
      ENDDO

      DO N=K-2, L, -1
         N1=N + 1
         N2=N + 2
         F =AA(N+1,N)
         IF (F.NE.ZERO) THEN
!MOD
!MOD        PROTECT AGAINST F*WK(N+1+M) BEING 0.
!MOD        F =F * WK(N+1+M)
            F1 = WK(N+1+M)
            DO I=N+2, K
               WK(I+M)=AA(I,N)
            ENDDO
            IF (N+1.LE.K) THEN
               DO J=1, M
                  G=ZERO
                  DO I=N+1, K
                     G=G + WK(I+M) * EVEC(I,J)
                  ENDDO
                  G=G / F
!MOD
!MOD              INCORPORATE F1 INTO THE DENOMINATOR NOW.
                  G=G / F1
                  DO I=N+1, K
                     EVEC(I,J)=EVEC(I,J) + G * WK(I+M)
                  ENDDO
               ENDDO
            END IF
         END IF
      ENDDO

350   CONTINUE
      N=1
      DO I=1, M
         DO J=N, M
            RNORM=RNORM + ABS(AA(I,J))
         ENDDO
         N=I
         IF (I.LT.L .OR. I.GT.K) EVAL(I)=AA(I,I)
      ENDDO
      N=K
      T=ZERO
!                                         ** SEARCH FOR NEXT EIGENVALUES
380   IF (N.LT.L) GO TO 530
      IN=0
      N1=N - 1
      N2=N - 2
!                          ** LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
390   CONTINUE
      DO I=L, N
         LB=N+L - I
         IF (LB.EQ.L) GO TO 410
         S=ABS(AA(LB-1,LB-1)) + ABS(AA(LB,LB))
         IF (S.EQ.ZERO) S=RNORM
         IF (ABS(AA(LB,LB-1)) .LE. RRIGHT * S) GO TO 410
      ENDDO
!                                                      ** ONE EVAL FOUND
410   X=AA(N,N)
      IF (LB.EQ.N ) THEN
         AA(N,N) =X + T
         EVAL(N)=AA(N,N)
         N=N1
         GO TO 380
      END IF
!                                                     ** TWO EVALS FOUND
      Y=AA(N1,N1)
      W=AA(N,N1) * AA(N1,N)
      IF (LB.EQ.N1) THEN
         P=(Y-X) * C2
         Q=P * P + W
         Z=SQRT(ABS(Q))
         AA(N,N)=X + T
         X=AA(N,N)
         AA(N1,N1)=Y + T
!                                                           ** REAL PAIR
         Z=P + SIGN(Z,P)
         EVAL(N1)=X + Z
         EVAL(N) =EVAL(N1)
         IF (Z.NE.ZERO) EVAL(N)=X - W / Z
         X=AA(N,N1)
!                                  ** EMPLOY SCALE FACTOR IN CASE
!                                  ** X AND Z ARE VERY SMALL
         R=SQRT(X * X + Z * Z)
         P=X / R
         Q=Z / R
!                                                    ** ROW MODIFICATION
         DO J=N1, M
            Z=AA(N1,J)
            AA(N1,J)=Q * Z + P * AA(N,J)
            AA(N,J) =Q * AA(N,J) - P * Z
         ENDDO
!                                                 ** COLUMN MODIFICATION
         DO I=1, N
            Z=AA(I,N1)
            AA(I,N1)=Q * Z + P * AA(I,N)
            AA(I,N) =Q * AA(I,N) - P * Z
         ENDDO
!                                          ** ACCUMULATE TRANSFORMATIONS
         DO I=L, K
            Z=EVEC(I,N1)
            EVEC(I,N1)=Q * Z + P * EVEC(I,N)
            EVEC(I,N) =Q * EVEC(I,N) - P * Z
         ENDDO
         N=N2
         GO TO 380
      END IF
!                    ** NO CONVERGENCE AFTER 30 ITERATIONS; SET ERROR
!                    ** INDICATOR TO THE INDEX OF THE CURRENT EIGENVALUE

      IF (IN.EQ.30) THEN
         IER=128 + N
         RETURN
      END IF
!                                                          ** FORM SHIFT
      IF (IN.EQ.10 .OR. IN.EQ.20) THEN
         T=T + X
         DO I=L, N
            AA(I,I)=AA(I,I) - X
         ENDDO
         S=ABS(AA(N,N1)) + ABS(AA(N1,N2))
         X=C3 * S
         Y=X
         W=-C1 * S * S
      END IF

      IN=IN + 1
!                ** LOOK FOR TWO CONSECUTIVE SMALL SUB-DIAGONAL ELEMENTS

      DO J=LB, N2
         I=N2+LB - J
         Z=AA(I,I)
         R=X - Z
         S=Y - Z
         P=(R * S-W) / AA(I+1,I) + AA(I,I+1)
         Q=AA(I+1,I+1) - Z - R - S
         R=AA(I+2,I+1)
         S=ABS(P) + ABS(Q) + ABS(R)
         P=P / S
         Q=Q / S
         R=R / S
         IF (I.EQ.LB) GO TO 470
         UU=ABS(AA(I,I-1)) * (ABS(Q) + ABS(R))
         VV=ABS(P) * (ABS(AA(I-1,I-1)) + ABS(Z) + ABS(AA(I+1,I+1)))
         IF (UU .LE. RRIGHT*VV) GO TO 470
      ENDDO

470   CONTINUE
      AA(I+2,I)=ZERO
      DO J=I+3, N
         AA(J,J-2)=ZERO
         AA(J,J-3)=ZERO
      ENDDO

!             ** DOUBLE QR STEP INVOLVING ROWS K TO N AND COLUMNS M TO N

      DO KA=I, N1
         NOTLAS=KA.NE.N1
         IF (KA.EQ.I) THEN
            S=SIGN(SQRT(P*P + Q*Q + R*R),P)
            IF (LB.NE.I) AA(KA,KA-1)=- AA(KA,KA-1)
         ELSE
            P=AA(KA,KA-1)
            Q=AA(KA+1,KA-1)
            R=ZERO
            IF (NOTLAS) R=AA(KA+2,KA-1)
            X=ABS(P) + ABS(Q) + ABS(R)
            IF (X.EQ.ZERO) GO TO 520
            P=P / X
            Q=Q / X
            R=R / X
            S=SIGN(SQRT(P*P + Q*Q + R*R),P)
            AA(KA,KA-1)=- S * X
         END IF
         P=P + S
         X=P / S
         Y=Q / S
         Z=R / S
         Q=Q / P
         R=R / P
!                                                    ** ROW MODIFICATION
         DO J=KA, M
            P=AA(KA,J) + Q * AA(KA+1,J)
            IF (NOTLAS) THEN
               P=P + R * AA(KA+2,J)
               AA(KA+2,J)=AA(KA+2,J) - P * Z
            END IF
            AA(KA+1,J)=AA(KA+1,J) - P * Y
            AA(KA,J)  =AA(KA,J)   - P * X
         ENDDO
!                                                 ** COLUMN MODIFICATION
         IILIM=MIN0(N,KA+3)
         DO II=1, IILIM
            P=X * AA(II,KA) + Y * AA(II,KA+1)
            IF (NOTLAS) THEN
               P=P + Z * AA(II,KA+2)
               AA(II,KA+2)=AA(II,KA+2) - P * R
            END IF
            AA(II,KA+1)=AA(II,KA+1) - P * Q
            AA(II,KA)  =AA(II,KA) - P
         ENDDO
!                                          ** ACCUMULATE TRANSFORMATIONS
         DO II=L, K
            P=X * EVEC(II,KA) + Y * EVEC(II,KA+1)
            IF (NOTLAS) THEN
               P=P + Z * EVEC(II,KA+2)
               EVEC(II,KA+2)=EVEC(II,KA+2) - P * R
            END IF
            EVEC(II,KA+1)=EVEC(II,KA+1) - P * Q
            EVEC(II,KA)  =EVEC(II,KA) - P
         ENDDO
520   CONTINUE
      ENDDO
      GO TO 390
!                     ** ALL EVALS FOUND, NOW BACKSUBSTITUTE REAL VECTOR
530   CONTINUE
      IF (RNORM.NE.ZERO) THEN
         DO N=M, 1, -1
            N2=N
            AA(N,N)=ONE
            DO I=N-1, 1, -1
               W=AA(I,I) - EVAL(N)
               IF (W.EQ.ZERO) W=RRIGHT * RNORM
               R=AA(I,N)
               DO J=N2, N-1
                  R=R + AA(I,J) * AA(J,N)
               ENDDO
               AA(I,N)=-R / W
               N2=I
            ENDDO
         ENDDO
!                      ** END BACKSUBSTITUTION VECTORS OF ISOLATED EVALS

         DO I=1, M
            IF (I.LT.L .OR. I.GT.K) THEN
               DO J=I, M
                  EVEC(I,J)=AA(I,J)
               ENDDO
            END IF
         ENDDO
!                                   ** MULTIPLY BY TRANSFORMATION MATRIX
         IF (K.NE.0) THEN
            DO J=M, L, -1
               NLIM=MIN0(J,K)
               DO I=L, K
                  Z=ZERO
                  DO N=L, NLIM
                     Z=Z + EVEC(I,N) * AA(N,J)
                  ENDDO
                  EVEC(I,J)=Z
               ENDDO
            ENDDO
         END IF

      END IF

      DO I=L, K
         DO J=1, M
            EVEC(I,J)=EVEC(I,J) * WK(I)
         ENDDO
      ENDDO
!                           ** INTERCHANGE ROWS IF PERMUTATIONS OCCURRED
      DO I=L-1, 1, -1
         J=INT(WK(I))
         IF (I.NE.J) THEN
            DO N=1, M
               REPL     =EVEC(I,N)
               EVEC(I,N)=EVEC(J,N)
               EVEC(J,N)=REPL
            ENDDO
         END IF
      ENDDO

      DO I=K+1, M
         J=INT(WK(I))
         IF (I.NE.J) THEN
            DO N=1, M
               REPL     =EVEC(I,N)
               EVEC(I,N)=EVEC(J,N)
               EVEC(J,N)=REPL
            ENDDO
         END IF
      ENDDO

      RETURN
      END
