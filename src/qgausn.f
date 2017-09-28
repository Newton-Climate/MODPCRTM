      SUBROUTINE  QGAUSN(M, GMU, GWT)

!         COMPUTE WEIGHTS AND ABSCISSAE FOR ORDINARY GAUSSIAN QUADRATURE
!             (NO WEIGHT FUNCTION INSIDE INTEGRAL) ON THE INTERVAL (0,1)

!   INPUT :    M                     ORDER OF QUADRATURE RULE

!   OUTPUT :  GMU(I)  I=1 TO M,    ARRAY OF ABSCISSAE
!             GWT(I)  I=1 TO M,    ARRAY OF WEIGHTS

!   REFERENCE:  DAVIS, P.J. AND P. RABINOWITZ, METHODS OF NUMERICAL
!                   INTEGRATION, ACADEMIC PRESS, NEW YORK, PP. 87, 1975.

!   METHOD:  COMPUTE THE ABSCISSAE AS ROOTS OF THE LEGENDRE
!            POLYNOMIAL P-SUB-M USING A CUBICALLY CONVERGENT
!            REFINEMENT OF NEWTON's method.  Compute the
!            weights from EQ. 2.7.3.8 of Davis/Rabinowitz.  Note
!            that Newton'S METHOD CAN VERY EASILY DIVERGE; ONLY A
!            VERY GOOD INITIAL GUESS CAN GUARANTEE CONVERGENCE.
!            THE INITIAL GUESS USED HERE HAS NEVER LED TO DIVERGENCE
!            EVEN FOR M UP TO 1000.

!   ACCURACY:  AT LEAST 13 SIGNIFICANT DIGITS

!   INTERNAL VARIABLES:

!    ITER      : NUMBER OF NEWTON METHOD ITERATIONS
!    MAXIT     : MAXIMUM ALLOWED ITERATIONS OF NEWTON METHOD
!    PM2,PM1,P : 3 SUCCESSIVE LEGENDRE POLYNOMIALS
!    PPR       : DERIVATIVE OF LEGENDRE POLYNOMIAL
!    P2PRI     : 2ND DERIVATIVE OF LEGENDRE POLYNOMIAL
!    TOL       : CONVERGENCE CRITERION FOR LEGENDRE POLY ROOT ITERATION
!    X,XI      : SUCCESSIVE ITERATES IN CUBICALLY-CONVERGENT VERSION
!                OF NEWTONS METHOD (SEEKING ROOTS OF LEGENDRE POLY.)
!+---------------------------------------------------------------------+
      REAL     CONA, GMU(*), GWT(*), T
      INTEGER  ITER, LIM, M, MAXIT, NP1
      REAL     EN, NNP1, ONE, P, PM1, PM2, PPR, P2PRI, PROD,            &
     &         TMP, TOL, TWO, X, XI

!     /RCNSTN/
!       PI       THE CONSTANT PI
!       DEG      NUMBER OF DEGREES IN ONE RADIAN.
!       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.
!       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.
!       RRIGHT   SMALLEST SINGLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
      REAL PI,DEG,BIGNUM,BIGEXP,RRIGHT
      COMMON/RCNSTN/PI,DEG,BIGNUM,BIGEXP,RRIGHT
      SAVE     TOL
      DATA     TOL / 0.0 /,  MAXIT / 1000 /,  ONE / 1.0 /,  TWO / 2.0 /

      IF(TOL.EQ.0.0)TOL=10.*RRIGHT

      IF (M.LT.1)  CALL ERRMSG('QGAUSN--Bad value of M', .TRUE.)
      IF (M.EQ.1)  THEN
         GMU(1)=0.5
         GWT(1)=1.0
         RETURN
      END IF

      EN  =M
      NP1 =M + 1
      NNP1=M * NP1
      CONA=FLOAT(M-1) / (8 * M**3)

      LIM =M / 2
      DO 30  K=1, LIM
!                                INITIAL GUESS FOR K-TH ROOT OF LEGENDRE
!                           POLYNOMIAL, FROM DAVIS/RABINOWITZ (2.7.3.3A)

         T=(4*K - 1) * PI / (4*M + 2)
         X=COS (T + CONA / TAN(T))
         ITER=0
!                             UPWARD RECURRENCE FOR LEGENDRE POLYNOMIALS
10       ITER=ITER + 1
         PM2=ONE
         PM1=X
         DO 20 NN=2, M
            P  =((2*NN - 1) * X * PM1 - (NN-1) * PM2) / NN
            PM2=PM1
            PM1=P
20       CONTINUE
!                                                          NEWTON METHOD
         TMP  =ONE / (ONE - X**2)
         PPR  =EN * (PM2 - X * P) * TMP
         P2PRI=(TWO * X * PPR - NNP1 * P) * TMP
         XI   =X - (P / PPR) * (ONE +                                   &
     &               (P / PPR) * P2PRI / (TWO * PPR))

!                                                  CHECK FOR CONVERGENCE
         IF (ABS(XI-X) .GT. TOL) THEN
            IF(ITER.GT.MAXIT)                                           &
     &          CALL ERRMSG('QGAUSN--MAX ITERATION COUNT', .TRUE.)
            X=XI
            GO TO 10
         END IF

!            ITERATION FINISHED--CALCULATE WEIGHTS, ABSCISSAE FOR (-1,1)

         GMU(K)=- X
         GWT(K)=TWO / (TMP * (EN * PM2)**2)
         GMU(NP1 - K)=- GMU(K)
         GWT(NP1 - K)=  GWT(K)
30    CONTINUE
!                  SET MIDDLE ABSCISSA AND WEIGHT FOR RULES OF ODD ORDER

      IF (MOD(M,2) .NE. 0)  THEN
         GMU(LIM + 1)=0.0
         PROD=ONE
         DO 40 K=3, M, 2
            PROD=PROD * K / (K-1)
40       CONTINUE
         GWT(LIM + 1)=TWO / PROD**2
      END IF
!                                           CONVERT FROM (-1,1) TO (0,1)
      DO 50  K=1, M
         GMU(K)=0.5 * GMU(K) + 0.5
         GWT(K)=0.5 * GWT(K)
50    CONTINUE

      RETURN
      END
