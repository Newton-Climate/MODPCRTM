      REAL FUNCTION  PLKAVG(WNUMLO, WNUMHI, T)

!        COMPUTES PLANCK FUNCTION INTEGRATED BETWEEN TWO WAVENUMBERS

!  NOTE ** CHANGE 'R1MACH' TO 'D1MACH' TO RUN IN DOUBLE PRECISION

!  I N P U T :  WNUMLO : LOWER WAVENUMBER (INV CM) OF SPECTRAL
!                           INTERVAL
!               WNUMHI : UPPER WAVENUMBER
!               T      : TEMPERATURE (K)

!  O U T P U T :  PLKAVG : INTEGRATED PLANCK FUNCTION (WATTS/SQ M)
!                          =INTEGRAL (WNUMLO TO WNUMHI) OF
!                              2H C**2  NU**3 / (EXP(HC NU/KT) - 1)
!                              (WHERE H=PLANCKS CONSTANT, C=SPEED OF
!                              LIGHT, NU=WAVENUMBER, T=TEMPERATURE,
!                              AND K=BOLTZMANN CONSTANT)

!  REFERENCE : SPECIFICATIONS OF THE PHYSICAL WORLD: NEW VALUE
!                 OF THE FUNDAMENTAL CONSTANTS, DIMENSIONS/N.B.S.,
!                 JAN. 1974

!  METHOD :  FOR  -WNUMLO-  CLOSE TO  -WNUMHI-, A SIMPSON-RULE
!            QUADRATURE IS DONE TO AVOID ILL-CONDITIONING; OTHERWISE

!            (1)  FOR 'WNUMLO' OR 'WNUMHI' SMALL,
!                 INTEGRAL(0 TO WNUMLO/HI) IS CALCULATED BY EXPANDING
!                 THE INTEGRAND IN A POWER SERIES AND INTEGRATING
!                 TERM BY TERM;

!            (2)  OTHERWISE, INTEGRAL(WNUMLO/HI TO INFINITY) IS
!                 CALCULATED BY EXPANDING THE DENOMINATOR OF THE
!                 INTEGRAND IN POWERS OF THE EXPONENTIAL AND
!                 INTEGRATING TERM BY TERM.

!  ACCURACY :  AT LEAST 6 SIGNIFICANT DIGITS, ASSUMING THE
!              PHYSICAL CONSTANTS ARE INFINITELY ACCURATE

!  ERRORS WHICH ARE NOT TRAPPED:

!      * POWER OR EXPONENTIAL SERIES MAY UNDERFLOW, GIVING NO
!        SIGNIFICANT DIGITS.  THIS MAY OR MAY NOT BE OF CONCERN,
!        DEPENDING ON THE APPLICATION.

!      * SIMPSON-RULE SPECIAL CASE IS SKIPPED WHEN DENOMINATOR OF
!        INTEGRAND WILL CAUSE OVERFLOW.  IN THAT CASE THE NORMAL
!        PROCEDURE IS USED, WHICH MAY BE INACCURATE IF THE
!        WAVENUMBER LIMITS (WNUMLO, WNUMHI) ARE CLOSE TOGETHER.
! ----------------------------------------------------------------------
!                                   *** ARGUMENTS
      REAL     T, WNUMLO, WNUMHI

!     /RCNSTN/
!       PI       THE CONSTANT PI
!       DEG      NUMBER OF DEGREES IN ONE RADIAN.
!       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.
!       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.
!       RRIGHT   SMALLEST SINGLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
      REAL PI,DEG,BIGNUM,BIGEXP,RRIGHT
      COMMON/RCNSTN/PI,DEG,BIGNUM,BIGEXP,RRIGHT
!                                   *** LOCAL VARIABLES

!        A1,2,... :  POWER SERIES COEFFICIENTS
!        C2       :  H * C / K, IN UNITS CM*K (H=PLANCK'S CONSTANT,
!                      C=SPEED OF LIGHT, K=BOLTZMANN CONSTANT)
!        D(I)     :  EXPONENTIAL SERIES EXPANSION OF INTEGRAL OF
!                       PLANCK FUNCTION FROM WNUMLO (I=1) OR WNUMHI
!                       (I=2) TO INFINITY
!        EPSIL    :  SMALLEST NUMBER SUCH THAT 1+EPSIL .GT. 1 ON
!                       COMPUTER
!        EX       :  EXP(- V(I))
!        EXM      :  EX**M
!        MMAX     :  NO. OF TERMS TO TAKE IN EXPONENTIAL SERIES
!        MV       :  MULTIPLES OF 'V(I)'
!        P(I)     :  POWER SERIES EXPANSION OF INTEGRAL OF
!                       PLANCK FUNCTION FROM ZERO TO WNUMLO (I=1) OR
!                       WNUMHI (I=2)
!        SIGMA    :  STEFAN-BOLTZMANN CONSTANT (W/M**2/K**4)
!        SIGDPI   :  SIGMA / PI
!        SMALLV   :  NUMBER OF TIMES THE POWER SERIES IS USED (0,1,2)
!        V(I)     :  C2 * (WNUMLO(I=1) OR WNUMHI(I=2)) / TEMPERATURE
!        VCUT     :  POWER-SERIES CUTOFF POINT
!        VCP      :  EXPONENTIAL SERIES CUTOFF POINTS
!        VMAX     :  LARGEST ALLOWABLE ARGUMENT OF 'EXP' FUNCTION

      REAL C2,SIGMA,VCUT
      PARAMETER(C2=1.438786,SIGMA=5.67032E-8,VCUT=1.5,A1=1./3,          &
     &  A2=-1./8,A3=1./60,A4=-1./5040,A5=1./272160,A6=-1./13305600)
      INTEGER  SMALLV
      REAL     CONC,D(2),EPSIL,EX,MV,P(2),SIGDPI,V(2),VCP(7),VSQ
      SAVE     CONC, VMAX, EPSIL, SIGDPI
      DATA     EPSIL/ 0. / , VCP/ 10.25, 5.7, 3.9, 2.9, 2.3, 1.9, 0.0 /
      F(X)=X**3 / (EXP(X) - 1)

      IF (EPSIL.EQ.0.0)  THEN
         VMAX=BIGEXP
         EPSIL=RRIGHT
         SIGDPI=SIGMA / PI
         CONC=15. / PI**4
      END IF

      IF(T.LT.0.0 .OR. WNUMHI.LE.WNUMLO .OR. WNUMLO.LT.0.)              &
     &    CALL ERRMSG('PLKAVG--TEMPERATURE OR WAVENUMS. WRONG', .TRUE.)

      IF (T.LT.1.E-4)  THEN
         PLKAVG=0.0
         RETURN
      ENDIF

      V(1)=C2 * WNUMLO / T
      V(2)=C2 * WNUMHI / T
      IF (V(1).GT.EPSIL .AND. V(2).LT.VMAX .AND.                        &
     &     (WNUMHI-WNUMLO)/WNUMHI .LT. 1.E-2)  THEN

!                          ** WAVENUMBERS ARE VERY CLOSE.  GET INTEGRAL
!                          ** BY ITERATING SIMPSON RULE TO CONVERGENCE.
         HH=V(2) - V(1)
         OLDVAL=0.0
         VAL0=F(V(1)) + F(V(2))

         DO  2  N=1, 10
            DEL=HH / (2*N)
            VAL=VAL0
            DO  1  K=1, 2*N-1
               VAL=VAL + 2*(1+MOD(K,2)) * F(V(1) + K*DEL)
    1       CONTINUE
            VAL=DEL/3. * VAL
            IF (ABS((VAL-OLDVAL)/VAL) .LE. 1.E-6)  GO TO 3
            OLDVAL=VAL
    2    CONTINUE
         CALL ERRMSG('PLKAVG--SIMPSON RULE DIDNT CONVERGE', .FALSE.)

    3    PLKAVG=SIGDPI * T**4 * CONC * VAL
         RETURN
      END IF

      SMALLV=0
      DO  50  I=1, 2

         IF(V(I).LT.VCUT)  THEN
!                                   ** USE POWER SERIES
            SMALLV=SMALLV + 1
            VSQ=V(I)**2
            P(I)= CONC * VSQ * V(I) * (A1 + V(I) * (A2 + V(I) *         &
     &                (A3 + VSQ * (A4 + VSQ * (A5 + VSQ*A6)))))
         ELSE
!                    ** USE EXPONENTIAL SERIES
            MMAX=0
!                                ** FIND UPPER LIMIT OF SERIES
   20       MMAX=MMAX + 1
               IF (V(I).LT.VCP(MMAX))  GO TO 20

            EX=EXP(- V(I))
            EXM=1.0
            D(I)=0.0

            DO  30  M=1, MMAX
               MV=M * V(I)
               EXM=EX * EXM
               D(I)=D(I) +                                              &
     &                EXM * (6. + MV*(6. + MV*(3. + MV))) / M**4
   30       CONTINUE

            D(I)=CONC * D(I)
         END IF

   50 CONTINUE

      IF (SMALLV .EQ. 2) THEN
!                                    ** WNUMLO AND WNUMHI BOTH SMALL
         PLKAVG=P(2) - P(1)

      ELSE IF (SMALLV .EQ. 1) THEN
!                                    ** WNUMLO SMALL, WNUMHI LARGE
         PLKAVG=1. - P(1) - D(2)

      ELSE
!                                    ** WNUMLO AND WNUMHI BOTH LARGE
         PLKAVG=D(1) - D(2)

      END IF

      PLKAVG=SIGDPI * T**4 * PLKAVG
      IF(PLKAVG.EQ.0.0)                                                 &
     &    CALL ERRMSG('PLKAVG--RETURNS ZERO; POSSIBLE UNDERFLOW',       &
     &                 .FALSE.)

      RETURN
      END
