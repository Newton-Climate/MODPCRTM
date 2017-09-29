      REAL FUNCTION  SASUM( N, SX, INCX )

!  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR TO BE SUMMED
!            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
!          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'

! --OUTPUT-- SASUM   SUM FROM 0 TO N-1 OF  ABS(SX(1+I*INCX))

      REAL SX(*)

      SASUM = 0.0
      IF( N.LE.0 )  RETURN
      IF( INCX.NE.1 ) THEN
!                                          ** NON-UNIT INCREMENTS
          DO 10 I = 1, 1+(N-1)*INCX, INCX
             SASUM = SASUM + ABS(SX(I))
   10     CONTINUE
      ELSE
!                                          ** UNIT INCREMENTS
         M = MOD(N,6)
         IF( M.NE.0 ) THEN
!                             ** CLEAN-UP LOOP SO REMAINING VECTOR
!                             ** LENGTH IS A MULTIPLE OF 6.
            DO 30  I = 1, M
              SASUM = SASUM + ABS(SX(I))
   30       CONTINUE
         ENDIF
!                              ** UNROLL LOOP FOR SPEED
         DO 50  I = M+1, N, 6
           SASUM = SASUM + ABS(SX(I))   + ABS(SX(I+1)) + ABS(SX(I+2))   &
     &                   + ABS(SX(I+3)) + ABS(SX(I+4)) + ABS(SX(I+5))
   50    CONTINUE
      ENDIF

      RETURN
      END
