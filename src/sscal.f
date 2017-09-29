      SUBROUTINE     SSCAL( N, SA, SX, INCX )

!         CALCULATE  X = A*X  (X = VECTOR, A = SCALAR)

!  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR
!            SA  SINGLE PRECISION SCALE FACTOR
!            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
!          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'

! --OUTPUT-- SX  REPLACE  SX(1+I*INCX)  WITH  SA * SX(1+I*INCX)
!                FOR I = 0 TO N-1

      REAL SA, SX(*)

      IF( N.LE.0 ) RETURN

      IF( INCX.NE.1 ) THEN

          DO 10  I = 1, 1+(N-1)*INCX, INCX
             SX(I) = SA * SX(I)
   10     CONTINUE

      ELSE

         M = MOD(N,5)
         IF( M.NE.0 ) THEN
!                           ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
!                           ** IS A MULTIPLE OF 5.
            DO 30  I = 1, M
               SX(I) = SA * SX(I)
   30       CONTINUE
         ENDIF
!                             ** UNROLL LOOP FOR SPEED
         DO 50  I = M+1, N, 5
            SX(I)   = SA * SX(I)
            SX(I+1) = SA * SX(I+1)
            SX(I+2) = SA * SX(I+2)
            SX(I+3) = SA * SX(I+3)
            SX(I+4) = SA * SX(I+4)
   50    CONTINUE

      ENDIF

      RETURN
      END
