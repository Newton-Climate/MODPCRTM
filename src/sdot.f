      REAL FUNCTION  SDOT( N, SX, INCX, SY, INCY )

!          S.P. DOT PRODUCT OF VECTORS  'X'  AND  'Y'

!  --INPUT--
!        N  NUMBER OF ELEMENTS IN INPUT VECTORS 'X' AND 'Y'
!       SX  SING-PREC ARRAY CONTAINING VECTOR 'X'
!     INCX  SPACING OF ELEMENTS OF VECTOR 'X' IN 'SX'
!       SY  SING-PREC ARRAY CONTAINING VECTOR 'Y'
!     INCY  SPACING OF ELEMENTS OF VECTOR 'Y' IN 'SY'

! --OUTPUT--
!     SDOT   SUM FOR I = 0 TO N-1 OF  SX(LX+I*INCX) * SY(LY+I*INCY),
!            WHERE  LX = 1          IF INCX .GE. 0,
!                      = (-INCX)*N  IF INCX .LT. 0,
!            AND LY IS DEFINED IN A SIMILAR WAY USING INCY.

      REAL SX(*), SY(*)

      SDOT = 0.0
      IF( N.LE.0 )  RETURN

      IF ( INCX.EQ.INCY .AND. INCX.GT.1 )  THEN

          DO 10  I = 1, 1+(N-1)*INCX, INCX
             SDOT = SDOT + SX(I) * SY(I)
   10     CONTINUE

      ELSE IF ( INCX.EQ.INCY .AND. INCX.EQ.1 )  THEN

!                                        ** EQUAL, UNIT INCREMENTS
         M = MOD(N,5)
         IF( M .NE. 0 ) THEN
!                            ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
!                            ** IS A MULTIPLE OF 4.
            DO 20  I = 1, M
               SDOT = SDOT + SX(I) * SY(I)
   20       CONTINUE
         ENDIF
!                              ** UNROLL LOOP FOR SPEED
         DO 30  I = M+1, N, 5
            SDOT = SDOT + SX(I)*SY(I)     + SX(I+1)*SY(I+1)             &
     &                  + SX(I+2)*SY(I+2) + SX(I+3)*SY(I+3)             &
     &                  + SX(I+4)*SY(I+4)
   30    CONTINUE

      ELSE
!               ** NONEQUAL OR NONPOSITIVE INCREMENTS.
         IX = 1
         IY = 1
         IF( INCX.LT.0 )  IX = 1 + (N-1)*(-INCX)
         IF( INCY.LT.0 )  IY = 1 + (N-1)*(-INCY)
         DO 40  I = 1, N
            SDOT = SDOT + SX(IX) * SY(IY)
            IX = IX + INCX
            IY = IY + INCY
   40    CONTINUE

      ENDIF

      RETURN
      END
