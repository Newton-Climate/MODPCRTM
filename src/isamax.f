      INTEGER FUNCTION ISAMAX(N,SX,INCX)

!  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR OF INTEREST
!            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
!          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'

! --OUTPUT-- ISAMAX   FIRST I, I = 1 TO N, TO MAXIMIZE
!                         ABS(SX(1+(I-1)*INCX))

      REAL SX(*),SMAX,ABS_SX

!     ASSUMES N>0:
      ISAMAX=1
      SMAX=ABS(SX(1))
      I=1
      DO II=2,N
          I=I+INCX
          ABS_SX=ABS(SX(I))
          IF(ABS_SX.GT.SMAX)THEN
              ISAMAX=II
              SMAX=ABS_SX
          ENDIF
      ENDDO
      RETURN
      END
