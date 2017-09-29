      INTEGER FUNCTION PSLCT(ANGLE,HRANGE,BETA)

!     THIS ROUTINE RETURNS THE AN INTEGER VALUE INDICATING THE TYPE
!     OF SLANT (ITYPE=2) PATH.  THE FOLLOWING VALUES ARE RETURNED:
!          PSLCT = 1 FOR CASE 2A (H1,H2,ANGLE)
!          PSLCT = 2 FOR CASE 2B (H1,ANGLE,HRANGE)
!          PSLCT = 3 FOR CASE 2C (H1,H2,HRANGE)
!          PSLCT = 4 FOR CASE 2D (H1,H2,BETA)

!     H1     H2     ANGLE     HRANGE     BETA        CASE  PSLCT
!--------------------------------------------     ----------------
!     X      X      X                                 2A     1

!     X             X         X                       2B     2

!     X      X                X                       2C     3

!     X      X                           X            2D     4
      IMPLICIT NONE

!     INPUT ARGUMENTS:
      DOUBLE PRECISION ANGLE,HRANGE,BETA
      IF(BETA.GT.0.D0)THEN

!         BETA > 0 IMPLIES CASE 2D (H1,H2,BETA)
          PSLCT=4
      ELSEIF(HRANGE.GT.0.D0)THEN

!         HRANGE > 0 AND BETA = 0 IMPLIES CASE 2B OR 2C.
          IF(ANGLE.EQ.0.D0)THEN

!             HRANGE > 0, BETA = 0 AND ANGLE = 0 IMPLIES CASE 2C
              PSLCT=3
          ELSE

!             HRANGE > 0, BETA = 0 AND ANGLE NON-ZERO IMPLIES CASE 2B
              PSLCT=2
          ENDIF
      ELSE

!         HRANGE =0 AND BETA = 0 IMPLIES CASE 2A
          PSLCT=1
      ENDIF
      RETURN
      END
