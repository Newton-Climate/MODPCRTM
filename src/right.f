      SUBROUTINE RIGHT(RRIGHT,DRIGHT)

!     RIGHT RETURNS THE SMALLEST SINGLE AND DOUBLE PRECISION REALS THAT
!     CAN BE ADDED TO ONE AND PRODUCE A NUMBER GREATER THAN ONE.

!     ARGUMENTS:
!       RRIGHT   SMALLEST SINGLE PRECISION REAL THAT CAN BE ADDED TO 1.
!       DRIGHT   SMALLEST DOUBLE PRECISION REAL THAT CAN BE ADDED TO 1.
      REAL RRIGHT
      DOUBLE PRECISION DRIGHT

!     LOCAL:
      REAL RRTLOW,R1PLUS,RRTMID
      DOUBLE PRECISION DRTLOW,D1PLUS,DRTMID

!     DATA:
      REAL R1,RRTSAV
      DOUBLE PRECISION D1,DRTSAV
      DATA R1,RRTSAV/2*1.E0/,D1,DRTSAV/2*1.D0/

!     DETERMINE LOWER AND UPPER BOUNDS ON RRIGHT:
      RRTLOW=R1
   10 CONTINUE
      RRIGHT=RRTLOW
      RRTLOW=RRTLOW/2
      R1PLUS=R1+RRTLOW
      IF(R1PLUS.GT.R1)GOTO 10

!     USE BISECTION METHOD TO FIND RRIGHT:
   20 CONTINUE
      RRTMID=(RRTLOW+RRIGHT)/2
      IF(RRTSAV.NE.RRTMID .AND. RRIGHT.GT.1.E-12)THEN
          R1PLUS=R1+RRTMID
          IF(R1PLUS.GT.R1)THEN
              RRIGHT=RRTMID
          ELSE
              RRTLOW=RRTMID
          ENDIF
          RRTSAV=RRTMID
          GOTO 20
      ENDIF

!     FOR SOME MACHINES, THE ABOVE PROCEDURE PRODUCES AN "RRIGHT" VALUE
!     NEARER TO 1.E-16 THAN TO 1.E-08.  PROTECT AGAINST THIS PROBLEM:
      IF(RRIGHT.LT.1.E-12)RRIGHT=6.E-08

!     DETERMINE LOWER AND UPPER BOUNDS ON DRIGHT:
      DRTLOW=DBLE(RRTLOW)
   30 CONTINUE
      DRIGHT=DRTLOW
      DRTLOW=DRTLOW/2
      D1PLUS=D1+DRTLOW
      IF(D1PLUS.GT.D1)GOTO 30

!     USE BISECTION METHOD TO FIND DRIGHT
   40 CONTINUE
      DRTMID=(DRTLOW+DRIGHT)/2
      IF(DRTSAV.NE.DRTMID)THEN
          D1PLUS=D1+DRTMID
          IF(D1PLUS.GT.D1)THEN
              DRIGHT=DRTMID
          ELSE
              DRTLOW=DRTMID
          ENDIF
          DRTSAV=DRTMID
          GOTO 40
      ENDIF
      RETURN
      END
