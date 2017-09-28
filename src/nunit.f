      INTEGER FUNCTION NUNIT()

!     THIS INTEGER FUNCTION DETERMINES AN UNUSED
!     FILE UNIT NUMBER LESS THAN OR EQUAL TO 99.

!     COMMONS:

!     COMMON/IFIL/:  FILE UNIT NUMBERS.
      INCLUDE 'IFIL.h'

!     LOCAL VARIABLES
      LOGICAL LOPEN
      INTEGER MUNIT

!     FIND UNUSED UNIT NUMBER AND RETURN (1 THROUGH
!     35 ARE ASSIGNED IN THE devcbd.f ROUTINE).
      DO MUNIT=99,36,-1
          INQUIRE(UNIT=MUNIT,OPENED=LOPEN)
          IF(.NOT.LOPEN)THEN
              NUNIT=MUNIT
              RETURN
          ENDIF
      ENDDO

!     NO FILE UNIT NUMBER AVAILABLE
      WRITE(IPR ,'(/A)')                                                &
     &  ' Error in routine NUNIT:  No file unit numbers available?'
      STOP ' Error in routine NUNIT:  No file unit numbers available?'
      END
