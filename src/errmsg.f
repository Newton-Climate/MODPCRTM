      SUBROUTINE ERRMSG(MESSAG,LFATAL)

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     COMMONS:
      INCLUDE 'IFIL.h'

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD
!        PRINT OUT A WARNING OR ERROR MESSAGE;  ABORT IF ERROR

      LOGICAL       LFATAL, ONCE
      CHARACTER     MESSAG*(*)
      INTEGER       NUMMSG
      SAVE          NUMMSG, ONCE
      DATA NUMMSG / 0 /,  ONCE / .FALSE. /

      IF ( LFATAL )  THEN
          IF(LJMASS)THEN
              WRITE ( IPR, '(/,2A)' )  ' ******* ERROR >>>>>>  ', MESSAG
              CALL WRTBUF(FATAL)
          ELSE
              WRITE ( *, '(/,2A)' )  ' ******* ERROR >>>>>>  ', MESSAG
          ENDIF
          STOP
      END IF

      NUMMSG = NUMMSG + 1

      IF ( NUMMSG.GT.MAXMSG )  THEN
         IF ( .NOT.ONCE )  THEN
             IF(LJMASS)THEN
                 WRITE (IPR,100)
             ELSE
                 WRITE ( *,100)
             ENDIF
             ONCE = .TRUE.
         ENDIF
      ELSE
         IF(LJMASS)THEN
             WRITE (IPR, '(/,2A)' )  ' ******* WARNING >>>>>>  ', MESSAG
         ELSE
             WRITE ( *, '(/,2A)' )  ' ******* WARNING >>>>>>  ', MESSAG
         ENDIF
      ENDIF

      RETURN

  100 FORMAT( ///,' >>>>>>  TOO MANY WARNING MESSAGES --  ',            &
     &   'THEY WILL NO LONGER BE PRINTED  <<<<<<<', /// )
      END
