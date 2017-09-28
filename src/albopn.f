      SUBROUTINE ALBOPN(LNFLRT,FLRT,LABEL)

!     OPENS ATM CORRECTION (SPHERICAL ALBEDO & DIFFUSE TRANSM) FILE.
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     COMMONS:
      INCLUDE 'IFIL.h'

!     ARGUMENTS:
!       LNFLRT   LENGTH OF FILE ROOT NAME, 0 IF NO MOD5ROOT FILE OPEN.
!       FLRT     FILE ROOT NAME.
      INTEGER LNFLRT
      CHARACTER FLRT*(NAMLEN-4) , LABEL*6

!     LOCAL VARIABLES:
!       LOPEN    LOGICAL FLAG, TRUE IF FILE IS OPEN.
      LOGICAL LOPEN

      IF(BINOUT)THEN

!         ATM COR DATA IS WRITTEN TO THE .tp8 BINARY OUTPUT FILES:
          INQUIRE(JPR1,OPENED=LOPEN)
          IF(.NOT.LOPEN)THEN
              IF(LNFLRT.GT.0)THEN
                  CALL OPNFL(JPR1,0,FLRT(1:LNFLRT)//'b.tp8',            &
     &              'UNKNOWN','UNFORMATTED','ALBOPN')
              ELSE
                  CALL OPNFL(JPR1,0,'tape8b',                           &
     &              'UNKNOWN','UNFORMATTED','ALBOPN')
              ENDIF
          ENDIF
      ELSE
          INQUIRE(IALB,OPENED=LOPEN)
          IF(.NOT.LOPEN)THEN
              IF(LNFLRT.GT.0)THEN

!                 OPEN *.acd ASCII FILE FOR ATMOSPHERIC CORRECTION DATA:
                  CALL OPNFL(IALB,0,FLRT(1:LNFLRT)//LABEL//'.acd',      &
     &              'UNKNOWN','FORMATTED','ALBOPN')
              ELSE

!                 OPEN ASCII atmcor.dat FOR ATMOSPHERIC CORRECTION DATA:
                  CALL OPNFL(IALB,0,'atmcor.dat',                       &
     &              'UNKNOWN','FORMATTED','ALBOPN')
              ENDIF
          ENDIF

!         WRITE HEADER EVEN IF THE ASCII FILE WAS ALREADY OPENED:
          WRITE(IALB,'(/(A))')                                          &
     &      ' SPECTRAL LOS  k           k         SUN->GND'//           &
     &      '   SUN-GND-OBS      OBS->GND       OBS-GND     SPHERICAL', &
     &      'FREQUENCY  #  INT       WEIGHT        DIFFUSE'//           &
     &      '        DIRECT      EMBEDDED        DIRECT        ALBEDO', &
     &      '   [CM-1]                              TRANSM'//           &
     &      '        TRANSM    DIF TRANSM        TRANSM      FROM GND', &
     &      '--------- --- ---     --------    -----------'//           &
     &      '   -----------   -----------   -----------   -----------'
      ENDIF

!     RETURN TO GEODRV:
      RETURN
      END
