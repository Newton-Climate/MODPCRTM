      SUBROUTINE OPNFL(IUNIT,LENREC,FNAME,STAT,FORMT,ROUTIN)

!     OPNFL CHECKS & OPENS A FILE, PROVIDING LOCAL CONTROL FOR FAILURES.
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     INPUT ARGUMENTS:
!       IUNIT    FILE UNIT NUMBER.
!       LENREC   RECORD LENGTH (>0 FOR DIRECT ACCESS FILES).
!       FNAME    NAME OF FILE THAT IS TO BE OPENED.
!       STAT     FILE STATUS ('OLD', 'NEW', 'UNKNOWN' OR 'SCRATCH').
!       FORMT    FILE FORMAT ('FORMATTED' OR 'UNFORMATTED').
!       ROUTIN   CALLING ROUTINE NAME.
      INTEGER IUNIT,LENREC
      CHARACTER FNAME*(*),STAT*(*),FORMT*(*),ROUTIN*(*)

!     COMMONS:
      INCLUDE 'IFIL.h'

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     LOCAL VARIABLES:
!       LEXIST   LOGICAL FLAG, TRUE IF FILE EXISTS.
!       LOPEN    LOGICAL FLAG, TRUE IF FILE IS OPEN.
!       LENNAM   LENGTH OF FILE NAME.
!       NUM      FILE UNIT NUMBER OF CONNECTED FILE.
!       FILNAM   FILE NAME.
      LOGICAL LEXIST,LOPEN
      INTEGER LENNAM,NUM
      CHARACTER FILNAM*(NAMLEN)

!     FUNCTIONS:
!       LENSTR   RETURNS STRING LENGTH AFTER TRIMMING LEADING BLANKS.
      INTEGER LENSTR

!     CHECK FOR PROPER FILE STATUS STRING:
      CALL UPCASE(STAT)
      IF(STAT.NE.'OLD' .AND. STAT.NE.'NEW' .AND.                        &
     &  STAT.NE.'UNKNOWN' .AND. STAT.NE.'SCRATCH')THEN

!         FAILURE:  File STATUS has an unallowed value.
          INQUIRE(UNIT=IPR,OPENED=LOPEN)
          IF(LOPEN)WRITE(IPR,'(/4A)')'OPNFL error from ',ROUTIN,        &
     &      ': File STATUS =',STAT
          IF(LJMASS)THEN
              CALL WRTBUF(FATAL)
          ELSE
              WRITE(*,'(4A)')'OPNFL error from ',ROUTIN,                &
     &          ': File STATUS =',STAT
          ENDIF
          STOP 'OPNFL ERROR: File STATUS string has an unallowed value.'
      ENDIF

!     CHECK FOR PROPER FILE FORMAT STRING:
      CALL UPCASE(FORMT)
      IF(FORMT.NE.'FORMATTED' .AND. FORMT.NE.'UNFORMATTED')THEN

!         FAILURE:  File FORMAT string has an unallowed value.
          INQUIRE(UNIT=IPR,OPENED=LOPEN)
          IF(LOPEN)WRITE(IPR,'(/4A)')'OPNFL error from ',ROUTIN,        &
     &      ': File FORMAT =',FORMT
          IF(LJMASS)THEN
              CALL WRTBUF(FATAL)
          ELSE
              WRITE(*,'(4A)')'OPNFL error from ',ROUTIN,                &
     &          ': File FORMAT =',FORMT
          ENDIF
          STOP 'OPNFL ERROR: File FORMAT string has an unallowed value.'
      ENDIF

!     IS UNIT NUMBER ALREADY CONNECTED?
      INQUIRE(UNIT=IUNIT,OPENED=LOPEN,NAME=FILNAM)
      IF(LOPEN)THEN

!         FAILURE:  UNIT NUMBER ALREADY CONNECTED.
          LENNAM=MAX(LENSTR(FILNAM),1)
          INQUIRE(UNIT=IPR,OPENED=LOPEN)
          IF(LOPEN)WRITE(IPR,'(/(3A,I3,A))')                            &
     &      'OPNFL error from ',ROUTIN,': Unit number',IUNIT,           &
     &      ' is already connected to file',FILNAM(1:LENNAM)
          IF(LJMASS)THEN
              CALL WRTBUF(FATAL)
          ELSE
              WRITE(*,'(/(3A,I3,A))')                                   &
     &          'OPNFL error from ',ROUTIN,': Unit number',IUNIT,       &
     &          ' is already connected to file',FILNAM(1:LENNAM)
          ENDIF
          STOP 'OPNFL ERROR:  Unit number is already connected.'
      ENDIF

!     IF FILE STATUS IS SCRATCH, OPEN IT NOW:
      IF(STAT.EQ.'SCRATCH')THEN

!         CHECK ACCESS:
          IF(LENREC.LE.0)THEN

!             OPEN SEQUENTIAL ACCESS SCRATCH FILE:
              OPEN(UNIT=IUNIT,STATUS='SCRATCH',FORM=FORMT)
          ELSE

!             OPEN DIRECT ACCESS SCRATCH FILE:
              OPEN(UNIT=IUNIT,STATUS='SCRATCH',                         &
     &          ACCESS='DIRECT',RECL=LENREC,FORM=FORMT)
          ENDIF
          RETURN
      ENDIF

!     CHECK FILE NAME LENGTH:
      LENNAM=LENSTR(FNAME)
      IF(LENNAM.LE.0)THEN

!         FAILURE:  FILE NAME IS BLANK:
          INQUIRE(UNIT=IPR,OPENED=LOPEN)
          IF(LOPEN)WRITE(IPR,'(/3A)')'OPNFL error from ',ROUTIN,        &
     &      ':  File name is blank.'
          IF(LJMASS)THEN
              CALL WRTBUF(FATAL)
          ELSE
              WRITE(*,'(3A)')'OPNFL error from ',ROUTIN,                &
     &          ':  File name is blank and status is not SCRATCH.'
          ENDIF
          STOP 'OPNFL ERROR:  File name is blank.'
      ENDIF

!     CHECK FILE STATUS:
      INQUIRE(FILE=FNAME(1:LENNAM),EXIST=LEXIST,OPENED=LOPEN,NUMBER=NUM)

!     DOES 'NEW' FILE ALREADY EXIST?
      IF(STAT.EQ.'NEW' .AND. LEXIST)THEN

!         FAILURE:  FILE ALREADY EXISTS.
          INQUIRE(UNIT=IPR,OPENED=LOPEN)
          IF(LOPEN)WRITE(IPR,'(5A)')'OPNFL error from ',ROUTIN,         &
     &      ':  "NEW" File "',FNAME(1:LENNAM),'" already exists.'
          IF(LJMASS)THEN
              CALL WRTBUF(FATAL)
          ELSE
              WRITE(*,'(5A)')'OPNFL error from ',ROUTIN,                &
     &          ':  "NEW" File "',FNAME(1:LENNAM),'" already exists.'
          ENDIF
          STOP 'OPNFL ERROR:  File being opened as NEW already exists.'

!     IF FILE IS ALREADY OPEN, JUST REWIND IT:
      ELSEIF(LOPEN)THEN
          IF(IUNIT.NE.NUM)THEN

!             FAILURE: FILE IS CONNECTED TO A DIFFERENT UNIT NUMBER.
              INQUIRE(UNIT=IPR,OPENED=LOPEN)
              IF(LOPEN)WRITE(IPR,'(/4A,2(A,I3))')'OPNFL error from ',   &
     &          ROUTIN,':  File "',FNAME(1:LENNAM),                     &
     &          '" is connected to unit no.',NUM,'instead of',IUNIT
              IF(LJMASS)THEN
                  CALL WRTBUF(FATAL)
              ELSE
                  WRITE(*,'(4A,2(A,I3))')'OPNFL error from ',           &
     &              ROUTIN,':  File "',FNAME(1:LENNAM),                 &
     &              '" is connected to unit no.',NUM,'instead of',IUNIT
              ENDIF
              STOP                                                      &
     &          'OPNFL ERROR: Open file connected to wrong unit number.'
          ENDIF
          REWIND(IUNIT)
          RETURN
      ENDIF

!     DOES 'OLD' FILE NOT EXIST?
      IF(STAT.EQ.'OLD' .AND. .NOT.LEXIST)THEN

!         FAILURE:  FILE DOES NOT EXIST.
          INQUIRE(UNIT=IPR,OPENED=LOPEN)
          IF(LOPEN)WRITE(IPR,'(5A)')'OPNFL error from ',ROUTIN,         &
     &      ':  "OLD" File "',FNAME(1:LENNAM),'" does not exist.'
          IF(LJMASS)THEN
              CALL WRTBUF(FATAL)
          ELSE
              WRITE(*,'(5A)')'OPNFL error from ',ROUTIN,                &
     &          ':  "OLD" File "',FNAME(1:LENNAM),'" does not exist.'
          ENDIF
          STOP 'OPNFL ERROR:  File being opened as OLD does not exist.'
      ENDIF

!     CHECK ACCESS:
      IF(LENREC.LE.0)THEN

!         OPEN SEQUENTIAL ACCESS FILE:
          OPEN(UNIT=IUNIT,FILE=FNAME(1:LENNAM),STATUS=STAT,FORM=FORMT)

!         REMOVE POSSIBLE OLD 'UNKNOWN' FILE.
          IF(STAT.EQ.'UNKNOWN')THEN
              CLOSE(IUNIT,STATUS='DELETE')
              OPEN(UNIT=IUNIT,FILE=FNAME(1:LENNAM),                     &
     &          STATUS='NEW',FORM=FORMT)
          ENDIF
      ELSEIF(STAT.EQ.'UNKNOWN')THEN

!         REMOVE POSSIBLE OLD 'UNKNOWN' FILE.
          OPEN(UNIT=IUNIT,FILE=FNAME(1:LENNAM),STATUS='UNKNOWN')
          CLOSE(IUNIT,STATUS='DELETE')
          OPEN(UNIT=IUNIT,FILE=FNAME(1:LENNAM),STATUS='NEW',            &
     &      ACCESS='DIRECT',RECL=LENREC,FORM=FORMT)
      ELSE

!         OPEN 'NEW' OR 'OLD' DIRECT ACCESS FILE:
          OPEN(UNIT=IUNIT,FILE=FNAME(1:LENNAM),STATUS=STAT,             &
     &      ACCESS='DIRECT',RECL=LENREC,FORM=FORMT)
      ENDIF
      RETURN
      END
