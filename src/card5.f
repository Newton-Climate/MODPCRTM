      BLOCK DATA CRD5BD

!     COMMONS:

!     /JM5/
!       IRPT     REPEAT INPUT FLAG (0=NONE, 1=ALL, 3=GEOM, 4=SPEC).
!       IFAC     CURRENT COLUMN SCALING FACTOR INDEX.
!       NFACMN   NUMBER OF COLUMN SCALING FACTOR LESS THAN 1.
!       NFACMX   NUMBER OF COLUMN SCALING FACTOR GREATER THAN 1.
!       FACMC    CURRENT COLUMN SCALING FACTOR.
!       SCALMN   MINIMUM COLUMN SCALING FACTOR.
!       SCALMX   MAXIMUM COLUMN SCALING FACTOR.
!       LBMRES   LOGICAL FLAG, .TRUE. FOR BAND MODEL AND .FALSE.
!                FOR 1 NM SPECTRAL RESOLUTION OUTPUT.
      INTEGER IRPT,IFAC,NFACMN,NFACMX
      REAL FACMC
      DOUBLE PRECISION SCALMN,SCALMX
      LOGICAL LBMRES
      COMMON/JM5/SCALMN,SCALMX,IRPT,IFAC,NFACMN,NFACMX,FACMC,LBMRES

!     /CJM5/
!       AMOD3D   FLAG INDICATING OUTPUT DATABASE FILE TYPE:
      CHARACTER AMOD3D*1
      COMMON/CJM5/AMOD3D

!     DATA:
      DATA IRPT,IFAC,NFACMN,NFACMX/4*0/,FACMC/1./,SCALMN,SCALMX/2*1.D0/
      DATA AMOD3D/' '/
      END

      SUBROUTINE CARD5(LNFLRT,FLRT)

!     ROUTINE TO READ IN AND PROCESS CARD5 INPUTS:
      IMPLICIT NONE

!     PARAMETERS:
!       ONE      THE NUMBER ONE IN DOUBLE PRECISION.
      DOUBLE PRECISION ONE
      PARAMETER(ONE=1.D0)
      INCLUDE 'PARAMS.h'

!     INPUT ARGUMENTS:
!       LNFLRT   LENGTH OF ROOT NAME FOR ALL I/O FILES.
!       FLRT     ROOT NAME FOR ALL I/O FILES.
      INTEGER LNFLRT
      CHARACTER FLRT*(NAMLEN-4)

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /JM5/
!       IRPT     REPEAT INPUT FLAG (0=NONE, 1=ALL, 3=GEOM, 4=SPEC).
!       IFAC     CURRENT COLUMN SCALING FACTOR INDEX.
!       NFACMN   NUMBER OF COLUMN SCALING FACTOR LESS THAN 1.
!       NFACMX   NUMBER OF COLUMN SCALING FACTOR GREATER THAN 1.
!       SCALMN   MINIMUM COLUMN SCALING FACTOR.
!       FACMC    CURRENT COLUMN SCALING FACTOR.
!       SCALMX   MAXIMUM COLUMN SCALING FACTOR.
!       LBMRES   LOGICAL FLAG, .TRUE. FOR BAND MODEL AND .FALSE.
!                FOR 1 NM SPECTRAL RESOLUTION OUTPUT.
      INTEGER IRPT,IFAC,NFACMN,NFACMX
      REAL FACMC
      DOUBLE PRECISION SCALMN,SCALMX
      LOGICAL LBMRES
      COMMON/JM5/SCALMN,SCALMX,IRPT,IFAC,NFACMN,NFACMX,FACMC,LBMRES

!     /CJM5/
!       AMOD3D   FLAG INDICATING OUTPUT DATABASE FILE TYPE:
      CHARACTER AMOD3D*1
      COMMON/CJM5/AMOD3D

!     /CNTRL/
!       NSEG     NUMBER OF PATH SEGMENTS ALONG LINE-OF-SIGHT.
!       ML       NUMBER OF ATMOSPHERIC PROFILE LEVELS.
!       MLFLX    NUMBER OF LEVELS FOR WHICH FLUX VALUES ARE WRITTEN.
!       IMULT    MULTIPLE SCATTERING FLAG
!                  (0=NONE, 1=AT SENSOR, -1=AT FINAL OR TANGENT POINT).
!       THERML   FLAG TO CALCULATE THERMAL SCATTER.
      INTEGER NSEG,ML,MLFLX,IMULT
      LOGICAL THERML
      COMMON/CNTRL/NSEG(0:MLOSP1),ML,MLFLX,IMULT,THERML

!     FUNCTIONS:
!       IRECLN   RETURNS THE RECORD LENGTH FOR A DIRECT ACCESS FILES
!                CONTAINING A KNOWN NUMBER OF VARIABLES IN EACH RECORD.
!       NUNIT    RETURNS AN UNUSED FILE UNIT NUMBER.
      INTEGER IRECLN,NUNIT

!     LOCAL VARIABLES:
!       IOS      RESULT OF IOSTAT CHECK.
!       LRECLN   RECORD LENGTH OF DIRECT ACCESS FILES.
!       LENDIR   LENGTH OF OUTPUT DIRECTORY NAME.
!       LROOTS   UNIT NUMBER FOR ROOT NAME OUTPUT FILE.
!       LOPEN    LOGICAL FLAG, TRUE IF FILE IS OPEN.
!       EXTEN    FILE EXTENSION NAME.
!       CHAR80   INPUT STRING.
!       FULNAM   FULL PATH NAME.
      INTEGER IOS,LRECLN,LENDIR,LROOTS
      LOGICAL LOPEN
      CHARACTER EXTEN*8,CHAR80*80,FULNAM*(NAMLEN+5)

!     SAVED VARIABLES:
!       SNUMER   NUMERATOR FACTOR USED TO DEFINE COLUMN SCALING.
!       SDENOM   DENOMINATOR FACTOR USED TO DEFINE COLUMN SCALING.
      DOUBLE PRECISION SNUMER,SDENOM
      SAVE SNUMER,SDENOM

!     CHECK FOR NEW READ:
      IF(IFAC.LT.NFACMX)THEN
          IFAC=IFAC+1
      ELSE
          READ(IRD,'(A80)')CHAR80
          READ(CHAR80,'(I5,50X,A1,2(I3,F9.0))',IOSTAT=IOS)              &
     &      IRPT,AMOD3D,NFACMN,SCALMN,NFACMX,SCALMX
          IF(IOS.NE.0)THEN

!             OLD FORMAT:
              READ(CHAR80(1:5),'(I5)',IOSTAT=IOS)IRPT
              IF(IOS.NE.0)THEN

!                 UNABLE TO READ CARD5:
                  WRITE(IPR,'(/3A)')' Error in CARD5: '                 &
     &              //' Unable to read CARD5 "',CHAR80(1:5),'"'
                  STOP 'Error in CARD5:  Unable to read CARD5!'
              ENDIF
              AMOD3D=' '
              FACMC=1.
              IFAC=0
          ELSEIF(AMOD3D.EQ.'T' .OR. AMOD3D.EQ.'t')THEN

!             MOD3D DATABASE:
              AMOD3D='T'
              FACMC=1.
              NFACMX=0
              IFAC=0
          ELSEIF(AMOD3D.EQ.'C' .OR. AMOD3D.EQ.'c')THEN

!             MC CONTINUA DATABASE:
              AMOD3D='C'
              LRECLN=IRECLN(10*NSEG(1))
              IF(LNFLRT.LE.0)THEN
                  CALL OPNFL(IDBOUT,LRECLN,'mc.bin','UNKNOWN',          &
     &              'UNFORMATTED','CARD5')
              ELSE
                  CALL OPNFL(IDBOUT,LRECLN,FLRT(1:LNFLRT)//'.mcb',      &
     &              'UNKNOWN','UNFORMATTED','CARD5')
                  DO LENDIR=LNFLRT-1,1,-1

!                     CHAR(92) IS A BACKSLASH, THE PC DIRECTORY SYMBOL.
                      IF(FLRT(LENDIR:LENDIR).EQ.CHAR(92) .OR.           &
     &                   FLRT(LENDIR:LENDIR).EQ.'/')THEN

!                         CHARACTER STRING "FULNAM" IS DEFINED BECAUSE
!                         STRING CONCATENATION FOR FILE SPECIFICATION
!                         IN OPEN STATEMENTS IS NOT SUPPORTED BY THE
!                         G77 COMPILER!
                          FULNAM=FLRT(1:LENDIR)//'roots.txt'
                          INQUIRE(FILE=FULNAM,                          &
     &                      OPENED=LOPEN,NUMBER=LROOTS)
                          IF(.NOT.LOPEN)THEN
                              LROOTS=NUNIT()
                              CALL OPNFL(LROOTS,0,FULNAM,'UNKNOWN',     &
     &                          'FORMATTED','CARD5')
                          ENDIF
                          WRITE(LROOTS,'(A)')FLRT(LENDIR+1:LNFLRT)
                          GOTO 10
                      ENDIF
                  ENDDO
              ENDIF
   10         CONTINUE
              NFACMX=0
              IFAC=0

!             WRITE NM BIN DATA IF NFACMN=0; ELSE, BAND MODEL BIN DATA:
              LBMRES=NFACMN.NE.0
          ELSEIF(AMOD3D.EQ.'M' .OR. AMOD3D.EQ.'m')THEN

!             MC MOLECULAR TRANSMITTANCE DATABASE:
              AMOD3D='M'
              IFAC=-NFACMN
              IF(NFACMN.LE.0 .OR. SCALMN.GE.ONE .OR.                    &
     &           NFACMX.LE.0 .OR. SCALMX.LE.ONE)STOP 'Bad CARD5 input'
              SDENOM=NFACMN*NFACMX*(SCALMX-SCALMN)
              SNUMER=(SCALMX*(ONE-SCALMN)*NFACMX                        &
     &               -SCALMN*(SCALMX-ONE)*NFACMN)/SDENOM
              SDENOM=((ONE-SCALMN)*NFACMX-(SCALMX-ONE)*NFACMN)/SDENOM
              INQUIRE(IDBOUT,OPENED=LOPEN)
              IF(.NOT.LOPEN)THEN
                  IF(LNFLRT.LE.0)THEN
                      CALL OPNFL(IDBOUT,0,'mc.dat',                     &
     &                  'UNKNOWN','FORMATTED','CD5')
                  ELSE
                      CALL OPNFL(IDBOUT,0,FLRT(1:LNFLRT)//'.mc',        &
     &                  'UNKNOWN','FORMATTED','CD5')
                  ENDIF
              ENDIF
              WRITE(IDBOUT,'(/(I13,A,/F13.4,A))')                       &
     &          NFACMN,' NUMBER OF SCALE FACTORS LESS THAN 1.',         &
     &          SCALMN,' MINIMUM SCALE FACTOR',                         &
     &          NFACMX,' NUMBER OF SCALE FACTORS GREATER THAN 1.',      &
     &          SCALMX,' MAXIMUM SCALE FACTOR'
          ELSE
              AMOD3D=' '
              FACMC=1.
              IFAC=0
          ENDIF
      ENDIF
      IF(AMOD3D.EQ.'M')THEN
          FACMC=REAL((ONE+IFAC*SNUMER)/(ONE+IFAC*SDENOM))
          REWIND(IPR)
          REWIND(IPU)
          CLOSE(IDBOUT)
          EXTEN(1:2)='.M'
          WRITE(EXTEN(3:8),'(I6.6)')INT(100*FACMC+.5)
          LRECLN=IRECLN(8*NSEG(1)+8)
          IF(LNFLRT.LE.0)THEN
              CALL OPNFL(IDBOUT,LRECLN,EXTEN(2:8),'UNKNOWN',            &
     &          'UNFORMATTED','CARD5')
          ELSE

!                 CHARACTER STRING "FULNAM" IS DEFINED BECAUSE STRING
!                 CONCATENATION FOR FILE SPECIFICATION IN OPEN
!                 STATEMENTS IS NOT SUPPORTED BY THE G77 COMPILER!
              FULNAM=FLRT(1:LNFLRT)//EXTEN
              CALL OPNFL(IDBOUT,LRECLN,FULNAM,'UNKNOWN','UNFORMATTED',  &
     &          'CARD5')
          ENDIF
      ENDIF
      RETURN
      END
