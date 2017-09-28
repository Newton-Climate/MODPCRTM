      SUBROUTINE OPENBM(LENDAT,DATDIR,BMROOT,RESCHR,LABEL)

!     OPENBM OPENS THE BAND MODEL DATA FILES, AND READS HEADER DATA.
!     THE BAND MODEL FILE NAMES, WHICH ARE IN BINARY, ARE:
!       <BMROOT>c.bin FOR LINE-CENTER PARAMETERS
!       <BMROOT>t.bin FOR LINE-TAIL PARAMETERS

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     ARGUMENTS:
!       LENDAT   LENGTH OF DATDIR STRING.
!       DATDIR   NAME OF MODTRAN DATA DIRECTORY.
!       BMROOT   PREFIX OF THE BAND MODEL PARAMETER FILE NAME.
!       RESCHR   BAND MODEL RESOLUTION CHARACTER STRING.
      INTEGER LENDAT
      CHARACTER BMROOT*(NAMLEN),DATDIR*(*),RESCHR*2,LABEL*6

!     COMMONS:
      INCLUDE 'IFIL.h'
      INCLUDE 'BMHEAD.h'

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     FUNCTIONS:
!       IRECLN   RETURNS THE RECORD LENGTH FOR A DIRECT ACCESS FILES
!                CONTAINING A KNOWN NUMBER OF VARIABLES IN EACH RECORD.
      INTEGER IRECLN,LENSTR

!     LOCAL VARIABLES:
!       LOPEN    FILE OPENED LOGICAL.
!       LEXIST   FILE EXISTS LOGICAL.
!       ITEMP    TEMPERATURE GRID INDEX
!       IPRES    PRESSURE GRID INDEX (1 FOR CENTERS, UP TO 2 FOR TAILS)
!       LENFUL   FILE FULL NAME LENGTH.
      LOGICAL LOPEN,LEXIST
      INTEGER ITEMP,IPRES,MAXTL,MAXCN,CNRECL,TLRECL,                    &
     &  LNROOT,NTEMP0,ITLSUB,LENFUL
      REAL TBAND0(MTEMP),EDGEN0,BNDWD0,SLOPE
      CHARACTER FULNAM*(NAMLEN),NAMECN*(NAMLEN),NAMETL*(NAMLEN)

!     DATA:
!       EDGEMN   MINIMUM ACCEPTABLE VALUE FOR DEDGE [CM-1].
!       TOL      BAND MODEL RESOLUTION TOLERANCE [CM-1].
      REAL EDGEMN,TOL
      DATA EDGEMN/.0/,TOL/.01/

!     CLOSE FILE IF UNIT ALREADY CONNECTED:
      INQUIRE(ITBCN,OPENED=LOPEN)
      IF(LOPEN)CLOSE(ITBCN,STATUS='KEEP')
      INQUIRE(ITBTL,OPENED=LOPEN)
      IF(LOPEN)CLOSE(ITBTL,STATUS='KEEP')

!     CHECK THAT FILE EXISTS:
      CALL LCTRIM(BMROOT)
      LNROOT=LENSTR(BMROOT)
      IF(LNROOT.EQ.0)THEN
          NAMECN=DATDIR//'01_2008c.bin'
          NAMETL=DATDIR//'01_2008t.bin'
          LENFUL=LENDAT+12
      ELSE
          NAMECN=BMROOT(1:LNROOT)//'c.bin'
          NAMETL=BMROOT(1:LNROOT)//'t.bin'
          NAMECN=DATDIR//NAMECN(1:LNROOT+5)
          NAMETL=DATDIR//NAMETL(1:LNROOT+5)
          LENFUL=LENDAT+LNROOT+5
      ENDIF
      IF(.NOT.LJMASS)WRITE(IPR,'(/(1X,A))')                             &
     &  'MOLECULAR BAND MODEL DATA FILES',                              &
     &  '-------------------------------',                              &
     &  NAMECN(1:LENFUL),NAMETL(1:LENFUL)
      INQUIRE(FILE=NAMECN,EXIST=LEXIST)
      IF(.NOT.LEXIST)THEN
          WRITE(IPR,'(/3A,/18X,A)')' Error in OPENBM:  The molecular',  &
     &      ' band model data file ',NAMECN(1:LENFUL),' does not exist.'
          STOP 'Error:  Molecular band model data file not found.'
      ENDIF
      INQUIRE(FILE=NAMETL,EXIST=LEXIST)
      IF(.NOT.LEXIST)THEN
          WRITE(IPR,'(/3A,/18X,A)')' Error in OPENBM:  The molecular',  &
     &      ' band model data file ',NAMETL(1:LENFUL),' does not exist.'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP 'OPENBM error: Molecular band model data file not found'
      ENDIF

!     INITIALLY OPEN THE BAND MODEL LINE CENTER DATA FILES WITH
!     DUMMY RECORD LENGTH (8) AND READ THE FIRST VARIABLE (NTEMP0)
!     SO THAT THE REQUIRED RECORD LENGTH CAN BE COMPUTED.  EACH
!     RECORD CONTAINS NTEMP0 S/d's, NTEMP0 1/d's, AND 4 ADDITIONAL
!     VARIABLES: A SPECTRAL FREQUENCY INDEX, A MOLECULAR LABEL
!     AND BOTH SELF AND FOREIGN BROADENED LORENTZ HALF-WIDTHS.
      CALL OPNFL(ITBCN,8,NAMECN,'OLD','UNFORMATTED','OPENBM')
      READ(ITBCN,REC=1)NTEMP0
      CLOSE(ITBCN)
      CNRECL=IRECLN(2*NTEMP0+4,LABEL)
      CALL OPNFL(ITBCN,CNRECL,NAMECN,'OLD','UNFORMATTED','OPENBM')
      READ(ITBCN,REC=1)NTEMP0,NPRESS,BNDWD0,MAXCN,LRECCN,EDGEN0,        &
     &  (TBAND0(ITEMP),ITEMP=1,NTEMP0),(PBAND(IPRES),IPRES=1,NPRESS)

!     INITIALLY OPEN THE BAND MODEL LINE TAIL DATA FILES WITH DUMMY
!     RECORD LENGTH (12) AND READ THE FIRST 2 VARIABLES (NTEMP AND
!     NPRESS) SO THAT THE REQUIRED RECORD LENGTH CAN BE COMPUTED.
!     EACH RECORD CONTAINS NTEMP0*NPRESS SETS OF PADE APPROXIMATE
!     COEFFICIENTS (5 TERMS IN EACH SET), AND 2 ADDITIONAL
!     VARIABLES: A SPECTRAL FREQUENCY INDEX AND A MOLECULAR LABEL.
      CALL OPNFL(ITBTL,12,NAMETL,'OLD','UNFORMATTED','OPENBM')
      READ(ITBTL,REC=1)NTEMP,NPRESS
      CLOSE(ITBTL)
      TLRECL=IRECLN(5*NTEMP*NPRESS+2,LABEL)
      CALL OPNFL(ITBTL,TLRECL,NAMETL,'OLD','UNFORMATTED','OPENBM')
      READ(ITBTL,REC=1)NTEMP,NPRESS,BNDWID,MAXTL,LRECTL,EDGENR,         &
     &  (TBAND(ITEMP),ITEMP=1,NTEMP),(PBAND(IPRES),IPRES=1,NPRESS)
      HBNDWD=BNDWID/2

!     DEFINE UNIFORM GRID BETWEEN -1 AND 1:
      DELSUB(0)=-1.
      SLOPE=2./NTLSUB
      DO ITLSUB=1,NTLSUB-1
          DELSUB(ITLSUB)=ITLSUB*SLOPE-1
      ENDDO
      DELSUB(NTLSUB)=1.

!     CHECK FILE COMPATIBILITY:
      IF(NTEMP0.NE.NTEMP .OR. BNDWD0.NE.BNDWID                          &
     &  .OR. EDGEN0.NE.EDGENR)THEN
          WRITE(IPR,'(/A)')' Error in OPENBM: '//                       &
     &      ' LINE CENTER & TAIL BAND MODEL DATA ARE INCOMPATIBLE'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP 'LINE CENTER & TAIL BAND MODEL DATA ARE INCOMPATIBLE'
      ENDIF
      DO ITEMP=1,NTEMP
          IF(TBAND0(ITEMP).NE.TBAND(ITEMP))THEN
              WRITE(IPR,'(/A)')' Error in OPENBM: '//                   &
     &          ' LINE CENTER & TAIL BAND MODEL DATA ARE INCOMPATIBLE'
              IF(LJMASS)CALL WRTBUF(FATAL)
              STOP 'LINE CENTER & TAIL BAND MODEL DATA ARE INCOMPATIBLE'
          ENDIF
      ENDDO

      IF(BNDWID.LT..2)THEN

!         SHIFT IS AN OFFSET FOR OUTPUT ONLY WITH 0.1 CM-1 BAND MODEL
          OSHIFT=HBNDWD
      ELSE
          OSHIFT=0.
      ENDIF
      MBINBM=MAX(MAXTL,MAXCN)
      MXFREQ=.1*NINT(10*(MBINBM*BNDWID))

!2BIN
!2BIN DOUBLE BIN ALGORITHM:
!2BIN IF(EDGENR.LE.0.)THEN
!2BIN     EDGENR=0.
!2BIN     ADTEST=.2
!2BIN ENDIF

!     SINGLE BIN ALGORITHM:
      EDGENR=EDGENR/BNDWID
      IF(EDGENR.LT.EDGEMN)THEN
          EDGENR=EDGEMN
          EDGEFR=1.-EDGENR
          EDGEDF=EDGEFR-EDGENR
      ELSEIF(EDGENR.GT..5)THEN
          EDGENR=.5
          EDGEFR=.5
          EDGEDF=.0
      ELSE
          EDGEFR=1.-EDGENR
          EDGEDF=EDGEFR-EDGENR
      ENDIF
      ADTEST=.2*EDGENR

!     OPEN THE FORMATTED BAND MODEL FILE, UNIT ITBX.
      INQUIRE(ITBX,OPENED=LOPEN)
      IF(LOPEN)CLOSE(ITBX,STATUS='KEEP')
      IF(ABS(BNDWID-15.).LT.TOL)THEN

!         CHECK THAT FILE EXISTS:
          RESCHR='15'
          FULNAM=DATDIR//'CFC04_15.ASC'
      ELSEIF(ABS(BNDWID-5.).LT.TOL)THEN

!         CHECK THAT FILE EXISTS:
          RESCHR='05'
          FULNAM=DATDIR//'CFC04_05.ASC'
      ELSEIF(ABS(BNDWID-1.).LT.TOL)THEN

!         CHECK THAT FILE EXISTS:
          RESCHR='01'
          FULNAM=DATDIR//'CFC04_01.ASC'
      ELSEIF(ABS(BNDWID-.1).LT.TOL)THEN

!         CHECK THAT FILE EXISTS:
          RESCHR='p1'
          FULNAM=DATDIR//'CFC04_p1.ASC'
      ENDIF
      LENFUL=LENSTR(FULNAM)
      INQUIRE(FILE=FULNAM,EXIST=LEXIST)
      IF(.NOT.LEXIST)THEN
          WRITE(IPR,'(/2A,/18X,3A)')' Error in OPENBM: ',               &
     &      ' The CFC molecular band model data file',                  &
     &      ' "',FULNAM(1:LENFUL),'" does not exist.'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP 'OPENBM error: CFC molecular band model file not found.'
      ENDIF
      IF(.NOT.LJMASS)WRITE(IPR,'(/A,15X,A)')                            &
     &  ' CFC BAND MODEL DATA FILE:',FULNAM
      CALL OPNFL(ITBX,0,FULNAM,'OLD','FORMATTED','OPENBM')

!     WRITE WATER CONTINUUM INFORMATION:
      IF(.NOT.LJMASS)WRITE(IPR,'(/A)')' Version 2.4 of the'//           &
     &  ' Clough-Kneizys Water Continuum Data from LBLRTM (24mar2000).'

!     RETURN TO DRIVER:
      RETURN
      END
