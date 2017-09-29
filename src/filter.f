      LOGICAL FUNCTION FILTER(FILTNM,LNFILT,NMFLRT,NLOS)

!     THIS LOGICAL FUNCTION READS THE NAME OF A FILTER
!     FUNCTION FILE, OPENS THE FILE, AND PROCESSES THE FILTER
!     FOR EACH BAND TO ENABLE RAPID PROCESSING.  IF AN ERROR
!     OCCURS, FILTER IS RETURNED WITH A VALUE OF .FALSE.
      IMPLICIT NONE

!     INPUT ARGUMENTS:
!       FILTNM   FILTER FUNCTION FILE NAME.
!       LNFILT   LENGTH OF FILTER FUNCTION FILE NAME.
!       NMFLRT   CASE NUMBER.
!       NLOS     NUMBER OF LINE-OF-SIGHT PATHS.
      CHARACTER FILTNM*(*)
      INTEGER LNFILT,NMFLRT,NLOS

!     PARAMETERS:
      INCLUDE 'PARAMS.h'
      INCLUDE 'CHANLS.h'

!     COMMONS:
      INCLUDE 'IFIL.h'
      INCLUDE 'BMHEAD.h'

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     FUNCTIONS:
      INTEGER LENSTR,NUNIT
      LOGICAL WTCHAN

!     LOCAL VARIABLES:
!       FRMT     FORMAT STRING.
!       ICOMMA   LOCATION OF COMMA IN CHARACTER STRING.
!       IEND     LOCATION OF END CHARACTER IN STRING.
!       IBLANK   LOCATION OF BLANK IN CHARACTER STRING.
!       IOUT     CHANNEL SPECTRAL OUTPUT INDEX.
!       CONVRT   WAVELENGTH TO CM-1 CONVERSION FACTOR OVER THE
!                CALCULATION BIN WIDTH [NM OR MICRONS].
!       NFLRTL   LAST FILE NUMBER PROCESSED.
!       ILOS     LOOP INDEX FOR LINE-OF-SIGHT PATHS.
!       SPEC_P   PREVIOUS SPECTRAL POINT.
!       SPEC_C   CURRENT SPECTRAL POINT.
!       WGT_C    CHANNEL RESPONSE VALUE AT SPEC_C.
!       SPEC_N   NEXT SPECTRAL POINT.
!       WGT_N    CHANNEL RESPONSE VALUE AT SPEC_N.
!       WGT      NORMALIZATION WEIGHT.
!       NORM     NORMALIZATION SUM.
!       MOM2ND   CHANNEL SECOND SPECTRAL MOMENT [CM-2, UM^2 or NM^2].
!       MOM3RD   CHANNEL THIRD SPECTRAL MOMENT [CM-2, UM^2 or NM^2].
!       COEF     COEFFICIENT USED IN SPECTRAL MOMENT CALCULATIONS.
!       COEF2    COEFFICIENT USED IN SPECTRAL MOMENT CALCULATIONS.
      CHARACTER INLINE*80,FRMT*7
      INTEGER IOTEST,IFILTR,LNTEST,LSTBIN,NSPEC,IBIN,NBNDWD,            &
     &  ICOMMA,IEND,IBLANK,IOUT,NFLRTL,ILOS
      LOGICAL LOPEN,LEXIST
      REAL SPEC_P,SPEC_C,WGT_C,SPEC_N,WGT_N,CONVRT,WGT,NORM,            &
     &  MOM2ND,MOM3RD,COEF,COEF2

!     DATA:
!       DATTST   CHARACTER STRING CONTAINING POSSIBLE INITIAL NON-BLANK
!                CHARACTERS FOR FREQUENCY/WAVELENGTH AND WEIGHT LINES.
      CHARACTER DATTST*11,FLTNMS*(NAMLEN)
      LOGICAL NOT0LO,NOT0HI
      SAVE NOT0LO,NOT0HI,NBNDWD,FLTNMS,NFLRTL
      DATA DATTST/'.0123456789'/,NOT0LO,NOT0HI/2*.TRUE./,NBNDWD/0/,     &
     &  FLTNMS/' '/,NFLRTL/0/

      IF(NMFLRT.NE.NFLRTL)THEN
          NFLRTL=NMFLRT
          NOT0LO=.TRUE.
          NOT0HI=.TRUE.
      ENDIF

!     INITIALIZE FILTER TO .FALSE.:
      FILTER=.FALSE.

!     CHECK IF FILTER FILE HAS ALREADY BEEN PROCESSED
!     WITH CURRENT SPECTRAL FREQUENCY STEP SIZE.
      IF(FILTNM.EQ.FLTNMS .AND. BNDWID.EQ.NBNDWD)THEN
          WRITE(IPR,'(2(/2A),F6.3,A)')' FILTER FUNCTION FILE  ',        &
     &      FILTNM(1:LNFILT),' HAS ALREADY BEEN PROCESSED',             &
     &      ' WITH FREQUENCY STEP SIZE',BNDWID,' CM-1.'
          RETURN
      ENDIF
      FLTNMS=' '

!     CHECK FILTER FILE STATUS:
      IF(LNFILT.LE.0)THEN
          WRITE(IPR,'(/A)')                                             &
     &      ' WARNING:  The input filter file name is BLANK.'
          RETURN
      ENDIF
      INQUIRE(FILE=FILTNM(1:LNFILT),EXIST=LEXIST,OPENED=LOPEN)
      IF(.NOT.LEXIST)THEN
          WRITE(IPR,'(/A,/(11X,A))')' WARNING:  The filter file',       &
     &      FILTNM(1:LNFILT),'does not exist.'
          LNFILT=0
          RETURN
      ELSEIF(LOPEN)THEN
          WRITE(IPR,'(/A,/(11X,A))')' WARNING:  The filter file',       &
     &      FILTNM(1:LNFILT),'is already opened.'
          LNFILT=0
          RETURN
      ENDIF

!     OPEN FILTER FUNCTION:
      IFILTR=NUNIT()
      CALL OPNFL(IFILTR,0,FILTNM(1:LNFILT),'OLD','FORMATTED','FILTER')
      WRITE(IPR,'(/2A)')' Opened filter file:  ',FILTNM(1:LNFILT)

!     READ UNIT FLAG:
      READ(IFILTR,'(A80)',IOSTAT=IOTEST)INLINE
      IF(IOTEST.GT.0)THEN
          WRITE(IPR,'(/A)')                                             &
     &      ' WARNING:  Error reading FILTER FILE UNIT FLAG.'
          CLOSE(IFILTR,STATUS='KEEP')
          LNFILT=0
          RETURN
      ELSEIF(IOTEST.LT.0)THEN
          WRITE(IPR,'(/2A)')' WARNING:  End-Of-File',                   &
     &      ' occurred during read of FILTER FILE UNIT FLAG.'
          CLOSE(IFILTR,STATUS='KEEP')
          LNFILT=0
          RETURN
      ENDIF
      LNTEST=LENSTR(INLINE)
      IF(INLINE(1:1).EQ.'M' .OR. INLINE(1:1).EQ.'m')THEN

!         MICRONS:
          UNTFLG='M'
          CONVRT=1.E4/BNDWID
      ELSEIF(INLINE(1:1).EQ.'N' .OR. INLINE(1:1).EQ.'n')THEN

!         NANOMETERS:
          UNTFLG='N'
          CONVRT=1.E7/BNDWID
      ELSE

!         WAVENUMBERS:
          IF(INLINE(1:1).NE.'W' .AND. INLINE(1:1).NE.'w')               &
     &      WRITE(IPR,'(/(A))')                                         &
     &        ' WARNING:  Filter file unit flag does not begin with',   &
     &        '           "W", but WAVENUMBERS is being assumed.'
          UNTFLG='W'
      ENDIF

!     INITIALIZE NUMCHN ARRAY TO ZERO.
      DO IBIN=MNBIN,MXBIN
          NUMCHN(IBIN)=0
      ENDDO

!     READ TITLE OF FIRST CHANNEL AND THE FIRST TWO FREQUENCY/WAVELENGTH
!     AND WEIGHT PAIRS, AND THEN PROCESS THOSE WEIGHTS.
      NCHAN=1
      READ(IFILTR,'(A80)',IOSTAT=IOTEST)NMCHAN(1)
      IF(IOTEST.GT.0)THEN
          WRITE(IPR,'(/A)')' WARNING:  Error reading CHANNEL   1 NAME'
          CALL RMCHAN(IFILTR)
          RETURN
      ELSEIF(IOTEST.LT.0)THEN
          WRITE(IPR,'(/2A)')' WARNING:  End-Of-File',                   &
     &      ' occurred during read of CHANNEL   1 NAME'
          CALL RMCHAN(IFILTR)
          RETURN
      ENDIF
      LNCHAN(1)=LENSTR(NMCHAN(1))
      READ(IFILTR,*,IOSTAT=IOTEST)SPEC_C,WGT_C
      IF(IOTEST.GT.0)THEN
          WRITE(IPR,'(/2A)')' WARNING:  Error reading data',            &
     &      ' pair   1 for channel   1 in filter function file.'
          CALL RMCHAN(IFILTR)
          RETURN
      ELSEIF(IOTEST.LT.0)THEN
          WRITE(IPR,'(/2A)')' WARNING:  End-Of-File - data',            &
     &      ' pair   1 for channel   1 in filter function file.'
          CALL RMCHAN(IFILTR)
          RETURN
      ELSEIF(WGT_C.NE.0. .AND. NOT0LO)THEN
          WRITE(IPR,'(/2A,/(11X,2A))')' WARNING: ',                     &
     &      ' First filter weight for channel   1,',                    &
     &      NMCHAN(1)(1:LNCHAN(1)),',','is non-zero.'
          NOT0LO=.FALSE.
      ENDIF
      READ(IFILTR,*,IOSTAT=IOTEST)SPEC_N,WGT_N
      IF(IOTEST.GT.0)THEN
          WRITE(IPR,'(/2A)')' WARNING:  Error reading data',            &
     &      ' pair   2 for channel   1 in filter function file.'
          CALL RMCHAN(IFILTR)
          RETURN
      ELSEIF(IOTEST.LT.0)THEN
          WRITE(IPR,'(/2A)')' WARNING:  End-Of-File - data',            &
     &      ' pair   2 for channel   1 in filter function file.'
          CALL RMCHAN(IFILTR)
          RETURN
      ENDIF
      NSPEC=2
      LSTBIN=-1
      WDFREQ(1)=0.
      WDWAVE(1)=0.
      IF(.NOT.WTCHAN(LSTBIN,SPEC_C,WGT_C,SPEC_N,WGT_N))THEN
          WRITE(IPR,'(/(A))')                                           &
     &      ' WARNING:  Routine WTCHAN error for pairs   1 and   2',    &
     &      '           of channel   1 in filter function file.'
          CALL RMCHAN(IFILTR)
          RETURN
      ENDIF
      DO ILOS=1,NLOS
          DO IOUT=1,MOUT
              VLCHAN(IOUT,ILOS,1)=0.D0
          ENDDO
      ENDDO
      SPECLO(1)=SPEC_C
      IF(UNTFLG.EQ.'W')THEN
          NFRQLO(1)=INT(BNDWID*NINT(SPEC_C/BNDWID))
      ELSE
          NFRQHI(1)=INT(BNDWID*NINT(CONVRT/SPEC_C))
      ENDIF

!     INITIALIZE SPECTRAL MOMENTS:
      WGT=(SPEC_N-SPEC_C)*WGT_C
      NORM=WGT
      COEF=SPEC_N+2*SPEC_C
      MOM1ST(1)=COEF*WGT
      COEF=SPEC_N*COEF+3*SPEC_C**2
      MOM2ND=COEF*WGT
      COEF=SPEC_N*COEF+4*SPEC_C**3
      MOM3RD=COEF*WGT

!     READ AND PARSE NEXT LINE OF FILTER FILE.  THE FIRST NON-BLANK
!     CHARACTER OF A LINE CONTAINING A FREQUENCY/WAVELENGTH AND WEIGHT
!     MUST BE A "-", A ".", OR A NUMERIC CHARACTER (0, 1, ..., 9).
   10 CONTINUE
          READ(IFILTR,'(A80)',IOSTAT=IOTEST)INLINE
          IF(IOTEST.GT.0)THEN

!             ERROR READING FILTER RESPONSE FUNCTION FILE.
              WRITE(IPR,'(/2A,/(11X,A,I5,A))')' WARNING:  Error',       &
     &          ' reading filter function file;  the last line',        &
     &          'successfully read was from CHANNEL',NCHAN,             &
     &          ' entitled',NMCHAN(NCHAN)(1:LNCHAN(NCHAN))
              CALL RMCHAN(IFILTR)
              RETURN
          ELSEIF(IOTEST.LT.0)THEN

!             END-OF-FILE - NORMAL TERMINATION.
              IF(WGT_N.NE.0. .AND. NOT0HI)THEN
                  WRITE(IPR,'(/2A,I5,A,/(11X,2A))')' WARNING: ',        &
     &              ' Last filter weight for channel',NCHAN,',',        &
     &              NMCHAN(NCHAN)(1:LNCHAN(NCHAN)),',','is non-zero.'
                  NOT0HI=.FALSE.
              ENDIF
              SPECHI(NCHAN)=SPEC_N
              IF(UNTFLG.EQ.'W')THEN
                  NFRQHI(NCHAN)=INT(BNDWID*NINT(SPEC_N/BNDWID))
              ELSE
                  NFRQLO(NCHAN)=INT(BNDWID*NINT(CONVRT/SPEC_N))
              ENDIF

!             FINALIZE SPECTRAL MOMENTS:
              WGT=(SPEC_N-SPEC_C)*WGT_N
              NORM=NORM+WGT
              COEF=2*SPEC_N+SPEC_C
              MOM1ST(NCHAN)=(MOM1ST(NCHAN)+COEF*WGT)/(3*NORM)
              COEF=3*SPEC_N**2+SPEC_C*COEF
              MOM2ND=(MOM2ND+COEF*WGT)/(6*NORM)
              COEF=4*SPEC_N**3+SPEC_C*COEF
              MOM3RD=(MOM3RD+COEF*WGT)/(10*NORM)

!             TEST BRIGHTNESS TEMPERATURE CALCULATION:
              M2DIFF(NCHAN)=MOM2ND/MOM1ST(NCHAN)**2-1
              IF(ABS(MOM3RD/MOM1ST(NCHAN)**3-3*M2DIFF(NCHAN)-1).GT..01  &
     &          .OR. ABS(M2DIFF(NCHAN)).GT..02)THEN
                  WRITE(IPR,'(/A,I5,A)')' Warning from routine'//       &
     &              ' FILTER:  Brightness temperature for channel',     &
     &              NCHAN,' may be inaccurate!'
                  C_STAR(NCHAN)='*'
              ELSE
                  C_STAR(NCHAN)=' '
              ENDIF

!             FINISH UP:
              FILTER=.TRUE.
              FLTNMS=FILTNM
              NBNDWD=INT(BNDWID)
              WRITE(IPR,'(/I5,(A))')NCHAN,' channels were read from'    &
     &          //' the filter function file',FILTNM(1:LNFILT)
              CLOSE(IFILTR,STATUS='KEEP')
              RETURN
          ENDIF
          LNTEST=LENSTR(INLINE)
          IF(INDEX(DATTST,INLINE(1:1)).GT.0)THEN

!             INLINE CONTAINS FREQUENCY/WAVELENGTH AND WEIGHT.
              SPEC_P=SPEC_C
              SPEC_C=SPEC_N
              WGT_C=WGT_N
              NSPEC=NSPEC+1

!             SINCE READING FROM A CHARACTER STRING USING A FREE FORMAT,
!             I.E.  "READ(INLINE,FMT='(*)',IOSTAT=IOTEST)SPEC_N,WGT_N"
!             IS NOT ASCII STANDARD AND NOT PERMITTED FOR SOME FORTRAN
!             SOME COMPILERS, INLINE IS PARSED FOR TWO NUMBERS.
              ICOMMA=INDEX(INLINE,',')
              IF(ICOMMA.GT.0)INLINE(ICOMMA:ICOMMA)=' '
              IEND=INDEX(INLINE,' ')-1
              IF(IEND.LE.0)THEN
                  WRITE(IPR,'(/3(A,I5))')                               &
     &              ' WARNING:  Error reading data pair',NSPEC,         &
     &              ' for channel',NCHAN,' in filter function file.'
                  CALL RMCHAN(IFILTR)
                  RETURN
              ENDIF
              FRMT='(F00.0)'
              WRITE(FRMT(3:4),'(I2.2)')IEND
              READ(INLINE(1:IEND),FMT=FRMT,IOSTAT=IOTEST)SPEC_N
              IF(IOTEST.GT.0)THEN
                  WRITE(IPR,'(/3(A,I5))')                               &
     &              ' WARNING:  Error reading data pair',NSPEC,         &
     &              ' for channel',NCHAN,' in filter function file.'
                  CALL RMCHAN(IFILTR)
                  RETURN
              ELSEIF(IOTEST.LT.0)THEN
                  WRITE(IPR,'(/3(A,I5))')                               &
     &              ' WARNING:  End-Of-File - data pair',NSPEC,         &
     &              ' for channel',NCHAN,' in filter function file.'
                  CALL RMCHAN(IFILTR)
                  RETURN
              ENDIF

!             INCREMENT SPECTRAL MOMENTS:
              WGT=(SPEC_N-SPEC_P)*WGT_C
              NORM=NORM+WGT
              COEF=SPEC_N+SPEC_C+SPEC_P
              MOM1ST(NCHAN)=MOM1ST(NCHAN)+COEF*WGT
              COEF2=SPEC_N**2+SPEC_C**2+SPEC_P**2
              MOM2ND                                                    &
     &          =MOM2ND+(COEF2+SPEC_N*(SPEC_C+SPEC_P)+SPEC_C*SPEC_P)*WGT
              MOM3RD=MOM3RD+(COEF*COEF2+SPEC_N*SPEC_C*SPEC_P)*WGT

!             ELIMINATE SPECTRAL VALUE FROM INLINE CHARACTER STRING:
              DO IBLANK=1,IEND
                  INLINE(IBLANK:IBLANK)=' '
              ENDDO
              LNTEST=LENSTR(INLINE)
              IF(LNTEST.LE.0)THEN
                  WRITE(IPR,'(/3(A,I5))')                               &
     &              ' WARNING:  Error reading data pair',NSPEC,         &
     &              ' for channel',NCHAN,' in filter function file.'
                  CALL RMCHAN(IFILTR)
                  RETURN
              ENDIF
              ICOMMA=INDEX(INLINE,',')
              IF(ICOMMA.GT.0)INLINE(ICOMMA:ICOMMA)=' '
              IEND=INDEX(INLINE,' ')-1
              IF(IEND.LE.0)THEN
                  WRITE(IPR,'(/3(A,I5))')                               &
     &              ' WARNING:  Error reading data pair',NSPEC,         &
     &              ' for channel',NCHAN,' in filter function file.'
                  CALL RMCHAN(IFILTR)
                  RETURN
              ENDIF
              WRITE(FRMT(3:4),'(I2.2)')IEND
              READ(INLINE(1:IEND),FMT=FRMT,IOSTAT=IOTEST)WGT_N
              IF(IOTEST.GT.0)THEN
                  WRITE(IPR,'(/3(A,I5))')                               &
     &              ' WARNING:  Error reading data pair',NSPEC,         &
     &              ' for channel',NCHAN,' in filter function file.'
                  CALL RMCHAN(IFILTR)
                  RETURN
              ELSEIF(IOTEST.LT.0)THEN
                  WRITE(IPR,'(/3(A,I5))')                               &
     &              ' WARNING:  End-Of-File - data pair',NSPEC,         &
     &              ' for channel',NCHAN,' in filter function file.'
                  CALL RMCHAN(IFILTR)
                  RETURN
              ENDIF
          ELSEIF(NCHAN.LT.MXCHAN)THEN

!             SAVE PREVIOUS CHANNEL MAXIMUM SPECTRAL POINT,
!             INCREMENT CHANNEL, ASSIGN TITLE, AND READ AND CHECK
!             THE FIRST TWO FREQUENCY/WAVELENGTH AND WEIGHT PAIRS.
              IF(WGT_N.NE.0. .AND. NOT0HI)THEN
                  WRITE(IPR,'(/2A,I5,A,/(11X,2A))')' WARNING: ',        &
     &              ' Last filter weight for channel',NCHAN,',',        &
     &              NMCHAN(NCHAN)(1:LNCHAN(NCHAN)),',','is non-zero.'
                  NOT0HI=.FALSE.
              ENDIF
              SPECHI(NCHAN)=SPEC_N
              IF(UNTFLG.EQ.'W')THEN
                  NFRQHI(NCHAN)=INT(BNDWID*NINT(SPEC_N/BNDWID))
              ELSE
                  NFRQLO(NCHAN)=INT(BNDWID*NINT(CONVRT/SPEC_N))
              ENDIF

!             FINALIZE SPECTRAL MOMENTS:
              WGT=(SPEC_N-SPEC_C)*WGT_N
              NORM=NORM+WGT
              COEF=2*SPEC_N+SPEC_C
              MOM1ST(NCHAN)=(MOM1ST(NCHAN)+COEF*WGT)/(3*NORM)
              COEF=3*SPEC_N**2+SPEC_C*COEF
              MOM2ND=(MOM2ND+COEF*WGT)/(6*NORM)
              COEF=4*SPEC_N**3+SPEC_C*COEF
              MOM3RD=(MOM3RD+COEF*WGT)/(10*NORM)

!             TEST BRIGHTNESS TEMPERATURE CALCULATION:
              M2DIFF(NCHAN)=MOM2ND/MOM1ST(NCHAN)**2-1
              IF(ABS(MOM3RD/MOM1ST(NCHAN)**3-3*M2DIFF(NCHAN)-1).GT..01  &
     &          .OR. ABS(M2DIFF(NCHAN)).GT..02)THEN
                  WRITE(IPR,'(/A,I5,A)')' Warning from routine'//       &
     &              ' FILTER:  Brightness temperature for channel',     &
     &              NCHAN,' may be inaccurate!'
                  C_STAR(NCHAN)='*'
              ELSE
                  C_STAR(NCHAN)=' '
              ENDIF

!             NEXT CHANNEL:
              NCHAN=NCHAN+1
              NMCHAN(NCHAN)=INLINE
              LNCHAN(NCHAN)=LNTEST
              READ(IFILTR,*,IOSTAT=IOTEST)SPEC_C,WGT_C
              IF(IOTEST.GT.0)THEN
                  WRITE(IPR,'(/2A,I5,A)')' WARNING: ',                  &
     &              ' Error reading data pair   1 for channel',         &
     &              NCHAN,' in filter function file.'
                  CALL RMCHAN(IFILTR)
                  RETURN
              ELSEIF(IOTEST.LT.0)THEN
                  WRITE(IPR,'(/2A,I5,A)')' WARNING: ',                  &
     &              ' End-Of-File - data pair   1 for channel',         &
     &              NCHAN,' in filter function file.'
                  CALL RMCHAN(IFILTR)
                  RETURN
              ELSEIF(WGT_C.NE.0. .AND. NOT0LO)THEN
                  WRITE(IPR,'(/2A,I5,A,/(11X,2A))')' WARNING: ',        &
     &              ' First filter weight for channel',NCHAN,',',       &
     &              NMCHAN(NCHAN)(1:LNCHAN(NCHAN)),',','is non-zero.'
                  NOT0LO=.FALSE.
              ENDIF
              READ(IFILTR,*,IOSTAT=IOTEST)SPEC_N,WGT_N
              IF(IOTEST.GT.0)THEN
                  WRITE(IPR,'(/2A,I5,A)')' WARNING: ',                  &
     &              ' Error reading data pair   2 for channel',         &
     &              NCHAN,' in filter function file.'
                  CALL RMCHAN(IFILTR)
                  RETURN
              ELSEIF(IOTEST.LT.0)THEN
                  WRITE(IPR,'(/2A,I5,A)')' WARNING: ',                  &
     &              ' End-Of-File - data pair   2 for channel',         &
     &              NCHAN,' in filter function file.'
                  CALL RMCHAN(IFILTR)
                  RETURN
              ENDIF
              NSPEC=2
              LSTBIN=-1
              SPECLO(NCHAN)=SPEC_C
              IF(UNTFLG.EQ.'W')THEN
                  NFRQLO(NCHAN)=INT(BNDWID*NINT(SPEC_C/BNDWID))
              ELSE
                  NFRQHI(NCHAN)=INT(BNDWID*NINT(CONVRT/SPEC_C))
              ENDIF
              WDFREQ(NCHAN)=0.
              WDWAVE(NCHAN)=0.
              DO ILOS=1,NLOS
                  DO IOUT=1,MOUT
                      VLCHAN(IOUT,ILOS,NCHAN)=0.D0
                  ENDDO
              ENDDO

!             INITIALIZE SPECTRAL MOMENTS:
              WGT=(SPEC_N-SPEC_C)*WGT_C
              NORM=WGT
              COEF=SPEC_N+2*SPEC_C
              MOM1ST(NCHAN)=COEF*WGT
              COEF=SPEC_N*COEF+3*SPEC_C**2
              MOM2ND=COEF*WGT
              COEF=SPEC_N*COEF+4*SPEC_C**3
              MOM3RD=COEF*WGT
          ELSE
              WRITE(IPR,'(/A,I5,A,/(10X,A))')                           &
     &          ' WARNING:  Only the first',NCHAN,                      &
     &          ' channels in the filter function',                     &
     &          ' file are being used; increase parameter MXCHAN',      &
     &          ' to accommodate all the channels in the file.'
              RETURN
          ENDIF

!     PROCESS DATA AND THEN RETURN TO 10 TO READ NEXT LINE.
      IF(WTCHAN(LSTBIN,SPEC_C,WGT_C,SPEC_N,WGT_N))GOTO 10

!     PROBLEM WITH CURRENT CHANNEL:
      WRITE(IPR,'(/2(A,I5),/21X,A,I5,A)')'Warning from FILTER: '//      &
     &  ' Routine WTCHAN error for pairs',NSPEC-1,' and',NSPEC,         &
     &  ' of channel',NCHAN,' in filter function file.'
      CALL RMCHAN(IFILTR)
      RETURN
      END
