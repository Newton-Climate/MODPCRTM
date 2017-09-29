      SUBROUTINE RDSUN(SFWHM,DFTSUN,USRSUN,SOLCON,LNFLRT,FLRT)

!     THIS ROUTINE READS 'DFTSUN', AND POSSIBLY, 'USRSUN' FILES.
!     THESE SOLAR IRRADIANCE FILES ARE MERGED.  FINALLY THIS
!     ROUTINE PASSES A TRIANGULAR SLIT WITH SFWHM=|SFWHM|.
!     THE RESULTING DATA IS STORED IN ARRAY SUN OF /SOL01/.
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     INPUT ARGUMENTS:
!       SFWHM    FULL-WIDTH-AT-HALF-MAXIMUM USED IN TRIANGULAR SLIT
!                FUNCTION [IF NEGATIVE, ABS(SFWHM) IS USED IN THE
!                SLIT FUNCTION AND THE DATA IS OUTPUT TO A FILE].
!       DFTSUN   DEFAULT TOP-OF-ATMOSPHERE SOLAR IRRADIANCE FILE.
!       USRSUN   USER-SPECIFIED TOP-OF-ATMOSPHERE SOLAR IRRADIANCE FILE.
!       SOLCON   SOLAR CONSTANT (WATTS/M^2) IF > 0.0; SCALE FACTOR IF <
!       LNFLRT   LENGTH OF I/O FILE ROOT NAME, 0 IF NO mod5root.in FILE.
!       FLRT     ROOT NAME FOR ALL I/O FILES.
      INTEGER LNFLRT
      REAL SOLCON,SFWHM
      CHARACTER DFTSUN*(*),USRSUN*(*),FLRT*(*)

!     COMMONS:
      INCLUDE 'IFIL.h'
      REAL SUN01
      COMMON/SOL01/SUN01(0:MSUN01)

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      SAVE /SOL01/
      EXTERNAL DEVCBD,S01BD

!     FUNCTIONS:
!       NUNIT    RETURNS AN UNUSED FILE UNIT NUMBER.
!       LENSTR   RETURNS STRING LENGTH AFTER TRIMMING LEADING BLANKS.
      INTEGER NUNIT,LENSTR

!     LOCAL VARIABLES:
!       LENFL1   LENGTH OF DFTSUN FILE NAME.
      INTEGER NSUN,IVLO,IVHI,IV,IRES,IRES2,IRESM1,IVOFF,I,J,LENFL1
      LOGICAL LEXIST,LOPEN
      CHARACTER SUNCHR*12,FULLNM*(NAMLEN+8)
      REAL SUNIV,SUNSUM,SSCALE

!     DATA:
!       A        COEFFICIENT FOR LOW FREQUENCY POWER LAW APPROXIMATION.
!       B        EXPONENT FOR LOW FREQUENCY POWER LAW APPROXIMATION.
!       LBLKDT   LOGICAL FLAG, TRUE TO CREATE sXXbd.f FILE.
!       FILNAM   NAME OF BLOCK DATA OUTPUT FILE WHEN LBLKDT IS TRUE.
      REAL A,B
      LOGICAL LBLKDT
      CHARACTER FILNAM*7
      DATA A,B/3.50187E-13,1.93281/,LBLKDT/.FALSE./,FILNAM/'sXXbd.f'/

!     OPEN INPUT DATA FILE.
      INQUIRE(FILE=DFTSUN,EXIST=LEXIST,OPENED=LOPEN,NUMBER=NSUN)
      LENFL1=LENSTR(DFTSUN)
      IF(.NOT.LEXIST)THEN
          WRITE(IPR,'(/3A)')                                            &
     &      ' Error in RDSUN:  File ',DFTSUN(1:LENFL1),' was not found.'

!         IF(LJMASS)CALL WRTBUF(FATAL)
          STOP ' Error in RDSUN:  Solar irradiance file was not found.'
      ENDIF
      IF(LOPEN)THEN
          REWIND(NSUN)
      ELSE
          NSUN=NUNIT()
          CALL OPNFL(NSUN,0,DFTSUN(1:LENFL1),'OLD','FORMATTED','RDSUN')
      ENDIF

!     READ FIRST LINE OF SOLAR DATA
!       IVLO     INITIAL FREQUENCY [CM-1]
!       SUN01    SOLAR IRRADIANCE AT IVLO [W CM-2 / CM-1]
      READ(NSUN,*,ERR=40)
      READ(NSUN,*,ERR=40)
      READ(NSUN,*,ERR=40)IVLO,SUNIV
      IF(IVLO.LT.0 .OR. IVLO.GT.MSUN01)GOTO 40
      IVHI=IVLO
      SUN01(IVHI)=SUNIV

!     READ REMAINING SOLAR DATA (1 CM-1 INCREMENTS):
   10 CONTINUE
      READ(NSUN,*,END=20,ERR=50)IV,SUNIV
      IF(IV.NE.IVHI+1)GOTO 50
      IVHI=IV
      SUN01(IVHI)=SUNIV
      IF(IVHI.LT.MSUN01)GOTO 10
   20 CONTINUE
      CLOSE(NSUN)

!     IF REQUESTED, READ ADDITIONAL "USRSUN" AND MERGE WITH SUN ARRAY:
      IF(USRSUN.NE.'1' .AND. USRSUN.NE.'2' .AND. USRSUN.NE.'3' .AND.    &
     &   USRSUN.NE.'4' .AND. USRSUN.NE.'5' .AND. USRSUN.NE.'6' .AND.    &
     &   USRSUN.NE.'7')CALL RDUSRS(USRSUN,.FALSE.,MSUN01,SUN01(0))

!     SCALE OR RENORMALIZE IF SOLCON (SOLAR CONSTANT) IS NOT ZERO:
      IF(SOLCON.LT.0.)THEN

!         TREAT ABSOLUTE VALUE AS A SCALE FACTOR.  NORMALLY, THE
!         ABSOLUTE VALUE WILL BE CLOSE TO UNITY.
          SSCALE=ABS(SOLCON)
      ELSEIF(SOLCON.GT.0.)THEN

!         RENORMALIZE SOLAR DATA TO INPUT SOLAR CONSTANT [WATTS/M*2].
!         FOR REFERENCE, SOLAR "CONSTANTS" FOR THE 4 DATA FILES ARE:
!           1368.0 W/M**2 FOR NEWKUR.DAT      1359.8 FOR CHKUR.DAT
!           1362.1 W/M**2 FOR CEBCHKUR.DAT    1376.2 FOR THKUR.DAT
          SUNSUM=0.
          DO I=IVLO,IVHI
              SUNSUM=SUN01(I)+SUNSUM
          ENDDO

!         MULTIPLY BY 1.E4 TO CONVERT FROM W/CM^2 TO WATTS/M^2:
          SUNSUM=SUNSUM*1.E4
          SSCALE=SOLCON/SUNSUM
      ELSEIF(SOLCON.EQ.0.)THEN

!         DO NOTHING
          SSCALE=1.
      ENDIF

!     PASS TRIANGULAR SLIT FUNCTION OVER DATA.
      IRES=NINT(ABS(SFWHM))
      IF(IRES.GT.1)THEN
          IRES2=IRES*IRES
          IRESM1=IRES-1
          DO IV=IVLO+IRESM1,IVHI-IRESM1
              SUNSUM=IRES*SUN01(IV)
              DO IVOFF=1,IRESM1
                  SUNSUM=SUNSUM                                         &
     &              +(IRES-IVOFF)*(SUN01(IV+IVOFF)+SUN01(IV-IVOFF))
              ENDDO

!             AVOID OVERWRITING BY STORING DATA OFFSET BY IRESM1 CM-1.
              SUN01(IV-IRESM1)=SUNSUM/IRES2
          ENDDO
      ELSE
          IRESM1=0
      ENDIF

!     MOVE DATA TO PROPER LOCATION (FREQUENCY).
!     ALSO, MULTIPLY BY SCALE HERE.
      DO IV=IVHI-IRESM1,IVLO+IRESM1,-1
          SUN01(IV)=SUN01(IV-IRESM1)*SSCALE
      ENDDO

!     DEFINE LOW FREQUENCY DATA [W CM-2 / CM-1] USING
!     POWER LAW APPROXIMATION.
      SUN01(0)=0.
      DO IV=1,IVLO+IRESM1-1
          SUN01(IV)=A*FLOAT(IV)**B
      ENDDO

!     INSERT A CONSTANT VALUE AT HIGH FREQUENCIES.
      DO IV=IVHI-IRESM1+1,MSUN01
          SUN01(IV)=SUN01(IVHI-IRESM1)
      ENDDO

!     THE FOLLOWING CODING IS USED TO CREATE sXXbd.f BLOCK DATA files.
!     THE LOGICAL FLAG LBLKDT IS HARD-WIRED TO FALSE, BUT THE CODING IS
!     LEFT HERE TO AID IN INCORPORATION OF NEW SOLAR IRRADIANCE DATA.
      IF(LBLKDT)THEN

!         OPEN sXXbd.f FILE.
          WRITE(FILNAM(2:3),'(I2.2)')IRES
          CALL OPNFL(NSUN,0,FILNAM,'UNKNOWN','FORMATTED','RDSUN')

!         WRITE HEADER.
          WRITE(NSUN,'(A,I2.2,A)')'      BLOCK DATA S',IRES,'BD'
          WRITE(NSUN,'(2(/A,I3,A),/(A))')                               &
     &      '!     SOLAR IRRADIANCES [W CM-2 / CM-1]'                   &
     &      //' TABULATED EACH',IRES,' CM-1',                           &
     &      '!     FROM 0 TO 50,000 CM-1 AT',IRES,' CM-1'               &
     &      //' SPECTRAL RESOLUTION (FWHM',                             &
     &      '!     TRIANGULAR SLIT).  DATA DERIVED'                     &
     &      //' FROM THE 1 CM-1 SQUARE SLIT',                           &
     &      '!     TABULATION OF KURUCZ HIGH SPECTRAL'                  &
     &      //' RESOLUTION CALCULATIONS.'
          WRITE(NSUN,'(/(A))')                                          &
     &      '!     PARAMETERS:',                                        &
     &      '      INCLUDE ''PARAMS.h'''
          WRITE(NSUN,'(/A)')'!     COMMONS:'
          IF(IRES.EQ.1)THEN
              WRITE(NSUN,'(/(A))')                                      &
     &          '!     /SUNFLG/',                                       &
     &          '!       LNFLRT   LENGTH OF I/O FILE ROOT'              &
     &          //' NAME, 0 IF NO mod5root.in FILE.',                   &
     &          '!       LRDSUN   SOLAR DATA FLAG, TRUE IF'             &
     &          //' IRRADIANCES IN COMMON BLOCK',                       &
     &          '!                /SOL01/ HAVE BEEN MODIFIED'           &
     &          //' FROM THE BLOCK DATA.',                              &
     &          '      INTEGER LNFLRT',                                 &
     &          '      LOGICAL LRDSUN',                                 &
     &          '      COMMON/SUNFLG/LNFLRT,LRDSUN'
              WRITE(NSUN,'(/(A))')                                      &
     &          '!     /SUNNAM/',                                       &
     &          '!       SUNFIL   NAME OF DEFAULT FILE CONTAINING'      &
     &          //' SOLAR IRRADIANCE DATA.',                            &
     &          '!       FLRT     ROOT NAME FOR ALL I/O FILES.',        &
     &          '      CHARACTER FLRT*(NAMLEN-4),SUNFIL*(LENSUN)',      &
     &          '      COMMON/SUNNAM/FLRT,SUNFIL'
          ENDIF
          WRITE(NSUN,'(/(A,I2.2,A))')                                   &
     &      '!     /SOL',IRES,'/',                                      &
     &      '!       SUN',IRES,'    TOP-OF-ATMOSPHERE SOLAR'            &
     &      //' IRRADIANCES [W CM-2 / CM-1].',                          &
     &      '      REAL SUN',IRES
          WRITE(NSUN,'(3(A,I2.2),A)')                                   &
     &      '      COMMON/SOL',IRES,'/SUN',IRES,'(0:MSUN',IRES,')'
          IF(IRES.EQ.1)THEN
              WRITE(NSUN,'(/A)')'      SAVE /SUNFLG/,/SUNNAM/,/SOL01/'
          ELSE
              WRITE(NSUN,'(/A,I2.2,A)')'      SAVE /SOL',IRES,'/'
          ENDIF
          WRITE(NSUN,'(/A,/(A,I2.2,A))')                                &
     &      '!     LOCAL VARIABLES:',                                   &
     &      '!       WAV',IRES,'    FREQUENCY INDEX [CM-1].',           &
     &      '      INTEGER WAV',IRES
          IF(IRES.EQ.1)WRITE(NSUN,'(/(A))')                             &
     &      '!     DATA:',                                              &
     &      '      DATA LRDSUN/.FALSE./,SUNFIL/''SUN01kurucz1997.dat''/'

!         WRITE 0 TO 50 CM-1 DATA.
          WRITE(NSUN,'(/A,I6,A,/A,3(I2.2,A),26X,A,                      &
     &      /5X,A,16X,1P,3(E15.6,A),A,/(5X,A,4(E15.6,A),A))')           &
     &      '!         0 TO',50*IRES,' CM-1 DATA:','      DATA (SUN',   &
     &      IRES,'(WAV',IRES,'),WAV',IRES,'=     0,    50)/','&',       &
     &       '&',(DBLE(SUN01(I*IRES)),',',I=  0, 2),'  &',              &
     &      ('&',(DBLE(SUN01(I*IRES)),',',I=J-3, J),'  &',J=6,46,4),    &
     &       '&',(DBLE(SUN01(I*IRES)),',',I= 47,49),                    &
     &      DBLE(SUN01(50*IRES)),'/'

!         WRITE 51 TO 50*INT(MSUN01/50) CM-1 DATA.
          DO IV=100,MSUN01/IRES,50
              WRITE(NSUN,'(/A,2(I6,A),/A,3(I2.2,A),2(I6,A),26X,A,       &
     &          /5X,A,32X,1P,2(E15.6,A),A,/(5X,A,4(E15.6,A),A))')       &
     &          '!    ',(IV-49)*IRES,' TO',IV*IRES,' CM-1 DATA:',       &
     &          '      DATA (SUN',IRES,'(WAV',IRES,'),WAV',IRES,'=',    &
     &          IV-49,',',IV,')/','&','&',(DBLE(SUN01(I*IRES)),',',     &
     &          I=IV-49,IV-48),'  &',('&',(DBLE(SUN01(I*IRES)),         &
     &          ',',I=J-3,J),'  &',J=IV-44,IV-4,4),'&',(DBLE(SUN01      &
     &          (I*IRES)),',',I=IV-3,IV-1),DBLE(SUN01(IV*IRES)),'/'
          ENDDO

!         WRITE REMAINING CM-1 DATA.
          IF(J.LT.MSUN01/IRES)THEN
              WRITE(NSUN,'(/A,2(I6,A),/A,3(I2.2,A),2(I6,A),26X,A)')     &
     &          '!    ',(J+1)*IRES,' TO',(MSUN01/IRES)*IRES,            &
     &          ' CM-1 DATA:','      DATA (SUN',IRES,'(WAV',IRES,       &
     &          '),WAV',IRES,'=',J+1,',',MSUN01/IRES,')/','&'
   30         CONTINUE
              IF(J+4.LT.MSUN01/IRES)THEN
                  WRITE(NSUN,'(5X,A,1P,4(E15.6,A),A)')                  &
     &              '&',(DBLE(SUN01((J+I)*IRES)),',',I=1,4),'  &'
                  J=J+4
                  GOTO 30
              ENDIF
              IF(J+1.EQ.MSUN01/IRES)THEN
                  WRITE(NSUN,'(5X,A,1P,E15.6,A)')                       &
     &              '&',DBLE(SUN01((MSUN01/IRES)*IRES)),'/'
              ELSE
                  WRITE(NSUN,'(5X,1P,5(A,E12.6))')                      &
     &              '&',(DBLE(SUN01(I*IRES)),',',I=J+1,MSUN01-1),       &
     &              DBLE(SUN01((MSUN01/IRES)*IRES)),'/'
              ENDIF
          ENDIF
          WRITE(NSUN,'(A)')'      END'
          CLOSE(NSUN)
      ENDIF

!     RETURN UNLESS SFWHM IS NEGATIVE
      IF(SFWHM.GE.0.)RETURN

!     OPEN OUTPUT DATA FILE.
      IF(USRSUN.EQ.'1' .OR. USRSUN.EQ.'2' .OR. USRSUN.EQ.'3' .OR.       &
     &   USRSUN.EQ.'4' .OR. USRSUN.EQ.'5' .OR. USRSUN.EQ.'6' .OR.       &
     &   USRSUN.EQ.'7')THEN
          SUNCHR='_sn_0000.dat'
          SUNCHR(3:3)=USRSUN
      ELSE
          SUNCHR='_s0_0000.dat'
      ENDIF
      WRITE(SUNCHR(5:8),'(I4.4)')IRES
      IF(LNFLRT.LE.0)THEN
          CALL OPNFL(NSUN,0,SUNCHR,'UNKNOWN','FORMATTED','RDSUN')
      ELSE

!         FULLNM MUST BE DEFINED BECAUSE OF GNU COMPILER PROBLEMS:
          FULLNM(1:LNFLRT+12)=FLRT(1:LNFLRT)//SUNCHR
          CALL OPNFL(NSUN,0,FULLNM(1:LNFLRT+12),                        &
     &      'UNKNOWN','FORMATTED','RDSUN')
      ENDIF

!     WRITE DATA TO SUNCHR IN ORIGINAL UNITS [W CM-2 / CM-1].
      WRITE(NSUN,'(A,/A,/(I7,1P,E13.3))')'  FREQ   SOLAR IRRADIANCE',   &
     &  ' (CM-1)  (W CM-2 / CM-1)',(IV,SUN01(IV),IV=1,MSUN01)
      CLOSE(NSUN)

!     RETURN TO DRIVER.
      RETURN

!     WRITE OUT ERROR MESSAGE AND STOP.
   40 CONTINUE
      WRITE(IPR,'(/3A)')                                                &
     &  ' Error in RDSUN:  Problem reading/using file ',DFTSUN,' data.'
!     IF(LJMASS)CALL WRTBUF(FATAL)
      STOP 'Error in RDSUN:  Problem reading solar irradiance data'

!     WRITE OUT ERROR MESSAGE AND STOP.
   50 CONTINUE
      WRITE(IPR,'(/2A,/17X,A,I6,A)')                                    &
     &  ' Error in RDSUN:  Problem reading/using file ',DFTSUN,         &
     &  ' data.  Last successfully frequency was',IVHI,' cm-1'
!     IF(LJMASS)CALL WRTBUF(FATAL)
      STOP 'Error in RDSUN:  Problem reading solar irradiance data'
      END

      SUBROUTINE RDUSRS(USRSUN,L_P1,MSUNDT,SUNDAT)

!     RDUSRS READS USER-SPECIFIED SOLAR IRRADIANCE FILE AND INTERPOLATES
!     DATA ONTO EITHER A 1.0 CM-1 OR 0.1 CM-1 SPECTRAL FREQUENCY GRID.

!     FILE USRSUN CONTAINS HEADER INFORMATION ON THE FIRST LINE FOLLOWED
!     BY ONE DATA PAIR PER LINE.  HEADER HAS TWO INTEGERS, EACH OF WHICH
!     IS 1, 2 OR 3.  THE FIRST DESIGNATES THE SPECTRAL GRID UNIT:
!       1 FOR SPECTRAL FREQUENCY IN CM-1
!       2 FOR SPECTRAL WAVELENTH IN NM
!       3 FOR SPECTRAL WAVELENTH IN MICRONS

!     THE SECOND INTEGER ON THE FIRST LINE DESIGNATES THE SOLAR
!     IRRADIANCE UNITS.
!       1 FOR SOLAR IRRADIANCE IN W CM-2 / CM-1
!       2 FOR SOLAR IRRADIANCE IN PHOTONS SEC-1 CM-2 / NM
!       3 FOR SOLAR IRRADIANCE IN mW M-2 / NM = W M-2 / MICRON

!     E.G., FOR FREQUENCY IN CM-1 AND IRRADIANCE IN W CM-2 / CM-1,
!     THE USER-CHOSEN FILE SHOULD LOOK LIKE THIS:
!       1   1
!       51    7.453E-10
!       52    7.711E-10
!       53    7.974E-10
!       54    8.243E-10
!       55    8.516E-10
!       ...
!       49982    2.800E-09
!       49983    2.603E-09
      IMPLICIT NONE

!     PARAMETERS:
!       HC    PLANCK CONSTANT TIMES THE SPEED OF LIGHT [W SEC NM].
      REAL HC
      PARAMETER(HC=1.9864455E-16)

!     ARGUMENTS:
!       USRSUN   USER-SPECIFIED TOP-OF-ATMOSPHERE SOLAR IRRADIANCE FILE.
!       L_P1     LOGICAL, .TRUE. FOR 0.1 CM-1 SOLAR IRRADIANCE DATA.
!       MSUNDT   MAXIMUM INDEX OF SUNDAT ARRAY.
!       SUNDAT   TOP-OF-ATMOSPHERE SOLAR IRRADIANCES [W CM-2 / CM-1].
      CHARACTER USRSUN*(*)
      LOGICAL L_P1
      INTEGER MSUNDT
      REAL SUNDAT(0:MSUNDT)

!     COMMONS:
      INCLUDE 'IFIL.h'

!     FUNCTIONS:
!       NUNIT    RETURNS AN UNUSED FILE UNIT NUMBER.
!       LENSTR   RETURNS STRING LENGTH AFTER TRIMMING LEADING BLANKS.
      INTEGER NUNIT,LENSTR

!     LOCAL VARIABLES:
!       IUSR     FILE UNIT NUMBER FOR THE USRSUN FILE.
!       IOS      IOSTAT RESULT.
!       JSPEC    SPECTRAL GRID UNIT INDEX (1=CM-1, 2=NM, 3=MICRONS).
!       JIRRAD   SOLAR IRRADIANCE UNIT INDEX.
!       IFRQHI   UPPER INDEX FOR SPECTRAL FREQUENCY INTERPOLATION LOOP.
!       IFRQLO   LOWER INDEX FOR SPECTRAL FREQUENCY INTERPOLATION LOOP.
!       IFRQ     SPECTRAL INTERPOLATION INDEX.
!       SPECHI   CURRENT SPECTRAL GRID VALUE.
!       SUNHI    SOLAR IRRADIANCE AT SPECHI.
!       SPECLO   PREVIOUS SPECTRAL GRID VALUE.
!       SUNLO    SOLAR IRRADIANCE AT SPECLO.
!       SPCHI    SPECTRAL INTERPOLATION NUMERATOR TERM.
!       SPCDEN   SPECTRAL INTERPOLATION DENOMINATOR.
!       WAVMIN   MINIMUM SPECTRAL WAVELENGTH OF SUNDAT ARRAY [um].
      INTEGER IUSR,IOS,JSPEC,JIRRAD,IFRQHI,IFRQLO,IFRQ
      REAL SPECHI,SUNHI,SPECLO,SUNLO,SPCHI,SPCDEN,WAVMIN

!     OPEN USER SOLAR IRRADIANCE DATA FILE:
      IUSR=NUNIT()
      CALL OPNFL(IUSR,0,USRSUN(1:LENSTR(USRSUN)),                       &
     &  'OLD','FORMATTED','RDUSRS')

!     READ AND CHECK HEADER INFORMATION:
      READ(IUSR,*,IOSTAT=IOS)JSPEC,JIRRAD
      IF(IOS.NE.0 .OR. ABS(JSPEC-2).GT.1 .OR. ABS(JIRRAD-2).GT.1)THEN
          WRITE(IPR,'(/2A,/22X,A)')' Warning from RDUSRS:  Unable to '//&
     &      'read solar irradiance data file ',USRSUN(1:LENSTR(USRSUN)),&
     &      ' Default solar irradiance data will be used!'
      ELSEIF(JSPEC.EQ.1)THEN

!         SPECTRAL UNIT IS CM-1: READ AND CHECK FIRST PAIR OF DATA.
          READ(IUSR,*,IOSTAT=IOS)SPECHI,SUNHI
          IF(IOS.NE.0)THEN
              WRITE(IPR,'(/2A,/22X,A)')' Warning from RDUSRS: '         &
     &          //' Unable to read solar irradiance data file ',        &
     &          USRSUN(1:LENSTR(USRSUN)),                               &
     &          ' Default solar irradiance data will be used!'
              CLOSE(IUSR,STATUS='KEEP')
              RETURN
          ENDIF
          IF(L_P1)THEN
              IFRQHI=INT(10*SPECHI)
          ELSE
              IFRQHI=INT(SPECHI)
          ENDIF
          IF(IFRQHI.GE.MSUNDT .OR. SPECHI.LT.0.)THEN
              WRITE(IPR,'(/2A,/22X,A)')' Warning from RDUSRS: '         &
     &          //' Spectral incompatibility with user solar',          &
     &          ' irradiance data file ',USRSUN(1:LENSTR(USRSUN)),      &
     &          ' Default solar irradiance data will be used!'
              CLOSE(IUSR,STATUS='KEEP')
              RETURN
          ENDIF

!         ENTRY POINT FOR NEXT PAIR OF DATA:
   10     CONTINUE

!         READ AND CHECK SOLAR IRRADIANCE DATA:
          SPECLO=SPECHI
          SUNLO=SUNHI
          READ(IUSR,*,IOSTAT=IOS)SPECHI,SUNHI
          IF(IOS.LT.0)THEN

!             SUCCESSFULLY REACHED END-OF-FILE.
              CLOSE(IUSR,STATUS='KEEP')
              RETURN
          ELSEIF(IOS.GT.0)THEN
              WRITE(IPR,'(/2A,/22X,A,F10.2,A)')' Warning from RDUSRS: ' &
     &          //' Error reading solar irradiance data file ',         &
     &          USRSUN(1:LENSTR(USRSUN)),                               &
     &          ' Only processed data up to',SPECLO,' cm-1.'
              CLOSE(IUSR,STATUS='KEEP')
              RETURN
          ELSEIF(SPECHI.LE.SPECLO)THEN
              WRITE(IPR,'(/2A,/22X,A,F10.2,A)')' Warning from RDUSRS: ' &
     &          //' Spectral grid not monotonically increasing in',     &
     &          ' solar irradiance data file ',USRSUN(1:LENSTR(USRSUN)),&
     &          ' Only processed data up to',SPECLO,' cm-1.'
              CLOSE(IUSR,STATUS='KEEP')
              RETURN
          ENDIF

!         INTERPOLATE IN FREQUENCY AND THEN CONVERT IRRADIANCE UNIT:
          SPCHI=SPECHI
          SPCDEN=SPECHI-SPECLO
          IFRQLO=IFRQHI+1
          IF(L_P1)THEN
              SPCHI=10*SPCHI
              SPCDEN=10*SPCDEN
              IFRQHI=MIN(INT(SPCHI),MSUNDT)
              IF(JIRRAD.EQ.2)THEN

!                 INPUT SOLAR IRRADIANCE IN PHOTONS SEC-1 CM-2 / NM:
                  DO IFRQ=IFRQLO,IFRQHI
                      SUNDAT(IFRQ)=10*HC*(SUNHI+(SPCHI-IFRQ)            &
     &                  *(SUNLO-SUNHI)/SPCDEN)/IFRQ
                  ENDDO
              ELSEIF(JIRRAD.EQ.3)THEN

!                 INPUT SOLAR IRRADIANCE IN W M-2 / MICRON
                  DO IFRQ=IFRQLO,IFRQHI
                      SUNDAT(IFRQ)=100*(SUNHI+(SPCHI-IFRQ)              &
     &                  *(SUNLO-SUNHI)/SPCDEN)/IFRQ**2
                  ENDDO
              ELSE
                  DO IFRQ=IFRQLO,IFRQHI
                      SUNDAT(IFRQ)=SUNHI+(SPCHI-IFRQ)                   &
     &                  *(SUNLO-SUNHI)/SPCDEN
                  ENDDO
              ENDIF
          ELSE
              IFRQHI=MIN(INT(SPCHI),MSUNDT)
              IF(JIRRAD.EQ.2)THEN

!                 INPUT SOLAR IRRADIANCE IN PHOTONS SEC-1 CM-2 / NM:
                  DO IFRQ=IFRQLO,IFRQHI
                      SUNDAT(IFRQ)=HC*(SUNHI+(SPCHI-IFRQ)               &
     &                  *(SUNLO-SUNHI)/SPCDEN)/IFRQ
                  ENDDO
              ELSEIF(JIRRAD.EQ.3)THEN

!                 INPUT SOLAR IRRADIANCE IN W M-2 / MICRON
                  DO IFRQ=IFRQLO,IFRQHI
                      SUNDAT(IFRQ)=(SUNHI+(SPCHI-IFRQ)                  &
     &                  *(SUNLO-SUNHI)/SPCDEN)/IFRQ**2
                  ENDDO
              ELSE
                  DO IFRQ=IFRQLO,IFRQHI
                      SUNDAT(IFRQ)=SUNHI+(SPCHI-IFRQ)                   &
     &                  *(SUNLO-SUNHI)/SPCDEN
                  ENDDO
              ENDIF
          ENDIF

!         MORE DATA?
          IF(IFRQHI.LT.MSUNDT)GOTO 10
      ELSE

!         SPECTRAL UNIT IS NM OR uM: READ AND CHECK FIRST PAIR OF DATA.
          READ(IUSR,*,IOSTAT=IOS)SPECHI,SUNHI
          IF(IOS.NE.0)THEN
              WRITE(IPR,'(/2A,/22X,A)')' Warning from RDUSRS: '         &
     &          //' Unable to read solar irradiance data file ',        &
     &          USRSUN(1:LENSTR(USRSUN)),                               &
     &          ' Default solar irradiance data will be used!'
              CLOSE(IUSR,STATUS='KEEP')
              RETURN
          ELSEIF(SPECHI.LE.0.)THEN
              WRITE(IPR,'(/2A,/22X,A)')' Warning from RDUSRS: '         &
     &          //' Non-positive spectral wavelength in user solar',    &
     &          ' irradiance data file ',USRSUN(1:LENSTR(USRSUN)),      &
     &          ' Default solar irradiance data will be used!'
              CLOSE(IUSR,STATUS='KEEP')
              RETURN
          ENDIF

!         IF NM, CONVERT TO MICRONS:
          IF(JSPEC.EQ.2)SPECHI=SPECHI/1000

!         DETERMINE MINIMUM WAVELENGTH [um] AND FREQUENCY LOOP MINIMUM:
          IF(L_P1)THEN
              WAVMIN=100000./MSUNDT
              IFRQLO=INT(100000/SPECHI+1)
          ELSE
              WAVMIN=10000./MSUNDT
              IFRQLO=INT(10000/SPECHI+1)
          ENDIF

!         ENTRY POINT FOR NEXT PAIR OF DATA:
   20     CONTINUE

!         READ AND CHECK SOLAR IRRADIANCE DATA:
          SPECLO=SPECHI
          SUNLO=SUNHI
          READ(IUSR,*,IOSTAT=IOS)SPECHI,SUNHI
          IF(IOS.LT.0)THEN

!             SUCCESSFULLY REACHED END-OF-FILE.
              CLOSE(IUSR,STATUS='KEEP')
              RETURN
          ELSEIF(IOS.GT.0)THEN
              WRITE(IPR,'(/2A,/22X,A,F10.2,A)')' Warning from RDUSRS: ' &
     &          //' Error reading solar irradiance data file ',         &
     &          USRSUN(1:LENSTR(USRSUN)),                               &
     &          ' Only processed data up to',SPECLO,' microns.'
              CLOSE(IUSR,STATUS='KEEP')
              RETURN
          ENDIF
          IF(JSPEC.EQ.2)SPECHI=SPECHI/1000
          IF(SPECHI.LE.SPECLO)THEN
              WRITE(IPR,'(/2A,/22X,A,F10.2,A)')' Warning from RDUSRS: ' &
     &          //' Spectral grid not monotonically increasing in',     &
     &          ' solar irradiance data file ',USRSUN(1:LENSTR(USRSUN)),&
     &          ' Only processed data up to',SPECLO,' microns.'
              CLOSE(IUSR,STATUS='KEEP')
              RETURN
          ELSEIF(SPECHI.LE.WAVMIN)THEN

!             CURRENT DATA IS OUT OF THE SPECTRAL RANGE OF SUNDAT:
              GOTO 20
          ENDIF

!         INTERPOLATE IN WAVELENGTH [um] & THEN CONVERT IRRADIANCE UNIT:
          SPCDEN=SPECHI-SPECLO
          IFRQHI=IFRQLO-1
          IF(L_P1)THEN
              IFRQLO=INT(100000/SPECHI+1)
              IF(JIRRAD.EQ.2)THEN

!                 INPUT SOLAR IRRADIANCE IN PHOTONS SEC-1 CM-2 / NM:
                  DO IFRQ=IFRQLO,IFRQHI
                      SUNDAT(IFRQ)=10*HC*(SUNHI+(SPECHI-1.E5/IFRQ)      &
     &                  *(SUNLO-SUNHI)/SPCDEN)/IFRQ
                  ENDDO
              ELSEIF(JIRRAD.EQ.3)THEN

!                 INPUT SOLAR IRRADIANCE IN W M-2 / MICRON
                  DO IFRQ=IFRQLO,IFRQHI
                      SUNDAT(IFRQ)=100*(SUNHI+(SPECHI-1.E5/IFRQ)        &
     &                  *(SUNLO-SUNHI)/SPCDEN)/FLOAT(IFRQ)**2
                  ENDDO
              ELSE
                  DO IFRQ=IFRQLO,IFRQHI
                      SUNDAT(IFRQ)=SUNHI+(SPECHI-1.E5/IFRQ)             &
     &                  *(SUNLO-SUNHI)/SPCDEN
                  ENDDO
              ENDIF
          ELSE
              IFRQLO=INT(10000/SPECHI+1)
              IF(JIRRAD.EQ.2)THEN

!                 INPUT SOLAR IRRADIANCE IN PHOTONS SEC-1 CM-2 / NM:
                  DO IFRQ=IFRQLO,IFRQHI
                      SUNDAT(IFRQ)=HC*(SUNHI+(SPECHI-1.E4/IFRQ)         &
     &                  *(SUNLO-SUNHI)/SPCDEN)/IFRQ
                  ENDDO
              ELSEIF(JIRRAD.EQ.3)THEN

!                 INPUT SOLAR IRRADIANCE IN W M-2 / MICRON
                  DO IFRQ=IFRQLO,IFRQHI
                      SUNDAT(IFRQ)=(SUNHI+(SPECHI-1.E4/IFRQ)            &
     &                  *(SUNLO-SUNHI)/SPCDEN)/FLOAT(IFRQ)**2
                  ENDDO
              ELSE
                  DO IFRQ=IFRQLO,IFRQHI
                      SUNDAT(IFRQ)=SUNHI+(SPECHI-1.E4/IFRQ)             &
     &                  *(SUNLO-SUNHI)/SPCDEN
                  ENDDO
              ENDIF
          ENDIF

!         MORE DATA?
          IF(IFRQLO.GT.1)GOTO 20
      ENDIF
      CLOSE(IUSR,STATUS='KEEP')
      RETURN
      END

      SUBROUTINE HUNT(XA,N,X,JLO,JHI)

!     MODIFIED FROM NUMERICAL RECIPES' ROUTINE HUNT.
!     LIKE LOCATE (FROM NR) BUT FASTER IF JLO IS A GOOD GUESS.
      IMPLICIT NONE

!     ARGUMENTS:
!       XA(1), XA(2), XA(3), ..., XA(N) ARE THE MONOTONIC GRID PTS.
!       RETURNS JLO=J, IF X IS IN BETWEEN XA(J) & XA(J+1).
!       IF X=XA(J), RETURNS JLO=J-1.  IF OUT OF RANGE, RETURNS 0 OR N.
      INTEGER N,JLO,JHI
      REAL XA(N),X

!     LOCAL VARIABLES:
      INTEGER INC,JM
      LOGICAL ASCEND

      IF(X.EQ.XA(1))THEN
          JLO=1
          JHI=2
          RETURN
      ELSEIF(X.EQ.XA(N))THEN
          JLO=N-1
          JHI=N
          RETURN
      ENDIF

      ASCEND=XA(N).GT.XA(1)
      IF(JLO.LE.0 .OR. JLO.GT.N)THEN
          JLO=0
          JHI=N+1
      ELSE
          INC=1
          IF(X.GE.XA(JLO) .EQV. ASCEND)THEN
   10         CONTINUE
              JHI=JLO+INC
              IF(JHI.GT.N)THEN
                  JHI=N+1
              ELSEIF(X.GE.XA(JHI) .EQV. ASCEND)THEN
                  JLO=JHI
                  INC=INC+INC
                  GOTO 10
              ENDIF
          ELSE
              JHI=JLO
   20         CONTINUE
              JLO=JHI-INC
              IF(JLO.LT.1)THEN
                  JLO=0
              ELSEIF(X.LT.XA(JLO) .EQV. ASCEND)THEN
                  JHI=JLO
                  INC=INC+INC
                  GOTO 20
              ENDIF
          ENDIF
      ENDIF

   30 CONTINUE
      IF(JHI-JLO.EQ.1)RETURN

      JM=(JHI+JLO)/2
      IF(X.GT.XA(JM) .EQV. ASCEND)THEN
          JLO=JM
      ELSE
          JHI=JM
      ENDIF
      GOTO 30
      END

      SUBROUTINE DPHUNT(XA,N,X,JLO,JHI)

!     MODIFIED FROM NUMERICAL RECIPES' ROUTINE HUNT.
!     LIKE LOCATE (FROM NR) BUT FASTER IF JLO IS A GOOD GUESS.
      IMPLICIT NONE

!     ARGUMENTS:
!       XA(1), XA(2), XA(3), ..., XA(N) ARE THE MONOTONIC GRID PTS.
!       RETURNS JLO=J, IF X IS IN BETWEEN XA(J) & XA(J+1).
!       IF X=XA(J), RETURNS JLO=J-1.  IF OUT OF RANGE, RETURNS 0 OR N.
      INTEGER N,JLO,JHI
      DOUBLE PRECISION XA(N),X

!     LOCAL VARIABLES:
      INTEGER INC,JM
      LOGICAL ASCEND

      IF(X.EQ.XA(1))THEN
          JLO=1
          JHI=2
          RETURN
      ELSEIF(X.EQ.XA(N))THEN
          JLO=N-1
          JHI=N
          RETURN
      ENDIF

      ASCEND=XA(N).GT.XA(1)
      IF(JLO.LE.0 .OR. JLO.GT.N)THEN
          JLO=0
          JHI=N+1
      ELSE
          INC=1
          IF(X.GE.XA(JLO) .EQV. ASCEND)THEN
   10         CONTINUE
              JHI=JLO+INC
              IF(JHI.GT.N)THEN
                  JHI=N+1
              ELSEIF(X.GE.XA(JHI) .EQV. ASCEND)THEN
                  JLO=JHI
                  INC=INC+INC
                  GOTO 10
              ENDIF
          ELSE
              JHI=JLO
   20         CONTINUE
              JLO=JHI-INC
              IF(JLO.LT.1)THEN
                  JLO=0
              ELSEIF(X.LT.XA(JLO) .EQV. ASCEND)THEN
                  JHI=JLO
                  INC=INC+INC
                  GOTO 20
              ENDIF
          ENDIF
      ENDIF
   30 CONTINUE
      IF(JHI-JLO.EQ.1)RETURN
      JM=(JHI+JLO)/2
      IF(X.GT.XA(JM) .EQV. ASCEND)THEN
          JLO=JM
      ELSE
          JHI=JM
      ENDIF
      GOTO 30
      END

      SUBROUTINE RSUNP1(SFWHM,DFTSUN,USRSUN,SOLCON,LNFLRT,FLRT)

!     THIS IS A MODIFIED COPY OF 'RDSUN' FOR READING 0.1 1/CM SOLAR
!     IRRADIANCE BINARY FILE, WHICH IS A BLUE-FORCED PATCH AND SHOULD
!     BE REVISED. A NEW BLOCK OF DATA '/SOLP1/SUNP1(0:MSUNP1)'
!     WAS CREATED. RECOMMENDATION:
!       (1) TWO BINARY FILES, 0.1 AND 1 1/CM.
!       (2) THESE FILES HAVE 500,000 AND 50,000 4BYTE REAL DATA,
!           RESPECTIVELY.  THEY TAKE 0.025 AND 0.010 SECOND TO READ
!           USING A 1 GH PC, RESPECTIVELY.
!       (3) THEIR RANGE IS FROM 0.05 TO 49999.95 AND 0.5 TO
!           49999.5 1/CM, RESPECTIVELY.
!       (4) DEGRADE THEM RESPECTIVELY TO OBTAIN USER'S RESOLUTION,
!           0.2, 5, ETC.
!       (5) RESAMPLE 'SUNP1' TO USER'S GRID AND RANGE.
!       (6) REVISE 'SOURCE', 'S01BD', 'TRANS', 'DRIVER',
!           'CD1A', 'RSUNP1'.

!     by JL/SSi_4/22/03

!     THIS ROUTINE READS 'DFTSUN', AND POSSIBLY, 'USRSUN' FILES.
!     THESE SOLAR IRRADIANCE FILES ARE MERGED.  FINALLY THIS
!     ROUTINE PASSES A TRIANGULAR SLIT WITH SFWHM=|SFWHM|.
!     THE RESULTING DATA IS STORED IN ARRAY SUN OF /SOLP1/.
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     ARGUMENTS:
!       SFWHM    FULL-WIDTH-AT-HALF-MAXIMUM USED IN TRIANGULAR SLIT
!                FUNCTION [IF NEGATIVE, ABS(SFWHM) IS USED IN THE
!                SLIT FUNCTION AND THE DATA IS OUTPUT TO A FILE].
!       DFTSUN   DEFAULT TOP-OF-ATMOSPHERE SOLAR IRRADIANCE FILE.
!       USRSUN   USER-SPECIFIED TOP-OF-ATMOSPHERE SOLAR IRRADIANCE FILE.
!       SOLCON   SOLAR CONSTANT (WATTS/M^2) IF > 0.0; SCALE FACTOR IF <
!       FLRT     ROOT NAME FOR ALL I/O FILES.
!       LNFLRT   LENGTH OF I/O FILE ROOT NAME, 0 IF NO mod5root.in FILE.
      INTEGER LNFLRT
      REAL SOLCON,SSCALE,SFWHM
      CHARACTER DFTSUN*(*),USRSUN*(*),FLRT*(*)

!     COMMONS:
      INCLUDE 'IFIL.h'
      REAL SUNP1
      COMMON/SOLP1/SUNP1(0:MSUNP1)

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      SAVE /SOLP1/
      EXTERNAL DEVCBD,S01BD

!     FUNCTIONS:
!       NUNIT    RETURNS AN UNUSED FILE UNIT NUMBER.
!       LENSTR   RETURNS STRING LENGTH AFTER TRIMMING LEADING BLANKS.
      INTEGER NUNIT,LENSTR

!     LOCAL VARIABLES:
!       LENFL1   LENGTH OF DFTSUN FILE NAME.
      INTEGER NSUN,IVLO,IVHI,IV,IRES,IRES2,IRESM1,IVOFF,I,J,IVOLD,LINE, &
     &  LENFL1
      LOGICAL LEXIST,LOPEN
      CHARACTER SUNCHR*12,FULLNM*(NAMLEN+8)
      REAL SUNIV,SUNSUM

!     DATA:
!       A        COEFFICIENT FOR LOW FREQUENCY POWER LAW APPROXIMATION.
!       B        EXPONENT FOR LOW FREQUENCY POWER LAW APPROXIMATION.
!       LBLKDT   LOGICAL FLAG, TRUE TO CREATE S01BD.F FILE.
      REAL A,B
      LOGICAL LBLKDT
      DATA A,B/3.50187E-13,1.93281/,LBLKDT/.FALSE./

!     OPEN INPUT DATA FILE:
      INQUIRE(FILE=DFTSUN,EXIST=LEXIST,OPENED=LOPEN,NUMBER=NSUN)
      LENFL1=LENSTR(DFTSUN)
      IF(.NOT.LEXIST)THEN
          WRITE(IPR,'(/3A)')' Error in RSUNP1:  File ',                 &
     &      DFTSUN(1:LENFL1),' was not found.'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP ' Error in RSUNP1:  Solar irradiance file was not found.'
      ENDIF
      IF(LOPEN)THEN
          REWIND(NSUN)
      ELSE
          NSUN=NUNIT()

!jl+      CHECK IF THE FILE IS A BINARY OR ASCII FILE:
          IF(INDEX(DFTSUN,'.B').GT.0 .OR. INDEX(DFTSUN,'.b').GT.0)THEN
              CALL OPNFL(NSUN,0,DFTSUN(1:LENFL1),                       &
     &          'OLD','UNFORMATTED','RDUSRS')
          ELSE
              CALL OPNFL(NSUN,0,DFTSUN(1:LENFL1),                       &
     &          'OLD','FORMATTED','RDUSRS')
          ENDIF
      ENDIF

!     READ 0.1 1/CM BINARY FILE. OR 0.1 1/CM BINARY FILE
!       SUNP1      SOLAR IRRADIANCE [W CM-2 / CM-1] FROM 0 TO 500000 1/CM
      IF(INDEX(DFTSUN,'.B').GT.0 .OR. INDEX(DFTSUN,'.b').GT.0)THEN
          DO I=0,MSUNP1
              READ(NSUN,END=10,ERR=60)SUNP1(I)
          ENDDO
   10     CONTINUE
          IVHI=I-1
          IVLO=0
      ELSE
!         READ FIRST LINE OF SOLAR DATA
!           IVLO     INITIAL FREQUENCY [CM-1]
!           SUN      SOLAR IRRADIANCE AT IVLO [W CM-2 / CM-1]
          READ(NSUN,*,ERR=50)
          READ(NSUN,*,ERR=50)
          READ(NSUN,*,ERR=50)IVLO,SUNIV
          IF(IVLO.LT.0 .OR. IVLO.GT.MSUNP1)GOTO 50
          IVHI=IVLO
          SUNP1(IVHI)=SUNIV

!         READ REMAINING SOLAR DATA (1 CM-1 INCREMENTS):
   20     CONTINUE
          READ(NSUN,*,END=30,ERR=50)IV,SUNIV
          IF(IV.NE.IVHI+1)GOTO 50
          IVHI=IV
          SUNP1(IVHI)=SUNIV
          IF(IVHI.LT.MSUNP1)GOTO 20
      ENDIF
   30 CONTINUE
      CLOSE(NSUN)

!     IF REQUESTED, READ ADDITIONAL "USRSUN" AND MERGE WITH SUN ARRAY:
      IF(USRSUN.NE.'1' .AND. USRSUN.NE.'5' .AND. USRSUN.NE.'7')         &
     &  CALL RDUSRS(USRSUN,.TRUE.,MSUNP1,SUNP1(0))

!     SCALE OR RENORMALIZE IF SOLCON (SOLAR CONSTANT) IS NOT ZERO
      IF(SOLCON.LT.0.)THEN

!         ABSOLUTE VALUE IS SCALE FACTOR, NORMALLY NEAR UNITY.
          SSCALE=ABS(SOLCON)
      ELSEIF(SOLCON.GT.0.)THEN

!         RENORMALIZE SOLAR DATA TO INPUT SOLAR CONSTANT [WATTS/M*2].
!         FOR REFERENCE, SOLAR "CONSTANTS" FOR THE 4 DATA FILES ARE:
!           1368.0 W/M**2 FOR NEWKUR.DAT      1359.8 FOR CHKUR.DAT
!           1362.1 W/M**2 FOR CEBCHKUR.DAT    1376.2 FOR THKUR.DAT
          SUNSUM=0.
          DO I=IVLO,IVHI
              SUNSUM=SUNP1(I)+SUNSUM
          ENDDO

!         MULTIPLY BY 1.0E4 TO CONVERT FROM W/CM^2 TO WATTS/M^2
          SUNSUM=SUNSUM*1.E4
          SSCALE=SOLCON/SUNSUM
      ELSE

!         DO NOTHING
          SSCALE=1.
      ENDIF

!     PASS TRIANGULAR SLIT FUNCTION OVER DATA.
!      MAKE 0.1 1/CM AS DEFAULT GRID. RESAMPLE 'SUNP1' SHOULD APPLY TO
!      THE NEW GRID, 0.2, 1, 2, 5, etc. IN FUTURE. jl+/-
      IRES=NINT(ABS(10*SFWHM))
      IF(IRES.GT.1)THEN
          IRES2=IRES*IRES
          IRESM1=IRES-1
          DO IV=IVLO+IRESM1,IVHI-IRESM1
              SUNSUM=IRES*SUNP1(IV)
              DO IVOFF=1,IRESM1
                  SUNSUM=SUNSUM                                         &
     &              +(IRES-IVOFF)*(SUNP1(IV+IVOFF)+SUNP1(IV-IVOFF))
              ENDDO

!             STORE DATA OFFSET BY IRESM1 CM-1 TO
!             AVOID OVERWRITING REQUIRED DATA.
              SUNP1(IV-IRESM1)=SUNSUM/IRES2
          ENDDO
      ELSE
          IRESM1=0
      ENDIF

!     MOVE DATA TO PROPER LOCATION (FREQUENCY).
!     ALSO, MULTIPLY BY SCALE HERE.
      DO IV=IVHI-IRESM1,IVLO+IRESM1,-1
          SUNP1(IV)=SUNP1(IV-IRESM1)*SSCALE
      ENDDO

!     DEFINE LOW FREQUENCY DATA [W CM-2 / CM-1] USING
!     POWER LAW APPROXIMATION.
!     UNNECESSARY IF USING BINARY INPUT FILE WHICH COVERS FROM 0
!     TO 50000 1/CM
      IF(INDEX(DFTSUN,'.BIN').GT.0 .OR. INDEX(DFTSUN,'.bin').GT.0)THEN
          SUNP1(0)=0.
          DO IV=1,IVLO+IRESM1-1
              SUNP1(IV)=A*FLOAT(IV)**B
          ENDDO

!         INSERT A CONSTANT VALUE AT HIGH FREQUENCIES.
          DO IV=IVHI-IRESM1+1,MSUNP1
              SUNP1(IV)=SUNP1(IVHI-IRESM1)
          ENDDO
      ENDIF

!     THE FOLLOWING CODING WAS USED TO CREATE THE BLOCK DATA S01BD.
!     THE LOGICAL FLAG LBLKDT IS HARD-WIRED TO FALSE, BUT THE CODING IS
!     LEFT HERE TO AID IN INCORPORATION OF NEW SOLAR IRRADIANCE DATA:
      IF(LBLKDT)THEN

!         OPEN S01BD.F FILE:
          CALL OPNFL(NSUN,0,'s01bd.f','UNKNOWN','FORMATTED','RDUSRS')

!         WRITE HEADER:
          WRITE(NSUN,'((A))')                                           &
     &      '      BLOCK DATA S01BD',                                   &
     &      ' ',                                                        &
     &      '!     SOLAR IRRADIANCES [W CM-2 / CM-1] TABULATED EACH',   &
     &      '!     TENTH WAVENUMBER FROM 0 TO 50,000 CM-1 AT 0.1 CM-1', &
     &      '!     SPECTRAL RESOLUTION (SFWHM RECTANGULAR SLIT).  THIS',&
     &      '!     DATA IS DERIVED FROM KURUCZ 0.01 CM-1 TABULATION',   &
     &      '!     OF TOP-OF-ATMOSPHERE SOLAR IRRADIANCES.',            &
     &      '      INCLUDE ''PARAMS.h''',                               &
     &      '      REAL SUNP1',                                         &
     &      '      COMMON/SOLP1/SUNP1(0:MSUNP1)',                       &
     &      '      INTEGER I'

!         WRITE 0 TO 100 CM-1 DATA:
          WRITE(NSUN,'(A,2(/A),25X,1P,E10.3,A)')                        &
     &      '!','!         0 TO   100 CM-1 DATA.',                      &
     &      '      DATA (SUNP1(I),I=     0,   100)/',SUNP1(0),','
          WRITE(NSUN,'((9(I6,6(1P,E10.3,A),/),A,6(1P,E10.3,A)))')       &
     &      (J,      (SUNP1(I),',',I=6*J-5,6*J),J=1,9),                 &
     &      '     &',(SUNP1(I),',',I=55   ,60 ),                        &
     &      (J-10,   (SUNP1(I),',',I=6*J-5,6*J),J=11,16),               &
     &       J-10,   (SUNP1(I),',',I=97,99),SUNP1(100),'/'

!         WRITE 101 TO 100*INT(MSUNP1/100) CM-1 DATA:
          IVOLD=100
          DO IV=200,MSUNP1,100
              WRITE(NSUN,'(A,2(/A,2(I6,A)))')                           &
     &          '!','!    ',IVOLD+1,' TO',IV,' CM-1 DATA.',             &
     &          '      DATA (SUNP1(I),I=',IVOLD+1,',',IV,')/'
              WRITE(NSUN,'((9(I6,6(1P,E10.3,A),/),A,6(1P,E10.3,A)))')   &
     &          (J      ,(SUNP1(IVOLD+I),',',I=6*J-5,6*J),J=1,9),       &
     &          '     &',(SUNP1(IVOLD+I),',',I=55,60),                  &
     &          (J-10   ,(SUNP1(IVOLD+I),',',I=6*J-5,6*J),J=11,16),     &
     &           J-10   ,(SUNP1(IVOLD+I),',',I=97,99),SUNP1(IV),'/'
              IVOLD=IV
          ENDDO

!         WRITE REMAINING CM-1 DATA:
          IF(IVOLD.LT.MSUNP1)THEN
              WRITE(NSUN,'(A,2(/A,2(I6,A)))')                           &
     &          '!','!    ',IVOLD+1,' TO',MSUNP1,' CM-1 DATA.',         &
     &          '      DATA (SUNP1(I),I=',IVOLD+1,',',MSUNP1,')/'
              LINE=1
   40         CONTINUE
              IF(IVOLD+6.LT.MSUNP1)THEN
                  IF(LINE.LT.10)THEN
                      WRITE(NSUN,'(I6,6(1P,E10.3,A))')LINE,             &
     &                  (SUNP1(IVOLD+I),',',I=1,6)
                      LINE=LINE+1
                  ELSE
                      WRITE(NSUN,'(A,6(1P,E10.3,A))')'     &',          &
     &                  (SUNP1(IVOLD+I),',',I=1,6)
                      LINE=1
                  ENDIF
                  IVOLD=IVOLD+6
                  GOTO 40
              ENDIF
              IF(LINE.LT.10)THEN
                  WRITE(NSUN,'(I6,6(1P,E10.3,A))')LINE,                 &
     &              (SUNP1(I),',',I=IVOLD+1,MSUNP1-1),SUNP1(MSUNP1),'/'
              ELSE
                  WRITE(NSUN,'(A,6(1P,E10.3,A))')'     &',              &
     &              (SUNP1(I),',',I=IVOLD+1,MSUNP1-1),SUNP1(MSUNP1),'/'
              ENDIF
          ENDIF
          WRITE(NSUN,'(A)')'      END'
          CLOSE(NSUN)
      ENDIF

!     RETURN UNLESS SFWHM IS NEGATIVE:
      IF(SFWHM.GE.0.)RETURN

!     OPEN OUTPUT DATA FILE:
      IF(USRSUN.EQ.'1' .OR. USRSUN.EQ.'5' .OR. USRSUN.EQ.'7')THEN
          SUNCHR='_sn_00p0.dat'
          SUNCHR(3:3)=USRSUN
      ELSE
          SUNCHR='_s0_00p0.dat'
      ENDIF
      WRITE(SUNCHR(5:6),'(I2.2)')IRES/10
      WRITE(SUNCHR(8:8),'(I1)')IRES-10*(IRES/10)
      IF(LNFLRT.LE.0)THEN
          CALL OPNFL(NSUN,0,SUNCHR,'UNKNOWN','FORMATTED','RDUSRS')
      ELSE

!         FULLNM MUST BE DEFINED BECAUSE OF GNU COMPILER PROBLEMS:
          FULLNM(1:LNFLRT+12)=FLRT(1:LNFLRT)//SUNCHR
          CALL OPNFL(NSUN,0,FULLNM(1:LNFLRT+12),                        &
     &      'UNKNOWN','FORMATTED','RDUSRS')
      ENDIF

!     WRITE DATA TO SUNCHR IN ORIGINAL UNITS [W CM-2 / CM-1]:
      WRITE(NSUN,'(A,/A,/(I7,1P,E13.3))')                               &
     &  '10XFREQ  SOLAR IRRADIANCE',                                    &
     &  ' (CM-1)  (W CM-2 / CM-1)',(IV,SUNP1(IV),IV=1,MSUNP1)
      CLOSE(NSUN)

!     RETURN TO DRIVER:
      RETURN

!     WRITE OUT ERROR MESSAGE AND STOP:
   50 CONTINUE
      WRITE(IPR,'(/A,/18X,A)')' Error in RSUNP1:  Problem'              &
     &  //' reading/using ASCII solar irradiance data file ',DFTSUN
      IF(LJMASS)CALL WRTBUF(FATAL)
      STOP 'Error in RSUNP1:  Problem reading/using ASCII solar data.'

!     WRITE OUT ERROR MESSAGE AND STOP:
   60 CONTINUE
      WRITE(IPR,'(/A,/18X,A)')' Error in RSUNP1:  Problem'              &
     &  //' reading/using BINARY solar irradiance data file ',DFTSUN
      IF(LJMASS)CALL WRTBUF(FATAL)
      STOP 'Error in RSUNP1:  Problem reading/using BINARY solar data.'
      END
