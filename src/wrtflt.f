      SUBROUTINE WRTFLT(IEMSCT,NMOLY,NLOS,DIS)

!     THIS ROUTINE PERFORMS THE FILTER RESPONSE SUMS FOR EACH CHANNEL.
      IMPLICIT NONE

!     ARGUMENTS:
!       IEMSCT   MODTRAN RADIATION TRANSPORT MODE FLAG
!                  0 FOR TRANSMITTANCE ONLY
!                  1 FOR THERMAL RADIANCE
!                  2 FOR THERMAL + SOLAR RADIANCE
!                  3 FOR TRANSMITTED SOLAR IRRADIANCE
!       NMOLY    NUMBER OF ACTIVE Y-SPECIES.
!       NLOS     NUMBER OF LINE-OF-SIGHT PATHS.
!       DIS      LOGICAL FLAG, TRUE FOR DISORT MULTIPLE SCATTERING.
      INTEGER IEMSCT,NMOLY,NLOS
      LOGICAL DIS

!     PARAMETERS:
!       ZERO     THE NUMBER 0.
      REAL ZERO
      PARAMETER(ZERO=0.)
      INCLUDE 'PARAMS.h'

!     COMMONS:
      INCLUDE 'CHANLS.h'
      INCLUDE 'IFIL.h'
      INCLUDE 'YPROPC.h'

!     /NAMEX/
!       CNAMEX   NAME OF CROSS-SECTION (X) SPECIES.
      CHARACTER CNAMEX*8
      COMMON/NAMEX/CNAMEX(NMOLX)
      SAVE /NAMEX/

!     /VRANGE/
!       VBNDMN   COMPUTATIONAL BANDPASS MINIMUM FREQUENCY [CM-1].
!       VBNDMX   COMPUTATIONAL BANDPASS MAXIMUM FREQUENCY [CM-1].
      REAL VBNDMN,VBNDMX
      COMMON/VRANGE/VBNDMN,VBNDMX

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     LOCAL VARIABLES:
!       UNITS    UNITS USED FOR WIDTH AND SPECTRAL MINIMUM AND MAXIMUM.
!       FRMT     FORMAT FOR TRANSMITTANCE DATA.
!       ICHAN    CHANNEL INDEX.
!       IMOL     MOLECULAR INDEX.
!       IOUT     CHANNEL OUTPUT VARIABLE INDEX.
!       NOUTR    NUMBER OF CHANNEL RADIANCE OUTPUT VARIABLES.
!       ILOS     LOOP INDEX FOR LINE-OF-SIGHT PATH.
!       LOPEN    FILE OPEN FLAG.
!       L_STAR   BRIGHTNESS TEMPERATURE LOW ACCURACY FLAG.
!       EXTINC   ONE MINUS TRANSMITTANCE.
      CHARACTER UNITS*121,FRMT*77
      INTEGER ICHAN,IMOL,IOUT,NOUTR,ILOS
      LOGICAL LOPEN,L_STAR
      REAL EXTINC

!     FUNCTIONS:
!       BT_wn    BRIGHTNESS TEMPERATURE FOR [W CM-1 SR-1] RADIANCE.
!       BT_uM    BRIGHTNESS TEMPERATURE FOR [W CM-2 SR-1 / uM] RADIANCE.
!       BT_NM    BRIGHTNESS TEMPERATURE FOR [W CM-2 SR-1 / NM] RADIANCE.
      REAL BT_WN,BT_UM,BT_NM

!     CHECK UNTFLG AND DEFINE UNITS:
      IF(UNTFLG.EQ.'M')THEN
          UNITS(1:37)='  (MICRONS)        NO.     (PER CM-1)'
          UNITS(38:121)='   (PER MICRON)  TEMP (K)  '//                 &
     &      '(W SR-1 CM-2)    (CM-1)   (MICRONS)  (MICRONS)  (MICRONS)'
      ELSEIF(UNTFLG.EQ.'N')THEN
          UNITS(1:37)='(NANOMETERS)       NO.     (PER CM-1)'
          UNITS(38:121)='       (PER NM)  TEMP (K)  '//                 &
     &      '(W SR-1 CM-2)    (CM-1)      (NM)       (NM)       (NM)  '
      ELSEIF(UNTFLG.EQ.'W')THEN
          UNITS(1:37)='(WAVENUMBERS)      NO.     (PER CM-1)'
          UNITS(38:121)='   (PER MICRON)  TEMP (K)  '//                 &
     &      '(W SR-1 CM-2)    (CM-1)   (MICRONS)    (CM-1)     (CM-1) '
      ELSE
          WRITE(IPR,'(/2A,/10X,A)')' WARNING:  No write'//              &
     &      ' to file "',CHNOUT(1:LENCHN)//'" because UNTFLG',          &
     &      ' in "CHANLS.h" does not equal "M", "N" or "W".'
          RETURN
      ENDIF

!     OPEN CHANNEL OUTPUT FILE IF NOT ALREADY OPENED.
      INQUIRE(ICHNUN,OPENED=LOPEN)
      IF(.NOT.LOPEN)                                                    &
     &  CALL OPNFL(ICHNUN,0,CHNOUT,'UNKNOWN','FORMATTED','WRTFLT')

!     HEADER AND RESULTS:
!       IF NPR=-1, PRINT ALL CHANNELS;
!       IF NPR= 0, ONLY PRINT CHANNELS WITH A POSITIVE VLCHAN; AND
!       IF NPR= 1, ONLY PRINT FULLY INTEGRATED CHANNELS.
      IF(IEMSCT.EQ.0)THEN
          FRMT='(/,000A)                             '                  &
     &      //'                                        '
          WRITE(FRMT(4:6),'(I3)')NMOLX+NMOLY+5
          WRITE(ICHNUN,FMT=FRMT)                                        &
     &      '1ST SPECTRAL CHAN      AVERAGE        CHANNEL  '/          &
     &      /'      FULL CHANNEL      SPECTRAL   SPECTRAL ',            &
     &      '      H2O      UNIFORMLY                 TRACE  '/         &
     &      /'      N2          H2O      MOLECULAR    AEROSOL ',        &
     &      '                AEROSOL                         '/         &
     &      /'                                                ',        &
     &      '                                                ',         &
     &      ('            ',IMOL=1,NMOLX+NMOLY),'  CHANNEL'
          FRMT(2:3)='  '
          IF(NMOLY.EQ.0)THEN
              WRITE(ICHNUN,FMT=FRMT)                                    &
     &          '   MOMENT     NEL    EXTINCTION     EXTINCTION '/      &
     &          /'    EQUIVALENT WIDTH     MINIMUM    MAXIMUM ',        &
     &          '  (NO CONTM)   MIX GASES      O3         GASES  '/     &
     &          /'   CONTINUUM   CONTINUUM  (RAYLEIGH)  PLUS CLOUD',    &
     &          '     HNO3     PLUS CLOUD      CO2         CO    '/     &
     &          /'      CH4         N2O         O2          NH3   ',    &
     &          '      NO          NO2         SO2        CLOUD  ',     &
     &          ('    '//CNAMEX(IMOL),IMOL=1,NMOLX),'  DESCRIPTION'
          ELSE
              WRITE(ICHNUN,FMT=FRMT)                                    &
     &          '   MOMENT     NEL    EXTINCTION     EXTINCTION '/      &
     &          /'    EQUIVALENT WIDTH     MINIMUM    MAXIMUM ',        &
     &          '  (NO CONTM)   MIX GASES      O3         GASES  '/     &
     &          /'   CONTINUUM   CONTINUUM  (RAYLEIGH)  PLUS CLOUD',    &
     &          '     HNO3     PLUS CLOUD      CO2         CO    '/     &
     &          /'      CH4         N2O         O2          NH3   ',    &
     &          '      NO          NO2         SO2        CLOUD  ',     &
     &          ('    '//CNAMEX(IMOL),IMOL=1,NMOLX),                    &
     &          ('    '//CNAMEY(IMOL),IMOL=1,NMOLY),'  DESCRIPTION'
          ENDIF
          WRITE(ICHNUN,FMT=FRMT)UNITS(1:13),                            &
     &      ' NO.   (1 - TRANSM)      (CM-1)   '//UNITS(78:121),        &
     &      '  ABSORBANCE  ABSORBANCE  ABSORBANCE  ABSORBANCE'/         &
     &      /'  ABSORBANCE  ABSORBANCE  SCATTERING  EXTINCTION',        &
     &      '  ABSORBANCE  ABSORBANCE  ABSORBANCE  ABSORBANCE'/         &
     &      /'  ABSORBANCE  ABSORBANCE  ABSORBANCE  ABSORBANCE',        &
     &      '  ABSORBANCE  ABSORBANCE  ABSORBANCE  EXTINCTION',         &
     &      ('  ABSORBANCE',IMOL=1,NMOLX+NMOLY)
          WRITE(ICHNUN,FMT=FRMT)                                        &
     &      '------------  ---  -------------  -------------'/          &
     &      /'  ---------  ---------  ---------  ---------',            &
     &      '  ----------  ----------  ----------  ----------'/         &
     &      /'  ----------  ----------  ----------  ----------',        &
     &      '  ----------  ----------  ----------  ----------'/         &
     &      /'  ----------  ----------  ----------  ----------',        &
     &      '  ----------  ----------  ----------  ----------',         &
     &      ('  ----------',IMOL=1,NMOLX+NMOLY),'  -----------'
          FRMT(1:42)='(F12.5,I5,F15.5,F15.4,4F11.4,00F12.5,2X,A)'
          WRITE(FRMT(30:31),'(I2.2)')NMOLX+NMOLY+20
          DO ICHAN=1,NCHAN
              IF((VLCHAN(1,1,ICHAN).GT.0.D0 .OR. NPR.LT.0)              &
     &          .AND. WDFREQ(ICHAN).GT.0.)THEN
                  IF(REAL(NFRQLO(ICHAN)).LT.VBNDMN .OR.                 &
     &              REAL(NFRQHI(ICHAN)).GT.VBNDMX)THEN
                      IF(NPR.LE.0)WRITE(ICHNUN,FMT=FRMT)MOM1ST(ICHAN),  &
     &                  -ICHAN,ZERO,VLCHAN(1,1,ICHAN),WDFREQ(ICHAN),    &
     &                  WDWAVE(ICHAN),SPECLO(ICHAN),SPECHI(ICHAN),      &
     &                  (VLCHAN(IOUT,1,ICHAN)/DBLE(WDFREQ(ICHAN)),IOUT= &
     &                  2,NMOLX+NMOLY+21),NMCHAN(ICHAN)(1:LNCHAN(ICHAN))
                  ELSE
                      EXTINC=REAL(VLCHAN(1,1,ICHAN))/WDFREQ(ICHAN)
                      IF(EXTINC.GT..999995)                             &
     &                  VLCHAN(1,1,ICHAN)=DBLE(WDFREQ(ICHAN))
                      WRITE(ICHNUN,FMT=FRMT)MOM1ST(ICHAN),              &
     &                  ICHAN,EXTINC,VLCHAN(1,1,ICHAN),WDFREQ(ICHAN),   &
     &                  WDWAVE(ICHAN),SPECLO(ICHAN),SPECHI(ICHAN),      &
     &                  (VLCHAN(IOUT,1,ICHAN)/DBLE(WDFREQ(ICHAN)),IOUT= &
     &                  2,NMOLX+NMOLY+21),NMCHAN(ICHAN)(1:LNCHAN(ICHAN))
                  ENDIF
              ENDIF
              DO IOUT=1,NMOLX+NMOLY+21
                  VLCHAN(IOUT,1,ICHAN)=0.D0
              ENDDO
          ENDDO
      ELSEIF(IEMSCT.EQ.3)THEN
          UNITS(65:77)='   (W CM-2)  '
          WRITE(ICHNUN,'(/(2A))')                                       &
     &      '1ST SPECTRAL CHAN    TRANSMITTED SPECTRAL SOLAR'           &
     &      //'   TRANSMITTED       FULL CHANNEL      SPECTRAL ',       &
     &      '  SPECTRAL    TOP OF ATMOS     SOLAR    CHANNEL',          &
     &      '   MOMENT     NEL    IRRADIANCE (W CM-2 / XXXX)'           &
     &      //'   SOLAR IRRAD.    EQUIVALENT WIDTH     MINIMUM ',       &
     &      '   MAXIMUM    SOLAR IRRAD.     PATH     DESCRIPTION',      &
     &      UNITS(1:13)//UNITS(19:54),                                  &
     &      UNITS(65:121)//'     (W CM-2)      TRANSM',                 &
     &      '------------  ---  -------------  -------------'           &
     &      //'  -------------  ---------  ---------  ---------',       &
     &      '  ---------  -------------  ----------  -----------'
          DO ICHAN=1,NCHAN
              IF((VLCHAN(2,1,ICHAN).GT.0.D0 .OR. NPR.LT.0)              &
     &          .AND. WDFREQ(ICHAN).GT.0.)THEN
                  IF(REAL(NFRQLO(ICHAN)).LT.VBNDMN .OR.                 &
     &              REAL(NFRQHI(ICHAN)).GT.VBNDMX)THEN
                      IF(NPR.LE.0)WRITE(ICHNUN,'(F12.5,I5,1P,3E15.6,    &
     &                  0P,4F11.4,1P,E15.6,0P,F12.7,2X,A)')             &
     &                  MOM1ST(ICHAN),-ICHAN,ZERO,ZERO,                 &
     &                  VLCHAN(1,1,ICHAN),WDFREQ(ICHAN),WDWAVE(ICHAN),  &
     &                  SPECLO(ICHAN),SPECHI(ICHAN),VLCHAN(2,1,ICHAN),  &
     &                  1-VLCHAN(3,1,ICHAN)/DBLE(WDFREQ(ICHAN)),        &
     &                  NMCHAN(ICHAN)(1:LNCHAN(ICHAN))
                  ELSE
                      WRITE(ICHNUN,'(F12.5,I5,1P,3E15.6,0P,4F11.4,      &
     &                  1P,E15.6,0P,F12.7,2X,A)')MOM1ST(ICHAN),ICHAN,   &
     &                  VLCHAN(1,1,ICHAN)/DBLE(WDFREQ(ICHAN)),          &
     &                  VLCHAN(1,1,ICHAN)/DBLE(WDWAVE(ICHAN)),          &
     &                  VLCHAN(1,1,ICHAN),WDFREQ(ICHAN),WDWAVE(ICHAN),  &
     &                  SPECLO(ICHAN),SPECHI(ICHAN),VLCHAN(2,1,ICHAN),  &
     &                  1-VLCHAN(3,1,ICHAN)/DBLE(WDFREQ(ICHAN)),        &
     &                  NMCHAN(ICHAN)(1:LNCHAN(ICHAN))
                  ENDIF
              ENDIF
              VLCHAN(1,1,ICHAN)=0.D0
              VLCHAN(2,1,ICHAN)=0.D0
              VLCHAN(3,1,ICHAN)=0.D0
          ENDDO
      ELSE
          FRMT='(F12.5,2I5,1P,2E15.6,0P,F10.3,A1,1P,E14.6,'             &
     &      //'0P,4F11.4,1P,3E15.6,0P,3F12.7,2X,A)'
          IF(IEMSCT.EQ.1)THEN
              IF(DIS)THEN
                  WRITE(ICHNUN,'(/(3A))')                               &
     &              '1ST SPECTRAL  LOS CHAN        SPECTRAL '//         &
     &              ' RADIANCE       BRIGHT-      CHANNEL        FULL ',&
     &              'CHANNEL      SPECTRAL   SPECTRAL    PATH + SCAT '//&
     &              '  TRANSM GROUND   TOTAL TRANSM    SENSOR  ',       &
     &              '     SURFACE  CHANNEL',                            &
     &              '   MOMENT     NO.  NEL       (W SR-1 CM'//         &
     &              '-2 / XXXX)       NESS       RADIANCE      EQUIVAL',&
     &              'ENT WIDTH     MINIMUM    MAXIMUM      EMISSION  '//&
     &              '     EMISSION     GRND REFLECT     PATH   ',       &
     &              ' DIRECTIONAL  DESCRIPTION',                        &
     &              UNITS,'  (W SR-1 CM-2)  (W SR-1 CM-2)  ',           &
     &              '(W SR-1 CM-2)    TRANSM    EMISSIVITY',            &
     &              '------------  ---  ---  -------------  '//         &
     &              '-------------  --------  -------------  ---------',&
     &              '  ---------  ---------  ---------  -------------'//&
     &              '  -------------  -------------  ----------',       &
     &              '  ----------  -----------'
                  NOUTR=4
                  FRMT(53:71)='1P,3E15.6,0P,2F12.7'
              ELSE
                  WRITE(ICHNUN,'(/(3A))')                               &
     &              '1ST SPECTRAL  LOS CHAN        SPECTRAL '//         &
     &              ' RADIANCE       BRIGHT-      CHANNEL        FULL ',&
     &              'CHANNEL      SPECTRAL   SPECTRAL    PATH + SCAT '//&
     &              '    SCATTERED    TRANSM GROUND   TOTAL TRANSM',    &
     &              '    SENSOR       SURFACE  CHANNEL',                &
     &              '   MOMENT     NO.  NEL       (W SR-1 CM'//         &
     &              '-2 / XXXX)       NESS       RADIANCE      EQUIVAL',&
     &              'ENT WIDTH     MINIMUM    MAXIMUM      EMISSION  '//&
     &              '     EMISSION       EMISSION     GRND REFLECT',    &
     &              '     PATH    DIRECTIONAL  DESCRIPTION',            &
     &              UNITS,'  (W SR-1 CM-2)  (W SR-1 CM-2)  (W SR-1',    &
     &              ' CM-2)  (W SR-1 CM-2)    TRANSM    EMISSIVITY',    &
     &              '------------  ---  ---  -------------  '//         &
     &              '-------------  --------  -------------  ---------',&
     &              '  ---------  ---------  ---------  -------------'//&
     &              '  -------------  -------------  -------------',    &
     &              '  ----------  ----------  -----------'
                  NOUTR=5
                  FRMT(53:71)='1P,4E15.6,0P,2F12.7'
              ENDIF
          ELSEIF(DIS)THEN
              WRITE(ICHNUN,'(/(3A))')                                   &
     &          '1ST SPECTRAL  LOS CHAN        SPECTRAL  RADIANCE    '  &
     &          //'   BRIGHT-      CHANNEL        FULL CHANNEL    ',    &
     &          '  SPECTRAL   SPECTRAL    PATH + SCAT   TRANSM GROUND'  &
     &          //'    PATH TOTAL    PATH SINGLE    TOTAL TRANSM',      &
     &          '  DIRECT TRANSM   TRANSM SOLAR   TRANSM SOLAR'         &
     &          //'    SENSOR       SURFACE  CHANNEL',                  &
     &          '   MOMENT     NO.  NEL       (W SR-1 CM-2 / XXXX)   '  &
     &          //'    NESS       RADIANCE      EQUIVALENT WIDTH  ',    &
     &          '   MINIMUM    MAXIMUM      EMISSION       EMISSION  '  &
     &          //'    SCAT SOLAR     SCAT SOLAR    GRND REFLECT',      &
     &          '   GRND REFLECT   LOS+SUN PATH     TO SENSOR '         &
     &          //'     PATH    DIRECTIONAL  DESCRIPTION'
              WRITE(ICHNUN,'((3A))')UNITS,                              &
     &          '  (W SR-1 CM-2)  (W SR-1 CM-2)'                        &
     &          //'  (W SR-1 CM-2)  (W SR-1 CM-2)  (W SR-1 CM-2)',      &
     &          '  (W SR-1 CM-2)     (W CM-2)       (W CM-2)  '         &
     &          //'    TRANSM    EMISSIVITY',                           &
     &          '------------  ---  ---  -------------  -------------'  &
     &          //'  --------  -------------  ---------  ---------',    &
     &          '  ---------  ---------  -------------  -------------'  &
     &          //'  -------------  -------------  -------------',      &
     &          '  -------------  -------------  -------------'         &
     &          //'  ----------  ----------  -----------'
              NOUTR=9
              FRMT(53:71)='1P,8E15.6,0P,2F12.7'
          ELSE
              WRITE(ICHNUN,'(/(3A))')                                   &
     &          '1ST SPECTRAL  LOS CHAN        SPECTRAL  RADIANCE    '  &
     &          //'   BRIGHT-      CHANNEL        FULL CHANNEL    ',    &
     &          '  SPECTRAL   SPECTRAL    PATH + SCAT     SCATTERED  '  &
     &          //'  TRANSM GROUND    PATH TOTAL    PATH SINGLE ',      &
     &          '   TOTAL TRANSM  DIRECT TRANSM   TRANSM SOLAR  '       &
     &          //' TRANSM SOLAR    SENSOR       SURFACE  CHANNEL',     &
     &          '   MOMENT     NO.  NEL       (W SR-1 CM-2 / XXXX)   '  &
     &          //'    NESS       RADIANCE      EQUIVALENT WIDTH  ',    &
     &          '   MINIMUM    MAXIMUM      EMISSION       EMISSION  '  &
     &          //'     EMISSION      SCAT SOLAR     SCAT SOLAR ',      &
     &          '   GRND REFLECT   GRND REFLECT   LOS+SUN PATH  '       &
     &          //'   TO SENSOR      PATH    DIRECTIONAL  DESCRIPTION'
              WRITE(ICHNUN,'((3A))')UNITS,                              &
     &          '  (W SR-1 CM-2)  (W SR-1 CM-2)  (W SR-1 CM-2)'         &
     &          //'  (W SR-1 CM-2)  (W SR-1 CM-2)  (W SR-1 CM-2)',      &
     &          '  (W SR-1 CM-2)     (W CM-2)       (W CM-2)  '         &
     &          //'    TRANSM    EMISSIVITY',                           &
     &          '------------  ---  ---  -------------  -------------'  &
     &          //'  --------  -------------  ---------  ---------',    &
     &          '  ---------  ---------  -------------  -------------'  &
     &          //'  -------------  -------------  -------------',      &
     &          '  -------------  -------------  -------------  '       &
     &          //'-------------  ----------  ----------  -----------'
              NOUTR=10
              FRMT(53:71)='1P,9E15.6,0P,2F12.7'
          ENDIF
          L_STAR=.FALSE.
          DO ILOS=1,NLOS
              DO ICHAN=1,NCHAN
                  IF((VLCHAN(1,ILOS,ICHAN).GT.0.D0 .OR. NPR.LT.0)       &
     &              .AND. WDFREQ(ICHAN).GT.0.)THEN
                      IF(REAL(NFRQLO(ICHAN)).LT.VBNDMN .OR.             &
     &                  REAL(NFRQHI(ICHAN)).GT.VBNDMX)THEN
                          IF(NPR.LE.0)WRITE(ICHNUN,FMT=FRMT)            &
     &                      MOM1ST(ICHAN),ILOS,-ICHAN,ZERO,ZERO,ZERO,   &
     &                      ' ',REAL(VLCHAN(1,ILOS,ICHAN)),             &
     &                      WDFREQ(ICHAN),WDWAVE(ICHAN),                &
     &                      SPECLO(ICHAN),SPECHI(ICHAN),                &
     &                      (REAL(VLCHAN(IOUT,ILOS,ICHAN)),             &
     &                                                  IOUT=2,NOUTR),  &
     &                      MAX(0.D0,1-VLCHAN(NOUTR+1,ILOS,ICHAN)/      &
     &                      DBLE(WDFREQ(ICHAN))),                       &
     &                      REAL(VLCHAN(NOUTR+2,ILOS,ICHAN))/           &
     &                      WDFREQ(ICHAN),NMCHAN(ICHAN)(1:LNCHAN(ICHAN))
                      ELSEIF(UNTFLG.EQ.'W')THEN
                          WRITE(ICHNUN,FMT=FRMT)                        &
     &                      MOM1ST(ICHAN),ILOS,ICHAN,                   &
     &                      REAL(VLCHAN(1,ILOS,ICHAN))/WDFREQ(ICHAN),   &
     &                      REAL(VLCHAN(1,ILOS,ICHAN))/WDWAVE(ICHAN),   &
     &                      BT_WN(REAL(VLCHAN(1,ILOS,ICHAN))            &
     &                                                /WDFREQ(ICHAN),   &
     &                      MOM1ST(ICHAN),M2DIFF(ICHAN)),C_STAR(ICHAN), &
     &                      REAL(VLCHAN(1,ILOS,ICHAN)),WDFREQ(ICHAN),   &
     &                      WDWAVE(ICHAN),SPECLO(ICHAN),SPECHI(ICHAN),  &
     &                      (REAL(VLCHAN(IOUT,ILOS,ICHAN)),             &
     &                                                  IOUT=2,NOUTR),  &
     &                      MAX(1-VLCHAN(NOUTR+1,ILOS,ICHAN)            &
     &                                     /DBLE(WDFREQ(ICHAN)),0.D0),  &
     &                      REAL(VLCHAN(NOUTR+2,ILOS,ICHAN))            &
     &                                                 /WDFREQ(ICHAN),  &
     &                      NMCHAN(ICHAN)(1:LNCHAN(ICHAN))
                          IF(C_STAR(ICHAN).EQ.'*')L_STAR=.TRUE.
                      ELSEIF(UNTFLG.EQ.'M')THEN
                          WRITE(ICHNUN,FMT=FRMT)                        &
     &                      MOM1ST(ICHAN),ILOS,ICHAN,                   &
     &                      REAL(VLCHAN(1,ILOS,ICHAN))/WDFREQ(ICHAN),   &
     &                      REAL(VLCHAN(1,ILOS,ICHAN))/WDWAVE(ICHAN),   &
     &                      BT_UM(REAL(VLCHAN(1,ILOS,ICHAN))            &
     &                                                /WDWAVE(ICHAN),   &
     &                      MOM1ST(ICHAN),M2DIFF(ICHAN)),C_STAR(ICHAN), &
     &                      REAL(VLCHAN(1,ILOS,ICHAN)),WDFREQ(ICHAN),   &
     &                      WDWAVE(ICHAN),SPECLO(ICHAN),SPECHI(ICHAN),  &
     &                      (REAL(VLCHAN(IOUT,ILOS,ICHAN)),             &
     &                                                  IOUT=2,NOUTR),  &
     &                      MAX(1-VLCHAN(NOUTR+1,ILOS,ICHAN)            &
     &                                     /DBLE(WDFREQ(ICHAN)),0.D0),  &
     &                      REAL(VLCHAN(NOUTR+2,ILOS,ICHAN))            &
     &                                                 /WDFREQ(ICHAN),  &
     &                      NMCHAN(ICHAN)(1:LNCHAN(ICHAN))
                          IF(C_STAR(ICHAN).EQ.'*')L_STAR=.TRUE.
                      ELSEIF(UNTFLG.EQ.'N')THEN
                          WRITE(ICHNUN,FMT=FRMT)                        &
     &                      MOM1ST(ICHAN),ILOS,ICHAN,                   &
     &                      REAL(VLCHAN(1,ILOS,ICHAN))/WDFREQ(ICHAN),   &
     &                      REAL(VLCHAN(1,ILOS,ICHAN))/WDWAVE(ICHAN),   &
     &                      BT_NM(REAL(VLCHAN(1,ILOS,ICHAN))            &
     &                                                /WDWAVE(ICHAN),   &
     &                      MOM1ST(ICHAN),M2DIFF(ICHAN)),C_STAR(ICHAN), &
     &                      REAL(VLCHAN(1,ILOS,ICHAN)),WDFREQ(ICHAN),   &
     &                      WDWAVE(ICHAN),SPECLO(ICHAN),SPECHI(ICHAN),  &
     &                      (REAL(VLCHAN(IOUT,ILOS,ICHAN)),             &
     &                                                  IOUT=2,NOUTR),  &
     &                      MAX(1-VLCHAN(NOUTR+1,ILOS,ICHAN)            &
     &                                     /DBLE(WDFREQ(ICHAN)),0.D0),  &
     &                      REAL(VLCHAN(NOUTR+2,ILOS,ICHAN))            &
     &                                                 /WDFREQ(ICHAN),  &
     &                      NMCHAN(ICHAN)(1:LNCHAN(ICHAN))
                          IF(C_STAR(ICHAN).EQ.'*')L_STAR=.TRUE.
                      ENDIF
                  ENDIF
                  DO IOUT=1,NOUTR+2
                      VLCHAN(IOUT,ILOS,ICHAN)=0.D0
                  ENDDO
              ENDDO
          ENDDO
          IF(L_STAR)WRITE(ICHNUN,'(/A)')'Brightness temperatures with'  &
     &      //' a "*" suffix have an uncertainty greater than 0.1%'
      ENDIF

!     TABLE COMPLETE:
      RETURN
      END
