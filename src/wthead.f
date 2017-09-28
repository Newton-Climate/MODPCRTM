      SUBROUTINE WTHEAD(IMULT,NLOS,DIS,DISAZM,TRANSM)

!     WRITE HEADER FOR SPECTRAL TABLES.
      IMPLICIT NONE

!     INPUT ARGUMENTS:
!       IMULT    MULTIPLE SCATTERING FLAG
!                  (0=NONE, 1=AT SENSOR, -1=AT FINAL OR TANGENT POINT).
!       NLOS     NUMBER OF LINE-OF-SIGHT PATHS.
!       DIS      LOGICAL FLAG, TRUE FOR DISORT MULTIPLE SCATTERING.
!       DISAZM   LOGICAL FLAG, TRUE FOR DISORT WITH AZIMUTH DEPENDENCE.
!       TRANSM   FLAG [TRUE FOR TRANSMITTANCE ONLY CALCULATIONS].
      INTEGER IMULT,NLOS
      LOGICAL DIS,DISAZM,TRANSM

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     COMMONS:
      INCLUDE 'IFIL.h'
      INCLUDE 'YPROP.h'
      INCLUDE 'YPROPC.h'

!     /CARD1/
!       MODEL    MODEL ATMOSPHERE INDEX.
!       ITYPE    SLANT PATH TYPE.
!       IEMSCT   RADIATIVE TRANSFER MODE.
!                  0 FOR TRANSMITTANCE
!                  1 FOR THERMAL EMISSION ONLY
!                  2 FOR THERMAL EMISSION PLUS SOLAR SCATTER
!                  3 FOR TRANSMITTED SOLAR IRRADIANCE
!                  4 FOR SOLAR SCATTER ONLY
!       M1       MODEL ATMOSPHERE FOR PRESSURE & TEMPERATURE PROFILES.
!       M2       MODEL ATMOSPHERE FOR H2O PROFILE.
!       M3       MODEL ATMOSPHERE FOR O3 PROFILE.
!       I_RD2C   READ CARD 2C, 2C1, ... IF EQUAL 1; SKIP IF EQUAL TO 0.
!       NOPRNT   PRINT FLAG.
!       MODTRN   MODTRAN BAND MODEL FLAG.
      INTEGER MODEL,ITYPE,IEMSCT,M1,M2,M3,I_RD2C,NOPRNT
      LOGICAL MODTRN
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,I_RD2C,NOPRNT,MODTRN

!     /NAMEX/
!       CNAMEX   NAME OF CROSS-SECTION (X) SPECIES.
      CHARACTER CNAMEX*8
      COMMON/NAMEX/CNAMEX(NMOLX)
      SAVE /NAMEX/

!     /CSCAN/
!       CHUNIT   UNIT FLAG ('W'=WAVENUMBERS;'M'=MICRONS;'N'=NANOMETERS).
!       RELABS   SPECTRAL RESOLUTION FLAG('A'=ABSOLUTE;'R'=RELATIVE[%]).
!       LNFEED   LINE FEED FLAG FOR .FLX FILE ('T' FOR 80 CHARACTER
!                  LINES, 'F' FOR LONG LINES, ' ' FOR NO .FLX FILE).
      CHARACTER CHUNIT*1,RELABS*1,LNFEED*1
      COMMON/CSCAN/CHUNIT,RELABS,LNFEED

!     /SPCCTL/
!       DODRGD   .TRUE. FOR NON-DEFAULT SPECTRAL ANALYSIS.
      LOGICAL DODGRD
      COMMON/SPCCTL/DODGRD

!     /CFLAGS/
!       YFLAG    Y COORDINATE FLAG FOR PLOT.DAT FILE
!                  = "T" FOR TRANSMITTANCE
!                  = "R" FOR RADIANCE (IRRADIANCE FOR IEMSCT=3)
!                  = "N" FOR NO PLOT.DAT OUTPUT
!       XFLAG    X COORDINATE FLAG FOR PLOT.DAT FILE
!                  = "W" FOR FREQUENCY IN WAVENUMBERS (CM-1) AND
!                        RADIANCE IN W SR-1 CM-2 / CM-1
!                  = "M" FOR WAVELENGTH IN MICRONS AND
!                        RADIANCE IN W SR-1 CM-2 / MICRON
!                  = "N" FOR WAVELENGTH IN NANOMETERS AND
!                        RADIANCE IN MICRO-WATTS SR-1 CM-2 / NANOMETER
!       DLIMIT   DELIMITER CHARACTER STRING BETWEEN MODTRAN RUNS
!       FLAGS    SCANNING FUNCTION FLAGS.
      CHARACTER YFLAG*1,XFLAG*1,DLIMIT*8,FLAGS*7
      COMMON/CFLAGS/YFLAG,XFLAG,DLIMIT,FLAGS

!     LOCAL VARIABLES:
!       K        SPECIES LOOP INDEX
      INTEGER K
      IF(DODGRD)THEN
          WRITE(JPUSCR)BINSPR
          WRITE(JPUSCR)BINSPR
          WRITE(JPUSCR)'FREQ',IEMSCT,NMOLX,NMOLY,IMULT,DIS
      ENDIF

!     BRANCH BASED ON RADIATIVE TRANSFER MODE:
      IF(IEMSCT.EQ.0)THEN
          IF(NMOLY.EQ.0)THEN
              WRITE(IPU,'(100A)')'    FREQ COMBIN    H2O   UMIX',       &
     &          '     O3  TRACE     N2    H2O MOLEC AER+CLD',           &
     &          '  HNO3 AER+CLD    -LOG    CO2     CO    CH4    N2O',   &
     &          '     O2    NH3     NO    NO2    SO2  CLOUD',           &
     &          (CNAMEX(K)(2:8),K=1,NMOLX)
          ELSE
              WRITE(IPU,'(100A)')'    FREQ COMBIN    H2O   UMIX',       &
     &          '     O3  TRACE     N2    H2O MOLEC AER+CLD',           &
     &          '  HNO3 AER+CLD    -LOG    CO2     CO    CH4    N2O',   &
     &          '     O2    NH3     NO    NO2    SO2  CLOUD',           &
     &          (CNAMEX(K)(2:8),K=1,NMOLX),(CNAMEY(K)(2:8),K=1,NMOLY)
          ENDIF
          WRITE(IPU,'(100A)')'    CM-1  TRANS  TRANS  TRANS',           &
     &      '  TRANS  TRANS   CONT   CONT   SCAT  TRANS',               &
     &      '  TRANS abTRNS  COMBIN',('  TRANS',K=1,NMOLX+NMOLY+10)

          IF(DODGRD)THEN
              IF(CHUNIT.EQ.'W')THEN

!                 WAVENUMBERS:
                  WRITE(IPUSCN,'(103A)')' FREQ (CM-1)    COMBIN    H2O',&
     &              '   UMIX     O3  TRACE     N2 H2Ocon MOLEC AER+CLD',&
     &              '  HNO3 AER+CLD    -LOG    CO2     CO    CH4    ',  &
     &              'N2O     O2    NH3     NO    NO2    SO2  CLOUD',    &
     &              (CNAMEX(K)(2:8),K=1,NMOLX),                         &
     &              (CNAMEY(K)(2:8),K=1,NMOLY)
              ELSEIF(CHUNIT.EQ.'M')THEN

!                 MICRONS:
                  WRITE(IPUSCN,'(103A)')' WAVLEN (uM)    COMBIN    H2O',&
     &              '   UMIX     O3  TRACE     N2 H2Ocon MOLEC AER+CLD',&
     &              '  HNO3 AER+CLD    -LOG    CO2     CO    CH4    ',  &
     &              'N2O     O2    NH3     NO    NO2    SO2  CLOUD',    &
     &              (CNAMEX(K)(2:8),K=1,NMOLX),                         &
     &              (CNAMEY(K)(2:8),K=1,NMOLY)
              ELSE

!                 NANOMETERS:
                  WRITE(IPUSCN,'(103A)')' WAVLEN (NM)    COMBIN    H2O',&
     &              '   UMIX     O3  TRACE     N2 H2Ocon MOLEC AER+CLD',&
     &              '  HNO3 AER+CLD    -LOG    CO2     CO    CH4    ',  &
     &              'N2O     O2    NH3     NO    NO2    SO2  CLOUD',    &
     &              (CNAMEX(K)(2:8),K=1,NMOLX),                         &
     &              (CNAMEY(K)(2:8),K=1,NMOLY)

              ENDIF

              DO K=1,NMOLX
                  WRITE(JPUSCR)CNAMEX(K)(2:8)
              ENDDO
              DO K=1,NMOLY
                  WRITE(JPUSCR)CNAMEY(K)(2:8)
              ENDDO
          ENDIF
      ELSEIF(IEMSCT.EQ.3)THEN
          WRITE(IPU,'(A)')'    FREQ   TRANS     SOL TR  SOLAR'
          IF(DODGRD)THEN
              IF(CHUNIT.EQ.'W')THEN
                  WRITE(IPUSCN,'(2A)')'  FREQ(CM-1)   TRAN  SOL TR  ',  &
     &              'SOLAR'
              ELSEIF(CHUNIT.EQ.'M')THEN
                  WRITE(IPUSCN,'(2A)')'WAVLEN(MCRN)   TRAN  SOL TR  ',  &
     &              'SOLAR'
              ELSE
                  WRITE(IPUSCN,'(2A)')'  WAVLEN(NM)   TRAN  SOL TR  ',  &
     &              'SOLAR'
              ENDIF
          ENDIF
      ELSE
          IF(NLOS.EQ.1)THEN
              WRITE(IPU,'(2A)')'    FREQ  TOT_TRANS  PTH_THRML'         &
     &          //'  THRML_SCT  SURF_EMIS   SOL_SCAT  SING_SCAT',       &
     &          '  GRND_RFLT  DRCT_RFLT  TOTAL_RAD  REF_SOL'            &
     &          //'  SOL@OBS   DEPTH DIR_EM    TOA_SUN BBODY_T[K]'
          ELSE
              WRITE(IPU,'(2A)')'    FREQ LOS  TOT_TRANS  PTH_THRML'     &
     &          //'  THRML_SCT  SURF_EMIS   SOL_SCAT  SING_SCAT',       &
     &          '  GRND_RFLT  DRCT_RFLT  TOTAL_RAD  REF_SOL'            &
     &          //'  SOL@OBS   DEPTH DIR_EM    TOA_SUN BBODY_T[K]'
          ENDIF
          IF(DODGRD)THEN
              IF(CHUNIT.EQ.'W')THEN
                  WRITE(IPUSCN,'(4A)')                                  &
     &              '  FREQ(CM-1)      TRAN  PTH_THRML  THRML_SCT  ',   &
     &              'SURF_EMIS   SOL_SCAT  SING_SCAT  GRND_RFLT  ',     &
     &              'DRCT_RFLT  TOTAL_RAD  REF_SOL  SOL@OBS   DEPTH ',  &
     &              'DIR_EM    TOA_SUN BBODY_T[K]'
              ELSEIF(CHUNIT.EQ.'M')THEN
                  WRITE(IPUSCN,'(4A)')                                  &
     &              'WAVLEN(MCRN)      TRAN  PTH_THRML  THRML_SCT  ',   &
     &              'SURF_EMIS   SOL_SCAT  SING_SCAT  GRND_RFLT  ',     &
     &              'DRCT_RFLT  TOTAL_RAD  REF_SOL  SOL@OBS   DEPTH ',  &
     &              'DIR_EM    TOA_SUN BBODY_T[K]'
              ELSE
                  WRITE(IPUSCN,'(4A)')                                  &
     &              '  WAVLEN(NM)      TRAN  PTH_THRML  THRML_SCT  ',   &
     &              'SURF_EMIS   SOL_SCAT  SING_SCAT  GRND_RFLT  ',     &
     &              'DRCT_RFLT  TOTAL_RAD  REF_SOL  SOL@OBS   DEPTH ',  &
     &              'DIR_EM    TOA_SUN BBODY_T[K]'
              ENDIF
          ENDIF
      ENDIF

      IF(BINOUT)THEN
          WRITE(JPU)                                                    &
     &      -9998,IEMSCT,NMOLX,NMOLY,IMULT,DIS,DISAZM,NLOS,YFLAG,XFLAG

!     BINARY FILES DO NOT OUTPUT TO IPR1
      ELSEIF(NOPRNT.LE.-1)THEN
          IF(NLOS.GT.1)THEN
              WRITE(IPR1,'((3A))')                                      &
     &          ' FREQ    LOS  ALT     TOTAL    DELTA   UP-DN-DIR',     &
     &          '   THRML_UP   THRML_DN  THRML_SRC  THRML_SUM',         &
     &          '   SOLAR_UP   SOLAR_DN  SOLAR_SRC  SOLAR_SUM',         &
     &          '(CM-1)       (KM)     TRANS    TRANS (          ',     &
     &          '                                 W CM-2 / CM',         &
     &          '-1                                         )'
          ELSEIF(IMULT.NE.0)THEN
              WRITE(IPR1,'((3A))')                                      &
     &          ' FREQ     ALT     TOTAL    DELTA   UP-DN-DIR',         &
     &          '   THRML_UP   THRML_DN  THRML_SRC  THRML_SUM',         &
     &          '   SOLAR_UP   SOLAR_DN  SOLAR_SRC  SOLAR_SUM',         &
     &          '(CM-1)   (KM)     TRANS    TRANS (          ',         &
     &          '                                 W CM-2 / CM',         &
     &          '-1                                         )'
          ELSEIF(.NOT.TRANSM)THEN
              WRITE(IPR1,'((1X,2A))')                                   &
     &          '             ALTITUDES                 B(V,T)  ',      &
     &          '           TRANSMISSION            RADIANCE',          &
     &          '  FREQ  BEGINNING  ENDING  INT   LAYER     BOUNDARY  ',&
     &          ' TO BEGIN   IN LAYER    LAYER      TOTAL',             &
     &          '(CM-1)    (KM)     (KM)         (W SR-1 CM-2 / CM-1) ',&
     &          '                      (W SR-1 CM-2 / CM-1)'
          ENDIF
      ENDIF
      RETURN
      END

      SUBROUTINE T6HEAD(IPR,IEMSCT,IMULT,NLOS,DIS)

!     T6HEAD WRITES OUT SPECTRAL TABLE HEADER FOR .TP6 FILE.
      IMPLICIT NONE

!     INPUT ARGUMENTS:
!       IPR      UNIT NUMBER OF STANDARD OUTPUT FILE, "<rootname>.tp6".
!       IEMSCT   RADIATIVE TRANSFER MODE
!                  0 FOR TRANSMITTANCE
!                  1 FOR THERMAL EMISSION
!                  2 FOR THERMAL EMISSION PLUS SOLAR SCATTER
!                  3 FOR TRANSMITTED SOLAR IRRADIANCE
!       IMULT    MULTIPLE SCATTERING FLAG
!                  (0=NONE, 1=AT SENSOR, -1=AT FINAL OR TANGENT POINT).
!       NLOS     NUMBER OF LINE-OF-SIGHT PATHS.
!       DIS      LOGICAL FLAG, TRUE FOR DISORT MULTIPLE SCATTERING.
      INTEGER IPR,IEMSCT,IMULT,NLOS
      LOGICAL DIS
      IF(IEMSCT.EQ.0)THEN
          WRITE(IPR,'(//(3A))')                                         &
     &      '   FREQ WAVELENGTH    TOTAL   H2O     CO2+  ',             &
     &      '   OZONE    TRACE  N2 CONT  H2O CONT MOL SCAT',            &
     &      '  AER&CLD  HNO3    AER&CLD  INTEGRATED',                   &
     &      '   1/CM   MICRONS     TRANS  TRANS    TRANS ',             &
     &      '   TRANS    TRANS   TRANS    TRANS    TRANS  ',            &
     &      '   TRANS   TRANS   TAU-ABS  ABSORPTION'
      ELSEIF(IEMSCT.EQ.1)THEN
          IF(NLOS.EQ.1)THEN
              WRITE(IPR,'(/52X,A,3(/3A))')                              &
     &          'RADIANCE(WATTS/CM2-STER-XXX)',                         &
     &          '   FREQ  WAVELENGTH  DIREC      PATH_THERMAL   ',      &
     &          'SCAT_PART    SURFACE_EMISSION   SURFACE_REFLECTED',    &
     &          '    TOTAL_RADIANCE  INTEGRAL   TOTAL',                 &
     &          '  (CM-1) (MICRON)     EMIS    (CM-1)   (MICRN) ',      &
     &          '   (CM-1)    (CM-1)   (MICRN)    (CM-1)   (MICRN)',    &
     &          '   (CM-1)   (MICRN)            TRANS'
          ELSE
              WRITE(IPR,'(/52X,A,3(/3A))')                              &
     &          'RADIANCE(WATTS/CM2-STER-XXX)',                         &
     &          '   FREQ  WAVELENGTH  LOS  DIREC      PATH_THERMAL   ', &
     &          'SCAT_PART    SURFACE_EMISSION   SURFACE_REFLECTED',    &
     &          '    TOTAL_RADIANCE  INTEGRAL   TOTAL',                 &
     &          '  (CM-1) (MICRON)          EMIS    (CM-1)   (MICRN) ', &
     &          '   (CM-1)    (CM-1)   (MICRN)    (CM-1)   (MICRN)',    &
     &          '   (CM-1)   (MICRN)            TRANS'
          ENDIF
      ELSEIF(IEMSCT.EQ.3)THEN
          WRITE(IPR,'(/22X,A,//(2A))')                                  &
     &      'IRRADIANCE (WATTS/CM2-XXXX)',                              &
     &      '   FREQ  WAVELENGTH      TRANSMITTED     ',                &
     &      '      SOLAR           INTEGRATED         TOTAL',           &
     &      '  (CM-1) (MICRON)     (CM-1)    (MICRN)  ',                &
     &      ' (CM-1)    (MICRN)   TRANS.    SOLAR     TRANS'
      ELSEIF(IMULT.EQ.0)THEN
          WRITE(IPR,'(/27X,2A,/3(/3A))')                                &
     &      'RADIANCE (WATTS/CM2-STER-XXX)',                            &
     &      '   [No STER in TOA SOLAR IRRADIANCE]',                     &
     &      '  FREQ   WAVELENGTH    PATH_THERMAL    SURFACE_',          &
     &      'EMISSION  SOL_IRR    SINGLE_SCATTER  GROUND_',             &
     &      'REFLECTED   TOTAL_RADIANCE  INTEGRAL    TOTAL',            &
     &      ' (CM-1) (MICRON)      (CM-1)  (MICRN)   (CM-1) ',          &
     &      ' (MICRN)   (CM-1)   (CM-1)  (MICRN)   (CM-1)',             &
     &      '  (MICRN)   (CM-1)  (MICRN)             TRANS'
      ELSEIF(DIS)THEN
          IF(NLOS.EQ.1)THEN
              WRITE(IPR,'(/45X,A,/4(/3A))')                             &
     &          'RADIANCE(WATTS/CM2-STER-XXX)',                         &
     &          '  FREQ  WAVELENGTH     PATH_THERMAL    SURFACE   ',    &
     &          '   PATH_SCATTERED_SOLAR   GROUND_REFLECTED_',          &
     &          'RADIANCE   TOTAL_RADIANCE   INTEGRAL    TOTAL',        &
     &          '                                       EMISSION  ',    &
     &          '   TOTAL RAD      SINGLE        TOTAL      ',          &
     &          '  DIRECT                                TRANS',        &
     &          ' (CM-1) (MICRON)       (CM-1)  (MICRN)   (CM-1)  ',    &
     &          ' (CM-1)  (MICRN)   (CM-1)   (CM-1)  (MICRN)',          &
     &          '   (CM-1)   (CM-1)  (MICRN)'
          ELSE
              WRITE(IPR,'(/45X,A,/4(/3A))')                             &
     &          'RADIANCE(WATTS/CM2-STER-XXX)',                         &
     &          '  FREQ  WAVELENGTH    LOS   PATH_THERMAL    SURFACE ', &
     &          '     PATH_SCATTERED_SOLAR   GROUND_REFLECTED_',        &
     &          'RADIANCE   TOTAL_RADIANCE   INTEGRAL    TOTAL',        &
     &          '                                            EMISSION', &
     &          '     TOTAL RAD      SINGLE        TOTAL      ',        &
     &          '  DIRECT                                TRANS',        &
     &          ' (CM-1) (MICRON)            (CM-1)  (MICRN)   (CM-1)', &
     &          '   (CM-1)  (MICRN)   (CM-1)   (CM-1)  (MICRN)',        &
     &          '   (CM-1)   (CM-1)  (MICRN)'
          ENDIF
      ELSEIF(NLOS.EQ.1)THEN
          WRITE(IPR,'(/52X,A,/4(/3A))')                                 &
     &      'RADIANCE(WATTS/CM2-STER-XXX)',                             &
     &      '  FREQ  WAVELENGTH   DIREC          PATH_THERMAL     ',    &
     &      ' SURFACE   PATH_SCAT_SOLAR   GROUND_REFLECTED',            &
     &      '     TOTAL_RADIANCE   INTEGRAL   TOTAL',                   &
     &      '                      EMIS       TOTAL      SCATTERED',    &
     &      ' EMISSION    TOTAL   SINGLE    TOTAL   DIRECT',            &
     &      '                                 TRANS',                   &
     &      ' (CM-1) (MICRON)             (CM-1)  (MICRN)   (CM-1)',    &
     &      '   (CM-1)   (CM-1)   (CM-1)   (CM-1)   (CM-1)',            &
     &      '    (CM-1)   (MICRN)'
      ELSE
          WRITE(IPR,'(/52X,A,/4(/3A))')                                 &
     &      'RADIANCE(WATTS/CM2-STER-XXX)',                             &
     &      '  FREQ  WAVELENGTH    LOS DIREC          PATH_THERMA',     &
     &      'L      SURFACE   PATH_SCAT_SOLAR   GROUND_REFLECTED',      &
     &      '     TOTAL_RADIANCE   INTEGRAL   TOTAL',                   &
     &      '                           EMIS       TOTAL      SCA',     &
     &      'TTERED EMISSION    TOTAL   SINGLE    TOTAL   DIRECT',      &
     &      '                                 TRANS',                   &
     &      ' (CM-1) (MICRON)                  (CM-1)  (MICRN)   ',     &
     &      '(CM-1)   (CM-1)   (CM-1)   (CM-1)   (CM-1)   (CM-1)',      &
     &      '    (CM-1)   (MICRN)'
      ENDIF

!     RETURN TO TRANS:
      RETURN
      END

      SUBROUTINE TWRITE(NEXTXY,NOPRNT,NREC,                             &
     &  VCEN,ALAM,UNIF,TRACE,SUMA,LN_TRN,TX)

!     WRITE SPECTRAL TRANSMITTANCE DATA:
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     ARGUMENTS:
!       NEXTXY   NUMBER OF ACTIVE SPECIES.
!       NOPRNT   PRINT FLAG.
!       NREC     NUMBER OF RECORDS IN BINARY FILE.
!       VCEN     SPECTRAL BAND CENTRAL FREQUENCY [CM-1].
!       ALAM     SPECTRAL WAVELENGTH [MICRONS].
!       UNIF     UNIFORMLY MIXED GASES SPECTRAL TRANSMITTANCE.
!       TRACE    TRACE GAS SPECTRAL TRANSMITTANCE.
!       SUMA     CUMULATIVE ABSORPTION SPECTRAL INTEGRAL.
!       LN_TRN   LOGARITHM OF COMBINED SPECIES TRANSMITTANCE.
      INTEGER NEXTXY,NOPRNT,NREC
      REAL VCEN,ALAM,UNIF,TRACE,SUMA,LN_TRN,TX(MEXTXY)

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /SPCCTL/
!       DODRGD   .TRUE. FOR NON-DEFAULT SPECTRAL ANALYSIS.
      LOGICAL DODGRD
      COMMON/SPCCTL/DODGRD

!     LOCAL VARIABLES:
!       IEXT     EXTINCTION SPECIES INDEX.
      INTEGER IEXT

!     WRITE SPECTRAL DATA:
      IF(BINOUT)THEN
          NREC=NREC+1
          WRITE(JPU)VCEN,TX(9),TX(17),UNIF,TX(31),TRACE,TX(4),TX(5),    &
     &      TX(6),TX(7),TX(11),TX(10),LN_TRN,TX(36),TX(44),TX(46),      &
     &      TX(47),TX(50),TX(52),TX(54),TX(64),TX(65),TX(16),           &
     &      (TX(IEXT),IEXT=MEXT+1,NEXTXY)
      ELSE
          IF(NOPRNT.LE.1)WRITE(IPR,'(F8.2,1X,G12.7,F6.4,F7.4,9F9.4,     &
     &      F12.3)')VCEN,ALAM,TX(9),TX(17),UNIF,TX(31),TRACE,TX(4),     &
     &      TX(5),TX(6),TX(7),TX(11),TX(10),SUMA
          WRITE(IPU,'(F8.2,11(1X,F6.4),1X,F7.3,99(1X,F6.4))')VCEN,TX(9),&
     &      TX(17),UNIF,TX(31),TRACE,TX(4),TX(5),TX(6),TX(7),TX(11),    &
     &      TX(10),LN_TRN,TX(36),TX(44),TX(46),TX(47),TX(50),TX(52),    &
     &      TX(54),TX(64),TX(65),TX(16),(TX(IEXT),IEXT=MEXT+1,NEXTXY)
      ENDIF
      IF(DODGRD)WRITE(JPUSCR)                                           &
     &  VCEN,TX(9),TX(17),UNIF,TX(31),TRACE,TX(4),TX(5),TX(6),TX(7),    &
     &  TX(11),TX(10),LN_TRN,TX(36),TX(44),TX(46),TX(47),TX(50),TX(52), &
     &  TX(54),TX(64),TX(65),TX(16),(TX(IEXT),IEXT=MEXT+1,NEXTXY)

!     RETURN TO TRANS:
      RETURN
      END

      SUBROUTINE EWRITE(DIS,IMULT,NOPRNT,NREC,VCEN,ALAM,NLOS,ILOS,DEMIS,&
     &  CONVRT,THMCUM,SUMTMS,BBG,RFSURF,SUMT,RADSUM,TRNLOS,LN_TRN,IFWHM)

!     WRITE PATH, SURFACE AND SCATTERED SPECTRAL EMISSION DATA:
      IMPLICIT NONE

!     ARGUMENTS:
!       DIS      LOGICAL FLAG, TRUE FOR DISORT MULTIPLE SCATTERING.
!       IMULT    MULTIPLE SCATTERING FLAG
!                  (0=NONE, 1=AT SENSOR, -1=AT FINAL OR TANGENT POINT).
!       NOPRNT   PRINT FLAG.
!       NREC     NUMBER OF RECORDS IN BINARY FILE.
!       VCEN     SPECTRAL BAND CENTRAL FREQUENCY [CM-1].
!       NLOS     NUMBER OF LINE-OF-SIGHT PATHS.
!       ILOS     LOOP INDEX FOR LINE-OF-SIGHT PATHS.
!       ALAM     SPECTRAL WAVELENGTH [MICRONS].
!       DEMIS    DIRECTIONAL EMISSIVITY CONVOLVED WITH SLIT FUNCTION.
!       CONVRT   CONVERSION FROM (1 / CM-1) TO (1 / MICRON).
!       THMCUM   PATH & SCATTERED THERMAL EMISSION [W SR-1 CM-2 / CM-1].
!       SUMTMS   SCATTERED THERMAL EMISSION [W SR-1 CM-2 / CM-1].
!       BBG      TRANSMITTED GROUND EMISSION [W SR-1 CM-2 / CM-1].
!       RFSURF   TRANSMITTED SURFACE REFLECTANCE [W SR-1 CM-2 / CM-1].
!       SUMT     TOTAL SPECTRAL RADIANCE AT SENSOR [W SR-1 CM-2 / CM-1].
!       RADSUM   SPECTRALLY INTEGRATED RADIANCE AT SENSOR [W SR-1 CM-2].
!       TRNLOS   LINE-OF-SIGHT COMBINED SPECIES TRANSMITTANCE [TX(9)].
!       LN_TRN   LOGARITHM OF COMBINED SPECIES PATH TRANSMITTANCE.
!       IFWHM    TRIANGULAR SLIT FULL-WIDTH-HALF-MAXIMUM [CM-1].
      LOGICAL DIS
      INTEGER IMULT,NOPRNT,NREC,NLOS,ILOS
      REAL VCEN,ALAM,DEMIS,CONVRT,THMCUM,SUMTMS,                        &
     &  BBG,RFSURF,SUMT,RADSUM,TRNLOS,LN_TRN,IFWHM

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /SPCCTL/
!       DODRGD   .TRUE. FOR NON-DEFAULT SPECTRAL ANALYSIS.
      LOGICAL DODGRD
      COMMON/SPCCTL/DODGRD

!     LOCAL VARIABLES:
!       M2DIFF   CHANNEL SECOND SPECTRAL MOMENT DIFFERENCE INTEGRAL.
      REAL M2DIFF

!     FUNCTIONS:
!       BT_wn    BRIGHTNESS TEMPERATURE FOR [W CM-1 / SR] RADIANCE.
      REAL BT_WN

!     WRITE SPECTRAL DATA:
      M2DIFF=(IFWHM/VCEN)**2/12
      IF(DIS)THEN

!         DISORT RUN, SCATTERED EMISSION NOT SEPARATED OUT:
          IF(BINOUT)THEN
              NREC=NREC+1
              WRITE(JPU)VCEN,TRNLOS,THMCUM,BBG,RFSURF,SUMT,             &
     &          LN_TRN,DEMIS,BT_WN(SUMT,VCEN,M2DIFF)
          ELSEIF(NLOS.EQ.1)THEN
              IF(NOPRNT.LE.1)WRITE(IPR,                                 &
     &          '(F8.2,1X,G12.7,F5.3,1P,2E10.2,10X,7E10.2,0P,F8.5)')    &
     &          VCEN,ALAM,DEMIS,THMCUM,CONVRT*THMCUM,BBG,CONVRT*BBG,    &
     &          RFSURF,CONVRT*RFSURF,SUMT,CONVRT*SUMT,RADSUM,TRNLOS
              WRITE(IPU,'(F8.2,F11.8,2(1P,E11.4,11X),2E22.4,            &
     &          0P,F26.3,F7.4,F22.3)')VCEN,TRNLOS,THMCUM,BBG,           &
     &          RFSURF,SUMT,LN_TRN,DEMIS,BT_WN(SUMT,VCEN,M2DIFF)
          ELSEIF(ILOS.EQ.1)THEN
              IF(NOPRNT.LE.1)WRITE(IPR,'(/F8.2,1X,G12.7,I4,F6.3,        &
     &          1P,2E10.2,10X,7E10.2,0P,F8.5)')VCEN,ALAM,ILOS,          &
     &          DEMIS,THMCUM,CONVRT*THMCUM,BBG,CONVRT*BBG,RFSURF,       &
     &          CONVRT*RFSURF,SUMT,CONVRT*SUMT,RADSUM,TRNLOS
              WRITE(IPU,'(/F8.2,I4,F11.8,1P,2(E11.4,11X),2E22.4,        &
     &          0P,F26.3,F7.4,F22.3)')VCEN,ILOS,TRNLOS,THMCUM,BBG,      &
     &          RFSURF,SUMT,LN_TRN,DEMIS,BT_WN(SUMT,VCEN,M2DIFF)
          ELSE
              IF(NOPRNT.LE.1)                                           &
     &          WRITE(IPR,'(21X,I4,F6.3,1P,2E10.2,10X,7E10.2,0P,F8.5)') &
     &          ILOS,DEMIS,THMCUM,CONVRT*THMCUM,BBG,CONVRT*BBG,         &
     &          RFSURF,CONVRT*RFSURF,SUMT,CONVRT*SUMT,RADSUM,TRNLOS
              WRITE(IPU,'(I12,F11.8,2(1P,E11.4,11X),2E22.4,             &
     &          0P,F26.3,F7.4,F22.3)')ILOS,TRNLOS,THMCUM,BBG,           &
     &          RFSURF,SUMT,LN_TRN,DEMIS,BT_WN(SUMT,VCEN,M2DIFF)
          ENDIF
          IF(DODGRD)                                                    &
     &      WRITE(JPUSCR)VCEN,TRNLOS,THMCUM,BBG,RFSURF,SUMT,LN_TRN,DEMIS

!         RETURN TO TRANS:
          RETURN
      ENDIF
      IF(BINOUT)THEN

!         BINARY OUTPUT WITH NO DISORT.
          NREC=NREC+1
          WRITE(JPU)VCEN,TRNLOS,THMCUM,SUMTMS,BBG,                      &
     &      RFSURF,SUMT,LN_TRN,DEMIS,BT_WN(SUMT,VCEN,M2DIFF)
      ELSEIF(NLOS.EQ.1)THEN

!         NO MULTIPLE SCATTERING:
          IF(NOPRNT.LE.1)WRITE(IPR,                                     &
     &      '(F8.2,1X,G12.7,F5.3,1P,6E10.2,2(E9.2,E10.2),0P,F8.5)')     &
     &      VCEN,ALAM,DEMIS,THMCUM,CONVRT*THMCUM,SUMTMS,BBG,CONVRT*BBG, &
     &      RFSURF,CONVRT*RFSURF,SUMT,CONVRT*SUMT,RADSUM,TRNLOS
          WRITE(IPU,'(F8.2,F11.8,1P,3E11.4,11X,2E22.4,                  &
     &      0P,F26.3,F7.4,F22.3)')VCEN,TRNLOS,THMCUM,SUMTMS,BBG,        &
     &      RFSURF,SUMT,LN_TRN,DEMIS,BT_WN(SUMT,VCEN,M2DIFF)
      ELSEIF(ILOS.EQ.1)THEN

!         ISAACS RUN, SCATTERED EMISSION SEPARATED OUT (FIRST LOS):
          IF(NOPRNT.LE.1)WRITE(IPR,                                     &
     &      '(/F8.2,1X,G12.7,I4,F6.3,1P,6E10.2,2(E9.2,E10.2),0P,F8.5)') &
     &      VCEN,ALAM,ILOS,DEMIS,THMCUM,CONVRT*THMCUM,SUMTMS,BBG,CONVRT &
     &      *BBG,RFSURF,CONVRT*RFSURF,SUMT,CONVRT*SUMT,RADSUM,TRNLOS
          WRITE(IPU,'(/F8.2,I4,F11.8,1P,3E11.4,11X,2E22.4,              &
     &      0P,F26.3,F7.4,F22.3)')VCEN,ILOS,TRNLOS,THMCUM,SUMTMS,BBG,   &
     &      RFSURF,SUMT,LN_TRN,DEMIS,BT_WN(SUMT,VCEN,M2DIFF)
      ELSE

!         ISAACS RUN, SCATTERED EMISSION SEPARATED OUT (MULTIPLE LOS):
          IF(NOPRNT.LE.1)WRITE(IPR,'(I25,F6.3,1P,6E10.2,2(E9.2,E10.2),  &
     &      0P,F8.5)')ILOS,DEMIS,THMCUM,CONVRT*THMCUM,SUMTMS,BBG,CONVRT &
     &      *BBG,RFSURF,CONVRT*RFSURF,SUMT,CONVRT*SUMT,RADSUM,TRNLOS
          WRITE(IPU,'(I12,F11.8,1P,3E11.4,11X,2E22.4,                   &
     &      0P,F26.3,F7.4,F22.3)')ILOS,TRNLOS,THMCUM,SUMTMS,BBG,        &
     &      RFSURF,SUMT,LN_TRN,DEMIS,BT_WN(SUMT,VCEN,M2DIFF)
      ENDIF

!     SPECTRAL DEGRADE POST-PROCESSING?
      IF(DODGRD)WRITE(JPUSCR)                                           &
     &  VCEN,TRNLOS,THMCUM,SUMTMS,BBG,RFSURF,SUMT,LN_TRN,DEMIS

!     RETURN TO TRANS:
      RETURN
      END

      SUBROUTINE RWRITE(DIS,IMULT,NOPRNT,NREC,VCEN,ALAM,NLOS,ILOS,      &
     &  DEMIS,CONVRT,THMCUM,SUMTMS,BBG,S0,SUMSSR,SUMSSS,RFSURF,         &
     &  RFLSS,SUMT,TSNREF,TSNOBS,RADSUM,TRNLOS,LN_TRN,IFWHM)

!     WRITE PATH, SURFACE & SCATTERED SPECTRAL EMISSION AND SOLAR DATA:
      IMPLICIT NONE

!     ARGUMENTS:
!       DIS      LOGICAL FLAG, TRUE FOR DISORT MULTIPLE SCATTERING.
!       IMULT    MULTIPLE SCATTERING FLAG
!                  (0=NONE, 1=AT SENSOR, -1=AT FINAL OR TANGENT POINT).
!       NOPRNT   PRINT FLAG.
!       NREC     NUMBER OF RECORDS IN BINARY FILE.
!       VCEN     SPECTRAL BAND CENTRAL FREQUENCY [CM-1].
!       ALAM     SPECTRAL WAVELENGTH [MICRONS].
!       NLOS     NUMBER OF LINE-OF-SIGHT PATHS.
!       ILOS     LOOP INDEX FOR LINE-OF-SIGHT PATHS.
!       DEMIS    DIRECTIONAL EMISSIVITY CONVOLVED WITH SLIT FUNCTION.
!       CONVRT   CONVERSION FROM (1 / CM-1) TO (1 / MICRON).
!       THMCUM   PATH & SCATTERED THERMAL EMISSION [W SR-1 CM-2 / CM-1].
!       SUMTMS   SCATTERED THERMAL EMISSION [W SR-1 CM-2 / CM-1].
!       BBG      TRANSMITTED GROUND EMISSION [W SR-1 CM-2 / CM-1].
!       S0       TOP-OF-ATMOSPHERE SOLAR IRRADIANCE [W CM-2 / CM-1].
!       SUMSSR   PATH SCATTERED SOLAR [W SR-1 CM-2 / CM-1].
!       SUMSSS   PATH SINGLE SCATTERED SOLAR [W SR-1 CM-2 / CM-1].
!       RFSURF   TRANSMITTED SURFACE REFLECTANCE [W SR-1 CM-2 / CM-1].
!       RFLSS    SUN-SURFACE-SENSOR DIRECT RADIANCE[W SR-1 CM-2 / CM-1].
!       SUMT     TOTAL SPECTRAL RADIANCE AT SENSOR [W SR-1 CM-2 / CM-1].
!       RADSUM   SPECTRALLY INTEGRATED RADIANCE AT SENSOR [W SR-1 CM-2].
!       TSNREF   SENSOR-FINAL_ALTITUDE-SUN TRANSMITTANCE CONVOLVED WITH
!                THE TOP-OF-ATMOSPHERE SOLAR IRRADIANCE [W CM-2 / CM-1].
!       TSNOBS   TRANSMITTED SOLAR IRRADIANCE AT SENSOR [W CM-2 / CM-1].
!       TRNLOS   LINE-OF-SIGHT COMBINED SPECIES TRANSMITTANCE [TX(9)].
!       LN_TRN   LOGARITHM OF COMBINED SPECIES PATH TRANSMITTANCE.
!       IFWHM    TRIANGULAR SLIT FULL-WIDTH-HALF-MAXIMUM [CM-1].
      LOGICAL DIS
      INTEGER IMULT,NOPRNT,NREC,NLOS,ILOS
      REAL VCEN,ALAM,DEMIS,CONVRT,THMCUM,SUMTMS,BBG,S0,SUMSSR,SUMSSS,   &
     &  RFSURF,RFLSS,SUMT,TSNREF,TSNOBS,RADSUM,TRNLOS,LN_TRN,IFWHM

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /SPCCTL/
!       DODRGD   .TRUE. FOR NON-DEFAULT SPECTRAL ANALYSIS.
      LOGICAL DODGRD
      COMMON/SPCCTL/DODGRD

!     LOCAL VARIABLES:
!       M2DIFF   CHANNEL SECOND SPECTRAL MOMENT DIFFERENCE INTEGRAL.
      REAL M2DIFF

!     FUNCTIONS:
!       BT_wn    BRIGHTNESS TEMPERATURE FOR [W CM-1 / SR] RADIANCE.
      REAL BT_WN

!     TRIANGULAR CHANNEL SECOND SPECTRAL MOMENT DIFFERENCE INTEGRAL:
      M2DIFF=(IFWHM/VCEN)**2/12

!     WRITE SPECTRAL DATA:
      IF(DIS)THEN

!         DISORT MULTIPLE SCATTERING (NO SUMTMS):
          IF(BINOUT)THEN
              NREC=NREC+1
              WRITE(JPU)VCEN,TRNLOS,THMCUM,BBG,SUMSSR,                  &
     &          SUMSSS,RFSURF,RFLSS,SUMT,TSNREF,TSNOBS,                 &
     &          LN_TRN,DEMIS,S0,BT_WN(SUMT,VCEN,M2DIFF)
          ELSEIF(NLOS.EQ.1)THEN
              IF(NOPRNT.LE.1)WRITE(IPR,'(F8.2,1X,G12.7,1P,E8.2,10E9.2,  &
     &          E10.2,0P,F8.5)')VCEN,ALAM,THMCUM,CONVRT*THMCUM,BBG,     &
     &          SUMSSR,CONVRT*SUMSSR,SUMSSS,RFSURF,CONVRT*RFSURF,       &
     &          RFLSS,SUMT,CONVRT*SUMT,RADSUM,TRNLOS
              WRITE(IPU,'(F8.2,F11.8,1P,E11.4,11X,6E11.4,2E9.2,         &
     &          0P,F8.3,F7.4,1P,E11.4,0P,F11.3)')VCEN,TRNLOS,           &
     &          THMCUM,BBG,SUMSSR,SUMSSS,RFSURF,RFLSS,SUMT,TSNREF,      &
     &          TSNOBS,LN_TRN,DEMIS,S0,BT_WN(SUMT,VCEN,M2DIFF)
          ELSEIF(ILOS.EQ.1)THEN
              IF(NOPRNT.LE.1)WRITE(IPR,'(/F8.2,1X,G12.7,I4,1P,          &
     &          11E9.2,E10.2,0P,F8.5)')VCEN,ALAM,ILOS,THMCUM,           &
     &          CONVRT*THMCUM,BBG,SUMSSR,CONVRT*SUMSSR,SUMSSS,RFSURF,   &
     &          CONVRT*RFSURF,RFLSS,SUMT,CONVRT*SUMT,RADSUM,TRNLOS
              WRITE(IPU,'(/F8.2,I4,F11.8,1P,E11.4,11X,6E11.4,2E9.2,     &
     &          0P,F8.3,F7.4,1P,E11.4,0P,F11.3)')VCEN,ILOS,TRNLOS,      &
     &          THMCUM,BBG,SUMSSR,SUMSSS,RFSURF,RFLSS,SUMT,TSNREF,      &
     &          TSNOBS,LN_TRN,DEMIS,S0,BT_WN(SUMT,VCEN,M2DIFF)
          ELSE
              IF(NOPRNT.LE.1)                                           &
     &          WRITE(IPR,'(I25,1P,11E9.2,E10.2,0P,F8.5)')ILOS,THMCUM,  &
     &          CONVRT*THMCUM,BBG,SUMSSR,CONVRT*SUMSSR,SUMSSS,RFSURF,   &
     &          CONVRT*RFSURF,RFLSS,SUMT,CONVRT*SUMT,RADSUM,TRNLOS
              WRITE(IPU,'(I12,F11.8,1P,E11.4,11X,6E11.4,2E9.2,          &
     &          0P,F8.3,F7.4,1P,E11.4,0P,F11.3)')ILOS,TRNLOS,           &
     &          THMCUM,BBG,SUMSSR,SUMSSS,RFSURF,RFLSS,SUMT,TSNREF,      &
     &          TSNOBS,LN_TRN,DEMIS,S0,BT_WN(SUMT,VCEN,M2DIFF)
          ENDIF
          IF(DODGRD)WRITE(JPUSCR)VCEN,TRNLOS,THMCUM,BBG,SUMSSR,         &
     &      SUMSSS,RFSURF,RFLSS,SUMT,TSNREF,TSNOBS,LN_TRN,DEMIS,S0

!         RETURN TO TRANS:
          RETURN
      ENDIF
      IF(BINOUT)THEN

!         BINARY OUTPUT:
          NREC=NREC+1
          WRITE(JPU)VCEN,TRNLOS,THMCUM,SUMTMS,BBG,SUMSSR,               &
     &      SUMSSS,RFSURF,RFLSS,SUMT,TSNREF,TSNOBS,                     &
     &      LN_TRN,DEMIS,S0,BT_WN(SUMT,VCEN,M2DIFF)
      ELSEIF(IMULT.EQ.0)THEN

!         NO MULTIPLE SCATTERING:
          WRITE(IPU,'(F8.2,F11.8,1P,8E11.4,2E9.2,0P,F8.3,F7.4,          &
     &      1P,E11.4,0P,F11.3)')VCEN,TRNLOS,THMCUM,SUMTMS,              &
     &      BBG,SUMSSR,SUMSSS,RFSURF,RFLSS,SUMT,TSNREF,                 &
     &      TSNOBS,LN_TRN,DEMIS,S0,BT_WN(SUMT,VCEN,M2DIFF)
          IF(NOPRNT.LE.1)WRITE(IPR,                                     &
     &      '(F8.2,1X,G12.7,1P,12E9.2,0P,F9.5)')VCEN,ALAM,THMCUM,       &
     &      CONVRT*THMCUM,BBG,CONVRT*BBG,S0,SUMSSR,CONVRT*SUMSSR,       &
     &      RFSURF,CONVRT*RFSURF,SUMT,CONVRT*SUMT,RADSUM,TRNLOS
      ELSEIF(NLOS.EQ.1)THEN

!         SINGLE LINE-OF-SIGHT:
          WRITE(IPU,'(F8.2,F11.8,1P,8E11.4,2E9.2,0P,F8.3,F7.4,          &
     &      1P,E11.4,0P,F11.3)')VCEN,TRNLOS,THMCUM,SUMTMS,              &
     &      BBG,SUMSSR,SUMSSS,RFSURF,RFLSS,SUMT,TSNREF,                 &
     &      TSNOBS,LN_TRN,DEMIS,S0,BT_WN(SUMT,VCEN,M2DIFF)
          IF(NOPRNT.LE.1)WRITE(IPR,'(F8.2,1X,G12.7,F5.3,1P,8E9.2,3E10.3,&
     &      0P,F8.5)')VCEN,ALAM,DEMIS,THMCUM,CONVRT*THMCUM,SUMTMS,BBG,  &
     &      SUMSSR,SUMSSS,RFSURF,RFLSS,SUMT,CONVRT*SUMT,RADSUM,TRNLOS
      ELSEIF(ILOS.EQ.1)THEN

!         ISAACS MULTIPLE SCATTERING (FIRST LINE-OF-SIGHT):
          WRITE(IPU,'(/F8.2,I4,F11.8,1P,8E11.4,2E9.2,0P,F8.3,F7.4,      &
     &      1P,E11.4,0P,F11.3)')VCEN,ILOS,TRNLOS,THMCUM,SUMTMS,         &
     &      BBG,SUMSSR,SUMSSS,RFSURF,RFLSS,SUMT,TSNREF,                 &
     &      TSNOBS,LN_TRN,DEMIS,S0,BT_WN(SUMT,VCEN,M2DIFF)
          IF(NOPRNT.LE.1)WRITE(IPR,                                     &
     &      '(/F8.2,1X,G12.7,I4,F6.3,1P,8E9.2,3E10.3,0P,F8.5)')         &
     &      VCEN,ALAM,ILOS,DEMIS,THMCUM,CONVRT*THMCUM,SUMTMS,BBG,       &
     &      SUMSSR,SUMSSS,RFSURF,RFLSS,SUMT,CONVRT*SUMT,RADSUM,TRNLOS
      ELSE

!         ISAACS MULTIPLE SCATTERING (SECOND LINE-OF-SIGHT):
          WRITE(IPU,'(I12,F11.8,1P,8E11.4,2E9.2,0P,F8.3,F7.4,           &
     &      1P,E11.4,0P,F11.3)')ILOS,TRNLOS,THMCUM,SUMTMS,              &
     &      BBG,SUMSSR,SUMSSS,RFSURF,RFLSS,SUMT,TSNREF,                 &
     &      TSNOBS,LN_TRN,DEMIS,S0,BT_WN(SUMT,VCEN,M2DIFF)
          IF(NOPRNT.LE.1)WRITE(IPR,'(I25,F6.3,1P,8E9.2,3E10.3,0P,F8.5)')&
     &      ILOS,DEMIS,THMCUM,CONVRT*THMCUM,SUMTMS,BBG,SUMSSR,          &
     &      SUMSSS,RFSURF,RFLSS,SUMT,CONVRT*SUMT,RADSUM,TRNLOS
      ENDIF

!     SPECTRAL DEGRADE POST-PROCESSING?
      IF(DODGRD)WRITE(JPUSCR)VCEN,TRNLOS,THMCUM,SUMTMS,BBG,SUMSSR,      &
     &  SUMSSS,RFSURF,RFLSS,SUMT,TSNREF,TSNOBS,LN_TRN,DEMIS,S0

!     RETURN TO TRANS:
      RETURN
      END

      SUBROUTINE SWRITE(NOPRNT,NREC,VCEN,ALAM,CONVRT,                   &
     &  TS0,S0,RADSUM,SSOL,TRNLOS,LN_TRN)

!     WRITE TRANSMITTED SOLAR IRRADIANCE SPECTRAL DATA:
      IMPLICIT NONE

!     ARGUMENTS:
!       NOPRNT   PRINT FLAG.
!       NREC     NUMBER OF RECORDS IN BINARY FILE.
!       VCEN     SPECTRAL BAND CENTRAL FREQUENCY [CM-1].
!       ALAM     SPECTRAL WAVELENGTH [MICRONS].
!       CONVRT   CONVERSION FROM (1 / CM-1) TO (1 / MICRON).
!       TS0      TRANSMITTED SOLAR IRRADIANCE [W CM-2 / CM-1].
!       S0       TOP-OF-ATMOSPHERE SOLAR IRRADIANCE [W CM-2 / CM-1].
!       RADSUM   SPECTRALLY INTEGRATED TRANSMITTED IRRADIANCE [W CM-2].
!       SSOL     INTEGRATED TOP-OF-ATMOSPHERE SOLAR IRRADIANCE [W CM-2].
!       TRNLOS   LINE-OF-SIGHT COMBINED SPECIES TRANSMITTANCE [TX(9)].
!       LN_TRN   LOGARITHM OF COMBINED SPECIES PATH TRANSMITTANCE.
      INTEGER NOPRNT,NREC
      REAL VCEN,ALAM,CONVRT,TS0,S0,RADSUM,SSOL,TRNLOS,LN_TRN

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /SPCCTL/
!       DODRGD   .TRUE. FOR NON-DEFAULT SPECTRAL ANALYSIS.
      LOGICAL DODGRD
      COMMON/SPCCTL/DODGRD

!     WRITE SPECTRAL DATA:
      IF(BINOUT)THEN
          NREC=NREC+1
          WRITE(JPU)VCEN,TRNLOS,TS0,S0,LN_TRN
      ELSE
          IF(NOPRNT.LE.1)                                               &
     &      WRITE(IPR,'(F8.2,1X,G12.7,1P,E8.2,5E10.2,0P,F9.4)')         &
     &      VCEN,ALAM,TS0,CONVRT*TS0,S0,CONVRT*S0,RADSUM,SSOL,TRNLOS
          WRITE(IPU,'(F8.2,F8.4,1P,2E9.2,T96,E10.3)')                   &
     &      VCEN,TRNLOS,TS0,S0,LN_TRN
      ENDIF
      IF(DODGRD)WRITE(JPUSCR)VCEN,TRNLOS,TS0,S0,LN_TRN

!     RETURN TO TRANS:
      RETURN
      END

      SUBROUTINE ENDWRT(NREC)

!     ENDWRT writes down total number of spectral data,
!     recorded in binary part of tape7
      IMPLICIT NONE
      INTEGER NREC
      INCLUDE 'IFIL.h'
      WRITE(IPU,'(A)')'TO RESTORE SPECTRAL DATA FROM'                   &
     &  //' THE BINARY FILE EXECUTE M5BINRESTORE.EXE'
      WRITE(IPU, '(A,I6)')                                              &
     &  'NUMBER OF SPEC.  RECORDS IN THE BINARY FILE :  ',NREC
      RETURN
      END
