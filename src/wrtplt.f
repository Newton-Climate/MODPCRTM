      SUBROUTINE WRTPLT(IPLOT,VCEN,ALAM,NLOS,YPLTMX,TRNLOS,CONVRT,SUMT)

!     WRITE <rootname>.plt SPECTRAL DATA:
      IMPLICIT NONE

!     INPUT ARGUMENTS:
!       IPLOT    UNIT NUMBER FOR <rootname>.plt FILE.
!       VCEN     SPECTRAL BAND CENTRAL FREQUENCY [CM-1].
!       ALAM     SPECTRAL WAVELENGTH [MICRONS].
!       NLOS     NUMBER OF LINE-OF-SIGHT PATHS.
!       TRNLOS   LINE-OF-SIGHT COMBINED SPECIES TRANSMITTANCE [TX(9)].
!       CONVRT   CONVERSION FROM (1 / CM-1) TO (1 / MICRON).
!       SUMT     TOTAL SPECTRAL RADIANCE AT SENSOR [W SR-1 CM-2 / CM-1].
      INTEGER IPLOT,NLOS
      REAL VCEN,ALAM,TRNLOS(NLOS),CONVRT,SUMT(NLOS)

!     INPUT/OUTPUT ARGUMENTS:
!       YPLTMX   MAXIMUM ORDINATE IN <rootname>.plt FILE.
      REAL YPLTMX(NLOS)

!     COMMONS:

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
!       FRMT     OUTPUT FORMAT.
!       ILOS     LOOP INDEX FOR LINE-OF-SIGHT PATHS.
!       CNVERT   LOCAL VALUE OF CONVRT.
      CHARACTER*22 FRMT
      INTEGER ILOS
      REAL CNVERT

!     BRANCH BASED ON OUTPUT FLAGS:
      IF(YFLAG.EQ.'T')THEN

!         TRANSMITTANCE OUTPUT:
          IF(XFLAG.EQ.'N')THEN
              FRMT='(0P,F16.4,   F16.8)   '
              WRITE(FRMT(11:13),'(I3)')NLOS
              WRITE(IPLOT,FRMT)1000*ALAM,(TRNLOS(ILOS),ILOS=1,NLOS)
          ELSEIF(XFLAG.EQ.'M')THEN
              FRMT='(0P,F16.7,   F16.8)   '
              WRITE(FRMT(11:13),'(I3)')NLOS
              WRITE(IPLOT,FRMT)ALAM,(TRNLOS(ILOS),ILOS=1,NLOS)
          ELSE
              FRMT='(0P,F16.2,   F16.8)   '
              WRITE(FRMT(11:13),'(I3)')NLOS
              WRITE(IPLOT,FRMT)VCEN,(TRNLOS(ILOS),ILOS=1,NLOS)
          ENDIF
          DO ILOS=1,NLOS
              IF(TRNLOS(ILOS).GT.YPLTMX(ILOS))YPLTMX(ILOS)=TRNLOS(ILOS)
          ENDDO
      ELSEIF(YFLAG.EQ.'R')THEN

!         TOTAL RADIANCE OR ATTENUATED SOLAR IRRADIANCE OUTPUT:
          IF(XFLAG.EQ.'N')THEN
              CNVERT=1000*CONVRT
              FRMT='(0P,F16.4,1P,   E16.5)'
              WRITE(FRMT(14:16),'(I3)')NLOS
              WRITE(IPLOT,FRMT)1000*ALAM,(CONVRT*SUMT(ILOS),ILOS=1,NLOS)
          ELSEIF(XFLAG.EQ.'M')THEN
              CNVERT=CONVRT
              FRMT='(0P,F16.7,1P,   E16.5)'
              WRITE(FRMT(14:16),'(I3)')NLOS
              WRITE(IPLOT,FRMT)ALAM,(CONVRT*SUMT(ILOS),ILOS=1,NLOS)
          ELSE
              CNVERT=1.
              FRMT='(0P,F16.2,1P,   E16.5)'
              WRITE(FRMT(14:16),'(I3)')NLOS
              WRITE(IPLOT,FRMT)VCEN,(SUMT(ILOS),ILOS=1,NLOS)
          ENDIF
          DO ILOS=1,NLOS
              IF(CNVERT*SUMT(ILOS).GT.YPLTMX(ILOS))                     &
     &          YPLTMX(ILOS)=CNVERT*SUMT(ILOS)
          ENDDO
      ENDIF

!     RETURN TO TRANS:
      RETURN
      END

      SUBROUTINE WPLTBN(JPLOT,VCEN,ALAM,NLOS,YPLTMX,TRNLOS,CONVRT,SUMT)

!     BINARY FORM Of WRTPLT.  WRITES TO BINARY FILE.  SHOULD BE CALLED
!     IF GLOBAL FLAG BINOUT IS .TRUE.  SEE COMMENTS IN WRTPLT.
      IMPLICIT NONE

!     INPUT ARGUMENTS:
!       JPLOT    UNIT NUMBER FOR <rootname>b.plt BINARY FILE.
!       VCEN     SPECTRAL BAND CENTRAL FREQUENCY [CM-1].
!       ALAM     SPECTRAL WAVELENGTH [MICRONS].
!       NLOS     NUMBER OF LINE-OF-SIGHT PATHS.
!       TRNLOS   LINE-OF-SIGHT COMBINED SPECIES TRANSMITTANCE [TX(9)].
!       CONVRT   CONVERSION FROM (1 / CM-1) TO (1 / MICRON).
!       SUMT     TOTAL SPECTRAL RADIANCE AT SENSOR [W SR-1 CM-2 / CM-1].
      INTEGER JPLOT,NLOS
      REAL VCEN,ALAM,TRNLOS(NLOS),CONVRT,SUMT(NLOS)

!     INPUT/OUTPUT ARGUMENTS:
!       YPLTMX   MAXIMUM ORDINATE IN <rootname>.plt FILE.
      REAL YPLTMX(NLOS)

!     COMMONS:
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
!       ILOS     LOOP INDEX FOR LINE-OF-SIGHT PATHS.
!       CNVERT   LOCAL VALUE OF CONVRT.
      INTEGER ILOS
      REAL CNVERT

!     BRANCH BASED ON OUTPUT FLAGS:
      IF(YFLAG.EQ.'T')THEN

!         TRANSMITTANCE OUTPUT:
          IF(XFLAG.EQ.'N')THEN
              WRITE(JPLOT)1000*ALAM,(TRNLOS(ILOS),ILOS=1,NLOS)
          ELSEIF(XFLAG.EQ.'M')THEN
              WRITE(JPLOT)ALAM,(TRNLOS(ILOS),ILOS=1,NLOS)
          ELSE
              WRITE(JPLOT)VCEN,(TRNLOS(ILOS),ILOS=1,NLOS)
          ENDIF
          DO ILOS=1,NLOS
              IF(TRNLOS(ILOS).GT.YPLTMX(ILOS))YPLTMX(ILOS)=TRNLOS(ILOS)
          ENDDO
      ELSEIF(YFLAG.EQ.'R')THEN

!         TOTAL RADIANCE OR ATTENUATED SOLAR IRRADIANCE OUTPUT:
          IF(XFLAG.EQ.'N')THEN
              CNVERT=1000*CONVRT
              WRITE(JPLOT)1000*ALAM,(CNVERT*SUMT(ILOS),ILOS=1,NLOS)
          ELSEIF(XFLAG.EQ.'M')THEN
!**************** VINCENT ROSS ADDED TO SOLVE BUG ************
              CNVERT=CONVRT
!*********************** END VINCENT ROSS ********************
              WRITE(JPLOT)ALAM,(CNVERT*SUMT(ILOS),ILOS=1,NLOS)
          ELSE
!**************** VINCENT ROSS ADDED TO SOLVE BUG ************
              CNVERT=1.
!*********************** END VINCENT ROSS ********************
              WRITE(JPLOT)VCEN,(SUMT(ILOS),ILOS=1,NLOS)
          ENDIF
          DO ILOS=1,NLOS
              IF(CNVERT*SUMT(ILOS).GT.YPLTMX(ILOS))                     &
     &          YPLTMX(ILOS)=CNVERT*SUMT(ILOS)
          ENDDO
      ENDIF

!     RETURN TO TRANS:
      RETURN
      END

      SUBROUTINE ENDPLT(IPLOT,IEMSCT,NLOS,YPLTMX)

!     WRITE DELIMITER FOR <rootname>.plt FILE.
      IMPLICIT NONE

!     INPUT ARGUMENTS:
!       IPLOT    UNIT NUMBER FOR <rootname>.plt FILE.
!       IEMSCT   RADIATIVE TRANSFER MODE
!                  0 FOR TRANSMITTANCE
!                  1 FOR THERMAL EMISSION
!                  2 FOR THERMAL EMISSION PLUS SOLAR SCATTER
!                  3 FOR TRANSMITTED SOLAR IRRADIANCE
!       NLOS     NUMBER OF LINE-OF-SIGHT PATHS.
!       YPLTMX   MAXIMUM ORDINATE IN <rootname>.plt FILE.
      INTEGER IPLOT,IEMSCT,NLOS
      REAL YPLTMX(NLOS)

!     COMMONS:

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
!       FRMT     OUTPUT FORMAT.
!       ILOS     LOOP INDEX FOR LINE-OF-SIGHT PATHS.
      CHARACTER*18 FRMT
      INTEGER ILOS

!     BRANCH BASED ON OUTPUT FLAGS:
      IF(YFLAG.EQ.'T')THEN

!         BLANK DELIMITER?
          IF(DLIMIT.EQ.'        ')THEN
              WRITE(IPLOT,*)

!             RETURN TO TRANS:
              RETURN
          ENDIF

!         FORMAT:
          FRMT='(2A,1P,   E16.5,A)'
          WRITE(FRMT(8:10),'(I3)')NLOS

!         TRANSMITTANCE OUTPUT:
          IF(XFLAG.EQ.'N')THEN
              WRITE(IPLOT,FRMT)DLIMIT,'    MAX:',(YPLTMX(ILOS),         &
     &          ILOS=1,NLOS),' (WAVELENGTHS IN NANOMETERS)'
          ELSEIF(XFLAG.EQ.'M')THEN
              WRITE(IPLOT,FRMT)DLIMIT,'    MAX:',(YPLTMX(ILOS),         &
     &          ILOS=1,NLOS),' (WAVELENGTHS IN MICRONS)'
          ELSE
              WRITE(IPLOT,FRMT)DLIMIT,'    MAX:',(YPLTMX(ILOS),         &
     &          ILOS=1,NLOS),' (FREQUENCIES IN CM-1)'
          ENDIF
      ELSEIF(YFLAG.EQ.'R')THEN

!         BLANK DELIMITER?
          IF(DLIMIT.EQ.'        ')THEN
              WRITE(IPLOT,*)

!             RETURN TO TRANS:
              RETURN
          ENDIF

!         FORMAT:
          FRMT='(2A,1P,   E16.5,A)'
          WRITE(FRMT(8:10),'(I3)')NLOS

!         IRRADIANCE OR RADIANCE?
          IF(IEMSCT.EQ.3)THEN

!             ATTENUATED SOLAR IRRADIANCE OUTPUT:
              IF(XFLAG.EQ.'N')THEN
                  WRITE(IPLOT,FRMT)DLIMIT,'    MAX:',                   &
     &              (YPLTMX(ILOS),ILOS=1,NLOS),                         &
     &              ' (WAVELENGTHS IN NANOMETERS; TRANSMITTED SOLAR'    &
     &              //' IRRADIANCE IN MICRO-WATTS CM-2 / NANOMETER)'
              ELSEIF(XFLAG.EQ.'M')THEN
                  WRITE(IPLOT,FRMT)DLIMIT,'    MAX:',                   &
     &              (YPLTMX(ILOS),ILOS=1,NLOS),                         &
     &              ' (WAVELENGTHS IN MICRONS; TRANSMITTED SOLAR'       &
     &              //' IRRADIANCE IN WATTS CM-2 / MICRON)'
              ELSE
                  WRITE(IPLOT,FRMT)DLIMIT,'    MAX:',                   &
     &              (YPLTMX(ILOS),ILOS=1,NLOS),                         &
     &              ' (FREQUENCY IN CM-1; TRANSMITTED SOLAR'            &
     &              //' IRRADIANCE IN WATTS CM-2 / CM-1)'
              ENDIF
          ELSE

!             TOTAL RADIANCE OUTPUT:
              IF(XFLAG.EQ.'N')THEN
                  WRITE(IPLOT,FRMT)DLIMIT,'    MAX:',                   &
     &              (YPLTMX(ILOS),ILOS=1,NLOS),                         &
     &              ' (WAVELENGTHS IN NANOMETERS; TOTAL SPECTRAL'       &
     &              //' RADIANCE IN MICRO-WATTS SR-1 CM-2 / NANOMETER)'
              ELSEIF(XFLAG.EQ.'M')THEN
                  WRITE(IPLOT,FRMT)DLIMIT,'    MAX:',                   &
     &              (YPLTMX(ILOS),ILOS=1,NLOS),                         &
     &              ' (WAVELENGTHS IN MICRONS; TOTAL SPECTRAL'          &
     &              //' RADIANCE IN WATTS SR-1 CM-2 / MICRON)'
              ELSE
                  WRITE(IPLOT,FRMT)DLIMIT,'    MAX:',                   &
     &              (YPLTMX(ILOS),ILOS=1,NLOS),                         &
     &              ' (FREQUENCY IN CM-1; TOTAL SPECTRAL'               &
     &              //' RADIANCE IN WATTS SR-1 CM-2 / CM-1)'
              ENDIF
          ENDIF
      ENDIF

!     RETURN TO TRANS
      RETURN
      END

      SUBROUTINE ENDPTB(JPLOT,IEMSCT,NLOS,YPLTMX)

!     WRITE DELIMITER FOR <rootname>.plt FILE (SAME AS ENDPLT BUT BINARY
!     WRITE; TO MARK A DELIMITER IT WRITES FIRST -1.0 -1.0 TO FILE).
      IMPLICIT NONE

!     INPUT ARGUMENTS:
!       JPLOT    UNIT NUMBER FOR <rootname>b.plt BINARY FILE.
!       IEMSCT   RADIATIVE TRANSFER MODE
!                  0 FOR TRANSMITTANCE
!                  1 FOR THERMAL EMISSION
!                  2 FOR THERMAL EMISSION PLUS SOLAR SCATTER
!                  3 FOR TRANSMITTED SOLAR IRRADIANCE
!       NLOS     NUMBER OF LINE-OF-SIGHT PATHS.
!       YPLTMX   MAXIMUM ORDINATE IN <rootname>.plt FILE.
      INTEGER JPLOT,IEMSCT,NLOS
      REAL YPLTMX(NLOS)

!     COMMONS:

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
      CHARACTER YFLAG*1,XFLAG*1,DLIMIT*8,FLAGS*7,BUF*128
      COMMON/CFLAGS/YFLAG,XFLAG,DLIMIT,FLAGS

!     LOCAL VARIABLES:
!       FRMT     OUTPUT FORMAT.
!       ILOS     LOOP INDEX FOR LINE-OF-SIGHT PATHS.
      CHARACTER*18 FRMT
      INTEGER ILOS

!     BRANCH BASED ON OUTPUT FLAGS:
      IF(YFLAG.EQ.'T')THEN

!         BLANK DELIMITER?
          IF(DLIMIT.EQ.'        ')THEN
              WRITE(JPLOT) -1.0, -1.0
              BUF=' '
              WRITE(JPLOT)BUF

!             RETURN TO TRANS:
              RETURN
          ENDIF

!         FORMAT:
          FRMT='(2A,1P,   E16.5,A)'
          WRITE(FRMT(8:10),'(I3)')NLOS

!         TRANSMITTANCE OUTPUT:
          IF(XFLAG.EQ.'N')THEN
              WRITE(BUF,FRMT)DLIMIT,'    MAX:',(YPLTMX(ILOS),           &
     &          ILOS=1,NLOS),' (WAVELENGTHS IN NANOMETERS)'
          ELSEIF(XFLAG.EQ.'M')THEN
              WRITE(BUF,FRMT)DLIMIT,'    MAX:',(YPLTMX(ILOS),           &
     &          ILOS=1,NLOS),' (WAVELENGTHS IN MICRONS)'
          ELSE
              WRITE(BUF,FRMT)DLIMIT,'    MAX:',(YPLTMX(ILOS),           &
     &          ILOS=1,NLOS),' (FREQUENCIES IN CM-1)'
          ENDIF
      ELSEIF(YFLAG.EQ.'R')THEN

!         BLANK DELIMITER?
          IF(DLIMIT.EQ.'        ')THEN
              WRITE(JPLOT) -1.0, -1.0
              BUF=' '
              WRITE(JPLOT)BUF

!             RETURN TO TRANS:
              RETURN
          ENDIF

!         FORMAT:
          FRMT='(2A,1P,   E16.5,A)'
          WRITE(FRMT(8:10),'(I3)')NLOS

!         IRRADIANCE OR RADIANCE?
          IF(IEMSCT.EQ.3)THEN

!             ATTENUATED SOLAR IRRADIANCE OUTPUT:
              IF(XFLAG.EQ.'N')THEN
                  WRITE(BUF,FRMT)DLIMIT,'    MAX:',                     &
     &              (YPLTMX(ILOS),ILOS=1,NLOS),                         &
     &              ' (WAVELENGTHS IN NANOMETERS; TRANSMITTED SOLAR'    &
     &              //' IRRADIANCE IN MICRO-WATTS CM-2 / NANOMETER)'
              ELSEIF(XFLAG.EQ.'M')THEN
                  WRITE(BUF,FRMT)DLIMIT,'    MAX:',                     &
     &              (YPLTMX(ILOS),ILOS=1,NLOS),                         &
     &              ' (WAVELENGTHS IN MICRONS; TRANSMITTED SOLAR'       &
     &              //' IRRADIANCE IN WATTS CM-2 / MICRON)'
              ELSE
                  WRITE(BUF,FRMT)DLIMIT,'    MAX:',                     &
     &              (YPLTMX(ILOS),ILOS=1,NLOS),                         &
     &              ' (FREQUENCY IN CM-1; TRANSMITTED SOLAR'            &
     &              //' IRRADIANCE IN WATTS CM-2 / CM-1)'
              ENDIF
          ELSE

!             TOTAL RADIANCE OUTPUT:
              IF(XFLAG.EQ.'N')THEN
                  WRITE(BUF,FRMT)DLIMIT,'    MAX:',                     &
     &              (YPLTMX(ILOS),ILOS=1,NLOS),                         &
     &              ' (WAVELENGTHS IN NANOMETERS; TOTAL SPECTRAL'       &
     &              //' RADIANCE IN MICRO-WATTS SR-1 CM-2 / NANOMETER)'
              ELSEIF(XFLAG.EQ.'M')THEN
                  WRITE(BUF,FRMT)DLIMIT,'    MAX:',                     &
     &              (YPLTMX(ILOS),ILOS=1,NLOS),                         &
     &              ' (WAVELENGTHS IN MICRONS; TOTAL SPECTRAL'          &
     &              //' RADIANCE IN WATTS SR-1 CM-2 / MICRON)'
              ELSE
                  WRITE(BUF,FRMT)DLIMIT,'    MAX:',                     &
     &              (YPLTMX(ILOS),ILOS=1,NLOS),                         &
     &              ' (FREQUENCY IN CM-1; TOTAL SPECTRAL'               &
     &              //' RADIANCE IN WATTS SR-1 CM-2 / CM-1)'
              ENDIF
          ENDIF
      ELSE
          RETURN
      ENDIF
      WRITE(JPLOT) -1.0, -1.0
      WRITE(JPLOT)BUF

!     RETURN TO TRANS
      RETURN
      END
