      SUBROUTINE AEREXT(V,MNAER)

!     SPECTRALLY INTERPOLATES AEROSOL EXTINCTION, ABSORPTION
!     AND ASYMMETRY (JAN 1986, A.E.R. INC.) COEFFICIENTS.
      IMPLICIT NONE

!     INPUT ARGUMENTS
!       V        SPECTRAL FREQUENCY [CM-1].
!       MNAER    MINIMUM AEROSOL INDEX TO BE USED
!                (=1 FOR ALL AEROSOLS, =6 FOR CLOUDS ONLY).
      REAL V
      INTEGER MNAER

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     COMMONS:
      INCLUDE 'BASE.h'
      INCLUDE 'IFIL.h'

!     /CARD2A/
      INTEGER NCRALT,NCRSPC
      DOUBLE PRECISION CTHIK,CALT
      REAL CEXT,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP
      COMMON/CARD2A/CTHIK,CALT,CEXT,NCRALT,NCRSPC,CWAVLN,CCOLWD,        &
     &  CCOLIP,CHUMID,ASYMWD,ASYMIP
      INTEGER IPH
      REAL G
      COMMON/CARD3A/IPH,G

!     /DISRT/
!       UMU      MONOTONICALLY INCREASING LIST OF DISTINCT USER-PATH
!                COSINE POLAR ANGLES.
!       PHI      MONOTONICALLY INCREASING LIST OF DISTINCT RELATIVE
!                SOLAR AZIMUTH ANGLES [0 TO 180 DEG].
!       NSTR     NUMBER OF DISCRETE ORDINATE STREAMS.
!       NAZ      NUMBER OF DISORT AZIMUTH COMPONENTS.
!       N2GAUS   ORDER OF DOUBLE-GAUSS QUADRATURES.
!       NUMU     NUMBER OF DISTINCT USER LINE-OF-SIGHT POLAR ANGLES.
!       MAPUMU   MAPPING FROM LINE-OF-SIGHT INDEX TO UMU ARRAY ENTRY.
!       NPHI     NUMBER OF DISTINCT RELATIVE SOLAR AZIMUTH ANGLES.
!       MAPPHI   MAPPING FROM LINE-OF-SIGHT INDEX TO PHI ARRAY ENTRY.
!       DIS      LOGICAL FLAG, TRUE FOR DISORT MULTIPLE SCATTERING.
!       DISAZM   LOGICAL FLAG, TRUE FOR DISORT WITH AZIMUTH DEPENDENCE.
!       DISALB   LOGICAL FLAG, TRUE FOR DISORT SPHERICAL ALBEDO OPTION.
!       LDISCL   LOGICAL FLAG, TRUE FOR ISAACS SCALED TO DISORT.
      REAL UMU,PHI
      INTEGER NSTR,NAZ,N2GAUS,NUMU,MAPUMU,NPHI,MAPPHI
      LOGICAL DIS,DISAZM,DISALB,LDISCL
      COMMON/DISRT/UMU(MXUMU),PHI(MAXPHI),NSTR,NAZ,N2GAUS,NUMU,         &
     &  MAPUMU(MLOS),NPHI,MAPPHI(MLOS),DIS,DISAZM,DISALB,LDISCL
      SAVE /DISRT/

!     /EXTWAV/
!       VX0      AEROSOL SPECTRAL GRID [CM-1].
      REAL VX0
      COMMON/EXTWAV/VX0(NWAVLN)

!     /AER/ THERE ARE "MAER=17" PARTICULATE COMPONENTS:
!           1       AEROSOL 1 (NOMINALLY, BOUNDARY LAYER AEROSOL).
!           2       AEROSOL 2 (NOMINALLY, TROPOSPHERIC AEROSOL).
!           3       AEROSOL 3 (NOMINALLY, STRATOSPHERIC AEROSOL).
!           4       AEROSOL 4 (NOMINALLY, VOLCANIC AEROSOL).
!           5       CIRRUS CLOUD.
!           6       CLOUD 1 (NOMINALLY, WATER CLOUD).
!           7       CLOUD 2 (NOMINALLY, ICE CLOUD).
!           8-17    NOVAM (NAVY OCEANIC VERTICAL AEROSOL MODEL) LAYERS.
!       NAER     NUMBER OF ACTIVE AEROSOLS.
!       EXTV     SPECTRAL EXTINCTION (NORMALIZED TO 1 AT 550 NM).
!       ABSV     SPECTRAL ABSORPTION (1-ABSV/EXTV=SCATTERING ALBEDO).
!       ASYV     SPECTRAL LEGENDRE MOMENT (DIVIDED BY 2N+1).
!       FRAC5    5 CM-1 GRID SPECTRAL INTERPOLATION FRACTION.
!       ASYVLO   ASYMMETRY FACTOR FROM PREVIOUS SPECTRAL FREQUENCY.
      INTEGER NAER
      REAL EXTV,ABSV,ASYV,FRAC5,ASYVLO
      COMMON/AER/NAER,EXTV(MAER),ABSV(MAER),ASYV(MXCMU,MAER),           &
     &  FRAC5,ASYVLO(MAER)

      REAL TAER
      COMMON/AERTM/TAER(MAER)

!     /COMNOV/
      LOGICAL LNOVAM
      REAL EXTNOV(MNOV,MXWVLN),ABSNOV(MNOV,MXWVLN),                     &
     &     ASMNOV(MNOV,MXWVLN),WLNOV(MXWVLN)
      INTEGER NNOV,NWLNOV
      COMMON/COMNOV/LNOVAM,EXTNOV,ABSNOV,ASMNOV,WLNOV,NNOV,NWLNOV

!     /USSPC/
      INTEGER NARSPC
      REAL VARSPC
      LOGICAL LARUSS
      COMMON/USSPC/NARSPC(4),VARSPC(4,NWAVLN),LARUSS
      INTEGER NUSS

!     /aeralb/
!       lssalb   index of final aerosol albedo spectral grid point.
!       acoalb   aerosol single scattering co-albedo scaling factor.
!       awavln   aerosol albedo spectral grid points [microns].
!       assalb   boundary layer and troposphere aerosol single
!                scattering albedos.
!       aastrm   wavelength interpolation exponents (equivalent to an
!                angstrom coefficient for single scattering albedo)
!       asclmn   aerosol (1=boundary layer, 2=tropospheric)
!                single scattering albedo scaling factors
!                for wavelengths below awavln(0).
!       asclmx   aerosol (1=boundary layer, 2=tropospheric)
!                single scattering albedo scaling factors
!                for wavelengths above awavln(lssalb).
      INTEGER LSSALB
      REAL ACOALB,AWAVLN,ASSALB,AASTRM,ASCLMN,ASCLMX
      COMMON/AERALB/LSSALB,ACOALB,AWAVLN(0:MSSALB),ASSALB(0:MSSALB),    &
     &     ASCLMN(2),AASTRM(MSSALB),ASCLMX(2)

!     /aerc/
!       wavlen   cloud/rain data wavelength grid [microns].
!       awccon   water column content.
!       extc     aerosol/cloud extinction coefficients.
!       absc     aerosol/cloud absorption coefficients.
!       asym     aerosol/cloud asymmetry factors.
!       astrom   aerosol angstrom law exponents.
!       astmo    aerosol angstrom law offset.
      REAL ASTROM,ASTMO
      COMMON/AERC/ASTROM(MAER,MXWVLN),ASTMO

!     /LSTART/
!       TOLABS   OUTPUT WARNING THRESHOLD FOR EXCESS ABSORPTION.
!       LLORMN   OUTPUT WARNING FLAG FOR SMALL LORENTZ HALF WIDTH.
!       LLORMX   OUTPUT WARNING FLAG FOR LARGE LORENTZ HALF WIDTH.
!       LDOPMN   OUTPUT WARNING FLAG FOR SMALL DOPPLER HALF WIDTH.
!       LDOPMX   OUTPUT WARNING FLAG FOR LARGE DOPPLER HALF WIDTH.
!       LLINMN   OUTPUT WARNING FLAG FOR SMALL EFFECTIVE LINE NUMBER.
!       LSSA     OUTPUT WARNING FLAG FOR SCATTERING ALBEDO NEAR 1.
!       LAEREX   OUTPUT WARNING FLAG FOR ROUTINE AEREXT.
!       LBMCRK   OUTPUT WARNING FLAG FOR ROUTINE BMCRKS.
!       LFLUXS   OUTPUT WARNING FLAG FOR ROUTINE FLUXES.
!       LSSRAD   OUTPUT WARNING FLAG FOR ROUTINE SSRAD.
!       LTRLAY   OUTPUT WARNING FLAG FOR ROUTINE TRLAY.
!       LNVERS   OUTPUT WARNING FLAG FOR ROUTINE DENFAC.
!       LGEOM    OUTPUT WARNING FLAG FOR INDEX OF REFRACTION GRADIENT.
!       LO3TRN   OUTPUT WARNING FLAG FOR OZONE CURTIS-GODSON PROBLEM.
      REAL TOLABS
      LOGICAL LLORMN,LLORMX,LDOPMN,LDOPMX,LLINMN,LSSA,                  &
     &  LAEREX,LBMCRK,LFLUXS,LSSRAD,LTRLAY,LNVERS,LGEOM,LO3TRN
      COMMON/LSTART/TOLABS,LLORMN,LLORMX,LDOPMN,LDOPMX,LLINMN,LSSA,     &
     &  LAEREX,LBMCRK,LFLUXS,LSSRAD,LTRLAY,LNVERS,LGEOM,LO3TRN
      SAVE /LSTART/

!     LOCAL VARIABLES FOR ANGSTROM COEFFICIENTS:
      REAL ALAM,FACTOR,GAMFOG,GMFOG,VMN,VMX
      INTEGER IAER,IWAV,IWAVM1,MXAER,NWAVM1,JAER
      INTEGER ISSALB
      REAL COALB,RATLAM

!     LOCAL VARIABLES:
!       ISTR     LOOP INDEX FOR MULTIPLE SCATTERING STREAMS.
!       ISTRM1   ISTR MINUS ONE.
      INTEGER ISTR,ISTRM1

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      SAVE /EXTWAV/
      EXTERNAL DEVCBD,EXTDT1,EXTDT2

!     NO EXTINCTION ASSUMED FOR V=0.
      MXAER=7
      IF(V.LE.1.E-5)THEN
          DO IAER=MNAER,MXAER
              EXTV(IAER)=0.
              ABSV(IAER)=0.
              DO ISTR=1,NSTR
                  ASYV(ISTR,IAER)=0.
              ENDDO
          ENDDO
          DO IAER=8,17
              EXTV(IAER)=0.
              ABSV(IAER)=0.
              DO ISTR=1,NSTR
                  ASYV(ISTR,IAER)=0.
              ENDDO
          ENDDO
          IF(IPH.EQ.0)GOTO 90
          RETURN
      ELSEIF(.NOT.LNOVAM)THEN

!         FOR REPEAT RUNS NEED TO ZERO OUT PREVIOUS VALUES
          DO IAER=8,17
              EXTV(IAER)=0.
              ABSV(IAER)=0.
              DO ISTR=1,NSTR
                  ASYV(ISTR,IAER)=0.
              ENDDO
          ENDDO
      ENDIF

!     CHECK FOR USER-DEFINED CLOUD SPECTRAL DATA
      ALAM=10000./V
      IF(NCRSPC.GE.2)THEN

!         USER-DEFINED CLOUD SPECTRAL DATA (IAER=6 AND IAER=7)
!         IS BASED ON "WAVLEN" WAVELENGTHS, NOT "VX0" WAVELENGTHS.
          MXAER=5
          IF(ALAM.LT.WAVLEN(1).OR. ALAM.GT.WAVLEN(NCRSPC))THEN
              IF(LAEREX)THEN
                  WRITE(IPR,'(A,F8.4,A,/10X,2(A,F8.4),A,                &
     &              /10X,A,F8.4,2A,/10X,2A,/)')                         &
     &              ' WARNING:  CLOUD SPECTRAL DATA IS REQUIRED AT',    &
     &              ALAM,' MICRONS, BUT DATA WAS ONLY',                 &
     &              ' INPUT FROM',WAVLEN(1),' TO', WAVLEN(NCRSPC),      &
     &              ' MICRONS.  END POINT DATA WILL BE USED',           &
     &              ' FOR',ALAM,' MICRONS.  THIS WARNING IS',           &
     &              ' NOT REPEATED IF ADDITIONAL DATA',                 &
     &              ' IS REQUIRED THAT FALLS OUTSIDE',                  &
     &              ' THE INPUT SPECTRAL RANGE.'
                  LAEREX=.FALSE.
              ENDIF
              IWAV=1
              IF(ALAM.GT.WAVLEN(NCRSPC))IWAV=NCRSPC
              EXTV(6)=EXTC(6,IWAV)
              ABSV(6)=ABSC(6,IWAV)
              ASYV(1,6)=ASYM(6,IWAV)
              EXTV(7)=EXTC(7,IWAV)
              ABSV(7)=ABSC(7,IWAV)
              ASYV(1,7)=ASYM(7,IWAV)
          ELSE
              IWAVM1=1
              DO IWAV=2,NCRSPC-1
                  IF(ALAM.LE.WAVLEN(IWAV))GOTO 10
                  IWAVM1=IWAV
              ENDDO
              IWAV=NCRSPC
   10         CONTINUE
              FACTOR=(WAVLEN(IWAV)-ALAM)/(WAVLEN(IWAV)-WAVLEN(IWAVM1))
              EXTV(6)=EXTC(6,IWAV)+FACTOR*(EXTC(6,IWAVM1)-EXTC(6,IWAV))
              ABSV(6)=ABSC(6,IWAV)+FACTOR*(ABSC(6,IWAVM1)-ABSC(6,IWAV))
              ASYV(1,6)                                                 &
     &          =ASYM(6,IWAV)+FACTOR*(ASYM(6,IWAVM1)-ASYM(6,IWAV))
              EXTV(7)=EXTC(7,IWAV)+FACTOR*(EXTC(7,IWAVM1)-EXTC(7,IWAV))
              ABSV(7)=ABSC(7,IWAV)+FACTOR*(ABSC(7,IWAVM1)-ABSC(7,IWAV))
              ASYV(1,7)                                                 &
     &          =ASYM(7,IWAV)+FACTOR*(ASYM(7,IWAVM1)-ASYM(7,IWAV))
          ENDIF
          ISTRM1=1
          DO ISTR=2,NSTR
              ASYV(ISTR,6)=ASYV(ISTRM1,6)*ASYV(1,6)
              ASYV(ISTR,7)=ASYV(ISTRM1,7)*ASYV(1,7)
              ISTRM1=ISTR
          ENDDO
          IF(MNAER.GE.6)RETURN
      ELSEIF(NCRSPC.EQ.1)THEN
          MXAER=5
          CALL GTCLDV(ALAM)
          IF(MNAER.GE.6)RETURN
      ENDIF

!     IF NOVAM, ....
      IF(LNOVAM)THEN

!         DETERMINE BRACKETING WAVELENGTHS:
          IF(ALAM.LT.WLNOV(1))THEN
              DO IAER=8,(8-1)+(NNOV-1)
                  JAER=IAER-7
                  EXTV(IAER)=EXTNOV(JAER,1)
                  ABSV(IAER)=ABSNOV(JAER,1)
                  ASYV(1,IAER)=ASMNOV(JAER,1)
                  ISTRM1=1
                  DO ISTR=2,NSTR
                      ASYV(ISTR,IAER)=ASYV(ISTRM1,IAER)*ASYV(1,IAER)
                      ISTRM1=ISTR
                  ENDDO
              ENDDO
          ELSEIF(ALAM.GT.WLNOV(NWLNOV))THEN
              DO IAER=8,(8-1)+(NNOV-1)
                  JAER=IAER-7
                  EXTV(IAER)=EXTNOV(JAER,NWLNOV)
                  ABSV(IAER)=ABSNOV(JAER,NWLNOV)
                  ASYV(1,IAER)=ASMNOV(JAER,NWLNOV)
                  ISTRM1=1
                  DO ISTR=2,NSTR
                      ASYV(ISTR,IAER)=ASYV(ISTRM1,IAER)*ASYV(1,IAER)
                      ISTRM1=ISTR
                  ENDDO
              ENDDO
          ELSE
              IWAVM1=1
              DO IWAV=2,NWLNOV
                  IF(ALAM.LE.WLNOV(IWAV))GOTO 20
                  IWAVM1=IWAV
              ENDDO
   20         CONTINUE
              FACTOR=(WLNOV(IWAV)-ALAM)/(WLNOV(IWAV)-WLNOV(IWAVM1))
              DO IAER=8,(8-1)+(NNOV-1)
                  JAER=IAER-7
                  EXTV(IAER)=EXTNOV(JAER,IWAV)                          &
     &              +FACTOR*(EXTNOV(JAER,IWAVM1)-EXTNOV(JAER,IWAV))
                  ABSV(IAER)=ABSNOV(JAER,IWAV)                          &
     &              +FACTOR*(ABSNOV(JAER,IWAVM1)-ABSNOV(JAER,IWAV))
                  ASYV(1,IAER)=ASMNOV(JAER,IWAV)                        &
     &              +FACTOR*(ASMNOV(JAER,IWAVM1)-ASMNOV(JAER,IWAV))
                  ISTRM1=1
                  DO ISTR=2,NSTR
                      ASYV(ISTR,IAER)=ASYV(ISTRM1,IAER)*ASYV(1,IAER)
                      ISTRM1=ISTR
                  ENDDO
              ENDDO
          ENDIF
      ENDIF

!     NOW BACK TO THE REGULAR AEROSOLS (EXCLUDING NOVAM)

!     FIRST CHECK TO SEE IF USER-SUPPLIED SPECTRA EXISTS
!     FOR THE 4 AEROSOLS.  IF SO, TAKE CARE OF IT HERE.
      DO IAER=1,4
          IF(NARSPC(IAER).GT.0)THEN
              NUSS=NARSPC(IAER)
              VMX=10000./VARSPC(IAER,1)
              VMN=10000./VARSPC(IAER,NUSS)
              IF(V.GE.VMX)THEN
                  EXTV(IAER)=EXTC(IAER,1)
                  ABSV(IAER)=ABSC(IAER,1)
                  ASYV(1,IAER)=ASYM(IAER,1)
              ELSEIF(V.LE.VMN)THEN
                  EXTV(IAER)=EXTC(IAER,NUSS)
                  ABSV(IAER)=ABSC(IAER,NUSS)
                  ASYV(1,IAER)=ASYM(IAER,NUSS)
              ELSE

!                 DETERMINE BRACKETING WAVELENGTH
                  NWAVM1=NUSS-1
                  IWAVM1=1
                  DO IWAV=2,NWAVM1
                      IF(ALAM.LE.VARSPC(IAER,IWAV))GOTO 30
                      IWAVM1=IWAV
                  ENDDO
   30             CONTINUE
                  FACTOR=(VARSPC(IAER,IWAV)-ALAM)                       &
     &                  /(VARSPC(IAER,IWAV)-VARSPC(IAER,IWAVM1))
                  EXTV(IAER)=EXTC(IAER,IWAV)                            &
     &                      +FACTOR*(EXTC(IAER,IWAVM1)-EXTC(IAER,IWAV))
                  ABSV(IAER)=ABSC(IAER,IWAV)                            &
     &                      +FACTOR*(ABSC(IAER,IWAVM1)-ABSC(IAER,IWAV))
                  ASYV(1,IAER)=ASYM(IAER,IWAV)                          &
     &                      +FACTOR*(ASYM(IAER,IWAVM1)-ASYM(IAER,IWAV))
              ENDIF
              IF(IPH.EQ.0)ASYV(1,IAER)=G
              ISTRM1=1
              DO ISTR=2,NSTR
                  ASYV(ISTR,IAER)=ASYV(ISTRM1,IAER)*ASYV(1,IAER)
                  ISTRM1=ISTR
              ENDDO
          ENDIF
      ENDDO

!     COMPUTE MICROWAVE ATTENUATION COEFFICIENTS (THE AWCCON FACTOR IS
!     INCLUDED IN CALLS TO GAMFOG BECAUSE EXTV AND ABSV ARE MULTIPLIED
!     BY COLUMN DENSITIES WHICH HAVE BEEN DIVIDED BY THIS FACTOR).
      VMN=10000./VX0(NWAVLN)
      IF(V.LE.VMN)THEN
          DO 40 IAER=MNAER,MXAER
              IF(IAER.LE.4)THEN
                  IF(NARSPC(IAER).GT.0)GOTO 40
              ENDIF
              EXTV(IAER)=GAMFOG(IAER,V,TAER(IAER),AWCCON(IAER))
              ABSV(IAER)=EXTV(IAER)
              DO ISTR=1,NSTR
                  ASYV(ISTR,IAER)=0.
              ENDDO
   40     CONTINUE
          IF(IPH.EQ.0)GOTO 90
          RETURN
      ENDIF

!     DETERMINE BRACKETING WAVELENGTHS:
      NWAVM1=NWAVLN-1
      IWAVM1=1
      DO IWAV=2,NWAVM1
          IF(ALAM.LE.VX0(IWAV))GOTO 60
          IWAVM1=IWAV
      ENDDO

!     ALAM IS BETWEEN 289.7 AND 299.9 MICRONS.  DEFINE THE
!     299.9 MICRON SPECTRAL DATA BEFORE INTERPOLATION
      FACTOR=(VX0(NWAVLN)-ALAM)/(VX0(NWAVLN)-VX0(NWAVM1))
      DO 50 IAER=MNAER,MXAER
          IF(IAER.LE.4)THEN
              IF(NARSPC(IAER).GT.0)GOTO 50
          ENDIF
          GMFOG=GAMFOG(IAER,VMN,TAER(IAER),AWCCON(IAER))
          EXTV(IAER)=GMFOG+FACTOR*(EXTC(IAER,NWAVM1)-GMFOG)
          ABSV(IAER)=GMFOG+FACTOR*(ABSC(IAER,NWAVM1)-GMFOG)
          ASYV(1,IAER)=FACTOR*ASYM(IAER,NWAVM1)
          ISTRM1=1
          DO ISTR=2,NSTR
              ASYV(ISTR,IAER)=ASYV(ISTRM1,IAER)*ASYV(1,IAER)
              ISTRM1=ISTR
          ENDDO
   50 CONTINUE
      IF(IPH.EQ.0)GOTO 90
      RETURN

!     COMPUTE INFRARED SPECTRAL DATA
   60 CONTINUE
      FACTOR=(VX0(IWAV)-ALAM)/(VX0(IWAV)-VX0(IWAVM1))
      DO 80 IAER=MNAER,MXAER
          IF(IAER.LE.4)THEN
              IF(NARSPC(IAER).GT.0)GOTO 80
          ENDIF
          EXTV(IAER)=EXTC(IAER,IWAV)                                    &
     &        +FACTOR*(EXTC(IAER,IWAVM1)-EXTC(IAER,IWAV))
          ASYV(1,IAER)=ASYM(IAER,IWAV)                                  &
     &        +FACTOR*(ASYM(IAER,IWAVM1)-ASYM(IAER,IWAV))
          RATLAM=VX0(IWAVM1)/ALAM
          IF(IAER.LE.2)THEN
             IF(LSSALB.LT.0)THEN

!               no user-defined aerosol single scattering albedos.
                ABSV(IAER) = ABSC(IAER,IWAV)                            &
     &               +FACTOR*(ABSC(IAER,IWAVM1)-ABSC(IAER,IWAV))
                IF(LSSALB.LT.-1)THEN
                   ABSV(IAER)=ACOALB*ABSV(IAER)
                   IF(ABSV(IAER).GT.EXTV(IAER))ABSV(IAER)=EXTV(IAER)
                ENDIF
             ELSE
                IF(ALAM.LT.AWAVLN(0))THEN

!                  wavelength < single scattering albedo spectral grid.
                   ABSV(IAER)=ABSC(IAER,IWAV)+                          &
     &                  FACTOR*(ABSC(IAER,IWAVM1)-ABSC(IAER,IWAV))
                   ABSV(IAER)=                                          &
     &                  EXTV(IAER)-ASCLMN(IAER)*(EXTV(IAER)-ABSV(IAER))
                   IF(ABSV(IAER).LT.0.)ABSV(IAER)=0.
                ELSEIF(ALAM.GT.AWAVLN(LSSALB))THEN

!                  wavelength > single scattering albedo spectral grid.
                   ABSV(IAER)=ABSC(IAER,IWAV)+                          &
     &                  FACTOR*(ABSC(IAER,IWAVM1)-ABSC(IAER,IWAV))
                   ABSV(IAER)=                                          &
     &                  EXTV(IAER)-ASCLMX(IAER)*(EXTV(IAER)-ABSV(IAER))
                   IF(ABSV(IAER).LT.0.)ABSV(IAER)=0.
                ELSE
                   IF(IAER.EQ.1)THEN
                      DO ISSALB=1,LSSALB
                         IF(ALAM.LE.AWAVLN(ISSALB))GOTO 70
                      ENDDO
   70                 CONTINUE
                      COALB=1-ASSALB(ISSALB)                            &
     &                  *(AWAVLN(ISSALB)/ALAM)**AASTRM(ISSALB)
                   ENDIF
                   ABSV(IAER)=COALB*EXTV(IAER)
                ENDIF
                ASYV(1,IAER) = ASYM(IAER,IWAV)                          &
     &               +FACTOR*(ASYM(IAER,IWAVM1)-ASYM(IAER,IWAV))
                IF(LSSALB.NE.-1 .AND. ASTROM(IAER,IWAV).GT.0. .AND.     &
     &            EXTC(IAER,IWAV).GT.0.)THEN
                   COALB=ABSV(IAER)/EXTV(IAER)
                   EXTV(IAER)=ASTMO+(EXTC(IAER,IWAVM1)-ASTMO)           &
     &                  *RATLAM**ASTROM(IAER,IWAVM1)
                   ABSV(IAER)=COALB*EXTV(IAER)
                ENDIF
             ENDIF
          ELSE
             ABSV(IAER)=ABSC(IAER,IWAV)                                 &
     &            +FACTOR*(ABSC(IAER,IWAVM1)-ABSC(IAER,IWAV))
          ENDIF
          IF(LSSALB.NE.-1 .AND. ASTROM(IAER,IWAV).GT.0.)THEN
             IF(IAER.LE.2 .AND. EXTC(IAER,IWAV).GT.0.)THEN
                COALB=ABSV(IAER)/EXTV(IAER)
                EXTV(IAER)=ASTMO+(EXTC(IAER,IWAVM1)-ASTMO)              &
     &               *RATLAM**ASTROM(IAER,IWAVM1)
                ABSV(IAER)=COALB*EXTV(IAER)
             ELSEIF(IAER.GT.2)THEN
                COALB=ABSV(IAER)/EXTV(IAER)
                EXTV(IAER)=EXTC(IAER,IWAVM1)*RATLAM**ASTROM(IAER,IWAVM1)
                ABSV(IAER)=COALB*EXTV(IAER)
             ENDIF
          ENDIF
          ISTRM1=1
          DO ISTR=2,NSTR
              ASYV(ISTR,IAER)=ASYV(ISTRM1,IAER)*ASYV(1,IAER)
              ISTRM1=ISTR
          ENDDO
   80 CONTINUE
      IF(IPH.NE.0)RETURN

!     OVERWRITE AEROSOL ASYMMETRY FACTORS WITH USER-DEFINED VALUE
   90 CONTINUE
      DO IAER=MNAER,5
          ASYV(1,IAER)=G
          ISTRM1=1
          DO ISTR=2,NSTR
              ASYV(ISTR,IAER)=ASYV(ISTRM1,IAER)*G
              ISTRM1=ISTR
          ENDDO
      ENDDO
      IF(LNOVAM)THEN
          DO IAER=8,(8-1)+(NNOV-1)
              ASYV(1,IAER)=G
              ISTRM1=1
              DO ISTR=2,NSTR
                  ASYV(ISTR,IAER)=ASYV(ISTRM1,IAER)*G
                  ISTRM1=ISTR
              ENDDO
          ENDDO
      ENDIF
      RETURN
      END
