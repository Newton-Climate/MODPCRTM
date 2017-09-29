      LOGICAL FUNCTION LRDSAP(LEVEL)

!     LRDSAP READS IN AND CHECKS DATA FROM THE SpecAerProf.dat FILE.
!     RETURNS TRUE IF LEVEL CONTAINS SPECTRAL AEROSOL DATA.
      IMPLICIT NONE

!     PARAMETERS:
!       TWOPI    TWO TIMES PI.
!       FOURPI   FOUR TIMES PI.
!       ALTTOL   ATMOSPHERIC LEVEL TOLERANCE (m).
      DOUBLE PRECISION TWOPI,FOURPI,ALTTOL
!     PARAMETER(TWOPI=6.283185307179586D0,FOURPI=2*TWOPI,ALTTOL=.001D0)
      PARAMETER(TWOPI=6.283185307179586D0,FOURPI=2*TWOPI,ALTTOL=.005D0)
      INCLUDE 'PARAMS.h'

!     INPUT ARGUMENTS:
!       LEVEL    ATMOSPHERIC LEVEL INDEX.
      INTEGER LEVEL

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /SAP/
!       NWVSAP   NUMBER OF AEROSOL SPECTRAL GRID POINTS.
!       NLGSAP   HIGHEST AEROSOL PHASE FUNCTION LEGENDRE MOMENT.
!       NANSAP   NUMBER OF AEROSOL PHASE FUNCTION ANGULAR POINTS.
!       LEVSAP   NUMBER OF AEROSOL ATMOSPHERIC LEVELS.
!       ANGSAP   AEROSOL PHASE FUNCTION ANGULAR GRID [DEGREES].
!       COSSAP   COSINE OF THE AEROSOL PHASE FUNCTION ANGLES.
!       WAVSAP   AEROSOL SPECTRAL GRID [MICRONS].
!       LEGSAP   AEROSOL LEGENDRE COEFFICIENTS.
!       PFSAP    AEROSOL PHASE FUNCTION [SR-1].
!       LOSSAP   LOGICAL FLAG, TRUE FOR LINE-OF-SIGHT PATH.
      INTEGER NWVSAP,NLGSAP,NANSAP,LEVSAP
      DOUBLE PRECISION ANGSAP,COSSAP,PFSAP
      REAL WAVSAP,LEGSAP
      LOGICAL LOSSAP
      COMMON/SAP/NWVSAP,NLGSAP,NANSAP,LEVSAP,ANGSAP(MANSAP),            &
     &  COSSAP(MANSAP),WAVSAP(MWVSAP),LEGSAP(MLGSAP,MWVSAP,LAYDIM),     &
     &  PFSAP(MANSAP,MWVSAP,LAYDIM),LOSSAP

!     /SAPDAT/
!       BSAPX    AEROSOL EXTINCTION COEFFICIENTS [KM-1].
!       BSAPA    AEROSOL ABSORPTION COEFFICIENTS [KM-1].
      REAL BSAPX,BSAPA
      COMMON/SAPDAT/BSAPX(MWVSAP,LAYDIM),BSAPA(MWVSAP,LAYDIM)

!     /MPROF/
!       ZM       PROFILE LEVEL ALTITUDES [KM].
!       PM       PROFILE LEVEL PRESSURES [MBAR].
!       TM       PROFILE LEVEL TEMPERATURES [K].
!       RFNDX    PROFILE LEVEL REFRACTIVITIES.
!       LRHSET   FLAG, .TRUE. IF RELATIVE HUMIDITY IS NOT TO BE SCALED.
      DOUBLE PRECISION ZM
      REAL PM,TM,RFNDX
      LOGICAL LRHSET
      COMMON/MPROF/ZM(LAYDIM),PM(LAYDIM),TM(LAYDIM),                    &
     &  RFNDX(LAYDIM),LRHSET(LAYDIM)

!     LOCAL VARIABLES:
!       IOS      RESULT OF IOSTAT CHECK.
!       IWVSAP   SPECTRAL GRID INDEX FOR SPECTRAL AEROSOL PROFILES.
!       IWAVM1   IWVSAP MINUS ONE.
!       ILGSAP   LEGENDRE MOMENT INDEX FOR SPECTRAL AEROSOL PROFILES.
!       IANSAP   AEROSOL PHASE FUNCTION ANGULAR GRID INDEX.
!       IANGM1   IANSAP MINUS ONE.
!       LEVM1    PREVIOUS (LOWER) ATMOSPHERIC LEVEL INDEX.
!       ICOEF    COEFFICIENT USED IN COMPUTING FORWARD SCATTER PEAK.
!       ZLEVEL   ATMOSPHERIC LEVEL [M].
!       PFNORM   PHASE FUNCTION NORMALIZATION [SR-1].
!       PFRAT    RATIO OF PHASE FUNCTION VALUES AN NEIGHBORING ANGLES.
!       FORWAR   FORWARD SCATTER PHASE FUNCTION VALUE.
!       BSAPS    AEROSOL SCATTERING COEFFICIENT [KM-1].
!       LG0SAP   FIRST LEGENDRE MOMENT (NORMALIZATION).
!       WAVLEN   AEROSOL SPECTRAL GRID [MICRONS].
      INTEGER IOS,IWVSAP,IWAVM1,ILGSAP,IANSAP,IANGM1,LEVM1,ICOEF
      DOUBLE PRECISION ZLEVEL,PFNORM,PFRAT
      REAL FORWAR,BSAPS(MWVSAP),LG0SAP(MWVSAP),WAVLEN(MWVSAP)

!     SAVED VARIABLES:
!       PFWARN   WARNS IF PHASE FUNCTION IS NOT PROPERLY NORMALIZED.
!       LGWARN   WARNS IF LEGENDRE SUM DOES NOT AGREE WITH FORWARD PEAK.
      LOGICAL PFWARN,LGWARN
      SAVE PFWARN,LGWARN

!     FIRST READ OF SPECTRAL DATA?
      IF(LEVEL.EQ.1)THEN

!         INITIALIZE PFWARN
          PFWARN=.TRUE.
          LGWARN=.TRUE.

!         READ LEVEL DATA:
          READ(ISAP,*,IOSTAT=IOS)(ZLEVEL,WAVSAP(IWVSAP),                &
     &      BSAPX(IWVSAP,LEVEL),BSAPS(IWVSAP),LG0SAP(IWVSAP),           &
     &      (LEGSAP(ILGSAP,IWVSAP,LEVEL),ILGSAP=1,NLGSAP),              &
     &      (PFSAP(IANSAP,IWVSAP,LEVEL),IANSAP=1,NANSAP),               &
     &      IWVSAP=1,NWVSAP)
          IF(IOS.NE.0)THEN
              WRITE(IPR,'(/A)')' Error in LRDSAP:  Unable to read'//    &
     &          ' atmosphere level   1 data from SpecAerProf.dat file.'
              STOP 'Unable to read data from SpecAerProf.dat file.'
          ENDIF

!         CHECK ALTITUDE:
          IF(ABS(1000*ZM(LEVEL)-ZLEVEL).GT.ALTTOL)THEN
              WRITE(IPR,'(/A,/(18X,A,F12.6))')                          &
     &          ' Error in LRDSAP:  Altitude mismatch.',                &
     &          '  profile altitude    (km):',ZM(LEVEL),                &
     &          'SpecAerProf.dat altitude (km):',ZLEVEL/1000
              STOP  'Error in LRDSAP:  Altitude mismatch.'
          ENDIF

!         CHECK SPECTRAL GRID:
          IF(WAVSAP(1).LT.0.)THEN
              WRITE(IPR,'(/2A,F12.6,A)')' Error in LRDSAP: ',           &
     &          ' First spectral wavelength in SpecAerProf.dat'//       &
     &          ' file is negative:',WAVSAP(1),' microns.'
              STOP 'Error in LRDSAP:  First wavelength is negative.'
          ENDIF
          IWAVM1=1
          DO IWVSAP=2,NWVSAP
              IF(WAVSAP(IWVSAP).LE.WAVSAP(IWAVM1))THEN
                  WRITE(IPR,'(/2A,/(17X,A,I4,A,F12.6))')' Error',       &
     &              ' in LRDSAP:  Spectral grid in SpecAerProf.dat'     &
     &              //' file is not monotonically increasing',          &
     &              ' Wavelength',IWAVM1,' (microns):',WAVSAP(IWAVM1),  &
     &              ' Wavelength',IWVSAP,' (microns):',WAVSAP(IWVSAP)
                  STOP 'Error in LRDSAP: Wavelength grid not monotonic.'
              ENDIF
              IWAVM1=IWVSAP
          ENDDO

!         SET LEVSAP TO ZERO TO FLAG ADDITIONAL DATA:
          LEVSAP=0
      ELSEIF(LEVSAP.GT.0)THEN

!         NO MORE DATA TO READ SO JUST RETURN:
          LRDSAP=.FALSE.
          RETURN
      ELSE

!         READ LEVEL DATA:
          READ(ISAP,*,IOSTAT=IOS)(ZLEVEL,WAVLEN(IWVSAP),                &
     &      BSAPX(IWVSAP,LEVEL),BSAPS(IWVSAP),LG0SAP(IWVSAP),           &
     &      (LEGSAP(ILGSAP,IWVSAP,LEVEL),ILGSAP=1,NLGSAP),              &
     &      (PFSAP(IANSAP,IWVSAP,LEVEL),IANSAP=1,NANSAP),               &
     &      IWVSAP=1,NWVSAP)
          IF(IOS.NE.0)THEN
              WRITE(IPR,'(/A,I4,A)')                                    &
     &          ' Error in LRDSAP:  Unable to read atmosphere level',   &
     &          LEVEL,' data from SpecAerProf.dat file'
              STOP 'Unable to read data from SpecAerProf.dat file.'
          ENDIF

!         CHECK ALTITUDE:
          IF(ABS(1000*ZM(LEVEL)-ZLEVEL).GT.ALTTOL)THEN
              WRITE(IPR,'(/A,/(18X,A,F12.6))')                          &
     &          ' Error in LRDSAP:  Altitude mismatch.',                &
     &          '  profile altitude    (km):',ZM(LEVEL),                &
     &          'SpecAerProf.dat altitude (km):',ZLEVEL/1000
              STOP  'Error in LRDSAP:  Altitude mismatch.'
          ENDIF

!         CHECK SPECTRAL GRID:
          DO IWVSAP=1,NWVSAP
              IF(WAVSAP(IWVSAP).NE.WAVLEN(IWVSAP))THEN
                  WRITE(IPR,'(/A,/(17X,A,I4,A,2F12.6))')                &
     &              ' Error in LRDSAP:  Spectral grid mismatch.'        &
     &              //' Wavelength',IWVSAP,' (microns):',               &
     &              WAVSAP(IWVSAP),WAVLEN(IWVSAP)
                  STOP 'Error in LRDSAP:  Spectral grid mismatch.'
              ENDIF
          ENDDO

!         CHECK FOR FINAL LEVEL:
          IF(LG0SAP(1).EQ.0.)THEN

!             LG0SAP(1) EQUALS ZERO FLAGS END OF NON-ZERO DATA:
              LEVM1=LEVEL-1
              WRITE(IPR,'(/I4,A)')LEVM1,' levels of non-zero data'      &
     &          //' was read from SpecAerProf.dat file.'
              LEVSAP=LEVEL
              DO IWVSAP=1,NWVSAP

!                 SET COEFFICIENTS TO ZERO:
                  BSAPX(IWVSAP,LEVEL)=0.
                  BSAPA(IWVSAP,LEVEL)=0.

!                 SET PHASE FUNCTIONS TO LOWER LEVEL VALUES:
                  DO ILGSAP=1,NLGSAP
                      LEGSAP(ILGSAP,IWVSAP,LEVEL)                       &
     &                  =LEGSAP(ILGSAP,IWVSAP,LEVM1)
                  ENDDO
                  DO IANSAP=1,NANSAP
                      PFSAP(IANSAP,IWVSAP,LEVEL)                        &
     &                  =PFSAP(IANSAP,IWVSAP,LEVM1)
                  ENDDO
              ENDDO
              LRDSAP=.FALSE.
              RETURN
          ENDIF
      ENDIF

!     CHECK SPECTRAL DATA:
      DO IWVSAP=1,NWVSAP

!         CHECK EXTINCTION COEFFICIENTS:
          IF(BSAPX(IWVSAP,LEVEL).LT.0.)THEN
              WRITE(IPR,'(/A,2(F12.6,A),/21X,F12.6,A)')                 &
     &          ' Error in LRDSAP:  Negative extinction coefficient at',&
     &          ZLEVEL,'m altitude and',WAVSAP(IWVSAP),' microns:',     &
     &          BSAPX(IWVSAP,LEVEL),' km-1.'
              STOP 'Error in LRDSAP:  Negative extinction coefficient'
          ENDIF

!         CHECK SCATTERING COEFFICIENTS:
          IF(BSAPS(IWVSAP).LT.0. .OR.                                   &
     &       BSAPS(IWVSAP).GT.BSAPX(IWVSAP,LEVEL))THEN
              WRITE(IPR,'(/A,2(F12.6,A),/(21X,A,F12.6))')' Error in'//  &
     &          ' LRDSAP:  Scattering albedo not between 0 and 1 at',   &
     &          ZLEVEL,'m altitude and',WAVSAP(IWVSAP),' microns:',     &
     &          'Scattering coefficient (km-1):',BSAPS(IWVSAP),         &
     &          'Extinction coefficient (km-1):',BSAPX(IWVSAP,LEVEL)
              STOP                                                      &
     &          'Error in LRDSAP: Scattering albedo not between 0 and 1'
          ENDIF
          BSAPA(IWVSAP,LEVEL)=BSAPX(IWVSAP,LEVEL)-BSAPS(IWVSAP)

!         NORMALIZE LEGENDRE EXPANSION AND COMPUTE FORWARD SCATTER:
          ICOEF=1
          FORWAR=1.
          DO ILGSAP=1,NLGSAP
              LEGSAP(ILGSAP,IWVSAP,LEVEL)                               &
     &          =LEGSAP(ILGSAP,IWVSAP,LEVEL)/LG0SAP(IWVSAP)
              ICOEF=ICOEF+2
              FORWAR=FORWAR+ICOEF*LEGSAP(ILGSAP,IWVSAP,LEVEL)
          ENDDO

!         CHECK PHASE FUNCTION VALUES AND NORMALIZATION (ASSUME
!         EXPONENTIAL VARIATION BETWEEN COSINE ANGLE VALUES):
          PFNORM=0.D0
          IF(PFSAP(1,IWVSAP,LEVEL).LE.0.D0)THEN
              WRITE(IPR,'(/A,3(F12.6,A),/21X,1P,E15.6,A))')' Error in'  &
     &          //' LRDSAP:  Non-positive phase function value at',     &
     &          ZLEVEL,'m altitude,',WAVSAP(IWVSAP),' microns and',     &
     &          ANGSAP(1),' deg.',PFSAP(1,IWVSAP,LEVEL),' sr-1'
              STOP 'Error in LRDSAP:  Negative phase function value.'
          ENDIF
          IANGM1=1
          DO IANSAP=2,NANSAP
              IF(PFSAP(IANSAP,IWVSAP,LEVEL).LE.0.D0)THEN
                  WRITE(IPR,'(/A,3(F12.6,A),/21X,1P,E15.6,A))')         &
     &              ' Error in LRDSAP:  Non-positive phase function'    &
     &              //' value at',ZLEVEL,'m altitude,',WAVSAP(IWVSAP),  &
     &              ' microns and',ANGSAP(IANSAP),' deg.',              &
     &              PFSAP(IANSAP,IWVSAP,LEVEL),' sr-1'
                  STOP 'Error in LRDSAP: Negative phase function value.'
              ENDIF
              PFRAT                                                     &
     &          =PFSAP(IANGM1,IWVSAP,LEVEL)/PFSAP(IANSAP,IWVSAP,LEVEL)
              IF(ABS(PFRAT-1).GT..00001D0)THEN
                  PFNORM=PFNORM+(COSSAP(IANGM1)-COSSAP(IANSAP))*        &
     &              (PFSAP(IANGM1,IWVSAP,LEVEL)                         &
     &              -PFSAP(IANSAP,IWVSAP,LEVEL))/LOG(PFRAT)
              ELSE
                  PFNORM=PFNORM+2*(COSSAP(IANGM1)-COSSAP(IANSAP))*      &
     &              PFSAP(IANGM1,IWVSAP,LEVEL)/(3-PFRAT)
              ENDIF
              IANGM1=IANSAP
          ENDDO
          PFNORM=TWOPI*PFNORM
          IF(ABS(PFNORM-1).GT..00001D0)THEN
              IF(PFWARN .AND. ABS(PFNORM-1).GT..01D0)THEN
                  WRITE(IPR,'(/A,2(F12.6,A),/21X,F12.6,A)')             &
     &              ' Warning from LRDSAP:  At',ZLEVEL,'m altitude and',&
     &              WAVSAP(IWVSAP),' microns, the aerosol phase'        &
     &              //' function was normalized to',PFNORM,             &
     &              ' instead of 1; it is being renormalized.'          &
     &              //'  This warning will not be repeated.'
                  PFWARN=.FALSE.
              ENDIF
              DO IANSAP=1,NANSAP
                  PFSAP(IANSAP,IWVSAP,LEVEL)                            &
     &              =PFSAP(IANSAP,IWVSAP,LEVEL)/PFNORM
              ENDDO
          ENDIF

!         CHECK THE FORWARD PEAK:
          IF(ABS(DBLE(FORWAR)/PFSAP(1,IWVSAP,LEVEL)-FOURPI).GT..001D0   &
     &      .AND. LGWARN)THEN
              WRITE(IPR,'(/A,F12.6,A,/(21X,2(F12.6,A)))')               &
     &          ' Warning from LRDSAP:  The aerosol phase function'     &
     &          //' forward peak at',ZLEVEL,'m altitude and',           &
     &          WAVSAP(IWVSAP),' microns computed from the'             &
     &          //' Legendre expansion is',DBLE(FORWAR)/FOURPI,         &
     &          ' sr-1; however, the tabulated value is',               &
     &          PFSAP(1,IWVSAP,LEVEL),' sr-1.  Before normalization,'// &
     &          ' the tabulated value was',PFNORM*PFSAP(1,IWVSAP,LEVEL)
              LGWARN=.FALSE.
          ENDIF
      ENDDO

!     RETURN TO AERNSM:
      LRDSAP=.TRUE.
      RETURN
      END
