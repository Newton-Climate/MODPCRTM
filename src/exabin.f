      SUBROUTINE EXABIN(ICH,ASTMC,ASTMX,LDASTM,RHASYM,AERRH,LASTM)

!     LOADS EXTINCTION, ABSORPTION AND ASYMMETRY COEFFICIENTS
!     FOR THE FOUR AEROSOL ALTITUDE REGIONS
!     MODIFIED FOR ASYMMETRY - JAN 1986 (A.E.R. INC.)
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     ARGUMENTS:
      INTEGER ICH(4)
      REAL ASTMC,ASTMX,RHASYM,AERRH
      LOGICAL LASTM,LDASTM(2)

!     COMMONS:
      INCLUDE 'BASE.h'
      INCLUDE 'SEGDAT.h'

!     /CARD2/
!       IHAZE    BOUNDARY LAYER AEROSOL MODEL NUMBER.
!       ISEASN   SEASON NUMBER (1=SPRING-SUMMER, 2=FALL-WINTER).
!       IVULCN   VOLCANIC AEROSOL MODEL NUMBER.
!       ICSTL    COASTAL AIRMASS MODEL NUMBER.
!       ICLD     CLOUD MODEL NUMBER.
!       IVSA     VERTICAL STRUCTURE ALGORITHM (0=OFF, 1=ON).
!       VIS      SURFACE VISIBILITY (GROUND METEOROLOGICAL RANGE) [KM].
!       WSS      CURRENT WIND SPEED (M/S).
!       WHH      24-HOUR WIND SPEED (M/S).
!       RAINRT   RAIN RATE (MM/HR)
!       LSAP     LOGICAL FLAG FOR SPECTRAL AEROSOL PROFILES INPUT.
      INTEGER IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA
      REAL VIS,WSS,WHH,RAINRT
      LOGICAL LSAP
      COMMON/CARD2/IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,                 &
     &  VIS,WSS,WHH,RAINRT,LSAP

!     /CARD2D/
      INTEGER IREG,IREGC
      DOUBLE PRECISION ALTB
      COMMON/CARD2D/IREG(4),ALTB(4),IREGC(4)

!           0-2KM
!             RUREXT=RURAL EXTINCTION   RURABS=RURAL ABSORPTION
!             RURSYM=RURAL ASYMMETRY FACTORS
!             URBEXT=URBAN EXTINCTION   URBABS=URBAN ABSORPTION
!             URBSYM=URBAN ASYMMETRY FACTORS
!             OCNEXT=MARITIME EXTINCTION  OCNABS=MARITIME ABSORPTION
!             OCNSYM=MARITIME ASYMMETRY FACTORS
!             TROEXT=TROPSPHER EXTINCTION  TROABS=TROPOSPHER ABSORPTION
!             TROSYM=TROPSPHERIC ASYMMETRY FACTORS
!             FG1EXT=FOG1 .2KM VIS EXTINCTION  FG1ABS=FOG1 ABSORPTION
!             FG1SYM=FOG1 ASYMMETRY FACTORS
!             FG2EXT=FOG2 .5KM VIS EXTINCTION  FG2ABS=FOG2 ABSORPTION
!             FG2SYM=FOG2 ASYMMETRY FACTORS
!           >2-10KM
!             TROEXT=TROPOSPHERIC EXTINCTION
!             TROABS=TROPOSPHERIC ABSORPTION
!             TROSYM=TROPOSPHERIC ASYMMETRY FACTORS
!           >10-30KM
!             BSTEXT=BACKGROUND STRATOSPHERIC EXTINCTION
!             BSTABS=BACKGROUND STRATOSPHERIC ABSORPTION
!             BSTSYM=BACKGROUND STRATOSPHERIC ASYMMETRY FACTORS
!             AVOEXT=AGED VOLCANIC EXTINCTION
!             AVOABS=AGED VOLCANIC ABSORPTION
!             AVOSYM=AGED VOLCANIC ASYMMETRY FACTORS
!             FVOEXT=FRESH VOLCANIC EXTINCTION
!             FVOABS=FRESH VOLCANIC ABSORPTION
!             FVOSYM=FRESH VOLCANIC ASYMMETRY FACTORS
!           >30-100KM
!             DMEEXT=METEORIC DUST EXTINCTION
!             DMEABS=METEORIC DUST ABSORPTION
!             DMESYM=METEORIC DUST ASYMMETRY FACTORS

!     AEROSOL EXTINCTION AND ABSORPTION DATA:

!     MODIFIED TO INCLUDE ASYMMETRY DATA - JAN 1986 (A.E.R. INC.)
      REAL RUREXT,RURABS,RURSYM,URBEXT,URBABS,URBSYM,                   &
     &  OCNEXT,OCNABS,OCNSYM
      COMMON/EXTD1/RUREXT(NWAVLN,4),RURABS(NWAVLN,4),RURSYM(NWAVLN,4),  &
     &             URBEXT(NWAVLN,4),URBABS(NWAVLN,4),URBSYM(NWAVLN,4),  &
     &             OCNEXT(NWAVLN,4),OCNABS(NWAVLN,4),OCNSYM(NWAVLN,4)
      REAL TROEXT,TROABS,TROSYM,                                        &
     &  FG1EXT,FG1ABS,FG1SYM,FG2EXT,FG2ABS,FG2SYM,BSTEXT,BSTABS,BSTSYM, &
     &  AVOEXT,AVOABS,AVOSYM,FVOEXT,FVOABS,FVOSYM,DMEEXT,DMEABS,DMESYM
      COMMON/EXTD2/TROEXT(NWAVLN,4),TROABS(NWAVLN,4),TROSYM(NWAVLN,4),  &
     &             FG1EXT(NWAVLN),FG1ABS(NWAVLN),FG1SYM(NWAVLN),        &
     &             FG2EXT(NWAVLN),FG2ABS(NWAVLN),FG2SYM(NWAVLN),        &
     &             BSTEXT(NWAVLN),BSTABS(NWAVLN),BSTSYM(NWAVLN),        &
     &             AVOEXT(NWAVLN),AVOABS(NWAVLN),AVOSYM(NWAVLN),        &
     &             FVOEXT(NWAVLN),FVOABS(NWAVLN),FVOSYM(NWAVLN),        &
     &             DMEEXT(NWAVLN),DMEABS(NWAVLN),DMESYM(NWAVLN)
      REAL CLDSPC
      COMMON/CLDDAT/CLDSPC(NWAVLN,3,NCLDS)
      REAL CIRSPC
      COMMON/CIRR/CIRSPC(NWAVLN,3,NCIRS)

      REAL VX0
      COMMON/EXTWAV/VX0(NWAVLN)

      REAL COALB,WAVCOR

!     LOADS EXTINCTION, ABSORPTION AND ASYMMETRY COEFFICIENTS
!     FOR THE FOUR AEROSOL ALTITUDE REGIONS

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      SAVE /CLDDAT/,/CIRR/,/EXTD1/,/EXTD2/,/EXTWAV/
      EXTERNAL CLDDTA,EXTDT1,EXTDT2

!     /AERALB/
!       LSSALB   INDEX OF FINAL AEROSOL ALBEDO SPECTRAL GRID POINT.
!       ACOALB   AEROSOL SINGLE SCATTERING CO-ALBEDO SCALING FACTOR.
!       AWAVLN   AEROSOL ALBEDO SPECTRAL GRID POINTS [MICRONS].
!       ASSALB   BOUNDARY LAYER AND TROPOSPHERE AEROSOL SINGLE
!                SCATTERING ALBEDOS.
!       AASTRM   WAVELENGTH INTERPOLATION EXPONENTS (EQUIVALENT TO AN
!                ANGSTROM COEFFICIENT FOR SINGLE SCATTERING ALBEDO)
!       ASCLMN   AEROSOL (1=BOUNDARY LAYER, 2=TROPOSPHERIC)
!                SINGLE SCATTERING ALBEDO SCALING FACTORS
!                FOR WAVELENGTHS BELOW AWAVLN(0).
!       ASCLMX   AEROSOL (1=BOUNDARY LAYER, 2=TROPOSPHERIC)
!                SINGLE SCATTERING ALBEDO SCALING FACTORS
!                FOR WAVELENGTHS ABOVE AWAVLN(LSSALB).
      INTEGER LSSALB
      REAL ACOALB,AWAVLN,ASSALB,AASTRM,ASCLMN,ASCLMX
      COMMON/AERALB/LSSALB,ACOALB,AWAVLN(0:MSSALB),ASSALB(0:MSSALB),    &
     &     ASCLMN(2),AASTRM(MSSALB),ASCLMX(2)

!     /AERC/
!       WAVLEN   CLOUD/RAIN DATA WAVELENGTH GRID [MICRONS].
!       AWCCON   WATER COLUMN CONTENT.
!       EXTC     AEROSOL/CLOUD EXTINCTION COEFFICIENTS.
!       ABSC     AEROSOL/CLOUD ABSORPTION COEFFICIENTS.
!       ASYM     AEROSOL/CLOUD ASYMMETRY FACTORS.
!       ASTROM   AEROSOL ANGSTROM LAW EXPONENTS.
!       ASTMO    AEROSOL ANGSTROM LAW OFFSET.
      REAL ASTROM,ASTMO
      COMMON/AERC/ASTROM(MAER,MXWVLN),ASTMO

      REAL FACTOR,EXTV,SSALB,ABSV,DENOM
      INTEGER IWAVP1
      REAL A,A1,A2,E1,E2,EC,WRH,X,X1,X2,Y,Y1,Y2,Z1,Z2,ZK
      INTEGER IRH,IM1RH,ITA,ITAS,ITC,IAER,IWAV

!     DATA:
      REAL AFLWC,ASLWC,AVLWC,BSLWC,CULWC,FVLWC,RFLWC,SCLWC,SNLWC,       &
     &  STLWC,MDLWC,RHZONE(4),ELWCR(4),ELWCU(4),ELWCM(4),ELWCT(4)
      DATA RHZONE/0., 70., 80., 99./
      DATA ELWCR/3.517E-04, 3.740E-04, 4.439E-04, 9.529E-04/
      DATA ELWCM/4.675E-04, 6.543E-04, 1.166E-03, 3.154E-03/
      DATA ELWCU/3.102E-04, 3.802E-04, 4.463E-04, 9.745E-04/
      DATA ELWCT/1.735E-04, 1.820E-04, 2.020E-04, 2.408E-04/
      DATA AFLWC/1.295E-02/, RFLWC/1.804E-03/, CULWC/7.683E-03/
      DATA ASLWC/4.509E-03/, STLWC/5.272E-03/, SCLWC/4.177E-03/
      DATA SNLWC/7.518E-03/, BSLWC/1.567E-04/, FVLWC/5.922E-04/
      DATA AVLWC/1.675E-04/, MDLWC/4.775E-04/

!     LOOP OVER AEROSOL REGIONS:
      DO IAER=1,4
          AWCCON(IAER)=0.
          IF(IREG(IAER).EQ.0 .OR. (ICH(1).EQ.3 .AND. IAER.EQ.1))THEN
              ITA=ICH(IAER)
              ITC=ICH(IAER)-7
              ITAS=ITA
              IF(IREGC(IAER).NE.0)THEN

!                 SECTION TO LOAD EXTINCTION AND ABSORPTION
!                 COEFFICIENTS FOR CLOUD AND OR RAIN MODELS
                  DO IWAV=1,NWAVLN-1
                      IF(ICLD.EQ.2)THEN
                          ABSC(IAER,IWAV)=CLDSPC(IWAV,2,2)
                          EXTC(IAER,IWAV)=CLDSPC(IWAV,1,2)
                          ASYM(IAER,IWAV)=CLDSPC(IWAV,3,2)
                          IF(IWAV.EQ.1)AWCCON(IAER)=ASLWC
                      ELSEIF(ICLD.EQ.3 .OR. ICLD.EQ.6)THEN
                          ABSC(IAER,IWAV)=CLDSPC(IWAV,2,3)
                          EXTC(IAER,IWAV)=CLDSPC(IWAV,1,3)
                          ASYM(IAER,IWAV)=CLDSPC(IWAV,3,3)
                          IF(IWAV.EQ.1)AWCCON(IAER)=STLWC
                      ELSEIF(ICLD.EQ.4)THEN
                          ABSC(IAER,IWAV)=CLDSPC(IWAV,2,4)
                          EXTC(IAER,IWAV)=CLDSPC(IWAV,1,4)
                          ASYM(IAER,IWAV)=CLDSPC(IWAV,3,4)
                          IF(IWAV.EQ.1)AWCCON(IAER)=SCLWC
                      ELSEIF(ICLD.EQ.5 .OR. ICLD.EQ.7                   &
     &                                 .OR. ICLD.EQ.8)THEN
                          ABSC(IAER,IWAV)=CLDSPC(IWAV,2,5)
                          EXTC(IAER,IWAV)=CLDSPC(IWAV,1,5)
                          ASYM(IAER,IWAV)=CLDSPC(IWAV,3,5)
                          IF(IWAV.EQ.1)AWCCON(IAER)=SNLWC
                      ELSE
                          ABSC(IAER,IWAV)=CLDSPC(IWAV,2,1)
                          EXTC(IAER,IWAV)=CLDSPC(IWAV,1,1)
                          ASYM(IAER,IWAV)=CLDSPC(IWAV,3,1)
                          IF(IWAV.EQ.1)AWCCON(IAER)=CULWC
                      ENDIF
                  ENDDO
              ELSE
                  IF(RHASYM.GT.0. .AND. RHASYM.LE.99.)THEN
                      WRH=RHASYM
                  ELSEIF(ITA.EQ.6 .AND. IAER.GT.1)THEN

!                     DOES NOT ALLOW TROP RH DEPENDENT ABOVE EH(7,I)
!                     DEFAULTS TO TROPOSPHERIC AT 70. PERCENT
                      WRH=70.
                  ELSEIF(AERRH.GT.0.)THEN
                      WRH=AERRH
                  ELSE
                      WRH=WTOTAL(15)
                  ENDIF
                  IF(WRH.LT.RHZONE(2))THEN
                      IM1RH=1
                      IRH=2
                  ELSEIF(WRH.LT.RHZONE(3))THEN
                      IM1RH=2
                      IRH=3
                  ELSE
                      IM1RH=3
                      IRH=4
                  ENDIF
                  IF(WRH.GT.0. .AND. WRH.LT.99.)X=ALOG(100-WRH)
                  X1=ALOG(100-RHZONE(IM1RH))
                  X2=ALOG(100-RHZONE(IRH))
                  IF(WRH.GE.99.)X=X2
                  IF(WRH.LE.0.)X=X1
                  DO 10 IWAV=1,NWAVLN-1
                      ITA=ITAS
                      IF(ITA.NE.3 .OR. IAER.NE.1)THEN
                          ABSC(IAER,IWAV)=0.
                          EXTC(IAER,IWAV)=0.
                          ASYM(IAER,IWAV)=0.
                          IF(ITA.GT.6)THEN
                              IF(ITA.GT.19)THEN
                                  ABSC(IAER,IWAV)=DMEABS(IWAV)
                                  EXTC(IAER,IWAV)=DMEEXT(IWAV)
                                  ASYM(IAER,IWAV)=DMESYM(IWAV)
                                  IF(IWAV.EQ.1)AWCCON(IAER)=MDLWC
                                  GOTO 10
                              ENDIF
                              IF(ITC.GE.1)THEN
                                  IF(ITC.EQ.2)THEN
                                      ABSC(IAER,IWAV)=FG2ABS(IWAV)
                                      EXTC(IAER,IWAV)=FG2EXT(IWAV)
                                      ASYM(IAER,IWAV)=FG2SYM(IWAV)
                                      IF(IWAV.EQ.1)AWCCON(IAER)=RFLWC
                                  ELSEIF(ITC.EQ.4 .OR. ITC.EQ.9         &
     &                                            .OR. ITC.EQ.10)THEN
                                      ABSC(IAER,IWAV)=BSTABS(IWAV)
                                      EXTC(IAER,IWAV)=BSTEXT(IWAV)
                                      ASYM(IAER,IWAV)=BSTSYM(IWAV)
                                      IF(IWAV.EQ.1)AWCCON(IAER)=BSLWC
                                  ELSEIF(ITC.EQ.5 .OR. ITC.EQ.7)THEN
                                      ABSC(IAER,IWAV)=AVOABS(IWAV)
                                      EXTC(IAER,IWAV)=AVOEXT(IWAV)
                                      ASYM(IAER,IWAV)=AVOSYM(IWAV)
                                      IF(IWAV.EQ.1)AWCCON(IAER)=AVLWC
                                  ELSEIF(ITC.EQ.6 .OR. ITC.EQ.8         &
     &                                           .OR. ITC.EQ.11)THEN
                                      ABSC(IAER,IWAV)=FVOABS(IWAV)
                                      EXTC(IAER,IWAV)=FVOEXT(IWAV)
                                      ASYM(IAER,IWAV)=FVOSYM(IWAV)
                                      IF(IWAV.EQ.1)AWCCON(IAER)=FVLWC
                                  ELSEIF(ITC.EQ.12)THEN
                                      ABSC(IAER,IWAV)=DMEABS(IWAV)
                                      EXTC(IAER,IWAV)=DMEEXT(IWAV)
                                      ASYM(IAER,IWAV)=DMESYM(IWAV)
                                      IF(IWAV.EQ.1)AWCCON(IAER)=MDLWC
                                      GOTO 10
                                  ELSEIF(ITC.NE.3)THEN
                                      ABSC(IAER,IWAV)=FG1ABS(IWAV)
                                      EXTC(IAER,IWAV)=FG1EXT(IWAV)
                                      ASYM(IAER,IWAV)=FG1SYM(IWAV)
                                      IF(IWAV.EQ.1)AWCCON(IAER)=AFLWC
                                  ENDIF
                              ENDIF
                              GOTO 10
                          ELSEIF(ITA.LE.0)THEN
                              GOTO 10
                          ENDIF
                      ENDIF

!                     NAVY MARITIME AEROSOL BECOMES MARINE IN MICROWAVE
!                     IWAV=743 FOR 50.00 MICRONS.
                      IF(IWAV.GE.743 .AND. ITA.EQ.3)ITA=4

!                     RH DEPENDENT AEROSOLS
                      IF(ITA.EQ.3 .AND. IAER.EQ.1)THEN
                          A2=ALOG(OCNSYM(IWAV,IRH))
                          A1=ALOG(OCNSYM(IWAV,IM1RH))
                          A=A1+(A2-A1)*(X-X1)/(X2-X1)
                          ASYM(IAER,IWAV)=EXP(A)
                          E2=ALOG(ELWCM(IRH))
                          E1=ALOG(ELWCM(IM1RH))
                          GOTO 10
                      ELSEIF(ITA.EQ.5)THEN
                          Y2=ALOG(URBEXT(IWAV,IRH))
                          Y1=ALOG(URBEXT(IWAV,IM1RH))
                          Z2=ALOG(URBABS(IWAV,IRH))
                          Z1=ALOG(URBABS(IWAV,IM1RH))
                          A2=ALOG(URBSYM(IWAV,IRH))
                          A1=ALOG(URBSYM(IWAV,IM1RH))
                          E2=ALOG(ELWCU(IRH))
                          E1=ALOG(ELWCU(IM1RH))
                      ELSEIF(ITA.EQ.6)THEN
                          Y2=ALOG(TROEXT(IWAV,IRH))
                          Y1=ALOG(TROEXT(IWAV,IM1RH))
                          Z2=ALOG(TROABS(IWAV,IRH))
                          Z1=ALOG(TROABS(IWAV,IM1RH))
                          A2=ALOG(TROSYM(IWAV,IRH))
                          A1=ALOG(TROSYM(IWAV,IM1RH))
                          E2=ALOG(ELWCT(IRH))
                          E1=ALOG(ELWCT(IM1RH))
                      ELSEIF(ITA.EQ.3 .OR. ITA.EQ.4)THEN
                          Y2=ALOG(OCNEXT(IWAV,IRH))
                          Y1=ALOG(OCNEXT(IWAV,IM1RH))
                          Z2=ALOG(OCNABS(IWAV,IRH))
                          Z1=ALOG(OCNABS(IWAV,IM1RH))
                          A2=ALOG(OCNSYM(IWAV,IRH))
                          A1=ALOG(OCNSYM(IWAV,IM1RH))
                          E2=ALOG(ELWCM(IRH))
                          E1=ALOG(ELWCM(IM1RH))
                      ELSE
                          Y2=ALOG(RUREXT(IWAV,IRH))
                          Y1=ALOG(RUREXT(IWAV,IM1RH))
                          Z2=ALOG(RURABS(IWAV,IRH))
                          Z1=ALOG(RURABS(IWAV,IM1RH))
                          A2=ALOG(RURSYM(IWAV,IRH))
                          A1=ALOG(RURSYM(IWAV,IM1RH))
                          E2=ALOG(ELWCR(IRH))
                          E1=ALOG(ELWCR(IM1RH))
                      ENDIF
                      Y=Y1+(Y2-Y1)*(X-X1)/(X2-X1)
                      ZK=Z1+(Z2-Z1)*(X-X1)/(X2-X1)
                      A=A1+(A2-A1)*(X-X1)/(X2-X1)
                      ABSC(IAER,IWAV)=EXP(ZK)
                      EXTC(IAER,IWAV)=EXP(Y)
                      ASYM(IAER,IWAV)=EXP(A)
                      IF(IWAV.EQ.1)THEN
                          EC=E1+(E2-E1)*(X-X1)/(X2-X1)
                          AWCCON(IAER)=EXP(EC)
                      ENDIF
   10             CONTINUE
              ENDIF
          ENDIF
      ENDDO
      DO IWAV=1,NWAVLN
          IF(ICLD.EQ.18)THEN
              EXTC(5,IWAV)=CIRSPC(IWAV,1,1)
              ABSC(5,IWAV)=CIRSPC(IWAV,2,1)
              ASYM(5,IWAV)=CIRSPC(IWAV,3,1)
              AWCCON(5)=5.811E-2
          ELSEIF(ICLD.EQ.19)THEN
              EXTC(5,IWAV)=CIRSPC(IWAV,1,2)
              ABSC(5,IWAV)=CIRSPC(IWAV,2,2)
              ASYM(5,IWAV)=CIRSPC(IWAV,3,2)
              AWCCON(5)=3.446E-3
          ELSE
              ABSC(5,IWAV)=0.
              EXTC(5,IWAV)=0.
              ASYM(5,IWAV)=0.
              AWCCON(5)=0.
          ENDIF
      ENDDO

!     DEFINE ANGSTROM COEFFICIENTS FOR SPECTRAL INTERPOLATION:
      IWAV=1
      DO IWAVP1=2,NWAVLN
         DENOM=LOG(VX0(IWAV)/VX0(IWAVP1))
         DO IAER=1,5
            IF(EXTC(IAER,IWAV).GT.0. .AND. EXTC(IAER,IWAVP1).GT.0.)THEN
               ASTROM(IAER,IWAV)=                                       &
     &           LOG(EXTC(IAER,IWAVP1)/EXTC(IAER,IWAV))/DENOM
            ENDIF
         ENDDO
         IWAV=IWAVP1
      ENDDO

!     INCORPORATE USER-DEFINED ANGSTROM LAW FOR THE BOUNDARY
!     LAYER AND TROPOSPHERIC AEROSOLS.
!     NOTE:  VX0(71)=0.5495 MICRON.
      IF(LASTM)THEN

          DO IAER=1,2
              IF(EXTC(IAER,71).GT.0.)THEN
                  DO IWAV=1,788
                      IF(EXTC(IAER,IWAV).GT.0.)THEN

!                         REPLACE SPECTRAL DEPENDENCE OF EXTINCTION WITH
!                         THE ANGSTROM LAW, BUT CONSERVE CO-ALBEDO.
                          COALB=ABSC(IAER,IWAV)/EXTC(IAER,IWAV)
                          WAVCOR=(VX0(71)/VX0(IWAV))**ASTMX
                          IF(IAER.EQ.1)THEN
                              IF(LDASTM(1))THEN

!                                 IF(LDASTM), TREAT ANGSTROM EXPONENT AS
!                                 A DELTA VALUE & IGNORE THE COEFFICIENT
                                  EXTC(1,IWAV)=EXTC(1,IWAV)*WAVCOR
                              ELSE

!                                 FULL ANGSTROM LAW
                                  EXTC(1,IWAV)=ASTMO+ASTMC*WAVCOR
                              ENDIF
                          ELSE
                              IF(LDASTM(2))THEN
                                  EXTC(2,IWAV)=EXTC(2,IWAV)*WAVCOR
                              ELSE

!                                 FULL ANGSTROM LAW
                                  EXTC(2,IWAV)=ASTMO+ASTMC*WAVCOR
                              ENDIF
                          ENDIF
                          IF(EXTC(IAER,IWAV).LT.0.)EXTC(IAER,IWAV)=0.
                          ABSC(IAER,IWAV)=COALB*EXTC(IAER,IWAV)
                      ENDIF
                  ENDDO
              ENDIF
          ENDDO
      ENDIF

      IF(LSSALB.GE.0)THEN

!         DETERMINE AEROSOL (1=BOUNDARY LAYER, 2=TROPOSPHERIC) SINGLE
!         SCATTERING ALBEDO SCALING FACTORS AT SPECTRAL GRID END POINTS.
          ASCLMN(1)=1.
          ASCLMN(2)=1.
          ASCLMX(1)=1.
          ASCLMX(2)=1.
          IF(AWAVLN(0).GE.VX0(1))THEN
              IWAV=1
              DO IWAVP1=2,NWAVLN
                  IF(AWAVLN(0).LE.VX0(IWAVP1))THEN
                      FACTOR                                            &
     &                  =(AWAVLN(0)-VX0(IWAV))/(VX0(IWAVP1)-VX0(IWAV))
                      EXTV=EXTC(1,IWAV)                                 &
     &                  +FACTOR*(EXTC(1,IWAVP1)-EXTC(1,IWAV))
                      IF(EXTV.GT.0.)THEN
                          ABSV=ABSC(1,IWAV)                             &
     &                      +FACTOR*(ABSC(1,IWAVP1)-ABSC(1,IWAV))
                          SSALB=(EXTV-ABSV)/EXTV
                          IF(SSALB.GT.0.)ASCLMN(1)=ASSALB(0)/SSALB
                      ENDIF
                      EXTV=EXTC(2,IWAV)                                 &
     &                  +FACTOR*(EXTC(2,IWAVP1)-EXTC(2,IWAV))
                      IF(EXTV.GT.0.)THEN
                          ABSV=ABSC(2,IWAV)                             &
     &                      +FACTOR*(ABSC(2,IWAVP1)-ABSC(2,IWAV))
                          SSALB=(EXTV-ABSV)/EXTV
                          IF(SSALB.GT.0.)ASCLMN(2)=ASSALB(0)/SSALB
                      ENDIF
                      GOTO 20
                  ENDIF
                  IWAV=IWAVP1
              ENDDO
              RETURN
          ENDIF
   20     CONTINUE
          IF(AWAVLN(LSSALB).GE.VX0(1))THEN
              IWAV=1
              DO IWAVP1=2,NWAVLN
                  IF(AWAVLN(LSSALB).LE.VX0(IWAVP1))THEN                 &
                      FACTOR=(AWAVLN(LSSALB)-VX0(IWAV))                 &
     &                  /(VX0(IWAVP1)-VX0(IWAV))
                      EXTV=EXTC(1,IWAV)                                 &
     &                  +FACTOR*(EXTC(1,IWAVP1)-EXTC(1,IWAV))
                      IF(EXTV.GT.0.)THEN
                          ABSV=ABSC(1,IWAV)                             &
     &                  +FACTOR*(ABSC(1,IWAVP1)-ABSC(1,IWAV))
                          SSALB=(EXTV-ABSV)/EXTV
                          IF(SSALB.GT.0.)ASCLMX(1)=ASSALB(LSSALB)/SSALB
                      ENDIF
                      EXTV=EXTC(2,IWAV)                                 &
     &                  +FACTOR*(EXTC(2,IWAVP1)-EXTC(2,IWAV))
                      IF(EXTV.GT.0.)THEN
                          ABSV=ABSC(2,IWAV)                             &
     &                  +FACTOR*(ABSC(2,IWAVP1)-ABSC(2,IWAV))
                          SSALB=(EXTV-ABSV)/EXTV
                          IF(SSALB.GT.0.)ASCLMX(2)=ASSALB(LSSALB)/SSALB
                      ENDIF
                      RETURN
                  ENDIF
                  IWAV=IWAVP1
              ENDDO
              RETURN
          ENDIF
      ENDIF

      RETURN
      END

      SUBROUTINE RSSALB(NSSALB,RHASYM)

!     RSSALB READS IN A BOUNDARY LAYER + TROPOSPHERIC
!     AEROSOL SPECTRAL SINGLE SCATTERING ALBEDO GRID AND
!     SETS UP INTERPOLATION. IF NSSALB IS LESS THAN ZERO, A
!     CONSTANT COALBEDO SCALING FACTOR IS READ IN INSTEAD.
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     INPUT ARGUMENTS:
!       NSSALB   NUMBER OF AEROSOL ALBEDO SPECTRAL GRID POINTS.
!       RHASYM   RELATIVE HUMIDITY TO DEFINE AEROSOL ASYMMETRY FACTORS.
      INTEGER NSSALB
      REAL RHASYM

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /AERALB/
!       LSSALB   INDEX OF FINAL AEROSOL ALBEDO SPECTRAL GRID POINT.
!       ACOALB   AEROSOL SINGLE SCATTERING COALBEDO SCALING FACTOR.
!       AWAVLN   AEROSOL ALBEDO SPECTRAL GRID POINTS [MICRONS].
!       ASSALB   BOUNDARY LAYER AND TROPOSPHERE AEROSOL SINGLE
!                SCATTERING ALBEDOS.
!       AASTRM   WAVELENGTH INTERPOLATION EXPONENTS (EQUIVALENT TO AN
!                ANGSTROM COEFFICIENT FOR SINGLE SCATTERING ALBEDO)
!       ASCLMN   AEROSOL (1=BOUNDARY LAYER, 2=TROPOSPHERIC)
!                SINGLE SCATTERING ALBEDO SCALING FACTORS
!                FOR WAVELENGTHS BELOW AWAVLN(0).
!       ASCLMX   AEROSOL (1=BOUNDARY LAYER, 2=TROPOSPHERIC)
!                SINGLE SCATTERING ALBEDO SCALING FACTORS
!                FOR WAVELENGTHS ABOVE AWAVLN(LSSALB).
      INTEGER LSSALB
      REAL ACOALB,AWAVLN,ASSALB,AASTRM,ASCLMN,ASCLMX
      COMMON/AERALB/LSSALB,ACOALB,AWAVLN(0:MSSALB),ASSALB(0:MSSALB),    &
     &     ASCLMN(2),AASTRM(MSSALB),ASCLMX(2)

!     LOCAL VARIABLES:
!       IOTEST   OUTPUT OF IOSTAT I/O STATUS SPECIFIER.
!       NEXT     INDEX FOR AEROSOL SINGLE SCATTERING ALBEDO GRID.
!       LAST     PREVIOUS AEROSOL SINGLE SCATTERING ALBEDO GRID INDEX.
      INTEGER IOTEST,NEXT,LAST

!     DEFINE LAST (FINAL) INDEX:
      LSSALB=NSSALB-1
      RHASYM=0.
      IF(LSSALB.EQ.-1)RETURN
      IF(LSSALB.LT.-1)THEN

!         READ IN AEROSOL SINGLE SCATTERING COALBEDO FACTOR AND
!         RELATIVE HUMIDITY USED TO DEFINE ASYMMETRY FACTORS:
          READ(IRD,'(2F10.5)',IOSTAT=IOTEST)ACOALB,RHASYM
          IF(IOTEST.NE.0)STOP                                           &
     &      'error reading aerosol coalbedo factor and asymmetry rh.'
          IF(ACOALB.LT.0.)ACOALB=0.
          IF(RHASYM.LT.0.)THEN
              RHASYM=0.
          ELSEIF(RHASYM.GT.99.)THEN
              RHASYM=99.
          ENDIF
          WRITE(IPR,'(/(F10.5,3X,A))')                                  &
     &      ACOALB,'aerosol single scattering coalbedo factor',         &
     &      RHASYM,'relative humidity for aerosol asymmetry factors'
          RETURN
      ENDIF
      IF(LSSALB.GT.MSSALB)STOP 'parameter mssalb must be increased.'

!     READ THE CARD 1B INPUTS:
      READ(IRD,'((8F10.0))',IOSTAT=IOTEST)                              &
     &  (AWAVLN(NEXT),ASSALB(NEXT),NEXT=0,LSSALB)
      IF(IOTEST.NE.0)                                                   &
     &  STOP 'error reading aerosol single scattering albedos.'
      WRITE(IPR,'(//A,/2(/A),/(4(F10.5,F10.7)))')                       &
     &  ' aerosol single scattering albedo table:',                     &
     &  ' wavln(um) ss albedo wavln(um) ss albedo'                      &
     &  //' wavln(um) ss albedo wavln(um) ss albedo',                   &
     &  ' --------- --------- --------- ---------'                      &
     &  //' --------- --------- --------- ---------',                   &
     &  (AWAVLN(NEXT),ASSALB(NEXT),NEXT=0,LSSALB)
      IF(AWAVLN(0).LE.0.)                                               &
     &  STOP 'first aerosol wavelength not greater than zero.'
      IF(ASSALB(0).LE.0. .OR. ASSALB(0).GT.1.)                          &
     &  STOP 'first aerosol coalbedo is out-of-range.'
      LAST=0
      DO NEXT=1,LSSALB
          IF(AWAVLN(NEXT).LE.AWAVLN(LAST))                              &
     &      STOP 'aerosol wavelength grid is out of order.'
          IF(ASSALB(NEXT).LE.0. .OR. ASSALB(NEXT).GT.1.)                &
     &      STOP 'an aerosol coalbedo is out-of-range.'
          AASTRM(NEXT)=LOG(ASSALB(LAST)/ASSALB(NEXT))                   &
     &      /LOG(AWAVLN(NEXT)/AWAVLN(LAST))
          LAST=NEXT
      ENDDO
      RETURN
      END
