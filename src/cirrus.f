      SUBROUTINE CIRRUS(CTHIK,CALT,CPROB,CEXT)
!*********************************************************************
!*  ROUTINE TO GENERATE ALTITUDE PROFILES OF CIRRUS DENSITY         **
!*  PROGRAMMED BY   M.J. POST                                       **
!*                  R.A. RICHTER        NOAA/WPL                    **
!*                                      BOULDER, COLORADO           **
!*                                      01/27/1981                  **
!*                                                                  **
!*  INPUTS!                                                         **
!*           CTHIK    -  CIRRUS THICKNESS (KM)                      **
!*                       0 = USE THICKNESS STATISTICS               **
!*                       .NE. 0 = USER DEFINES THICKNESS            **
!*                                                                  **
!*           CALT     -  CIRRUS BASE ALTITUDE (KM)                  **
!*                       0 = USE CALCULATED VALUE                   **
!*                       .NE. 0 = USER DEFINES BASE ALTITUDE        **
!*                                                                  **
!*           ICIR     -  CIRRUS PRESENCE FLAG                       **
!*                       0 = NO CIRRUS                              **
!*                       .NE. 0 = USE CIRRUS PROFILE                **
!*                                                                  **
!*           MODEL    -  ATMOSPHERIC MODEL                          **
!*                       1-5  AS IN MAIN PROGRAM                    **
!*                       MODEL = 0,6,7 NOT USED SET TO 2            **
!*                                                                  **
!*  OUTPUTS!                                                        **
!*         CTHIK        -  CIRRUS THICKNESS (KM)                    **
!*         CALT         -  CIRRUS BASE ALTITUDE (KM)                **
!*         DENSTY(16,I) -  ARRAY, ALTITUDE PROFILE OF CIRRUS DENSITY**
!*         CPROB        -  CIRRUS PROBABILITY                       **
!*                                                                  **
!*********************************************************************
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     ARGUMENTS:
      DOUBLE PRECISION CTHIK,CALT
      REAL CPROB,CEXT

!     COMMONS:

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

!     /DEN/
!       DENSTY   PROFILE LEVEL DENSITIES [ATM CM / KM FOR MOST SPECIES].
      REAL DENSTY
      COMMON/DEN/DENSTY(0:MEXTXY,1:LAYDIM)

!     LOCAL VARIABLES:
      INTEGER MDL,I,IBOT,ITOP,IHGT1,IHGT2,IML
      REAL DENCIR
      DOUBLE PRECISION ZMIN,ZMAX,TOP

!     DATA:
      DOUBLE PRECISION CAMEAN(5)
      REAL PTAB(5)
      DATA CAMEAN/1.1D1,1.D1,8.D0,7.D0,5.D0/,PTAB/.8,.4,.5,.45,.4/

!     SET CIRRUS PROBABILITY AND PROFILE TO ALL ZEROES
      CPROB=0.
      MDL=MODEL
      DO I=1,ML
          DENSTY(16,I)=0.
      ENDDO

!     CHECK IF USER WANTS TO USE A THICKNESS VALUE HE PROVIDES OR USE A
!     MEAN THICKNESS. DEFAULTED MEAN CIRRUS THICKNESS IS 1.0 KM.
      IF(CTHIK.LE.0.D0)CTHIK=1.D0

!     DENCIR IS CIRRUS DENSITY IN KM-1
      IF(CEXT.GT.0.)THEN
          DENCIR=CEXT/2
      ELSE
          DENCIR=.07*SNGL(CTHIK)
      ENDIF

!     BASE HEIGHT CALCULATIONS
      IF(MODEL.LT.1 .OR. MODEL.GT.5)MDL=2
      CPROB=100*PTAB(MDL)
      IF(CALT.LE.0.)CALT=CAMEAN(MDL)

!     PUT CIRRUS DENSITY IN CORRECT ALTITUDE BINS. IF MODEL=7,
!     INTERPOLATE EH(16,I) FOR NON-STANDARD ALTITUDE BOUNDARIES.
      IF(MODEL.LT.7)THEN
          IBOT=INT(CALT)
          ITOP=INT(CALT+CTHIK)
          DO I=2,16
              IF(I.GE.IBOT .AND. I.LE.ITOP)DENSTY(16,I+1)=DENCIR
          ENDDO

!         ADJUST FIRST AND LAST CIRRUS LEVEL IF CLOUD DOES NOT ENTIRELY
!         FILL EACH LEVEL.
          IHGT1=INT(CALT)
          IHGT2=INT(CALT+CTHIK)
          IF(IHGT1.EQ.IHGT2)THEN
              DENSTY(16,IHGT1+1)=DENSTY(16,IHGT1+1)*REAL(CTHIK)
              RETURN
          ENDIF
          DENSTY(16,IHGT1+1)=DENSTY(16,IHGT1+1)*REAL(1-CALT+IHGT1)
          DENSTY(16,IHGT2+1)=DENSTY(16,IHGT2+1)*REAL(CALT+CTHIK-IHGT2)
          RETURN
      ENDIF

!     INTERPOLATE DENSTY(16,I) FOR USER SUPPLIED ALTITUDE BOUNDARIES
      TOP=CALT+CTHIK
      IF(TOP.LT.ZM(1) .OR. CALT.GT.ZM(ML))RETURN
      IML=ML-1
      DO I=1,IML
          ZMIN=ZM(I)
          ZMAX=ZM(I+1)
          IF(CALT.LE.ZMIN .AND. TOP.GE.ZMAX)DENSTY(16,I)=DENCIR
          IF(CALT.GE.ZMIN .AND. TOP.LT.ZMAX)                            &
     &      DENSTY(16,I)=DENCIR*REAL(CTHIK/(ZMAX-ZMIN))
          IF(CALT.GE.ZMIN .AND. TOP.GE.ZMAX .AND. CALT.LT.ZMAX)         &
     &      DENSTY(16,I)=DENCIR*REAL((ZMAX-CALT)/(ZMAX-ZMIN))
          IF(CALT.LT.ZMIN .AND. TOP.LE.ZMAX .AND. TOP.GT.ZMIN)          &
     &      DENSTY(16,I)=DENCIR*REAL((TOP-ZMIN)/(ZMAX-ZMIN))
      ENDDO
      RETURN
      END
