      SUBROUTINE CIRR18(ICLD,LCIRZ)
      IMPLICIT NONE
      INTEGER ICLD
      LOGICAL LCIRZ
!*********************************************************************
!*  ROUTINE TO SET CTHIK CALT CEXT  FOR  CIRRUS CLOUDS 18 19        **
!*  INPUTS]                                                         **
!*           CTHIK    -  CIRRUS THICKNESS (KM)                      **
!*                       0 = USE THICKNESS STATISTICS               **
!*                       .NE. 0 = USER DEFINES THICKNESS            **
!*                                                                  **
!*           CALT     -  CIRRUS BASE ALTITUDE (KM)                  **
!*                       0 = USE CALCULATED VALUE                   **
!*                       .NE. 0 = USER DEFINES BASE ALTITUDE        **
!*                                                                  **
!*           ICLD     -  CIRRUS PRESENCE FLAG                       **
!*                       0 = NO CIRRUS                              **
!*                       18  19 = USE CIRRUS PROFILE                **
!*                                                                  **
!*           MODEL    -  ATMOSPHERIC MODEL                          **
!*                       1-5  AS IN MAIN PROGRAM                    **
!*                       MODEL = 0,6,7 NOT USED SET TO 2            **
!*                                                                  **
!*  OUTPUTS]                                                        **
!*         CTHIK        -  CIRRUS THICKNESS (KM)                    **
!*         CALT         -  CIRRUS BASE ALTITUDE (KM)                **
!          CEXT IS THE EXTINCTION COEFFIENT(KM-1) AT 0.55
!               DEFAULT VALUE 0.14*CTHIK
!*                                                                  **
!*********************************************************************

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

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

!     /CARD2A/
      INTEGER NCRALT,NCRSPC
      DOUBLE PRECISION CTHIK,CALT
      REAL CEXT,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP
      COMMON/CARD2A/CTHIK,CALT,CEXT,NCRALT,NCRSPC,CWAVLN,CCOLWD,        &
     &  CCOLIP,CHUMID,ASYMWD,ASYMIP
      INCLUDE 'IFIL.h'

!     LOCAL VARIABLES:
      INTEGER MDL

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD
      REAL CAMEAN(5)
      DATA CAMEAN/11.,10.,8.,7.,5./
      MDL = MODEL

!     INITIAL LCIRZ TO TRUE.
      LCIRZ=.TRUE.

!  CHECK IF USER WANTS TO USE A THICKNESS VALUE HE PROVIDES
!  DEFAULTED MEAN CIRRUS THICKNESS IS 1.0KM  OR 0.2 KM.

      IF(CTHIK.GT.0.D0)GOTO 25
      IF(ICLD.EQ.18)CTHIK=1.D0
      IF(ICLD.EQ.19)CTHIK=.2D0
25    IF(CEXT.EQ.0.)CEXT=.14*SNGL(CTHIK)

!  BASE HEIGHT CALCULATIONS

      IF ( MODEL .LT. 1  .OR.  MODEL .GT. 5 ) MDL = 2

      IF(CALT.LE.0.D0)CALT=DBLE(CAMEAN(MDL))

      IF(LJMASS)RETURN
      IF(ICLD.EQ.18)THEN
          WRITE(IPR,'(15X,A)')                                          &
     &      'CIRRUS ATTENUATION INCLUDED   (STANDARD CIRRUS)'
      ELSEIF(ICLD.EQ.19)THEN
          WRITE(IPR,'(15X,A)')                                          &
     &      'CIRRUS ATTENUATION INCLUDED   (THIN CIRRUS)'
      ENDIF
      WRITE(IPR,'((15X,A,F10.3,A))')                                    &
     &  'CIRRUS THICKNESS          ',CTHIK,' KM',                       &
     &  'CIRRUS BASE ALTITUDE      ',CALT,' KM',                        &
     &  'CIRRUS PROFILE EXTINCTION ',CEXT,' KM-1'

!     END OF CIRRUS MODEL SET UP
      RETURN
      END
