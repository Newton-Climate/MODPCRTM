      BLOCK DATA TITLE

!     TITLE INFORMATION

!     COMMONS:

!     /TITL/
!       HHAZE    AEROSOL MODEL NAMES.
!       HSEASN   SEASON CHARACTER STRING.
!       HVULCN   VOLCANIC AEROSOL MODEL NAME.
!       HMET     METEOROLOGICAL REGION NAME.
!       HMODEL   ATMOSPHERIC PROFILE MODEL NAME.
      CHARACTER HHAZE(16)*20,HSEASN(2)*20,HVULCN(8)*20,                 &
     &  HMET(2)*20,HMODEL(0:8)*20
      COMMON/TITL/HHAZE,HSEASN,HVULCN,HMET,HMODEL

!     /VSBD/
!       VSB      SURFACE METEOROLOGICAL RANGE (VISIBILITY) GRID [KM].
      REAL VSB
      COMMON/VSBD/VSB(10)
      SAVE /TITL/,/VSBD/
      DATA HHAZE/              'RURAL               ',                  &
     &  'RURAL               ','NAVY_MARITIME       ',                  &
     &  'MARITIME            ','URBAN               ',                  &
     &  'TROPOSPHERIC        ','USER_DEFINED        ',                  &
     &  'FOG1 (ADVECTTION)   ','FOG2 (RADIATI0N)    ',                  &
     &  'DESERT_AEROSOL      ','BCKGD_STRATOSPHERIC ',                  &
     &  'AGED_VOLCANIC       ','FRESH_VOLCANIC      ',                  &
     &  'AGED_VOLCANIC       ','FRESH_VOLCANIC      ',                  &
     &  'METEORIC DUST       '/
      DATA HSEASN/                                                      &
     &  'SPRING-SUMMER       ','FALL-WINTER         ' /
      DATA HVULCN/                                                      &
     &  'BCKGD_STRATOSPHERIC ','MODERATE_VOLCANIC   ',                  &
     &  'HIGH     VOLCANIC   ','HIGH     VOLCANIC   ',                  &
     &  'MODERATE_VOLCANIC   ','MODERATE_VOLCANIC   ',                  &
     &  'HIGH     VOLCANIC   ','EXTREME  VOLCANIC   '/
      DATA HMET/                                                        &
     &  'NORMAL              ','TRANSITION          '/
      DATA HMODEL/             'MODEL=0 HORIZONTAL  ',                  &
     &  'TROPICAL MODEL      ','MID-LATITUDE SUMMER ',                  &
     &  'MID-LATITUDE WINTER ','SUB-ARCTIC   SUMMER ',                  &
     &  'SUB-ARCTIC   WINTER ','1976 U S STANDARD   ',                  &
     &  'USER PROFILES       ','HYDROSTATIC EQUATION'/
      DATA VSB/23.,5.,0.,23.,5.,50.,23.,0.2,0.5,0./
      END
