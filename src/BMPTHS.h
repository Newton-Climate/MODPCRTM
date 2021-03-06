
!     /BMPTHS/
!       PTLOSS   PT SCALING FOR SOLAR ILLUMINATION PATHS TO LOS[ATM].
!       PTMSS    PT SCALING FOR SOLAR PATH TO VERTICAL MS PATHS [ATM].
!       P2LOSS   SUN TO LOS PATH P-SQUARED INTERPOLATION FRACTION.
!       P2MSS    SUN TO VERTICAL PATH P-SQUARED INTERPOLATION FRACTION.
!       T5LOSS   SQRT(T/273.15K) FOR SOLAR ILLUMINATION PATHS TO LOS.
!       T5MSS    SQRT(T/273.15K) FOR SOLAR PATHS TO VERTICAL MS PATH.
!       JTLOSS   BAND MODEL TEMP UPPER INDEX FOR SOLAR PATHS TO LOS.
!       JTMSS    BAND MODEL TEMP UPPER INDEX FOR SOLAR PATHS TO MS PATH.
!       KTLOSS   X-SECTION TEMP UPPER INDEX FOR SOLAR PATHS TO LOS.
!       KTMSS    X-SECTION TEMP UPPER INDEX FOR SOLAR PATHS TO MS PATH.
!       FTLOSS   TEMPERATURE INTERPOLATION FRACTION USED WITH JTLOSS.
!       FTMSS    TEMPERATURE INTERPOLATION FRACTION USED WITH JTMSS.
!       GTLOSS   TEMPERATURE INTERPOLATION FRACTION USED WITH KTLOSS.
!       GTMSS    TEMPERATURE INTERPOLATION FRACTION USED WITH KTMSS.
      INTEGER JTLOSS,JTMSS,KTLOSS,KTMSS
      REAL PTLOSS,PTMSS,P2LOSS,P2MSS,T5LOSS,T5MSS,                      &
     &  FTLOSS,FTMSS,GTLOSS,GTMSS
      COMMON/BMPTHS/PTLOSS(MMOLT,LAYTWO,MLOS),PTMSS(MMOLT,LAYDIM),      &
     &  P2LOSS(MMOLT,LAYTWO,MLOS),P2MSS(MMOLT,LAYDIM),                  &
     &  T5LOSS(MMOLT,LAYTWO,MLOS),T5MSS(MMOLT,LAYDIM),                  &
     &  JTLOSS(MMOLT,LAYTWO,MLOS),JTMSS(MMOLT,LAYDIM),                  &
     &  KTLOSS(MMOLT,LAYTWO,MLOS),KTMSS(MMOLT,LAYDIM),                  &
     &  FTLOSS(MMOLT,LAYTWO,MLOS),FTMSS(MMOLT,LAYDIM),                  &
     &  GTLOSS(MMOLT,LAYTWO,MLOS),GTMSS(MMOLT,LAYDIM)
      SAVE /BMPTHS/
