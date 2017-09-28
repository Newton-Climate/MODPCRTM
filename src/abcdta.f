      SUBROUTINE ABCDTA(V)

!     LOWTRAN BAND MODEL DATA:
!        1   H2O (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)
!        2   CO2 (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)
!        3   O3  (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)
!        4   N2O (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)
!        5   CO  (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)
!        6   CH4 (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)
!        7   O2  (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)
!        8   NO  (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)
!        9   SO2 (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)
!       10   NO2 (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)
!       11   NH3 (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)
      IMPLICIT NONE

!     INPUT ARGUMENTS:
!       V        SPECTRAL FREQUENCY [CM-1].
      REAL V

!     COMMONS:
      INTEGER IBND
      REAL AA,BB,CC,QA,CPS
      COMMON/AABBCC/AA(11),BB(11),CC(11),IBND(11),QA(11),CPS(11)

!     LOCAL VARIABLES:
      INTEGER IW,IBAND

!     DATA:
      REAL AH2O(14),AAH2O(14),BBH2O(14),CCH2O(14),ACO2(10),             &
     &  AACO2(10),BBCO2(10),CCCO2(10),AO3(5),AAO3(5),BBO3(5),           &
     &  CCO3(5),AN2O(11),AAN2O(11),BBN2O(11),CCN2O(11),ACO(3),          &
     &  AACO(3),BBCO(3),CCCO(3),ACH4(4),AACH4(4),BBCH4(4),              &
     &  CCCH4(4),AO2(6),AAO2(6),BBO2(6),CCO2(6),ANO,AANO,BBNO,          &
     &  CCNO,ASO2(4),AASO2(4),BBSO2(4),CCSO2(4),ANO2(3),AANO2(3),       &
     &  BBNO2(3),CCNO2(3),ANH3(2),AANH3(2),BBNH3(2),CCNH3(2)
      DATA AH2O,AAH2O,BBH2O,CCH2O/.5274,.5299,.5416,.5479,.5495, .5464, &
     &   .5454,  .5474,  .5579,  .5621,  .5847,  .6076,  .6508,  .6570, &
     &                 .219312,.216415,.206349,.196196,.194540,.198500, &
     & .198500,.196196,.184148,.179360,.154120,.130095,.091341,.086549, &
     &                 .334884,.336904,.343272,.348610,.349810,.347498, &
     & .347498,.348610,.353429,.354864,.357640,.352497,.327526,.322898, &
     &                 21.8352,21.9588,22.4234,22.9517,23.0750,22.8262, &
     & 22.8262,22.9517,23.6654,23.9774,25.9207,28.2957,33.3998,34.1575/
      DATA ACO2,AACO2,BBCO2,CCCO2/                                      &
     &  .6176,  .6810,  .6033,  .6146,  .6513,  .6050,  .6160,  3*.7070,&
     &.120300,.069728,.134448,.123189,.090948,.132717,.121835,3*.054348,&
     &.348172,.303510,.354002,.349583,.327160,.353435,.348936,3*.280674,&
     &29.4277,37.0842,27.8241,29.0834,33.4608,28.0093,29.2436,3*40.1951/
      DATA AO3,AAO3,BBO3,CCO3/.8559,   .7593,   .7819,   .9175,   .7703,&
     &                      .006712, .030870, .023278, .000458, .027004,&
     &                      .138026, .231722, .209952, .078492, .221153,&
     &                      55.6442, 46.1189, 48.5155, 60.7802, 47.2982/
      DATA AN2O,AAN2O,BBN2O,CCN2O/                                      &
     &    .8997,   5*.7201,   5*.6933, .001679, 5*.047599, 5*.062106,   &
     &  .095621, 5*.268696, 5*.292891, 59.3660, 5*41.7251, 5*38.5667/
      DATA ACO,AACO,BBCO,CCCO/                    .6397,   2*.6133,     &
     &  .100401, 2*.124454, .335296, 2*.350165, 32.0496, 2*28.9354/
      DATA ACH4,AACH4,BBCH4,CCCH4/4*.5844,4*.154447,4*.357657,4*25.8920/
      DATA AO2,AAO2,BBO2,CCO2/                    .6011,   5*.5641,     &
     &  .136706, 5*.177087, .354683, 5*.355447, 27.5869, 5*24.1314/
      DATA ANO,AANO,BBNO,CCNO/.6613, .083336, .319585, 34.6834/
      DATA ASO2,AASO2,BBSO2,CCSO2/                .8907,   3*.8466,     &
     &  .002468, 3*.008192, .104307, 3*.147065, 58.6298, 3*54.8078/
      DATA ANO2,AANO2,BBNO2,CCNO2/3*.7249,3*.045281,3*.264248,3*42.2784/
      DATA ANH3,AANH3,BBNH3,CCNH3/            .4704,   .6035,           &
     &  .285772, .134244, .269839, .353937, 19.9507, 27.8458/
      SAVE AH2O,AAH2O,BBH2O,CCH2O,ACO2,AACO2,BBCO2,CCCO2,AO3,AAO3,BBO3, &
     &  CCO3,AN2O,AAN2O,BBN2O,CCN2O,ACO,AACO,BBCO,CCCO,ACH4,AACH4,      &
     &  BBCH4,CCCH4,AO2,AAO2,BBO2,CCO2,ANO,AANO,BBNO,CCNO,ASO2,AASO2,   &
     &  BBSO2,CCSO2,ANO2,AANO2,BBNO2,CCNO2,ANH3,AANH3,BBNH3,CCNH3

!  ---H2O
      IW=-1
      IF(V.GE.     0..AND.V.LE.   345.) IW=17
      IF(V.GE.   350..AND.V.LE.  1000.) IW=18
      IF(V.GE.  1005..AND.V.LE.  1640.) IW=19
      IF(V.GE.  1645..AND.V.LE.  2530.) IW=20
      IF(V.GE.  2535..AND.V.LE.  3420.) IW=21
      IF(V.GE.  3425..AND.V.LE.  4310.) IW=22
      IF(V.GE.  4315..AND.V.LE.  6150.) IW=23
      IF(V.GE.  6155..AND.V.LE.  8000.) IW=24
      IF(V.GE.  8005..AND.V.LE.  9615.) IW=25
      IF(V.GE.  9620..AND.V.LE. 11540.) IW=26
      IF(V.GE. 11545..AND.V.LE. 13070.) IW=27
      IF(V.GE. 13075..AND.V.LE. 14860.) IW=28
      IF(V.GE. 14865..AND.V.LE. 16045.) IW=29
      IF(V.GE. 16340..AND.V.LE. 17860.) IW=30
      IBAND=IW - 16
      IBND(1)=IW
      IF(IW .GT.  0) THEN
           QA(1)= AH2O(IBAND)
           AA(1) =AAH2O(IBAND)
           BB(1) =BBH2O(IBAND)
           CC(1) =CCH2O(IBAND)
      ENDIF
!  ---O3
      IW=-1
      IF (V .GE.     0. .AND. V .LE.   200.)  IW=31
      IF (V .GE.   515. .AND. V .LE.  1275.)  IW=32
      IF (V .GE.  1630. .AND. V .LE.  2295.)  IW=33
      IF (V .GE.  2670. .AND. V .LE.  2845.)  IW=34
      IF (V .GE.  2850. .AND. V .LE.  3260.)  IW=35
      IBAND     =IW - 30
      IBND(3)=IW
      IF(IW .GT.  0) THEN
           QA(3)=AO3(IBAND)
           AA(3)=AAO3(IBAND)
           BB(3)=BBO3(IBAND)
           CC(3)=CCO3(IBAND)
      ENDIF
!  ---CO2
      IW=-1
      IF (V .GE.   425. .AND. V .LE.   835.)  IW=36
      IF (V .GE.   840. .AND. V .LE.  1440.)  IW=37
      IF (V .GE.  1805. .AND. V .LE.  2855.)  IW=38
      IF (V .GE.  3070. .AND. V .LE.  3755.)  IW=39
      IF (V .GE.  3760. .AND. V .LE.  4065.)  IW=40
      IF (V .GE.  4530. .AND. V .LE.  5380.)  IW=41
      IF (V .GE.  5905. .AND. V .LE.  7025.)  IW=42
      IF((V .GE.  7395. .AND. V .LE.  7785.) .OR.                       &
     &   (V .GE.  8030. .AND. V .LE.  8335.) .OR.                       &
     &   (V .GE.  9340. .AND. V .LE.  9670.)) IW=43
      IBAND=IW - 35
      IBND(2)=IW
      IF(IW .GT.  0) THEN
           QA(2)=ACO2(IBAND)
           AA(2)=AACO2(IBAND)
           BB(2)=BBCO2(IBAND)
           CC(2)=CCCO2(IBAND)
      ENDIF
!  ---CO
      IW=-1
      IF (V .GE.     0. .AND. V .LE.   175.) IW=44
      IF((V .GE.  1940. .AND. V .LE.  2285.) .OR.                       &
     &   (V .GE.  4040. .AND. V .LE.  4370.)) IW=45
      IBAND=IW - 43
      IBND(5)=IW
      IF(IW .GT.  0) THEN
           QA(5)=ACO(IBAND)
           AA(5)=AACO(IBAND)
           BB(5)=BBCO(IBAND)
           CC(5)=CCCO(IBAND)
      ENDIF
!  ---CH4
      IW=-1
      IF((V .GE.  1065. .AND. V .LE.  1775.) .OR.                       &
     &   (V .GE.  2345. .AND. V .LE.  3230.) .OR.                       &
     &   (V .GE.  4110. .AND. V .LE.  4690.) .OR.                       &
     &   (V .GE.  5865. .AND. V .LE.  6135.))IW=46
      IBAND=IW - 45
      IBND(6)=IW
      IF(IW .GT.  0) THEN
           QA(6)=ACH4(IBAND)
           AA(6)=AACH4(IBAND)
           BB(6)=BBCH4(IBAND)
           CC(6)=CCCH4(IBAND)
      ENDIF
!  ---N2O
      IW=-1
      IF (V .GE.     0. .AND. V .LE.   120.)  IW=47
      IF((V .GE.   490. .AND. V .LE.   775.) .OR.                       &
     &   (V .GE.   865. .AND. V .LE.   995.) .OR.                       &
     &   (V .GE.  1065. .AND. V .LE.  1385.) .OR.                       &
     &   (V .GE.  1545. .AND. V .LE.  2040.) .OR.                       &
     &   (V .GE.  2090. .AND. V .LE.  2655.)) IW=48
      IF((V .GE.  2705. .AND. V .LE.  2865.) .OR.                       &
     &   (V .GE.  3245. .AND. V .LE.  3925.) .OR.                       &
     &   (V .GE.  4260. .AND. V .LE.  4470.) .OR.                       &
     &   (V .GE.  4540. .AND. V .LE.  4785.) .OR.                       &
     &   (V .GE.  4910. .AND. V .LE.  5165.)) IW=49
      IBAND=IW - 46
      IBND(4)=IW

      IF(IW .EQ. 49)IBAND=7

!     THIS CORRECTION IS ONLY FOR N2O AS CURRENTLY WRITTEN
      IF(IW .GT.  0) THEN
           QA(4)=AN2O(IBAND)
           AA(4)=AAN2O(IBAND)
           BB(4)=BBN2O(IBAND)
           CC(4)=CCN2O(IBAND)
      ENDIF
!  ---O2
      IW=-1
      IF (V .GE.     0. .AND. V .LE.   265.)  IW=50
      IF((V .GE.  7650. .AND. V .LE.  8080.) .OR.                       &
     &   (V .GE.  9235. .AND. V .LE.  9490.) .OR.                       &
     &   (V .GE. 12850. .AND. V .LE. 13220.) .OR.                       &
     &   (V .GE. 14300. .AND. V .LE. 14600.) .OR.                       &
     &   (V .GE. 15695. .AND. V .LE. 15955.)) IW=51
       IF(V .GE. 49600. .AND. V. LE. 52710.)  IW=51
      IBAND=IW - 49
      IBND(7)=IW
      IF(IW .GT.  0) THEN
           QA(7)=AO2(IBAND)
           IF(V .GE. 49600. .AND. V. LE. 52710.)  QA(7) =.4704
           AA(7)=AAO2(IBAND)
           BB(7)=BBO2(IBAND)
           CC(7)=CCO2(IBAND)
      ENDIF
!  ---NH3
      IW=-1
      IF (V .GE.     0. .AND. V .LE.   385.)  IW=52
      IF (V .GE.   390. .AND. V .LE.  2150.)  IW=53
      IBAND=IW - 51
      IBND(11)=IW
      IF(IW .GT.  0) THEN
           QA(11)=ANH3(IBAND)
           AA(11)=AANH3(IBAND)
           BB(11)=BBNH3(IBAND)
           CC(11)=CCNH3(IBAND)
      ENDIF
!  ---NO
      IW=-1
      IF (V .GE.  1700. .AND. V .LE.  2005.) IW =54
      IBAND=IW - 53
      IBND(8)=IW
      IF(IW .GT.  0) THEN
           QA(8)=ANO
           AA(8)=AANO
           BB(8)=BBNO
           CC(8)=CCNO
      ENDIF
!  ---NO2
      IW=-1
      IF((V .GE.   580. .AND. V .LE.   925.) .OR.                       &
     &   (V .GE.  1515. .AND. V .LE.  1695.) .OR.                       &
     &   (V .GE.  2800. .AND. V .LE.  2970.)) IW=55
      IBAND=IW - 54
      IBND(10)=IW
      IF(IW .GT.  0) THEN
           QA(10)=ANO2(IBAND)
           AA(10)=AANO2(IBAND)
           BB(10)=BBNO2(IBAND)
           CC(10)=CCNO2(IBAND)
      ENDIF
!  ---SO2
      IW=-1
      IF (V .GE.     0. .AND. V .LE.   185.)  IW=56
      IF((V .GE.   400. .AND. V .LE.   650.) .OR.                       &
     &   (V .GE.   950. .AND. V .LE.  1460.) .OR.                       &
     &   (V .GE.  2415. .AND. V .LE.  2580.)) IW=57
      IBAND=IW - 55
      IBND(9)=IW
      IF(IW .GT.  0) THEN
           QA(9)=ASO2(IBAND)
           AA(9)=AASO2(IBAND)
           BB(9)=BBSO2(IBAND)
           CC(9)=CCSO2(IBAND)
      ENDIF
      RETURN
      END
