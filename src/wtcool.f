      SUBROUTINE WTCOOL(BNDWID,IEMSCT)

!     ROUTINE TO WRITE THE IN-BAND COOLING RATES AND VERTICAL
!     FLUXES.  VALUES BETWEEN LAYER BOUNDARIES ARE DETERMINED
!     BY FITTING A CUBIC POLYNOMIAL IN LOG PRESSURE TO UPPER
!     AND LOWER LAYER BOUNDARY VALUES AND THEIR FIRST DERIVATIVES.
      IMPLICIT NONE

!     ARGUMENTS:
!       BNDWID   SPECTRAL BIN WIDTH [CM-1].
!       IEMSCT   MODTRAN RADIATION TRANSPORT MODE FLAG
!                  0 FOR TRANSMITTANCE ONLY
!                  1 FOR THERMAL RADIANCE
!                  2 FOR THERMAL + SOLAR RADIANCE
!                  3 FOR TRANSMITTED SOLAR IRRADIANCE
      REAL BNDWID
      INTEGER IEMSCT

!     PARAMETERS:
      INCLUDE 'PARAMS.h'
      INTEGER NGRID
      PARAMETER(NGRID=3)

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /COOLRT/
!       ISPCCR  UNIT NUMBER FOR SPECTRAL COOLING RATE FILE.
!       P3COEF  PRESSURE DEPENDENT COEFFICIENTS FOR CALCULATING COOLING
!               RATES WITH 3 POINT SPLINE FIT [(K/DAY) / (W/CM2)].
!       P5COEF  PRESSURE DEPENDENT COEFFICIENTS FOR CALCULATING COOLING
!               RATES WITH 5 POINT SPLINE FIT [(K/DAY) / (W/CM2)].
!       CLRTSM  LAYER BOUNDARY IN-BAND COOLING RATES [K/DAY].
!       CLRTS0  LAYER BOUNDARY IN-BAND COOLING RATES
!               EXCLUDING DIRECT SOLAR CONTRIBUTION [K/DAY].
!       NTFSM   LAYER BOUNDARY IN-BAND NET (THERMAL PLUS SCATTERED
!               SOLAR PLUS DIRECT SOLAR) UPWARD FLUX [W/CM2].
!       UPFSM   LAYER BOUNDARY IN-BAND UPWARD THERMAL (PLUS SCATTERED
!               SOLAR IF DISORT & NO AZIMUTH DEPENDENCE) FLUX [W/CM2].
!       DNFSM   LAYER BOUNDARY IN-BAND DOWNWARD THERMAL (PLUS SCATTERED
!               SOLAR IF DISORT & NO AZIMUTH DEPENDENCE) FLUX [W/CM2].
!       UPFSSM  LAYER BOUNDARY IN-BAND UPWARD SCATTERED SOLAR FLUX
!               (USED WITH DISORT ONLY IF AZIMUTH DEPENDENT) [W/CM2].
!       DNFSSM  LAYER BOUNDARY IN-BAND DOWNWARD SCATTERED SOLAR FLUX
!               (USED WITH DISORT ONLY IF AZIMUTH DEPENDENT) [W/CM2].
      DOUBLE PRECISION CLRTSM,CLRTS0,NTFSM,UPFSM,DNFSM,UPFSSM,DNFSSM
      REAL P3COEF,P5COEF
      INTEGER ISPCCR
      COMMON/COOLRT/CLRTSM(LAYDIM),CLRTS0(LAYDIM),NTFSM(LAYDIM),        &
     &  UPFSM(LAYDIM),DNFSM(LAYDIM),UPFSSM(LAYDIM),DNFSSM(LAYDIM),      &
     &  P3COEF(3,LAYDIM),P5COEF(5,LAYDIM),ISPCCR

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

!     /VRANGE/
!       VBNDMN   COMPUTATIONAL BANDPASS MINIMUM FREQUENCY [CM-1].
!       VBNDMX   COMPUTATIONAL BANDPASS MAXIMUM FREQUENCY [CM-1].
      REAL VBNDMN,VBNDMX
      COMMON/VRANGE/VBNDMN,VBNDMX

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     LOCAL VARIABLES:
      INTEGER I,L,LP1,LP2
      REAL CONVRT,PLOG1,PLOG2,PLOG3,PLOG21,PLOG31,PLOG32,DCOEF1,        &
     &  DCOEF2,DCOEF3,DZ1,DZ2,DCLRT1,DCLRT2,D0CRT1,D0CRT2,DUPF1,        &
     &  DUPF2,DDNF1,DDNF2,DUPFS1,DUPFS2,DDNFS1,DDNFS2,DNTF1,            &
     &  DNTF2,DENOM,RAT,RATM1,COEF1,COEF2,RAT2,P1LG21,P2LG21,           &
     &  TSTDIF,TSTMIN,TSTMAX,TST,COEF(2,NGRID),DCOEF(2,NGRID)
!9    INTEGER ICASE
!9    DATA ICASE/0/

!     THE IN-BAND FLUXES NEED TO BE MULTIPLIED BY THE
!     BAND WIDTH AND CONVERTED FROM W/CM2 TO W/M2.
      CONVRT=10000*BNDWID

!     WRITE HEADER AND GROUND DATA
      IF(IEMSCT.EQ.1)THEN
          WRITE(IPR,'(/A,2(F9.2,A),/3(/A),/4X,A1,F12.6,F13.6,4F14.6)')  &
     &      ' INTEGRATED COOLING RATES AND VERTICAL FLUXES FROM',       &
     &      VBNDMN,' TO',VBNDMX,' CM-1:',                               &
     &      ' CALC                                COOLING    '          &
     &      //'TOTAL FLUX          TOTAL THERMAL',                      &
     &      'LAYER    ALTITUDE     PRESSURE          RATE    '          &
     &      //' UP - DOWN       FLUX UP     FLUX DOWN',                 &
     &      '           (KM)           (MB)       (K/DAY)    '          &
     &      //'    (W/M2)        (W/M2)        (W/M2)','1',ZM(1),       &
     &      PM(1),-DBLE(BNDWID)*CLRTSM(1),DBLE(CONVRT)*NTFSM(1),        &
     &      DBLE(CONVRT)*UPFSM(1),DBLE(CONVRT)*DNFSM(1)
      ELSE
          WRITE(IPR,'(/A,2(F9.2,A),/3(/2A),/4X,A1,F12.6,F13.6,7F14.6)') &
     &      ' INTEGRATED COOLING RATES AND VERTICAL FLUXES FROM'        &
     &      ,VBNDMN,' TO',VBNDMX,' CM-1:',' CALC'//                     &
     &      '                                   COOLING RATE      ',    &
     &      '    TOTAL FLUX          TOTAL THERMAL     '//              &
     &      '         SCATTERED SOLAR','LAYER'//                        &
     &      '    ALTITUDE     PRESSURE        TOTAL  NO DIRECT SUN',    &
     &      '     UP-DN-DIR       FLUX UP     FLUX DOWN'//              &
     &      '       FLUX UP     FLUX DOWN','     '//                    &
     &      '      (KM)           (MB)       (K/DAY)       (K/DAY)',    &
     &      '        (W/M2)        (W/M2)        (W/M2)'//              &
     &      '        (W/M2)        (W/M2)',                             &
     &      '1',ZM(1),PM(1),-DBLE(BNDWID)*CLRTSM(1),                    &
     &      -DBLE(BNDWID)*CLRTS0(1),DBLE(CONVRT)*NTFSM(1),              &
     &      DBLE(CONVRT)*UPFSM(1),DBLE(CONVRT)*DNFSM(1),                &
     &      DBLE(CONVRT)*UPFSSM(1),DBLE(CONVRT)*DNFSSM(1)
      ENDIF

!     DETERMINE THE DERIVATIVES AT THE SECOND LAYER BOUNDARY.
      PLOG1=LOG(PM(1))
      PLOG2=LOG(PM(2))
      PLOG3=LOG(PM(3))
      PLOG21=PLOG2-PLOG1
      PLOG31=PLOG3-PLOG1
      PLOG32=PLOG3-PLOG2
      DCOEF1=-PLOG32/(PLOG21*PLOG31)
      DCOEF2=1/PLOG21-1/PLOG32
      DCOEF3=PLOG21/(PLOG31*PLOG32)
      DZ2=DCOEF1*REAL(ZM(1))+DCOEF2*REAL(ZM(2))+DCOEF3*REAL(ZM(3))
      DCLRT2=DCOEF1*REAL(CLRTSM(1))                                     &
     &      +DCOEF2*REAL(CLRTSM(2))+DCOEF3*REAL(CLRTSM(3))
      D0CRT2=DCOEF1*REAL(CLRTS0(1))                                     &
     &      +DCOEF2*REAL(CLRTS0(2))+DCOEF3*REAL(CLRTS0(3))
      DUPF2=DCOEF1*REAL(UPFSM(1))                                       &
     &     +DCOEF2*REAL(UPFSM(2))+DCOEF3*REAL(UPFSM(3))
      DDNF2=DCOEF1*REAL(DNFSM(1))                                       &
     &     +DCOEF2*REAL(DNFSM(2))+DCOEF3*REAL(DNFSM(3))
      DUPFS2=DCOEF1*REAL(UPFSSM(1))                                     &
     &      +DCOEF2*REAL(UPFSSM(2))+DCOEF3*REAL(UPFSSM(3))
      DDNFS2=DCOEF1*REAL(DNFSSM(1))                                     &
     &      +DCOEF2*REAL(DNFSSM(2))+DCOEF3*REAL(DNFSSM(3))
      DNTF2=DCOEF1*REAL(NTFSM(1))                                       &
     &     +DCOEF2*REAL(NTFSM(2))+DCOEF3*REAL(NTFSM(3))

!     DETERMINE VALUES WITHIN FIRST LAYER (BETWEEN PLOG1 AND PLOG2)
      DENOM=NGRID+1
      IF(CLRTSM(1).LT.CLRTSM(2))THEN
          TSTDIF=REAL(CLRTSM(2)-CLRTSM(1))/4
          TSTMIN=REAL(CLRTSM(1))-TSTDIF
          TSTMAX=REAL(CLRTSM(2))+TSTDIF
      ELSE
          TSTDIF=REAL(CLRTSM(1)-CLRTSM(2))/4
          TSTMIN=REAL(CLRTSM(2))-TSTDIF
          TSTMAX=REAL(CLRTSM(1))+TSTDIF
      ENDIF
      DO I=1,NGRID
          RAT=I/DENOM
          RATM1=RAT-1
          COEF1=RATM1**2
          COEF2=RAT*(2-RAT)
          DCOEF2=RAT*RATM1*PLOG21
          TST=COEF1*REAL(CLRTSM(1))+COEF2*REAL(CLRTSM(2))+DCOEF2*DCLRT2
          IF(TST.GT.TSTMIN .AND. TST.LT.TSTMAX)THEN
              IF(IEMSCT.EQ.1)THEN
                  WRITE(IPR,'(F17.6,F13.6,4F14.6)')                     &
     &              (COEF1*REAL(ZM(1))+COEF2*REAL(ZM(2))+DCOEF2*DZ2),   &
     &              EXP(PLOG1+RAT*PLOG21),-BNDWID*TST,                  &
     &               CONVRT*(COEF1*REAL(NTFSM(1))                       &
     &                      +COEF2*REAL(NTFSM(2))+DCOEF2*DNTF2 ),       &
     &               CONVRT*(COEF1*REAL(UPFSM(1))                       &
     &                      +COEF2*REAL(UPFSM(2))+DCOEF2*DUPF2 ),       &
     &               CONVRT*(COEF1*REAL(DNFSM(1))                       &
     &                      +COEF2*REAL(DNFSM(2))+DCOEF2*DDNF2 )
              ELSE
                  WRITE(IPR,'(F17.6,F13.6,7F14.6)')                     &
     &              (COEF1*REAL(ZM(1))+COEF2*REAL(ZM(2))+DCOEF2*DZ2),   &
     &              EXP(PLOG1+RAT*PLOG21),-BNDWID*TST,                  &
     &              -BNDWID*(COEF1*REAL(CLRTS0(1))                      &
     &                      +COEF2*REAL(CLRTS0(2))+DCOEF2*D0CRT2),      &
     &               CONVRT*(COEF1*REAL(NTFSM(1) )                      &
     &                      +COEF2*REAL(NTFSM(2) )+DCOEF2*DNTF2 ),      &
     &               CONVRT*(COEF1*REAL(UPFSM(1) )                      &
     &                      +COEF2*REAL(UPFSM(2) )+DCOEF2*DUPF2 ),      &
     &               CONVRT*(COEF1*REAL(DNFSM(1) )                      &
     &                      +COEF2*REAL(DNFSM(2) )+DCOEF2*DDNF2 ),      &
     &               CONVRT*(COEF1*REAL(UPFSSM(1))                      &
     &                      +COEF2*REAL(UPFSSM(2))+DCOEF2*DUPFS2),      &
     &               CONVRT*(COEF1*REAL(DNFSSM(1))                      &
     &                      +COEF2*REAL(DNFSSM(2))+DCOEF2*DDNFS2)
               ENDIF
           ENDIF

!         DETERMINE COEFFICIENTS FOR INTERMEDIATE LAYERS
          RAT2=RAT**2
          COEF(1,I)=COEF1*(1+2*RAT)
          COEF(2,I)=RAT2*(3-2*RAT)
          DCOEF(1,I)=COEF1*RAT
          DCOEF(2,I)=RAT2*RATM1
      ENDDO

!     WRITE VALUES FOR TOP OF THE FIRST LAYER.
      IF(IEMSCT.EQ.1)THEN
          WRITE(IPR,'(4X,A1,F12.6,F13.6,4F14.6)')'2',ZM(2),PM(2),       &
     &      -DBLE(BNDWID)*CLRTSM(2),DBLE(CONVRT)*NTFSM(2),              &
     &      DBLE(CONVRT)*UPFSM(2),DBLE(CONVRT)*DNFSM(2)
      ELSE
          WRITE(IPR,'(4X,A1,F12.6,F13.6,7F14.6)')                       &
     &      '2',ZM(2),PM(2),-DBLE(BNDWID)*CLRTSM(2),                    &
     &      -DBLE(BNDWID)*CLRTS0(2),DBLE(CONVRT)*NTFSM(2),              &
     &      DBLE(CONVRT)*UPFSM(2),DBLE(CONVRT)*DNFSM(2),                &
     &      DBLE(CONVRT)*UPFSSM(2),DBLE(CONVRT)*DNFSSM(2)
      ENDIF

!     LOOP OVER INTERMEDIATE LAYERS.
      L=2
      LP1=3
      DO LP2=4,ML

!         DETERMINE UPPER BOUNDARY DERIVATIVES.
          DZ1=DZ2
          DCLRT1=DCLRT2
          D0CRT1=D0CRT2
          DUPF1=DUPF2
          DDNF1=DDNF2
          DUPFS1=DUPFS2
          DDNFS1=DDNFS2
          DNTF1=DNTF2
          PLOG1=PLOG2
          PLOG2=PLOG3
          PLOG3=LOG(PM(LP2))
          PLOG21=PLOG2-PLOG1
          PLOG31=PLOG3-PLOG1
          PLOG32=PLOG3-PLOG2
          DCOEF1=-PLOG32/(PLOG21*PLOG31)
          DCOEF2=1/PLOG21-1/PLOG32
          DCOEF3=PLOG21/(PLOG31*PLOG32)
          DZ2=DCOEF1*REAL(ZM(L))+DCOEF2*REAL(ZM(LP1))                   &
     &                          +DCOEF3*REAL(ZM(LP2))
          DCLRT2=DCOEF1*REAL(CLRTSM(L))                                 &
     &          +DCOEF2*REAL(CLRTSM(LP1))+DCOEF3*REAL(CLRTSM(LP2))
          D0CRT2=DCOEF1*REAL(CLRTS0(L))                                 &
     &          +DCOEF2*REAL(CLRTS0(LP1))+DCOEF3*REAL(CLRTS0(LP2))
          DUPF2=DCOEF1*REAL(UPFSM(L))                                   &
     &         +DCOEF2*REAL(UPFSM(LP1))+DCOEF3*REAL(UPFSM(LP2))
          DDNF2=DCOEF1*REAL(DNFSM(L))                                   &
     &         +DCOEF2*REAL(DNFSM(LP1))+DCOEF3*REAL(DNFSM(LP2))
          DUPFS2=DCOEF1*REAL(UPFSSM(L))                                 &
     &          +DCOEF2*REAL(UPFSSM(LP1))+DCOEF3*REAL(UPFSSM(LP2))
          DDNFS2=DCOEF1*REAL(DNFSSM(L))                                 &
     &          +DCOEF2*REAL(DNFSSM(LP1))+DCOEF3*REAL(DNFSSM(LP2))
          DNTF2=DCOEF1*REAL(NTFSM(L))                                   &
     &         +DCOEF2*REAL(NTFSM(LP1))+DCOEF3*REAL(NTFSM(LP2))
          IF(CLRTSM(L).LT.CLRTSM(LP1))THEN
              TSTDIF=REAL(CLRTSM(LP1)-CLRTSM(L  ))/4
              TSTMIN=REAL(CLRTSM(L  ))-TSTDIF
              TSTMAX=REAL(CLRTSM(LP1))+TSTDIF
          ELSE
              TSTDIF=REAL(CLRTSM(L  )-CLRTSM(LP1))/4
              TSTMIN=REAL(CLRTSM(LP1))-TSTDIF
              TSTMAX=REAL(CLRTSM(L  ))+TSTDIF
          ENDIF

!         WRITE VALUES WITHIN LAYER.
          DO I=1,NGRID
              P1LG21=PLOG21*DCOEF(1,I)
              P2LG21=PLOG21*DCOEF(2,I)
              TST=COEF(1,I)*REAL(CLRTSM(L  ))+P1LG21*DCLRT1+            &
     &            COEF(2,I)*REAL(CLRTSM(LP1))+P2LG21*DCLRT2
              IF(TST.GT.TSTMIN .AND. TST.LT.TSTMAX)THEN
                  IF(IEMSCT.EQ.1)THEN
                      WRITE(IPR,'((F17.6,F13.6,4F14.6))')               &
     &                  (COEF(1,I)*REAL(ZM(L  ))+P1LG21*DZ1             &
     &                  +COEF(2,I)*REAL(ZM(LP1))+P2LG21*DZ2),           &
     &                  EXP(PLOG1+I*PLOG21/DENOM),-BNDWID*TST,          &
     &                  CONVRT*(COEF(1,I)*REAL(NTFSM(L))+P2LG21*DNTF2   &
     &                    +COEF(2,I)*REAL(NTFSM(LP1))+P1LG21*DNTF1),    &
     &                  CONVRT*(COEF(1,I)*REAL(UPFSM(L))+P1LG21*DUPF1   &
     &                    +COEF(2,I)*REAL(UPFSM(LP1))+P2LG21*DUPF2),    &
     &                  CONVRT*(COEF(1,I)*REAL(DNFSM(L))+P1LG21*DDNF1   &
     &                    +COEF(2,I)*REAL(DNFSM(LP1))+P2LG21*DDNF2)
                  ELSE
                      WRITE(IPR,'((F17.6,F13.6,7F14.6))')               &
     &                  (COEF(1,I)*REAL(ZM(L  ))+P1LG21*DZ1             &
     &                  +COEF(2,I)*REAL(ZM(LP1))+P2LG21*DZ2),           &
     &                  EXP(PLOG1+I*PLOG21/DENOM),-BNDWID*TST,          &
     &                  -BNDWID*(COEF(1,I)*REAL(CLRTS0(L  ))            &
     &                          +COEF(2,I)*REAL(CLRTS0(LP1))            &
     &                          +P1LG21*D0CRT1+P2LG21*D0CRT2),          &
     &                  CONVRT*(COEF(1,I)*REAL(NTFSM(L))+P1LG21*DNTF1   &
     &                    +COEF(2,I)*REAL(NTFSM(LP1))+P2LG21*DNTF2),    &
     &                  CONVRT*(COEF(1,I)*REAL(UPFSM(L))+P1LG21*DUPF1   &
     &                    +COEF(2,I)*REAL(UPFSM(LP1))+P2LG21*DUPF2),    &
     &                  CONVRT*(COEF(1,I)*REAL(DNFSM(L))+P1LG21*DDNF1   &
     &                    +COEF(2,I)*REAL(DNFSM(LP1))+P2LG21*DDNF2),    &
     &                  CONVRT*(COEF(1,I)*REAL(UPFSSM(L))+P1LG21*DUPFS1 &
     &                    +COEF(2,I)*REAL(UPFSSM(LP1))+P2LG21*DUPFS2),  &
     &                  CONVRT*(COEF(1,I)*REAL(DNFSSM(L))+P1LG21*DDNFS1 &
     &                    +COEF(2,I)*REAL(DNFSSM(LP1))+P2LG21*DDNFS2)
                  ENDIF
              ENDIF
          ENDDO

!         WRITE VALUES FOR TOP OF THE CURRENT LAYER.
          L=LP1
          LP1=LP2
          IF(IEMSCT.EQ.1)THEN
              WRITE(IPR,'(I5,F12.6,F13.6,4F14.6)')L,ZM(L),PM(L),        &
     &          -DBLE(BNDWID)*CLRTSM(L),DBLE(CONVRT)*NTFSM(L),          &
     &          DBLE(CONVRT)*UPFSM(L),DBLE(CONVRT)*DNFSM(L)
          ELSE
              WRITE(IPR,'(I5,F12.6,F13.6,7F14.6)')L,ZM(L),PM(L),        &
     &          -DBLE(BNDWID)*CLRTSM(L),-DBLE(BNDWID)*CLRTS0(L),        &
     &          DBLE(CONVRT)*NTFSM(L),DBLE(CONVRT)*UPFSM(L),            &
     &          DBLE(CONVRT)*DNFSM(L),DBLE(CONVRT)*UPFSSM(L),           &
     &          DBLE(CONVRT)*DNFSSM(L)
          ENDIF
      ENDDO

!     DETERMINE VALUES WITHIN LAST LAYER.
      IF(CLRTSM(L).LT.CLRTSM(ML))THEN
          TSTDIF=REAL(CLRTSM(ML)-CLRTSM(L ))/4
          TSTMIN=REAL(CLRTSM(L ))-TSTDIF
          TSTMAX=REAL(CLRTSM(ML))+TSTDIF
      ELSE
          TSTDIF=REAL(CLRTSM(L )-CLRTSM(ML))/4
          TSTMIN=REAL(CLRTSM(ML))-TSTDIF
          TSTMAX=REAL(CLRTSM(L ))+TSTDIF
      ENDIF
      DO I=1,NGRID
          RAT=I/DENOM
          COEF1=(1+RAT)*(1-RAT)
          DCOEF1=RAT*(1-RAT)*PLOG32
          COEF2=RAT**2
          TST=COEF1*REAL(CLRTSM(L))+DCOEF1*DCLRT2+COEF2*REAL(CLRTSM(ML))
          IF(TST.GT.TSTMIN .AND. TST.LT.TSTMAX)THEN
              IF(IEMSCT.EQ.1)THEN
                  WRITE(IPR,'(F17.6,F13.6,4F14.6)')                     &
     &              (COEF1*REAL(ZM(L))+DCOEF1*DZ2+COEF2*REAL(ZM(ML))),  &
     &              EXP(PLOG2+RAT*PLOG32),-BNDWID*TST,                  &
     &              CONVRT*(COEF1*REAL(NTFSM( L)) +DCOEF1*DNTF2         &
     &                     +COEF2*REAL(NTFSM(ML)) ),                    &
     &              CONVRT*(COEF1*REAL(UPFSM( L)) +DCOEF1*DUPF2         &
     &                     +COEF2*REAL(UPFSM(ML)) ),                    &
     &              CONVRT*(COEF1*REAL(DNFSM( L)) +DCOEF1*DDNF2         &
     &                     +COEF2*REAL(DNFSM(ML)) )
              ELSE
                  WRITE(IPR,'(F17.6,F13.6,7F14.6)')                     &
     &              (COEF1*REAL(ZM(L))+DCOEF1*DZ2+COEF2*REAL(ZM(ML))),  &
     &              EXP(PLOG2+RAT*PLOG32),-BNDWID*TST,                  &
     &              -BNDWID*(COEF1*REAL(CLRTS0( L))+DCOEF1*D0CRT2       &
     &                      +COEF2*REAL(CLRTS0(ML))),                   &
     &              CONVRT*(COEF1*REAL(NTFSM( L)) +DCOEF1*DNTF2         &
     &                     +COEF2*REAL(NTFSM(ML)) ),                    &
     &              CONVRT*(COEF1*REAL(UPFSM( L)) +DCOEF1*DUPF2         &
     &                     +COEF2*REAL(UPFSM(ML)) ),                    &
     &              CONVRT*(COEF1*REAL(DNFSM( L)) +DCOEF1*DDNF2         &
     &                     +COEF2*REAL(DNFSM(ML)) ),                    &
     &              CONVRT*(COEF1*REAL(UPFSSM( L))+DCOEF1*DUPFS2        &
     &                     +COEF2*REAL(UPFSSM(ML))),                    &
     &              CONVRT*(COEF1*REAL(DNFSSM( L))+DCOEF1*DDNFS2        &
     &                     +COEF2*REAL(DNFSSM(ML)))
              ENDIF
          ENDIF
      ENDDO

!     WRITE VALUES FOR TOP OF ATMOSPHERE.
      IF(IEMSCT.EQ.1)THEN
          WRITE(IPR,'(I5,F12.6,F13.6,4F14.6)')ML,ZM(ML),PM(ML),         &
     &      -DBLE(BNDWID)*CLRTSM(ML),DBLE(CONVRT)*NTFSM(ML),            &
     &      DBLE(CONVRT)*UPFSM(ML),DBLE(CONVRT)*DNFSM(ML)
      ELSE
          WRITE(IPR,'(I5,F12.6,F13.6,7F14.6)')ML,ZM(ML),PM(ML),         &
     &      -DBLE(BNDWID)*CLRTSM(ML),-DBLE(BNDWID)*CLRTS0(ML),          &
     &      DBLE(CONVRT)*NTFSM(ML),DBLE(CONVRT)*UPFSM(ML),              &
     &      DBLE(CONVRT)*DNFSM(ML),DBLE(CONVRT)*UPFSSM(ML),             &
     &      DBLE(CONVRT)*DNFSSM(ML)
      ENDIF
      RETURN
      END
