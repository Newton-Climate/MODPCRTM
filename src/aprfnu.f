      SUBROUTINE APRFNU(LMODEL,IHAZE,GNDALT)

!     THIS ROUTINE RELAYERS THE AEROSOLS
!     A PAIR OF ALTITUDES (TOP AND BASE) PER AEROSOL IS INPUT
!     E.G., FOR THE 1ST AEROSOL, ZAER11 (ZAER12) IS THE NEW BASE (TOP)
!     THE AEROSOL'S ORIGINAL BASE IS MOVED TO THE NEW BASE

!     THE TOP CAN BE LEFT BLANK OR ZERO, IN WHICH CASE THE AEROSOL IS
!     SIMPLY TRANSLATED TO THE NEW BASE
!     IF TOP .LE. BOTTOM, THE TRANSLATION TAKES PLACE

!     ZM IS THE MODTRAN ALTITUDES.
!     IF THE BASE AND THE TOP DO NOT COINCIDE A ZM ALTITUDE,
!     THEY ARE INSERTED INTO ZM.
!     OTHER MODTRAN SPECIES ARE GIVEN VALUES AT THE NEW ALTITUDES
      IMPLICIT NONE

!     PARAMETERS:
!       ZMICRO   LAYER WIDTH [KM]
      DOUBLE PRECISION ZMICRO
      PARAMETER(ZMICRO=.001D0)
      INCLUDE 'PARAMS.h'

!     ARGUMENTS:
      LOGICAL LMODEL
      INTEGER IHAZE
      DOUBLE PRECISION GNDALT

!     COMMONS:

!     /JM2APLUS/
      DOUBLE PRECISION ZAER11,ZAER12,ZAER21,ZAER22,                     &
     &  ZAER31,ZAER32,ZAER41,ZAER42
      REAL SCALE1,SCALE2,SCALE3,SCALE4
      COMMON/JM2APLUS/ZAER11,ZAER12,ZAER21,ZAER22,ZAER31,ZAER32,        &
     &  ZAER41,ZAER42,SCALE1,SCALE2,SCALE3,SCALE4

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

!     DECLARE LOCAL VARIABLES AND ARGUMENTS
      INTEGER INDX1,INDX2,NAER1,NAER2,NAER3,NAER4,NMRG
      INTEGER I,J,JLO,JHI,II,IIMIN1
      LOGICAL LAER1,LAER2,LAER3,LAER4,LALLZERO
      REAL FAC,DEN,AERNU1(LAYDIM),AERNU2(LAYDIM),AERNU3(LAYDIM),        &
     &  AERNU4(LAYDIM)
      DOUBLE PRECISION FACZ,ZDFLT1,ZDFLT2,DALT,ZMNU1(LAYDIM),           &
     &  ZMNU2(LAYDIM),ZMNU3(LAYDIM),ZMNU4(LAYDIM),ALTMRG(LAYDIM)

!     DECLARE FUNCTION NAMES
      REAL EXPINT

!     IF (IHAZE.EQ.7) AEROSOL ENHANCEMENTS DO NOT WORK
!     IF (IHAZE.EQ.7) USER-DEFINED AEROSOL PROFILE OPTION IS ACTIVATED
      NAER1=0
      NAER2=0
      NAER3=0
      NAER4=0
      DO I=1,LAYDIM
          ZMNU1(I)=0.D0
          AERNU1(I)=0.
          ZMNU2(I)=0.D0
          AERNU2(I)=0.
          ZMNU3(I)=0.D0
          AERNU3(I)=0.
          ZMNU4(I)=0.D0
          AERNU4(I)=0.
      ENDDO
      IF(SCALE1.EQ.0.)SCALE1=1.
      IF(SCALE2.EQ.0.)SCALE2=1.
      IF(SCALE3.EQ.0.)SCALE3=1.
      IF(SCALE4.EQ.0.)SCALE4=1.

      IF(IHAZE.LE.0 .OR. IHAZE.EQ.7)RETURN

      IF(ZAER11.NE.ZAER12 .OR. SCALE1.NE.1.)THEN

!        AEROSOL NUMBER 1 (BOUNDARY LAYER AEROSOL (0-3 KM))
         ZDFLT1=ZM(1)
         IF(IHAZE.EQ.6)THEN
!           IHAZE.EQ.6 COMBINES AEROSOLS 1 AND 2 INTO A SINGLE REGION.
            ZDFLT2=11.D0
         ELSE
            ZDFLT2=3.D0
         ENDIF

!        LIFTED FROM FLAYZ (FOR NON-ZERO GROUND)
         IF(GNDALT.GT.0.D0 .AND. IHAZE.NE.6)THEN
!           ZDFLT2 WILL NOT CHANGE FOR IHAZE=6 BECAUSE IT IS ABOVE 6 KM.
            DALT=(6-GNDALT)/6
            ZDFLT2=GNDALT+DALT*ZDFLT2
         ENDIF

         IF(.NOT.LMODEL)THEN
            LALLZERO=.TRUE.
!           find zdflt2
            DO I=ML,1,-1
               IF(DENSTY(7,I).GT.0.)THEN
                  LALLZERO=.FALSE.
                  ZDFLT2=ZM(I+1)
                  GOTO 10
               ENDIF
            ENDDO
   10       CONTINUE
            IF(LALLZERO)GOTO 20
         ENDIF

         ZAER11=MAX(ZAER11,ZDFLT1)
         IF(ZAER12.LE.ZAER11)ZAER12=ZAER11+(ZDFLT2-ZDFLT1)

         ZAER12=MIN(ZAER12,ZM(ML))
         IF(ZAER12.LE.ZAER11)THEN
            IF(LJMASS)CALL WRTBUF(FATAL)
            STOP'APRFNU: TOP OF AER 1 .LE. BASE'
         ENDIF
         ZMNU1(1)=ZAER11-ZMICRO
         AERNU1(1)=0
         CALL FINDEX(ZM,ML,ZDFLT2,INDX2,JLO)
         IF (INDX2.EQ.-1)THEN
            IF(LJMASS)CALL WRTBUF(FATAL)
            STOP'APRFNU: PROBLEM W AEROSOL 1'
         ENDIF
         NAER1=INDX2+1
         FACZ=(ZAER11-ZAER12)/(ZM(1)-ZM(INDX2))
         J=1
         DO I=1,INDX2
            J=J+1
            ZMNU1(J)=FACZ*(ZM(I)-ZM(INDX2))+ZAER12
            AERNU1(J)=(DENSTY(7,I)/SNGL(FACZ))*SCALE1
         ENDDO
         IF(ZMNU1(1).LT.ZDFLT1)THEN

!           DROP THE BOTTOM "MICRO LAYER"
            NAER1=NAER1-1
            DO I=1,NAER1
               ZMNU1(I)=ZMNU1(I+1)
               AERNU1(I)=AERNU1(I+1)
            ENDDO
         ENDIF

!        ZMNU1 IS A MAPPING OF THE OLD ALTITUDES ZM WHERE AEROSOL WAS
!        DEFINED.  THERE WILL PROBABLY BE SOME ZM VALUES INTERSPERSED
!        BETWEEN ZMNU1 WHERE AERNU1 IS NOT DEFINED,
!        FOR EXAMPLE, SUPPOSE ZMNU1 IS AT 0, .5, 1.5, 2.5 AND 3.5
!        THIS ARRAY SHOULD EXTEND TO HAVE 1, 2 AND 3 KM, WHICH ARE IN ZM
!        SIMILARLY, SOME NEW ALTITUDES COULD EMERGE FROM OTHER AEROSOLS.
!        SO WE NEED TO KEEP A LIST OF MERGED ALTITUDES FOR LATER USE.

         CALL ZMRG(ZM,ML,ZMNU1,NAER1,ALTMRG,NMRG)
      ELSE
         NMRG=ML
         DO I=1,ML
            ALTMRG(I)=ZM(I)
         ENDDO
      ENDIF

   20 CONTINUE
      IF((ZAER21.NE.ZAER22 .OR. SCALE2.NE.1.) .AND. IHAZE.NE.6)THEN
!        IHAZE.EQ.6 DOES NOT NEED THIS SECTION;
!        THIS REGION IS MERGED INTO REGION 1.

!        AEROSOL NUMBER 2 (TROPOSPHERE, 2-11 KM)
         ZDFLT1=2.D0
         ZDFLT2=11.D0

!        LIFTED FROM FLAYZ (FOR NON-ZERO GROUND)
         IF(GNDALT.GT.0.D0)THEN
            DALT=(6-GNDALT)/6
            ZDFLT1=GNDALT+DALT*ZDFLT1
         ENDIF

         IF(.NOT.LMODEL)THEN
            LALLZERO=.TRUE.
!           find zdflt1
            DO I=1,ML
               IF(DENSTY(12,I).GT.0.)THEN
                  LALLZERO=.FALSE.
                  ZDFLT1=ZM(I-1)
                  GOTO 30
               ENDIF
            ENDDO
   30       CONTINUE
!           find zdflt2
            DO I=ML,1,-1
               IF(DENSTY(12,I).GT.0.)THEN
                  ZDFLT2=ZM(I+1)
                  GOTO 40
               ENDIF
            ENDDO
   40       CONTINUE

            IF(LALLZERO)GOTO 50
         ENDIF

         ZAER21=MAX(ZAER21,ZM(1))
         IF(ZAER22.LE.ZAER21)ZAER22=ZAER21+(ZDFLT2-ZDFLT1)

         ZAER22=MIN(ZAER22,ZM(ML))
         IF(ZAER22.LE.ZAER21)THEN
            IF(LJMASS)CALL WRTBUF(FATAL)
            STOP'APRFNU: TOP OF AER 2 .LE. BASE'
         ENDIF
         CALL FINDEX(ZM,ML,ZDFLT1,INDX1,JLO)
         CALL FINDEX(ZM,ML,ZDFLT2,INDX2,JLO)
         IF (INDX1.EQ.-1)THEN
            IF(LJMASS)CALL WRTBUF(FATAL)
            STOP'APRFNU: PROBLEM W INDX1 W AEROSOL2'
         ENDIF
         IF(INDX2.EQ.-1)THEN
            IF(LJMASS)CALL WRTBUF(FATAL)
            STOP'APRFNU: PROBLEM W INDX2 W AEROSOL2'
         ENDIF
         NAER2=INDX2-INDX1+1
         FACZ=(ZAER21-ZAER22)/(ZM(INDX1)-ZM(INDX2))
         J=0
         DO I=INDX1,INDX2
            J=J+1
            ZMNU2(J)=FACZ*(ZM(I)-ZM(INDX2))+ZAER22
            AERNU2(J)=(DENSTY(12,I)/SNGL(FACZ))*SCALE2
         ENDDO

!        ZMNU2 IS A MAPPING OF THE OLD ALTS ZM.
!        THERE WILL PROBABLY BE SOME ZM VALUES INTERSPERSED BETWEEN
!        ZMNU2 WHERE AERNU2 IS NOR DEFINED,
!        FOR EXAMPLE, SUPPOSE ZMNU2 IS AT 0, .5, 1.5, 2.5 AND 3.5
!        THIS ARRAY SHOULD EXTENDED TO INCLUDE 1, 2 AND 3 KM
!        SIMILARLY, SOME NEW ALTITUDES COULD EMERGE FROM OTHER AEROSOLS.
!        SO WE NEED TO KEEP A LIST OF MERGED ALTITUDES FOR LATER USE.

         CALL ZMRG2(ALTMRG,NMRG,ZMNU2,NAER2)
      ENDIF

   50 CONTINUE
      IF(ZAER31.NE.ZAER32 .OR. SCALE3.NE.1.)THEN

!        AEROSOL NUMBER 3 (STRATOSPHERE, 10-35 KM)
         ZDFLT1=10.D0
         ZDFLT2=35.D0

         IF(.NOT.LMODEL)THEN
            LALLZERO=.TRUE.
!           find zdflt1
            DO I=1,ML
               IF(IHAZE.NE.6)THEN
                  DEN=DENSTY(13,I)
               ELSE
                  DEN=DENSTY(12,I)
               ENDIF
               IF(DEN.GT.0.)THEN
                  LALLZERO=.FALSE.
                  ZDFLT1=ZM(I-1)
                  GOTO 60
               ENDIF
            ENDDO
   60       CONTINUE

!           find zdflt2
            DO I=ML,1,-1
               IF(IHAZE.NE.6)THEN
                  DEN=DENSTY(13,I)
               ELSE
                  DEN=DENSTY(12,I)
               ENDIF
               IF(DEN.GT.0.)THEN
                  ZDFLT2=ZM(I+1)
                  GOTO 70
               ENDIF
            ENDDO
   70       CONTINUE

            IF(LALLZERO)GOTO 80
         ENDIF

         ZAER31=MAX(ZAER31,ZM(1))
         IF(ZAER32.LE.ZAER31)ZAER32=ZAER31+(ZDFLT2-ZDFLT1)

         ZAER32=MIN(ZAER32,ZM(ML))
         IF(ZAER32.LE.ZAER31)STOP'APRFNU: TOP OF AER 3 .LE. BASE'
         CALL FINDEX(ZM,ML,ZDFLT1,INDX1,JLO)
         CALL FINDEX(ZM,ML,ZDFLT2,INDX2,JLO)
         IF(INDX1.EQ.-1)THEN
            IF(LJMASS)CALL WRTBUF(FATAL)
            STOP 'APRFNU:  PROBLEM W INDX1 W AEROSOL3'
         ENDIF
         IF(INDX2.EQ.-1)THEN
            IF(LJMASS)CALL WRTBUF(FATAL)
            STOP 'APRFNU:  PROBLEM W INDX2 W AEROSOL3'
         ENDIF
         NAER3=INDX2-INDX1+1
         FACZ=(ZAER31-ZAER32)/(ZM(INDX1)-ZM(INDX2))
         J=0
         DO I=INDX1,INDX2
            J=J+1
            ZMNU3(J)=FACZ*(ZM(I)-ZM(INDX2))+ZAER32
            IF(IHAZE.NE.6)THEN
               AERNU3(J)=(DENSTY(13,I)/SNGL(FACZ))*SCALE3
            ELSE
!              IHAZE.EQ.6 PUTS AEROSOL 3 IN THE SECOND REGION.
               AERNU3(J)=(DENSTY(12,I)/SNGL(FACZ))*SCALE3
            ENDIF
         ENDDO

!        ZMNU3 IS A MAPPING OF THE OLD ALTS ZM.
!        THERE WILL PROBABLY BE SOME ZM VALUES INTERSPERSED BETWEEN
!        ZMNU3 WHERE AERNU3 IS NOR DEFINED,
!        FOR EXAMPLE, SUPPOSE ZMNU3 IS AT 0, .5, 1.5, 2.5 AND 3.5
!        THIS ARRAY SHOULD EXTENDED TO INCLUDE 1, 2 AND 3 KM
!        SIMILARLY, SOME NEW ALTITUDES COULD EMERGE FROM OTHER AEROSOLS.
!        SO WE NEED TO KEEP A LIST OF MERGED ALTITUDES FOR LATER USE.

         CALL ZMRG2(ALTMRG,NMRG,ZMNU3,NAER3)
      ENDIF

   80 CONTINUE
      IF(ZAER42.NE.ZAER41 .OR. SCALE4.NE.1.)THEN

!        AEROSOL NUMBER 4 (STRATOSPHERE AND BEYOND, 30-100 KM)
         ZDFLT1=30.D0
         ZAER41=MAX(ZAER41,ZM(1))
         ZDFLT2=MIN(100.D0,ZM(ML))

         IF(.NOT.LMODEL)THEN
            LALLZERO=.TRUE.
!           find zdflt1
            DO I=1,ML
               IF(IHAZE.NE.6)THEN
                  DEN=DENSTY(14,I)
               ELSE
                  DEN=DENSTY(13,I)
               ENDIF
               IF(DEN.GT.0.)THEN
                  LALLZERO=.FALSE.
                  ZDFLT1=ZM(I-1)
                  GOTO 90
               ENDIF
            ENDDO
   90       CONTINUE

            IF(LALLZERO)GOTO 100
         ENDIF

         IF(ZAER42.LE.ZAER41)ZAER42=ZAER41+(ZDFLT2-ZDFLT1)

         ZAER42=MIN(ZAER42,ZM(ML))
         IF(ZAER42.LE.ZAER41)THEN
            IF(LJMASS)CALL WRTBUF(FATAL)
            STOP 'APRFNU: TOP OF AER 4 .LE. BASE'
         ENDIF
         CALL FINDEX(ZM,ML,ZDFLT1,INDX1,JLO)
         CALL FINDEX(ZM,ML,ZDFLT2,INDX2,JLO)
         IF(INDX1.EQ.-1)THEN
            IF(LJMASS)CALL WRTBUF(FATAL)
            STOP'APRFNU: PROBLEM W INDX1 W AEROSOL 4'
         ENDIF
         IF(INDX2.EQ.-1)THEN
            IF(LJMASS)CALL WRTBUF(FATAL)
            STOP'APRFNU: PROBLEM W INDX2 W AEROSOL 4'
         ENDIF
         NAER4=INDX2-INDX1+1
         FACZ=(ZAER41-ZAER42)/(ZM(INDX1)-ZM(INDX2))
         J=0
         DO I=INDX1,INDX2
            J=J+1
            ZMNU4(J)=FACZ*(ZM(I)-ZM(INDX2))+ZAER42
            IF(IHAZE.NE.6)THEN
               AERNU4(J)=(DENSTY(14,I)/SNGL(FACZ))*SCALE4
            ELSE
!              IHAZE.EQ.6 PUTS AEROSOL 4 IN THE THIRD REGION.
               AERNU4(J)=(DENSTY(13,I)/SNGL(FACZ))*SCALE4
            ENDIF
         ENDDO
         IF(ZAER42+ZMICRO.LE.ZDFLT2)THEN

!           ADD A MICRO LAYER
            NAER4=NAER4+1
            ZMNU4(NAER4)=ZAER42+ZMICRO
            AERNU4(NAER4)=0.
         ENDIF

!        ZMNU4 IS A MAPPING OF THE OLD ALTS ZM.
!        THERE WILL PROBABLY BE SOME ZM VALUES INTERSPERSED BETWEEN
!        ZMNU4 WHERE AERNU4 IS NOR DEFINED,
!        FOR EXAMPLE, SUPPOSE ZMNU4 IS AT 0, .5, 1.5, 2.5 AND 3.5
!        THIS ARRAY SHOULD EXTENDED TO INCLUDE 1, 2 AND 3 KM.
!        SIMILARLY, SOME NEW ALTITUDES COULD EMERGE FROM OTHER AEROSOLS.
!        SO WE NEED TO KEEP A LIST OF MERGED ALTITUDES FOR LATER USE.

         CALL ZMRG2(ALTMRG,NMRG,ZMNU4,NAER4)
      ENDIF

  100 CONTINUE
      IF(NAER1.GT.0)THEN

          DO I=1,NMRG
              IF(ALTMRG(I).GE.ZMNU1(NAER1))GOTO 110
              IF(ALTMRG(I).GT.ZMNU1(1))THEN
                  CALL FINDEX(ZMNU1,NAER1,ALTMRG(I),INDX1,JLO)
                  IF(INDX1.EQ.-1)THEN

!                     NO MATCH W/ ZMNU; SO INSERT.
                      NAER1=NAER1+1
                      J=JLO+1
                      JHI=JLO+2
                      DO II=NAER1,JHI,-1
                          IIMIN1=II-1
                          ZMNU1(II)=ZMNU1(IIMIN1)
                          AERNU1(II)=AERNU1(IIMIN1)
                      ENDDO
                      ZMNU1(J)=ALTMRG(I)
                      FAC=SNGL((ZMNU1(J)-ZMNU1(JHI))                    &
     &                        /(ZMNU1(JLO)-ZMNU1(JHI)))
                      AERNU1(J)=EXPINT(AERNU1(JHI),AERNU1(JLO),FAC)
                  ENDIF
              ENDIF
          ENDDO
 110      CONTINUE

          DO I=1,LAYDIM
              DENSTY(7,I)=0.
           ENDDO
          LAER1=.TRUE.
          LAER2=.FALSE.
          LAER3=.FALSE.
          LAER4=.FALSE.
          CALL AERMRG(NAER1,ZMNU1,AERNU1,ML,LAER1,LAER2,LAER3,LAER4)
      ENDIF
      IF(NAER2.GT.0.AND.IHAZE.NE.6)THEN
!         IHAZE.EQ.6 MEANS SKIP THIS REGION.

          DO I=1,NMRG
              IF(ALTMRG(I).GE.ZMNU2(NAER2))GOTO 120
              IF(ALTMRG(I).GT.ZMNU2(1))THEN
                  CALL FINDEX(ZMNU2,NAER2,ALTMRG(I),INDX1,JLO)
                  IF(INDX1.EQ.-1)THEN

!                     NO MATCH W/ ZMNU; SO INSERT.
                      NAER2=NAER2+1
                      J=JLO+1
                      JHI=JLO+2
                      DO II=NAER2,JHI,-1
                          IIMIN1=II-1
                          ZMNU2(II)=ZMNU2(IIMIN1)
                          AERNU2(II)=AERNU2(IIMIN1)
                      ENDDO
                      ZMNU2(J)=ALTMRG(I)
                      FAC=SNGL((ZMNU2(J)-ZMNU2(JHI))                    &
     &                        /(ZMNU2(JLO)-ZMNU2(JHI)))
                      AERNU2(J)=EXPINT(AERNU2(JHI),AERNU2(JLO),FAC)
                  ENDIF
              ENDIF
          ENDDO
 120      CONTINUE

          DO I=1,LAYDIM
              DENSTY(12,I)=0.
          ENDDO
          LAER1=.FALSE.
          LAER2=.TRUE.
          LAER3=.FALSE.
          LAER4=.FALSE.
          CALL AERMRG(NAER2,ZMNU2,AERNU2,ML,LAER1,LAER2,LAER3,LAER4)
      ENDIF
      IF(NAER3.GT.0)THEN

          DO I=1,NMRG
              IF(ALTMRG(I).GE.ZMNU3(NAER3))GOTO 130
              IF(ALTMRG(I).GT.ZMNU3(1))THEN
                  CALL FINDEX(ZMNU3,NAER3,ALTMRG(I),INDX1,JLO)
                  IF(INDX1.EQ.-1)THEN

!                     NO MATCH W/ ZMNU; SO INSERT.
                      NAER3=NAER3+1
                      J=JLO+1
                      JHI=JLO+2
                      DO II=NAER3,JHI,-1
                          IIMIN1=II-1
                          ZMNU3(II)=ZMNU3(IIMIN1)
                          AERNU3(II)=AERNU3(IIMIN1)
                      ENDDO
                      ZMNU3(J)=ALTMRG(I)
                      FAC=SNGL((ZMNU3(J)-ZMNU3(JHI))                    &
     &                  /(ZMNU3(JLO)-ZMNU3(JHI)))
                      AERNU3(J)=EXPINT(AERNU3(JHI),AERNU3(JLO),FAC)
                  ENDIF
              ENDIF
          ENDDO
 130      CONTINUE

          LAER1=.FALSE.
          LAER4=.FALSE.
          IF(IHAZE.NE.6)THEN
              DO I=1,LAYDIM
                  DENSTY(13,I)=0.
              ENDDO
              LAER2=.FALSE.
              LAER3=.TRUE.
          ELSE
!             IHAZE.EQ.6 MOVES WHAT IS NORMALLY AEROSOL 3
!             TO AEROSOL 2 REGION.
              DO I=1,LAYDIM
                  DENSTY(12,I)=0.
              ENDDO
              LAER2=.TRUE.
              LAER3=.FALSE.
          ENDIF
          CALL AERMRG(NAER3,ZMNU3,AERNU3,ML,LAER1,LAER2,LAER3,LAER4)
      ENDIF
      IF(NAER4.GT.0)THEN

          DO I=1,NMRG
              IF(ALTMRG(I).GE.ZMNU4(NAER4))GOTO 140
              IF(ALTMRG(I).GT.ZMNU4(1))THEN
                  CALL FINDEX(ZMNU4,NAER4,ALTMRG(I),INDX1,JLO)
                  IF(INDX1.EQ.-1)THEN

!                     NO MATCH W/ ZMNU; SO INSERT.
                      NAER4=NAER4+1
                      J=JLO+1
                      JHI=JLO+2
                      DO II=NAER4,JHI,-1
                          IIMIN1=II-1
                          ZMNU4(II)=ZMNU4(IIMIN1)
                          AERNU4(II)=AERNU4(IIMIN1)
                      ENDDO
                      ZMNU4(J)=ALTMRG(I)
                      FAC=SNGL((ZMNU4(J)-ZMNU4(JHI))                    &
     &                  /(ZMNU4(JLO)-ZMNU4(JHI)))
                      AERNU4(J)=EXPINT(AERNU4(JHI),AERNU4(JLO),FAC)
                  ENDIF
              ENDIF
          ENDDO
  140     CONTINUE

          LAER1=.FALSE.
          LAER2=.FALSE.
          IF(IHAZE.NE.6)THEN
              DO I=1,LAYDIM
                  DENSTY(14,I)=0.
              ENDDO
              LAER3=.FALSE.
              LAER4=.TRUE.
          ELSE
!             IHAZE.EQ.6 MOVES WHAT IS NORMALLY AEROSOL 4
!             TO AEROSOL 3 REGION
              DO I=1,LAYDIM
                  DENSTY(13,I)=0.
              ENDDO
              LAER3=.TRUE.
              LAER4=.FALSE.
          ENDIF
          CALL AERMRG(NAER4,ZMNU4,AERNU4,ML,LAER1,LAER2,LAER3,LAER4)
      ENDIF
      END

      SUBROUTINE ZMRG(A,NA,B,NB,ABMRG,NABMRG)
      IMPLICIT NONE
      INTEGER NA,NB,NABMRG,I,J,JLO,INDX
      DOUBLE PRECISION A(NA),B(NB),ABMRG(*)

!     INPUT:  A(*) CONTAINS REALS; ITS INPUT LENGTH IS NA
!     INPUT:  B(*) CONTAINS REALS; INPUT LENGTH NB
!     BOTH INCREASE W/ ARRAY INDEX

!     MERGE A AND B TO CREATE THE MERGED ARRAY ABMRG OF SIZE
!     NABMRG BY ELIMINATING REDUNDANTS AND PRESERVING ORDER

      NABMRG=NA
      DO I=1,NA
          ABMRG(I)=A(I)
      ENDDO

      DO I=1,NB
          CALL FINDEX(ABMRG,NABMRG,B(I),INDX,JLO)
          IF(INDX.EQ.-1)THEN

!             NO MATCH, SO INSERT BUT FIRST MOVE
              NABMRG=NABMRG+1
              DO J=NABMRG,JLO+2,-1
                  ABMRG(J)=ABMRG(J-1)
              ENDDO
              ABMRG(JLO+1)=B(I)
          ENDIF
      ENDDO
      RETURN
      END

      SUBROUTINE ZMRG2(A,NA,B,NB)
      IMPLICIT NONE
      INTEGER NA,NB,I,J,JLO,INDX
      DOUBLE PRECISION A(*),B(NB)

!     INPUT:  A(*) CONTAINS REALS; ITS INPUT LENGTH IS NA
!     INPUT:  B(*) CONTAINS REALS; INPUT LENGTH NB
!     BOTH INCREASE W/ ARRAY INDEX

!     MERGE A AND B TO CREATE THE MERGED ARRAY A OF SIZE NA BY
!     ELIMINATING REDUNDANTS AND PRESERVING ORDER

      DO I=1,NB
          CALL FINDEX(A,NA,B(I),INDX,JLO)
          IF(INDX.EQ.-1)THEN

!             NO MATCH, SO INSERT BUT FIRST MOVE
              NA=NA+1
              DO J=NA,JLO+2,-1
                  A(J)=A(J-1)
              ENDDO
              A(JLO+1)=B(I)
          ENDIF
      ENDDO
      END
