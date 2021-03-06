      SUBROUTINE FLAYZ(ML,ICLD,GNDALT,IVSA)

!     THIS ROUTINE DETERMINES THE LAYER BOUNDARIES WHEN
!     ONE OF THE MODEL ATMOSPHERES IS SELECTED.
      IMPLICIT NONE

!     INCLUDE PARAMETERS.
      INCLUDE 'PARAMS.h'
      INTEGER NDFLT
      PARAMETER(NDFLT=36)

!     INPUT ARGUMENT:
      INTEGER ML,ICLD,IVSA
      DOUBLE PRECISION GNDALT

!     COMMONS:

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

!     /CARD2A/
      INTEGER NCRALT,NCRSPC
      DOUBLE PRECISION CTHIK,CALT
      REAL CEXT,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP
      COMMON/CARD2A/CTHIK,CALT,CEXT,NCRALT,NCRSPC,CWAVLN,CCOLWD,        &
     &  CCOLIP,CHUMID,ASYMWD,ASYMIP

!     /ZVSALY/
      DOUBLE PRECISION ZVSA
      REAL RHVSA,AHVSA
      INTEGER IHVSA
      COMMON/ZVSALY/ZVSA(10),RHVSA(10),AHVSA(10),IHVSA(10)
      INCLUDE 'IFIL.h'

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     DECLARE LOCAL VARIABLES.
      INTEGER I,ILAST,INEXT,NEXT,IHI,J,IM1,ICIR(4)
      DOUBLE PRECISION DALT,CLDD,TOL,ZCIR(4)

!     DATA:
      DOUBLE PRECISION ZDFLT(NDFLT)
      DATA ZDFLT/ 0.D0, 1.D0, 2.D0, 3.D0, 4.D0, 5.D0, 6.D0, 7.D0, 8.D0, &
     &            9.D0,10.D0,11.D0,12.D0,13.D0,14.D0,15.D0,16.D0,17.D0, &
     &           18.D0,19.D0,20.D0,21.D0,22.D0,23.D0,24.D0,25.D0,30.D0, &
     &           35.D0,40.D0,45.D0,50.D0,55.D0,60.D0,70.D0,80.D0,100.D0/
      IF(IVSA.EQ.1)THEN

!         VERTICAL STRUCTURE ALGORITHM (VSA) OPTION IS ON.
!         FIRST 9 BOUNDARY ALTITUDES WERE DEFINED IN ROUTINE VSA.
          DO I=1,9
              ZM(I)=ZVSA(I)
          ENDDO

!         ADD A BOUNDARY ALTITUDE 10 METERS ABOVE ZVSA(9).
          ML=10
          ZM(10)=ZM(9)+.01D0

!         FIND THE FIRST DEFAULT ALTITUDE HALF A METER ABOVE ZM(10).
          ILAST=0
          DO INEXT=1,NDFLT
              IF(ZDFLT(INEXT).GE.ZM(10)+.0005D0)GOTO 10
              ILAST=INEXT
          ENDDO

!         NO DEFAULT ALTITUDES ABOVE ZM(10)+.0005D0.
          RETURN

!         ADD DEFAULT ALTITUDES.
   10     CONTINUE
          IF(ML+NDFLT-ILAST.GT.LAYDIM)THEN
              WRITE(IPR,'(2A,I3,/8X,2A,I3,A)')                          &
     &          ' ERROR:  VERTICAL STRUCTURE ALGORITHM INCREASES',      &
     &          ' LAYER BOUNDARY NUMBER ABOVE',LAYDIM,' INCREASE',      &
     &          ' PARAMETER LAYDIM TO',ML+NDFLT-ILAST,' IN PARAMS.h'
              IF(LJMASS)CALL WRTBUF(FATAL)
              STOP 'INCREASE PARAMETER LAYDIM'
          ENDIF
          DO I=INEXT,NDFLT
              ML=ML+1
              ZM(ML)=ZDFLT(I)
          ENDDO
          RETURN
      ENDIF

!     ASSIGN DEFAULT ALTITUDES
      ML=NDFLT
      DO I=1,NDFLT
          ZM(I)=ZDFLT(I)
      ENDDO

!     IF THE INPUT GROUND ALTITUDE IS LESS THAN 6 KM ABOVE
!     SEA LEVEL, SCALE BOTTOM 5 LAYER BOUNDARIES.
      IF(GNDALT.LT.ZM(7))THEN
          DALT=(ZM(7)-GNDALT)/(ZM(7)-ZM(1))
          DO I=6,1,-1
              ZM(I)=GNDALT+DALT*(ZM(I)-ZM(1))
          ENDDO
      ENDIF
      IF(ICLD.NE.18 .AND. ICLD.NE.19)RETURN

!     CIRRUS CLOUD.  DEFINE CLOUD BOUNDARIES.
      CLDD=CTHIK/10
      ZCIR(1)=CALT-CLDD/2
      IF(ZCIR(1).LE.GNDALT)ZCIR(1)=GNDALT
      ZCIR(2)=ZCIR(1)+CLDD
      ZCIR(3)=ZCIR(1)+CTHIK
      ZCIR(4)=ZCIR(3)+CLDD
      IF(ZCIR(4).GE.ZM(NDFLT))THEN
          WRITE(IPR,'(2A,/8X,2A)')' ERROR:  CIRRUS LAYER',              &
     &      ' IS ABOVE THE TOP OF THE MODEL ATMOSPHERE.',               &
     &      ' IF YOU REALLY WANT A CLOUD AT THIS ALTITUDE,',            &
     &      ' INPUT A USER-DEFINED ATMOSPHERE.'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP 'CIRRUS CLOUD ABOVE TOP OF ATMOSPHERE'
      ENDIF

!     CHECK ZM DIMENSION
      ML=NDFLT+4
      IF(LAYDIM.LT.ML)THEN
          WRITE(IPR,'(3A,I3,/14X,A,I3,A)')' ERROR:  CIRRUS',            &
     &      ' LAYER INCREASES LAYER BOUNDARY NUMBER ABOVE',LAYDIM,      &
     &      ' INCREASE PARAMETER LAYDIM TO',ML,' IN PARAMS.h'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP 'INCREASE PARAMETER LAYDIM'
      ENDIF

!     COMBINE DEFAULT AND CIRRUS LAYERS STARTING FROM THE TOP.
      NEXT=ML
      IHI=NDFLT
      DO J=4,1,-1
          DO I=IHI,1,-1
              IF(ZCIR(J).GE.ZM(I))GOTO 20
              ZM(NEXT)=ZM(I)
              NEXT=NEXT-1
          ENDDO
   20     CONTINUE
          ZM(NEXT)=ZCIR(J)
          ICIR(J)=NEXT
          NEXT=NEXT-1
          IHI=I
      ENDDO

!     SET TOLERANCE FOR MERGING LAYERS
      TOL=MIN(CLDD/2,.0005D0)

!     CHECK FOR MERGING
      DO J=4,1,-1
          IF(ZM(ICIR(J)+1)-ZM(ICIR(J)).LT.TOL)THEN
              IM1=ICIR(J)+1
              DO I=ICIR(J)+2,ML
                  ZM(IM1)=ZM(I)
                  IM1=I
              ENDDO
              ML=ML-1
          ELSEIF(ZM(ICIR(J))-ZM(ICIR(J)-1).LT.TOL)THEN
              IM1=ICIR(J)-1
              DO I=ICIR(J),ML
                  ZM(IM1)=ZM(I)
                  IM1=I
              ENDDO
              ML=ML-1
          ENDIF
      ENDDO
      RETURN
      END
