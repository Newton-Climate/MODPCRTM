      SUBROUTINE RDEXA

!     READ IN USER DEFINED EXTINCTION ABSORPTION COEFFICIENTS AND
!     ASYMMETRY PARAMETERS
      IMPLICIT NONE

      INCLUDE 'PARAMS.h'

!     COMMONS:
      INCLUDE 'IFIL.h'
      INCLUDE 'BASE.h'
      INTEGER IREG,IREGC
      DOUBLE PRECISION ALTB
      COMMON/CARD2D/IREG(4),ALTB(4),IREGC(4)
      INTEGER IK,IHC,JLO,JHI
      CHARACTER AERNAM(18)*4
      REAL VX(47),EXTEMP(4,47),ABTEMP(4,47),ASTEMP(4,47)
      REAL Y1,Y2,X1,X2,X
      REAL XNTRP,LNTRP
      LOGICAL LOGINT

      INTEGER ISPC,ISPC3,LIMSPC,IOS
      REAL VX0
      COMMON/EXTWAV/VX0(NWAVLN)

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      SAVE /EXTWAV/
      EXTERNAL DEVCBD,EXTDT1,EXTDT2

      DATA VX /                                                         &
     &     .2000,   .3000,   .3371,   .5500,   .6943,  1.0600,  1.5360, &
     &     2.0000,  2.2500,  2.5000,  2.7000,  3.0000,  3.3923,  3.7500,&
     &     4.5000,  5.0000,  5.5000,  6.0000,  6.2000,  6.5000,  7.2000,&
     &     7.9000,  8.2000,  8.7000,  9.0000,  9.2000, 10.0000, 10.5910,&
     &     11.0000, 11.5000, 12.5000, 14.8000, 15.0000, 16.4000,17.2000,&
     &     18.5000, 21.3000, 25.0000, 30.0000, 40.0000, 50.0000,60.0000,&
     &     80.0000, 100.000, 150.000, 200.000, 300.000/

      IF(LJMASS)THEN
          CALL INITCD( 'CARD2D' )
      ELSE
          READ(IRD,'(4I5)')(IREG(IK),IK=1,4)
          WRITE(IPR,'(A,4I5)')'1 CARD 2D *****',(IREG(IK),IK=1,4)
      ENDIF

      LOGINT=.FALSE.
      DO IHC=1,4

          IF(IREG(IHC).NE.0)THEN
              LOGINT=.TRUE.
              IF(LJMASS)THEN
                  CALL INITCD( 'CARD2D1' )
              ELSE
                  READ(IRD,'(E10.3,18A4)') AWCCON(IHC),AERNAM
                  WRITE(IPR,'(/A,1P,E10.3,18A4)')' CARD 2D1 ****'//     &
     &              ' EQUIVALENT WATER = ',AWCCON(IHC),AERNAM
                  WRITE(IPR,'(/A)')'0 CARD 2D2 ****'

!                 NWAVLN WAS REPLACED BY 47 FOR BACKWARD COMPATIBILITY
!                 AS NWAVLN IS NOW 788 (SEE PARAMS.h).
                  DO ISPC3=1,47,3
                      LIMSPC=MIN(47,ISPC3+2)
                      READ(IRD,'(3(6X,2F7.5,F6.4))',IOSTAT=IOS)         &
     &                  (EXTEMP(IHC,ISPC),ABTEMP(IHC,ISPC),             &
     &                  ASTEMP(IHC,ISPC),ISPC=ISPC3,LIMSPC)
                      IF(IOS.NE.0)THEN
                          WRITE(IPR,'(/A,I3,A)')'Error in RDEXA:  '//   &
     &                    'Unable to read line',1+ISPC3/3,' of CARD 2D2'
                          STOP 'Error in RDEXA: Unable to read CARD 2D2'
                      ENDIF
                      WRITE(IPR,'(3(F6.2,1X,2(F7.5,1X),F6.4,1X))')      &
     &                  (VX(ISPC),EXTEMP(IHC,ISPC),ABTEMP(IHC,ISPC),    &
     &                  ASTEMP(IHC,ISPC),ISPC=ISPC3,LIMSPC)
                      DO ISPC=ISPC3,LIMSPC
                          IF(EXTEMP(IHC,ISPC).LT.0.)THEN
                              WRITE(IPR,'(/A)')'Error in RDEXA:  Neg'// &
     &                          'ative CARD 2D2 extinction coefficient.'
                              STOP 'CARD 2D2 extinction coef < 0.'
                          ELSEIF(ABTEMP(IHC,ISPC).LT.0.)THEN
                              WRITE(IPR,'(/A)')'Error in RDEXA:  Neg'// &
     &                          'ative CARD 2D2 absorption coefficient.'
                              STOP 'CARD 2D2 absorption coef < 0.'
                          ELSEIF(ABTEMP(IHC,ISPC)                       &
     &                       .GT.EXTEMP(IHC,ISPC))THEN
                              WRITE(IPR,'(/A)')'Error in RDEXA:  CARD'  &
     &                          //' 2D2 absorption exceeds extinction.'
                              STOP 'CARD 2D2 absorption > extinction.'
                          ELSEIF(ABS(ASTEMP(IHC,ISPC)).GE.1.)THEN
                              WRITE(IPR,'(/A)')'Error in RDEXA:  As'//  &
     &                          'ymmetry factor magnitude exceeds one.'
                              STOP 'CARD 2D2 |asymmetry factor| > 1.'
                          ENDIF
                      ENDDO
                  ENDDO
              ENDIF

!             AFTER THE 47 POINTS ARE READ, INTERPOLATE FOR VX0 WITH
!             788 POINTS.  BEFORE THE MODIFICATIONS, VX WAS READ BUT
!             NOT USED BECAUSE THE INFORMATION WAS IN COMMON ELSEWHERE.
!             FOR BACKWARD COMPATIBILITY, THE 47 POINTS WILL CORRESPOND
!             TO THE OLD VX ARRAY, EXPLICITLY STATED IN THE ABOVE DATA
!             STATEMENT.  BECAUSE THERE IS A NEED FOR VX0 IN THIS
!             ROUTINE, THE COMMON BLOCK /EXTWAV/ IS NOW IN THIS ROUTINE.
          ENDIF
      ENDDO
      IF(.NOT.LOGINT)RETURN

!     INTERPOLATE
      DO ISPC=1,NWAVLN
          CALL HUNT(VX,47,VX0(ISPC),JLO,JHI)
          X=1.E4/VX0(ISPC)
          DO IHC=1,4
              IF(IREG(IHC).NE.0)THEN
                  IF(JLO.EQ.0 )THEN
                      EXTC(IHC,1)=EXTEMP(IHC,1)
                      ABSC(IHC,1)=ABTEMP(IHC,1)
                      ASYM(IHC,1)=ASTEMP(IHC,1)
                  ELSEIF(JLO.EQ.47)THEN
                      EXTC(IHC,NWAVLN)=EXTEMP(IHC,47)
                      ABSC(IHC,NWAVLN)=ABTEMP(IHC,47)
                      ASYM(IHC,NWAVLN)=ASTEMP(IHC,47)
                  ELSE
                      Y1=EXTEMP(IHC,JLO)
                      Y2=EXTEMP(IHC,JHI)
                      X1=1.E4/VX(JLO)
                      X2=1.E4/VX(JHI)
                      EXTC(IHC,ISPC)=XNTRP(Y1,Y2,X1,X2,X)
                      Y1=ABTEMP(IHC,JLO)
                      Y2=ABTEMP(IHC,JHI)
                      ABSC(IHC,ISPC)=LNTRP(Y1,Y2,X1,X2,X)
                      Y1=ASTEMP(IHC,JLO)
                      Y2=ASTEMP(IHC,JHI)
                      ASYM(IHC,ISPC)=LNTRP(Y1,Y2,X1,X2,X)
                  ENDIF
              ENDIF
          ENDDO
      ENDDO
      RETURN
      END

      REAL FUNCTION XNTRP(Y1,Y2,X1,X2,X)
      IMPLICIT NONE
      REAL Y1,Y2,X1,X2,X,LNTRP

!     EXPONENTIAL INTERPOLATION BETWEEN (X1,Y1) AND (X2,Y2).
!     FIND VALUE AT X.

      IF(Y1.EQ.0. .OR. Y2.EQ.0.0)THEN
          XNTRP=LNTRP(Y1,Y2,X1,X2,X)
      ELSEIF(Y1.EQ.Y2)THEN
          XNTRP=Y1
      ELSE
          XNTRP=Y1*EXP((X-X1)*LOG(Y2/Y1)/(X2-X1))
      ENDIF
      RETURN
      END

      REAL FUNCTION LNTRP(YA,YB,XA,XB,X)

!     LINEAR INTERPOLATION BETWEEN (XA,YA) AND (XB,YB); FIND VALUE @ X.
      IMPLICIT NONE

      REAL YA,YB,XA,XB,X
      LNTRP=YA+(YB-YA)*(X-XA)/(XB-XA)
      RETURN
      END
