      SUBROUTINE CD1(NLOS,TPTMPS,APPREF,LMODEL,LYMOL,INITC1,KPRINT,CLRT)
!GLUT SUBROUTINE CD1(NLOS,TPTMPS,APPREF,LMODEL,LYMOL,INITC1,KPRINT,CLRT,&
!    &  DODB)

!     PROCESS CARD1 INPUTS:
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     INPUT ARGUMENTS:
!       INITC1   LOGICAL FLAG, .TRUE. FOR INITIAL PASS THROUGH CD1.
!       CLRT     NAME OF COOLING RATE FILE.
      CHARACTER CLRT*(NAMLEN)

!     OUTPUT ARGUMENTS:
!       NLOS     NUMBER OF OBSERVER ZENITH ANGLES.
!       TPTMPS   INPUT TARGET-PIXEL SURFACE TEMPERATURE [K].
!       APPREF   FLAG, TRUE FOR APPREF REFLECTANCE INSTEAD OF RADIANCE.
!       LMODEL   FLAG, .TRUE. IF MODEL ATMOSPHERE OPTION IS USED.
!       LYMOL    FLAG, .TRUE. IF INTERNAL Y-SPECIES OPTION IS USED.
!       INITC1   LOGICAL FLAG, .TRUE. FOR INITIAL PASS THROUGH CD1.
!       KPRINT   LOGICAL FLAG DICTATING K-DISTRIBUTION OUTPUT.
      INTEGER NLOS
      REAL TPTMPS
      LOGICAL APPREF,LMODEL,LYMOL,INITC1,KPRINT
!GLUT CHARACTER DODB*7

!     COMMONS:
      INCLUDE 'IFIL.h'

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

!     /JM1/
!       CODE     CODE INPUT FLAG: 'F' OR 'L' FOR LOWTRAN BAND MODEL.
!                                 'T' OR 'M' FOR MODTRAN BAND MODEL.
!                                 'C' OR 'K' FOR CORRELATED-K MODEL.
!       SPEED    COMPUTATIONAL SPEED FLAG USED WITH CK MODEL ['F' OR '2'
!                FOR FAST (12 K'S), 'M' OR '1' FOR MODERATE (17 K'S),
!                AND 'S' OR '0' FOR SLOW (33 K'S); 'S' IS DEFAULT]
!       SURREF   SURFACE REFLECTANCE CHARACTER STRING.
!                IF FIRST NON-BLANK CHARACTER IS "B",
!                  BI-DIRECTIONAL REFLECTANCE DISTRIBUTION
!                  FUNCTION (BRDF) DATA IS READ IN.
!                IF FIRST NON-BLANK CHARACTER IS "L" OR "-"
!                  SURFACE IS MODELED AS A LAMBERTIAN REFLECTOR AND
!                  SPECTRAL ALBEDO IS READ FROM FILE "DATA/spec_alb.dat"
!                OTHERWISE, THE CHARACTER STRING IS ASSUMED TO CONTAIN
!                  A SPECTRALLY INDEPENDENT VALUE FOR SURFACE ALBEDO.
      CHARACTER CODE*1,SPEED*1,SURREF*7
      COMMON/JM1/CODE,SPEED,SURREF

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

!     /CARD1A/
      INTEGER M4,M5,M6,MDEF,IRD1,IRD2
      COMMON/CARD1A/M4,M5,M6,MDEF,IRD1,IRD2

!     /TITL/
      CHARACTER HHAZE(16)*20,HSEASN(2)*20,HVULCN(8)*20,                 &
     &  HMET(2)*20,HMODEL(0:8)*20
      COMMON/TITL/HHAZE,HSEASN,HVULCN,HMET,HMODEL

!     SAVED COMMONS:
      SAVE /TITL/

!     FUNCTIONS:
!       LENSTR   RETURNS STRING LENGTH AFTER TRIMMING LEADING BLANKS.
      INTEGER LENSTR

!     LOCAL VARIABLES:
!       LOPEN    LOGICAL FLAG, .TRUE. IF FILE IS OPEN.
!       LYMOLC   EQUALS '+' IF INTERNAL Y-SPECIES OPTION IS USED.
!       BINARY   EQUALS 't' OR 'T' FOR BINARY SPECTRAL OUTPUT.
!       CKPRNT   EQUALS 't' OR 'T' FOR K-DISTRIBUTION OUTPUT.
!       LNCLRT   LENGTH OF COOLING RATE FILE NAME.
      LOGICAL LOPEN
      CHARACTER LYMOLC*1,BINARY*1,CKPRNT*1
      INTEGER LNCLRT

!     INPUT CARD1:
      IF(LJMASS)THEN
          CALL INITCD('CARD1')
!GLUT     DODB='       '
      ELSE
          WRITE(IPR,'(/(A))')                                           &
     &      ' ********************************************************',&
     &      ' *                                                      *',&
     &      ' *  MODTRAN5:  Developmental Version 3.5  February 2008 *',&
     &      ' *                                                      *',&
     &      ' *  Developed collaboratively by SPECTRAL SCIENCES,     *',&
     &      ' *  INC. (www.spectral.com) and the AIR FORCE RESEARCH  *',&
     &      ' *  LABORATORY (www.vs.afrl.af.mil/Division/VSBYB)      *',&
     &      ' *                                                      *',&
     &      ' ********************************************************'
          READ(IRD,'(4A1,I1,11I5,A1,I4,F8.0,2A7)')                      &
     &      CODE,SPEED,BINARY,LYMOLC,MODEL,ITYPE,IEMSCT,IMULT,M1,       &
     &      M2,M3,M4,M5,M6,MDEF,I_RD2C,CKPRNT,NOPRNT,TPTMPS,SURREF

      !WRITE(*,*) 'CODE = ',CODE
      !WRITE(*,*) 'SPEED = ',SPEED
      !WRITE(*,*) 'BINARY = ',BINARY
      !WRITE(*,*) 'LYMOLC = ',LYMOLC
      !WRITE(*,*) 'MODEL = ',MODEL
      !WRITE(*,*) 'ITYPE = ',ITYPE
      !WRITE(*,*) 'IEMSCT = ',IEMSCT
      !WRITE(*,*) 'IMULT = ',IMULT
      !WRITE(*,*) 'M1 = ',M1
      !WRITE(*,*) 'M2 = ',M2
      !WRITE(*,*) 'M3 = ',M3
      !WRITE(*,*) 'M4 = ',M4
      !WRITE(*,*) 'M5 = ',M5
      !WRITE(*,*) 'M6 = ',M6
      !WRITE(*,*) 'MDEF = ',MDEF
      !WRITE(*,*) 'I_RD2C = ',I_RD2C
      !WRITE(*,*) 'CKPRNT = ',CKPRNT
      !WRITE(*,*) 'NPRNT = ',NOPRNT
      !WRITE(*,*) 'TPTMPS = ',TPTMPS
      !WRITE(*,*) 'SURREF = ',SURREF
      

!GLUT&      M2,M3,M4,M5,M6,MDEF,I_RD2C,CKPRNT,NOPRNT,TPTMPS,SURREF,DODB
      ENDIF
      IF(TPTMPS.LT.0.)TPTMPS=0.
      BINOUT=BINARY.EQ.'T' .or. BINARY.EQ.'t'
      THERML=IEMSCT.NE.4
      KPRINT=CKPRNT.EQ.'T' .or. CKPRNT.EQ.'t'
      APPREF=IEMSCT.EQ.-2 .OR. IEMSCT.EQ.-3 .OR. IEMSCT.EQ.-4
      IF(APPREF)THEN
          IEMSCT=ABS(IEMSCT)
          WRITE(IPR,'(/A)')'*** Input IEMSCT value is negative:  APPAR' &
     &      //'ENT REFLECTANCE will be output instead of RADIANCE! ***'
      ENDIF
      IF(IEMSCT.EQ.4)IEMSCT=2
      IF(IMULT.NE.0)THEN
          IF(IEMSCT.EQ.0 .OR. IEMSCT.EQ.3 .OR. ITYPE.EQ.1)THEN
              IF(.NOT.LJMASS)WRITE(IPR,'(/A)')' Warning from CD1: '     &
     &          //' Multiple scattering has been turned off.'
              IMULT=0
          ELSEIF(.NOT.LJMASS)THEN
              WRITE(IPR,'(/A)')                                         &
     &          ' Calculations will be done using multiple scattering.'
          ENDIF
      ENDIF
      CALL UPCASE(CODE)
      MODTRN=CODE.NE.'L' .AND. CODE.NE.'F'
      LYMOL=LYMOLC.EQ.'+'

!     CHECK FOR MULTIPLE PATHS:
      IF(ITYPE.LT.0)THEN
          NLOS=-ITYPE
          IF(NLOS.GT.1 .AND. IMULT.EQ.0)THEN
              WRITE(IPR,'(/A,/15X,A)')'Error in CD1:  The number'       &
     &          //' of Lines-Of-Sight exceeds 1 (ITYPE < -1),',         &
     &          'but multiple scattering is off (IMULT=0)!'
              STOP 'CD1 Error: Multiple LOS with no multiple scattering'
          ELSEIF(NLOS.GT.MLOS)THEN
              WRITE(IPR,'(/A,I5)')'Error in CD1:  Too many lines-of'//  &
     &          '-sight.  Increase parameter MLOS in PARAMS.h to',NLOS
              STOP 'CD1 Error: Too many lines-of-sight to be processed!'
          ENDIF
          ITYPE=2
      ELSE
          NLOS=1
      ENDIF

!     CHECK RANGE FOR M1 THROUGH M6 INPUTS:
      IF(M1.LT.0 .OR. M1.GT.6)THEN
          WRITE(IPR,'(/A,I5,A,/18X,A)')'Warning from CD1:  Improper'//  &
     &      ' input value for M1 (=',M1,');',' M1 is being reset to 0.'
          M1=0
      ENDIF
      IF(M2.LT.-6 .OR. M2.GT.6)THEN
          WRITE(IPR,'(/A,I5,A,/18X,A)')'Warning from CD1:  Improper'//  &
     &      ' input value for M2 (=',M2,');',' M2 is being reset to 0.'
          M2=0
      ENDIF
      IF(M3.LT.0 .OR. M3.GT.6)THEN
          WRITE(IPR,'(/A,I5,A,/18X,A)')'Warning from CD1:  Improper'//  &
     &      ' input value for M3 (=',M3,');',' M3 is being reset to 0.'
          M3=0
      ENDIF
      IF(M4.LT.0 .OR. M4.GT.6)THEN
          WRITE(IPR,'(/A,I5,A,/18X,A)')'Warning from CD1:  Improper'//  &
     &      ' input value for M4 (=',M4,');',' M4 is being reset to 0.'
          M4=0
      ENDIF
      IF(M5.LT.0 .OR. M5.GT.6)THEN
          WRITE(IPR,'(/A,I5,A,/18X,A)')'Warning from CD1:  Improper'//  &
     &      ' input value for M5 (=',M5,');',' M5 is being reset to 0.'
          M5=0
      ENDIF
      IF(M6.LT.0 .OR. M6.GT.6)THEN
          WRITE(IPR,'(/A,I5,A,/18X,A)')'Warning from CD1:  Improper'//  &
     &      ' input value for M6 (=',M6,');',' M6 is being reset to 0.'
          M6=0
      ENDIF
      LMODEL=MODEL.GE.1 .AND. MODEL.LE.6
      IF(MDEF.LT.-1 .OR. MDEF.GT.2 .OR. (MDEF.EQ.-1 .AND. LMODEL))THEN
          WRITE(IPR,'(/A,I5,A,/18X,A)')                                 &
     &      'Warning from CD1:  Improper input value for MDEF (=',      &
     &      MDEF,');',' MDEF is being reset to 0.'
          MDEF=0
      ENDIF
      IF(LMODEL)THEN
          IF(M1.EQ.0)M1=MODEL
          IF(M2.EQ.0)M2=MODEL
          IF(M3.EQ.0)M3=MODEL
          IF(M4.EQ.0)M4=MODEL
          IF(M5.EQ.0)M5=MODEL
          IF(M6.EQ.0)M6=MODEL
          IF(MDEF.EQ.0)MDEF=1
      ELSEIF(I_RD2C.NE.1)THEN
          IF(INITC1)THEN
              WRITE(IPR,'(/A,/18X,A)')'Warning from CD1:  The user-'//  &
     &          'specified atmosphere option has been selected, but',   &
     &          'I_RD2C on CARD 1 is not equal to 1.'                   &
     &          //'  It is being set to 1.'
              I_RD2C=1
          ELSE
              WRITE(IPR,'(/A,/18X,A)')'Warning from CD1:  The user-'//  &
     &          'specified atmosphere option has been selected, and',   &
     &          'I_RD2C on CARD 1 is not equal to 1.  The atmosphere'// &
     &          ' from the previous run will be reused.'
          ENDIF
      ENDIF
      INITC1=.FALSE.
      IF(LJMASS)THEN
          NOPRNT=3
          NPR=3
          RETURN
      ENDIF
      IF(NOPRNT.GT.3)THEN
          WRITE(IPR,'(/A)')' Input NOPRNT reduced to 3 (minimum output)'
          NOPRNT=3
          NPR=3
      ELSEIF(NOPRNT.LE.-2)THEN
          IF(NOPRNT.LT.-2)THEN
              WRITE(IPR,'(/A)')                                         &
     &          ' Input NOPRNT increased to -2 (maximum output)'
              NOPRNT=-2
          ENDIF

!         OPEN COOLING RATE FILE IF MULTIPLE SCATTERING
!         IS ON AND THE FILE IS NOT ALREADY OPEN:
          IF(IMULT.EQ.0)THEN
              NPR=-1
              WRITE(IPR,'(/A)')' Since Multiple Scattering is off '//   &
     &          '(IMULT=0), no cooling rate data will be generated.'
          ELSE
              INQUIRE(ICR,OPENED=LOPEN)
              IF(.NOT.LOPEN)THEN
                  LNCLRT=LENSTR(CLRT)
                  CALL OPNFL(ICR,0,                                     &
     &              CLRT(1:LNCLRT),'UNKNOWN','FORMATTED','CD1')
              ENDIF
              NPR=-2
          ENDIF
      ELSE
          NPR=NOPRNT
      ENDIF
      WRITE(IPR,'(/A,3A1,I2,11I5,A1,I4,F8.3,A7)')                       &
     &  ' CARD 1  *****',CODE,SPEED,BINARY,MODEL,ITYPE,IEMSCT,          &
     &  IMULT,M1,M2,M3,M4,M5,M6,MDEF,I_RD2C,CKPRNT,NOPRNT,TPTMPS,SURREF
      IF(IEMSCT.EQ.0)THEN
          WRITE(IPR,'(/A)')' PROGRAM WILL COMPUTE TRANSMITTANCE'
      ELSEIF(IEMSCT.EQ.1)THEN
          WRITE(IPR,'(/A)')' PROGRAM WILL COMPUTE RADIANCE'
      ELSEIF(IEMSCT.EQ.2)THEN
          WRITE(IPR,'(/A)')                                             &
     &      ' PROGRAM WILL COMPUTE RADIANCE + SOLAR SCATTER'
      ELSEIF(IEMSCT.EQ.3)THEN
          WRITE(IPR,'(/A)')' PROGRAM WILL COMPUTE'//                    &
     &      ' DIRECTLY TRANSMITTED SOLAR IRRADIANCE'
      ENDIF
      IF(MODEL.EQ.0)RETURN
      HMODEL(7)='USER PROFILES       '
      HMODEL(8)='HYDROSTATIC EQUATION'
      WRITE(IPR,'(/A)')' ATMOSPHERIC MODEL'
      IF(M1.EQ.0)THEN
          WRITE(IPR,'(9X,A,I5,5X,A)')'TEMPERATURE =',MODEL,HMODEL(MODEL)
      ELSE
          WRITE(IPR,'(9X,A,I5,5X,A)')'TEMPERATURE =',M1,HMODEL(M1)
      ENDIF
      IF(M2.EQ.0)THEN
          WRITE(IPR,'(9X,A,I5,5X,A)')'WATER VAPOR =',MODEL,HMODEL(MODEL)
      ELSE
          WRITE(IPR,'(9X,A,I5,5X,A)')'WATER VAPOR =',M2,HMODEL(ABS(M2))
      ENDIF
      IF(M3.EQ.0)THEN
          WRITE(IPR,'(9X,A,I5,5X,A)')'OZONE       =',MODEL,HMODEL(MODEL)
      ELSE
          WRITE(IPR,'(9X,A,I5,5X,A)')'OZONE       =',M3,HMODEL(M3)
      ENDIF
      WRITE(IPR,'(22X,4(A,I6))')                                        &
     &  'M4 =',M4,' M5 =',M5,' M6 =',M6,' MDEF =' ,MDEF

!     RETURN TO DRIVER
      RETURN
      END
