      SUBROUTINE CD2C(LMODEL,LAPLUS,LENDAT,MARIC1,MARK,JPRT,ICH,        &
     &  H2OCOL,O3COL,LYMOL,REARTH)

!     PROCESS CARD2C INPUTS:
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     ARGUMENTS:
!       LMODEL   FLAG, .TRUE. IF MODEL ATMOSPHERE IS USED.
!       LAPLUS   LOGICAL FLAG FOR THE AEROSOL A+ OPTION.
!       LENDAT   LENGTH OF DATDIR STRING.
!       ICH      NUMERIC LABEL FOR AEROSOL MODEL.
!       H2OCOL   WATER COLUMN [KM GM /M3].
!       O3COL    OZONE COLUMN [KM GM /M3].
!       LYMOL    FLAG, .TRUE. IF INTERNAL Y-SPECIES IS USED.
!       REARTH   RADIUS OF THE EARTH [KM].
      LOGICAL LMODEL,LAPLUS,LYMOL
      INTEGER LENDAT,MARIC1,MARK,JPRT,ICH(4)
      REAL H2OCOL,O3COL
      DOUBLE PRECISION REARTH

!     COMMONS:
      INCLUDE 'IFIL.h'
      INCLUDE 'YPROP.h'
      INCLUDE 'YPROPC.h'

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

!     /JM1A1/
!       RESCHR   BAND MODEL RESOLUTION CHARACTER STRING.
!       DISSTR   CHARACTER STRING USED TO READ IN DISORT LOGICALS.
!       H2OSTR   VERTICAL WATER COLUMN CHARACTER STRING (IF THE
!                FIRST NON-BLANK CHARACTER IS A "G", THE WATER
!                COLUMN IN GM/CM2 FOLLOWS "G"; IF THE FIRST
!                NON-BLANK CHARACTER IS AN "A", THE WATER COLUMN
!                IN ATM-CM AT 273.15K FOLLOWS "A"; OTHERWISE THE
!                STRING CONTAINS A WATER COLUMN SCALING FACTOR).
!       O3STR    VERTICAL OZONE COLUMN CHARACTER STRING (IF THE
!                FIRST NON-BLANK CHARACTER IS A "G", THE OZONE
!                COLUMN IN GM/CM2 FOLLOWS "G"; IF THE FIRST
!                NON-BLANK CHARACTER IS AN "A", THE OZONE COLUMN
!                IN ATM-CM AT 273.15K FOLLOWS "A"; OTHERWISE THE
!                STRING CONTAINS AN OZONE COLUMN SCALING FACTOR).
!       USRSUN   USER-SPECIFIED TOP-OF-ATMOSPHERE SOLAR IRRADIANCE FILE.
!       BMROOT   PREFIX OF MOLECULAR BAND MODEL PARAMETERS FILE.
!       FILTNM   NAME OF FILTER RESPONSE FUNCTION FILE.
!       H2OAER   FLAG, TRUE IF DEFAULT AEROSOL PROPERTIES ARE REVISED
!                BASED ON WATER COLUMN SCALING.
!       DATDIR   NAME OF THE MODTRAN DATA DIRECTORY.
      CHARACTER RESCHR*2,DISSTR*3,H2OSTR*10,O3STR*10,USRSUN*(NAMLEN),   &
     &  FILTNM*(NAMLEN),BMROOT*(NAMLEN),H2OAER*1,DATDIR*(NAMLEN-LENSUN)
      COMMON/JM1A1/RESCHR,DISSTR,H2OSTR,O3STR,USRSUN,                   &
     &  FILTNM,BMROOT,H2OAER,DATDIR

!     /CARD2/
!       IHAZE    BOUNDARY LAYER AEROSOL MODEL NUMBER.
!       ISEASN   SEASON NUMBER (1=SPRING-SUMMER, 2=FALL-WINTER).
!       IVULCN   VOLCANIC AEROSOL MODEL NUMBER.
!       ICSTL    COASTAL AIRMASS MODEL NUMBER.
!       ICLD     CLOUD MODEL NUMBER.
!       IVSA     VERTICAL STRUCTURE ALGORITHM (0=OFF, 1=ON).
!       VIS      SURFACE VISIBILITY (GROUND METEOROLOGICAL RANGE) [KM].
!       WSS      CURRENT WIND SPEED (M/S).
!       WHH      24-HOUR WIND SPEED (M/S).
!       RAINRT   RAIN RATE (MM/HR)
!       LSAP     LOGICAL FLAG FOR SPECTRAL AEROSOL PROFILES INPUT.
      INTEGER IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA
      REAL VIS,WSS,WHH,RAINRT
      LOGICAL LSAP
      COMMON/CARD2/IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,                 &
     &  VIS,WSS,WHH,RAINRT,LSAP

!     /JM2/
      CHARACTER APLUS*2,CNOVAM*1,ARUSS*3
      COMMON/JM2/APLUS,CNOVAM,ARUSS

!     /GRND/
!       GNDALT   GROUND ALTITUDE, ABOVE SEA LEVEL [KM].
      DOUBLE PRECISION GNDALT
      COMMON/GRND/GNDALT

!     /COMNOV/
!       LNOVAM   LOGICAL FLAG, .TRUE. IF NOVAM AEROSOLS ARE USED.
      LOGICAL LNOVAM
      REAL EXTNOV(MNOV,MXWVLN),ABSNOV(MNOV,MXWVLN),                     &
     &  ASMNOV(MNOV,MXWVLN),WLNOV(MXWVLN)
      INTEGER NNOV,NWLNOV
      COMMON/COMNOV/LNOVAM,EXTNOV,ABSNOV,ASMNOV,WLNOV,NNOV,NWLNOV

!     /NOVDAT/
      DOUBLE PRECISION ALTNOV
      INTEGER NLNOV
      REAL RHNOV,PNOV,TNOV,DENNOV
      COMMON/NOVDAT/ALTNOV(MLNOV),NLNOV,RHNOV(MLNOV),                   &
     &  PNOV(MLNOV),TNOV(MLNOV),DENNOV(MNOV,MLNOV)

!     /CARD2A/
      INTEGER NCRALT,NCRSPC
      DOUBLE PRECISION CTHIK,CALT
      REAL CEXT,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP
      COMMON/CARD2A/CTHIK,CALT,CEXT,NCRALT,NCRSPC,CWAVLN,CCOLWD,        &
     &  CCOLIP,CHUMID,ASYMWD,ASYMIP

!     SAVED COMMONS:
      SAVE /TITL/

!     LOCAL VARIABLES:
!       YFILE1   FULL PATH Y-SPECIES LINE CENTER OR X-SECTION FILE NAME.
!       YFILE2   FULL PATH Y-SPECIES LINE TAIL FILE NAME.
!       ICLDV    SAVED VALUE OF ICLD.
!       RAINSV   SAVED VALUE OF RAINRT.
!       CPROB    PROBABILITY OF CLOUD OCCURRING [%].
!       LEXIST   LOGICAL FLAG, .TRUE. IF FILE EXISTS.
!       LOPEN    LOGICAL FLAG, .TRUE. IF FILE IS OPEN.
!       LYBMOD   LOGICAL FLAG, .TRUE. IF Y-SPECIES USES BAND MODEL.
!       IOS      STATUS OF CARD2C READ.
!       JMOLY    INDEX USED TO LOOK OVER BUILT-IN Y-SPECIES PROFILES.
      CHARACTER YNAME(MMOLY)*10,YNAM16(16)*10,YFILE1*28,YFILE2*28,CDUM*8
      REAL DUM,RAINSV,CPROB
      LOGICAL LEXIST,LOPEN,LYBMOD(MMOLY)
      INTEGER FUNIT,IDUM,NUNIT,IMOLY,ICLDSV,IHMET,IHVUL,IOS,JMOLY
      DATA YNAM16/ '       OH*','       HF*','      HCl*','      HBr*', &
     &             '       HI*','      ClO*','      OCS*','     H2CO*', &
     &             '     HOCl*','       N2*','      HCN*','    CH3Cl*', &
     &             '     H2O2*','     C2H2*','     C2H6*','      PH3*'/

      IF(LMODEL)THEN
          IRD1=0
          IRD2=0

!         IF Y FILES ARE OPEN, CLOSE THEM:
          IF(NMOLY.GT.0)THEN
              DO IMOLY=1,NMOLY
                  CLOSE(ITBY(IMOLY,1),STATUS='KEEP')
                  IF(ITBY(IMOLY,2).GT.0)                                &
     &              CLOSE(ITBY(IMOLY,2),STATUS='KEEP')
              ENDDO
          ENDIF
          IF(LYMOL)THEN

!             USING INTERNAL DEFAULT ATMOSPHERIC PROFILE TRAC(ALT,6:21)
              NMOLY=16
              NMOLYS=NMOLY
              NMOLT=NMOLXT+NMOLY
              DO IMOLY=1,NMOLY
                  CALL YFILE(IPR,YNAM16(IMOLY),YFILE1,YFILE2,RESCHR)
                  MAPY_F(IMOLY)=IMOLY
                  MAPY_D(IMOLY)=IMOLY+5
                  LYBMOD(IMOLY)=.TRUE.
                  CALL MVRGT(YNAM16(IMOLY))
                  CNAMEY(IMOLY)=YNAM16(IMOLY)(3:10)

!                 CHECK THE STATUS OF THE Y-SPECIES DATA FILE(S):
                  INQUIRE(FILE=DATDIR(1:LENDAT)//YFILE1,                &
     &              EXIST=LEXIST,OPENED=LOPEN,NUMBER=FUNIT)
                  IF(LOPEN)THEN
                      ITBY(IMOLY,1)=FUNIT
                      REWIND(FUNIT)
                  ELSEIF(LEXIST)THEN
                      ITBY(IMOLY,1)=NUNIT()
                      CALL OPNFL(ITBY(IMOLY,1),0,DATDIR(1:LENDAT)       &
     &                  //YFILE1,'OLD','FORMATTED','CD2C')
                      WRITE(IPR,'(/(A))')'Successfully opened file',    &
     &                  DATDIR(1:LENDAT)//YFILE1
                  ELSE
                      WRITE(IPR,'(/(A))')'Error in routine CD2C: '//    &
     &                  ' Band model file ',DATDIR(1:LENDAT)//YFILE1,   &
     &                  ' does not exist.  Create it!'
                      STOP 'Band Model File does not exist.'
                  ENDIF
                  INQUIRE(FILE=DATDIR(1:LENDAT)//YFILE2,                &
     &              EXIST=LEXIST,OPENED=LOPEN,NUMBER=FUNIT)
                  IF(LOPEN)THEN
                      ITBY(IMOLY,2)=FUNIT
                      REWIND(FUNIT)
                  ELSEIF(LEXIST)THEN
                      ITBY(IMOLY,2)=NUNIT()
                      CALL OPNFL(ITBY(IMOLY,2),0,DATDIR(1:LENDAT)       &
     &                  //YFILE2,'OLD','FORMATTED','CD2C')
                      WRITE(IPR,'(/(A))')'Successfully opened file',    &
     &                  DATDIR(1:LENDAT)//YFILE2
                  ELSE
                      WRITE(IPR,'(/(A))')'Error in routine CD2C: '//    &
     &                  ' band model file ',DATDIR(1:LENDAT)//YFILE2,   &
     &                  ' does not exist.  Create it!'
                      STOP 'Band Model File does not exist.'
                  ENDIF

!                 READ MOLECULAR WEIGHT IN BOTH TAIL AND CENTER FILES:
                  READ(ITBY(IMOLY,1),'(48X,F8.0)')AMWTY(IMOLY)
                  READ(ITBY(IMOLY,2),'(48X,F8.0)')AMWTY(IMOLY)

!                 Compute DOPY:
                  DOPY(IMOLY)=5.91861E-6/SQRT(AMWTY(IMOLY))
              ENDDO
          ELSE

!             THE NUMBER OF USER-DEFINED SPECIES IS INITIALIZED TO ZERO.
              NMOLY=0
              NMOLYS=0
              NMOLT=NMOLXT
          ENDIF
      ELSEIF(I_RD2C.NE.1)THEN

!         REUSE INPUT ATMOSPHERE FROM PREVIOUS RUN:
          IRD1=0
          IRD2=0
      ELSE

!         IF Y FILES ARE OPEN, CLOSE THEM:
          IF(NMOLY.GT.0)THEN
              DO IMOLY=1,NMOLY
                  CLOSE(ITBY(IMOLY,1),STATUS='KEEP')
                  IF(ITBY(IMOLY,2).GT.0)                                &
     &              CLOSE(ITBY(IMOLY,2),STATUS='KEEP')
              ENDDO
          ENDIF

!         THE NUMBER OF USER-DEFINED SPECIES IS INITIALIZED TO ZERO.
          NMOLY=0
          NMOLYS=0
          NMOLT=NMOLXT

!         CARD 2C:  USER SUPPLIED ATMOSPHERIC PROFILE
          IF(LJMASS)THEN
              CALL INITCD('CARD2C')
          ELSE
              READ(IRD,'(3I5,A20,F10.0,I5)',IOSTAT=IOS)                 &
     &          ML,IRD1,IRD2,HMODEL(MODEL),REARTH,NMOLYC
              IF(IOS.NE.0)THEN
                  WRITE(IPR,'(/A,/(16X,A))')'Error in CD2C:  Unable to' &
     &              //' read CARD 2C.  Verify that profile title does', &
     &              'not extend beyond columns 16 through 35, inclusive'
                  STOP 'Error in CD2C:  Unable to read CARD 2C!'
              ELSEIF(NMOLYC.GT.MMOLY)THEN
                  WRITE(IPR,'(/A,I3,A,/15X,A)')'Error in CD2C:  At'//   &
     &              'tempting to read in too many Y species (NMOLYC =', &
     &              NMOLYC,')',' Increase parameter MMOLY in'           &
     &              //' YPROP.h to at least the NMOLYC value.'
                  STOP 'Error in CD2C:  Too many Y specied on CARD 2C!'
              ENDIF
              IF(REARTH.LE.0.D0)THEN
                  IF(M1.EQ.1)THEN
                      REARTH=6378.39D0
                  ELSEIF(M1.EQ.4)THEN
                      REARTH=6356.91D0
                  ELSEIF(M1.EQ.5)THEN
                      REARTH=6356.91D0
                  ELSE
                      REARTH=6371.23D0
                  ENDIF
              ENDIF
              WRITE(IPR,'(/A,3I5,A20,F10.0,I5)')' CARD 2C *****',       &
     &          ML,IRD1,IRD2,HMODEL(MODEL),REARTH,NMOLYC
               !WRITE(*,*) 'ML = ',ML 
               !WRITE(*,*) ' IRD1 = ',IRD1
               !WRITE(*,*) ' IRD2 = ',IRD2
               !WRITE(*,*) ' HMODEL = ',HMODEL(MODEL)
               !WRITE(*,*) ' REARTH = ',REARTH
               !WRITE(*,*) ' NMOLYC = ',NMOLYC

              IF(NMOLYC.GT.0 .OR. LYMOL)THEN
                  IF(NMOLYC.GT.0)THEN

!                     READ IN "Y-SPECIES" NAMES (CARD 2CY):
                      READ(IRD,'((8A10))')(YNAME(IMOLY),IMOLY=1,NMOLYC)
                  ELSE

!                     ASSIGN DEFAULT "Y-SPECIES" NAMES:
                      NMOLYC=16
                      DO IMOLY=1,NMOLYC
                          YNAME(IMOLY)=YNAM16(IMOLY)
                      ENDDO
                  ENDIF
                  WRITE(IPR,'((14X,8A10))')(YNAME(IMOLY),IMOLY=1,NMOLYC)
                  DO IMOLY=1,NMOLYC
                      CALL YFILE(IPR,YNAME(IMOLY),YFILE1,YFILE2,RESCHR)
                      IF(YNAME(IMOLY)(1:1).NE.'-')THEN

                          NMOLY=NMOLY+1
                          MAPY_F(NMOLY)=IMOLY
                          CALL MVRGT(YNAME(IMOLY))
                          CNAMEY(NMOLY)=YNAME(IMOLY)(3:10)

!                         CHECK STATUS OF Y-SPECIES LINE CENTER FILE:
                          INQUIRE(FILE=DATDIR(1:LENDAT)//YFILE1,        &
     &                      EXIST=LEXIST,OPENED=LOPEN,NUMBER=FUNIT)
                          IF(LOPEN)THEN
                              ITBY(NMOLY,1)=FUNIT
                              REWIND(FUNIT)
                          ELSEIF(LEXIST)THEN
                              ITBY(NMOLY,1)=NUNIT()
                              CALL OPNFL(ITBY(NMOLY,1),0,               &
     &                          DATDIR(1:LENDAT)//YFILE1,               &
     &                          'OLD','FORMATTED','CD2C')
                              WRITE(IPR,'(/(A))')'Successfully opened'  &
     &                          //' file',DATDIR(1:LENDAT)//YFILE1
                          ELSE
                              WRITE(IPR,'(/(A))')'Error in routine CD2C'&
     &                          //':  BANDMODEL FILE ',DATDIR(1:LENDAT) &
     &                          //YFILE1,'does not exist.  Create it'   &
     &                          //' or skip species by adding a minus', &
     &                          'sign prefix to name in tape5.'
                              STOP 'Band Model file does not exist!'
                          ENDIF

!                         READ MOLECULAR WEIGHT FROM FILE HEADER:
                          READ(ITBY(NMOLY,1),'(48X,F8.0)')AMWTY(NMOLY)

!                         LINE TAIL DATA IF YFILE2 IS NOT BLANK:
                          LYBMOD(NMOLY)=YNAME(IMOLY)(10:10).EQ.'*'
                          IF(LYBMOD(NMOLY))THEN

!                             CHECK STATUS OF Y-SPECIES LINE TAIL FILE:
                              INQUIRE(FILE=DATDIR(1:LENDAT)//YFILE2,    &
     &                          EXIST=LEXIST,OPENED=LOPEN,NUMBER=FUNIT)
                              IF(LOPEN)THEN
                                  ITBY(NMOLY,2)=FUNIT
                                  REWIND(FUNIT)
                              ELSEIF(LEXIST)THEN
                                  ITBY(NMOLY,2)=NUNIT()
                                  CALL OPNFL(ITBY(NMOLY,2),0,           &
     &                              DATDIR(1:LENDAT)//YFILE2,           &
     &                              'OLD','FORMATTED','CD2C')
                                  WRITE(IPR,'(/(A))')'Successfully open'&
     &                              //'ed file',DATDIR(1:LENDAT)//YFILE2
                              ELSE
                                  WRITE(IPR,'(/(A))')'Error in'//       &
     &                              ' routine CD2C:  BANDMODEL FILE ',  &
     &                              DATDIR(1:LENDAT)//YFILE2,           &
     &                              'does not exist.  Create it or'//   &
     &                              ' skip species by adding a minus',  &
     &                              ' sign prefix to name in tape5.'
                                  STOP 'Band Model file does not exist!'
                              ENDIF

!                             READ MOLECULAR WEIGHT FROM FILE HEADER:
                              READ(ITBY(NMOLY,2),'(48X,F8.0)')          &
     &                          AMWTY(NMOLY)

!                             DOES A DEFAULT PROFILE EXIST FOR NMOLY?
                              DO JMOLY=1,16
                                  IF(YNAME(IMOLY).EQ.YNAM16(JMOLY))THEN

!                                     MATCH TO TRAC(JMOLY+5,*) PROFILE:
                                      MAPY_D(NMOLY)=JMOLY+5
                                      GOTO 10
                                  ENDIF
                              ENDDO

!                             NO MATCH TO A DEFAULT PROFILE:
                              MAPY_D(NMOLY)=0
                          ELSE

!                             DOES A DEFAULT PROFILE EXIST FOR NMOLY?
                              DO JMOLY=1,16
                                  IF(YNAME(IMOLY)(2:10).EQ.             &
     &                              YNAM16(JMOLY)(1:9))THEN

!                                     MATCH TO TRAC(JMOLY+5,*) PROFILE:
                                      MAPY_D(NMOLY)=JMOLY+5
                                      GOTO 10
                                  ENDIF
                              ENDDO

!                             NO MATCH TO A DEFAULT PROFILE:
                              MAPY_D(NMOLY)=0
                          ENDIF

!                         Compute DOPY (also done twice also):
   10                     CONTINUE
                          DOPY(NMOLY)=5.91861E-6/SQRT(AMWTY(NMOLY))
                      ENDIF
                  ENDDO
                  NMOLT=NMOLXT+NMOLY
              ENDIF
          ENDIF

!         REORDER SPECIES SO BAND MODEL SPECIES ARE FIRST:
          NMOLYS=0
          DO IMOLY=1,NMOLY
              IF(LYBMOD(IMOLY))THEN
                  NMOLYS=NMOLYS+1
                  IF(IMOLY.NE.NMOLYS)THEN

!                     SWAP MAPY_D, MAPY_F, ITBY, CNAMEY, DOPY AND AMWTY:
                      IDUM=MAPY_D(NMOLYS)
                      MAPY_D(NMOLYS)=MAPY_D(IMOLY)
                      MAPY_D(IMOLY)=IDUM
                      IDUM=MAPY_F(NMOLYS)
                      MAPY_F(NMOLYS)=MAPY_F(IMOLY)
                      MAPY_F(IMOLY)=IDUM
                      DUM=DOPY(NMOLYS)
                      DOPY(NMOLYS)=DOPY(IMOLY)
                      DOPY(IMOLY)=DUM
                      DUM=AMWTY(NMOLYS)
                      AMWTY(NMOLYS)=AMWTY(IMOLY)
                      AMWTY(IMOLY)=DUM
                      IDUM=ITBY(NMOLYS,1)
                      ITBY(NMOLYS,1)=ITBY(IMOLY,1)
                      ITBY(IMOLY,1)=IDUM
                      IDUM=ITBY(NMOLYS,2)
                      ITBY(NMOLYS,2)=ITBY(IMOLY,2)
                      ITBY(IMOLY,2)=IDUM
                      CDUM=CNAMEY(NMOLYS)
                      CNAMEY(NMOLYS)=CNAMEY(IMOLY)
                      CNAMEY(IMOLY)=CDUM
                  ENDIF
              ENDIF
          ENDDO
          IF(LAPLUS .AND. (IRD2.EQ.1.OR.IRD2.EQ.2))THEN
              WRITE(IPR,'(2A,/(10X,A))')' WARNING: ',                   &
     &          ' WHEN APLUS ="A+", IRD2 CANNOT BE 1 OR 2 (FOR',        &
     &          ' READING AEROSOL PROFILES WITH MODEL=0/7).',           &
     &          ' *****  APLUS OPTION WILL BE IGNORED  *****'
              LAPLUS=.FALSE.
              APLUS='  '
          ENDIF
          IF(IVSA.EQ.1)CALL RDNSM
      ENDIF
      MARIC1=0
      MARK=0
      JPRT=0

      IF(ICLD.GE.1 .AND. ICLD.LE.10)THEN

!         CLOUD/RAIN MODELS 1-10 ARE NOW SET UP IN ROUTINE CRDRIV, NOT
!         ROUTINE AERNSM; TEMPORARILY SET ICLD AND RAINRT TO ZERO.
          ICLDSV=ICLD
          RAINSV=RAINRT
          ICLD=0
          RAINRT=0.
          CALL AERNSM(JPRT,GNDALT,MARIC1,MARK,ICH,LMODEL)
          ICLD=ICLDSV
          RAINRT=RAINSV
          CALL CRDRIV
      ELSE
          CALL AERNSM(JPRT,GNDALT,MARIC1,MARK,ICH,LMODEL)
      ENDIF

!     IF NOVAM IS CALLED MERGE NOVAM LAYERS WITH OTHER.
      IF(LNOVAM)THEN
          NLNOV=2*(NNOV+1)
          CALL NOVMRG(ALTNOV,RHNOV,PNOV,TNOV,DENNOV,NLNOV)
      ENDIF

!     AERNSM IS USED AS BEFORE BUT NEW A+ SCHEME NEEDS EXTRA
!     HANDLING.  THE NEW SCHEME GOES INTO EFFECT AFTER AERNSM.
      IF(LAPLUS)CALL APRFNU(LMODEL,IHAZE,GNDALT)
      IF(ICLD.EQ.20)THEN

!         SET UP CIRRUS MODEL
          CALL CIRRUS(CTHIK,CALT,CPROB,CEXT)
          IF(.NOT.LJMASS)THEN
              WRITE(IPR,'(15X,A)')                                      &
     &          'CIRRUS ATTENUATION INCLUDED (N O A A CIRRUS)'
              WRITE(IPR,'((15X,2A,F10.5,A))')' CIRRUS THICKNESS ',      &
     &          'DEFAULTED TO MEAN VALUE OF',CTHIK,'KM','CIRRUS ',      &
     &          'BASE ALTITUDE DEFAULTED TO MEAN VALUE OF',CALT,'KM'
              WRITE(IPR,'(15X,A,F7.1,A)')                               &
     &          'PROBABILITY OF CLOUD OCCURRING IS',CPROB,' PERCENT'
          ENDIF
      ENDIF

!     SCALE OZONE AND/OR WATER PROFILE:
      IF(H2OSTR.NE.'          ' .OR. O3STR.NE.'          ')             &
     &  CALL SCLCOL(H2OSTR,O3STR,H2OCOL,O3COL)

      IF(LJMASS)THEN

!         SET DEFAULT AEROSOLS:
          IF(IHAZE.GT.0 .AND. JPRT.NE.0)THEN
              IF(ISEASN.EQ.0)ISEASN=1
              IF(IVULCN.LE.0)IVULCN=1
          ENDIF
          RETURN
      ENDIF

      IF(IHAZE.LE.0 .OR. JPRT.EQ.0)RETURN
      IF(ISEASN.EQ.0)ISEASN=1
      IF(IVULCN.LE.1)THEN
          IHMET=1
          IVULCN=1
      ELSE
          IHMET=2
      ENDIF
      IF(IVULCN.EQ.6 .OR. IVULCN.EQ.7)THEN
          IHVUL=11
      ELSEIF(IVULCN.EQ.8)THEN
          IHVUL=13
      ELSE
          IHVUL=IVULCN+10
      ENDIF
      WRITE(IPR,'(/A,3(/5X,A),F10.2,A,/(4(5X,A)))')' AEROSOL MODEL',    &
     &  'REGIME                      AEROSOL TYPE             '//       &
     &     'PROFILE                  SEASON',                           &
     &  '-----------------------     --------------------     '//       &
     &     '--------------------     --------------------',             &
     &  'BOUNDARY LAYER (0-2KM)      '//HHAZE(IHAZE),                   &
     &     VIS,' KM METEOROLOGICAL_RANGE_AT_SEA_LEVEL',                 &
     &  'TROPOSPHERE  (2-10KM)  ',HHAZE(6),HHAZE(6),HSEASN(ISEASN),     &
     &  'STRATOSPHERE (10-30KM) ',HHAZE(IHVUL),HVULCN(IVULCN),          &
     &  HSEASN(ISEASN),'UPPER ATMOS (30-100KM) ',HHAZE(16),HMET(IHMET)

!     RETURN TO DRIVER:
      RETURN
      END

      SUBROUTINE YFILE(IPR,YNAME,YFILE1,YFILE2,RESCHR)

!     ARGUMENTS:
!       RESCHR   BAND MODEL RESOLUTION CHARACTER STRING.
!       IPR      UNIT NUMBER OF PRIMARY OUTPUT FILE (tape6 or *.tp6).
      INTEGER IPR
      CHARACTER YNAME*(*),YFILE1*(28),YFILE2*(28),RESCHR*2

!     LOCAL VARIABLES:
!       LNYNAM   LENGTH OF YNAME CHARACTER STRING.
      INTEGER LNYNAM

!     FUNCTIONS:
!       LENSTR   RETURNS STRING LENGTH AFTER TRIMMING LEADING BLANKS.
      INTEGER LENSTR
      CALL LCTRIM(YNAME)
      LNYNAM=LENSTR(YNAME)
      YFILE2='                            '
      IF(LNYNAM.LE.0)THEN
          WRITE(IPR,'(/A)')' Error in YFILE:  Y-Species name is blank.'
          STOP ' Error in YFILE:  Y-Species name is blank.'
      ELSEIF(YNAME(LNYNAM:LNYNAM).EQ.'*')THEN

!         '*' => BAND MODEL MOLECULE WITH LINE TAIL & CENTER PARAMETERS:
          YFILE2(1:LNYNAM+14)='HITRANtrace/p1_'//YNAME(1:LNYNAM-1)
          YFILE2(13:14)=RESCHR
          YFILE1=YFILE2(1:LNYNAM+14)//'.tBM'
          YFILE2=YFILE2(1:LNYNAM+14)//'.cBM'
      ELSE

!         CROSS-SECTION MOLECULE:
          YFILE1='                            '
          YFILE1(1:3)=RESCHR//'_'
          YFILE1(4:LNYNAM+6)=YNAME(1:LNYNAM)//'.BM'
      ENDIF
      RETURN
      END
