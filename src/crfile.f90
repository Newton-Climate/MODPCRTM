      SUBROUTINE CRFILE(REFWAV,REFDEP,CWDCOL,CIPCOL)

!     THIS ROUTINE DEFINES USER-DEFINED CLOUD SPECTRAL DATA.

!     MODULES:
      USE CLDMOD
      USE CIRMOD

!     INPUT ARGUMENTS:
!       REFWAV   REFERENCE WAVELENGTH [MICRONS]
!       REFDEP   REFERENCE VERTICAL OPTICAL DEPTH AT REFWAV
!       CWDCOL   CLOUD WATER DROPLET VERTICAL COLUMN DENSITY [KM GM/M3]
!       CIPCOL   CLOUD ICE PARTICLE VERTICAL COLUMN DENSITY [KM GM/M3]
      REAL REFWAV,REFDEP,CWDCOL,CIPCOL

!     PARAMETERS:
!       DEG      NUMBER OF RADIANS IN ONE DEGREE.
      REAL DEG,TWOPI,FOURPI
      INCLUDE 'PARAMS.h'
      PARAMETER(DEG=57.29578,TWOPI=6.283185,FOURPI=12.56637)

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /CRD2E2/
!       CFILE    CLOUD DATA FILE NAME.
!       CLDTYP   WATER CLOUD TYPE OR NAME.
!       CIRTYP   CIRRUS CLOUD TYPE OR NAME.
      CHARACTER CFILE*(NAMLEN),CLDTYP*(NAMLEN),CIRTYP*(NAMLEN)
      COMMON/CRD2E2/CFILE,CLDTYP,CIRTYP

!     /DISRT/
!       UMU      MONOTONICALLY INCREASING LIST OF DISTINCT USER-PATH
!                COSINE POLAR ANGLES.
!       PHI      MONOTONICALLY INCREASING LIST OF DISTINCT RELATIVE
!                SOLAR AZIMUTH ANGLES [0 TO 180 DEG].
!       NSTR     NUMBER OF DISCRETE ORDINATE STREAMS.
!       NAZ      NUMBER OF DISORT AZIMUTH COMPONENTS.
!       N2GAUS   ORDER OF DOUBLE-GAUSS QUADRATURES.
!       NUMU     NUMBER OF DISTINCT USER LINE-OF-SIGHT POLAR ANGLES.
!       MAPUMU   MAPPING FROM LINE-OF-SIGHT INDEX TO UMU ARRAY ENTRY.
!       NPHI     NUMBER OF DISTINCT RELATIVE SOLAR AZIMUTH ANGLES.
!       MAPPHI   MAPPING FROM LINE-OF-SIGHT INDEX TO PHI ARRAY ENTRY.
!       DIS      LOGICAL FLAG, TRUE FOR DISORT MULTIPLE SCATTERING.
!       DISAZM   LOGICAL FLAG, TRUE FOR DISORT WITH AZIMUTH DEPENDENCE.
!       DISALB   LOGICAL FLAG, TRUE FOR DISORT SPHERICAL ALBEDO OPTION.
!       LDISCL   LOGICAL FLAG, TRUE FOR ISAACS SCALED TO DISORT.
      REAL UMU,PHI
      INTEGER NSTR,NAZ,N2GAUS,NUMU,MAPUMU,NPHI,MAPPHI
      LOGICAL DIS,DISAZM,DISALB,LDISCL
      COMMON/DISRT/UMU(MXUMU),PHI(MAXPHI),NSTR,NAZ,N2GAUS,NUMU,         &
     &  MAPUMU(MLOS),NPHI,MAPPHI(MLOS),DIS,DISAZM,DISALB,LDISCL
      SAVE /DISRT/

!     LOCAL VARIABLES:
!       LEXIST   FILE EXISTS FLAG.
!       PFWARN   PHASE FUNCTION NORMALIZATION PROBLEM WARNING FLAG.
!       LEGWRN   LEGENDRE EXPANSION PROBLEM WARNING FLAG.
!       LENFIL   FILE NAME LENGTH.
!       IUNIT    FILE UNIT NUMBER.
!       LENTYP   LENGTH OF CLOUD TYPE (TITLE) NAME.
!       IOS      IOSTAT VALUE.
!       IANGM1   PREVIOUS ANGULAR GRID INDEX.
!       INPSTR   INPUT CHARACTER STRING.
!       PFSUM    SCATTERING PHASE FUNCTION NORMALIZATION.
!       PFRAT    RATIO OF PHASE FUNCTION VALUES AT ANGULAR GRID POINTS.
!       RENORM   RENORMALIZATION FACTOR.
!       LEGSUM   SUM OF LEGENDRE COEFFICIENTS.
!       EXTWD    WATER CLOUD EXTINCTION COEFFICIENT AT
!                REFERENCE WAVELENGTH [KM-1 M3 / G].
!       EXTIP    CIRRUS CLOUD EXTINCTION COEFFICIENT AT
!                REFERENCE WAVELENGTH [KM-1 M3 / G].
      LOGICAL LEXIST,PFWARN,LEGWRN
      INTEGER LENFIL,IUNIT,LENTYP,IOS,IANGM1,ICWVM1
      CHARACTER INPSTR*80
      REAL PFRAT,RENORM,LEGSUM,EXTWD,EXTIP
      DOUBLE PRECISION PFSUM

!     FUNCTIONS:
!       LENSTR   RETURNS STRING LENGTH AFTER TRIMMING LEADING BLANKS.
!       NUNIT    RETURNS AN UNUSED FILE UNIT NUMBER.
      INTEGER LENSTR,NUNIT
      EXTERNAL LENSTR,NUNIT

!     DATA:
!       FRMT     CHARACTER STRING FORMAT.
      CHARACTER FRMT*7
      DATA FRMT/'(A0256)'/

!     OPEN CLOUD DATA FILE NAME:
      WRITE(FRMT(3:6),'(I4.4)')NAMLEN
      READ(IRD,FRMT)CFILE
      LENFIL=LENSTR(CFILE)
      IF(LENFIL.LE.0)THEN
          WRITE(IPR,'(/A)')                                             &
     &      ' Error in CRFILE:  Water cloud data file name is blank.'
          STOP 'Error in CRFILE:  Water cloud data file name is blank.'
      ENDIF
      INQUIRE(FILE=CFILE(1:LENFIL),EXIST=LEXIST)
      IF(.NOT.LEXIST)THEN
          WRITE(IPR,'(/3A)')' Error in CRFILE:  Water cloud data file ',&
     &      CFILE(1:LENFIL),' does not exist.'
          STOP 'Error in CRFILE:  Water cloud data file does not exist.'
      ENDIF
      IUNIT=NUNIT()
      CALL OPNFL(IUNIT,0,CFILE(1:LENFIL),'OLD','FORMATTED','CRFILE')

!     READ WATER CLOUD TYPE AND FIND IT IN DATA FILE:
      READ(IRD,'(A80)')CLDTYP
      LENTYP=LENSTR(CLDTYP)
      IF(LENTYP.LE.0)THEN
          WRITE(IPR,'(/A)')                                             &
     &      ' Error in CRFILE:  Water cloud type is blank.'
          STOP 'Error in CRFILE:  Water cloud type is blank.'
      ENDIF
      CALL UPCASE(CLDTYP)
   10 CONTINUE
          READ(IUNIT,'(A80)',IOSTAT=IOS)INPSTR
          IF(IOS.LT.0)THEN
              WRITE(IPR,'(/3A,/(18X,2A))')' Error in CRFILE: ',         &
     &          ' End-of-File detected in',CFILE(1:LENFIL),             &
     &          ' and cloud type ',CLDTYP,' was not found.'
              STOP 'Error in CRFILE:  Water cloud type not found.'
          ENDIF
          IF(LENSTR(INPSTR).LE.0)GOTO 10
          CALL UPCASE(INPSTR)
      IF(INPSTR.NE.CLDTYP)GOTO 10

!     READ NUMBER OF ANGULAR, LEGENDRE EXPANSION & SPECTRAL GRID POINTS:
      READ(IUNIT,*,IOSTAT=IOS)NCLDAN,NCLDLG,NCLDWV
      IF(IOS.NE.0 .OR. NCLDAN.LT.2 .OR.                                 &
     &  NCLDLG.LT.NSTR .OR. NCLDWV.LT.2)THEN
          WRITE(IPR,'(/2A,/(18X,A))')' Error in CRFILE:  Problem',      &
     &      ' occurred while reading NUMBER OF ANGULAR LEGENDRE,',      &
     &      ' EXPANSION & SPECTRAL GRID POINTS for cloud type ',        &
     &      CLDTYP(1:LENTYP)
          STOP 'Error in CRFILE:  I/O Status reading cloud data file'
      ENDIF

!     ALLOCATE ARRAYS:
      ALLOCATE(CLDANG(NCLDAN),CLDCOS(NCLDAN),CLDWAV(0:NCLDWV),          &
     &         CLDEXT(NCLDWV),CLDABS(NCLDWV),CLDPF(NCLDAN,NCLDWV),      &
     &         CLDLEG(0:NCLDLG,1:NCLDWV),CLDPFV(NCLDAN))

!     READ ANGULAR GRID (INCREASING FROM 0 TO 180):
      READ(IUNIT,'(A80)')INPSTR
      READ(IUNIT,*,IOSTAT=IOS)CLDANG
      IF(IOS.NE.0)THEN
          WRITE(IPR,'(/2A,/19X,A)')                                     &
     &      ' Error in CRFILE:  Problem occurred while reading',        &
     &      ' ANGULAR GRID for cloud type,',CLDTYP(1:LENTYP)
          STOP 'Error in CRFILE:  Unable to read ANGULAR GRID.'
      ENDIF
      IF(CLDANG(1).NE.0.)THEN
          WRITE(IPR,'(/2A,/19X,A)')                                     &
     &      ' Error in CRFILE:  First angle of ANGULAR GRID is',        &
     &      ' not 0. for cloud type,',CLDTYP(1:LENTYP)
          STOP 'Error in CRFILE:  First angle is not 0.'
      ENDIF
      CLDCOS(1)=1.
      IANGM1=1
      DO ICLDAN=2,NCLDAN
          IF(CLDANG(ICLDAN).LE.CLDANG(IANGM1))THEN
              WRITE(IPR,'(/2A,/19X,A)')                                 &
     &          ' Error in CRFILE:  ANGULAR GRID is not monotonically', &
     &          ' increasing for cloud type,',CLDTYP(1:LENTYP)
              STOP 'Error in CRFILE:  Non-monotonic ANGULAR GRID.'
          ENDIF
          CLDCOS(ICLDAN)=COS(CLDANG(ICLDAN)/DEG)
          IANGM1=ICLDAN
      ENDDO
      IF(ABS(CLDANG(NCLDAN)-180.).GT..001)THEN
          WRITE(IPR,'(/2A,/19X,A)')                                     &
     &      ' Error in CRFILE:  Last angle of ANGULAR GRID is',         &
     &      ' not 180. for cloud type,',CLDTYP(1:LENTYP)
          STOP 'Error in CRFILE:  Last angle is not 180.'
      ENDIF
      CLDCOS(NCLDAN)=-1.

!     LOOP OVER WAVELENGTH:
      CLDWAV(0)=0.
      PFWARN=.TRUE.
      LEGWRN=.TRUE.
      DO ICLDWV=1,NCLDWV

!         READ WAVELENGTH (IN INCREASING MICRONS):
          READ(IUNIT,*,IOSTAT=IOS)                                      &
     &      CLDWAV(ICLDWV),CLDEXT(ICLDWV),CLDABS(ICLDWV)
          IF(IOS.NE.0)THEN
              WRITE(IPR,'(/3A,/19X,A)')' Error in CRFILE:  Problem',    &
     &          ' occurred while reading WAVELENGTH, EXTINCTION COEF',  &
     &          ' or ABSORPTION COEF for cloud type,',CLDTYP(1:LENTYP)
              STOP 'Error in CRFILE:  Unable to read WAVELENGTH GRID.'
          ELSEIF(CLDWAV(ICLDWV).LT.CLDWAV(ICLDWV-1))THEN
              WRITE(IPR,'(/3A,/19X,A)')' Error in CRFILE: ',            &
     &          ' WAVELENGTH GRID is not monotonically',                &
     &          ' increasing for cloud type,',CLDTYP(1:LENTYP)
              STOP 'Error in CRFILE:  Non-monotonic WAVELENGTH GRID.'
          ELSEIF(CLDEXT(ICLDWV).LT.CLDABS(ICLDWV)                       &
     &      .OR. CLDABS(ICLDWV).LT.0.)THEN
              WRITE(IPR,'(/3A,/19X,A)')' Error in CRFILE: ',            &
     &          ' A scattering albedo not between 0 and 1',             &
     &          ' detected for cloud type,',CLDTYP(1:LENTYP)
              STOP 'Error in CRFILE:  Non-monotonic WAVELENGTH GRID.'
          ENDIF

!         READ SCATTERING PHASE FUNCTION VALUE:
          READ(IUNIT,'(A80)')INPSTR
          READ(IUNIT,*,IOSTAT=IOS)(CLDPF(ICLDAN,ICLDWV),ICLDAN=1,NCLDAN)
          IF(IOS.NE.0)THEN
              WRITE(IPR,'(/2A,/19X,A)')                                 &
     &          ' Error in CRFILE:  Problem occurred while reading',    &
     &          ' PHASE FUNCTION DATA for cloud type,',CLDTYP(1:LENTYP)
              STOP 'Error in CRFILE:  Unable to read PHASE FUNCTION.'
          ENDIF

!         CHECK PHASE FUNCTION VALUES AND NORMALIZATION:
          IF(CLDPF(1,ICLDWV).LT.0.)THEN
              WRITE(IPR,'(/2A,/19X,A)')                                 &
     &          ' Error in CRFILE:  A negative PHASE FUNCTION value',   &
     &          ' was read in for cloud type,',CLDTYP(1:LENTYP)
              STOP 'Error in CRFILE:  Negative PHASE FUNCTION value.'
          ENDIF
          PFSUM=0.D0
          IANGM1=1

!         ZENITH ANGLE INTEGRATION (ASSUME EXPONENTIAL VARIATION):
          DO ICLDAN=2,NCLDAN
              IF(CLDPF(ICLDAN,ICLDWV).LT.0.)THEN
                  WRITE(IPR,'(/3A,/19X,A)')' Error in CRFILE: ',        &
     &              ' A negative PHASE FUNCTION value',                 &
     &              ' was read in for cloud type,',CLDTYP(1:LENTYP)
                  STOP 'Error in CRFILE:  Negative PHASE FUNCTION value'
              ENDIF
              PFRAT=CLDPF(IANGM1,ICLDWV)/CLDPF(ICLDAN,ICLDWV)
              IF(ABS(PFRAT-1).GT..001)THEN
                  PFSUM=PFSUM+DBLE((CLDCOS(IANGM1)-CLDCOS(ICLDAN))/LOG( &
     &              PFRAT)*(CLDPF(IANGM1,ICLDWV)-CLDPF(ICLDAN,ICLDWV)))
              ELSE
                  PFSUM=PFSUM+DBLE(2*(CLDCOS(IANGM1)-CLDCOS(ICLDAN))    &
     &              *CLDPF(IANGM1,ICLDWV)/(3-PFRAT))
              ENDIF
              IANGM1=ICLDAN

!         END ANGULAR GRID
          ENDDO

!         RENORMALIZE IF NOT CORRECT:
          RENORM=TWOPI*REAL(PFSUM)
          IF(ABS(RENORM-1).GT..00001)THEN
              IF(PFWARN)THEN
                  WRITE(IPR,'(/A,F10.5,2A,2(/23X,A),F10.5,A,/23X,2A)')  &
     &              ' Warning from CRFILE:  At',CLDWAV(ICLDWV),         &
     &              ' Microns, the SCATTERING PHASE FUNCTION for',      &
     &              ' cloud type',CLDTYP(1:LENTYP),'was normalized to', &
     &              RENORM,' instead of 1.  It is being renormalized.', &
     &              'This warning will not be repeated if the problem', &
     &              ' occurs for other wavelengths.'
                  PFWARN=.FALSE.
              ENDIF
              DO ICLDAN=1,NCLDAN
                  CLDPF(ICLDAN,ICLDWV)=CLDPF(ICLDAN,ICLDWV)/RENORM
              ENDDO
          ENDIF

!         READ PHASE FUNCTION LEGENDRE EXPANSION COEFFICIENTS:
          READ(IUNIT,'(A80)')INPSTR
          READ(IUNIT,*,IOSTAT=IOS)                                      &
     &      (CLDLEG(ICLDLG,ICLDWV),ICLDLG=0,NCLDLG)
          IF(IOS.NE.0)THEN
              WRITE(IPR,'(/3A,/19X,A)')' Error in CRFILE:  Problem',    &
     &          ' occurred while reading PHASE FUNCTION LEGENDRE',      &
     &          ' EXPANSION COEFS for cloud type,',CLDTYP(1:LENTYP)
              STOP 'Error in CRFILE: Unable to read LEGENDRE EXPANSION.'
          ENDIF

!         RENORMALIZE AND CHECK FORWARD SCATTER VALUE:
          RENORM=CLDLEG(0,ICLDWV)
          CLDLEG(0,ICLDWV)=1.
          LEGSUM=1.
          DO ICLDLG=1,NCLDLG
              CLDLEG(ICLDLG,ICLDWV)=CLDLEG(ICLDLG,ICLDWV)/RENORM
              LEGSUM=LEGSUM+(2*ICLDLG+1)*CLDLEG(ICLDLG,ICLDWV)
          ENDDO
          IF(LEGWRN .AND.                                               &
     &      ABS(LEGSUM/CLDPF(1,ICLDWV)-FOURPI).GT..001)THEN
              WRITE(IPR,'(/A,F10.5,2A,2(/23X,A,F12.3,A),2(/23X,A),A)')  &
     &          ' Warning from CRFILE:  At',CLDWAV(ICLDWV),' Microns,', &
     &          ' the forward scatter value of the unit normalized',    &
     &          'phase function (=',CLDPF(1,ICLDWV),                    &
     &          ' sr-1) does not equal the sum of the unit',            &
     &          'normalized LEGENDRE EXPANSION COEFFICIENTS (=',        &
     &          LEGSUM/FOURPI,' sr-1) for cloud type',                  &
     &          CLDTYP(1:LENTYP),'This warning will not be repeated',   &
     &          ' if the problem occurs for other wavelengths.'
              LEGWRN=.FALSE.
          ENDIF

!     END WAVELENGTH GRID:
      ENDDO

!     ANNOUNCE SUCCESS:
      WRITE(IPR,'(/(A))')                                               &
     &  ' Successfully read data for cloud type',CLDTYP(1:LENTYP)

!     READ CIRRUS CLOUD TYPE AND FIND IT IN DATA FILE:
      READ(IRD,'(A80)')CIRTYP
      LENTYP=LENSTR(CIRTYP)
      IF(LENTYP.LE.0)THEN
          WRITE(IPR,'(/A)')                                             &
     &      ' Error in CRFILE:  Cirrus cloud type is blank.'
          STOP 'Error in CRFILE:  Cirrus cloud type is blank.'
      ENDIF
      CALL UPCASE(CIRTYP)
      REWIND(IUNIT)
   20 CONTINUE
          READ(IUNIT,'(A80)',IOSTAT=IOS)INPSTR
          IF(IOS.LT.0)THEN
              WRITE(IPR,'(/3A,/(18X,2A))')' Error in CRFILE: ',         &
     &          ' End-of-File detected in',CFILE(1:LENFIL),             &
     &          ' and cloud type ',CIRTYP,' was not found.'
              STOP 'Error in CRFILE:  Cirrus cloud type not found.'
          ENDIF
          IF(LENSTR(INPSTR).LE.0)GOTO 20
          CALL UPCASE(INPSTR)
      IF(INPSTR.NE.CIRTYP)GOTO 20

!     READ NUMBER OF ANGULAR, LEGENDRE EXPANSION & SPECTRAL GRID POINTS:
      READ(IUNIT,*,IOSTAT=IOS)NCIRAN,NCIRLG,NCIRWV
      IF(IOS.NE.0 .OR. NCIRAN.LT.2 .OR.                                 &
     &  NCIRLG.LT.NSTR .OR. NCIRWV.LT.2)THEN
          WRITE(IPR,'(/2A,/(18X,A))')' Error in CRFILE:  Problem',      &
     &      ' occurred while reading NUMBER OF ANGULAR LEGENDRE,',      &
     &      ' EXPANSION & SPECTRAL GRID POINTS for cloud type ',        &
     &      CIRTYP(1:LENTYP)
          STOP 'Error in CRFILE:  I/O Status reading cloud data file'
      ENDIF

!     ALLOCATE ARRAYS:
      ALLOCATE(CIRANG(NCIRAN),CIRCOS(NCIRAN),CIRWAV(0:NCIRWV),          &
     &         CIREXT(NCIRWV),CIRABS(NCIRWV),CIRPF(NCIRAN,NCIRWV),      &
     &         CIRLEG(0:NCIRLG,1:NCIRWV),CIRPFV(NCIRAN))

!     READ ANGULAR GRID (INCREASING FROM 0 TO 180):
      READ(IUNIT,'(A80)')INPSTR
      READ(IUNIT,*,IOSTAT=IOS)CIRANG
      IF(IOS.NE.0)THEN
          WRITE(IPR,'(/2A,/19X,A)')                                     &
     &      ' Error in CRFILE:  Problem occurred while reading',        &
     &      ' ANGULAR GRID for cloud type,',CIRTYP(1:LENTYP)
          STOP 'Error in CRFILE:  Unable to read ANGULAR GRID.'
      ENDIF
      IF(CIRANG(1).NE.0.)THEN
          WRITE(IPR,'(/2A,/19X,A)')                                     &
     &      ' Error in CRFILE:  First angle of ANGULAR GRID is',        &
     &      ' not 0. for cloud type,',CIRTYP(1:LENTYP)
          STOP 'Error in CRFILE:  First angle is not 0.'
      ENDIF
      CIRCOS(1)=1.
      IANGM1=1
      DO ICIRAN=2,NCIRAN
          IF(CIRANG(ICIRAN).LE.CIRANG(IANGM1))THEN
              WRITE(IPR,'(/2A,/19X,A)')                                 &
     &          ' Error in CRFILE:  ANGULAR GRID is not monotonically', &
     &          ' increasing for cloud type,',CIRTYP(1:LENTYP)
              STOP 'Error in CRFILE:  Non-monotonic ANGULAR GRID.'
          ENDIF
          CIRCOS(ICIRAN)=COS(CIRANG(ICIRAN)/DEG)
          IANGM1=ICIRAN
      ENDDO
      IF(ABS(CIRANG(NCIRAN)-180.).GT..001)THEN
          WRITE(IPR,'(/2A,/19X,A)')                                     &
     &      ' Error in CRFILE:  Last angle of ANGULAR GRID is',         &
     &      ' not 180. for cloud type,',CIRTYP(1:LENTYP)
          STOP 'Error in CRFILE:  Last angle is not 180.'
      ENDIF
      CIRCOS(NCIRAN)=-1.

!     LOOP OVER WAVELENGTH:
      CIRWAV(0)=0.
      PFWARN=.TRUE.
      LEGWRN=.TRUE.
      DO ICIRWV=1,NCIRWV

!         READ WAVELENGTH (IN INCREASING MICRONS):
          READ(IUNIT,*,IOSTAT=IOS)                                      &
     &      CIRWAV(ICIRWV),CIREXT(ICIRWV),CIRABS(ICIRWV)
          IF(IOS.NE.0)THEN
              WRITE(IPR,'(/3A,/19X,A)')' Error in CRFILE:  Problem',    &
     &          ' occurred while reading WAVELENGTH, EXTINCTION COEF',  &
     &          ' or ABSORPTION COEF for cloud type,',CIRTYP(1:LENTYP)
              STOP 'Error in CRFILE:  Unable to read WAVELENGTH GRID.'
          ELSEIF(CIRWAV(ICIRWV).LT.CIRWAV(ICIRWV-1))THEN
              WRITE(IPR,'(/3A,/19X,A)')' Error in CRFILE: ',            &
     &          ' WAVELENGTH GRID is not monotonically',                &
     &          ' increasing for cloud type,',CIRTYP(1:LENTYP)
              STOP 'Error in CRFILE:  Non-monotonic WAVELENGTH GRID.'
          ELSEIF(CIREXT(ICIRWV).LT.CIRABS(ICIRWV)                       &
     &      .OR. CIRABS(ICIRWV).LT.0.)THEN
              WRITE(IPR,'(/3A,/19X,A)')' Error in CRFILE: ',            &
     &          ' A scattering albedo not between 0 and 1',             &
     &          ' detected for cloud type,',CIRTYP(1:LENTYP)
              STOP 'Error in CRFILE:  Non-monotonic WAVELENGTH GRID.'
          ENDIF

!         READ SCATTERING PHASE FUNCTION VALUE:
          READ(IUNIT,'(A80)')INPSTR
          READ(IUNIT,*,IOSTAT=IOS)(CIRPF(ICIRAN,ICIRWV),ICIRAN=1,NCIRAN)
          IF(IOS.NE.0)THEN
              WRITE(IPR,'(/2A,/19X,A)')                                 &
     &          ' Error in CRFILE:  Problem occurred while reading',    &
     &          ' PHASE FUNCTION DATA for cloud type,',CIRTYP(1:LENTYP)
              STOP 'Error in CRFILE:  Unable to read PHASE FUNCTION.'
          ENDIF

!         CHECK PHASE FUNCTION VALUES AND NORMALIZATION:
          IF(CIRPF(1,ICIRWV).LT.0.)THEN
              WRITE(IPR,'(/2A,/19X,A)')                                 &
     &          ' Error in CRFILE:  A negative PHASE FUNCTION value',   &
     &          ' was read in for cloud type,',CIRTYP(1:LENTYP)
              STOP 'Error in CRFILE:  Negative PHASE FUNCTION value.'
          ENDIF
          PFSUM=0.D0
          IANGM1=1

!         ZENITH ANGLE INTEGRATION (ASSUME EXPONENTIAL VARIATION):
          DO ICIRAN=2,NCIRAN
              IF(CIRPF(ICIRAN,ICIRWV).LT.0.)THEN
                  WRITE(IPR,'(/3A,/19X,A)')' Error in CRFILE: ',        &
     &              ' A negative PHASE FUNCTION value',                 &
     &              ' was read in for cloud type,',CIRTYP(1:LENTYP)
                  STOP 'Error in CRFILE:  Negative PHASE FUNCTION value'
              ENDIF
              PFRAT=CIRPF(IANGM1,ICIRWV)/CIRPF(ICIRAN,ICIRWV)
              IF(ABS(PFRAT-1).GT..001)THEN
                  PFSUM=PFSUM+DBLE((CIRCOS(IANGM1)-CIRCOS(ICIRAN))/LOG  &
     &              (PFRAT)*(CIRPF(IANGM1,ICIRWV)-CIRPF(ICIRAN,ICIRWV)))
              ELSE
                  PFSUM=PFSUM+DBLE(2*(CIRCOS(IANGM1)-CIRCOS(ICIRAN))    &
     &              *CIRPF(IANGM1,ICIRWV)/(3-PFRAT))
              ENDIF
              IANGM1=ICIRAN

!         END ANGULAR GRID
          ENDDO

!         RENORMALIZE IF NOT CORRECT:
          RENORM=TWOPI*REAL(PFSUM)
          IF(ABS(RENORM-1).GT..00001)THEN
              IF(PFWARN)THEN
                  WRITE(IPR,'(/A,F10.5,2A,2(/23X,A),F10.5,A,/23X,2A)')  &
     &              ' Warning from CRFILE:  At',CIRWAV(ICIRWV),         &
     &              ' Microns, the SCATTERING PHASE FUNCTION for',      &
     &              ' cloud type',CIRTYP(1:LENTYP),'was normalized to', &
     &              RENORM,' instead of 1.  It is being renormalized.', &
     &              'This warning will not be repeated if the problem', &
     &              ' occurs for other wavelengths.'
                  PFWARN=.FALSE.
              ENDIF
              DO ICIRAN=1,NCIRAN
                  CIRPF(ICIRAN,ICIRWV)=CIRPF(ICIRAN,ICIRWV)/RENORM
              ENDDO
          ENDIF

!         READ PHASE FUNCTION LEGENDRE EXPANSION COEFFICIENTS:
          READ(IUNIT,'(A80)')INPSTR
          READ(IUNIT,*,IOSTAT=IOS)                                      &
     &      (CIRLEG(ICIRLG,ICIRWV),ICIRLG=0,NCIRLG)
          IF(IOS.NE.0)THEN
              WRITE(IPR,'(/3A,/19X,A)')' Error in CRFILE:  Problem',    &
     &          ' occurred while reading PHASE FUNCTION LEGENDRE',      &
     &          ' EXPANSION COEFS for cloud type,',CIRTYP(1:LENTYP)
              STOP 'Error in CRFILE: Unable to read LEGENDRE EXPANSION.'
          ENDIF

!         RENORMALIZE AND CHECK FORWARD SCATTER VALUE:
          LEGSUM=0.
          RENORM=CIRLEG(0,ICIRWV)
          DO ICIRLG=0,NCIRLG
              CIRLEG(ICIRLG,ICIRWV)=CIRLEG(ICIRLG,ICIRWV)/RENORM
              LEGSUM=LEGSUM+(2*ICIRLG+1)*CIRLEG(ICIRLG,ICIRWV)
          ENDDO
          IF(LEGWRN .AND.                                               &
     &      ABS(LEGSUM/CIRPF(1,ICIRWV)-FOURPI).GT..001)THEN
              WRITE(IPR,'(/A,F10.5,2A,2(/23X,A,F12.3,A),2(/23X,A),A)')  &
     &          ' Warning from CRFILE:  At',CIRWAV(ICIRWV),' Microns,', &
     &          ' the forward scatter value of the unit normalized',    &
     &          'phase function (=',CIRPF(1,ICIRWV),                    &
     &          ' sr-1) does not equal the sum of the unit',            &
     &          'normalized LEGENDRE EXPANSION COEFFICIENTS (=',        &
     &          LEGSUM/FOURPI,' sr-1) for cloud type',                  &
     &          CIRTYP(1:LENTYP),'This warning will not be repeated',   &
     &          ' if the problem occurs for other wavelengths.'
              LEGWRN=.FALSE.
          ENDIF

!     END WAVELENGTH GRID:
      ENDDO

!     ANNOUNCE SUCCESS:
      WRITE(IPR,'(/(A))')                                               &
     &  ' Successfully read data for cloud type',CIRTYP(1:LENTYP)
      CLOSE(IUNIT)

!     RENORMALIZE SPECTRAL DATA IF REFERENCE OPTICAL DEPTH WAS INPUT
      IF(REFDEP.GT.0.)THEN

!         DETERMINE BRACKETING WAVELENGTHS:
          IF(REFWAV.LT.CLDWAV(1) .OR. REFWAV.GT.CLDWAV(NCLDWV))THEN
              WRITE(IPR,'(/A,F10.6,A,/14X,2A,/14X,A,2(F10.6,A))')       &
     &          ' Error in CRFILE:  The wavelength (',REFWAV,           &
     &          ' microns) used to define the cloud',                   &
     &          ' vertical optical depth is outside the range',         &
     &          ' of the user-defined',' water cloud spectral data (',  &
     &          CLDWAV(1),' to',CLDWAV(NCLDWV),' microns).'
              STOP 'Input spectral cloud depth outside spectral range.'
          ELSEIF(REFWAV.LT.CIRWAV(1) .OR. REFWAV.GT.CIRWAV(NCIRWV))THEN
              WRITE(IPR,'(/A,F10.6,A,/14X,2A,/14X,A,2(F10.6,A))')       &
     &          ' Error in CRFILE:  The wavelength (',REFWAV,           &
     &          ' microns) used to define the cloud',                   &
     &          ' vertical optical depth is outside the range',         &
     &          ' of the user-defined',' cirrus cloud spectral data (', &
     &          CIRWAV(1),' to',CIRWAV(NCIRWV),' microns).'
              STOP 'Input spectral cloud depth outside spectral range.'
          ENDIF

!         DETERMINE REFERENCE WATER CLOUD EXTINCTION COEFFICIENT:
          ICWVM1=1
          DO ICLDWV=2,NCLDWV-1
              IF(REFWAV.LE.CLDWAV(ICLDWV))GOTO 30
              ICWVM1=ICLDWV
          ENDDO
          ICLDWV=NCLDWV
   30     CONTINUE
          EXTWD=CLDEXT(ICWVM1)+(CLDEXT(ICLDWV)-CLDEXT(ICWVM1))          &
     &      *(REFWAV-CLDWAV(ICWVM1))/(CLDWAV(ICLDWV)-CLDWAV(ICWVM1))

!         DETERMINE REFERENCE CIRRUS CLOUD EXTINCTION COEFFICIENT:
          ICWVM1=1
          DO ICIRWV=2,NCIRWV-1
              IF(REFWAV.LE.CIRWAV(ICIRWV))GOTO 40
              ICWVM1=ICIRWV
          ENDDO
          ICIRWV=NCIRWV
   40     CONTINUE
          EXTIP=CIREXT(ICWVM1)+(CIREXT(ICIRWV)-CIREXT(ICWVM1))          &
     &      *(REFWAV-CIRWAV(ICWVM1))/(CIRWAV(ICIRWV)-CIRWAV(ICWVM1))

!         DETERMINE RATIO OF INPUT TO CURRENT CLOUD DEPTH
          RENORM=REFDEP/(EXTWD*CWDCOL+EXTIP*CIPCOL)

!         SCALE THE CLOUD PARTICLE SPECTRAL DATA
          DO ICLDWV=1,NCLDWV
              CLDEXT(ICLDWV)=RENORM*CLDEXT(ICLDWV)
              CLDABS(ICLDWV)=RENORM*CLDABS(ICLDWV)
          ENDDO
          DO ICIRWV=1,NCIRWV
              CIREXT(ICIRWV)=RENORM*CIREXT(ICIRWV)
              CIRABS(ICIRWV)=RENORM*CIRABS(ICIRWV)
          ENDDO
      ENDIF

!     RETURN TO CRUSPC IF SPECTRAL DATA IS NOT TO BE OUTPUT.
      IF(NPR.GE.0 .OR. LJMASS)RETURN

!     WATER CLOUD SPECTRAL DATA HEADER
      WRITE(IPR,'(/A,/(2A))')' WATER CLOUD SPECTRAL DATA',              &
     &  ' IWAV   WAVLEN       FREQ',                                    &
     &  '    EXT COEF    ABS COEF    SCT COEF     ASYM  SCT ALB',       &
     &  '      (MICRON)     (CM-1)',                                    &
     &  '  (           KM-1 M3/GM           )'

!     WATER CLOUD SPECTRAL DATA
      WRITE(IPR,'((I4,F10.4,F11.3,3(1X,F11.5),2F9.5))')                 &
     &  (ICLDWV,CLDWAV(ICLDWV),10000./CLDWAV(ICLDWV),CLDEXT(ICLDWV),    &
     &  CLDABS(ICLDWV),CLDEXT(ICLDWV)-CLDABS(ICLDWV),CLDLEG(1,ICLDWV),  &
     &  1.-CLDABS(ICLDWV)/CLDEXT(ICLDWV),ICLDWV=1,NCLDWV)
      WRITE(IPR,'(/A,/)')' END OF WATER CLOUD PARTICLE SPECTRAL DATA'

!     CIRRUS CLOUD SPECTRAL DATA HEADER
      WRITE(IPR,'(/A,/(2A))')' CIRRUS CLOUD SPECTRAL DATA',             &
     &  ' IWAV   WAVLEN       FREQ',                                    &
     &  '    EXT COEF    ABS COEF    SCT COEF     ASYM  SCT ALB',       &
     &  '      (MICRON)     (CM-1)',                                    &
     &  '  (           KM-1 M3/GM           )'

!     CIRRUS CLOUD SPECTRAL DATA
      WRITE(IPR,'((I4,F10.4,F11.3,3(1X,F11.5),2F9.5))')                 &
     &  (ICIRWV,CIRWAV(ICIRWV),10000./CIRWAV(ICIRWV),CIREXT(ICIRWV),    &
     &  CIRABS(ICIRWV),CIREXT(ICIRWV)-CIRABS(ICIRWV),CIRLEG(1,ICIRWV),  &
     &  1.-CIRABS(ICIRWV)/CIREXT(ICIRWV),ICIRWV=1,NCIRWV)
      WRITE(IPR,'(/A,/)')' END OF CIRRUS CLOUD PARTICLE SPECTRAL DATA'

!     RETURN TO ROUTINE CRSPEC
      RETURN
      END
