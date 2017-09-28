      SUBROUTINE RDCORK(CKNAME,NSPEED,KNTRVL)

!     ROUTINE TO INITIALIZE CORRELATED-K DISTRIBUTION DATA.

!     DECLARE ROUTINE ARGUMENT INPUTS.
!       CKNAME   NAME OF CORRELATED-K DISTRIBUTIONS FILE.
!       NSPEED   COMPUTATIONAL SPEED FLAG.
      CHARACTER CKNAME*(*)
      INTEGER NSPEED

!     DECLARE ROUTINE ARGUMENT OUTPUTS.
!       KNTRVL   NUMBER OF K-INTERVAL TO BE USED.
      INTEGER KNTRVL

!     PARAMETERS:
!       MXKSUB   DIMENSION OF K-DISTRIBUTION SUB-INTERVAL ARRAY.
!       MXGAML   DIMENSION OF LORENTZ HALF-WIDTH ARRAY.
!       MXGAMD   DIMENSION OF DOPPLER HALF-WIDTH ARRAY.
!       MXNUML   DIMENSION OF EFFECTIVE NUMBER OF LINES ARRAY.
      INCLUDE 'PARAMS.h'

!     COMMONS:

!     /CORKDT/
!       WTKSUB   SPECTRAL BIN SUB-INTERVAL FRACTIONAL WIDTHS.
!       WTKSAV   SAVED SPECTRAL BIN SUB-INTERVAL FRACTIONAL WIDTHS.
!       DEPLAY   INCREMENTAL OPTICAL DEPTHS.
!       TRNLAY   INCREMENTAL TRANSMITTANCES.
!       TRNCUM   CUMULATIVE TRANSMITTANCES.
!       K2TAIL   POINTER FROM K BIN TO LINE TAIL SUB-BIN
!                (=0 IF MULTIPLE LINE TAIL SUB-BINS CONTRIBUTE).
!       CONTWT   WEIGHTS FOR PARTITIONING LINE TAILS INTO K'S
!                (ONLY USED IF K2TAIL IS 0).
      REAL WTKSUB(MXKSUB),WTKSAV(NTLSUB),DEPLAY(MXKSUB),                &
     &  TRNLAY(MXKSUB),TRNCUM(MXKSUB),CONTWT(NTLSUB,MXKSUB)
      INTEGER K2TAIL(MXKSUB)
      COMMON/CORKDT/K2TAIL,WTKSUB,WTKSAV,DEPLAY,TRNLAY,TRNCUM,CONTWT
      SAVE /CORKDT/

!     COMMON /CORKTB/
!       GAMLIN   LORENTZ HALF-WIDTH BOUNDARY VALUES [CM-1].
!       GAMDIN   DOPPLER HALF-WIDTH BOUNDARY VALUES [CM-1].
!       RLININ   EFFECTIVE NUMBER OF LINES BOUNDARY VALUES.
!       VAL      ABSORPTION COEFFICIENT TABLE [ATM-1 CM-1].
      INTEGER NGAML,NGAMD,NNUML
      REAL GAMLIN,GAMDIN,RLININ,VAL
      COMMON/CORKTB/NGAML,NGAMD,NNUML,GAMLIN(MXGAML),GAMDIN(MXGAMD),    &
     &  RLININ(MXNUML),VAL(1:MXGAML,1:MXGAMD,1:MXNUML,0:MXKSUB)
      SAVE /CORKTB/
      INCLUDE 'IFIL.h'

!     BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     LOCAL VARIABLES:
!       LENCK    LENGTH OF CORRELATED-K DATA FILE NAME.
      INTEGER IFILE,IKSUB,IKSBP1,IGAML,IGAMD,INUML,NKSUB,INTRVL,        &
     &  ITLSUB,ITLSB0,NDIVM1,ILAST,INEXT,LENCK
      LOGICAL LTEST
      REAL WGT,TAILWT(NTLSUB)
!old  REAL SUMWTK,CUMWTK(MXKSUB)

!     FUNCTIONS:
!       NUNIT    RETURNS AN UNUSED FILE UNIT NUMBER.
!       LENSTR   RETURNS STRING LENGTH AFTER TRIMMING LEADING BLANKS.
      INTEGER NUNIT,LENSTR

!     DATA:
!       LSKIP    LIST OF K'S TO BE SKIPPED IF NSPEED IS 1 OR 2.
      LOGICAL LSKIP(MXKSUB,2)
      DATA (LSKIP(IKSUB,1),IKSUB=1,MXKSUB)     /.FALSE.,.FALSE.,.FALSE.,&
     &  .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.TRUE. ,.FALSE.,.TRUE. ,&
     &  .FALSE.,.TRUE. ,.FALSE.,.TRUE. ,.FALSE.,.TRUE. ,.FALSE.,.TRUE. ,&
     &  .FALSE.,.TRUE. ,.FALSE.,.TRUE. ,.FALSE.,.TRUE. ,.TRUE. ,        &
     &  .TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.FALSE./
      DATA (LSKIP(IKSUB,2),IKSUB=1,MXKSUB)     /.TRUE. ,.FALSE.,.TRUE. ,&
     &  .FALSE.,.TRUE. ,.FALSE.,.FALSE.,.TRUE. ,.FALSE.,.TRUE. ,.TRUE. ,&
     &  .FALSE.,.TRUE. ,.TRUE. ,.FALSE.,.TRUE. ,.TRUE. ,.FALSE.,.TRUE. ,&
     &  .TRUE. ,.FALSE.,.TRUE. ,.TRUE. ,.FALSE.,.TRUE. ,.TRUE. ,        &
     &  .TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.FALSE./

!     FIND AVAILABLE UNIT.
      IFILE=NUNIT()

!     TEST FOR EXISTENCE OF FILE.
      INQUIRE(FILE=CKNAME,EXIST=LTEST)
      LENCK=LENSTR(CKNAME)
      IF(.NOT.LTEST)THEN
          WRITE(IPR,'(/2A,/(11X,A))')                                   &
     &      ' WARNING:  The Correlated-k option was turned',            &
     &      ' off because file ',CKNAME(1:LENCK),'was not found.'
          RETURN
      ENDIF

!     OPEN CORRELATED-K DATA FILE.
      CALL OPNFL(IFILE,0,CKNAME,'OLD','UNFORMATTED','RDCORK')

!     READ DATA FROM BINARY FILE.
      READ(IFILE)NKSUB,NGAML,NGAMD,NNUML

!     CHECK ARRAY DIMENSIONS.
      IF(NKSUB.GT.MXKSUB .OR. NGAML.GT.MXGAML .OR.                      &
     &  NGAMD.GT.MXGAMD .OR. NNUML.GT.MXNUML)THEN
          WRITE(IPR,'(/2A,/2(10X,2A),/(10X,2(A,I3)))')' WARNING: ',     &
     &      ' The Correlated-k option was turned off.',                 &
     &      ' The Correlated-k distributions file ',CKNAME,             &
     &      ' requires minimally the following assignments',            &
     &      ' in file PARAMS.h.',                                       &
     &      '   MXKSUB =',NKSUB,'        MXGAML =',NGAML,               &
     &      '   MXGAMD =',NGAMD,'        MXNUML =',NNUML
          RETURN
      ENDIF
      READ(IFILE)                                                       &
     &  (WTKSUB(IKSUB),IKSUB=1,NKSUB),(GAMLIN(IGAML),IGAML=1,NGAML),    &
     &  (GAMDIN(IGAMD),IGAMD=1,NGAMD),(RLININ(INUML),INUML=1,NNUML)
      DO IKSUB=0,NKSUB
          READ(IFILE)(((VAL(IGAML,IGAMD,INUML,IKSUB),                   &
     &      IGAML=1,NGAML),IGAMD=1,NGAMD),INUML=1,NNUML)
      ENDDO

!     CLOSE CORRELATED-K DATA FILE.
      CLOSE(IFILE)
      IF(.NOT.LJMASS)WRITE(IPR,'(/2A)')                                 &
     &  ' SUCCESSFULLY READ CORRELATED-K DISTRIBUTIONS FILE:  ',CKNAME

!     TAILOR SELECTION OF K'S.
      IF(NSPEED.LE.0 .OR. NSPEED.GT.2)THEN

!         FOR HIGHEST RESOLUTION, USE ALL THE K'S.
          KNTRVL=NKSUB
      ELSE

!         THE FIRST AND LAST K'S CANNOT BE SKIPPED.
          KNTRVL=0
          IKSUB=1
          LSKIP(NKSUB,NSPEED)=.FALSE.
          DO IKSBP1=2,NKSUB+1
              IF(LSKIP(IKSUB,NSPEED))THEN

!                 SKIP THE "IKSUB" K'S BY ADDING WEIGHT TO "IKSUB+1".
                  WTKSUB(IKSBP1)=WTKSUB(IKSBP1)+WTKSUB(IKSUB)
              ELSE

!                 INCLUDE THE "IKSUB" K'S.
                  KNTRVL=KNTRVL+1
                  IF(KNTRVL.NE.IKSUB)THEN

!                     MOVE WEIGHTS AND K'S TO NEW LOCATION.
                      WTKSUB(KNTRVL)=WTKSUB(IKSUB)
                      DO INUML=1,NNUML
                          DO IGAMD=1,NGAMD
                              DO IGAML=1,NGAML
                                  VAL(IGAML,IGAMD,INUML,KNTRVL)         &
     &                              =VAL(IGAML,IGAMD,INUML,IKSUB)
                              ENDDO
                          ENDDO
                      ENDDO
                  ENDIF
              ENDIF
              IKSUB=IKSBP1
          ENDDO
      ENDIF

!     INITIALIZE TAIL WEIGHTS AND SAVE A COPY OF NTLSUB WTKSUB VALUES:
      DO ITLSUB=1,NTLSUB
          TAILWT(ITLSUB)=RTLSUB
          WTKSAV(ITLSUB)=WTKSUB(ITLSUB)
      ENDDO

!     FIRST DETERMINE WTKSUB FITTING IN A SINGLE INTERVAL:
      NDIVM1=NTLSUB-1
      ILAST=NTLSUB
      DO 20 INTRVL=1,KNTRVL

!         INITIALIZE K2TAIL AND CONTWT ARRAYS:
          K2TAIL(INTRVL)=0
          DO ITLSUB=1,NTLSUB
              CONTWT(ITLSUB,INTRVL)=0.
          ENDDO

!         CHECK IF WTKSUB FITS IN A SINGLE INTERVAL:
          ITLSB0=ILAST
          DO ITLSUB=ITLSB0,ITLSB0+NDIVM1
              INEXT=MOD(ITLSUB,NTLSUB)+1
              IF(WTKSUB(INTRVL).LE.TAILWT(INEXT))THEN

!                 K SUB-INTERVAL INTRVL WILL USE LINE TAIL INEXT:
                  TAILWT(INEXT)=TAILWT(INEXT)-WTKSUB(INTRVL)
                  K2TAIL(INTRVL)=INEXT
                  CONTWT(INEXT,INTRVL)=1.
                  ILAST=INEXT
                  GOTO 20
              ENDIF
          ENDDO
   20 CONTINUE

!     LOOP OVER K SUB-INTERVAL REQUIRING PARTITIONING OF WTKSUB:
!old  SUMWTK=0.
      DO 30 INTRVL=1,KNTRVL
!old      SUMWTK=SUMWTK+WTKSUB(INTRVL)
!old      CUMWTK(INTRVL)=SUMWTK
          IF(K2TAIL(INTRVL).GT.0)GOTO 30
          WGT=WTKSUB(INTRVL)
          DO ITLSUB=1,NTLSUB
              IF(WGT.GT.TAILWT(ITLSUB))THEN

!                 USE REMAINDER OF TAILWT(ITLSUB):
                  CONTWT(ITLSUB,INTRVL)=TAILWT(ITLSUB)/WTKSUB(INTRVL)
                  WGT=WGT-TAILWT(ITLSUB)
                  TAILWT(ITLSUB)=0.
              ELSE

!                 REMAINDER OF WEIGHT FITS IN ITLSUB
                  CONTWT(ITLSUB,INTRVL)=WGT/WTKSUB(INTRVL)
                  TAILWT(ITLSUB)=TAILWT(ITLSUB)-WGT
                  GOTO 30
              ENDIF
          ENDDO
   30 CONTINUE

!     WRITE WTKSUB INTO HEADER OF UNIT=IDBOUT FILE:
!old  WRITE(IDBOUT,'(4A,/90X,I2,33(1X,F13.7))')'    P[ATM]',
!old 1  '   AEROSOL_SCT   AEROSOL_EXT   AER_G     CLOUD_SCT',
!old 2  '     CLOUD_EXT   CLD_G      RAIN_SCT      RAIN_EXT',
!old 3  '  RAIN_G  RAYLEIGH_SCT  MOLECULAR_ABSORPTION...',
!old 4  KNTRVL,(CUMWTK(INTRVL),INTRVL=1,KNTRVL)

!     RETURN TO DRIVER.
      RETURN
      END