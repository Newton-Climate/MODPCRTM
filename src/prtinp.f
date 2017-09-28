      SUBROUTINE  PRTINP(HEADER, NLYR, DTAUC, SSALB, PMOM, TEMPER,      &
     &                    WVNMLO, WVNMHI, NTAU, UTAU, NSTR, NUMU, UMU,  &
     &                    NPHI, PHI, IBCND, FBEAM, UMU0, PHI0,          &
     &                    FISOT, LAMBER, ALBEDO, BTEMP, TTEMP,          &
     &                    TEMIS, DELTAM, PLANK, ONLYFL, ACCUR,          &
     &                    FLYR, LYRCUT, OPRIM, TAUC, TAUCPR, PRTMOM)

!                                        PRINT VALUES OF INPUT VARIABLES

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

      CHARACTER  HEADER*127
      LOGICAL  DELTAM, LAMBER, LYRCUT, PLANK, ONLYFL, PRTMOM
      REAL     UMU(*), FLYR(*), DTAUC(*), OPRIM(*), PHI(*),             &
     &         PMOM(0:MAXCOE,*), SSALB(*), UTAU(*), TAUC(0:*),          &
     &         TAUCPR(0:*), TEMPER(0:*)

!     NO NEED TO PRINT IF JMASS
      IF(LJMASS) RETURN

      WRITE(*,1000)  HEADER
      WRITE(*,1010)  NSTR, NLYR
      IF(IBCND.NE.1)WRITE(*,'(I4,A,1P,8G13.6,/(26X,8G13.6))')           &
     &  NTAU,' USER OPTICAL DEPTHS :',(UTAU(LU),LU=1,NTAU)
      IF (.NOT.ONLYFL)                                                  &
     &      WRITE(*,1030)  NUMU, (UMU(IU), IU=1, NUMU)
      IF (.NOT.ONLYFL .AND. IBCND.NE.1)                                 &
     &      WRITE(*,1040)  NPHI, (PHI(J), J=1, NPHI)
      IF (.NOT.PLANK .OR. IBCND.EQ.1)  WRITE(*,1100)
      WRITE(*,1050)  IBCND
      IF (IBCND.EQ.0)  THEN
         WRITE(*,1060) FBEAM, UMU0, PHI0, FISOT
         IF (LAMBER)   WRITE(*,1080) ALBEDO
!        IF (.NOT.LAMBER)  WRITE(*,1090) (HL(K), K=0, NSTR)
         IF (PLANK)  WRITE(*,1110) WVNMLO, WVNMHI, BTEMP,               &
     &                                 TTEMP, TEMIS
      ELSE IF (IBCND.EQ.1)  THEN
         WRITE(*,1070)
         WRITE(*,1080) ALBEDO
      ENDIF
      IF (DELTAM)      WRITE(*,1120)
      IF (.NOT.DELTAM) WRITE(*,1130)
      IF (IBCND.EQ.1)  THEN
         WRITE(*,1140)
      ELSE IF (ONLYFL)  THEN
         WRITE(*,1150)
      ELSE
         WRITE(*,1160)
      ENDIF
      WRITE(*,1170)  ACCUR
      IF (LYRCUT)  WRITE(*,1180)
      YESSCT=0.
      WRITE(*,'(/43X,A,/(A))')                                          &
     &      '<---------------- DELTA-M ------------------>',            &
     &      '                         TOTAL    SINGLE              '//  &
     &      '                   TOTAL    SINGLE',                       &
     &      '          OPTICAL      OPTICAL   SCATTER   TRUNCATED  '//  &
     &      '    OPTICAL      OPTICAL   SCATTER    ASYMM'
      IF(PLANK)THEN
          WRITE(*,'(A)')                                                &
     &      '            DEPTH        DEPTH    ALBEDO    FRACTION  '//  &
     &      '      DEPTH        DEPTH    ALBEDO   FACTOR   TEMPERATURE'
          DO LC=1,NLYR
             YESSCT=YESSCT + SSALB(LC)
             WRITE(*,'(I4,1P,2G13.6,0P,F10.5,F12.5,1P,2G13.6,           &
     &         0P,F10.5,F9.4,F14.3)')LC,DTAUC(LC),TAUC(LC),             &
     &         SSALB(LC),FLYR(LC),TAUCPR(LC)-TAUCPR(LC-1),              &
     &         TAUCPR(LC),OPRIM(LC),PMOM(1,LC),TEMPER(LC-1)
          ENDDO
          WRITE(*,'(97X,F14.3)')TEMPER(NLYR)
      ELSE
          WRITE(*,'(A)')                                                &
     &      '            DEPTH        DEPTH    ALBEDO    FRACTION  '//  &
     &      '      DEPTH        DEPTH    ALBEDO   FACTOR'
          DO LC=1,NLYR
             YESSCT=YESSCT + SSALB(LC)
             WRITE(*,'(I4,1P,2G13.6,0P,F10.5,F12.5,1P,2G13.6,           &
     &         0P,F10.5,F9.4,F14.3)')LC,DTAUC(LC),TAUC(LC),             &
     &         SSALB(LC),FLYR(LC),TAUCPR(LC)-TAUCPR(LC-1),              &
     &         TAUCPR(LC),OPRIM(LC),PMOM(1,LC)
          ENDDO
      ENDIF

      IF(PRTMOM .AND. YESSCT.GT.0.)THEN
         WRITE(*, '(/,A)')  ' LAYER   PHASE FUNCTION MOMENTS'
         DO 20 LC=1, NLYR
            IF(SSALB(LC).GT.0.0)                                        &
     &          WRITE(*,1230)  LC, (PMOM(K,LC), K=0, NSTR)
20       CONTINUE
      END IF

      RETURN

1000  FORMAT(////, 1X, 120('*'), /, 25X,                                &
     &  'DISCRETE ORDINATES RADIATIVE TRANSFER PROGRAM, VERSION 2.0',   &
     &  /, 1X, A, /, 1X, 120('*'))
1010  FORMAT(/, ' NO. STREAMS =', I4,                                   &
     &  '     NO. COMPUTATIONAL LAYERS =', I4)
1030  FORMAT(I4,' USER POLAR ANGLE COSINES :',10F9.5,/,(31X,10F9.5))
1040  FORMAT(I4,' USER AZIMUTHAL ANGLES :', 10F9.2, /, (28X,10F9.2))
1050  FORMAT(' BOUNDARY CONDITION FLAG: IBCND =', I2)
1060  FORMAT('    INCIDENT BEAM WITH INTENSITY =', 1P,E11.3, ' AND',    &
     & ' POLAR ANGLE COSINE=', 0P,F8.5,'  AND AZIMUTH ANGLE =', F7.2,   &
     & /,'    PLUS ISOTROPIC INCIDENT INTENSITY =', 1P,E11.3)
1070  FORMAT('    ISOTROPIC ILLUMINATION FROM TOP AND BOTTOM')
1080  FORMAT('    BOTTOM ALBEDO (LAMBERTIAN) =', 0P,F8.4)
!1090 FORMAT('    LEGENDRE COEFFS OF BOTTOM BIDIRECTIONAL',
!    $ ' REFLECTIVITY :', /, (10X,10F9.5))
1100  FORMAT(' NO THERMAL EMISSION')
1110  FORMAT('    THERMAL EMISSION IN WAVENUMBER INTERVAL :', 2F14.4,/, &
     &   '    BOTTOM TEMPERATURE =', F10.2, '     TOP TEMPERATURE =',   &
     &   F10.2,'    TOP EMISSIVITY =', F8.4)
1120  FORMAT(' USES DELTA-M METHOD')
1130  FORMAT(' DOES NOT USE DELTA-M METHOD')
1140  FORMAT(' CALCULATE ALBEDO AND TRANSMISSIVITY OF MEDIUM',          &
     &   ' VS. INCIDENT BEAM ANGLE')
1150  FORMAT(' CALCULATE FLUXES AND AZIM-AVERAGED INTENSITIES ONLY')
1160  FORMAT(' CALCULATE FLUXES AND INTENSITIES')
1170  FORMAT(' RELATIVE CONVERGENCE CRITERION FOR AZIMUTH SERIES =',    &
     &   1P,E11.2)
1180  FORMAT(' SETS RADIATION=0 BELOW ABSORPTION OPTICAL DEPTH 10')
1230  FORMAT(I6, 10F11.6, /, (6X,10F11.6))

      END
