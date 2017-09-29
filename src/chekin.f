      SUBROUTINE CHEKIN(ACCUR,ALBEDO,BTEMP,DTAUC,FBEAM,FISOT,IBCND,     &
     &  LAMBER,MAXCOE,MAXPHI,MAXULV,MXCLY,MXCMU,MXULV,MXUMU,NLYR,       &
     &  NPHI,NSTR,NTAU,NUMU,ONLYFL,PHI,PHI0,PLANK,PMOM,SSALB,TAUC,      &
     &  TEMIS,TEMPER,TTEMP,UMU,UMU0,USRANG,USRTAU,UTAU,WVNMHI,WVNMLO)

!     CHECKS THE INPUT DIMENSIONS AND VARIABLES
      LOGICAL  WRTBAD, WRTDIM
      LOGICAL  LAMBER, PLANK, ONLYFL, USRANG, USRTAU, INPERR
      INTEGER IBCND,MAXULV,MAXCOE,MAXPHI,NLYR,NUMU,NSTR,NPHI,NTAU,      &
     &  MXCMU,MXUMU,MXCLY,MXULV
      REAL     ACCUR, ALBEDO, BTEMP, DTAUC(*), FBEAM, FISOT,            &
     &         PHI(*), PMOM(0:MAXCOE,*), PHI0,                          &
     &         SSALB(*), TAUC(0:*), TEMPER(0:*), TEMIS,                 &
     &         TTEMP, WVNMLO, WVNMHI, UMU(*), UMU0, UTAU(*)
      LOGICAL MODTRN
      SAVE MODTRN
      DATA MODTRN/.TRUE./

      INPERR=.FALSE.
      IF (NLYR.LT.1) INPERR=WRTBAD('NLYR')
      IF (NLYR.GT.NTAU) INPERR=WRTBAD('NTAU')

      DO 10 LC=1, NLYR
           IF (DTAUC(LC).LT.0.0) INPERR=WRTBAD('DTAUC')
           IF (SSALB(LC).LT.0.0 .OR. SSALB(LC).GT.1.0)                  &
     &        INPERR=WRTBAD('SSALB')
           IF (PLANK .AND. IBCND.NE.1)  THEN
              IF(LC.EQ.1 .AND. TEMPER(0).LT.0.0)                        &
     &          INPERR=WRTBAD('TEMPER')
              IF(TEMPER(LC).LT.0.0) INPERR=WRTBAD('TEMPER')
           ENDIF
           DO 5 K=0, NSTR
              IF(PMOM(K,LC).LT.-1.0 .OR. PMOM(K,LC).GT.1.0)             &
     &          INPERR=WRTBAD('PMOM')
    5      CONTINUE
10    CONTINUE

      IF (IBCND.EQ.1)  THEN
           IF (MAXULV.LT.2) INPERR=WRTBAD('MAXULV')
      ELSE IF (USRTAU)  THEN
           IF (NTAU.LT.1) INPERR=WRTBAD('NTAU')
           IF (MAXULV.LT.NTAU) INPERR=WRTBAD('MAXULV')
           DO 20  LU=1, NTAU
              IF(ABS(UTAU(LU)-TAUC(NLYR)).LE.1.E-4) UTAU(LU)=TAUC(NLYR)
              IF(UTAU(LU).LT.0.0 .OR. UTAU(LU).GT.TAUC(NLYR))           &
     &          INPERR=WRTBAD('UTAU')
20       CONTINUE
      ELSE
           IF (MAXULV.LT.NLYR+1) INPERR=WRTBAD('MAXULV')
      END IF

      IF (NSTR.LT.2 .OR. MOD(NSTR,2).NE.0) INPERR=WRTBAD('NSTR')
      IF (NSTR.GT.MAXCOE) INPERR=WRTBAD('MAXCOE')

      IF (USRANG)  THEN
           IF (NUMU.LT.0) INPERR=WRTBAD('NUMU')
           IF (.NOT.ONLYFL .AND. NUMU.EQ.0) INPERR=WRTBAD('NUMU' )
           IF (NUMU.GT.MXUMU) INPERR=WRTBAD('MXUMU')
           IF (IBCND.EQ.1 .AND. 2*NUMU.GT.MXUMU)                        &
     &        INPERR=WRTBAD('MXUMU')
           DO 30 IU=1, NUMU
              IF(UMU(IU).LT.-1.0 .OR. UMU(IU).GT.1.0.OR.UMU(IU).EQ.0.0) &
     &           INPERR=WRTBAD('UMU')
              IF(IBCND.EQ.1 .AND. UMU(IU).LT.0.0)                       &
     &           INPERR=WRTBAD('UMU')
              IF(IU.GT.1) THEN
                 IF (UMU(IU).LT.UMU(IU-1)) INPERR=WRTBAD('UMU')
              ENDIF
 30        CONTINUE
      ELSE
           IF(MXUMU.LT.NSTR) INPERR=WRTBAD('MXUMU')
      END IF

      IF (.NOT.ONLYFL .AND. IBCND.NE.1)  THEN
           IF (NPHI.LE.0) INPERR=WRTBAD('NPHI')
           IF (NPHI.GT.MAXPHI) INPERR=WRTBAD('MAXPHI')
           DO 40 J=1, NPHI
              IF (PHI(J).LT.0.0 .OR. PHI(J).GT.360.0)                   &
     &           INPERR=WRTBAD('PHI')
40       CONTINUE
      END IF

      IF (IBCND.LT.0 .OR. IBCND.GT.1) INPERR=WRTBAD('IBCND')
      IF (IBCND.EQ.0)  THEN
           IF (FBEAM.LT.0.0) INPERR=WRTBAD('FBEAM')
           IF (FBEAM.GT.0.0 .AND. (UMU0.LE.0.0 .OR. UMU0.GT.1.0))       &
     &        INPERR=WRTBAD('UMU0')
           IF (FBEAM.GT.0.0 .AND. (PHI0.LT.0.0 .OR. PHI0.GT.360.0))     &
     &        INPERR=WRTBAD('PHI0')
           IF (FISOT.LT.0.0) INPERR=WRTBAD('FISOT')
           IF (LAMBER)  THEN
              IF (ALBEDO.LT.0.0 .OR. ALBEDO.GT.1.0)                     &
     &           INPERR=WRTBAD('ALBEDO')
!        ELSE
!                        MAKE SURE FLUX ALBEDO AT DENSE MESH OF INCIDENT
!                               ANGLES DOES NOT ASSUME UNPHYSICAL VALUES

!           DO 50 RMU=0.0, 1.0, 0.01
!              FLXALB=DREF(RMU, HL, NSTR)
!              IF (FLXALB.LT.0.0 .OR. FLXALB.GT.1.0)
!    $            INPERR=WRTBAD('HL')
!  50       CONTINUE
           ENDIF

      ELSE IF (IBCND.EQ.1)  THEN
           IF (ALBEDO.LT.0.0 .OR. ALBEDO.GT.1.0)                        &
     &        INPERR=WRTBAD('ALBEDO')
      END IF

      IF (PLANK .AND. IBCND.NE.1)  THEN
           IF (WVNMLO.LT.0.0 .OR. WVNMHI.LE.WVNMLO)                     &
     &        INPERR=WRTBAD('WVNMLO,HI')
           IF (TEMIS.LT.0.0 .OR. TEMIS.GT.1.0)                          &
     &        INPERR=WRTBAD('TEMIS')
           IF (BTEMP.LT.0.0) INPERR=WRTBAD('BTEMP')
           IF (TTEMP.LT.0.0) INPERR=WRTBAD('TTEMP')
      END IF

      IF (ACCUR.LT.0.0 .OR. ACCUR.GT.1.E-2)                             &
     &     INPERR=WRTBAD('ACCUR')

      IF (MXCLY.LT.NLYR) INPERR=WRTDIM('MXCLY', NLYR)
      IF (IBCND.NE.1)  THEN
           IF (USRTAU .AND. MXULV.LT.NTAU)                              &
     &        INPERR=WRTDIM('MXULV', NTAU)
           IF (.NOT.USRTAU .AND. MXULV.LT.NLYR+1)                       &
     &        INPERR=WRTDIM('MXULV', NLYR+1)
      ELSE
           IF (MXULV.LT.2) INPERR=WRTDIM('MXULV', 2)
      END IF

      IF (MXCMU.LT.NSTR) INPERR=WRTDIM('MXCMU', NSTR)
      IF (USRANG .AND. MXUMU.LT.NUMU)                                   &
     &     INPERR=WRTDIM('MXUMU', NUMU)
      IF (USRANG .AND. IBCND.EQ.1 .AND. MXUMU.LT.2*NUMU)                &
     &     INPERR=WRTDIM('MXUMU', NUMU)
      IF (.NOT.USRANG .AND. MXUMU.LT.NSTR)                              &
     &      INPERR=WRTDIM('MXUMU', NSTR)

      IF (INPERR)                                                       &
     &   CALL ERRMSG('DISORT--input and/or dimension errors', .TRUE.)
      IF(MODTRN)RETURN
      IF (PLANK)  THEN
        DO 60 LC=1, NLYR
           IF (ABS(TEMPER(LC)-TEMPER(LC-1)) .GT. 20.0)                  &
     &        CALL ERRMSG('CHEKIN--vertical temperature step may'       &
     &                  // ' be too large for good accuracy', .FALSE.)
 60     CONTINUE
      END IF

      RETURN
      END
