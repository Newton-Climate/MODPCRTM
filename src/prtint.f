      SUBROUTINE  PRTINT(UU, UTAU, NTAU, UMU, NUMU, PHI, NPHI, MXCLY)

!         PRINTS THE INTENSITY AT USER POLAR AND AZIMUTHAL ANGLES

!     ALL ARGUMENTS ARE DISORT INPUT OR OUTPUT VARIABLES

!+---------------------------------------------------------------------+

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

      REAL   PHI(*), UMU(*), UTAU(*), UU(MXUMU, MXCLY, *)

!     NO NEED TO PRINT IF JMASS
      IF(LJMASS) RETURN

      WRITE (*, '(//,A)')                                               &
     &         ' *********  I N T E N S I T I E S  *********'
      LENFMT=10
      NPASS=1 + NPHI / LENFMT
      IF (MOD(NPHI,LENFMT) .EQ. 0)  NPASS=NPASS - 1
      DO LU=1, NTAU
         DO NP=1, NPASS
            JMIN=1 + LENFMT * (NP-1)
            JMAX=MIN0(LENFMT*NP, NPHI)
            WRITE(*,'(13X,2(/A),10F26.2)')'   POLAR   AZIMUTH ANGLES'// &
     &        ' (DEGREES)','      OPTICAL   ANGLE',(PHI(J),J=JMIN,JMAX)
            WRITE(*,'(11A)')'        DEPTH   COSINE',                   &
     &        ('    INTENSITY     D(INTEN)',J=JMIN,JMAX)
            IF(UMU(1).GE.0.)THEN

!               DOWNLOOK:
                IF(LU.LT.NTAU)THEN
                    WRITE(*,'(1P,G13.6,0P,F8.4,1P,20E13.5)')UTAU(LU),   &
     &                UMU(1),(UU(1,LU,J),UU(1,LU,J)-UU(1,LU+1,J)        &
     &                *EXP((UTAU(LU)-UTAU(LU+1))/UMU(1)),J=JMIN,JMAX)
                ELSE
                    WRITE(*,'(1P,G13.6,0P,F8.4,1P,10(E13.5,13X))')      &
     &                UTAU(LU),UMU(1),(UU(1,LU,J),J=JMIN,JMAX)
                ENDIF
            ELSE

!               UPLOOK:
                IF(LU.GT.1)THEN
                    WRITE(*,'(1P,G13.6,0P,F8.4,1P,20E13.5)')UTAU(LU),   &
     &                UMU(1),(UU(1,LU,J),UU(1,LU,J)-UU(1,LU-1,J)        &
     &                *EXP((UTAU(LU)-UTAU(LU-1))/UMU(1)),J=JMIN,JMAX)
                ELSE
                    WRITE(*,'(1P,G13.6,0P,F8.4,1P,10(E13.5,13X))')      &
     &                UTAU(LU),UMU(1),(UU(1,LU,J),J=JMIN,JMAX)
                ENDIF
            ENDIF
            DO IU=2,NUMU
                IF(UMU(IU).GE.0.)THEN

!                   DOWNLOOK:
                    IF(LU.LT.NTAU)THEN
                        WRITE(*,'(13X,F8.4,1P,20E13.5)')UMU(IU),        &
     &                    (UU(IU,LU,J),UU(IU,LU,J)-UU(IU,LU+1,J)*EXP    &
     &                    ((UTAU(LU)-UTAU(LU+1))/UMU(IU)),J=JMIN,JMAX)
                    ELSE
                        WRITE(*,'(13X,F8.4,1P,10(E13.5,13X))')          &
     &                    UMU(IU),(UU(IU,LU,J),J=JMIN,JMAX)
                    ENDIF
                ELSE

!                   UPLOOK:
                    IF(LU.GT.1)THEN
                        WRITE(*,'(13X,F8.4,1P,20E13.5)')UMU(IU),        &
     &                    (UU(IU,LU,J),UU(IU,LU,J)-UU(IU,LU-1,J)*EXP    &
     &                    ((UTAU(LU)-UTAU(LU-1))/UMU(IU)),J=JMIN,JMAX)
                    ELSE
                        WRITE(*,'(13X,F8.4,1P,10(E13.5,13X))')          &
     &                    UMU(IU),(UU(IU,LU,J),J=JMIN,JMAX)
                    ENDIF
                ENDIF
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
