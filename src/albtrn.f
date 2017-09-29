      SUBROUTINE  ALBTRN( ALBEDO, AMB, APB, ARRAY, B, CBAND, CC,        &
     &                    CMU, CWT, EVAL, EVECC, GL, GC, GU, IPVT,      &
     &                    KK, LL, NLYR, NN, NSTR, NUMU, PRNT, TAUCPR,   &
     &                    UMU, U0U, WK, YLMC, YLMU, Z, MI9M2,           &
     &                    MAXULV, NNLYRI, WN0, ALBMED, TRNMED)

!        SPECIAL CASE TO GET ONLY ALBEDO AND TRANSMISSIVITY
!        OF ENTIRE MEDIUM AS A FUNCTION OF INCIDENT BEAM ANGLE
!        (MANY SIMPLIFICATIONS BECAUSE BOUNDARY CONDITION IS JUST
!        ISOTROPIC ILLUMINATION, THERE ARE NO THERMAL SOURCES, AND
!        PARTICULAR SOLUTIONS DO NOT NEED TO BE COMPUTED).  SEE
!        REF. S2 AND REFERENCES THEREIN FOR THEORY.

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!        ROUTINES CALLED:  ALTRIN, LEPOLY, PRALTR, SETMTX, SOLVE1,
!                          SOLEIG, ZEROIT

!     /RCNSTN/
!       PI       THE CONSTANT PI
!       DEG      NUMBER OF DEGREES IN ONE RADIAN.
!       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.
!       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.
!       RRIGHT   SMALLEST SINGLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
      REAL PI,DEG,BIGNUM,BIGEXP,RRIGHT
      COMMON/RCNSTN/PI,DEG,BIGNUM,BIGEXP,RRIGHT

      LOGICAL  PRNT(*)
      INTEGER  NLYR, NUMU, NSTR
      REAL       UMU(*), U0U( MXUMU,* ), TRNMD2(2)

      INTEGER IPVT(*)
      REAL    WN0, ALBMED(*), AMB(MI,*), APB(MI,*), ARRAY(MXCMU,*),     &
     &        B(*), BDR( MI,0:MI), CBAND( MI9M2,* ), CC( MXCMU,* ),     &
     &        CMU(*), CWT(*), EVAL(*), EVECC( MXCMU,* ),                &
     &  GL(0:MXCMU,*),GC(MXCMU,MXCMU,*),GU(0:MXUMU,MXCMU,*),            &
     &        KK( MXCMU,* ), LL( MXCMU,* ), TAUCPR( 0:* ), TRNMED(*),   &
     &  WK(*),YLMC(0:MXCMU,*),YLMU(0:MXCMU,0:*),Z(*)

      LOGICAL  LAMBER, LYRCUT

!     FUNCTIONS:
      INTEGER NUNIT

!     DATA:
!       NSPALB   FILE UNIT NUMBER FOR ATMOSPHERIC SPHERICAL ALBEDO.
!                (<0   OPEN A "SPHERALB.DAT" FILE             )
!                (>0   UNIT NUMBER FOR THE "SPHERALB.DAT" FILE)
!       LWRITE   WRITING FLAG TOGGLE.
      INTEGER NSPALB
      LOGICAL LWRITE
      SAVE NSPALB,LWRITE,UMUSUN,TMD2SV
      DATA NSPALB/-1/,LWRITE/.FALSE./

!     ATMOSPHERIC SPHERICAL ALBEDO FILE:
      IF(NSPALB.LT.0.)THEN
          NSPALB=NUNIT()
          CALL OPNFL(NSPALB,0,'spheralb.dat','UNKNOWN','FORMATTED',     &
     &      'ALBTRN')

!         HEADER:
          IF(.NOT.LJMASS)WRITE(NSPALB,'((2A))')                         &
     &      '   FREQ     COSINE    COSINE    SPHERICAL',                &
     &      '   TOTAL (DIRECT+DIFFUSE) TRANSMITANCES',                  &
     &      '  (CM-1)       SUN     NADIR      ALBEDO ',                &
     &      '      SUN-GND    GND-SPACE  SUN-GND-SPC',                  &
     &      '  ------   -------   -------  -----------',                &
     &      '  -----------  -----------  -----------'
      ENDIF

!                    ** SET DISORT VARIABLES THAT ARE IGNORED IN THIS
!                    ** SPECIAL CASE BUT ARE NEEDED BELOW IN ARGUMENT
!                    ** LISTS OF SUBROUTINES SHARED WITH GENERAL CASE
      NCUT = NLYR
      LYRCUT = .FALSE.
      FISOT = 1.0
      LAMBER = .TRUE.

      MAZIM = 0
      DELM0 = 1.0
!                          ** GET LEGENDRE POLYNOMIALS FOR COMPUTATIONAL
!                          ** AND USER POLAR ANGLE COSINES

      CALL LEPOLY(NUMU,MAZIM,MXCMU,NSTR-1,UMU,YLMU(0,1))
      CALL LEPOLY(NN,  MAZIM,MXCMU,NSTR-1,CMU,YLMC(0,1))

!                       ** EVALUATE LEGENDRE POLYNOMIALS WITH NEGATIVE
!                       ** -CMU- FROM THOSE WITH POSITIVE -CMU-;
!                       ** DAVE/ARMSTRONG EQ. (15)
      SGN  = -1.0
      DO  L = MAZIM, NSTR-1
         SGN = - SGN
         DO  IQ = NN+1, NSTR
            YLMC( L,IQ ) = SGN * YLMC( L,IQ-NN )
         ENDDO
      ENDDO
!                                  ** ZERO BOTTOM REFLECTIVITY
!                                  ** (-ALBEDO- IS USED ONLY IN ANALYTIC
!                                  ** FORMULAE INVOLVING ALBEDO = 0
!                                  ** SOLUTIONS; EQS 16-17 OF REF S2)
      CALL  ZEROIT( BDR(1,0), MI*(MI+1) )

! ===================  BEGIN LOOP ON COMPUTATIONAL LAYERS  =============

      DO LC = 1, NLYR

!                        ** SOLVE EIGENFUNCTION PROBLEM IN EQ. STWJ(8B)

         CALL  SOLEIG( AMB, APB, ARRAY, CMU, CWT, GL(0,LC), MAZIM,      &
     &     NN, NSTR, WK, YLMC, CC, EVECC, EVAL, KK(1,LC), GC(1,1,LC))

!                          ** INTERPOLATE EIGENVECTORS TO USER ANGLES

         CALL TERPEV(CWT,EVECC,GL(0,LC),GU(0,1,LC),MAZIM,MXCMU,         &
     &     MXUMU,NN,NSTR,1,NUMU,WK,YLMC,YLMU)
      ENDDO

! ===================  END LOOP ON COMPUTATIONAL LAYERS  ===============

!                      ** SET COEFFICIENT MATRIX OF EQUATIONS COMBINING
!                      ** BOUNDARY AND LAYER INTERFACE CONDITIONS

      CALL  SETMTX( BDR(1,0), CBAND(1,1), CMU, CWT, DELM0, GC, KK,      &
     &              LAMBER, LYRCUT, MI, MI9M2, MXCMU, NCOL, NCUT,       &
     &              NNLYRI, NN, NSTR, TAUCPR, WK )

      CALL  ZEROIT( U0U(1,1), MXUMU*MAXULV )

      NHOM = 2
      IF( NLYR.EQ.1 )  NHOM = 1
      SPHALB = 0.0
      SPHTRN = 0.0
      DO IHOM = 1, NHOM
!                             ** SOLVE FOR CONSTANTS OF INTEGRATION IN
!                             ** HOMOGENEOUS SOLUTION FOR ILLUMINATION
!                             ** FROM TOP (IHOM=1), THEN BOTTOM (IHOM=2)

         CALL  SOLVE1( B, CBAND(1,1), FISOT, IHOM, IPVT, LL, MI9M2,     &
     &                 MXCMU, NCOL, NLYR, NN, NNLYRI, NSTR, Z, .FALSE. )

!                             ** COMPUTE AZIMUTHALLY-AVERAGED INTENSITY
!                             ** AT USER ANGLES; GIVES ALBEDO IF MULTI-
!                             ** LAYER (EQ. 9 OF REF S2); GIVES BOTH
!                             ** ALBEDO AND TRANSMISSIVITY IF SINGLE
!                             ** LAYER (EQS. 3-4 OF REF S2)

         CALL  ALTRIN( GU, KK, LL, MXCMU, MXUMU, NLYR,                  &
     &                 NN, NSTR, NUMU, TAUCPR, UMU, U0U(1,1), WK )

         IF ( IHOM.EQ.1 )  THEN
!                                   ** SAVE ALBEDOS;  FLIP TRANSMISSIV.
!                                   ** END OVER END TO CORRESPOND TO
!                                   ** POSITIVE -UMU- INST. OF NEGATIVE
            DO IU = 1, NUMU/2
               ALBMED(IU) = U0U( IU + NUMU/2, 1 )
               IF( NLYR.EQ.1 )  TRNMED(IU) = U0U( NUMU/2+1-IU, 2 )      &
     &                          + EXP( - TAUCPR(NLYR) / UMU(IU+NUMU/2) )
               TRNMD2(IU) = U0U( NUMU/2+1-IU, 2 )                       &
     &                          + EXP( - TAUCPR(NLYR) / UMU(IU+NUMU/2) )
            ENDDO
!                                    ** GET SPHERICAL ALBEDO AND, FOR 1
!                                    ** LAYER, SPHERICAL TRANSMISSIVITY
            IF( ALBEDO.GT.0.0 )                                         &
     &          CALL SPALTR( CMU, CWT, GC, KK, LL, MXCMU, NLYR,         &
     &                       NN, NSTR, TAUCPR, SPHALB, SPHTRN )

         ELSE IF ( IHOM.EQ.2 )  THEN
!                                      ** SAVE TRANSMISSIVITIES
            DO IU = 1, NUMU/2
               TRNMED(IU) = U0U( IU + NUMU/2, 1 )                       &
     &                      + EXP( - TAUCPR(NLYR) / UMU(IU+NUMU/2) )
            ENDDO
!                             ** GET SPHERICAL ALBEDO AND TRANSMISSIVITY
            IF( ALBEDO.GT.0.0 )                                         &
     &          CALL SPALTR( CMU, CWT, GC, KK, LL, MXCMU, NLYR,         &
     &                       NN, NSTR, TAUCPR, SPHTRN, SPHALB )
         END IF
      ENDDO

      IF ( ALBEDO.GT.0.0 )  THEN
!                                ** REF. S2, EQS. 16-17 (THESE EQS. HAVE
!                                ** A SIMPLE PHYSICAL INTERPRETATION
!                                ** LIKE THAT OF THE DOUBLING EQS.)
      IF(LWRITE)THEN
         IF(.NOT.LJMASS)WRITE(NSPALB,'(F8.0,2F10.5,4F13.7)')            &
     &     WN0,UMUSUN,-UMU(1),SPHALB,TMD2SV,TRNMED(1),TMD2SV*TRNMED(1)
         LWRITE=.FALSE.
      ELSE
         UMUSUN=-UMU(1)
         TMD2SV=TRNMD2(1)
         LWRITE=.TRUE.
      ENDIF
         DO IU = 1, NUMU
            ALBMED(IU) = ALBMED(IU) + ( ALBEDO / (1.-ALBEDO*SPHALB) )   &
     &                                * SPHTRN * TRNMED(IU)
            TRNMED(IU) = TRNMED(IU) + ( ALBEDO / (1.-ALBEDO*SPHALB) )   &
     &                                * SPHALB * TRNMED(IU)
         ENDDO
      END IF
!                          ** RETURN -UMU- TO ALL POSITIVE VALUES, TO
!                          ** AGREE WITH ORDERING IN -ALBMED,TRNMED-
      NUMU = NUMU / 2
      DO IU = 1, NUMU
        UMU(IU) = UMU(IU+NUMU)
      ENDDO

      IF(.NOT.LJMASS .AND. PRNT(6))                                     &
     &   WRITE(*,'(/2(//A),/(0P,F13.4,F20.6,F12.5,1P,E17.4))')          &
     &  ' ***  FLUX ALBEDO AND/OR TRANSMISSIVITY OF ENTIRE MEDIUM  ***',&
     &  ' BEAM ZEN ANG  dcos(BEAM ZEN ANG)      ALBEDO  TRANSMISSIVITY',&
     &  (DEG*ACOS(UMU(IU)),UMU(IU),ALBMED(IU),TRNMED(IU),IU=1,NUMU)

      RETURN
      END
