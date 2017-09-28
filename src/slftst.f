      SUBROUTINE  SLFTST(ACCUR, ALBEDO, BPLANK, BTEMP, DELTAM, DTAUC,   &
     &                    FBEAM, FISOT, IBCND, CORINT, LAMBER, NLYR,    &
     &                    PLANK, NPHI, NUMU, NSTR, NTAU, ONLYFL, PHI,   &
     &                    PHI0, PKAG, PMOM, PRNT, SSALB, TEMIS, TEMPER, &
     &                    TPLANK, TTEMP, UMU, USRANG, USRTAU, UTAU,     &
     &                    UMU0, WVNMHI, WVNMLO, LTEST,                  &
     &                    FLUP, RFLDIR, RFLDN, UU)

!       IF  LTEST=FALSE, SAVE USER INPUT VALUES THAT WOULD OTHERWISE
!       BE DESTROYED AND REPLACE THEM WITH INPUT VALUES FOR SELF-TEST.
!       IF  LTEST=TRUE, COMPARE SELF-TEST CASE RESULTS WITH CORRECT
!       ANSWERS AND RESTORE USER INPUT VALUES IF TEST IS PASSED.

!       (SEE FILE 'DisORT.Doc' FOR VARIABLE DEFINITIONS.)

!                 I N T E R N A L    V A R I A B L E S:

!         ACC     RELATIVE ACCURACY REQUIRED FOR PASSING SELF-TEST
!         ERRORN  RELATIVE ERRORS IN 'DisORT' OUTPUT VARIABLES
!         OK      LOGICAL VARIABLE FOR DETERMINING FAILURE OF SELF-TEST
!         ALL VARIABLES ENDING IN 'S':  TEMPORARY 'S'TORAGE FOR INPUT
!+---------------------------------------------------------------------+
      REAL     PMOM(0:*), TEMPER(0:*), PKAG(0:*)
      LOGICAL  LTEST, CORINT, DELTAM, LAMBER, PLANK, OK, ONLYFL,        &
     &         PRNT(*), USRANG, USRTAU
      REAL     PMOMS(0:4), TEMPES (0:1), PKAGS(0:1)
      LOGICAL  CORINS, DELTAS, LAMBES, PLANKS, ONLYFS, PRNTS (7),       &
     &         USRANS, USRTAS, TSTBAD
      SAVE
      DATA     ACC / 1.E-4 /

      IF  (.NOT.LTEST)  THEN
!                                                 SAVE USER INPUT VALUES
         NLYRS =NLYR
         DTAUCS=DTAUC
         SSALBS=SSALB
         DO N=0, 4
            PMOMS(N)=PMOM(N)
         ENDDO
         NSTRS =NSTR
         USRANS=USRANG
         NUMUS =NUMU
         UMUS  =UMU
         USRTAS=USRTAU
         NTAUS =NTAU
         UTAUS =UTAU
         NPHIS =NPHI
         PHIS  =PHI
         IBCNDS=IBCND
         CORINS=CORINT
         FBEAMS=FBEAM
         UMU0S =UMU0
         PHI0S =PHI0
         FISOTS=FISOT
         LAMBES=LAMBER
         ALBEDS=ALBEDO
         DELTAS=DELTAM
         ONLYFS=ONLYFL
         ACCURS=ACCUR
         PLANKS=PLANK
         WVNMLS=WVNMLO
         WVNMHS=WVNMHI
         BTEMPS=BTEMP
         TTEMPS=TTEMP
         TEMISS=TEMIS
         BPLNKS=BPLANK
         TPLNKS=TPLANK
         TEMPES(0)=TEMPER(0)
         TEMPES(1)=TEMPER(1)
         PKAGS(0)=PKAG(0)
         PKAGS(1)=PKAG(1)
         DO I=1, 7
            PRNTS(I)=PRNT(I)
         ENDDO
!                                         SET INPUT VALUES FOR SELF-TEST
         NLYR =1
         DTAUC=1.0
         SSALB=0.9
!                                                         HAZE-L MOMENTS
         PMOM(0)=1.0
         PMOM(1)=0.8042
         PMOM(2)=0.646094
         PMOM(3)=0.481851
         PMOM(4)=0.359056
         NSTR  =4
         USRANG=.TRUE.
         NUMU  =1
         UMU   =0.5
         USRTAU=.TRUE.
         NTAU  =1
         UTAU  =0.5
         NPHI  =1
         PHI   =90.0
         IBCND =0
         CORINT=.FALSE.
         FBEAM =3.1415927
         UMU0  =0.866
         PHI0  =0.0
         FISOT =1.0
         LAMBER=.TRUE.
         ALBEDO=0.7
         DELTAM=.TRUE.
         ONLYFL=.FALSE.
         ACCUR =1.E-4
         PLANK =.TRUE.
         WVNMLO=0.0
         WVNMHI=50000.
         BTEMP =300.0
         TTEMP =100.0
         TEMIS =0.8
         TEMPER(0)=210.0
         TEMPER(1)=200.0
         DO I=1, 7
            PRNT(I)=.FALSE.
         ENDDO
         TPLANK=TEMIS * PLKAVG(WVNMLO, WVNMHI, TTEMP)
         BPLANK=        PLKAVG(WVNMLO, WVNMHI, BTEMP)
         PKAG(0)=PLKAVG(WVNMLO, WVNMHI, TEMPER(0))
         PKAG(1)=PLKAVG(WVNMLO, WVNMHI, TEMPER(1))

      ELSE

!        COMPARE TEST CASE RESULTS WITH CORRECT ANSWERS AND ABORT IF BAD

         OK=.TRUE.
         ERROR1=(UU  - 47.86005) / 47.86005
         ERROR2=(RFLDIR - 1.527286) / 1.527286
         ERROR3=(RFLDN - 28.37223) / 28.37223
         ERROR4=(FLUP   - 152.5853) / 152.5853
         IF(ABS(ERROR1).GT.ACC) OK=TSTBAD('UU',     ERROR1)
         IF(ABS(ERROR2).GT.ACC) OK=TSTBAD('RFLDIR', ERROR2)
         IF(ABS(ERROR3).GT.ACC) OK=TSTBAD('RFLDN',  ERROR3)
         IF(ABS(ERROR4).GT.ACC) OK=TSTBAD('FLUP',   ERROR4)

         IF(.NOT. OK)                                                   &
     &       CALL ERRMSG('DISORT--SELF-TEST FAILED', .TRUE.)

!                                              RESTORE USER INPUT VALUES
         NLYR =NLYRS
         DTAUC=DTAUCS
         SSALB=SSALBS
         DO N=0, 4
            PMOM(N)=PMOMS(N)
         ENDDO
         NSTR  =NSTRS
         USRANG=USRANS
         NUMU  =NUMUS
         UMU   =UMUS
         USRTAU=USRTAS
         NTAU  =NTAUS
         UTAU  =UTAUS
         NPHI  =NPHIS
         PHI   =PHIS
         IBCND =IBCNDS
         CORINT=CORINS
         FBEAM =FBEAMS
         UMU0  =UMU0S
         PHI0  =PHI0S
         FISOT =FISOTS
         LAMBER=LAMBES
         ALBEDO=ALBEDS
         DELTAM=DELTAS
         ONLYFL=ONLYFS
         ACCUR =ACCURS
         PLANK =PLANKS
         WVNMLO=WVNMLS
         WVNMHI=WVNMHS
         BTEMP =BTEMPS
         TTEMP =TTEMPS
         TEMIS =TEMISS
         BPLANK=BPLNKS
         TPLANK=TPLNKS
         TEMPER(0)=TEMPES(0)
         TEMPER(1)=TEMPES(1)
         PKAG(0)=PKAGS(0)
         PKAG(1)=PKAGS(1)
         DO I=1, 7
            PRNT(I)=PRNTS(I)
         ENDDO
      END IF

      RETURN
      END
