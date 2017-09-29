      SUBROUTINE  SOLVE1(B, CBAND, FISOT, IHOM, IPVT, LL, MI9M2,MXCMU,  &
     &                    NCOL, NCUT, NN, NNLYRI, NSTR, Z, LSGBCO)

!        CONSTRUCT RIGHT-HAND SIDE VECTOR -B- FOR ISOTROPIC INCIDENCE
!        (ONLY) ON EITHER TOP OR BOTTOM BOUNDARY AND SOLVE SYSTEM
!        OF EQUATIONS OBTAINED FROM THE BOUNDARY CONDITIONS AND THE
!        CONTINUITY-OF-INTENSITY-AT-LAYER-INTERFACE EQUATIONS

!     ROUTINES CALLED:  SGBCO, SGBSL, ZEROIT

!     I N P U T      V A R I A B L E S:

!       CBAND    :  LEFT-HAND SIDE MATRIX OF LINEAR SYSTEM EQ. SC(5),
!                   SCALED BY EQ. SC(12); IN BANDED FORM REQUIRED
!                   BY LINPACK SOLUTION ROUTINES
!       IHOM     :  DIRECTION OF ILLUMINATION FLAG
!       NCOL     :  COUNTS OF COLUMNS IN -CBAND-
!       NN       :  ORDER OF DOUBLE-GAUSS QUADRATURE (NSTR/2)
!       LSGBCO   :  IF .TRUE., CALL SGBCO FOR ILLUMINATION FROM BELOW.
!       (REMAINDER ARE 'DISORT' INPUT VARIABLES)

!    O U T P U T     V A R I A B L E S:

!       B        :  RIGHT-HAND SIDE VECTOR OF EQ. SC(5) GOING INTO
!                   *SGBSL*; RETURNS AS SOLUTION VECTOR OF EQ.
!                   SC(12), CONSTANTS OF INTEGRATION WITHOUT
!                   EXPONENTIAL TERM
!       LL      :   PERMANENT STORAGE FOR -B-, BUT RE-ORDERED

!   I N T E R N A L    V A R I A B L E S:

!       IPVT     :  INTEGER VECTOR OF PIVOT INDICES
!       NCD      :  NUMBER OF DIAGONALS BELOW OR ABOVE MAIN DIAGONAL
!       RCOND    :  INDICATOR OF SINGULARITY FOR -CBAND-
!       Z        :  SCRATCH ARRAY REQUIRED BY *SGBCO*
!----------------------------------------------------------------------+
      INTEGER  IPVT(*)
      REAL     B(NNLYRI), CBAND(MI9M2,NNLYRI), LL(MXCMU,*), Z(*)
      LOGICAL LSGBCO

      CALL  ZEROIT(B(1), NNLYRI)
      NCD=3*NN - 1

      IF (IHOM.EQ.1)  THEN
!                             ** BECAUSE THERE ARE NO BEAM OR EMISSION
!                             ** SOURCES, REMAINDER OF -B- ARRAY IS ZERO
         DO I=1, NN
            B(I)=FISOT
            B(NCOL-NN+I)=0.0
         ENDDO

         RCOND=0.0
         CALL  SGBCO(CBAND(1,1), MI9M2, NCOL, NCD, NCD, IPVT, RCOND, Z)
         TEST1=1.+RCOND
         IF(TEST1.EQ.1.)CALL ERRMSG                                     &
     &         ('SOLVE1--SGBCO SAYS MATRIX NEAR SINGULAR', .FALSE.)

      ELSE IF (IHOM.EQ.2)  THEN

         IF(LSGBCO)THEN
             RCOND=0.
             CALL SGBCO(CBAND(1,1),MI9M2,NCOL,NCD,NCD,IPVT,RCOND,Z)
             TEST1=1.+RCOND
             IF(TEST1.EQ.1.)CALL ERRMSG                                 &
     &         ('SOLVE1--SGBCO SAYS MATRIX NEAR SINGULAR', .FALSE.)
         ENDIF
         DO I=1, NN
            B(I)=0.0
            B(NCOL-NN+I)=FISOT
         ENDDO

      END IF

      CALL  SGBSL(CBAND(1,1), MI9M2, NCOL, NCD, NCD, IPVT, B(1), 0)

!                          ** ZERO -CBAND- TO GET RID OF 'FOREIGN'
!                          ** ELEMENTS PUT IN BY LINPACK
      DO LC=1, NCUT
         IPNT=LC*NSTR - NN
         DO IQ=1, NN
            LL(NN+1-IQ, LC)=B(IPNT+1-IQ)
            LL(IQ+NN,   LC)=B(IQ+IPNT)
         ENDDO
      ENDDO

      RETURN
      END
