      SUBROUTINE  SETMTX(BDR, CBAND, CMU, CWT, DELM0, GC, KK, LAMBER,   &
     &                    LYRCUT, MI, MI9M2, MXCMU, NCOL, NCUT, NNLYRI, &
     &                    NN, NSTR, TAUCPR, WK)

!        CALCULATE COEFFICIENT MATRIX FOR THE SET OF EQUATIONS
!        OBTAINED FROM THE BOUNDARY CONDITIONS AND THE CONTINUITY-
!        OF-INTENSITY-AT-LAYER-INTERFACE EQUATIONS;  STORE IN THE
!        SPECIAL BANDED-MATRIX FORMAT REQUIRED BY LINPACK ROUTINES

!     ROUTINES CALLED:  ZEROIT

!     I N P U T      V A R I A B L E S:

!       BDR      :  SURFACE BIDIRECTIONAL REFLECTIVITY
!       CMU      :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE
!       CWT      :  WEIGHTS FOR GAUSS QUADRATURE OVER ANGLE COSINE
!       DELM0    :  KRONECKER DELTA, DELTA-SUB-M0
!       GC       :  EIGENVECTORS AT POLAR QUADRATURE ANGLES, SC(1)
!       KK       :  EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)
!       LYRCUT   :  LOGICAL FLAG FOR TRUNCATION OF COMPUT. LAYER
!       NN       :  NUMBER OF STREAMS IN A HEMISPHERE (NSTR/2)
!       NCUT     :  TOTAL NUMBER OF COMPUTATIONAL LAYERS CONSIDERED
!       TAUCPR   :  CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED)
!       (REMAINDER ARE 'DISORT' INPUT VARIABLES)

!   O U T P U T     V A R I A B L E S:

!       CBAND    :  LEFT-HAND SIDE MATRIX OF LINEAR SYSTEM EQ. SC(5),
!                   SCALED BY EQ. SC(12); IN BANDED FORM REQUIRED
!                   BY LINPACK SOLUTION ROUTINES
!       NCOL     :  COUNTS OF COLUMNS IN -CBAND-

!   I N T E R N A L    V A R I A B L E S:

!       IROW     :  POINTS TO ROW IN  -CBAND-
!       JCOL     :  POINTS TO POSITION IN LAYER BLOCK
!       LDA      :  ROW DIMENSION OF -CBAND-
!       NCD      :  NUMBER OF DIAGONALS BELOW OR ABOVE MAIN DIAGONAL
!       NCOL     :  COUNTS OF COLUMNS IN -CBAND-
!       NSHIFT   :  FOR POSITIONING NUMBER OF ROWS IN BAND STORAGE
!       WK       :  TEMPORARY STORAGE FOR 'EXP' EVALUATIONS
! ---------------------------------------------------------------------+
      LOGICAL LAMBER, LYRCUT
      REAL    BDR(MI,0:*), CBAND(MI9M2,NNLYRI), CMU(*), CWT(*),         &
     &        GC(MXCMU,MXCMU,*), KK(MXCMU,*), TAUCPR(0:*), WK(*)

      CALL  ZEROIT(CBAND(1,1), MI9M2*NCUT*NSTR)
      NCD   =3*NN - 1
      LDA   =3*NCD + 1
      NSHIFT=LDA - 2*NSTR + 1
      NCOL  =0
!                         ** USE CONTINUITY CONDITIONS OF EQ. STWJ(17)
!                         ** TO FORM COEFFICIENT MATRIX IN STWJ(20);
!                         ** EMPLOY SCALING TRANSFORMATION STWJ(22)
      DO 30  LC=1, NCUT

         DO 4  IQ=1, NN
            WK(IQ)=EXP(KK(IQ,LC) * (TAUCPR(LC) - TAUCPR(LC-1)))
 4       CONTINUE

         JCOL=0
         DO 10  IQ=1, NN
            NCOL=NCOL + 1
            IROW=NSHIFT - JCOL
            DO 5  JQ=1, NSTR
               CBAND(IROW+NSTR,NCOL)=  GC(JQ,IQ,LC)
               CBAND(IROW,     NCOL)=- GC(JQ,IQ,LC) * WK(IQ)
               IROW=IROW + 1
 5          CONTINUE
            JCOL=JCOL + 1
10       CONTINUE

         DO 20  IQ=NN+1, NSTR
            NCOL=NCOL + 1
            IROW=NSHIFT - JCOL
            DO 15  JQ=1, NSTR
               CBAND(IROW+NSTR,NCOL)=  GC(JQ,IQ,LC) * WK(NSTR+1-IQ)
               CBAND(IROW,     NCOL)=- GC(JQ,IQ,LC)
               IROW=IROW + 1
15          CONTINUE
            JCOL=JCOL + 1
20       CONTINUE

30    CONTINUE
!                  ** USE TOP BOUNDARY CONDITION OF STWJ(20A) FOR
!                  ** FIRST LAYER
      JCOL=0
      DO 40  IQ=1, NN
         EXPA=EXP(KK(IQ,1) * TAUCPR(1))
         IROW=NSHIFT - JCOL + NN
         DO 35  JQ=NN, 1, -1
            CBAND(IROW,JCOL+1)=GC(JQ,IQ,1) * EXPA
            IROW=IROW+1
35       CONTINUE
         JCOL=JCOL+1
40    CONTINUE

      DO 50  IQ=NN+1, NSTR
         IROW=NSHIFT - JCOL + NN
         DO 45  JQ=NN, 1, -1
            CBAND(IROW,JCOL+1)=GC(JQ,IQ,1)
            IROW=IROW+1
45       CONTINUE
         JCOL=JCOL+1
50    CONTINUE
!                           ** USE BOTTOM BOUNDARY CONDITION OF
!                           ** STWJ(20C) FOR LAST LAYER
      NNCOL=NCOL - NSTR
      JCOL =0
      DO 70  IQ=1, NN
         NNCOL=NNCOL + 1
         IROW =NSHIFT - JCOL + NSTR

         DO 60  JQ=NN+1, NSTR
            IF (LYRCUT .OR. (LAMBER .AND. DELM0.EQ.0)) THEN

!                          ** NO AZIMUTHAL-DEPENDENT INTENSITY IF LAM-
!                          ** BERT SURFACE; NO INTENSITY COMPONENT IF
!                          ** TRUNCATED BOTTOM LAYER

               CBAND(IROW,NNCOL)=GC(JQ,IQ,NCUT)
            ELSE
               SUMVAL=0.0
               DO 55  K=1, NN
                  SUMVAL=SUMVAL + CWT(K) * CMU(K) * BDR(JQ-NN,K)        &
     &                        * GC(NN+1-K,IQ,NCUT)
55             CONTINUE
               CBAND(IROW,NNCOL)=GC(JQ,IQ,NCUT) - (1.+DELM0) * SUMVAL
            END IF

            IROW=IROW + 1
60       CONTINUE
         JCOL=JCOL + 1
70    CONTINUE

      DO 90  IQ=NN+1, NSTR
         NNCOL=NNCOL + 1
         IROW =NSHIFT - JCOL + NSTR
         EXPA=WK(NSTR+1-IQ)

         DO 80  JQ=NN+1, NSTR

            IF (LYRCUT .OR. (LAMBER .AND. DELM0.EQ.0)) THEN
               CBAND(IROW,NNCOL)=GC(JQ,IQ,NCUT) * EXPA
            ELSE
               SUMVAL=0.0
               DO 75  K=1, NN
                  SUMVAL=SUMVAL + CWT(K) * CMU(K) * BDR(JQ-NN,K)        &
     &                        * GC(NN+1-K,IQ,NCUT)
75             CONTINUE
               CBAND(IROW,NNCOL)=(GC(JQ,IQ,NCUT)                        &
     &                               - (1.+DELM0) * SUMVAL) * EXPA
            END IF

            IROW=IROW + 1
80       CONTINUE
         JCOL=JCOL + 1
90    CONTINUE

      RETURN
      END

      SUBROUTINE ST0MTX(CBAND,GC,KK,MI9M2,MXCMU,                        &
     &  NCOL,NCUT,NNLYRI,NN,NSTR,TAUCPR,WK)

!     VERSION OF SETMTX TAILORED FOR ZERO SURFACE REFLECTIVITY.
!     CALCULATE COEFFICIENT MATRIX FOR THE SET OF EQUATIONS
!     OBTAINED FROM THE BOUNDARY CONDITIONS AND THE CONTINUITY-
!     OF-INTENSITY-AT-LAYER-INTERFACE EQUATIONS;  STORE IN THE
!     SPECIAL BANDED-MATRIX FORMAT REQUIRED BY LINPACK ROUTINES
      IMPLICIT NONE

!     ROUTINES CALLED:  ZEROIT

!     I N P U T      V A R I A B L E S:

!       GC       :  EIGENVECTORS AT POLAR QUADRATURE ANGLES, SC(1)
!       KK       :  EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)
!       NN       :  NUMBER OF STREAMS IN A HEMISPHERE (NSTR/2)
!       NCUT     :  TOTAL NUMBER OF COMPUTATIONAL LAYERS CONSIDERED
!       TAUCPR   :  CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED)
!       (REMAINDER ARE 'DISORT' INPUT VARIABLES)

!   O U T P U T     V A R I A B L E S:

!       CBAND    :  LEFT-HAND SIDE MATRIX OF LINEAR SYSTEM EQ. SC(5),
!                   SCALED BY EQ. SC(12); IN BANDED FORM REQUIRED
!                   BY LINPACK SOLUTION ROUTINES
!       NCOL     :  COUNTS OF COLUMNS IN -CBAND-

!   I N T E R N A L    V A R I A B L E S:

!       IROW     :  POINTS TO ROW IN  -CBAND-
!       JCOL     :  POINTS TO POSITION IN LAYER BLOCK
!       LDA      :  ROW DIMENSION OF -CBAND-
!       NCD      :  NUMBER OF DIAGONALS BELOW OR ABOVE MAIN DIAGONAL
!       NCOL     :  COUNTS OF COLUMNS IN -CBAND-
!       NSHIFT   :  FOR POSITIONING NUMBER OF ROWS IN BAND STORAGE
!       WK       :  TEMPORARY STORAGE FOR 'EXP' EVALUATIONS
! ---------------------------------------------------------------------+

!     ARGUMENTS:
      INTEGER MI9M2,MXCMU,NCOL,NCUT,NNLYRI,NN,NSTR
      REAL CBAND(MI9M2,NNLYRI),GC(MXCMU,MXCMU,*),KK(MXCMU,*),           &
     &  TAUCPR(0:*),WK(*)

!     LOCAL VARIABLES:
      INTEGER NCD,LDA,NSHIFT,LC,IQ,JCOL,IROW,JQ,NNCOL
      REAL EXPA

!     INITIALIZATIONS
      CALL ZEROIT(CBAND(1,1),MI9M2*NCUT*NSTR)
      NCD=3*NN-1
      LDA=3*NCD+1
      NSHIFT=LDA-2*NSTR+1
      NCOL=0

!     USE CONTINUITY CONDITIONS OF EQ. STWJ(17) TO FORM COEFFICIENT
!     MATRIX IN STWJ(20); EMPLOY SCALING TRANSFORMATION STWJ(22)
      DO LC=1,NCUT
          DO IQ=1,NN
              WK(IQ)=EXP(KK(IQ,LC)*(TAUCPR(LC)-TAUCPR(LC-1)))
          ENDDO
          JCOL=0
          DO IQ=1,NN
              NCOL=NCOL+1
              IROW=NSHIFT-JCOL
              DO JQ=1,NSTR
                  CBAND(IROW+NSTR,NCOL)=GC(JQ,IQ,LC)
                  CBAND(IROW,NCOL)=-GC(JQ,IQ,LC)*WK(IQ)
                  IROW=IROW+1
              ENDDO
              JCOL=JCOL+1
          ENDDO
          DO IQ=NN+1,NSTR
              NCOL=NCOL+1
              IROW=NSHIFT-JCOL
              DO JQ=1,NSTR
                  CBAND(IROW+NSTR,NCOL)=GC(JQ,IQ,LC)*WK(NSTR+1-IQ)
                  CBAND(IROW,NCOL)=-GC(JQ,IQ,LC)
                  IROW=IROW+1
              ENDDO
              JCOL=JCOL+1
          ENDDO
      ENDDO

!     USE TOP BOUNDARY CONDITION OF STWJ(20A) FOR FIRST LAYER
      JCOL=0
      DO IQ=1,NN
          EXPA=EXP(KK(IQ,1)*TAUCPR(1))
          IROW=NSHIFT-JCOL+NN
          DO JQ=NN,1,-1
              CBAND(IROW,JCOL+1)=GC(JQ,IQ,1)*EXPA
              IROW=IROW+1
          ENDDO
          JCOL=JCOL+1
      ENDDO
      DO IQ=NN+1,NSTR
          IROW=NSHIFT-JCOL+NN
          DO JQ=NN,1,-1
              CBAND(IROW,JCOL+1)=GC(JQ,IQ,1)
              IROW=IROW+1
          ENDDO
          JCOL=JCOL+1
      ENDDO

!     USE BOTTOM BOUNDARY CONDITION OF STWJ(20C) FOR LAST LAYER
      NNCOL=NCOL-NSTR
      JCOL=0
      DO IQ=1,NN
          NNCOL=NNCOL+1
          IROW=NSHIFT-JCOL+NSTR
          DO JQ=NN+1,NSTR
              CBAND(IROW,NNCOL)=GC(JQ,IQ,NCUT)
              IROW=IROW+1
          ENDDO
          JCOL=JCOL+1
      ENDDO
      DO IQ=NN+1,NSTR
          NNCOL=NNCOL+1
          IROW=NSHIFT-JCOL+NSTR
          EXPA=WK(NSTR+1-IQ)
          DO JQ=NN+1,NSTR
              CBAND(IROW,NNCOL)=GC(JQ,IQ,NCUT)*EXPA
              IROW=IROW+1
          ENDDO
          JCOL=JCOL+1
      ENDDO
      RETURN
      END
