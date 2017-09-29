      SUBROUTINE TERPEV(CWT,EVECC,GL,GU,MAZIM,                          &
     &  MXCMU,MXUMU,NN,NSTR,IU_MIN,NUMU,WK,YLMC,YLMU)

!     INTERPOLATE EIGENVECTORS TO USER ANGLES; EQ SD(8)
      IMPLICIT NONE

!     INPUT ARGUMENTS:
      INTEGER MAZIM,MXCMU,MXUMU,NN,NSTR,IU_MIN,NUMU
      REAL CWT(*),EVECC(MXCMU,*),GL(0:*),                               &
     &  YLMC(0:MXCMU,*),YLMU(0:MXCMU,0:*)

!     OUTPUT ARGUMENTS:
      REAL GU(0:MXUMU,*),WK(*)

!     LOCAL VARIABLES:
      INTEGER IQ,L,JQ,IU
      REAL SUMVAL

!     LOOP OVER STREAMS:
      DO IQ=1,NSTR
          DO L=MAZIM,NSTR-1
!                                       ** INNER SUM IN SD(8) TIMES ALL
!                                   ** FACTORS IN OUTER SUM BUT PLM(MU)
              SUMVAL=0.
              DO JQ=1,NSTR
                  SUMVAL=SUMVAL+CWT(JQ)*YLMC(L,JQ)*EVECC(JQ,IQ)
              ENDDO
              WK(L+1)=GL(L)*SUMVAL/2
          ENDDO
!                                    ** FINISH OUTER SUM IN SD(8)
!                                    ** AND STORE EIGENVECTORS
          DO IU=IU_MIN,NUMU
              SUMVAL=0.
              DO L=MAZIM,NSTR-1
                  SUMVAL=SUMVAL+WK(L+1)*YLMU(L,IU)
              ENDDO
              IF(IQ.LE.NN)THEN
                  GU(IU,IQ+NN)=SUMVAL
              ELSE
                  GU(IU,NSTR+1-IQ)=SUMVAL
              ENDIF
          ENDDO
      ENDDO
      RETURN
      END
