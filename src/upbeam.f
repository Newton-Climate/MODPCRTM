      SUBROUTINE  UPBEAM(ARRAY, CC, CMU, DELM0, FBEAM, GL, IPVT,MAZIM,  &
     &  MXCMU, NN, NSTR, PI, UMU0, WK, YLM0, YLMC, ZJ, ZZ)

!         FINDS THE INCIDENT-BEAM PARTICULAR SOLUTION  OF SS(18)
!     DEC 2004:  Added dithering of cosine solar zenith if numerical
!                instability occurs.  Also toughened stability test.

!     ROUTINES CALLED:  SGECO, SGESL

!   I N P U T    V A R I A B L E S:

!       CC     :  CAPITAL-C-SUB-IJ IN EQ. SS(5)
!       CMU    :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE
!       DELM0  :  KRONECKER DELTA, DELTA-SUB-M0
!       GL     :  DELTA-M SCALED LEGENDRE COEFFICIENTS OF PHASE FUNCTION
!                    (INCLUDING FACTORS 2L+1 AND SINGLE-SCATTER ALBEDO)
!       MAZIM  :  ORDER OF AZIMUTHAL COMPONENT
!       YLM0   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL
!                 AT THE BEAM ANGLE
!       YLMC   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL
!                 AT THE QUADRATURE ANGLES
!       (REMAINDER ARE 'DISORT' INPUT VARIABLES)

!   O U T P U T    V A R I A B L E S:

!       ZJ     :  RIGHT-HAND SIDE VECTOR CAPITAL-X-SUB-ZERO IN SS(19);
!                 ALSO THE SOLUTION VECTOR CAPITAL-Z-SUB-ZERO
!                 AFTER SOLVING THAT SYSTEM
!       ZZ     :  PERMANENT STORAGE FOR -ZJ-, BUT RE-ORDERED

!   I N T E R N A L    V A R I A B L E S:

!       UMU0SV :  SAVED VALUE OF UMU0 IF DITHERING IS REQUIRED.
!       ARRAY  :  COEFFICIENT MATRIX IN LEFT-HAND SIDE OF EQ. SS(19)
!       IPVT   :  INTEGER VECTOR OF PIVOT INDICES REQUIRED BY *LINPACK*
!       WK     :  SCRATCH ARRAY REQUIRED BY *LINPACK*
!+---------------------------------------------------------------------+

      INTEGER  IPVT(*)
      REAL     UMU0SV, ARRAY(MXCMU,*), CC(MXCMU,*), CMU(*), GL(0:*),    &
     &         WK(*), YLM0(0:*), YLMC(0:MXCMU,*), ZJ(*), ZZ(*)

!     INITIALIZE UMU0SV TO ZERO TO INDICATE NO DITHERING (YET).
      UMU0SV=0.

!     REENTRY POINT IF DITHERING IS REQUIRED:
   10 CONTINUE
      DO IQ=1, NSTR

          DO JQ=1, NSTR
              ARRAY(IQ,JQ)=-CC(IQ,JQ)
          ENDDO
          ARRAY(IQ,IQ)=1 + CMU(IQ) / UMU0 + ARRAY(IQ,IQ)

          SUMVAL=0.
          DO K=MAZIM, NSTR-1
              SUMVAL=SUMVAL + GL(K) * YLMC(K,IQ) * YLM0(K)
          ENDDO
          ZJ(IQ)=(2 - DELM0) * FBEAM * SUMVAL / (4*PI)
      ENDDO
!                  ** FIND L-U (LOWER/UPPER TRIANGULAR) DECOMPOSITION
!                  ** OF -ARRAY- AND SEE IF IT IS NEARLY SINGULAR
!                  ** (NOTE:  -ARRAY- IS DESTROYED)
      RCOND=0.
      CALL SGECO(ARRAY, MXCMU, NSTR, IPVT, RCOND, WK)
      TEST1=1.+RCOND/8
      IF(TEST1.EQ.1.)THEN
          IF(UMU0SV.LE.0.)THEN

!             DITHER COSINE SOLAR ZENITH:
              CALL ERRMSG('UPBEAM--TEMPORARILY REDUCED COSINE SOLAR'//  &
     &          ' ZENITH BY 0.01% TO AVOID SINGULAR MATRIX',.FALSE.)
              UMU0SV=UMU0
              UMU0=.9999*UMU0
              GOTO 10
          ENDIF

!         DITHER OF COSINE SOLAR ZENITH FAILED:
          CALL ERRMSG('UPBEAM--SGECO SAYS MATRIX REMAINS NEAR'//        &
     &      ' SINGULAR AFTER DITHER OF COSINE SOLAR ZENITH',.TRUE.)
      ENDIF
!     IF(UMU0SV.GT.0.)UMU0=UMU0SV

!                ** SOLVE LINEAR SYSTEM WITH COEFF MATRIX -ARRAY-
!                ** (ASSUMED ALREADY L-U DECOMPOSED) AND R.H. SIDE(S)
!                ** -ZJ-;  RETURN SOLUTION(S) IN -ZJ-
      JOB=0
      CALL  SGESL(ARRAY, MXCMU, NSTR, IPVT, ZJ, JOB)

      DO IQ=1, NN
          ZZ(IQ+NN)  =ZJ(IQ)
          ZZ(NN+1-IQ)=ZJ(IQ+NN)
      ENDDO

      RETURN
      END
