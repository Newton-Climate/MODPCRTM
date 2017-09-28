      SUBROUTINE  UPISOT(ARRAY, CC, CMU, IPVT, MXCMU, NN, NSTR, OPRIM,  &
     &                    WK, XR0, XR1, Z0, Z1, ZPLK0, ZPLK1)

!       FINDS THE PARTICULAR SOLUTION OF THERMAL RADIATION OF SS(15)

!     ROUTINES CALLED:  SGECO, SGESL

!   I N P U T     V A R I A B L E S:

!       CC     :  CAPITAL-C-SUB-IJ IN EQ. SS(5)
!       CMU    :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE
!       OPRIM  :  DELTA-M SCALED SINGLE SCATTERING ALBEDO
!       XR0    :  EXPANSION OF THERMAL SOURCE FUNCTION
!       XR1    :  EXPANSION OF THERMAL SOURCE FUNCTION EQS. SS(14-16)
!       (REMAINDER ARE 'DISORT' INPUT VARIABLES)

!    O U T P U T    V A R I A B L E S:

!       Z0     :  SOLUTION VECTORS Z-SUB-ZERO OF EQ. SS(16)
!       Z1     :  SOLUTION VECTORS Z-SUB-ONE  OF EQ. SS(16)
!       ZPLK0, :  PERMANENT STORAGE FOR -Z0,Z1-, BUT RE-ORDERED
!        ZPLK1

!   I N T E R N A L    V A R I A B L E S:

!       ARRAY  :  COEFFICIENT MATRIX IN LEFT-HAND SIDE OF EQ. SS(16)
!       IPVT   :  INTEGER VECTOR OF PIVOT INDICES REQUIRED BY *LINPACK*
!       WK     :  SCRATCH ARRAY REQUIRED BY *LINPACK*
!+---------------------------------------------------------------------+

      INTEGER IPVT(*)
      REAL    ARRAY(MXCMU,*), CC(MXCMU,*), CMU(*), WK(*),               &
     &        Z0(*), Z1(*), ZPLK0(*), ZPLK1(*)

      DO 20 IQ=1, NSTR

         DO 10 JQ=1, NSTR
            ARRAY(IQ,JQ)=- CC(IQ,JQ)
10       CONTINUE
         ARRAY(IQ,IQ)=1.0 + ARRAY(IQ,IQ)

         Z1(IQ)=XR1
         Z0(IQ)=(1.-OPRIM) * XR0 + CMU(IQ) * Z1(IQ)
20    CONTINUE
!                       ** SOLVE LINEAR EQUATIONS: SAME AS IN *UPBEAM*,
!                       ** EXCEPT -ZJ- REPLACED BY -Z0-
      RCOND=0.0
      CALL  SGECO(ARRAY, MXCMU, NSTR, IPVT, RCOND, WK)
      TEST1=1.+RCOND
      IF(TEST1.EQ.1.)                                                   &
     &  CALL ERRMSG('UPISOT--SGECO SAYS MATRIX NEAR SINGULAR',.FALSE.)

      CALL  SGESL(ARRAY, MXCMU, NSTR, IPVT, Z0, 0)

      DO 30  IQ=1, NN
         ZPLK0(IQ+NN)  =Z0(IQ)
         ZPLK1(IQ+NN)  =Z1(IQ)
         ZPLK0(NN+1-IQ)=Z0(IQ+NN)
         ZPLK1(NN+1-IQ)=Z1(IQ+NN)
30    CONTINUE

      RETURN
      END
