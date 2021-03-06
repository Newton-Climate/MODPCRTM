      DOUBLE PRECISION FUNCTION RTBIS(REARTH,X1,X2,CPATH)

!     THIS FUNCTION FINDS THE TANGENT HEIGHT BETWEEN X1 AND X2 AS THE
!     ROOT OF FUNC(X) = X * (REFRACTIVE INDEX) - CPATH.  THIS ROUTINE
!     IS FROM "NUMERICAL RECIPES: THE ART OF SCIENTIFIC COMPUTING",
!     W.H. PRESS, B.P. FLANNERY, S.A. TEUKOLSKY AND W.T. VETTERING,
!     (CAMBRIDGE UNIVERSITY PRESS, CAMBRIDGE, 1986).
      IMPLICIT NONE

!     PARAMETERS:
      INTEGER JMX
      DOUBLE PRECISION XACC
      PARAMETER(JMX=40,XACC=1.D-5)
      INCLUDE 'PARAMS.h'

!     ARGUMENTS:
!       REARTH   RADIUS OF THE EARTH [KM].
      DOUBLE PRECISION REARTH,X1,X2,CPATH

!     LOCAL VARIABLES:
      DOUBLE PRECISION F,FMID,DX,XMID,F1,F2,RXRAT
      INTEGER J
      CALL IRFXN(X2,FMID,RXRAT)
      FMID=FMID*(X2+REARTH)-CPATH
      CALL IRFXN(X1,F,RXRAT)
      F=F*(X1+REARTH)-CPATH
      IF(F*FMID.GT.0.)THEN
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP 'Error in RTBIS:  Root must be bracketed for bisection.'
      ENDIF
      IF(F.LT.0.)THEN
          RTBIS=X1
          DX=X2-X1
      ELSE
          RTBIS=X2
          DX=X1-X2
      ENDIF
      DO J=1,JMX
          DX=DX/2
          XMID=RTBIS+DX
          CALL IRFXN(XMID,FMID,RXRAT)
          FMID=FMID*(XMID+REARTH)-CPATH
          IF(FMID.LE.0.)RTBIS=XMID
          IF(ABS(DX).LT.XACC .OR. FMID.EQ.0.)RETURN
      ENDDO

!     COMES HERE IF UNABLE TO SOLVE.
      CALL IRFXN(X2,F2,RXRAT)
      F2=F2*(X2+REARTH)-CPATH
      CALL IRFXN(X1,F1,RXRAT)
      F1=F1*(X1+REARTH)-CPATH
      IF(ABS(F2).LT.ABS(F1))THEN
          RTBIS=X2
      ELSE
          RTBIS=X1
      ENDIF
      END
