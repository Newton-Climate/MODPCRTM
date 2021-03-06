      REAL FUNCTION BETABS(COSANG,G)

!     FUNCTION BETABS SUPPLIES THE BACK SCATTER FRACTION FOR A
!     GIVEN ASYMMETRY FACTOR AND COSINE OF ANGLE (A.E.R. 1986).

!                         1
!                      1  /
!     BETA(COSANG)  =  -  |  P (-M,COSANG) dM
!                      2  /   G
!                         0

!     WHERE P (-M,COSANG) IS THE HENYEY-GREENSTEIN PHASE FUNCTION.
!            G
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     INPUTS:
      REAL COSANG,G

!     COMMONS:
      INCLUDE 'IFIL.h'

!     DECLARE BLOCK DATA ROUTINE EXTERNAL:
      EXTERNAL DEVCBD

!     DECLARE LOCAL VARIABLES:
      INTEGER N
      REAL ABSCOS,ABSG,COSN,BMAX,BMIN,COSINE(9),A(10,5)
      DATA COSINE/0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0/
      DATA A/10*.5,                                                     &
     &   0.0,         0.13979293, -0.12019000, -0.46017123, -0.40682880,&
     &  -0.30015417, -0.55347441, -0.62679466, -0.84678101, -0.40682368,&
     &   0.0,        -1.5989874,  -0.27242199,  1.1874739,   0.49409051,&
     &  -0.35928947,  0.37397958,  0.18057741,  0.50718036,  0.01406832,&
     &   0.0,         3.5184117,   0.63859601, -2.0812308,  -1.0144699, &
     &   0.15894758, -0.74761783, -0.37416958, -0.37404011,  0.10556077,&
     &   0.0,        -2.5592172,  -0.74598401,  0.85392041,  0.42720824,&
     &   0.00049606,  0.42711267,  0.32038684,  0.21367466, -0.21280542/

!     SPECIAL CASES
      ABSCOS=ABS(COSANG)
      ABSG=ABS(G)
      IF(ABSCOS.LT..000001 .OR. ABSG.LT..000001)THEN
          BETABS=.5
          RETURN
      ELSEIF(ABSG.GT..999999)THEN
          BETABS=0.
          IF(COSANG.LT.0.)BETABS=1.
          IF(G.LT.0.)BETABS=1.-BETABS
          RETURN
      ENDIF

!     BACKSCATTERING INTERPOLATION
      IF(ABSCOS.LT..7)THEN
          N=IFIX(10*ABSCOS+1)
      ELSEIF(ABSCOS.LT..8)THEN
          N=7
      ELSEIF(ABSCOS.LE.1.)THEN
          N=8
      ELSE
          WRITE(IPR,'(A,F15.7)')' Error in BETABS:  COSANG =',COSANG
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP ' Error in BETABS'
      ENDIF
      COSN=COSINE(N)
      BMAX=A(N,1)+ABSG*(A(N,2)+ABSG*(A(N,3)+ABSG*(A(N,4)+ABSG*A(N,5))))
      N=N+1
      BMIN=A(N,1)+ABSG*(A(N,2)+ABSG*(A(N,3)+ABSG*(A(N,4)+ABSG*A(N,5))))
      BETABS=BMAX+(BMIN-BMAX)*(ABSCOS-COSN)/(COSINE(N)-COSN)

!     CHECK FOR ROUND-OFF ERROR PROBLEMS
      IF(BETABS.GE..5)THEN
          BETABS=.5
          RETURN
      ENDIF
      IF(BETABS.LT.0.)BETABS=0.

!     IF G IS NEGATIVE, THE COMPLEMENT IS REQUIRED.
      IF(G.LT.0.)BETABS=1.-BETABS
      IF(COSANG.LT.0.)BETABS=1.-BETABS
      RETURN
      END
