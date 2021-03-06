      REAL FUNCTION DENLAY(DENBOT,DENTOP)

!     FUNCTION DENLAY RETURNS AN EXPONENTIALLY AVERAGED VERTICAL
!     PATH LAYER DENSITY GIVEN LAYER BOUNDARY DENSITY VALUES:
!                          1
!                        /                        X
!     DENLAY  =  DENBOT  |    ( DENTOP / DENBOT )     DX
!                        /
!                          0
!     THE MEAN VALUE IS RETURNED IF DENBOT OR DENTOP ARE NOT POSITIVE.
      IMPLICIT NONE

!     INPUT ARGUMENTS:
!       DENBOT   BOTTOM LAYER BOUNDARY DENSITY.
!       DENTOP   TOP LAYER BOUNDARY DENSITY.
      REAL DENBOT,DENTOP

!     LOCAL VARIABLES:
!       DENRAT   RATIO OF DENTOP TO DENBOT.
!       EPSILN   THE VARIATION OF DENRAT FROM ONE.
      REAL DENRAT,EPSILN
      IF(DENBOT.LE.0. .OR. DENTOP.LE.0.)THEN

!         ARITHMETIC MEAN:
          DENLAY=(DENBOT+DENTOP)/2
      ELSE
          DENRAT=DENTOP/DENBOT
          EPSILN=DENRAT-1
          IF(ABS(EPSILN).GT..01)THEN

!             ANALYTIC EXPRESSION:
              DENLAY=(DENTOP-DENBOT)/LOG(DENRAT)
          ELSE

!             TAYLOR EXPANSION:
              DENLAY=DENBOT/(1+EPSILN*(.5+EPSILN*(.3333333+EPSILN/4)))
          ENDIF
      ENDIF
      RETURN
      END
