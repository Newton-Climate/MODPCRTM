      REAL FUNCTION  XIFUNC(UMU1, UMU2, UMU3, TAU)

!                                            Xi function of EQ. STW (72)

!                    I N P U T   V A R I A B L E S

!        TAU         optical thickness of the layer
!        UMU1,2,3    cosine of zenith angle_1, _2, _3

      X1=1./UMU1 - 1./UMU2
      X2=1./UMU1 - 1./UMU3

      IF ((UMU2.EQ.UMU3) .AND. (UMU1.EQ.UMU2))  THEN
         XIFUNC=TAU**2*EXP(-TAU/UMU1) / (2.*UMU1*UMU2)

      ELSE IF ((UMU2.EQ.UMU3) .AND. (UMU1.NE.UMU2))  THEN
         XIFUNC=((TAU-1./X1)*EXP(-TAU/UMU2) + EXP(-TAU/UMU1)/X1) /      &
     &            (X1*UMU1*UMU2)

      ELSE IF ((UMU2.NE.UMU3) .AND. (UMU1.EQ.UMU2))  THEN
         XIFUNC=((EXP(-TAU/UMU3) - EXP(-TAU/UMU1))/X2                   &
     &              - TAU*EXP(-TAU/UMU1)) / (X2*UMU1*UMU2)

      ELSE IF ((UMU2.NE.UMU3) .AND. (UMU1.EQ.UMU3))  THEN
         XIFUNC=((EXP(-TAU/UMU2) - EXP(-TAU/UMU1))/X1                   &
     &              - TAU*EXP(-TAU/UMU1)) / (X1*UMU1*UMU2)

      ELSE
         XIFUNC=((EXP(-TAU/UMU3) - EXP(-TAU/UMU1))/X2 -                 &
     &              (EXP(-TAU/UMU2) - EXP(-TAU/UMU1))/X1) /             &
     &            (X2*UMU1*UMU2)
      END IF

      RETURN
      END
