      REAL FUNCTION RATIO(A,B)

!        CALCULATE RATIO A/B WITH OVER- AND UNDER-FLOW PROTECTION
         IF(ABS(A).LT.1.E-8 .AND. ABS(B).LT.1.E-8)THEN
             RATIO=1.
         ELSEIF(B.EQ.0.)THEN
             RATIO=1.E+20
         ELSE
             RATIO=A/B
         ENDIF

      RETURN
      END
