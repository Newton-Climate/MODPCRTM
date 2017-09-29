      SUBROUTINE BS(I,A,B,N,S)

!     THIS SUBROUTINE DOES THE BINARY SEARCH FOR THE INDEX I
!     SUCH THAT A IS IN BETWEEN B(I) AND B(I+1)
!     AND CALCULATES THE INTERPOLATION PARAMETER S
!     SUCH THAT A=S*B(I+1)+(1.-S)*B(I)
      REAL B(*),S,A
      INTEGER I,N,M,J
      I=1
      J=N
   10 CONTINUE
      M=(I+J)/2
      IF(A.LE.B(M)) THEN
          J=M
      ELSE
          I=M
      END IF
      IF(J.GT.I+1)GOTO 10
      S=(A-B(I))/(B(I+1)-B(I))
      RETURN
      END
