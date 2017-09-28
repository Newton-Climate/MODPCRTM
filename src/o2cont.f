      SUBROUTINE O2CONT(V,SIGMA,ALPHA,BETA)

!     THIS ROUTINE IS DRIVEN BY FREQUENCY, RETURNING ONLY THE
!     O2 COEFFICIENTS, INDEPENDENT OF TEMPERATURE.

!  *******************************************************************
!  *  THESE COMMENTS APPLY TO THE COLUMN ARRAYS FOR:                 *
!  *       PBAR*UBAR(O2)                                             *
!  *       PBAR*UBAR(O2)*DT                                          *
!  *   AND PBAR*UBAR(O2)*DT*DT    WHERE:  DT=TBAR-220.               *
!  *  THAT HAVE BEEN COMPILED IN OTHER PARTS OF THE LOWTRAN CODE     *
!  *                                                                 *
!  *  LOWTRAN7 COMPATIBLE:                                           *
!  *  O2 CONTINUUM SUBROUTINE FOR 1395-1760CM-1                      *
!  *  MODIFIED BY G.P. ANDERSON, APRIL '88                           *
!  *                                                                 *
!  *  THE EXPONENTIAL TEMPERATURE EMPLOYED IN THE FASCOD2 ALGORITHM  *
!  *  (SEE BELOW) IS NOT READILY SUITABLE FOR LOWTRAN.  THEREFORE    *
!  *  THE EXPONENTIALS HAVE BEEN LINEARLY EXPANDED, KEEPING ONLY THE *
!  *  LINEAR AND QUADRATIC TERMS:                                    *
!  *                                                                 *
!  *  EXP(A*DT)=1.+ A*DT + (A*DT)**2/2. + ....                       *
!  *                                                                 *
!  *     EXP(B*DT*DT)=1.+ B*DT*DT + (B*DT*DT)**2/2. + ....           *
!  *                                                                 *
!  *  THE PRODUCT OF THE TWO TERMS IS:                               *
!  *                                                                 *
!  *     (1. + A*DT + (A*A/2. + B)*DT*DT )                           *
!  *                                                                 *
!  *  THIS EXPANSION ONLY WORKS WELL FOR SMALL VALUES OF X IN EXP(X) *
!  *                                                                 *
!  *  SINCE DT = T-220., THE APPROXIMATION IS VERY GOOD UNTIL        *
!  *  T.GT.260. OR DT.GT.40.   AT T=280, THE MAXIMUM ERRORS ARE STILL*
!  *  LESS THAN 10% BUT AT T=300, THOSE ERRORS ARE AS LARGE AS 20%   *
!  *******************************************************************

!     THE FOLLOWING COMMENTS ARE EXCERPTED DIRECTLY FROM FASCOD2

!      THIS SUBROUTINE CONTAINS THE ROGERS AND WALSHAW
!      EQUIVALENT COEFFICIENTS DERIVED FROM THE THEORETICAL
!      VALUES SUPPLIED BY ROLAND DRAYSON. THESE VALUES USE
!      THE SAME DATA AS TIMOFEYEV AND AGREE WITH TIMOFEYEV'S RESULTS.
!      THE DATA ARE IN THE FORM OF STRENGTHS(O2SO) AND TWO
!      COEFFICIENTS (O2A & O2B),  WHICH ARE USED TO CORRECT FOR
!      TEMPERATURE. THE DEPENDENCY ON PRESSURE SQUARED
!      IS CONTAINED IN THE P*WO2 PART OF THE CONSTANT.
!      NOTE THAT SINCE THE COEFFICIENTS ARE FOR AIR, THE
!      THE STRENGTHS ARE DIVIDED BY THE O2 MIXING RATIO FOR
!      DRY AIR OF 0.20946 (THIS IS ASSUMED CONSTANT).
!      ORIGINAL FORMULATION OF THE COEFFICIENTS WAS BY LARRY GORDLEY.
!      THIS VERSION WRITTEN BY EARL THOMPSON, JULY 1984.
      IMPLICIT NONE

!     PARAMETERS:
!       O2FAC    RINSLAND FACTOR [JGR 94: 16,303 - 16,322 (1989)].
      REAL O2FAC
      PARAMETER(O2FAC=.78)

!     INPUT ARGUMENTS:
!       V        SPECTRAL FREQUENCY [CM-1].
      REAL V

!     OUTPUT ARGUMENTS:
      REAL SIGMA,ALPHA,BETA

!     LOCAL VARIABLES:
!       IND      SPECTRAL INDEX.
!       FRAC     SPECTRAL INTERPOLATION FRACTION.
      INTEGER IND
      REAL FRAC

!     DATA:
      REAL O2S0(74),O2A(74),O2B(74)
      DATA O2S0/      .0, .110E-8, .220E-8, .440E-8, .881E-8, .176E-7,  &
     &  .353E-7, .705E-7, .141E-6, .158E-6, .174E-6, .190E-6, .207E-6,  &
     &  .253E-6, .307E-6, .357E-6, .401E-6, .445E-6, .508E-6, .570E-6,  &
     &  .599E-6, .627E-6, .650E-6, .672E-6, .763E-6, .873E-6, .101E-5,  &
     &  .109E-5, .121E-5, .133E-5, .139E-5, .145E-5, .148E-5, .140E-5,  &
     &  .134E-5, .126E-5, .118E-5, .114E-5, .109E-5, .105E-5, .105E-5,  &
     &  .105E-5, .104E-5, .103E-5, .992E-6, .945E-6, .876E-6, .806E-6,  &
     &  .766E-6, .726E-6, .640E-6, .555E-6, .469E-6, .416E-6, .364E-6,  &
     &  .311E-6, .266E-6, .222E-6, .177E-6, .170E-6, .162E-6, .155E-6,  &
     &  .143E-6, .130E-6, .118E-6, .905E-7, .629E-7, .316E-7, .157E-7,  &
     &  .786E-8, .393E-8, .196E-8, .982E-9,      .0/
      DATA O2A /      .0, .147E-3, .147E-3, .147E-3, .147E-3, .147E-3,  &
     &  .147E-3, .147E-3, .147E-3, .122E-2, .204E-2, .217E-2, .226E-2,  &
     &  .126E-2, .362E-3,-.198E-2,-.545E-2,-.786E-2,-.624E-2,-.475E-2,  &
     & -.506E-2,-.533E-2,-.586E-2,-.635E-2,-.644E-2,-.679E-2,-.741E-2,  &
     & -.769E-2,-.780E-2,-.788E-2,-.844E-2,-.894E-2,-.899E-2,-.922E-2,  &
     & -.892E-2,-.857E-2,-.839E-2,-.854E-2,-.871E-2,-.889E-2,-.856E-2,  &
     & -.823E-2,-.796E-2,-.768E-2,-.715E-2,-.638E-2,-.570E-2,-.491E-2,  &
     & -.468E-2,-.443E-2,-.333E-2,-.184E-2, .313E-3,-.164E-4,-.417E-3,  &
     & -.916E-3,-.206E-2,-.343E-2,-.515E-2,-.365E-2,-.172E-2, .926E-3,  &
     &  .168E-2, .262E-2, .380E-2, .551E-2, .889E-2, .889E-2, .889E-2,  &
     &  .889E-2, .889E-2, .889E-2, .889E-2,      .0/
      DATA O2B  /     .0, .306E-4,-.306E-4,-.306E-4,-.306E-4,-.306E-4,  &
     & -.306E-4,-.306E-4,-.306E-4,-.218E-4,-.159E-4,-.346E-5, .642E-5,  &
     &  .360E-5,-.140E-5, .157E-4, .471E-4, .656E-4, .303E-4,-.192E-5,  &
     &  .705E-5, .149E-4, .200E-4, .245E-4, .158E-4, .841E-5, .201E-5,  &
     &  .555E-5, .108E-4, .150E-4, .193E-4, .230E-4, .243E-4, .226E-4,  &
     &  .184E-4, .157E-4, .169E-4, .197E-4, .226E-4, .258E-4, .235E-4,  &
     &  .212E-4, .185E-4, .156E-4, .125E-4, .872E-5, .760E-5, .577E-5,  &
     &  .334E-7,-.652E-5,-.977E-5,-.157E-4,-.273E-4,-.180E-4,-.641E-5,  &
     &  .817E-5, .326E-4, .626E-4, .101E-3, .755E-4, .430E-4,-.113E-5,  &
     & -.578E-5,-.120E-4,-.208E-4,-.235E-4,-.364E-4, .364E-4,-.364E-4,  &
     & -.364E-4,-.364E-4,-.364E-4,-.364E-4,      .0/
      SAVE O2S0,O2A,O2B

!     LOWTRAN7 FORMULATION (MODELED AFTER THE LOWTRAN UV-O3 BANDS)
      IF(V.LE.1395. .OR. V.GE.1760.)THEN
          SIGMA=0.
          ALPHA=0.
          BETA=0.
          RETURN
      ENDIF
      FRAC=(V-1390)/5
      IND=INT(FRAC)
      FRAC=FRAC-IND
      SIGMA=O2S0(IND)
      ALPHA=O2A(IND)
      BETA=O2B(IND)
      IND=IND+1
      SIGMA=SIGMA+FRAC*(O2S0(IND)-SIGMA)
      ALPHA=ALPHA+FRAC*(O2A(IND) -ALPHA)
      BETA =BETA +FRAC*(O2B(IND) -BETA )

!     OLD 'FASCOD2' TEMPERATURE DEPENDENCE USING BLOCK DATA ARRAYS

!     C(J)=O2S0(I)* EXP(O2A(I)*TD+O2B(I)*TD*TD) /(0.20946*VJ)

!     NEW COEFFICIENT DEFINITIONS FOR LOWTRAN FORMULATION

      BETA=ALPHA**2/2+BETA
      SIGMA=SIGMA/.20946

!     NEW 'LOWTRAN7' TEMPERATURE DEPENDENCE

!     THIS WOULD BE THE CODING FOR THE LOWTRAN7 FORMULATION, BUT
!       BECAUSE THE T-DEPENDENCE IS INCLUDED IN THE AMOUNTS, ONLY
!       THE COEFFICIENTS (SIGMA, ALPHA & BETA) ARE BEING RETURNED

!     C(J)=SIGMA*(1.+ALPHA*TD+BETA*TD*TD)

!     THE COEFFICIENTS FOR O2 HAVE BEEN MULTIPLIED BY A FACTOR
!     OF 0.78  [RINSLAND ET AL, 1989: JGR 94; 16,303 - 16,322.].
      SIGMA=O2FAC*SIGMA
      RETURN
      END
