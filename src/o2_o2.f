      REAL FUNCTION O2_O2(V)

!     RETURNS O2_O2 SPECTRAL ABSORBANCE BETWEEN 334.8 AND 666.7NM FOR
!     A 89.5 CM PATH AT 296K AND 55 ATM.  THE 10 CM-1 SPACED DATA ARE
!     FROM FIGURE 1 OF "Absorption Measurements of Oxygen Between 330
!     and 1140 nm," G.D. Greenblatt, J.J. Orlando, J.B. Burkholder,
!     and A.R. Ravishabkara, J. Geophys. Res., 95, 18577-18582, (1990).
      IMPLICIT NONE

!     INPUT ARGUMENTS:
!       V        SPECTRAL FREQUENCY [CM-1].
      REAL V

!     COMMONS:
!       O2VIS    O2_O2 ABSORBANCE FOR A 89.5 CM PATH AT 296K & 55ATM.
      REAL O2VIS
      COMMON/O2_VIS/O2VIS(1500:2988)

!     LOCAL VARIABLES:
!       FV10     DATA INTERPOLATION VARIABLE.
!       IV10     DATA INTERPOLATION LOWER INDEX.
      REAL FV10
      INTEGER IV10

      IF(V.LT.15140. .OR. V.GT.29870.)THEN
          O2_O2=0.
      ELSE
          FV10=V/10
          IV10=INT(FV10)
          O2_O2=O2VIS(IV10)
          O2_O2=O2_O2+(FV10-IV10)*(O2VIS(IV10+1)-O2_O2)
      ENDIF

!     RETURN TO FRQ5DT
      RETURN
      END

      BLOCK DATA BO2VIS

!     O2_O2 SPECTRAL ABSORBANCE BETWEEN 334.8 AND 666.7NM FOR A
!     89.5 CM PATH AT 296K AND 55 ATM.  THE 10 CM-1 SPACED DATA ARE
!     FROM FIGURE 1 OF "Absorption Measurements of Oxygen Between 330
!     and 1140 nm," G.D. Greenblatt, J.J. Orlando, J.B. Burkholder, and
!     A.R. Ravishabkara, J. Geophys. Res., 95, 18577-18582, (1990).

!     COMMONS:
!       O2VIS    O2_O2 ABSORBANCE FOR A 89.5 CM PATH AT 296K & 55ATM.
      REAL O2VIS
      COMMON/O2_VIS/O2VIS(1500:2988)

!     DATA AT 10 CM-1 SPACING:
      DATA (O2VIS(I10),I10=1500,1599)/  .000,   .000,   .000,   .000,   &
     &  .000,   .000,   .000,   .000,   .000,   .000,   .000,   .000,   &
     &  .000,   .000,   .000,   .001,   .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .000,   .000,   .000,   .000606,&
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .002,   &
     &  .002,   .002,   .002,   .002,   .00249, .003,   .003,   .003,   &
     &  .004,   .004,   .005,   .005,   .006,   .007,   .008,   .009,   &
     &  .010,   .011,   .0125,  .0146,  .016,   .018,   .020,   .0223,  &
     &  .025,   .0269,  .030,   .033,   .0363,  .0401,  .0442,  .0467,  &
     &  .0514,  .0555,  .0596,  .0643,  .0694,  .0737,  .0788,  .0838,  &
     &  .0886,  .0937,  .0989,  .103,   .107,   .110,   .114,   .116,   &
     &  .118,   .119,   .120,   .121,   .120,   .120,   .119,   .117,   &
     &  .116,   .113,   .110,   .107,   .103,   .0997,  .0958,  .0915/
      DATA (O2VIS(I10),I10=1600,1699)/  .088,   .0841,  .0794,  .0753,  &
     &  .0717,  .0683,  .0643,  .0608,  .0569,  .0531,  .0502,  .0477,  &
     &  .044,   .0423,  .0394,  .037,   .0351,  .033,   .031,   .029,   &
     &  .0279,  .026,   .025,   .0232,  .022,   .021,   .020,   .019,   &
     &  .018,   .017,   .0165,  .015,   .014,   .013,   .013,   .012,   &
     &  .011,   .011,   .010,   .010,   .009,   .009,   .009,   .008,   &
     &  .008,   .00701, .007,   .007,   .00698, .006,   .0058,  .005,   &
     &  .005,   .005,   .004,   .004,   .004,   .004,   .004,   .004,   &
     &  .004,   .004,   .004,   .004,   .003,   .003,   .003,   .003,   &
     &  .003,   .003,   .003,   .003,   .003,   .003,   .003,   .003,   &
     &  .003,   .003,   .003,   .003,   .003,   .003,   .003,   .003,   &
     &  .003,   .004,   .004,   .004,   .004,   .005,   .005,   .006,   &
     &  .006,   .007,   .00741, .00815, .009,   .0101,  .011,   .012/
      DATA (O2VIS(I10),I10=1700,1799)/  .014,   .015,   .017,   .0185,  &
     &  .0197,  .0224,  .0247,  .0274,  .0306,  .0336,  .037,   .0405,  &
     &  .0449,  .0493,  .0547,  .0601,  .0652,  .0723,  .0789,  .088,   &
     &  .0961,  .105,   .117,   .126,   .139,   .149,   .160,   .168,   &
     &  .174,   .179,   .182,   .184,   .185,   .184,   .183,   .181,   &
     &  .180,   .177,   .174,   .171,   .168,   .164,   .160,   .155,   &
     &  .151,   .146,   .140,   .136,   .130,   .125,   .120,   .114,   &
     &  .109,   .105,   .0993,  .093,   .0888,  .0838,  .0794,  .0751,  &
     &  .0708,  .0666,  .0632,  .0601,  .0555,  .0524,  .0493,  .0463,  &
     &  .0441,  .0415,  .039,   .0363,  .035,   .0326,  .0305,  .0294,  &
     &  .0273,  .0262,  .0246,  .0236,  .0225,  .021,   .020,   .019,   &
     &  .018,   .0176,  .017,   .016,   .015,   .0149,  .014,   .013,   &
     &  .013,   .0122,  .012,   .012,   .011,   .011,   .011,   .010/
      DATA (O2VIS(I10),I10=1800,1899)/  .010,   .010,   .010,   .00916, &
     &  .009,   .009,   .009,   .009,   .00849, .008,   .008,   .008,   &
     &  .008,   .008,   .008,   .008,   .007,   .008,   .007,   .007,   &
     &  .007,   .007,   .007,   .007,   .007,   .007,   .007,   .007,   &
     &  .007,   .007,   .007,   .007,   .007,   .007,   .007,   .007,   &
     &  .007,   .007,   .007,   .007,   .007,   .007,   .007,   .007,   &
     &  .007,   .007,   .007,   .008,   .008,   .008,   .008,   .008,   &
     &  .008,   .009,   .009,   .009,   .00907, .010,   .010,   .010,   &
     &  .011,   .011,   .012,   .0122,  .013,   .0131,  .014,   .015,   &
     &  .016,   .017,   .0182,  .020,   .0201,  .021,   .022,   .0228,  &
     &  .023,   .023,   .023,   .023,   .023,   .023,   .023,   .022,   &
     &  .022,   .022,   .021,   .021,   .020,   .020,   .019,   .019,   &
     &  .0182,  .018,   .0174,  .017,   .0163,  .016,   .015,   .0149/
      DATA (O2VIS(I10),I10=1900,1999)/  .014,   .0137,  .013,   .013,   &
     &  .0121,  .012,   .0113,  .0109,  .010,   .00934, .009,   .00843, &
     &  .008,   .00739, .007,   .006,   .006,   .00574, .005,   .005,   &
     &  .005,   .004,   .004,   .004,   .004,   .004,   .004,   .004,   &
     &  .004,   .004,   .00317, .003,   .003,   .003,   .003,   .003,   &
     &  .003,   .003,   .003,   .003,   .003,   .003,   .003,   .002,   &
     &  .002,   .002,   .002,   .002,   .002,   .002,   .002,   .002,   &
     &  .002,   .002,   .002,   .002,   .002,   .002,   .002,   .002,   &
     &  .00104, .001,   .001,   .001,   .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .001/
      DATA (O2VIS(I10),I10=2000,2099)/  .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .001,   &
     &  .00141, .002,   .002,   .002,   .002,   .002,   .002,   .002,   &
     &  .002,   .002,   .002,   .002,   .002,   .002,   .002,   .002,   &
     &  .00198, .00146, .00105, .001,   .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .001,   &
     &  .00171, .002,   .002,   .002,   .002,   .002,   .003,   .003,   &
     &  .00382, .004,   .00417, .005,   .006,   .007,   .00773, .00807, &
     &  .0097,  .0117,  .0131,  .0147,  .0164,  .0181,  .0207,  .0237,  &
     &  .027,   .0297,  .0327,  .037,   .0413,  .0449,  .0489,  .0538,  &
     &  .0598,  .0645,  .0694,  .0741,  .0801,  .0851,  .090,   .0949,  &
     &  .0988,  .101,   .104,   .107,   .107,   .106,   .103,   .100/
      DATA (O2VIS(I10),I10=2100,2199)/  .0966,  .0893,  .0835,  .0792,  &
     &  .0733,  .0684,  .064,   .0591,  .0557,  .0526,  .0503,  .0475,  &
     &  .0448,  .0426,  .0407,  .0383,  .0369,  .0347,  .0324,  .0311,  &
     &  .0285,  .0269,  .0255,  .0242,  .0221,  .0209,  .0193,  .0177,  &
     &  .0162,  .016,   .0144,  .0136,  .013,   .0116,  .011,   .010,   &
     &  .010,   .009,   .00827, .008,   .00745, .007,   .007,   .00618, &
     &  .006,   .006,   .005,   .005,   .005,   .005,   .004,   .004,   &
     &  .004,   .004,   .004,   .004,   .003,   .003,   .003,   .003,   &
     &  .003,   .003,   .003,   .003,   .00207, .002,   .002,   .002,   &
     &  .002,   .002,   .002,   .002,   .002,   .002,   .002,   .002,   &
     &  .002,   .002,   .002,   .002,   .002,   .002,   .002,   .002,   &
     &  .002,   .002,   .00128, .002,   .002,   .002,   .002,   .002,   &
     &  .002,   .002,   .002,   .002,   .002,   .002,   .002,   .002/
      DATA (O2VIS(I10),I10=2200,2299)/  .002,   .002,   .002,   .002,   &
     &  .002,   .002,   .002,   .002,   .002,   .002,   .002,   .002,   &
     &  .002,   .002,   .002,   .003,   .003,   .003,   .003,   .003,   &
     &  .004,   .004,   .004,   .00457, .005,   .005,   .00564, .006,   &
     &  .00667, .007,   .00735, .008,   .00836, .009,   .009,   .010,   &
     &  .010,   .010,   .010,   .010,   .010,   .010,   .00965, .009,   &
     &  .009,   .008,   .008,   .00769, .007,   .007,   .00644, .006,   &
     &  .006,   .006,   .005,   .005,   .005,   .005,   .005,   .004,   &
     &  .004,   .004,   .004,   .004,   .00398, .00301, .003,   .003,   &
     &  .003,   .003,   .003,   .003,   .003,   .003,   .00254, .002,   &
     &  .002,   .002,   .002,   .002,   .002,   .002,   .002,   .002,   &
     &  .002,   .002,   .002,   .002,   .002,   .002,   .002,   .002,   &
     &  .002,   .002,   .002,   .002,   .002,   .002,   .002,   .002/
      DATA (O2VIS(I10),I10=2300,2399)/  .002,   .002,   .002,   .002,   &
     &  .002,   .002,   .002,   .002,   .002,   .002,   .002,   .002,   &
     &  .002,   .002,   .002,   .00133, .00189, .00107, .00106, .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .001/
      DATA (O2VIS(I10),I10=2400,2499)/  .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .00055, .000,   .000,   &
     &  .001,   .001,   .000751,.000,   .000,   .000,   .000,   .000,   &
     &  .000,   .000,   .000,   .000,   .000,   .000,   .000,   .000,   &
     &  .000,   .000,   .000,   .000,   .000,   .000,   .000,   .000,   &
     &  .000,   .000,   .000,   .000,   .000,   .000,   .000,   .000,   &
     &  .000,   .000,   .000,   .000,   .000,   .000,   .000,   .000,   &
     &  .000,   .000,   .000,   .000,   .000,   .000,   .000134,.001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,  .0000765,&
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .001/
      DATA (O2VIS(I10),I10=2500,2599)/  .001,   .00012, .000,   .000,   &
     &  .000,   .000,   .000,   .000,   .000,   .000,   .000,   .000,   &
     &  .000,   .000,   .000,   .000,   .000,   .000,   .000,   .000,   &
     &  .000,   .000,   .000,   .000,   .000,   .000,   .000,   .000,   &
     &  .000,   .000,   .000,   .000,   .000,   .000,   .000,   .000,   &
     &  .000,   .000,   .000,   .000,   .000,   .000,   .000,   .000,   &
     &  .000,   .000,   .000,   .000,   .000,   .000,   .000,   .000,   &
     &  .000,   .000,   .000,   .000,   .000609,.000347,.000697,.00026, &
     &  .000781,.001,   .001,   .001,   .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .00168, &
     &  .002,   .002,   .002,   .002,   .00276, .003,   .003,   .003/
      DATA (O2VIS(I10),I10=2600,2699)/  .0038,  .004,   .00482, .005,   &
     &  .00584, .006,   .00685, .00785, .00886, .00986, .0109,  .0119,  &
     &  .0129,  .0147,  .0159,  .0177,  .0197,  .0209,  .0227,  .0247,  &
     &  .0267,  .0287,  .0307,  .0326,  .0338,  .0356,  .0368,  .0386,  &
     &  .039,   .0398,  .0407,  .041,   .041,   .0403,  .0393,  .0383,  &
     &  .0373,  .0364,  .0348,  .0334,  .0318,  .0299,  .0285,  .027,   &
     &  .025,   .0231,  .0211,  .0192,  .0176,  .0163,  .0147,  .0134,  &
     &  .0117,  .0107,  .00978, .00881, .00784, .00688, .006,   .00594, &
     &  .005,   .005,   .00405, .004,   .00313, .003,   .003,   .00224, &
     &  .002,   .002,   .002,   .002,   .002,   .002,   .00154, .00141, &
     &  .00164, .001,   .001,   .001,   .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .001/
      DATA (O2VIS(I10),I10=2700,2799)/  .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .001,   .001,   &
     &  .001,   .001,   .001,   .001,   .001,   .001,   .00115, .002,   &
     &  .002,   .002,   .002,   .002,   .002,   .00256, .003,   .003,   &
     &  .0033,  .004,   .004,   .00404, .00495, .00585, .006,   .00667, &
     &  .00758, .00848, .00939, .0103,  .0114,  .0131,  .014,   .0158,  &
     &  .0176,  .0194,  .0212,  .023,   .0256,  .0289,  .0316,  .0344,  &
     &  .038,   .0416,  .0452,  .0487,  .0523,  .0559,  .0591,  .062,   &
     &  .0653,  .0671,  .0689,  .0698,  .0707,  .071,   .071,   .0706,  &
     &  .0697,  .0689,  .068,   .0671,  .0654,  .0643,  .0629,  .0611,  &
     &  .0594,  .0574,  .0548,  .0531,  .0505,  .0486,  .0462,  .0441,  &
     &  .0423,  .0403,  .0378,  .0361,  .0343,  .0326,  .0308,  .0291/
      DATA (O2VIS(I10),I10=2800,2899)/  .0273,  .0258,  .0249,  .0231,  &
     &  .0222,  .0207,  .0195,  .0186,  .0177,  .0169,  .016,   .0151,  &
     &  .0143,  .014,   .0135,  .0127,  .0118,  .011,   .011,   .0102,  &
     &  .010,   .010,   .00967, .00881, .00805, .0089,  .00824, .008,   &
     &  .00753, .007,   .007,   .007,   .007,   .007,   .00642, .006,   &
     &  .006,   .006,   .006,   .00518, .005,   .005,   .005,   .0048,  &
     &  .00404, .00489, .00427, .004,   .004,   .004,   .004,   .004,   &
     &  .004,   .004,   .004,   .004,   .004,   .004,   .0032,  .003,   &
     &  .003,   .003,   .003,   .003,   .003,   .003,   .003,   .003,   &
     &  .003,   .003,   .003,   .003,   .003,   .003,   .003,   .003,   &
     &  .00375, .004,   .004,   .004,   .004,   .004,   .00469, .005,   &
     &  .005,   .00515, .00597, .006,   .00661, .00743, .008,   .00806, &
     &  .00888, .0097,  .0105,  .0113,  .0121,  .013,   .0138,  .0152/
      DATA (O2VIS(I10),I10=2900,2988)/  .0164,  .0172,  .018,   .0188,  &
     &  .0196,  .0204,  .021,   .021,   .021,   .021,   .021,   .021,   &
     &  .021,   .021,   .021,   .0205,  .020,   .0199,  .0191,  .019,   &
     &  .0185,  .018,   .0179,  .0171,  .0163,  .0155,  .0147,  .014,   &
     &  .014,   .0133,  .0125,  .012,   .0119,  .0111,  .0103,  .010,   &
     &  .00975, .009,   .009,   .00837, .008,   .008,   .008,   .00722, &
     &  .007,   .007,   .00686, .00607, .006,   .006,   .006,   .00593, &
     &  .00515, .005,   .005,   .005,   .005,   .005,   .005,   .005,   &
     &  .005,   .005,   .005,   .005,   .005,   .005,   .005,   .005,   &
     &  .00468, .004,   .004,   .004,   .004,   .004,   .004,   .004,   &
     &  .004,   .004,   .004,   .004,   .004,   .004,   .004,   .004,   &
     &  .004,   .001,   .0002,  .000,   .000/
      END
