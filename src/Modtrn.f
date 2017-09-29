      PROGRAM MDTRN5
      CALL DRIVER

!     AN ASCII VERSION OF SECTION 2 OF THE MODTRAN4 USER's
!     Manual is included below for reference.

!     2.  OVERVIEW OF INPUT DATA FORMAT

!        An attempt has been made in MODTRAN4 to make it easier for the
!     users to keep track of input and output (I/O) files.  The need
!     for easier file handling is evident to anyone who runs MODTRAN
!     using different tape5 input files and who wants to save the
!     corresponding output files (the tape6, pltout, tape7, and so on).
!     In the past, every MODTRAN input file had to have the name
!     'TAPE5' and previously generated I/O files had to be renamed to
!     avoid overwriting them with newer files.  The need for renaming
!     is now avoided by creating a new MODTRAN input file (referred to
!     as the root name file) called either 'MODROOT.IN' or 'MODROOT.IN'.
!     If 'MODROOT.IN' does not exist, MODTRAN checks for the existence
!     of a 'MODROOT.IN' file.  If neither of these files exists,
!     MODTRAN I/O files are the traditional ones: 'TAPE5', 'TAPE6',
!     'TAPE7', 'TAPE8', etc.  If a root name file exists and its very
!     first line contains a non-null string (maximum length is 80
!     characters), this string is treated as a prefix.  If the string
!     consists of all blanks or is a null string, the traditional
!     I/O file names are assumed.  The root name should contain no
!     embedded blanks; leading and trailing blanks are properly
!     ignored.  This string is used as a prefix for the I/O files
!     whose names have mnemonic suffixes.  As an example, if the
!     string is case1, the MODTRAN I/O files will have these names:

!           case1.tp5   (corresponding to tape5)
!           case1.tp6   (corresponding to tape6)
!           case1.tp7   (corresponding to tape7)
!           case1.tp8   (corresponding to tape8)
!           case1.7sc   (corresponding to tape7.scn)
!           case1.7sr   (corresponding to tape7.scr)
!           case1.plt   (corresponding to pltout)
!           case1.psc   (corresponding to pltout.scn)
!           case1.clr   (corresponding to clrates)
!           case1.chn   (corresponding to channels.out)
!           case1.flx   (corresponding to specflux)

!        MODTRAN is controlled by a single input file, 'TAPE5'
!     or 'ROOTNAME.TP5', which consists of a sequence of six
!     or more CARDS (inputs lines).  The input file format is
!     summarized below.  Except when specifying file names,
!     character inputs are case insensitive.  Also, blanks are
!     read as zeroes for numerical inputs, and as default values
!     otherwise.  Detailed descriptions of the card formats
!     and parameters are given in the following sections.

!     2.1  Listing of CARDs and Their Format

!        In the following, optional cards are indented.

!     CARD 1:   MODTRN, SPEED, MODEL, ITYPE, IEMSCT, IMULT, M1,  M2,
!               M3, M4, M5, M6, MDEF, IM, NOPRNT, TPTEMP, SURREF
!               FORMAT (2A1, I3, 12I5, F8.3, A7)

!     CARD 1A:  DIS, DISAZM, NSTR, LSUN, ISUN, CO2MX, H2OSTR, O3STR,
!               LSUNFL, LBMNAM, LFLTNM, H2OAER, SOLCON
!               FORMAT (2L1, I3, L1, I4, F10.5, 2A10,
!                      4(1X, A1), 2X, F10.3)

!         CARD 1A1:  USRSUN
!                    FORMAT (A80)   (If LSUNFL = True)

!         CARD 1A2:  BMROOT
!                    FORMAT (A80)   (If LBMNAM = True)

!         CARD 1A3:  FILTNM
!                    FORMAT (A80)   (If LFLTNM = True)

!     CARD 2:   APLUS, IHAZE, CNOVAM, ISEASN, ARUSS, IVULCN,
!               ICSTL,  ICLD, IVSA, VIS, WSS, WHH, RAINRT, GNDALT
!               FORMAT (A2, I3, A1, I4, A3, I2, 3I5, 5F10.5)

!         CARD 2A+:  ZAER11, ZAER12, SCALE1, ZAER21, ZAER22, SCALE2,
!                    ZAER31, ZAER32, SCALE3, ZAER41, ZAER42, SCALE4
!                    FORMAT ((3(1X, F9.0), 20X, 3(1X, F9.0)))
!                                                   (If APLUS = 'A+')

!         CARD 2A:   CTHIK, CALT, CEXT
!                    FORMAT (3F8.3)   (If ICLD = 18 or 19)

!         Alternate CARD 2A:   CTHIK, CALT, CEXT, NCRALT, NCRSPC,
!                    CWAVLN, CCOLWD, CCOLIP, CHUMID, ASYMWD, ASYMIP
!                    FORMAT (3F8.3, 2I4, 6F8.3)   (If ICLD = 1-10)

!         CARD 2B:   ZCVSA, ZTVSA, ZINVSA
!                    FORMAT (3F10.3)   (If IVSA = 1)

!         CARD 2C:   ML, IRD1, IRD2, TITLE
!                    FORMAT (3I5, A65)   (If MODEL = 0 or 7, and IM = 1)

!     CARDs 2C1, 2C2, 2C2X, and 2C3 (as required)
!     are each repeated ML times.

!         CARD 2C1:  ZM, P, T, WMOL(1), WMOL(2), WMOL(3),
!                    (JCHAR(J), J = 1, 14), JCHARX
!                    FORMAT (F10.3, 5E10.3, 14A1, 1X, A1)

!         CARD 2C2:  (WMOL(J), J = 4, 12)
!                    FORMAT (8E10.3, /E10.3)   (If IRD1 = 1)

!         CARD 2C2X: (WMOLX(J), J = 1, 15)
!                    FORMAT (8E10.3, /5E10.3)   (If MDEF = 2 & IRD1 = 1)

!         CARD 2C3:  AHAZE, EQLWCZ, RRATZ,
!                    IHA1, ICLD1, IVUL1, ISEA1, ICHR1
!                    FORMAT (10X, 3F10.3, 5I5)   (If IRD2 = 1)

!         CARD 2D:   (IREG(N), N = 1, 2, 3, 4)
!                    FORMAT (4I5)   (If IHAZE = 7 or ICLD = 11)

!         CARD 2D1:  AWCCON, TITLE
!                    FORMAT (E10.3, A70)

!     If ARUSS = 'USS' and IREG(N) > 1, then
!     Imax = IREG(N); Else Imax = 47

!         CARD 2D2:  (VX(N, I), EXTC(N, I), ABSC(N, I),
!                    ASYM(N, I), I = l, 2, ..., Imax)
!                    FORMAT ((3(F6.2, 2F7.5, F6.4)))

!         CARD 2E1:  (ZCLD(I, 0), CLD(I, 0), CLDICE(I, 0),
!                    RR(I, 0), I = 1, NCRALT)
!                    FORMAT((4F10.5))    (If ICLD = 1 - 10, NCRALT > 2)

!         CARD 2E2:  (WAVLEN(I), EXTC(6, I), ABSC(6, I), ASYM(6, I),
!                    EXTC(7, I), ABSC(7, I), ASYM(7, I), I = 1, NCRSPC)
!                    FORMAT((7F10.5))   (If ICLD = 1 - 10, NCRSPC > 1)

!     CARD 3:   H1, H2, ANGLE, HRANGE, BETA, RO, LENN, PHI
!               FORMAT (6F10.3, I5, 5X, F10.3)

!     Alternate CARD 3:   H1, H2, ANGLE, IDAY, RO, ISOURC, ANGLEM
!               FORMAT (3F10.3, I5, 5X, F10.3, I5, F10.3)
!                                                 (If IEMSCT = 3)

!         CARD 3A1:   IPARM, IPH, IDAY, ISOURC
!                     FORMAT (4I5)   (If IEMSCT = 2)

!         CARD 3A2:   PARM1, PARM2, SRCLAT, SRCLON, TIME,TRUEAZ,ANGLEM,G
!                     FORMAT (8F10.3)   (If IEMSCT = 2)

!         CARD 3B1:   NANGLS, NWLF
!                     FORMAT (2I5)   (If IPH = 1)

!         CARD 3B2:   (ANGF(I), F(1, I, 1), F(2, I, 1),
!                     F(3, I, 1), F(4, I, 1), I = l, NANGLS)
!                     FORMAT (8(1X, F9.0))   (If IPH = 1 and NWLF = 0)

!         CARD 3C1:   (ANGF(I), I = 1, NANGLS)
!                     FORMAT (8(1X, F9.0))   (If IPH = 1 and NWLF > 0)

!         CARD 3C2:   (WLF(J), J = 1, NWLF)
!                     FORMAT (8(1X, F9.0))   (If IPH = 1 and NWLF > 0)

!    In CARDs 3C3-3C6, 'I' is angle index as in CARD 3C1
!    and 'J' is the wavelength index as in CARD 3C2.

!         CARD 3C3:   (F(1, I, J), J = 1, NWLF)
!                     FORMAT (8(1X, E9.3))   (If IPH = 1 and NWLF > 0)

!         CARD 3C4:   (F(2, I, J), J = 1, NWLF)
!                     FORMAT (8(1X, E9.3))   (If IPH = 1 and NWLF > 0)

!         CARD 3C5:   (F(3, I, J), J = 1, NWLF)
!                     FORMAT (8(1X, E9.3))   (If IPH = 1 and NWLF > 0)

!         CARD 3C6:   (F(4, I, J), J = 1, NWLF)
!                     FORMAT (8(1X, E9.3))   (If IPH = 1 and NWLF > 0)

!     CARD 4:   V1, V2, DV, FWHM, YFLAG, XFLAG, DLIMIT, FLAGS
!               FORMAT (4F10.0, 2A1, A8, A7)

!         CARD 4A:    NSURF, AATEMP
!                     FORMAT (I1, F9.0)
!                     (If SURREF = 'BRDF' or 'LAMBER')

!     The set of CARD 4B1, 4B2, and 4B3 inputs is repeated NSURF times.

!         CARD 4B1:  CBRDF
!                    FORMAT (A80)    (If SURREF = 'BRDF')

!         CARD 4B2:   NWVSRF, SURFZN, SURFAZ
!                     FORMAT (*)   (If SURREF = 'BRDF')

!     CARD 4B3 is repeated NWVSRF times.

!         CARD 4B3:   WVSURF, (PARAMS(I), I = 1, NPARAM)
!                     FORMAT (*)   (If SURREF = 'BRDF')

!         CARD 4L1:   SALBFL
!                     FORMAT (A80)    (If SURREF = 'LAMBER')

!     CARD 4L2 is repeated NSURF times.

!         CARD 4L2:   CSALB
!                     FORMAT (A80)    (If SURREF = 'LAMBER')

!     CARD 5:   IRPT
!               FORMAT (I5)
      END
