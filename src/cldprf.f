      REAL FUNCTION CLDPRF(EQLWCZ,RELHUM,ICLD,IHA1,IC1,ICHIC1)

!     CLDPRF COMPUTES DENSITY PROFILES FOR CLOUDS / WATER AEROSOLS.
      IMPLICIT NONE

!     ARGUMENTS:
      REAL EQLWCZ,RELHUM
      INTEGER ICLD,IHA1,IC1,ICHIC1

!     COMMONS:
      INCLUDE 'IFIL.h'

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     LOCAL VARIABLES:
      INTEGER IM1,I
      REAL FAC

!     DATA:
      REAL LN20,LN30,LN100,AFLWC,RFLWC,CULWC,ASLWC,STLWC,               &
     &  SCLWC,SNLWC,BSLWC,FVLWC,AVLWC,MDLWC,TNLWC,TKLWC,                &
     &  ELWCR(4),ELWCM(4),ELWCU(4),ELWCT(4)
      DATA LN20,LN30,LN100/2.995732,3.401197,4.605170/,AFLWC/.01295/,   &
     &  RFLWC/.001804/,CULWC/.007683/,ASLWC/.004509/,STLWC/.005272/,    &
     &  SCLWC/.004177/,SNLWC/.007518/,BSLWC/.0001567/,FVLWC/.0005922/,  &
     &  AVLWC/.0001675/,MDLWC/.0004775/,TNLWC/.05811/,TKLWC/.003446/,   &
     &  ELWCR/.0003517,.0003740,.0004439,.0009529/,                     &
     &  ELWCM/.0004675,.0006543,.001166,.003154/,                       &
     &  ELWCU/.0003102,.0003802,.0004463,.0009745/,                     &
     &  ELWCT/.0001735,.0001820,.0002020,.0002408/
      CLDPRF=0.
      IF(ICLD.GT.0)THEN
          IF(ICLD.EQ.18)THEN
              CLDPRF=EQLWCZ/TNLWC
          ELSEIF(ICLD.EQ.19) THEN
              CLDPRF=EQLWCZ/TKLWC
          ELSEIF(ICLD.EQ.1 .OR. (ICLD.GE.9 .AND. ICLD.LE.11))THEN
              CLDPRF=EQLWCZ/CULWC
          ELSEIF(ICLD.EQ.2)THEN
              CLDPRF=EQLWCZ/ASLWC
          ELSEIF(ICLD.EQ.3 .OR. ICLD.EQ.6)THEN
              CLDPRF=EQLWCZ/STLWC
          ELSEIF(ICLD.EQ.4)THEN
              CLDPRF=EQLWCZ/SCLWC
          ELSEIF(ICLD.EQ.5 .OR. ICLD.EQ.7 .OR. ICLD.EQ.8)THEN
              CLDPRF=EQLWCZ/SNLWC
          ENDIF
      ELSEIF(IHA1.GT.0)THEN
          IF(ICHIC1.LE.6)THEN

!             RH DEPENDENT AEROSOLS
              IF(ICHIC1.EQ.6 .AND. IC1.GT.1)THEN

!                 FOR TROPOSPHERIC AEROSOL, USE 70% RH ABOVE EH(7,I).
                  IM1=3
                  I=3
                  FAC=0.
              ELSEIF(RELHUM.LT.70.)THEN
                  IM1=1
                  I=2
                  FAC=(LN100-LOG(100-RELHUM))/(LN100-LN30)
              ELSEIF(RELHUM.LT.80.)THEN
                  IM1=2
                  I=3
                  FAC=(LN30-LOG(100-RELHUM))/(LN30-LN20)
              ELSEIF(RELHUM.LT.99.)THEN
                  IM1=3
                  I=4
                  FAC=(LN20-LOG(100-RELHUM))/LN20
              ELSE
                  IM1=4
                  I=4
                  FAC=0.
              ENDIF
              IF(ICHIC1.EQ.6)THEN
                  CLDPRF=EQLWCZ/(ELWCT(IM1)*(ELWCT(I)/ELWCT(IM1))**FAC)
              ELSEIF(ICHIC1.EQ.5)THEN
                  CLDPRF=EQLWCZ/(ELWCU(IM1)*(ELWCU(I)/ELWCU(IM1))**FAC)
              ELSEIF(ICHIC1.GE.3)THEN
                  CLDPRF=EQLWCZ/(ELWCM(IM1)*(ELWCM(I)/ELWCM(IM1))**FAC)
              ELSE
                  CLDPRF=EQLWCZ/(ELWCR(IM1)*(ELWCR(I)/ELWCR(IM1))**FAC)
              ENDIF
          ELSEIF(ICHIC1.EQ.8)THEN
              CLDPRF=EQLWCZ/AFLWC
          ELSEIF(ICHIC1.EQ.9)THEN
              CLDPRF=EQLWCZ/RFLWC
          ELSEIF(ICHIC1.EQ.11 .OR. ICHIC1.EQ.16 .OR. ICHIC1.EQ.17)THEN
              CLDPRF=EQLWCZ/BSLWC
          ELSEIF(ICHIC1.EQ.12 .OR. ICHIC1.EQ.14)THEN
              CLDPRF=EQLWCZ/AVLWC
          ELSEIF(ICHIC1.EQ.13 .OR. ICHIC1.EQ.15 .OR. ICHIC1.EQ.18)THEN
              CLDPRF=EQLWCZ/FVLWC
          ELSE
              CLDPRF=EQLWCZ/MDLWC
          ENDIF
      ELSE
          WRITE(IPR,'(/A)')' WARNING:  ICLD not set in cldprf.f'
      ENDIF
      RETURN
      END
