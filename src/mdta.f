      BLOCK DATA MDTA

!     CLOUD AND RAIN DATA
      INCLUDE 'PARAMS.h'

!     /CLDRR/
!       ZCLD     INPUT (*,0) & MODEL (*,1) CLOUD ALTITUDE PROFILES [KM].
!       PCLD     CLOUD PROFILE PRESSURE [MB].
!       CLD      INPUT (*,0) & MODEL (*,>0) WATER CLOUD PROFILES [G/M3].
!       CLDICE   INPUT (*,0) & MODEL (*,>0) ICE CLOUD PROFILES [G/M3].
!       RR       INPUT (*,0) & MODEL (*,>0) RAIN RATE PROFILES [MM/HR].
      DOUBLE PRECISION ZCLD
      REAL PCLD,CLD,CLDICE,RR
      COMMON/CLDRR/ZCLD(1:NZCLD,0:1),PCLD,CLD(1:NZCLD,0:5),             &
     &  CLDICE(1:NZCLD,0:1),RR(1:NZCLD,0:5)
      SAVE /CLDRR/

!     PROFILE ALTITUDES [KM]
      DATA ZCLD/NZCLD*0.D0,                                             &
     &   .00D0,  .16D0,  .33D0,  .66D0, 1.00D0, 1.50D0, 2.00D0, 2.40D0, &
     &  2.70D0, 3.00D0, 3.50D0, 4.00D0, 4.50D0, 5.00D0, 5.50D0, 6.00D0/

!     WATER DROPLET DENSITIES [G/M3]
      DATA CLD/NZCLD*0,                                                 &
     &    .00,  .00,  .00,  .20,  .35, 1.00, 1.00,1.0, .3, .15, .0,5*.0,&
     &    .00,  .00,  .00,  .00,  .00,  .00,  .00, .3, .4, .30, .0,5*.0,&
     &    .00,  .00,  .15,  .30,  .15,  .00,  .00, .0, .0, .00, .0,5*.0,&
     &    .00,  .00,  .00,  .10,  .15,  .15,  .10, .0, .0, .00, .0,5*.0,&
     &    .00,  .30,  .65,  .40,  .00,  .00,  .00, .0, .0, .00, .0,5*.0/

!     ICE PARTICLE DENSITIES [G/M3]
      DATA CLDICE/NZCLD*0.,NZCLD*0./

!     RAIN RATES [MM/HR]
      DATA RR/NZCLD*0.,                                                 &
     &   2.00, 1.78, 1.43, 1.22,  .86,  .22,  .00, .0, .0, .00, .0,5*.0,&
     &   5.00, 4.00, 3.40, 2.60,  .80,  .20,  .00, .0, .0, .00, .0,5*.0,&
     &  12.50,10.50, 8.00, 6.00, 2.50,  .80,  .20, .0, .0, .00, .0,5*.0,&
     &  25.00,21.50,17.50,12.00, 7.50, 4.20, 2.50,1.0, .7, .20, .0,5*.0,&
     &  75.00,70.00,65.00,60.00,45.00,20.00,12.50,7.0,3.5,1.00, .2,5*.0/
      END
