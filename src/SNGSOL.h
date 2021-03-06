
!     /SNGSOL/
!       PHSANG   SOLAR SCATTERING PHASE ANGLE ALONG LINE-OF-SIGHT [DEG].
!       PHSCOS   COSINE OF SOLAR SCATTERING PHASE ANGLE.
!       IUPPHS   PHASE ANGLE UPPER INTERPOLATION INDEX FOR ANGF ARRAY.
!       PHSFAC   PHASE ANGLE INTERPOLATION FRACTION FOR ANGF ARRAY.
!       PF_HG    SPECTRALLY INDEPENDENT HENYEY-GREENSTEIN FNC [SR-1].
!       PF_WAT   SPECTRALLY INDEPENDENT WATER DROPLET PHASE FNC [SR-1].
!       PF_ICE   SPECTRALLY INDEPENDENT ICE PARTICLE PHASE FNC [SR-1].
      REAL PHSFAC,PHSANG,PHSCOS,PF_HG,PF_WAT,PF_ICE
      INTEGER IUPPHS
      COMMON/SNGSOL/PHSANG(0:LAYTWO,1:MLOS),PHSCOS(0:LAYTWO,1:MLOS),    &
     &              IUPPHS(0:LAYTWO,1:MLOS),PHSFAC(0:LAYTWO,1:MLOS),    &
     &               PF_HG(0:LAYTWO,1:MLOS),PF_WAT(0:LAYTWO,1:MLOS),    &
     &              PF_ICE(0:LAYTWO,1:MLOS)
