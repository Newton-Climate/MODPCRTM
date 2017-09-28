      SUBROUTINE  ZEROAL(ND1, EXPBEA, FLYR, OPRIM, PHASA, PHASE,PHASM,  &
     &                         TAUCPR, XR0, XR1,                        &
     &                    ND2, CMU, CWT, PSI, WK, Z0, Z1, ZJ,           &
     &                    ND3, YLM0,                                    &
     &                    ND4, ARRAY, CC, EVECC,                        &
     &                    ND5, GL,                                      &
     &                    ND6, YLMC,                                    &
     &                    ND7, YLMU,                                    &
     &                    ND8, KK, LL, ZZ, ZPLK0, ZPLK1,                &
     &                    ND9, GC,                                      &
     &                    ND10, LAYRU, UTAUPR,                          &
     &                    ND11, GU,                                     &
     &                    ND12, Z0U, Z1U, ZBEAM,                        &
     &                    ND13, EVAL,                                   &
     &                    ND14, AMB, APB,                               &
     &                    ND15, IPVT, Z)

!                                                            Zero arrays

      INTEGER  IPVT(*), LAYRU(*)
      REAL  AMB(*), APB(*), ARRAY(*), CC(*), CMU(*), CWT(*),            &
     &      EVAL(*), EVECC(*), EXPBEA(*), FLYR(*), GC(*),               &
     &      GL(*), GU(*), KK(*), LL(*), OPRIM(*),                       &
     &      PHASA(*), PHASE(*), PHASM(*), PSI(*),                       &
     &      TAUCPR(*), UTAUPR(*), WK(*), XR0(*), XR1(*),                &
     &      YLM0(*), YLMC(*), YLMU(*), Z(*), Z0(*),                     &
     &      Z0U(*), Z1(*), Z1U(*), ZBEAM(*), ZJ(*),                     &
     &      ZZ(*), ZPLK0(*), ZPLK1(*)

      DO N=1, ND1
         EXPBEA(N)= 0.0
         FLYR(N) =0.0
         OPRIM(N)=0.0
         PHASA(N)=0.0
         PHASE(N)=0.0
         PHASM(N)=0.0
         TAUCPR(N)= 0.0
         XR0(N)  =0.0
         XR1(N)  =0.0
      ENDDO

      DO N=1, ND2
         CMU(N)=0.0
         CWT(N)=0.0
         PSI(N)=0.0
         WK(N) =0.0
         Z0(N) =0.0
         Z1(N) =0.0
         ZJ(N) =0.0
      ENDDO

      DO N=1, ND3
         YLM0(N)=0.0
      ENDDO

      DO N=1, ND4
         ARRAY(N)=0.0
         CC(N)   =0.0
         EVECC(N)=0.0
      ENDDO

      DO N=1, ND5
         GL(N)=0.0
      ENDDO

      DO N=1, ND6
         YLMC(N)=0.0
      ENDDO

      DO N=1, ND7
         YLMU(N)=0.0
      ENDDO

      DO N=1, ND8
         KK(N)=0.0
         LL(N)=0.0
         ZZ(N)=0.0
         ZPLK0(N)=0.0
         ZPLK1(N)=0.0
      ENDDO

      DO N=1, ND9
         GC(N)=0.0
      ENDDO

      DO N=1, ND10
         LAYRU(N) =0
         UTAUPR(N)=0.0
      ENDDO

      DO N=1, ND11
         GU(N)=0.0
      ENDDO

      DO N=1, ND12
         Z0U(N)=0.0
         Z1U(N)=0.0
         ZBEAM(N)=0.0
      ENDDO

      DO N=1, ND13
         EVAL(N)=0.0
      ENDDO

      DO N=1, ND14
         AMB(N)=0.0
         APB(N)=0.0
      ENDDO

      DO N=1, ND15
         IPVT(N)=0
         Z(N)=0.0
      ENDDO

      RETURN
      END
