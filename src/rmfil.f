      SUBROUTINE RMFIL
      IMPLICIT NONE
      INCLUDE 'PARAMS.h'
!      REMOVES TEMPERARY FILES FROM DRIVER
!      FILES REMOVED INCLUDE:

      INTEGER LNFLRT
      LOGICAL LRDSUN
      COMMON/SUNFLG/LNFLRT,LRDSUN

      COMMON /PRONUM/LABEL

      CHARACTER LABEL*6
      INTEGER STAT


      COMMON/SUNNAM/FLRT,SUNFIL
      CHARACTER FLRT*(NAMLEN-4),SUNFIL*(LENSUN)

      open(unit=25, iostat=stat, file=FLRT(1:LNFLRT)//LABEL//'.acd',     &
     &       status='old')
      if (stat == 0) close(25, status='delete')
!      open(unit=2, iostat=stat, file=FLRT(1:LNFLRT)//LABEL//'.tp6',     & 
!     &     status='old')
!      if (stat == 0) close(2, status='delete')

      open(unit=7, iostat=stat, file=FLRT(1:LNFLRT)//LABEL//'.tp7',     & 
     &     status='old')
      if (stat == 0) close(7, status='delete')
      open(unit=29, iostat=stat, file=FLRT(1:LNFLRT)//LABEL//'b.7sr',   &
     &     status='old')
      if (stat == 0) close(29, status='delete')

      open(unit=15, iostat=stat, file=FLRT(1:LNFLRT)//LABEL//'.sc',     &
     &     status='old')
      if (stat == 0) close(15, status='delete')

      open(unit=15, iostat=stat, file=FLRT(1:LNFLRT)//LABEL//'.7sc',    &
     &     status='old')
      if (stat == 0) close(15, status='delete')

      RETURN
       END

