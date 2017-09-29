
!     COMMON BLOCK for PCRTM related calculations

      INTEGER MxMono, MxK
      PARAMETER(MxMono=2000,MxK=17)
      INTEGER MxPC, MxChnl, MxBnd,MxReg
      PARAMETER(MxPC=500,MxChnl=18000,MxBnd=5,MxReg=1000)
      integer*4 nb
      integer*4  nReg, nPcBnd, nch
      real*8 radStdCh,frqChBnd,radStd
      real*8 PC,RadMeanCh
      real*4 regCoef
      

      INTEGER ChanCheck
      INTEGER IKcheck
      INTEGER nMono, nK
      LOGICAL initflag

      common/PCparameter/radStdCh(MxChnl,MxBnd),
     & frqChBnd(MxChnl,MxBnd),radStd(mxReg,MxBnd),
     & PC(mxPc,MxChnl,MxBnd),
     & RadMeanCh(MxChnl,MxBnd),regCoef(mxReg,MxPC,MxBnd),
     & nReg(MxBnd), nPcBnd(MxBnd),nch(MxBnd),nb



      common/PCindex/ChanCheck(mxmono),IKcheck(mxmono, mxK),
     & nMono,nK,initflag
