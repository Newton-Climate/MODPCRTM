#!/bin/csh 
date
set n=$1                                              ## Pass in a unique run number
echo n = $n
echo $PWD
if ( -e  /nobackupp2/nhnguye3/modpcrtm/list_tape5_$n.asc) then   ## Look for existing list file.
   set PBS_O_WORKDIR = /nobackupp2/nhnguye3/modpcrtm/run_$n

   rm -fr /nobackupp2/nhnguye3/modpcrtm/run_$n                   ## Setup unique run environment
   mkdir -p /nobackupp2/nhnguye3/modpcrtm/run_$n
   cd /nobackupp2/nhnguye3/modpcrtm/run_$n
   cp    /nobackupp2/nhnguye3/modpcrtm/Mod5300alpha-ifort.exe mod5.exe
   ln -s /nobackupp2/nhnguye3/modpcrtm/DATA .
   ln -s /nobackupp2/nhnguye3/modpcrtm/support .
   cp    /nobackupp2/nhnguye3/modpcrtm/pcrtm.in pcrtm.in
   cp    /nobackupp2/nhnguye3/modpcrtm/output.tp5 .
   cp    /nobackupp2/nhnguye3/modpcrtm/list_tape5_$n.asc list_tape5.asc

   echo "Executing run $n on" `hostname` "in $PWD"

   set PBS_O_WORKDIR = /nobackupp2/nhnguye3/modpcrtm
   ./mod5.exe > output_$n
pwd
   rm -f DATA
   rm -f mod5.exe
   rm -f pcrtm.in
   rm -f pltout.asc
#  rm -f radiance.out
   rm -f tape?
   rm -f tape7.scn
   rm -f atmcor.dat
   rm -f pltout
   rm -f specflux
   rm -f support
   rm -f tape7b.scr
   #stat -c %s output_$n
   if ( -s output_$n ) then
   else
      rm -f output_$n
   endif

endif
