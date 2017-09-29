#!/bin/csh -fe
date
set n=$1                                              ## Pass in a unique run number
echo n = $n
if ( -e  /u/skizer/modpcrtm/list_tape5_$n.asc) then   ## Look for existing list file.
   set PBS_O_WORKDIR = /u/skizer/modpcrtm/run_$n

   rm -fr /u/skizer/modpcrtm/run_$n                   ## Setup unique run environment
   mkdir -p /u/skizer/modpcrtm/run_$n
   cd /u/skizer/modpcrtm/run_$n
   cp    /u/skizer/modpcrtm/Mod5300alpha-ifort.exe mod5.exe
   ln -s /u/skizer/modpcrtm/DATA .
   ln -s /u/skizer/modpcrtm/support .
   cp    /u/skizer/modpcrtm/pcrtm.in pcrtm.in
   cp    /u/skizer/modpcrtm/list_tape5_$n.asc list_tape5.asc

   echo "Executing run $n on" `hostname` "in $PWD"

   set PBS_O_WORKDIR = /u/skizer/modpcrtm

   ./mod5.exe > output_$n

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
