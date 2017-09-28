#!/bin/sh
cd obj-ifort 
ar rc libmodpcrtm.a *.o
mv libmodpcrtm.a ../modpcrtm_lib
cd ../



