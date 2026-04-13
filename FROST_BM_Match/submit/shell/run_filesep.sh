#!/bin/sh

cd ${NINJARECONDIR}/bin/FileSeparator

b2datadir=${YASUDIR}/Latest/testbench/data/reconstruction/physdata
b2datafile=${b2datadir}/neutrino_b2physics_fullsetup_fullstat_timedifcut0_loose_$1_$2_$3.root

ninjadatadir=${HOME}/data/ninja-root
if [ $1 = 2019 ]; then
    ninjadatafile=${ninjadatadir}/All_1573054729_conv.root
elif [ $1 = 2020 ]; then
    ninjadatafile=${ninjadatadir}/All_1578809458_conv.root 
fi

outputdir=${HOME}/data/filesep
outputfile=${outputdir}/ninja_rawdata_$1_$2_$3.root

#echo ${b2datafile} ${ninjadatafile} ${outputfile}
./FileSeparator ${b2datafile} ${ninjadatafile} ${outputfile}

