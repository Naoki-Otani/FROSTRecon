#!/bin/sh

cd ${NINJARECONDIR}/bin/HitConverter

b2datadir=${YASUDIR}/Latest/testbench/data/reconstruction/physdata
b2datafile=${b2datadir}/neutrino_b2physics_fullsetup_fullstat_timedifcut0_loose_$1_$2_$3.root

ninjadatadir=${HOME}/data/filesep
ninjadatafile=${ninjadatadir}/ninja_rawdata_$1_$2_$3.root

outputdir=${HOME}/data/hitconv
outputfile=${outputdir}/neutrino_b2physics_fullsetup_fullstat_timedifcut0_loose_ninja_$1_$2_$3.root

#echo ${b2datafile} ${ninjadatafile} ${outputfile}
./HitConverter ${b2datafile} ${ninjadatafile} ${outputfile}

