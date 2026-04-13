#!/bin/sh

cd ${NINJARECONDIR}/bin/TrackMatch

b2datadir=${HOME}/data/hitconv
b2datafile=${b2datadir}/neutrino_b2physics_fullsetup_fullstat_timedifcut0_loose_ninja_$1_$2_$3.root

matchdir=${HOME}/data/trackmatch
matchfile=${matchdir}/neutrino_b2physics_fullsetup_fullstat_timedifcut0_loose_ninjamatch_$1_$2_$3.root

./TrackMatch ${b2datafile} ${matchfile}
