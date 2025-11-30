#!/usr/bin/env bash

#source /cvmfs/sw.hsf.org/key4hep/setup.sh
export LIBGL_ALWAYS_SOFTWARE=1
#export LIBGL_ALWAYS_INDIRECT=1


# Absolute path to the directory that this script lives in
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Run the sibling executable "Druid"
"$SCRIPT_DIR/Druid" -coll.caloHit.filterOutSuffix Digi \
 -coll.simCaloHit.filterOutSuffix HitsEven \
 -coll.simCaloHit.filterOutSuffix HitsOdd \
 -coll.track.add MarlinTrkTracks \
 -coll.track.add ClupatraTrackSegments \
 -coll.cluster.add PandoraClusters \
 -tpc.innerRadius 330.0 \
 -tpc.outerRadius 1700.0 \
 -tpc.halfZ 2350.0 \
 -draw.tpc.cylinder \
 -logLevel 3 \
 "$@"


