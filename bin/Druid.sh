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
 -tpc.innerRadius 372.1 \
 -tpc.outerRadius 1692.1 \
 -tpc.halfZ 2225 \
 -draw.tpc.cylinder \
 -MCPtCut 0.0 \
 -logLevel 3 \
 "$@"


