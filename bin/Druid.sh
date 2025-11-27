#!/usr/bin/env bash

#source /cvmfs/sw.hsf.org/key4hep/setup.sh
export LIBGL_ALWAYS_SOFTWARE=1
#export LIBGL_ALWAYS_INDIRECT=1


# Absolute path to the directory that this script lives in
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Run the sibling executable "Druid"
"$SCRIPT_DIR/Druid" -BField 3.5 -coll.caloHit.filterOutSuffix Digi "$@"

