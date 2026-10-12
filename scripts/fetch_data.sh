#!/bin/sh
# Download the input graphs.
#   sh scripts/fetch_data.sh        the five SSTD 2005 road networks (required)
#   sh scripts/fetch_data.sh snap   additionally six graphs of the SNAP collection (about 200 MB)
# Sources: https://github.com/mehdi-ansari/Real-dataset-shortest-path-sample
#          (a mirror of https://users.cs.utah.edu/~lifeifei/SpatialDataset.htm)
#          https://snap.stanford.edu/data/
set -e
cd "$(dirname "$0")/.."
mkdir -p data
if [ ! -d data/sstd2005 ]; then
  git clone https://github.com/mehdi-ansari/Real-dataset-shortest-path-sample.git data/sstd2005
fi
# the commit the published measurements were taken from
git -C data/sstd2005 checkout -q 99cc5c06982d9f5345100c1d1f88d9ccdeb1adcf
ls "data/sstd2005/Datasets - SSTD 2005"

if [ "$1" = "snap" ]; then
  mkdir -p data/snap
  for f in web-Stanford web-Google roadNet-PA roadNet-TX roadNet-CA soc-pokec-relationships; do
    [ -f "data/snap/$f.txt.gz" ] || curl -fSL -o "data/snap/$f.txt.gz" "https://snap.stanford.edu/data/$f.txt.gz"
  done
  ls -l data/snap
fi
