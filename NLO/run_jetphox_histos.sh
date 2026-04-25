#!/bin/bash
cd "$(dirname "$0")"

root -b -q 'MakeJetPHOXhisto.C("05")' > logs/jetPHOX_05.log 2>&1 &
root -b -q 'MakeJetPHOXhisto.C("10")' > logs/jetPHOX_10.log 2>&1 &
root -b -q 'MakeJetPHOXhisto.C("20")' > logs/jetPHOX_20.log 2>&1 &

wait
echo "All done."
