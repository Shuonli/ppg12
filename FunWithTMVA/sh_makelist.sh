#!/usr/bin/env bash
for a in {0..0}; do
  for b in {0..0}; do
    for c in {0..3}; do
      echo "$a $b $c"
    done
  done
done > triples.txt
