#!/bin/bash
INPUT=pt_*
OUTPUT=combined.csv
cat $INPUT | head -n 1 > $OUTPUT
tail -n +2 -q $INPUT >> $OUTPUT