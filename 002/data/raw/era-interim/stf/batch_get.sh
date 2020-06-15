#!/bin/sh

for yr in {2000..2012}; do
    python get_interim_stf.py $yr
done
