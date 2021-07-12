#!/bin/bash

for wy in {1995..2009}; do 
        wy_p1=$((${wy}  + 1))
        start_date=${wy}-10-01
        end_date=${wy_p1}-09-30
        python inference.py $start_date $end_date
done
