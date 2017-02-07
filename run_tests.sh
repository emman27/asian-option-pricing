#!/bin/bash

formulae=( "asian_rs_fixed_call.py" "asian_vecer_fixed_call.py" "asian_dl_fixed_call.py" "asian_vecer_float_call.py" )

for i in "${formulae[@]}"
do
    echo $i
    python3 $i
done
