#!/bin/bash

formulae=( "asian_rs_fixed_call.py" "asian_vecer_fixed_call.py" "asian_dl_fixed_call.py" "asian_new_fixed_call.py" "asian_vecer_float_call.py" "asian_dl_float_call.py" "asian_new_float_call.py" )

printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
for i in "${formulae[@]}"
do
    echo $i
    python3 $i
    printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
done
