#!/bin/bash

cd $1;

# set the right observatory (environmental variables)
cd $sdesy/vts/ed/; 
source setEnv_ed.sh; 
cd -;

echo "sdesy $sdesy"
echo "EVNDISPSYS $EVNDISPSYS"
echo "PYTHONPATH $PYTHONPATH"

python $1/calcLC.py $2 $3