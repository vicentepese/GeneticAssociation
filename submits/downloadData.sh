#!/bin/bash

ssh vipese@login.scg.stanford.edu
cd labs/mignot/PLATES_114_AND_UP/
files=`ls CHR*`

for i in "${files[@]}"
do 
    echo "$i"
done 