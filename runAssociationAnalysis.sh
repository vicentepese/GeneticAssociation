#!/bin/bash

chrpre=`jq '.chrArray.CHRanalysis' options.json` 
if [ -z "$chrpre" ] ||  [[ $chrpre == *"all"* ]]
then 
    echo "No chromosomes specified. Association analysis will be performed in all the chromosomes"
    cd submits
    sbatch --array="1-22" submitAssociationAnalysis.sh
else
    chr=()
    for i in $chrpre ; do
        if [[ "$i" != *"["* ]] && [[ "$i" != *"]"* ]]
            then
            chr+=$i
        fi
    done
    echo "The Association Analysis will be performed for chromosomes ${chr}"
    cd submits
    sbatch --array="$chr" submitAssociationAnalysis.sh 
fi

