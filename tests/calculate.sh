#!/bin/bash
cd $1

jsons=$(ls *.json)
for file in $jsons
do 
    md5sums+=$(md5sum $file | awk '{print $1" "}')
done

md5sums+=$(echo $(awk '{$7=$8=$9=""}1' *.tsv | md5sum) | awk '{print $1}')

printf '%s\n' ${md5sums[*]} | sort


ls | echo $(grep -c "tsv") "tsv"
ls | echo $(grep -c "json") "json"