#!/bin/bash
cd $1

jsons=$(ls *.json)
for file in $jsons
do 
    echo $(md5sum $file) | awk '{print $1}'
done

echo $(awk '{$7=""}1' *.tsv | md5sum) | awk '{print $1}'

ls | echo $(grep -c "tsv") "tsv"
ls | echo $(grep -c "json") "json"