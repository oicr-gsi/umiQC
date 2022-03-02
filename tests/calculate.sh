#!/bin/bash
cd $1

for file in *
do 
    echo $(md5sum $file) | awk '{print $1}'
done

ls | echo $(grep -c "tsv") "tsv"
ls | echo $(grep -c "json") "json"