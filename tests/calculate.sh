#!/bin/bash
set -euo pipefail

cd $1

#find all json  files, return their md5sums to std out, list all file types
find . -name "*.json" -xtype f -exec sh -c "cat {} | md5sum" \;
ls | sed 's/.*\.//' | sort | uniq -c

