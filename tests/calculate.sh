#!/bin/bash
set -euo pipefail

cd $1

find . -name "output*" -xtype f -exec sh -c "cat {} | md5sum" \;
ls | sed 's/.*\.//' | sort | uniq -c

