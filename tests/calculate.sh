#!/bin/bash
set -euo pipefail

cd $1

find . -name "*.json" -xtype f -exec sh -c "cat {} | md5sum |sort" \;
ls | sed 's/.*\.//' | sort | uniq -c
