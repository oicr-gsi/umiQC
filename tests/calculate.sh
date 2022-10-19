#!/bin/bash
set -euo pipefail

cd $1

ls | sed 's/.*\.//' | sort | uniq -c
