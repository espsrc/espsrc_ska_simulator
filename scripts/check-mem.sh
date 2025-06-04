#!/bin/bash
echo "$(date '+%Y-%m-%d;%H:%M:%S');$(free | awk '/^Mem:/ {print $3}')"
