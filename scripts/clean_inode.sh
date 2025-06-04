#!/bin/bash
echo "Scanning for inotify watchers..."
for pid in $(find /proc/*/fd -lname anon_inode:inotify 2>/dev/null | cut -d/ -f3 | sort -u); do
    if ! ps -p $pid > /dev/null; then
        echo "PID $pid is not running. Skipping..."
        continue
    fi
    # Check if process is zombie or sleeping too long (optional logic)
    state=$(ps -o stat= -p $pid)
    if [[ $state == *Z* ]]; then
        echo "PID $pid is a zombie. Killing..."
        kill -9 $pid
    fi
done
