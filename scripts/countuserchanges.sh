#!/bin/bash
svn ls -R | egrep -v -e "\/$" | xargs svn blame | awk '{print $2}' | sort | uniq -c | sort -r

