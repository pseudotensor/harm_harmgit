#!/bin/bash

# http://stackoverflow.com/questions/67699/how-do-i-clone-all-remote-branches-with-git

#This will create tracking branches for all remote branches, except master (which you probably got from the original clone command).

for branch in `git branch -a | grep remotes | grep -v HEAD | grep -v master`; do
    git branch --track ${branch##*/} $branch
done


# I think you might still need to do a

git fetch --all
git pull --all
git fetch origin --tags
