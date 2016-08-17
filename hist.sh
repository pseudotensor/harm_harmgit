#!/bin/sh

# https://help.github.com/articles/changing-author-info/

#OLD_NAME="Jonathan C. McKinney"
#CORRECT_NAME="Jonathan C. McKinney"

git filter-branch -f --env-filter '
OLD_NAME="Jonathan McKinney"
CORRECT_NAME="Jonathan C. McKinney"
CORRECT_EMAIL="pseudotensor@gmail.com"
if [ "$GIT_COMMITTER_NAME" = "$OLD_NAME" ]
then
    export GIT_COMMITTER_NAME="$CORRECT_NAME"
    export GIT_COMMITTER_EMAIL="$CORRECT_EMAIL"
fi
if [ "$GIT_AUTHOR_NAME" = "$OLD_NAME" ]
then
    export GIT_AUTHOR_NAME="$CORRECT_NAME"
    export GIT_AUTHOR_EMAIL="$CORRECT_EMAIL"
fi
' --tag-name-filter cat -- --branches --tags


# git push --force --tags origin 'refs/heads/*'
