#!/bin/bash
files=`/bin/cat clist.txt`
echo $files
#splint +gnuextensions -weak -usedef -nestcomment -varuse -retvalother -fcnuse -fixedformalarray $files
#splint +gnuextensions -strict $files
splint $files
