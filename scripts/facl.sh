#!/bin/bash
# to get all in group to have rwx access for current and any future files in the tree
setfacl -Rm g:users:rwX,d:g:users:rwX . ; setfacl -Rm o:rwX . ;  chmod -R a+rwx .
