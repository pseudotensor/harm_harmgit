(ssh $1 "cd `pwd` ; ls ; scp bh:*.sh . ; rm *.fli ; sh mkflijon.sh $2 $3" > /dev/null 2>&1 &)
