#!/bin/sh
for fil in ls *.c *.h makefile; do
  echo $fil;
  sed -e 's/%g/%21.15g/g' -e 's/%21.15g/%21.15g/g' -e 's/%12.15g/%21.15g/g' -e 's/%15.10g/%21.15g/g' -e 's/%10.5g/%21.15g/g' -e 's/%26.20g/%21.15g/g' -e 's/%26.20e/%21.15g/g' -e 's/%20.10g/%21.15g/g' -e 's/%25.17g/%21.15g/g' -e 's/%15.8g/%21.15g/g'  $fil > $fil.new
  #
  mv $fil.new $fil
  rm $fil.new
done
