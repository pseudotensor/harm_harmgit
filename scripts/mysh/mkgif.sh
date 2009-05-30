mkdir rases ; cd rases
for fil in ../images/im0f1s0l[0-9]*.r8
  do
  echo converting $i
  r8torasjon 0  /home/jon/research/current/bin/i/genimages.pal $fil
done
mv ../images/*.ras .
convert im0f1s0l[0-9]*.ras im0f1s0l.gif
