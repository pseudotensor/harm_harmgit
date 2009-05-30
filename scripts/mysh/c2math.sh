#sed 's/\([0-9]\)e\([+,-]\)\([0-9]\+\)/\1*10^(\2\3)/g' $1 | sed 's/\*10^(\([+-][0-9]\+\))``20/``20*10^(\1)/g' > $2
#
# double to long double code
for fil in *.c *.h; do echo $fil; sed 's/\(%[0-9]*\.[0-9]*\)g/%31.25Lg/g' $fil | sed 's/M_PI /M_PIl /g' | sed 's/``20/``30/g' | sed 's/REALTYPE DOUBLETYPE/REALTYPE LONGDOUBLETYPE/g'| sed 's/SENSITIVE DOUBLETYPE/SENSITIVE LONGDOUBLETYPE/g'    > $1/$fil; done

