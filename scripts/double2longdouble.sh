#!/bin/sh
for fil in `ls *.c *.h makefile`; do
  echo $fil;
  #diff -bBdp $fil $1/$fil;
  sed -e 's/%g/%Lg/g' -e 's/%21.15g/%26.20Lg/g' -e 's/%31.26g/%31.26Lg/g' -e 's/%10.2g/%10.2Lg/g' -e 's/%12.15g/%26.20Lg/g' -e 's/%15.10g/%26.20Lg/g' -e 's/%10.5g/%26.20Lg/g' -e 's/%26.20g/%26.20Lg/g' -e 's/%26.20e/%26.20Lg/g' -e 's/%20.10g/%26.20Lg/g' -e 's/%25.17g/%26.20Lg/g' -e 's/%15.8g/%26.20Lg/g' -e 's/%lf/%Lf/g' -e 's/sqrt(/sqrtl(/g' -e 's/log10(/log10l(/g' -e 's/log(/logl(/g' -e 's/exp(/expl(/g' -e 's/fabs(/fabsl(/g' -e 's/sin(/sinl(/g' -e 's/ceil(/ceill(/g' -e 's/isfinite(/isfinitel(/g' -e 's/round(/roundl(/g' -e 's/floor(/floorl(/g' -e 's/lrint(/lrintl(/g' -e 's/fmod(/fmodl(/g' -e 's/modf(/modfl(/g' -e 's/cot(/cotl(/g' -e 's/cbrt(/cbrtl(/g' -e 's/cos(/cosl(/g' -e 's/tan(/tanl(/g' -e 's/atan(/atanl(/g' -e 's/atan2(/atan2l(/g' -e 's/pow(/powl(/g' -e 's/M_PI/M_PIl/g' -e 's/M_PIll/M_PIl/g' -e 's/#define CTYPE double/#define CTYPE long double/g' -e 's/#define MPI_CTYPE MPI_DOUBLE/#define MPI_CTYPE MPI_LONG_DOUBLE/g' -e 's/#define REALTYPE DOUBLETYPE/#define REALTYPE LONGDOUBLETYPE/g' -e 's/#define SENSITIVE DOUBLETYPE/#define SENSITIVE LONGDOUBLETYPE/g' -e 's/#define SUPERLONGDOUBLE 0/#define SUPERLONGDOUBLE 1/g'  $fil > $fil.new
  #
  sed -e 's/\([0-9.]\+[eE]\+[+-]*[0-9]\+\)/\1L/g' -e 's/\([0-9]\+\.[0-9]\+\)\([ ]\|[)]\|[\;]\|[\*]\|[+-]\|[=]\|$\|[/]\)/\1L\2/g' -e 's/\([0-9]\+\.\)\([ ]\|[)]\|[\;]\|[\*]\|[+-]\|[=]\|$\|[/]\)/\1L\2/g' -e 's/\([^0-9]\)\(\.[0-9]\+\)/\1\2L/g' -e 's/\(\.[0-9]\)L\([eE]\+[0-9]\+L\)/\1\2/g' $fil.new > $fil.new2
  mv $fil.new2 $fil
  rm $fil.new
done
#sed -e 's/2\.0/2\.0L/g' -e 's/1\.25/1\.25L/g' -e 's/3\.0/3\.0L/g' -e 's/6\.0/6\.0L/g'  -e 's/8\.0/8\.0L/g' -e 's/0\.5/0\.5L/g' interpbackup > interp.c

#define SUPERLONGDOUBLE 0

# 1.134e-19LL

# ptr8 and MPI_DOUBLE in init_mpi
# 

#  sed -e 's/\.\([0-9]\+\)/.\1L/g' -e 's/\([0-9]\+\.[0-9]\+\)/\1L/g' -e 's/\([0-9]\+\)\./\1.L/g'


#  sed -e  's/\([+-]?\([0-9]\+[.]?[0-9]*|[.][0-9]+\)\([eE][+-]?[0-9]+\)?\)/\1L/g'

#  sed -e  's/\([+-]?([0-9]\+[.]?[0-9]*|[.][0-9]\+)([eE][+-]?[0-9]\+)?\)/\1L/g' test.c


#  sed -e 's/\([0-9]\+\.[0-9]*\([eE][-+]?[0-9]+\)\)/\1L/g' -e 's/\([0-9]*\.[0-9]\+\)/\1L/g' 

#  sed -e 's/\([0-9]\+\.[0-9]*\)/\1L/g' -e 's/\([0-9]*\.[0-9]\+\)/\1L/g' 

#  sed -e 's/\([0-9]\+\.[0-9]*\([eE][-+]?[0-9]+\)\)/\1L/g' 

#  sed -e 's/\([0-9]\+\.*[+-]*[eE]\+[0-9]\+\)/\1L/g' -e 's/\([0-9]\+\.[0-9]\+\)/\1L/g' -e 's/\([0-9]\+\.\)\([+-\;]\)\+/\1L\2/g'

#  sed -e 's/\([0-9]\+\.*[+-]*[eE]\+[0-9]\+\)/\1L/g' -e 's/\([0-9]\+\.[0-9]\+\)/\1L/g' -e 's/\([0-9]\+\.\)/\1L/g'

#sed -e 's/\([0-9.]\+[eE]\+[+-]*[0-9]\+\)/\1L/g' -e 's/\([0-9]\+\.[0-9]\+\)/\1L/g' -e 's/\([0-9]\+\.\)\(;\|[+-]\|\*\|[/]\)\+/\1L\2/g' test.c


#sed -e 's/\([0-9.]\+[eE]\+[+-]*[0-9]\+\)/\1L/g' -e 's/\([0-9]\+\.[0-9]\+\)\([ ]\|[\;]\|[\*]\|[+-]\|$\|[/]\&[^0-9]\)\?/\1L\2/g' -e 's/\([0-9]\+\.\)\([ ]\|[\;]\|[\*]\|[+-]\|$\|[/]\&[^0-9]\)\?/\1L\2/g'  test.c

#sed -e 's/\([0-9.]\+[eE]\+[+-]*[0-9]\+\)/\1L/g' -e 's/\([0-9]\+\.[0-9]\+\)\([ ]\|[\;]\|[\*]\|[+-]\|$\|[/]\&[^0-9]\)\?/\1L\2/g' -e 's/\([0-9]\+\.\)\([ ]\|[\;]\|[\*]\|[+-]\|$\|[/]\&[^0-9]\)\?/\1L\2/g' -e 's/\([^0-9]\)\+\(\.[0-9]\+\)/\1\2L/g' test.c

#sed -e 's/\([0-9.]\+[eE]\+[+-]*[0-9]\+\)/\1L/g' -e 's/\([0-9]\+\.[0-9]\+\)\([ ]\|[\;]\|[\*]\|[+-]\|$\|[/]\&[^0-9]\)\?/\1L\2/g' -e 's/\([0-9]\+\.\)\([ ]\|[\;]\|[\*]\|[+-]\|$\|[/]\&[^0-9]\)\?/\1L\2/g' test.c

#sed -e 's/\([0-9]\+\.*[0-9]*[eE]\+[+-]\+[0-9]\+\)/\1L/g' test.c


#sed -e 's/\([0-9.]\+[eE]\+[+-]*[0-9]\+\)/\1L/g' -e 's/\([0-9]\+\.[0-9]\+\)\([ ]\|[\;]\|[\*]\|[+-]\|$\|[/]\)/\1L\2/g' -e 's/\([0-9]\+\.\)\([ ]\|[\;]\|[\*]\|[+-]\|$\|[/]\)/\1L\2/g' -e 's/\([^0-9]\)\+\(\.[0-9]\+\)/\1\2L/g' -e 's/\(\.[0-9]\)L\([eE]\+[0-9]\+L\)/\1\2/g' 


