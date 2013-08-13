#!/bin/bash
mymath=$1

echo " "
echo "pi: Good: `grep AGood $mymath |wc -l` Bad: `grep ABad $mymath |wc -l`  NEGS: `grep Aresultnegu $mymath |wc -l` `grep Aresultnegrho $mymath |wc -l` `grep AresultnegEr $mymath |wc -l`  COMPLX: `grep Aresultcomplexu $mymath |wc -l` `grep Aresultcomplexrho $mymath |wc -l` `grep AresultcomplexEr $mymath |wc -l`"

echo " "
echo "piS: Good: `grep ASGood $mymath |wc -l` Bad: `grep ASBad $mymath |wc -l`  NEGS: `grep ASresultnegu $mymath |wc -l` `grep ASresultnegrho $mymath |wc -l` `grep ASresultnegEr $mymath |wc -l`  COMPLX: `grep ASresultcomplexu $mymath |wc -l` `grep ASresultcomplexrho $mymath |wc -l` `grep ASresultcomplexEr $mymath |wc -l`"

echo " "
echo "Ui: Good: `grep 1Good $mymath |wc -l` Bad: `grep 1Bad $mymath |wc -l`  NEGS: `grep 1resultnegu $mymath |wc -l` `grep 1resultnegrho $mymath |wc -l` `grep 1resultnegEr $mymath |wc -l`  COMPLX: `grep 1resultcomplexu $mymath |wc -l` `grep 1resultcomplexrho $mymath |wc -l` `grep 1resultcomplexEr $mymath |wc -l`"

echo " "
echo "UiS: Good: `grep 1SGood $mymath |wc -l` Bad: `grep 1SBad $mymath |wc -l`  NEGS: `grep 1Sresultnegu $mymath |wc -l` `grep 1Sresultnegrho $mymath |wc -l` `grep 1SresultnegEr $mymath |wc -l`  COMPLX: `grep 1Sresultcomplexu $mymath |wc -l` `grep 1Sresultcomplexrho $mymath |wc -l` `grep 1SresultcomplexEr $mymath |wc -l`"

echo " "
echo "UU0: Good: `grep 2Good $mymath |wc -l` Bad: `grep 2Bad $mymath |wc -l`  NEGS: `grep 2resultnegu $mymath |wc -l` `grep 2resultnegrho $mymath |wc -l` `grep 2resultnegEr $mymath |wc -l`  COMPLX: `grep 2resultcomplexu $mymath |wc -l` `grep 2resultcomplexrho $mymath |wc -l` `grep 2resultcomplexEr $mymath |wc -l`"

echo " "
echo "UU0S: Good: `grep 0SnoGGood $mymath |wc -l` Bad: `grep 0SnoGBad $mymath |wc -l`  NEGS: `grep 0SnoGresultnegu $mymath |wc -l` `grep 0SnoGresultnegrho $mymath |wc -l` `grep 0SnoGresultnegEr $mymath |wc -l`  COMPLX: `grep 0SnoGresultcomplexu $mymath |wc -l` `grep 0SnoGresultcomplexrho $mymath |wc -l` `grep 0SnoGresultcomplexEr $mymath |wc -l`"


echo " "
echo "0: Good: `grep 0Good $mymath |wc -l` GoodW: `grep 0WGood $mymath |wc -l` Bad: `grep 0WBad $mymath |wc -l`  NEGS: `grep 0resultnegu $mymath |wc -l` `grep 0resultnegrho $mymath |wc -l` `grep 0resultnegEr $mymath |wc -l`  COMPLX: `grep 0resultcomplexu $mymath |wc -l` `grep 0resultcomplexrho $mymath |wc -l` `grep 0resultcomplexEr $mymath |wc -l`  NEGSW: `grep 0Wresultnegu $mymath |wc -l` `grep 0Wresultnegrho $mymath |wc -l` `grep 0WresultnegEr $mymath |wc -l`  COMPLXW: `grep 0Wresultcomplexu $mymath |wc -l` `grep 0Wresultcomplexrho $mymath |wc -l` `grep 0WresultcomplexEr $mymath |wc -l`"


echo " "
echo "0S: Good: `grep 0SGood $mymath |wc -l` GoodW: `grep 0WSGood $mymath |wc -l` Bad: `grep 0WSBad $mymath |wc -l`  NEGS: `grep 0Sresultnegu $mymath |wc -l` `grep 0Sresultnegrho $mymath |wc -l` `grep 0SresultnegEr $mymath |wc -l`  COMPLX: `grep 0Sresultcomplexu $mymath |wc -l` `grep 0Sresultcomplexrho $mymath |wc -l` `grep 0SresultcomplexEr $mymath |wc -l`   NEGSW: `grep 0WSresultnegu $mymath |wc -l` `grep 0WSresultnegrho $mymath |wc -l` `grep 0WSresultnegEr $mymath |wc -l`  COMPLXW: `grep 0WSresultcomplexu $mymath |wc -l` `grep 0WSresultcomplexrho $mymath |wc -l` `grep 0WSresultcomplexEr $mymath |wc -l`"


echo " "
echo "0M: Good: `grep 0MGood $mymath |wc -l` GoodW: `grep 0WSGood $mymath |wc -l` Bad: `grep 0WSBad $mymath |wc -l`  NEGS: `grep 0Mresultnegu $mymath |wc -l` `grep 0Mresultnegrho $mymath |wc -l` `grep 0MresultnegEr $mymath |wc -l`  COMPLX: `grep 0Mresultcomplexu $mymath |wc -l` `grep 0Mresultcomplexrho $mymath |wc -l` `grep 0MresultcomplexEr $mymath |wc -l`   NEGSW: `grep 0WSresultnegu $mymath |wc -l` `grep 0WSresultnegrho $mymath |wc -l` `grep 0WSresultnegEr $mymath |wc -l`  COMPLXW: `grep 0WSresultcomplexu $mymath |wc -l` `grep 0WSresultcomplexrho $mymath |wc -l` `grep 0WSresultcomplexEr $mymath |wc -l`"

echo " "
echo "0MS: Good: `grep 0MSGood $mymath |wc -l` GoodW: `grep 0WSGood $mymath |wc -l` Bad: `grep 0WSBad $mymath |wc -l`  NEGS: `grep 0MSresultnegu $mymath |wc -l` `grep 0MSresultnegrho $mymath |wc -l` `grep 0MSresultnegEr $mymath |wc -l`  COMPLX: `grep 0MSresultcomplexu $mymath |wc -l` `grep 0MSresultcomplexrho $mymath |wc -l` `grep 0MSresultcomplexEr $mymath |wc -l`   NEGSW: `grep 0WSresultnegu $mymath |wc -l` `grep 0WSresultnegrho $mymath |wc -l` `grep 0WSresultnegEr $mymath |wc -l`  COMPLXW: `grep 0WSresultcomplexu $mymath |wc -l` `grep 0WSresultcomplexrho $mymath |wc -l` `grep 0WSresultcomplexEr $mymath |wc -l`"

echo " "
echo "SHOULDS: `grep SHOULDUSEENTROPYNEGU $mymath | grep -v SHOULDUSEENTROPYNEGUBUTCANT |wc -l`  `grep SHOULDUSEENTROPYNEGRHO $mymath | grep -v SHOULDUSEENTROPYNEGRHOBUTCANT |wc -l`  `grep SHOULDUSEENTROPYNEGER $mymath | grep -v SHOULDUSEENTROPYNEGERBUTCANT |wc -l`"

echo " "
echo "SHOULDSBUTCANT: `grep SHOULDUSEENTROPYNEGUBUTCANT $mymath |wc -l`  `grep SHOULDUSEENTROPYNEGRHOBUTCANT $mymath |wc -l`  `grep SHOULDUSEENTROPYNEGERBUTCANT $mymath |wc -l`"

echo " "
echo "SHOULDCOND: `grep SHOULDUSEENTROPYCOND $mymath | grep -v SHOULDUSEENTROPYCONDBUTCANT |wc -l`"

echo " "
echo "SHOULDCONDBUTCANT: `grep SHOULDUSEENTROPYCONDBUTCANT $mymath |wc -l`"
