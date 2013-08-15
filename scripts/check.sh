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
echo "0: Good: `grep 0Good $mymath |wc -l` GoodW: `grep 0WGood $mymath |wc -l` Bad: `grep 0Bad $mymath |wc -l`  WBad: `grep 0WBad $mymath |wc -l`  NEGS: `grep 0resultnegu $mymath |wc -l` `grep 0resultnegrho $mymath |wc -l` `grep 0resultnegEr $mymath |wc -l`  COMPLX: `grep 0resultcomplexu $mymath |wc -l` `grep 0resultcomplexrho $mymath |wc -l` `grep 0resultcomplexEr $mymath |wc -l`  NEGSW: `grep 0Wresultnegu $mymath |wc -l` `grep 0Wresultnegrho $mymath |wc -l` `grep 0WresultnegEr $mymath |wc -l`  COMPLXW: `grep 0Wresultcomplexu $mymath |wc -l` `grep 0Wresultcomplexrho $mymath |wc -l` `grep 0WresultcomplexEr $mymath |wc -l`"


echo " "
echo "0S: Good: `grep 0SGood $mymath |wc -l` GoodW: `grep 0WSGood $mymath |wc -l` Bad: `grep 0SBad $mymath |wc -l`  WBad: `grep 0WSBad $mymath |wc -l`  NEGS: `grep 0Sresultnegu $mymath |wc -l` `grep 0Sresultnegrho $mymath |wc -l` `grep 0SresultnegEr $mymath |wc -l`  COMPLX: `grep 0Sresultcomplexu $mymath |wc -l` `grep 0Sresultcomplexrho $mymath |wc -l` `grep 0SresultcomplexEr $mymath |wc -l`   NEGSW: `grep 0WSresultnegu $mymath |wc -l` `grep 0WSresultnegrho $mymath |wc -l` `grep 0WSresultnegEr $mymath |wc -l`  COMPLXW: `grep 0WSresultcomplexu $mymath |wc -l` `grep 0WSresultcomplexrho $mymath |wc -l` `grep 0WSresultcomplexEr $mymath |wc -l`"


echo " "
echo "0M: Good: `grep 0MGood $mymath |wc -l` GoodW: `grep 0WMGood $mymath |wc -l` Bad: `grep 0MBad $mymath |wc -l`  NEGS: `grep 0Mresultnegu $mymath |wc -l` `grep 0Mresultnegrho $mymath |wc -l` `grep 0MresultnegEr $mymath |wc -l`  COMPLX: `grep 0Mresultcomplexu $mymath |wc -l` `grep 0Mresultcomplexrho $mymath |wc -l` `grep 0MresultcomplexEr $mymath |wc -l`   NEGSW: `grep 0WMresultnegu $mymath |wc -l` `grep 0WMresultnegrho $mymath |wc -l` `grep 0WMresultnegEr $mymath |wc -l`  COMPLXW: `grep 0WMresultcomplexu $mymath |wc -l` `grep 0WMresultcomplexrho $mymath |wc -l` `grep 0WMresultcomplexEr $mymath |wc -l`"

echo " "
echo "0MS: Good: `grep 0MSGood $mymath |wc -l` GoodW: `grep 0WMSGood $mymath |wc -l` Bad: `grep 0MSBad $mymath |wc -l`  NEGS: `grep 0MSresultnegu $mymath |wc -l` `grep 0MSresultnegrho $mymath |wc -l` `grep 0MSresultnegEr $mymath |wc -l`  COMPLX: `grep 0MSresultcomplexu $mymath |wc -l` `grep 0MSresultcomplexrho $mymath |wc -l` `grep 0MSresultcomplexEr $mymath |wc -l`   NEGSW: `grep 0WMSresultnegu $mymath |wc -l` `grep 0WMSresultnegrho $mymath |wc -l` `grep 0WMSresultnegEr $mymath |wc -l`  COMPLXW: `grep 0WMSresultcomplexu $mymath |wc -l` `grep 0WMSresultcomplexrho $mymath |wc -l` `grep 0WMSresultcomplexEr $mymath |wc -l`"


echo " "
echo "0Q: Good: `grep 0QGood $mymath |wc -l` GoodW: `grep 0WQGood $mymath |wc -l` Bad: `grep 0QBad $mymath |wc -l`  NEGS: `grep 0Qresultnegu $mymath |wc -l` `grep 0Qresultnegrho $mymath |wc -l` `grep 0QresultnegEr $mymath |wc -l`  COMPLX: `grep 0Qresultcomplexu $mymath |wc -l` `grep 0Qresultcomplexrho $mymath |wc -l` `grep 0QresultcomplexEr $mymath |wc -l`   NEGSW: `grep 0WQresultnegu $mymath |wc -l` `grep 0WQresultnegrho $mymath |wc -l` `grep 0WQresultnegEr $mymath |wc -l`  COMPLXW: `grep 0WQresultcomplexu $mymath |wc -l` `grep 0WQresultcomplexrho $mymath |wc -l` `grep 0WQresultcomplexEr $mymath |wc -l`"

echo " "
echo "0QS: Good: `grep 0QSGood $mymath |wc -l` GoodW: `grep 0WQSGood $mymath |wc -l` Bad: `grep 0QSBad $mymath |wc -l`  NEGS: `grep 0QSresultnegu $mymath |wc -l` `grep 0QSresultnegrho $mymath |wc -l` `grep 0QSresultnegEr $mymath |wc -l`  COMPLX: `grep 0QSresultcomplexu $mymath |wc -l` `grep 0QSresultcomplexrho $mymath |wc -l` `grep 0QSresultcomplexEr $mymath |wc -l`   NEGSW: `grep 0WQSresultnegu $mymath |wc -l` `grep 0WQSresultnegrho $mymath |wc -l` `grep 0WQSresultnegEr $mymath |wc -l`  COMPLXW: `grep 0WQSresultcomplexu $mymath |wc -l` `grep 0WQSresultcomplexrho $mymath |wc -l` `grep 0WQSresultcomplexEr $mymath |wc -l`"

echo " "
echo "SHOULDS: `grep SHOULDUSEENTROPYNEGU $mymath | grep -v SHOULDUSEENTROPYNEGUBUTCANT |wc -l`  `grep SHOULDUSEENTROPYNEGRHO $mymath | grep -v SHOULDUSEENTROPYNEGRHOBUTCANT |wc -l`  `grep SHOULDUSEENTROPYNEGER $mymath | grep -v SHOULDUSEENTROPYNEGERBUTCANT |wc -l`"

echo " "
echo "SHOULDSBUTCANT: `grep SHOULDUSEENTROPYNEGUBUTCANT $mymath |wc -l`  `grep SHOULDUSEENTROPYNEGRHOBUTCANT $mymath |wc -l`  `grep SHOULDUSEENTROPYNEGERBUTCANT $mymath |wc -l`"

echo " "
echo "SHOULDCOND: `grep SHOULDUSEENTROPYCOND $mymath | grep -v SHOULDUSEENTROPYCONDBUTCANT |wc -l`"

echo " "
echo "SHOULDCONDBUTCANT: `grep SHOULDUSEENTROPYCONDBUTCANT $mymath |wc -l`"

echo " "
echo "AllBad: `grep AllBad $mymath | grep -v NotAllBad |wc -l`"

echo " "
echo "OneGood: `grep OneGood $mymath | grep -v NotOneGood |wc -l`"

echo " "
echo "OneActualGood: `grep OneActualGood $mymath | grep -v NotOneActualGood |wc -l`"


echo " "
echo "NotAllBad: `grep NotAllBad $mymath |wc -l`"

echo " "
echo "NotOneGood: `grep NotOneGood $mymath |wc -l`"

echo " "
echo "NotOneActualGood: `grep NotOneActualGood $mymath |wc -l`"
