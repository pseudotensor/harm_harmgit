#!/bin/bash
# e.g. eps2gif.sh pic*.eps name.pic.gif
#
# deals with transparent background as white (useful for SM that inverts black-white)
#-density 200  (causes crash)
# -flatten  (causes no animation or crash)
#
#convert -flatten -density 200 -geometry 742x734 -delay 0 -loop 0 -background \#FFFFFF -dispose Background $1 $2
#
# e.g. 
# convert -flatten -density 200 -geometry 742x734 -background \#FFFFFF -dispose Background pic*.eps testnsbh_ffde_weight_oldbound.pic.gif

# ABOVE HAS ISSUES.  DO INSTEAD:

for fil in `ls pic*.eps`; do echo $fil ; convert -flatten -density 200 -geometry 742x734 -background \#FFFFFF -dispose Background $fil $fil.png ; done
convert *.png pic.gif

# convert *.png testnsbh_ffde_weight_oldbound.pic.gif
# convert *.png testnsbh_ffde_weight_fixbound_hasinsidefix_reverttest.pic.gif

