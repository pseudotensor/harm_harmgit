#!/bin/bash

# setting up for dependency list
# http://beige.ucs.indiana.edu/I590/node45.html
#

sed -e 's/firstrun=1/firstrun=0/g' -e 's/rada0.94 - a/rada0.94 a b/g' -e 's/rada0.94a/rada0.94b/g' batch.qsub.kraken.rada0.94a > batch.qsub.kraken.rada0.94b
                        
sed -e 's/rada0.94 a b/rada0.94 b c/g' -e 's/rada0.94b/rada0.94c/g' batch.qsub.kraken.rada0.94b > batch.qsub.kraken.rada0.94c
sed -e 's/rada0.94 b c/rada0.94 c d/g' -e 's/rada0.94c/rada0.94d/g' batch.qsub.kraken.rada0.94c > batch.qsub.kraken.rada0.94d
sed -e 's/rada0.94 c d/rada0.94 d e/g' -e 's/rada0.94d/rada0.94e/g' batch.qsub.kraken.rada0.94d > batch.qsub.kraken.rada0.94e
sed -e 's/rada0.94 d e/rada0.94 e f/g' -e 's/rada0.94e/rada0.94f/g' batch.qsub.kraken.rada0.94e > batch.qsub.kraken.rada0.94f
sed -e 's/rada0.94 e f/rada0.94 f g/g' -e 's/rada0.94f/rada0.94g/g' batch.qsub.kraken.rada0.94f > batch.qsub.kraken.rada0.94g
sed -e 's/rada0.94 f g/rada0.94 g h/g' -e 's/rada0.94g/rada0.94h/g' batch.qsub.kraken.rada0.94g > batch.qsub.kraken.rada0.94h
sed -e 's/rada0.94 g h/rada0.94 h i/g' -e 's/rada0.94h/rada0.94i/g' batch.qsub.kraken.rada0.94h > batch.qsub.kraken.rada0.94i
sed -e 's/rada0.94 h i/rada0.94 i j/g' -e 's/rada0.94i/rada0.94j/g' batch.qsub.kraken.rada0.94i > batch.qsub.kraken.rada0.94j
sed -e 's/rada0.94 i j/rada0.94 j k/g' -e 's/rada0.94j/rada0.94k/g' batch.qsub.kraken.rada0.94j > batch.qsub.kraken.rada0.94k
sed -e 's/rada0.94 j k/rada0.94 k l/g' -e 's/rada0.94k/rada0.94l/g' batch.qsub.kraken.rada0.94k > batch.qsub.kraken.rada0.94l
sed -e 's/rada0.94 k l/rada0.94 l m/g' -e 's/rada0.94l/rada0.94m/g' batch.qsub.kraken.rada0.94l > batch.qsub.kraken.rada0.94m
sed -e 's/rada0.94 l m/rada0.94 m n/g' -e 's/rada0.94m/rada0.94n/g' batch.qsub.kraken.rada0.94m > batch.qsub.kraken.rada0.94n
sed -e 's/rada0.94 m n/rada0.94 n o/g' -e 's/rada0.94n/rada0.94o/g' batch.qsub.kraken.rada0.94n > batch.qsub.kraken.rada0.94o
sed -e 's/rada0.94 n o/rada0.94 o p/g' -e 's/rada0.94o/rada0.94p/g' batch.qsub.kraken.rada0.94o > batch.qsub.kraken.rada0.94p
sed -e 's/rada0.94 o p/rada0.94 p q/g' -e 's/rada0.94p/rada0.94q/g' batch.qsub.kraken.rada0.94p > batch.qsub.kraken.rada0.94q
sed -e 's/rada0.94 p q/rada0.94 q r/g' -e 's/rada0.94q/rada0.94r/g' batch.qsub.kraken.rada0.94q > batch.qsub.kraken.rada0.94r
sed -e 's/rada0.94 q r/rada0.94 r s/g' -e 's/rada0.94r/rada0.94s/g' batch.qsub.kraken.rada0.94r > batch.qsub.kraken.rada0.94s
sed -e 's/rada0.94 r s/rada0.94 s t/g' -e 's/rada0.94s/rada0.94t/g' batch.qsub.kraken.rada0.94s > batch.qsub.kraken.rada0.94t
sed -e 's/rada0.94 s t/rada0.94 t u/g' -e 's/rada0.94t/rada0.94u/g' batch.qsub.kraken.rada0.94t > batch.qsub.kraken.rada0.94u
sed -e 's/rada0.94 t u/rada0.94 u v/g' -e 's/rada0.94u/rada0.94v/g' batch.qsub.kraken.rada0.94u > batch.qsub.kraken.rada0.94v
sed -e 's/rada0.94 u v/rada0.94 v w/g' -e 's/rada0.94v/rada0.94w/g' batch.qsub.kraken.rada0.94v > batch.qsub.kraken.rada0.94w
sed -e 's/rada0.94 v w/rada0.94 w x/g' -e 's/rada0.94w/rada0.94x/g' batch.qsub.kraken.rada0.94w > batch.qsub.kraken.rada0.94x
sed -e 's/rada0.94 w x/rada0.94 x y/g' -e 's/rada0.94x/rada0.94y/g' batch.qsub.kraken.rada0.94x > batch.qsub.kraken.rada0.94y
sed -e 's/rada0.94 x y/rada0.94 y z/g' -e 's/rada0.94y/rada0.94z/g' batch.qsub.kraken.rada0.94y > batch.qsub.kraken.rada0.94z
	    
sed -e 's/rada0.94 y z/rada0.94 z nexta/g' -e 's/rada0.94z/rada0.94nexta/g' batch.qsub.kraken.rada0.94z > batch.qsub.kraken.rada0.94nexta

sed -e 's/rada0.94 z nexta/rada0.94 nexta nextb/g' -e 's/rada0.94nexta/rada0.94nextb/g' batch.qsub.kraken.rada0.94nexta > batch.qsub.kraken.rada0.94nextb

sed -e 's/rada0.94 nexta nextb/rada0.94 nextb nextc/g' -e 's/rada0.94nextb/rada0.94nextc/g' batch.qsub.kraken.rada0.94nextb > batch.qsub.kraken.rada0.94nextc
sed -e 's/rada0.94 nextb nextc/rada0.94 nextc nextd/g' -e 's/rada0.94nextc/rada0.94nextd/g' batch.qsub.kraken.rada0.94nextc > batch.qsub.kraken.rada0.94nextd
sed -e 's/rada0.94 nextc nextd/rada0.94 nextd nexte/g' -e 's/rada0.94nextd/rada0.94nexte/g' batch.qsub.kraken.rada0.94nextd > batch.qsub.kraken.rada0.94nexte
sed -e 's/rada0.94 nextd nexte/rada0.94 nexte nextf/g' -e 's/rada0.94nexte/rada0.94nextf/g' batch.qsub.kraken.rada0.94nexte > batch.qsub.kraken.rada0.94nextf
sed -e 's/rada0.94 nexte nextf/rada0.94 nextf nextg/g' -e 's/rada0.94nextf/rada0.94nextg/g' batch.qsub.kraken.rada0.94nextf > batch.qsub.kraken.rada0.94nextg
sed -e 's/rada0.94 nextf nextg/rada0.94 nextg nexth/g' -e 's/rada0.94nextg/rada0.94nexth/g' batch.qsub.kraken.rada0.94nextg > batch.qsub.kraken.rada0.94nexth
sed -e 's/rada0.94 nextg nexth/rada0.94 nexth nexti/g' -e 's/rada0.94nexth/rada0.94nexti/g' batch.qsub.kraken.rada0.94nexth > batch.qsub.kraken.rada0.94nexti
sed -e 's/rada0.94 nexth nexti/rada0.94 nexti nextj/g' -e 's/rada0.94nexti/rada0.94nextj/g' batch.qsub.kraken.rada0.94nexti > batch.qsub.kraken.rada0.94nextj
sed -e 's/rada0.94 nexti nextj/rada0.94 nextj nextk/g' -e 's/rada0.94nextj/rada0.94nextk/g' batch.qsub.kraken.rada0.94nextj > batch.qsub.kraken.rada0.94nextk
sed -e 's/rada0.94 nextj nextk/rada0.94 nextk nextl/g' -e 's/rada0.94nextk/rada0.94nextl/g' batch.qsub.kraken.rada0.94nextk > batch.qsub.kraken.rada0.94nextl
sed -e 's/rada0.94 nextk nextl/rada0.94 nextl nextm/g' -e 's/rada0.94nextl/rada0.94nextm/g' batch.qsub.kraken.rada0.94nextl > batch.qsub.kraken.rada0.94nextm
sed -e 's/rada0.94 nextl nextm/rada0.94 nextm nextn/g' -e 's/rada0.94nextm/rada0.94nextn/g' batch.qsub.kraken.rada0.94nextm > batch.qsub.kraken.rada0.94nextn
sed -e 's/rada0.94 nextm nextn/rada0.94 nextn nexto/g' -e 's/rada0.94nextn/rada0.94nexto/g' batch.qsub.kraken.rada0.94nextn > batch.qsub.kraken.rada0.94nexto
sed -e 's/rada0.94 nextn nexto/rada0.94 nexto nextp/g' -e 's/rada0.94nexto/rada0.94nextp/g' batch.qsub.kraken.rada0.94nexto > batch.qsub.kraken.rada0.94nextp
sed -e 's/rada0.94 nexto nextp/rada0.94 nextp nextq/g' -e 's/rada0.94nextp/rada0.94nextq/g' batch.qsub.kraken.rada0.94nextp > batch.qsub.kraken.rada0.94nextq
sed -e 's/rada0.94 nextp nextq/rada0.94 nextq nextr/g' -e 's/rada0.94nextq/rada0.94nextr/g' batch.qsub.kraken.rada0.94nextq > batch.qsub.kraken.rada0.94nextr
sed -e 's/rada0.94 nextq nextr/rada0.94 nextr nexts/g' -e 's/rada0.94nextr/rada0.94nexts/g' batch.qsub.kraken.rada0.94nextr > batch.qsub.kraken.rada0.94nexts
sed -e 's/rada0.94 nextr nexts/rada0.94 nexts nextt/g' -e 's/rada0.94nexts/rada0.94nextt/g' batch.qsub.kraken.rada0.94nexts > batch.qsub.kraken.rada0.94nextt
sed -e 's/rada0.94 nexts nextt/rada0.94 nextt nextu/g' -e 's/rada0.94nextt/rada0.94nextu/g' batch.qsub.kraken.rada0.94nextt > batch.qsub.kraken.rada0.94nextu
sed -e 's/rada0.94 nextt nextu/rada0.94 nextu nextv/g' -e 's/rada0.94nextu/rada0.94nextv/g' batch.qsub.kraken.rada0.94nextu > batch.qsub.kraken.rada0.94nextv
sed -e 's/rada0.94 nextu nextv/rada0.94 nextv nextw/g' -e 's/rada0.94nextv/rada0.94nextw/g' batch.qsub.kraken.rada0.94nextv > batch.qsub.kraken.rada0.94nextw
sed -e 's/rada0.94 nextv nextw/rada0.94 nextw nextx/g' -e 's/rada0.94nextw/rada0.94nextx/g' batch.qsub.kraken.rada0.94nextw > batch.qsub.kraken.rada0.94nextx
sed -e 's/rada0.94 nextw nextx/rada0.94 nextx nexty/g' -e 's/rada0.94nextx/rada0.94nexty/g' batch.qsub.kraken.rada0.94nextx > batch.qsub.kraken.rada0.94nexty
sed -e 's/rada0.94 nextx nexty/rada0.94 nexty nextz/g' -e 's/rada0.94nexty/rada0.94nextz/g' batch.qsub.kraken.rada0.94nexty > batch.qsub.kraken.rada0.94nextz
