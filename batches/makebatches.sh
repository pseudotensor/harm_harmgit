#!/bin/bash

# setting up for dependency list
# http://beige.ucs.indiana.edu/I590/node45.html
#

#sed -e 's/rada0 - a/rada0 a b/g' -e 's/rada0a/rada0b/g' batch.qsub.kraken.rada0a > batch.qsub.kraken.rada0b
                         
sed -e 's/rada0 a b/rada0 b c/g' -e 's/rada0b/rada0c/g' batch.qsub.kraken.rada0b > batch.qsub.kraken.rada0c
sed -e 's/rada0 b c/rada0 c d/g' -e 's/rada0c/rada0d/g' batch.qsub.kraken.rada0c > batch.qsub.kraken.rada0d
sed -e 's/rada0 c d/rada0 d e/g' -e 's/rada0d/rada0e/g' batch.qsub.kraken.rada0d > batch.qsub.kraken.rada0e
sed -e 's/rada0 d e/rada0 e f/g' -e 's/rada0e/rada0f/g' batch.qsub.kraken.rada0e > batch.qsub.kraken.rada0f
sed -e 's/rada0 e f/rada0 f g/g' -e 's/rada0f/rada0g/g' batch.qsub.kraken.rada0f > batch.qsub.kraken.rada0g
sed -e 's/rada0 f g/rada0 g h/g' -e 's/rada0g/rada0h/g' batch.qsub.kraken.rada0g > batch.qsub.kraken.rada0h
sed -e 's/rada0 g h/rada0 h i/g' -e 's/rada0h/rada0i/g' batch.qsub.kraken.rada0h > batch.qsub.kraken.rada0i
sed -e 's/rada0 h i/rada0 i j/g' -e 's/rada0i/rada0j/g' batch.qsub.kraken.rada0i > batch.qsub.kraken.rada0j
sed -e 's/rada0 i j/rada0 j k/g' -e 's/rada0j/rada0k/g' batch.qsub.kraken.rada0j > batch.qsub.kraken.rada0k
sed -e 's/rada0 j k/rada0 k l/g' -e 's/rada0k/rada0l/g' batch.qsub.kraken.rada0k > batch.qsub.kraken.rada0l
sed -e 's/rada0 k l/rada0 l m/g' -e 's/rada0l/rada0m/g' batch.qsub.kraken.rada0l > batch.qsub.kraken.rada0m
sed -e 's/rada0 l m/rada0 m n/g' -e 's/rada0m/rada0n/g' batch.qsub.kraken.rada0m > batch.qsub.kraken.rada0n
sed -e 's/rada0 m n/rada0 n o/g' -e 's/rada0n/rada0o/g' batch.qsub.kraken.rada0n > batch.qsub.kraken.rada0o
sed -e 's/rada0 n o/rada0 o p/g' -e 's/rada0o/rada0p/g' batch.qsub.kraken.rada0o > batch.qsub.kraken.rada0p
sed -e 's/rada0 o p/rada0 p q/g' -e 's/rada0p/rada0q/g' batch.qsub.kraken.rada0p > batch.qsub.kraken.rada0q
sed -e 's/rada0 p q/rada0 q r/g' -e 's/rada0q/rada0r/g' batch.qsub.kraken.rada0q > batch.qsub.kraken.rada0r
sed -e 's/rada0 q r/rada0 r s/g' -e 's/rada0r/rada0s/g' batch.qsub.kraken.rada0r > batch.qsub.kraken.rada0s
sed -e 's/rada0 r s/rada0 s t/g' -e 's/rada0s/rada0t/g' batch.qsub.kraken.rada0s > batch.qsub.kraken.rada0t
sed -e 's/rada0 s t/rada0 t u/g' -e 's/rada0t/rada0u/g' batch.qsub.kraken.rada0t > batch.qsub.kraken.rada0u
sed -e 's/rada0 t u/rada0 u v/g' -e 's/rada0u/rada0v/g' batch.qsub.kraken.rada0u > batch.qsub.kraken.rada0v
sed -e 's/rada0 u v/rada0 v w/g' -e 's/rada0v/rada0w/g' batch.qsub.kraken.rada0v > batch.qsub.kraken.rada0w
sed -e 's/rada0 v w/rada0 w x/g' -e 's/rada0w/rada0x/g' batch.qsub.kraken.rada0w > batch.qsub.kraken.rada0x
sed -e 's/rada0 w x/rada0 x y/g' -e 's/rada0x/rada0y/g' batch.qsub.kraken.rada0x > batch.qsub.kraken.rada0y
sed -e 's/rada0 x y/rada0 y z/g' -e 's/rada0y/rada0z/g' batch.qsub.kraken.rada0y > batch.qsub.kraken.rada0z
	    
sed -e 's/rada0 y z/rada0 z nexta/g' -e 's/rada0z/rada0nexta/g' batch.qsub.kraken.rada0z > batch.qsub.kraken.rada0nexta

sed -e 's/rada0 z nexta/rada0 nexta nextb/g' -e 's/rada0nexta/rada0nextb/g' batch.qsub.kraken.rada0nexta > batch.qsub.kraken.rada0nextb

sed -e 's/rada0 nexta nextb/rada0 nextb nextc/g' -e 's/rada0nextb/rada0nextc/g' batch.qsub.kraken.rada0nextb > batch.qsub.kraken.rada0nextc
sed -e 's/rada0 nextb nextc/rada0 nextc nextd/g' -e 's/rada0nextc/rada0nextd/g' batch.qsub.kraken.rada0nextc > batch.qsub.kraken.rada0nextd
sed -e 's/rada0 nextc nextd/rada0 nextd nexte/g' -e 's/rada0nextd/rada0nexte/g' batch.qsub.kraken.rada0nextd > batch.qsub.kraken.rada0nexte
sed -e 's/rada0 nextd nexte/rada0 nexte nextf/g' -e 's/rada0nexte/rada0nextf/g' batch.qsub.kraken.rada0nexte > batch.qsub.kraken.rada0nextf
sed -e 's/rada0 nexte nextf/rada0 nextf nextg/g' -e 's/rada0nextf/rada0nextg/g' batch.qsub.kraken.rada0nextf > batch.qsub.kraken.rada0nextg
sed -e 's/rada0 nextf nextg/rada0 nextg nexth/g' -e 's/rada0nextg/rada0nexth/g' batch.qsub.kraken.rada0nextg > batch.qsub.kraken.rada0nexth
sed -e 's/rada0 nextg nexth/rada0 nexth nexti/g' -e 's/rada0nexth/rada0nexti/g' batch.qsub.kraken.rada0nexth > batch.qsub.kraken.rada0nexti
sed -e 's/rada0 nexth nexti/rada0 nexti nextj/g' -e 's/rada0nexti/rada0nextj/g' batch.qsub.kraken.rada0nexti > batch.qsub.kraken.rada0nextj
sed -e 's/rada0 nexti nextj/rada0 nextj nextk/g' -e 's/rada0nextj/rada0nextk/g' batch.qsub.kraken.rada0nextj > batch.qsub.kraken.rada0nextk
sed -e 's/rada0 nextj nextk/rada0 nextk nextl/g' -e 's/rada0nextk/rada0nextl/g' batch.qsub.kraken.rada0nextk > batch.qsub.kraken.rada0nextl
sed -e 's/rada0 nextk nextl/rada0 nextl nextm/g' -e 's/rada0nextl/rada0nextm/g' batch.qsub.kraken.rada0nextl > batch.qsub.kraken.rada0nextm
sed -e 's/rada0 nextl nextm/rada0 nextm nextn/g' -e 's/rada0nextm/rada0nextn/g' batch.qsub.kraken.rada0nextm > batch.qsub.kraken.rada0nextn
sed -e 's/rada0 nextm nextn/rada0 nextn nexto/g' -e 's/rada0nextn/rada0nexto/g' batch.qsub.kraken.rada0nextn > batch.qsub.kraken.rada0nexto
sed -e 's/rada0 nextn nexto/rada0 nexto nextp/g' -e 's/rada0nexto/rada0nextp/g' batch.qsub.kraken.rada0nexto > batch.qsub.kraken.rada0nextp
sed -e 's/rada0 nexto nextp/rada0 nextp nextq/g' -e 's/rada0nextp/rada0nextq/g' batch.qsub.kraken.rada0nextp > batch.qsub.kraken.rada0nextq
sed -e 's/rada0 nextp nextq/rada0 nextq nextr/g' -e 's/rada0nextq/rada0nextr/g' batch.qsub.kraken.rada0nextq > batch.qsub.kraken.rada0nextr
sed -e 's/rada0 nextq nextr/rada0 nextr nexts/g' -e 's/rada0nextr/rada0nexts/g' batch.qsub.kraken.rada0nextr > batch.qsub.kraken.rada0nexts
sed -e 's/rada0 nextr nexts/rada0 nexts nextt/g' -e 's/rada0nexts/rada0nextt/g' batch.qsub.kraken.rada0nexts > batch.qsub.kraken.rada0nextt
sed -e 's/rada0 nexts nextt/rada0 nextt nextu/g' -e 's/rada0nextt/rada0nextu/g' batch.qsub.kraken.rada0nextt > batch.qsub.kraken.rada0nextu
sed -e 's/rada0 nextt nextu/rada0 nextu nextv/g' -e 's/rada0nextu/rada0nextv/g' batch.qsub.kraken.rada0nextu > batch.qsub.kraken.rada0nextv
sed -e 's/rada0 nextu nextv/rada0 nextv nextw/g' -e 's/rada0nextv/rada0nextw/g' batch.qsub.kraken.rada0nextv > batch.qsub.kraken.rada0nextw
sed -e 's/rada0 nextv nextw/rada0 nextw nextx/g' -e 's/rada0nextw/rada0nextx/g' batch.qsub.kraken.rada0nextw > batch.qsub.kraken.rada0nextx
sed -e 's/rada0 nextw nextx/rada0 nextx nexty/g' -e 's/rada0nextx/rada0nexty/g' batch.qsub.kraken.rada0nextx > batch.qsub.kraken.rada0nexty
sed -e 's/rada0 nextx nexty/rada0 nexty nextz/g' -e 's/rada0nexty/rada0nextz/g' batch.qsub.kraken.rada0nexty > batch.qsub.kraken.rada0nextz
