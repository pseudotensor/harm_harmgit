#!/bin/bash

# setting up for dependency list
# http://beige.ucs.indiana.edu/I590/node45.html
#

sed -e 's/firstrun=1/firstrun=0/g' -e 's/radtma0.8 - a/radtma0.8 a b/g' -e 's/radtma0.8a/radtma0.8b/g' batch.slurm.kraken.radtma0.8a > batch.slurm.kraken.radtma0.8b
                        
sed -e 's/radtma0.8 a b/radtma0.8 b c/g' -e 's/radtma0.8b/radtma0.8c/g' batch.slurm.kraken.radtma0.8b > batch.slurm.kraken.radtma0.8c
sed -e 's/radtma0.8 b c/radtma0.8 c d/g' -e 's/radtma0.8c/radtma0.8d/g' batch.slurm.kraken.radtma0.8c > batch.slurm.kraken.radtma0.8d
sed -e 's/radtma0.8 c d/radtma0.8 d e/g' -e 's/radtma0.8d/radtma0.8e/g' batch.slurm.kraken.radtma0.8d > batch.slurm.kraken.radtma0.8e
sed -e 's/radtma0.8 d e/radtma0.8 e f/g' -e 's/radtma0.8e/radtma0.8f/g' batch.slurm.kraken.radtma0.8e > batch.slurm.kraken.radtma0.8f
sed -e 's/radtma0.8 e f/radtma0.8 f g/g' -e 's/radtma0.8f/radtma0.8g/g' batch.slurm.kraken.radtma0.8f > batch.slurm.kraken.radtma0.8g
sed -e 's/radtma0.8 f g/radtma0.8 g h/g' -e 's/radtma0.8g/radtma0.8h/g' batch.slurm.kraken.radtma0.8g > batch.slurm.kraken.radtma0.8h
sed -e 's/radtma0.8 g h/radtma0.8 h i/g' -e 's/radtma0.8h/radtma0.8i/g' batch.slurm.kraken.radtma0.8h > batch.slurm.kraken.radtma0.8i
sed -e 's/radtma0.8 h i/radtma0.8 i j/g' -e 's/radtma0.8i/radtma0.8j/g' batch.slurm.kraken.radtma0.8i > batch.slurm.kraken.radtma0.8j
sed -e 's/radtma0.8 i j/radtma0.8 j k/g' -e 's/radtma0.8j/radtma0.8k/g' batch.slurm.kraken.radtma0.8j > batch.slurm.kraken.radtma0.8k
sed -e 's/radtma0.8 j k/radtma0.8 k l/g' -e 's/radtma0.8k/radtma0.8l/g' batch.slurm.kraken.radtma0.8k > batch.slurm.kraken.radtma0.8l
sed -e 's/radtma0.8 k l/radtma0.8 l m/g' -e 's/radtma0.8l/radtma0.8m/g' batch.slurm.kraken.radtma0.8l > batch.slurm.kraken.radtma0.8m
sed -e 's/radtma0.8 l m/radtma0.8 m n/g' -e 's/radtma0.8m/radtma0.8n/g' batch.slurm.kraken.radtma0.8m > batch.slurm.kraken.radtma0.8n
sed -e 's/radtma0.8 m n/radtma0.8 n o/g' -e 's/radtma0.8n/radtma0.8o/g' batch.slurm.kraken.radtma0.8n > batch.slurm.kraken.radtma0.8o
sed -e 's/radtma0.8 n o/radtma0.8 o p/g' -e 's/radtma0.8o/radtma0.8p/g' batch.slurm.kraken.radtma0.8o > batch.slurm.kraken.radtma0.8p
sed -e 's/radtma0.8 o p/radtma0.8 p q/g' -e 's/radtma0.8p/radtma0.8q/g' batch.slurm.kraken.radtma0.8p > batch.slurm.kraken.radtma0.8q
sed -e 's/radtma0.8 p q/radtma0.8 q r/g' -e 's/radtma0.8q/radtma0.8r/g' batch.slurm.kraken.radtma0.8q > batch.slurm.kraken.radtma0.8r
sed -e 's/radtma0.8 q r/radtma0.8 r s/g' -e 's/radtma0.8r/radtma0.8s/g' batch.slurm.kraken.radtma0.8r > batch.slurm.kraken.radtma0.8s
sed -e 's/radtma0.8 r s/radtma0.8 s t/g' -e 's/radtma0.8s/radtma0.8t/g' batch.slurm.kraken.radtma0.8s > batch.slurm.kraken.radtma0.8t
sed -e 's/radtma0.8 s t/radtma0.8 t u/g' -e 's/radtma0.8t/radtma0.8u/g' batch.slurm.kraken.radtma0.8t > batch.slurm.kraken.radtma0.8u
sed -e 's/radtma0.8 t u/radtma0.8 u v/g' -e 's/radtma0.8u/radtma0.8v/g' batch.slurm.kraken.radtma0.8u > batch.slurm.kraken.radtma0.8v
sed -e 's/radtma0.8 u v/radtma0.8 v w/g' -e 's/radtma0.8v/radtma0.8w/g' batch.slurm.kraken.radtma0.8v > batch.slurm.kraken.radtma0.8w
sed -e 's/radtma0.8 v w/radtma0.8 w x/g' -e 's/radtma0.8w/radtma0.8x/g' batch.slurm.kraken.radtma0.8w > batch.slurm.kraken.radtma0.8x
sed -e 's/radtma0.8 w x/radtma0.8 x y/g' -e 's/radtma0.8x/radtma0.8y/g' batch.slurm.kraken.radtma0.8x > batch.slurm.kraken.radtma0.8y
sed -e 's/radtma0.8 x y/radtma0.8 y z/g' -e 's/radtma0.8y/radtma0.8z/g' batch.slurm.kraken.radtma0.8y > batch.slurm.kraken.radtma0.8z
	    
sed -e 's/radtma0.8 y z/radtma0.8 z nexta/g' -e 's/radtma0.8z/radtma0.8nexta/g' batch.slurm.kraken.radtma0.8z > batch.slurm.kraken.radtma0.8nexta

sed -e 's/radtma0.8 z nexta/radtma0.8 nexta nextb/g' -e 's/radtma0.8nexta/radtma0.8nextb/g' batch.slurm.kraken.radtma0.8nexta > batch.slurm.kraken.radtma0.8nextb

sed -e 's/radtma0.8 nexta nextb/radtma0.8 nextb nextc/g' -e 's/radtma0.8nextb/radtma0.8nextc/g' batch.slurm.kraken.radtma0.8nextb > batch.slurm.kraken.radtma0.8nextc
sed -e 's/radtma0.8 nextb nextc/radtma0.8 nextc nextd/g' -e 's/radtma0.8nextc/radtma0.8nextd/g' batch.slurm.kraken.radtma0.8nextc > batch.slurm.kraken.radtma0.8nextd
sed -e 's/radtma0.8 nextc nextd/radtma0.8 nextd nexte/g' -e 's/radtma0.8nextd/radtma0.8nexte/g' batch.slurm.kraken.radtma0.8nextd > batch.slurm.kraken.radtma0.8nexte
sed -e 's/radtma0.8 nextd nexte/radtma0.8 nexte nextf/g' -e 's/radtma0.8nexte/radtma0.8nextf/g' batch.slurm.kraken.radtma0.8nexte > batch.slurm.kraken.radtma0.8nextf
sed -e 's/radtma0.8 nexte nextf/radtma0.8 nextf nextg/g' -e 's/radtma0.8nextf/radtma0.8nextg/g' batch.slurm.kraken.radtma0.8nextf > batch.slurm.kraken.radtma0.8nextg
sed -e 's/radtma0.8 nextf nextg/radtma0.8 nextg nexth/g' -e 's/radtma0.8nextg/radtma0.8nexth/g' batch.slurm.kraken.radtma0.8nextg > batch.slurm.kraken.radtma0.8nexth
sed -e 's/radtma0.8 nextg nexth/radtma0.8 nexth nexti/g' -e 's/radtma0.8nexth/radtma0.8nexti/g' batch.slurm.kraken.radtma0.8nexth > batch.slurm.kraken.radtma0.8nexti
sed -e 's/radtma0.8 nexth nexti/radtma0.8 nexti nextj/g' -e 's/radtma0.8nexti/radtma0.8nextj/g' batch.slurm.kraken.radtma0.8nexti > batch.slurm.kraken.radtma0.8nextj
sed -e 's/radtma0.8 nexti nextj/radtma0.8 nextj nextk/g' -e 's/radtma0.8nextj/radtma0.8nextk/g' batch.slurm.kraken.radtma0.8nextj > batch.slurm.kraken.radtma0.8nextk
sed -e 's/radtma0.8 nextj nextk/radtma0.8 nextk nextl/g' -e 's/radtma0.8nextk/radtma0.8nextl/g' batch.slurm.kraken.radtma0.8nextk > batch.slurm.kraken.radtma0.8nextl
sed -e 's/radtma0.8 nextk nextl/radtma0.8 nextl nextm/g' -e 's/radtma0.8nextl/radtma0.8nextm/g' batch.slurm.kraken.radtma0.8nextl > batch.slurm.kraken.radtma0.8nextm
sed -e 's/radtma0.8 nextl nextm/radtma0.8 nextm nextn/g' -e 's/radtma0.8nextm/radtma0.8nextn/g' batch.slurm.kraken.radtma0.8nextm > batch.slurm.kraken.radtma0.8nextn
sed -e 's/radtma0.8 nextm nextn/radtma0.8 nextn nexto/g' -e 's/radtma0.8nextn/radtma0.8nexto/g' batch.slurm.kraken.radtma0.8nextn > batch.slurm.kraken.radtma0.8nexto
sed -e 's/radtma0.8 nextn nexto/radtma0.8 nexto nextp/g' -e 's/radtma0.8nexto/radtma0.8nextp/g' batch.slurm.kraken.radtma0.8nexto > batch.slurm.kraken.radtma0.8nextp
sed -e 's/radtma0.8 nexto nextp/radtma0.8 nextp nextq/g' -e 's/radtma0.8nextp/radtma0.8nextq/g' batch.slurm.kraken.radtma0.8nextp > batch.slurm.kraken.radtma0.8nextq
sed -e 's/radtma0.8 nextp nextq/radtma0.8 nextq nextr/g' -e 's/radtma0.8nextq/radtma0.8nextr/g' batch.slurm.kraken.radtma0.8nextq > batch.slurm.kraken.radtma0.8nextr
sed -e 's/radtma0.8 nextq nextr/radtma0.8 nextr nexts/g' -e 's/radtma0.8nextr/radtma0.8nexts/g' batch.slurm.kraken.radtma0.8nextr > batch.slurm.kraken.radtma0.8nexts
sed -e 's/radtma0.8 nextr nexts/radtma0.8 nexts nextt/g' -e 's/radtma0.8nexts/radtma0.8nextt/g' batch.slurm.kraken.radtma0.8nexts > batch.slurm.kraken.radtma0.8nextt
sed -e 's/radtma0.8 nexts nextt/radtma0.8 nextt nextu/g' -e 's/radtma0.8nextt/radtma0.8nextu/g' batch.slurm.kraken.radtma0.8nextt > batch.slurm.kraken.radtma0.8nextu
sed -e 's/radtma0.8 nextt nextu/radtma0.8 nextu nextv/g' -e 's/radtma0.8nextu/radtma0.8nextv/g' batch.slurm.kraken.radtma0.8nextu > batch.slurm.kraken.radtma0.8nextv
sed -e 's/radtma0.8 nextu nextv/radtma0.8 nextv nextw/g' -e 's/radtma0.8nextv/radtma0.8nextw/g' batch.slurm.kraken.radtma0.8nextv > batch.slurm.kraken.radtma0.8nextw
sed -e 's/radtma0.8 nextv nextw/radtma0.8 nextw nextx/g' -e 's/radtma0.8nextw/radtma0.8nextx/g' batch.slurm.kraken.radtma0.8nextw > batch.slurm.kraken.radtma0.8nextx
sed -e 's/radtma0.8 nextw nextx/radtma0.8 nextx nexty/g' -e 's/radtma0.8nextx/radtma0.8nexty/g' batch.slurm.kraken.radtma0.8nextx > batch.slurm.kraken.radtma0.8nexty
sed -e 's/radtma0.8 nextx nexty/radtma0.8 nexty nextz/g' -e 's/radtma0.8nexty/radtma0.8nextz/g' batch.slurm.kraken.radtma0.8nexty > batch.slurm.kraken.radtma0.8nextz
