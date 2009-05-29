#BEGIN:


# make directory, clean it, and go to it
mkdir run
rm -rf run/*
cd run

# copy over CODE and STELLAR MODELS
cp ../grmhd ~/research/grbmodel/stellar*.txt .


# latest monotonized w/ degen offset
export DIR=~/research/helm/100x50x20x50/
cp $DIR/eosnew.dat $DIR/eosnew.head $DIR/eosdegennew.dat .


# most recent HELM version is below (monotonized version with degen offset)
export DIR=~/research/helm/200x200x50x1/
cp $DIR/eosnew.dat eossimplenew.dat
cp $DIR/eosnew.head eossimplenew.head
cp $DIR/eosdegennew.dat eosdegensimplenew.dat

# SIMPLE ZOOM TABLE (not used now)


# RUN
nohup ./grmhd &


# DONE
