#BEGIN:


# make directory, clean it, and go to it
mkdir run
rm -rf run/*
cd run

# copy over CODE and STELLAR MODELS
cp ../grmhd .
# just link stellar models
ln -s ~/research/grbmodel/stellar*.txt .


# latest monotonized w/ degen offset
#export DIR=~/research/helm/100x50x20x50.badyeprec/
#export DIR=~/research/helm/100x50x20x50/
export DIR=~/research/helm/100rhox50tkx50yex20ynux1h/
# just link
ln -s $DIR/eosnew.dat $DIR/eosnew.head $DIR/eosdegennew.dat .


# most recent HELM version is below (monotonized version with degen offset)
#export DIR=~/research/helm/200x200x1x50.badyeprec/
#export DIR=~/research/helm/200x200x1x50/
#
#
#  TEMP::
export DIR=~/research/helm/200x200x50x1x1.ynumethod2/
#
# just link
ln -s $DIR/eosnew.dat eossimplenew.dat
ln -s $DIR/eosnew.head eossimplenew.head
ln -s $DIR/eosdegennew.dat eosdegensimplenew.dat

# SIMPLE ZOOM TABLE (not used now)


# RUN
nohup ./grmhd &


# DONE
