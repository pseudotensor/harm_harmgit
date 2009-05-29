#BEGIN:


# make directory, clean it, and go to it
mkdir run
rm -rf run/*
cd run

# copy over CODE and STELLAR MODELS
cp ../grmhd ~/research/grbmodel/stellar*.txt .

# LARGE TABLE
#cp ~/research/kazeos/eoslarge_corrxnuc_tkcolumn/eosnew.dat ~/research/kazeos/eoslarge_corrxnuc_tkcolumn/eosnew.head .
#cp ~/research/kazeos/run.4848163/"eosnew.dat eosnew.head eosnewdegen.dat" .

# latest monotonized w/ degen offset
export DIR=~/research/kazeos/4848163.allfixed.new/
cp $DIR/eosnew.dat $DIR/eosnew.head $DIR/eosdegennew.dat .


# SIMPLE TABLE
#cp ~/research/kazeos/eoslargesimple/eossimplenew.dat ~/research/kazeos/eoslargesimple/eossimplenew.head .
#cp ~/research/kazeos/eoslargesimple.200sq.1e40tdyn.1em2hcm/eossimplenew.dat ~/research/kazeos/eoslargesimple.200sq.1e40tdyn.1em2hcm/eossimplenew.head .
#cp ~/research/kazeos/run.200sq.1em15tdyn.1em2hcm/eossimplenew.dat ~/research/kazeos/run.200sq.1em15tdyn.1em2hcm/eossimplenew.head .
#cp ~/research/kazeos/run.200sq.1e15tdyn.1em15hcm/eossimplenew.dat ~/research/kazeos/run.200sq.1e15tdyn.1em15hcm/eossimplenew.head .
#
# most recent kaz version is below (monotonized version with degen offset)
#export DIR=~/research/kazeos/allfixed_200sq_new/
#cp $DIR/eossimplenew.dat $DIR/eossimplenew.head $DIR/eosdegensimplenew.dat .


# most recent HELM version is below (monotonized version with degen offset)
export DIR=~/research/helm/helm_mutot_linearmutot_yefit_200sq/
cp $DIR/eossimplenew.dat $DIR/eossimplenew.head $DIR/eosdegensimplenew.dat .

#
#
# test without offset
#export DIR=~/research/kazeos/testnooffset/
#cp $DIR/eossimplenew.dat $DIR/eossimplenew.head $DIR/eosdegensimplenew.dat .
#
# most recent HELM (OUTPUTTYPE==0 with entropy floor only) version (includes mutot/yefit)
#cp ~/research/helm/helm_mutot_yefit/eossimplenew.dat ~/research/helm/helm_mutot_yefit/eossimplenew.head .

# OLDER VERSIONS
#cp ~/research/kazeos/pwf99/eossimplenew.dat ~/research/kazeos/pwf99/eossimplenew.head .
#cp ~/research/helm/eoshelm/eossimplenew.dat ~/research/helm/eoshelm/eossimplenew.head .
#cp ~/research/helm/purehelmnew/eossimplenew.dat ~/research/helm/purehelmnew/eossimplenew.head .
#cp ~/research/helm/helmkaz/eossimplenew.dat ~/research/helm/helmkaz/eossimplenew.head .


# SIMPLE ZOOM TABLE
# below used bad original temperature range
#cp ~/research/helm/helm_mutot_yefit_refined50x800/eossimplezoomnew.dat ~/research/helm/helm_mutot_yefit_refined50x800/eossimplezoomnew.head .
# below uses extended temperature range, but limited u range
cp ~/research/helm/helm_mutot_yefit_refined100x800/eossimplezoomnew.dat ~/research/helm/helm_mutot_yefit_refined100x800/eossimplezoomnew.head .


# RUN
nohup ./grmhd &


# DONE
