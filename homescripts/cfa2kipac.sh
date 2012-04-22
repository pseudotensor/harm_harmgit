# Done on local computer
#
#
#
#
# Add ~/bin to PATH in .profile
# 
########### BEGIN CUT
ulimit -S -c 100000000 > /dev/null 2>&1
#umask 022
umask 077
PATH=~/bin:$PATH
export PATH
########### END CUT
cd /usr/src
#scp jon@relativity.cfa.harvard.edu:research/utils/bin2txt.c jon@relativity.cfa.harvard.edu:research/utils/sm jon@relativity.cfa.harvard.edu:research/utils/vis5d_anatoly jon@relativity.cfa.harvard.edu:research/utils/r8toras jon@relativity.cfa.harvard.edu:research/utils/mympeg2.tar.gz jon@relativity.cfa.harvard.edu:research/utils/ppm2fli jon@relativity.cfa.harvard.edu:research/utils/qslim-2.1 jon@relativity.cfa.harvard.edu:research/utils/gtkglarea-1.2.3/ .
scp -rp jon@relativity.cfa.harvard.edu:"research/utils/bin2txt.c research/utils/sm.tgz research/utils/vis5d_anatoly.tgz research/utils/r8toras research/utils/mympeg2.tar.gz research/utils/ppm2fli research/utils/qslim-2.1.tar.gz research/utils/gtkglarea*.gz" .
#
#
tar xvzf sm.tgz
cd sm/sm2_4_1
#
# Edit Makefile so that X11 libs point to lib64
emacs Makefile    
make clean ; make
cp src/makeyyl /u1/ki/usrlocal/      
cp src/sm /u1/ki/usrlocal/    
#
# Edit ~/.sm to move /usr/local to /u1/ki/usr/local
#
cd /u1/ki/usr/local/
mkdir lib
cd lib
scp jon@relativity.cfa.harvard.edu:/usr/local/lib/smlib.tgz.
tar xvzf smlib.tgz
#
# link to old name so don't have to change references to that directory
cd /home/
ln -s ~/ jon
#
# after this SM works
#
# make easier to get to my usr/local
ln -s /u1/ki/usr /usr2
#
cd ~/bin/   
mv * /usr2/local/bin/
#
#
# change KDE settings so focus strictly follows mouse
#
#
# install nvidia 3D graphics library/driver
#
#


