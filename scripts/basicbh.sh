sudo bash

hostname bh00.physics.umd.edu

cd /bin
rm -rf sh
ln -s bash sh

cd /mnt/data0/
mkdir user
chown -R user user
chgrp -R user user

cd /home
ln -s user jon
chown -R user user
chgrp -R user user

# run as 1 line
sudo apt-get -y install mpich2 python python2.7-scipy python-dateutil yasm ipython python-setuptools python-nose python2.7-scipy python2.7-numpy python-matplotlib python-matplotlib-data dvipng emacs make f2c

sudo aptitude install libhdf4-0 libhdf4-dev libhdf4-doc hdf4-tools

# run separately
cd /mnt/data0/
mkdir opt
scp user@bh01:/mnt/data0/opt.tar .

# run separately
tar xvf opt.tar
chmod a+rx opt/
cd /
rm -rf opt
ln -s /mnt/data0/opt/ .
cd /home/
ln -s /opt/ .


exit

cd ~/
scp user@bh02:.bash* .
scp user@bh02:.profile* .
scp user@bh02:.emacs .
scp user@bh02:.sm .
scp user@bh02:smstuff.tgz .
scp user@bh02:binstuff.tgz .
scp user@bh02:intel.tgz .
scp user@bh02:mpich-install.tgz .
tar xvzf smstuff.tgz 
ln -s harm_harmsmmacros sm
tar xvzf binstuff.tgz
tar xvzf intel.tgz 
tar xvzf mpich-install.tgz 



# on other computer:
tar cvzf waldnew4.tgz *.c *.P *.inc make* initboundcode *.f *.h
scp  waldnew4.tgz user@bh05:/mnt/data0/user/

exit and relogin to activate bash stuff.
