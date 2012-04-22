# .bashrc

# User specific aliases and functions

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

#alias be='emacs -geometry 150x100'
#alias emacs='emacs -geometry 80x50'
#alias mc='mc -a'
#alias mydiff='diff -bBdp'
#alias myps='ps -auxwf'
#alias myps2='ps -Heo pid,user,%cpu,%mem,bsdstart,cputime,etime,stime,nice,ni,pri,opri,stat,args,cmd'
#alias rm0='find -size 0 | xargs rm'
#alias myindent='indent -kr -lc72 -l72 -fca -fc1 -i2 -lp'
#alias gogrmhd='cd ~/research/grmhdcodes/grmhd3dweno/'
#alias lt='ls -alrt'
#alias lsdir='ls -la | egrep "^d"'
#alias lssdir='ls -ap | grep '/' | sed "s/\///"'
#alias lssdir2='ls -ap| grep '/' | tail -n +3 | sed "s/\///"'
#alias lsh='ls -Flagt $@ | head'
#export dolr='$'
#alias rubendir='ls -al | egrep "^d" | awk "{print "$dolr"9}"'
#alias rmindir='for fil in `lssdir2`; do cd $fil; rm -rf $@; cd ..; done'
#alias dudirs='for fil in `lssdir2`; do du -s $fil; done'
#alias dud='dudirs | sort -n'
#alias listruns='for fil in `find | grep mympirun`; do ls -alrt $fil; done'
#alias xdvi='xdvi -s 3'
#alias grmhdstatus='for fil in `find | grep grmhdoutput.txt`; do echo $fil; tail -5 $fil; done'
#alias xanim='xanim -WrT2'
#alias wgetfull='wget -b -m -k -o wget.log -e robots=off'

alias acroreadnew='acroread -openInNewInstance -openInNewWindow'
alias be='emacs -geometry 150x100'
alias cd..='cd ..'
alias cp='cp -i'
alias d='ls'
alias df='df -h -x supermount'
alias du='du -h'
alias dud='dudirs | sort -n'
alias dudirs='for fil in `lssdir2`; do du -s $fil; done'
alias emacs='emacs -geometry 80x50'
alias gogrmhd='cd ~/research/grmhdcodes/grmhd3dweno/'
alias grmhdstatus='for fil in `find | grep grmhdoutput.txt`; do echo $fil; tail -5 $fil; done'
alias kde='xinit /usr/bin/startkde'
alias l='ls'
alias la='ls -a'
alias listruns='for fil in `find | grep mympirun`; do ls -alrt $fil; done'
alias ll='ls -l'
alias ls='ls -F --color=auto'
alias lsd='ls -d */'
alias lsdir='ls -la | egrep "^d"'
alias lsh='ls -Flagt $@ | head'
alias lssdir='ls -ap | grep / | sed "s/\///"'
alias lssdir2='ls -ap| grep / | tail -n +3 | sed "s/\///"'
alias lt='ls -alrt'
alias mc='mc -a'
alias md='mkdir'
alias mv='mv -i'
alias mydiff='diff -bBdp'
alias mydiff2='diff -bBdpy -W 200'
alias myindent='indent -kr -lc72 -l72 -fca -fc1 -i2 -lp'
alias myps='ps -auxwf'
alias myps2='ps -Heo pid,user,%cpu,%mem,bsdstart,cputime,etime,stime,nice,ni,pri,opri,stat,args,cmd'
alias p='cd -'
alias rd='rmdir'
alias rm='rm -i'
alias rm0='find -size 0 | xargs rm'
alias rmindir='for fil in `lssdir2`; do cd $fil; rm -rf $@; cd ..; done'
alias rubendir='ls -al | egrep "^d" | awk "{print "$dolr"9}"'
alias s='cd ..'
alias wgetfull='wget -b -m -k -o wget.log -e robots=off'
alias xanim='xanim -WrT2'
alias xdvi='xdvi -s 3'
#alias findbyuser='find / -printf \"%u  %s\n\" | awk \''{user[$1]+=$2}; END{ for( i in user) print i \" \" user[i]}\\''
#find / -printf "%u  %s\n" | awk '{user[$1]+=$2}; END{ for( i in user) print i " " user[i]}'
alias rmallsvndots='rm -rf `find . -name .svn`'
alias bbcpjon='bbcp -k -a -P 5'
# below is correct even if emacs doesn't make it look correct visually
alias dujon='myfil=0 ; for fil in `ls -al | awk '\''{print $5}'\''` ; do myfil=`expr ${myfil} + ${fil}` ; done ; echo $myfil'
alias scpresume='rsync --partial --progress --rsh=ssh'
# for below, fill in SRC and DEST as $1 and $2
#-P is short for --partial and --progresss --partial-dir=.partial will put he partially downloaded file into directory called ".partial" for future downloading. It also tells rsync to look in this directory for a partially downloaded file to resume from. --exclude=.* tells rsync to stop any files for directories that start with a dot (eg .partial) 
# http://panela.blog-city.com/resume_scp_after_interrupted_downloads_use_rsync.htm
alias scpresume2='while true; do rsync -P --partial-dir=.partial --exclude=.* $1 $2 ; sleep 10; done'

export SVN_EDITOR=emacs

ulimit -S -c 100000000 > /dev/null 2>&1
#umask 022
umask 077

PATH=~/bin:/usr2/local/bin:/usr/local/bin/:$PATH
export PATH


source /opt/intel/mkl/10.0.3.020/tools/environment/mklvarsem64t.sh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/lib/

#echo "done bashrc"


#if [ "$PS1" ] ; then  
#           mkdir -m 0700 /sys/fs/cgroup/cpu/user/$$
#           echo $$ > /sys/fs/cgroup/cpu/user/$$/tasks
#fi
