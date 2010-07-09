# .bashrc

# User specific aliases and functions

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi


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




ulimit -S -c 100000000 > /dev/null 2>&1
#umask 022
umask 077

