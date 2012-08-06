# Example .cshrc
# $Header: /cvsroot/bcfg2/bcfg2/Cfg/etc/skel_NAS/.cshrc/.cshrc,v 1.1 2009/12/11 16:05:12 dtalcott Exp $

# Include the NAS standard .cshrc

source /usr/local/lib/init/global.cshrc

# Stop here unless interactive shell

if (! $?prompt ) exit

# Nothing above this line should produce any output.  If it does, it
# will interfere with scp, sftp, etc.

set noclobber				# Don't clobber file with redirect

set cdpath = (. $HOME)
set history=50                          # C shell history
set savehist=100                        # history saved between logins
set ignoreeof                           # prevent accidental exit
set filec				# file completion

set prompt="`uname -n`.$USER \!> "      # prompt

# Standard aliases
alias rm /bin/rm -i                     # prompting remove
alias mv /bin/mv -i                     # prompting move
alias cp /bin/cp -i                     # prompting copy

# This is is needed to set terminal type.
alias ts 'set noglob; eval `tset -sQ \!*`; unset noglob; set term=$TERM'

# PFE modules
module load comp-intel/2012.0.032
module load mpi-sgi/mpt.2.04.10789

# Add your own aliases here.
