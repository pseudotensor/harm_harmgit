#nohup x11vnc -rfbauth ~/.vnc/passwd -ncache 10 -forever -xinerama -bg -o logfile.x11vnc.txt &

# -ncache 10 # causes extra screens below for caching, and suggested SSVNC doesn't work even if -ycrop would work
#nohup x11vnc -rfbauth ~/.vnc/passwd  -forever -xinerama -bg -clip 5760x1200+0+0 -nodragging -solid  -progressive 100 -o logfile.x11vnc.txt &

dofull=0
dosplit=0
dofullnew=1
dosplitnew=0
tryx0vnc=0


# realvnc: still only serves up :0


# options to vnc module:
#http://www.realvnc.com/products/free/4.1/man/Xvnc.html


# doesn't change to other windows so can't type into windows
if [ $tryx0vnc -eq 1 ]
then
    nohup x0vncserver -rfbport 5905 -PasswordFile=/home/jon/.vnc/passwd &
fi


if [ $dofull -eq 1 ]
then

    nohup x11vnc -rfbauth ~/.vnc/passwd  -forever -xinerama -bg -clip 5760x1200+0+0 -nodragging -solid  -progressive 100 -rfbport 5904 -o logfile.x11vnc.txt &
fi


if [ $dosplit -eq 1 ]
then

############## 3 separate windows since full 3-screen window is fast on central window and super death slow on left-right windows

    nohup x11vnc -rfbauth ~/.vnc/passwd  -forever -xinerama -bg -clip xinerama0 -nodragging -solid  -progressive 100 -rfbport 5901 -o logfile0.x11vnc.txt &
    
    nohup x11vnc -rfbauth ~/.vnc/passwd  -forever -xinerama -bg -clip xinerama1 -nodragging -solid  -progressive 100 -rfbport 5902 -o logfile1.x11vnc.txt &
    
    nohup x11vnc -rfbauth ~/.vnc/passwd  -forever -xinerama -bg -clip xinerama2 -nodragging -solid  -progressive 100 -rfbport 5903 -o logfile2.x11vnc.txt &
    
fi

##################



if [ $dofullnew -eq 1 ]
then

#    nohup x11vnc -rfbauth ~/.vnc/passwd  -forever -bg -rfbport 5904 -deferupdate none -xtrap -wait none -defer none -o logfilenew.x11vnc.txt &

    # worked for a while, but then after some point crashes occurred and couldn't log-in and then got slow again on left-right screens even after restarting
#    nohup x11vnc -deferupdate none -wait none -defer none -rfbauth ~/.vnc/passwd  -forever -bg -rfbport 5904 -o logfilenew.x11vnc.txt &

# reports slow netrate: netrate: 102 KB/sec, latency: 10 ms
# reports: increased wireframe timeouts for slow network connection

# fast again with below.
    nohup x11vnc -noxdamage -deferupdate none -wait none -defer none -rfbauth ~/.vnc/passwd  -forever -bg -rfbport 5904 -o logfilenew.x11vnc.txt &
fi


if [ $dosplitnew -eq 1 ]
then

############## 3 separate windows since full 3-screen window is fast on central window and super death slow on left-right windows

#    nohup x11vnc -rfbauth ~/.vnc/passwd  -forever  -bg -clip xinerama0  -rfbport 5901 -o logfile0.x11vnc.txt &
    
    nohup x11vnc -rfbauth ~/.vnc/passwd  -forever  -bg -clip xinerama1  -rfbport 5902 -o logfile1.x11vnc.txt &
    
#    nohup x11vnc -rfbauth ~/.vnc/passwd  -forever  -bg -clip xinerama2  -rfbport 5903 -o logfile2.x11vnc.txt &
    
fi



#####


##
