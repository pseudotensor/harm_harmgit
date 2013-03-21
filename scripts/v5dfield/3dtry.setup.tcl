
set generalaud [expr 0]

set ixmin [expr -40]
set ixmax [expr 40]
set iymin [expr -40]
set iymax [expr 40]
set izmin [expr -40]
set izmax [expr 40]
#set bhspin [expr 0.9375]
set bhspin [expr 0.99]
set startsimtime [expr 0]
set simtime [expr 0]


set W  [expr 8]
set V  [expr 9]
set U  [expr 10]
set W2  [expr 11]
set V2  [expr 12]
set U2  [expr 13]
vis5d_set_wind_vars_and_owners $dtx $ctx $U $ctx $V $ctx $W $ctx $U2 $ctx $V2 $ctx $W2 $ctx $U2 $ctx $V2 $ctx $W2
set posr [expr 14]
set posh [expr 15]
set posph [expr 16]
# created
set rergo [expr 17]
set zeroergo [expr 18]



vis5d_graphics_mode $dtx VIS5D_CLOCK VIS5D_OFF
vis5d_graphics_mode $dtx VIS5D_PRETTY VIS5D_ON
# below actually turns off contour numbers
vis5d_graphics_mode $dtx VIS5D_CONTOUR_NUMBERS VIS5D_ON

vis5d_enable_graphics $ctx VIS5D_TRAJ 0 VIS5D_ON

# small-scale
set outputfilename "small3d"


# Add your own settings file from vis5d's "SAVE" button
#source 3dtry.set
#source 3dtry2.set
#source try11.SAVE
# put posr approx BH horizon directly in try12.SAVE
#source try12.SAVE
source try12nolines.SAVE


#set outputfilename "large3d"
#source largejet.set

# below bh.set must come after settingslatest3.set
puts "before bh.set.tcl"
#source bh.set.tcl
