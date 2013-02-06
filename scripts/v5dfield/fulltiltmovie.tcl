# movie.tcl


# set a unquotedstring
# set b "substitution processed string"
# set c {quoted literal string}
# set d [results of a command]


# A Vis5d Tcl script for saving a whole time series of images




#########################
#
# Source some default visuals
#
#########################

#source settings.set
#source settings2.set
#source settings3.set
#source settings4.set
#source settings6.set

#source settings7.set
source fulltilt8things.SAVE
#source settings8.set



#########################
#
# Set variable associations
#
#########################

# returns 18 things
#set getit [vis5d_get_wind_vars_and_owners $dtx]
#puts $dtx $getit
# 0 9 0 7 0 8 0 9 0 7 0 8 0 9 0 7 0 8
puts "here0"
#set U  [expr 10]
#set V  [expr 8]
#set W  [expr 9]
#set U2  [expr 13]
#set V2  [expr 11]
#set W2  [expr 12]

set W  [expr 8]
set V  [expr 9]
set U  [expr 10]
set W2  [expr 11]
set V2  [expr 12]
set U2  [expr 13]
vis5d_set_wind_vars_and_owners $dtx $ctx $U $ctx $V $ctx $W $ctx $U2 $ctx $V2 $ctx $W2 $ctx $U2 $ctx $V2 $ctx $W2
#vis5d_set_wind_vars_and_owners $dtx $getit

vis5d_graphics_mode $dtx VIS5D_CLOCK VIS5D_OFF
vis5d_graphics_mode $dtx VIS5D_PRETTY VIS5D_ON
# below actually turns off contour numbers
vis5d_graphics_mode $dtx VIS5D_CONTOUR_NUMBERS VIS5D_ON

vis5d_enable_graphics $ctx VIS5D_TRAJ 0 VIS5D_ON

#puts "here0.5"


#########################
#
# Get x,y,z position of inner and outer edges of box
#
#########################

# assume grid is constant in time
set time [expr 0]
# just pick first var since all vars have same size
set var [expr 0]

set rowmin [expr 0]
set colmin [expr 0]
set levmin [expr 0]
set triple "$rowmin $colmin $levmin"
set xyz [vis5d_grid_to_xyz $dtx $time $var $triple]
foreach {xmin ymin zmin} "$xyz" break
# -1.000000 1.000000 -0.968750

# max row,col,lev
set numrow [vis5d_get_dtx_grid_rows $dtx]
set numcol [vis5d_get_dtx_grid_columns $dtx]
set numlev [vis5d_get_dtx_grid_levels $dtx]

set rowmax [expr $numrow-1]
set colmax [expr $numcol-1]
set levmax [expr $numlev-1]
set triple "$rowmax $colmax $levmax"
set xyz [vis5d_grid_to_xyz $dtx $time $var $triple]
foreach {xmax ymax zmax} "$xyz" break
# 1.000000 -1.000000 0.968750

#puts "here0.75"

##############
#
# Set translation for user and stupid v5d grid
#
# Grid actually only goes from -1..1 in all directions due to make_box() in src/box.c
#
# Careful with floating point vs. integer.  Need to cast integer or inputted numbers into floats
#
# vis5d_xyz_to_grid data_context time var {x y z}
#    Convert a data graphics (x,y,z) coordinate to a data grid coordinate. The resulting grid coordinate is returned as a list of three values: {row column level}.
#
# NOTE!  x~col y~row z~lev according to proj.c
#
#
# vis5d_set_traj display_context step length ribbonflag
#    Set trajectory attributes.
#        * step - trajectory step size. 1.0 is default
#        * length - trajectory length scaling. 1.0 is default
#        * ribbonflag - if non-zero draw a ribbon trajectory, else draw as line segments. 
#
# vis5d_make_traj display_context row column level dtx_time trajset
#    Compute a window trajectory.
#        * row, column, level - starting coordinate of trajectory in dtx (virtual) grid coordinates
#        * dtx_time - display context timestep at which to start trajectory trace
#        * trajset - which trajectory set 
#
#
# Note!: bin2txt setup so that (x-z plane for \phi=y=0) so: row=yhat col=zhat lev=xhat
#        Normally x~col y~row z~lev according to proj.c
#        So jonxhat = v5dz
#           jonyhat = v5dy
#           jonzhat = v5dx
#
#
###############

set pi 3.14159265
#puts $pi
set fldph [expr 16.0]
set finalph [expr 2.0*$pi ]
set dph [expr $finalph/$fldph]
#set startz [expr -3.0]
#set finalz [expr 3.0]
set startz [expr -1.5]
set finalz [expr 1.5]
set numz [expr 2.0]
#set numz [expr 1.0]
set dz [expr $finalz/$numz]
#puts "here0.8"
puts $dph
for {set z $startz} {$z < $finalz} {set z [expr $z+$dz]} {
for {set ph 0.1} {$ph < $finalph} {set ph [expr $ph+$dph]} {

    #puts "here0.85"
    # cylindrical surface
    #set R [expr 10]
    set R [expr 1.2]
    #set z [expr 3]


    set userx [expr $R*cos($ph)]
    set usery [expr $R*sin($ph)]
    set userz [expr $z]

    #puts "here0.86"

    #puts $userx
    #puts $usery
    #puts $userz

    #set userx [expr -10.0]
    #set usery [expr 0.0]
    #set userz [expr 0.0]
    #
    set v5dorderx [expr $userz]
    set v5dordery [expr $usery]
    set v5dorderz [expr $userx]
    #
    #puts "here0.87"
    # translate
    set v5dx [expr ($v5dorderx-$ixmin)/(1.0*$ixmax-$ixmin)*(1.0*$xmax-$xmin) + $xmin]
    set v5dy [expr ($v5dordery-$iymin)/(1.0*$iymax-$iymin)*(1.0*$ymin-$ymax) + $ymax]
    set v5dz [expr ($v5dorderz-$izmin)/(1.0*$izmax-$izmin)*(1.0*$zmax-$zmin) + $zmin]
    set triple "$v5dx $v5dy $v5dz"
    # get internal position
    #puts "here0.88"
    set pos [vis5d_xyz_to_grid data_context $time $var $triple]
    #puts $pos
    foreach {v5drow v5dcol v5dlev} "$pos" break
    # 15.5 11.625 15.5
    #puts "here0.89"    
    #Trajectories
    vis5d_set_traj $dtx 1.0 1.0 0
    #8.922789 16.736794 13.349346
    #puts $v5drow
    #puts $v5dcol
    #puts $v5dlev
    #puts "here0.895"
    vis5d_make_traj $dtx $v5drow $v5dcol $v5dlev $time 0
    #vis5d_make_traj $dtx 8.922789 16.73679 15.5 $time 0
    #vis5d_make_traj $dtx 16.73679 8.922789 15.5 $time 0
    
    # row col lev
    # zhat
    #9.562092 15.500107 15.275778
}
}

#puts "here0.9"

set format VIS5D_PPM
set numtimes [ vis5d_get_dtx_numtimes $ctx ]

puts "$outputfilename"

for {set time 0} {$time<$numtimes} {set time [expr $time+1]} {

    puts "time $time"

    set newrow [expr 0]
    set newcol [expr 0]
    set newlev [expr 0]


    # V5D starts counting at $time=0
    puts "here1"
    # set which time to set as current context
    vis5d_set_dtx_timestep $ctx $time
    puts "here2"
    # redraw the window image for the current dtx timestep
    vis5d_draw_frame $ctx

    #vis5d_set_traj $ctx 1.0 1.0 $time
    #vis5d_make_traj $ctx $newrow $newcol $newlev $time 0

    #vis5d_draw_frame $ctx

    puts "$outputfilename $format"
    # Save the large window image to the named file. format is one of the file format identifiers returned by vis5d_get_image_formats.
    vis5d_save_window $outputfilename $format
    puts "here3"

}


# At this point you may want to invoke a Unix utility to convert a
# series of GIF files into an MPEG animation file.....
