#Vis5D 4.3 Tcl save file


source 3dtry.setup.tcl




#vis5d_draw_frame 0
#vis5d_draw_frame $ctx
#vis5d_set_traj $ctx 1.0 1.0 $time
#vis5d_make_traj $ctx $newrow $newcol $newlev $time 0
#vis5d_draw_frame $ctx



set format VIS5D_PPM
#set format VIS5D_PNG
#vis5d_save_window 3dtry.ppm $format
#puts "here3"
#vis5d_save_right_window 3dtry.right.ppm $format
#puts "here4"


#set dostereo [expr 2]
set dostereo [expr 0]
#set numtimes [ vis5d_get_dtx_numtimes $ctx ]
set numsteps [expr 1000]

#set skipupto [expr 917]
set skipupto [expr -1]

set skipafter [expr 0]

set myview [ vis5d_get_view $dtx ]
puts "myview $myview"
set view0 [lindex $myview 0]
puts "view0 $view0"
set view1 [lindex $myview 1]
puts "view1 $view1"
set view2 [lindex $myview 2]
puts "view2 $view2"
set view3 [lindex $myview 3]
puts "view3 $view3"
set view4 [lindex $myview 4]
puts "view4 $view4"
set view5 [lindex $myview 5]
puts "view5 $view5"
set view6 [lindex $myview 6]
puts "view6 $view6"
#


for {set step 0} {$step<$numsteps} {set step [expr $step+1]} {

    puts "step $step"





    ###################################3
    #Viewing matrix
    #vis5d_set_matrix $dtx { 0.60289 0.0781902 -0.793984 0 -0.0533052 0.99691 0.0576983 0 0.796042 0.00753765 0.605196 0 0 0 0 1 }
    #Camera
    #vis5d_set_camera $dtx 0 0 1

    # add to rotation around x-axis
    #set view1 [expr $view1 + $step*360.0/$numsteps]
    # fixed degree change
    set view1new [expr $view1 + $step*360.0/$numsteps]
    puts "view1new $view1new"
    #
    # set new view
    vis5d_set_view $dtx $view0 $view1new $view2 $view3 $view4 $view5 $view6



    if {$step>=$skipupto && $step<=$skipafter} {

	if {$dostereo==2} {

	    vis5d_stereo_on $ctx VIS5D_STEREO_LEFT
	    # redraw the window image for the current dtx timestep
	    puts "before vis5d_draw_frame1"
	    vis5d_draw_frame $ctx
	    vis5d_save_window $outputfilename.left.$step.$dumpnum.ppm $format

	    vis5d_sleep 4000

	    vis5d_stereo_on $ctx VIS5D_STEREO_RIGHT
	    # redraw the window image for the current dtx timestep
	    puts "before vis5d_draw_frame2"
	    vis5d_draw_frame $ctx
	    vis5d_save_window $outputfilename.right.$step.$dumpnum.ppm $format
	}
	if {$dostereo==0} {


	    #
	    # redraw the window image for the current dtx timestep
	    puts "before vis5d_draw_frame"
	    vis5d_draw_frame $ctx
	    #
	    puts "save_window with file: $outputfilename.$step.$dumpnum.ppm"
	    vis5d_save_window $outputfilename.$step.$dumpnum.ppm $format
	}
    }
}


# sleep at end to ensure X processes completed
vis5d_sleep 1000
puts "end after sleep"


# ~/bin/vis5d idumps/fieldline5736.cart.bin.box256x256x256.out17.modelsashaa99t1.5708.v5d -mbs 2802 -geometry 1600x1600 -verylarge 0 -offscreen -script 3dtry.tcl

# ffmpeg -y -fflags +genpts -r 30  -i 3dtry.%d.ppm -sameq -qmax 5 3dtry.mov

# ffmpeg -y -fflags +genpts -r 30  -i small3d.%d.ppm -sameq -qmax 5 small3d.mov

#
# ffmpeg -y -fflags +genpts -r 30  -i large3d.%d.ppm -sameq -qmax 5 large3d.mov
