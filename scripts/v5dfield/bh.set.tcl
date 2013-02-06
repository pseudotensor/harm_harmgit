puts "bh.set.tcl start"
########### setup black hole
set rhor [expr 1.0+sqrt(1.0-$bhspin*$bhspin)]
set rhorfake [expr $rhor*1.3]
#set rhor [expr 2.0]
puts "Horizon: $rhor"

puts "bh.set.tcl.1"
#Cloned or computed variables (must come before other settings)
vis5d_make_expr_var $dtx "rergo=1.0+sqrt(1.0-($bhspin*cos(posh))**2.0)"
vis5d_make_expr_var $dtx "zeroergo=posr-rergo"

puts "bh.set.tcl.2"
#Isosurfaces
# not enough resolution for the below -- generates continents on BH
#vis5d_set_isosurface $ctx "posr" $rhor
vis5d_set_isosurface $ctx "posr" $rhorfake
puts "bh.set.tcl.2b"
vis5d_make_isosurface $ctx VIS5D_ALL_TIMES "posr" 0
puts "bh.set.tcl.2c"
vis5d_enable_graphics $ctx VIS5D_ISOSURF $posr VIS5D_ON

puts "bh.set.tcl.3"
#vis5d_set_isosurface $ctx "zeroergo" 0.0
#vis5d_make_isosurface $ctx VIS5D_ALL_TIMES "zeroergo" 0
#vis5d_enable_graphics $ctx VIS5D_ISOSURF $zeroergo VIS5D_ON

puts "bh.set.tcl.4"
#Isosurface colors
vis5d_set_color $dtx VIS5D_ISOSURF $ctx "posr"               0               0               0               1
#vis5d_set_color $dtx VIS5D_ISOSURF $ctx "zeroergo"               0               1               0               0.2

puts "bh.set.tcl.5"
#Isosurface color tables
vis5d_set_color_table_params $dtx VIS5D_ISOSURF $ctx "posr"       1.4000000       1.0000000       2.0000000     255.0000000       0.0000000     255.0000000     255.0000000       1.2135766      37.2887840       0.0000000 
#vis5d_set_color_table_params $dtx VIS5D_ISOSURF $ctx "zeroergo"       1.4000000       1.0000000       2.0000000     255.0000000       0.0000000     255.0000000     255.0000000      -0.9673617      36.3595085       0.0000000 
puts "bh.set.tcl done"
