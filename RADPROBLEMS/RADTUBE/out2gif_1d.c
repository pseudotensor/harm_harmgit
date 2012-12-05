//KORAL - bc.c
//gives the gnuplot script used for plotting gif files on the go - 1 dimension
//**********************
//included in:
//int convert_out2gif_1d(char *fname,char *fname2,int niter,ldouble t)
//from misc.c
//**********************

/****************************************/
/****************************************/
/****************************************/
fprintf(fgnu,
	"set style line 1 lw 2 lc 1 lt 1\n"
	"set style line 2 lw 2 lc 2 lt 1\n"
	"set style line 3 lw 2 lc 3 lt 1\n"
	"set style line 4 lw 2 lc 4 lt 1\n"
	"set style line 10 lw 2 lc 3 lt 3\n"
	"set style line 11 lw 2 lc 4 lt 3\n"
	"set style line 12 lw 2 lc 5 lt 3\n"
	"set style line 13 lw 2 lc 0 lt 3\n"
	"set style line 14 lw 2 lc 7 lt 3\n"
	"set term gif large size 1200,600\n"
	"set output \"%s\"\n"
	"set size 1,1\n"
	"set origin 0,0\n"
	"set multiplot\n"
	"set autoscale\n" 
	"set label \"t=%.2Le (%.2Le s)\" at screen .48, .98\n"

	"set lmargin at screen 0.07\n"
	"set rmargin at screen 0.33\n"
	"set bmargin at screen .07\n"
	"set tmargin at screen .42\n"
	"set autoscale\n"
	"set xrange [%Lf:%Lf]\n"
	//	  "unset log y\n"
	"set format x \"%%.1f\"\n"
	"set format y \"%%.1e\"\n" 
	"set xlabel \"\"\n"
	"set ylabel \"\"\n"
	"plot \"%s\" u 1:($14) w lp ls 1 pt 7 ps .5 ti \"rho\"\n"

	"set lmargin at screen 0.40\n"
	"set rmargin at screen 0.66\n"
	"set bmargin at screen .07\n"
	"set tmargin at screen .42\n"
	//	  "unset log y\n"
	"set format x \"%%.1f\"\n"
	"set format y \"%%.1e\"\n" 
	"set xlabel \"\"\n"
	"set ylabel \"\"\n"
	//	  "plot \"%s\" u 1:27 w lp ls 2 pt 7 ps .5  ti \"tau_abs\", \"%s\" u 1:26 w lp ls 3 pt 7 ps .5  ti \"tau_tot\"\n"
	"plot \"%s\" u 1:20 w lp ls 2 pt 7 ps .5  ti \"E_rad\"\n"

	"set lmargin at screen 0.73\n"
	"set rmargin at screen 0.99\n"
	"set bmargin at screen .07\n"
	"set tmargin at screen .42\n"
	//	  "unset log y\n"
	"set format x \"%%.1f\"\n"
	"set format y \"%%.1e\"\n" 
	"set xlabel \"\"\n"
	"set ylabel \"\"\n"
	"plot \"%s\" u 1:($21+1.e-80) w lp ls 2 pt 7 ps .5  ti \"Fx\"\n"

	"set lmargin at screen 0.07\n"
	"set rmargin at screen 0.33\n"
	"set bmargin at screen .5\n"
	"set tmargin at screen .9\n"
	//	  "unset log y\n"
	"set format x \"\"\n"
	"set format y \"%%.1e\"\n" 
	"set xlabel \"\"\n"
	"set ylabel \"\"\n"
	//	  "set log y\n"
	"plot \"%s\" u 1:15 w lp ls 3 pt 7 ps .5 ti \"u_int\", \"%s\" u 1:20 w lp ls 7 pt 7 ps .5 ti \"E_rad\"\n"
	//	  "unset log y\n"

	"set lmargin at screen 0.40\n"
	"set rmargin at screen 0.66\n"
	"set bmargin at screen .5\n"
	"set tmargin at screen .9\n"
	"set format x \"\"\n"
	"set format y \"%%.1e\"\n" 
	"set xlabel \"\"\n"
	"set ylabel \"\"\n"
	"plot \"%s\" u 1:24 w p ls 4 pt 7 ti \"Tgas\", \"%s\" u 1:25 w l lc 9 lw 2 ti \"Trad\"\n"

	"set lmargin at screen 0.73\n"
	"set rmargin at screen 0.99\n"
	"set bmargin at screen .5\n"
	"set tmargin at screen .9\n"
	"set format x \"\"\n"
	"set format y \"%%.1e\"\n" 
	"set xlabel \"\"\n"
	"set ylabel \"\"\n"
	"plot \"%s\" u 1:($16) w lp ls 4 pt 7 ti \"vx\""
	,fname2,t,t/CCC,get_xb(-NG,0),get_xb(NX+NG,0),fname,fname,fname,fname,fname,fname,fname,fname);

/****************************************/
/****************************************/
/****************************************/

