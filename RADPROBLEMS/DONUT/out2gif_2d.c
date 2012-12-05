//int
//convert_out2gif_2d(char *fname,char *fname2,int niter,ldouble t)
//{
//  FILE *fgnu=fopen("plot.gp","w");
//  char bufor[50];
  fprintf(fgnu,
	  "set table \"table.gp\"\n"
	  "set contour base\n"
	  "unset surface\n"
	  "set log z\n"
	  "set cntrparam level discrete .01,.05,.1,.3,.5,.7,.9,1.1,1.3,1.5,1.7 \n"
	  "splot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):24 w l\n"
	  "unset dgrid3d\n"
	  "unset log z\n"
	  "unset table\n"
	  "unset contour\n"
	  "unset surface\n"

	  "set view map\n"
	  "set pm3d\n"
	  "unset surface\n"
	  "unset key\n"
	  "set style line 1 lt 1 lw 1 lc 3\n"
	  "set style line 11 lt 1 lw 2 lc 6\n"
	  "set style line 2 lt 1 lw 2 lc 2\n"
	  "set style line 3 lt 1 lw 2 lc 2\n"
	  "set style line 21 lt 3 lw 1 lc -1\n"
	  "set style arrow 1 head nofilled size screen 0.002,35 ls 3\n"

	  "set palette model RGB rgbformulae 7,5,15\n"
	  "set palette model RGB rgbformulae 30,31,32\n"
	  "set palette model RGB rgbformulae 21,22,23\n"
	  "set palette model RGB rgbformulae 23,28,3\n"
	  "set palette model RGB rgbformulae 7,8,9\n"
	  "set palette model RGB rgbformulae 6,3,21\n"
	  "set palette model RGB rgbformulae 35,3,9\n"

	  "unset log cb\n"
	  "set term gif large size 1000,700\n"
	  "set output \"%s\"\n"
	  "set size 1,1\n"
	  "set origin 0,0\n"
	  "set multiplot\n"
	  "set lmargin at screen 0.08\n"
	  "set rmargin at screen 0.43\n"
	  "set bmargin at screen .45\n"
	  "set tmargin at screen .92\n"
	  "set size .5,1\n"
	  "set autoscale\n"
	  "set xrange [%Lf:%Lf]\n"
	  "set yrange [%Lf:%Lf]\n"
	  "unset log cb\n"
	  //	  "set cbrange [0.005:.35]\n"
	  "set ylabel \"y\"\n"
	  "set cblabel \"\"\n"
	  //	  "set log cb\n"
	  "set title \"rho\"\n"
	  
	  "splot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):($14) ti \"\" w l ls 1\n"

	  "set ylabel \"\"\n"
	  "unset tics\n"
	  "unset border\n"
	  "unset pm3d\n"
	  "unset surface\n"
	  "plot \"table.gp\" w l ls 1\n"
	  "set pm3d\n"
	  "set tics\n"
	  "set border\n"

	  "set table \"table.gp\"\n"
	  "set contour base\n"
	  "set log z\n"
	  "set cntrparam level discrete .01,.05,.1,.3,.5,.7,.9,1.1,1.3,1.5,1.7 \n"
	  "splot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):14 w l\n"
	  "unset log z\n"
	  "unset table\n"
	  "unset contour\n"

	  "set ylabel \"\"\n"
	  "unset tics\n"
	  "unset border\n"
	  "unset pm3d\n"
	  "unset surface\n"
	  "plot \"table.gp\" w l ls 11\n"
	  "set pm3d\n"
	  "set tics\n"
	  "set border\n"	  

	  "set lmargin at screen 0.55\n"
	  "set rmargin at screen 0.9\n"
	  "set bmargin at screen .45\n"
	  "set tmargin at screen .92\n"
	  "set cblabel \"\"\n"
	  "set ylabel \"y\"\n"
	  "set title \"mass flux / poloidal velocity\"\n"
	  "unset log cb\n"
	  "set log cb\n"
	  "set cbrange [0.001:.5]\n"
	  "splot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):(($16*$16+$17*$17)**.5*$14) ti \"\" w l ls 1\n"
	  "unset pm3d\n"

	  "set isosam 10,10\n"
	  "set ylabel \"\"\n"
	  "unset tics\n"
	  "unset border\n"
	  "unset log cb\n"
	  "plot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):(($16*sin($2)+$17*cos($2))/(%Lf)):(($17*sin($2)+$16*cos($2))/(%Lf)) every %d:%d w vectors arrowstyle 1 ti \"\"\n"


	  "set table \"table.gp\"\n"
	  "set contour base\n"
	  "set log z\n"
	  "set cntrparam level discrete .01,.05,.1,.3,.5,.7,.9,1.1,1.3,1.5,1.7 \n"
	  "splot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):24 w l\n"
	  "unset log z\n"
	  "unset table\n"
	  "unset contour\n"

	  "set ylabel \"\"\n"
	  "unset tics\n"
	  "unset border\n"
	  "unset pm3d\n"
	  "unset surface\n"
	  "plot \"table.gp\" w l ls 1\n"
	  "set pm3d\n"
	  "set tics\n"
	  "set border\n"	  
	  
	  "unset log cb\n"
	  "unset colorbox\n"
	  "set lmargin at screen 0.08\n"
	  "set rmargin at screen 0.43\n"
	  "set bmargin at screen .10\n"
	  "set tmargin at screen .40\n"
	  "set title \"\"\n"
	  "set xlabel \"x\"\n"
	  "set ylabel \"density\"\n"
	  "set autoscale\n"
	  "unset log y\n"
	  "plot \"%s\" u 1:14 ti \"\" w l ls 11  \n"

	  "set lmargin at screen 0.55\n"
	  "set rmargin at screen 0.9\n"
	  "set bmargin at screen .10\n"
	  "set tmargin at screen .40\n"
	  "set title \"\"\n"
	  "set yrange [0.09:.5]\n"
	  "set ylabel \"poloidal velocity\"\n"
	  "set autoscale\n"
	  "unset log y\n"
	  "plot \"%s\" u 1:(($16*$16)**.5) ti \"\" w l ls 2 \n"
	  ,fname,fname2,
	  -.02*get_xb(NX,0),
	  1.02*get_xb(NX,0),
	  -.02*get_xb(NX,0),
	  1.02*get_xb(NX,0),
	  fname,fname,fname,
	  fname,1./(get_xb(NX,0)/11)*.5,1./(get_xb(NX,0)/11)*.5,NX/11+1,NY/11+1,
	  fname,fname,fname);  
	    

/*
  fprintf(fgnu,"\n");
  fclose(fgnu);   
  
  int i=system("gnuplot plot.gp ");
  return 0;
*/
