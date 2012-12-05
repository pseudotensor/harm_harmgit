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
	  "set cntrparam level discrete 0.01,.1,1,10\n"
	  "splot \"%s\" u 1:2:14 w l\n"
	  "unset dgrid3d\n"
	  "unset log z\n"
	  "unset table\n"
	  "unset contour\n"
	  "unset surface\n"

	  "set view map\n"
	  "set pm3d\n"
	  "unset key\n"
	  "set style line 1 lt 1 lw 1 lc 3\n"
	  "set style arrow 1 head nofilled size screen 0.002,35 ls 1\n"
	  "set palette model RGB rgbformulae 7,5,15\n"
	  "set palette model RGB rgbformulae 30,31,32\n"
	  "set palette model RGB rgbformulae 21,22,23\n"



	  "unset surface\n"
	  "set term gif large size 1200,150\n"
	  "set output \"%s\"\n"
	  "set size 1,1\n"
	  "set origin 0,0\n"
	  "set multiplot\n"
	  "set lmargin at screen 0.05\n"
	  "set rmargin at screen 0.4\n"
	  "set bmargin at screen .15\n"
	  "set tmargin at screen .85\n"
	  
	  "set size .5,1\n"
	  "set xrange [%Lf:%Lf]\n"
	  "set yrange [%Lf:%Lf]\n"
	  //	  "set log cb\n"
	  //"set cbrange [1.e6:5.e6]\n"
	  "set xlabel \"x\"\n"
	  "set ylabel \"y\"\n"
	  "set cblabel \"\"\n"
	  "set title \"T (colors), rho (contours)\"\n"
	  //	  "set dgrid3d 100,100,16\n"
	  
	  "splot \"%s\" u 1:2:20 ti \"\"\n"

	  "set format x \"\"\n"
	  "set format y \"\"\n"
	  "set xlabel \"\"\n"
	  "set ylabel \"\"\n"
	  "set title \"\"\n"
	  "unset pm3d\n"
	  "unset surface\n"
	  "plot \"table.gp\" w l lc 3 lw 3\n"
	  "set pm3d\n"
	  "set format x \"%%.2f\"\n"
	  "set format y \"%%.2f\"\n"
	  "set autoscale\n"
	  "set xrange [%Lf:%Lf]\n"
	  "set yrange [%Lf:%Lf]\n"

	  "set lmargin at screen 0.55\n"
	  "set rmargin at screen 0.9\n"
	  "set bmargin at screen .15\n"
	  "set tmargin at screen .85\n"
	  //	  "set dgrid3d 180,90,8\n"
	  "set cblabel \"\"\n"
	  "unset log cb\n"
	  "set title \"flux (colors,arrows), rho (contours)\"\n"

	  //	  "set cbrange [0:2]\n"	  
	  
	  "splot \"%s\" u 1:2:(($21*$21+$22*$22)**.5) ti \"\"\n"
	  //"splot \"%s\" u 1:2:(($11*$11)**.5) ti \"\"\n"
	  "unset dgrid3d\n"
	  
	  //"splot \"%s\" u 1:2:16 ti \"\"\n"


	  "unset pm3d\n"
	  "set isosam 10,10\n"
	  "set format x \"\"\n"
	  "set format y \"\"\n" 
	  "set xlabel \"\"\n"
	  "set ylabel \"\"\n"
	  "plot \"%s\" u 1:2:($21/(($21*$21+$22*$22)**.5)/35):($22/(($21*$21+$22*$22)**.5)/35) every %d:%d w vectors arrowstyle 1 ti \"\"\n"
	  "unset pm3d\n"
	  "unset surface\n"
	  "plot \"table.gp\" w l lc 3 lw 3\n"
	  "set pm3d\n"
	  //fname,fname2,get_xb(-2,0),get_xb(NX,0),get_xb(0,1),get_xb(NY,1),fname,get_xb(0,0),get_xb(NX,0),get_xb(0,1),get_xb(NY,1),fname,fname,(int)(NX/20),(int)(NY/8));
	  ,fname,fname2,get_x(-5,0),get_xb(NX,0),get_xb(0,1),get_xb(NY,1),fname,get_x(-5,0),get_xb(NX,0),get_xb(0,1),get_xb(NY,1),fname,fname,(int)(NX/20),(int)(NY/8));

	    

/*
  fprintf(fgnu,"\n");
  fclose(fgnu);   
  
  int i=system("gnuplot plot.gp ");
  return 0;
*/
