//int
//convert_out2gif_2d(char *fname,char *fname2,int niter,ldouble t)
//{
//  FILE *fgnu=fopen("plot.gp","w");
//  char bufor[50];
 fprintf(fgnu,
	  "set view map\n"
	  "set pm3d\n"
	  "unset surface\n"
	  //	  "set contour base\n"
	  "set cntrparam level discrete 2,4,6,8\n"
	  "unset key\n"
	  "set style line 1 lt 1 lw 1 lc 7\n"
	  "set style arrow 1 head nofilled size screen 0.002,35 ls 1\n"
	  "set palette model RGB rgbformulae 7,5,15\n"
	  "set palette model RGB rgbformulae 30,31,32\n"
	  "set palette model RGB rgbformulae 21,22,23\n"
	  "unset log cb\n"
	  "set term gif large size 1200,500\n"
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
	  "set cbrange [0.9:11]\n"
	  "set xlabel \"x\"\n"
	  "set ylabel \"y\"\n"
	  "set cblabel \"\"\n"
	  "set title \"gamma rho\"\n"
#ifdef CYLINDRICAL
	  "splot \"%s\" u (($1)*cos($3)):(($1)*sin($3)):14 ti \"\" w l ls 1\n"
#else
	  "splot \"%s\" u 1:3:14 ti \"\" w l\n"
#endif

	  
	  "set lmargin at screen 0.55\n"
	  "set rmargin at screen 0.9\n"
	  "set bmargin at screen .15\n"
	  "set tmargin at screen .85\n"
	  "set cbrange [0:.8]\n"
	  "set cblabel \"\"\n"
	  "set cntrparam level discrete 0.1,.3,.5,.7\n"
	  "set title \"velocity\"\n"
#ifdef CYLINDRICAL
	  "splot \"%s\" u (($1)*cos($3)):(($1)*sin($3)):(($16*sin($3)+$18*cos($3))/$14) ti \"\" w l \n"
#else
	  "splot \"%s\" u 1:3:($16/$14) ti \"\" w l ls 1\n"
#endif
	  "unset pm3d\n"
	  "set isosam 10,10\n"
	  "set format x \"\"\n"
	  "set format y \"\"\n"
	  "set xlabel \"\"\n"
	  "set ylabel \"\"\n"

#ifdef CYLINDRICAL
	  //	  ,fname2,get_xb(-NG,0),get_xb(NX+NG,0),sin(get_xb(-NG,2))*get_xb(-NG,0),sin(get_xb(NZ+NG,2))*get_xb(NX+NG,0),fname,fname);
	  ,fname2,
	  get_x(NX/2,0)-(sin(get_xb(-NG,2))*get_xb(-NG,0)),
	  get_x(NX/2,0)+(sin(get_xb(-NG,2))*get_xb(-NG,0)),
	  sin(get_xb(-NG,2))*get_xb(NX+NG,0),
	  sin(get_xb(NZ+NG,2))*get_xb(NX+NG,0),
	  fname,fname);
#else
	  ,fname2,get_xb(-NG,0),get_xb(NX+NG,0),get_xb(-NG,2),get_xb(NZ+NG,2),fname,fname);
#endif

	    

/*
  fprintf(fgnu,"\n");
  fclose(fgnu);   
  
  int i=system("gnuplot plot.gp ");
  return 0;
*/
