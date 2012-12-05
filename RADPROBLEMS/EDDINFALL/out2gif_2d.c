//int
//convert_out2gif_2d(char *fname,char *fname2,int niter,ldouble t)
//{
//  FILE *fgnu=fopen("plot.gp","w");
//  char bufor[50];
  fprintf(fgnu,
	  "set view map\n"
	  "set pm3d\n"
	  "unset surface\n"
	  "unset key\n"
	  "set style line 1 lt 1 lw 3 lc 3\n"
	  "set style line 2 lt 1 lw 3 lc 2\n"
	  "set style line 11 lt 3 lw 1 lc -1\n"
	  "set style line 21 lt 3 lw 1 lc -1\n"
	  "set style arrow 1 head nofilled size screen 0.002,35 ls 1\n"
	  "set palette model RGB rgbformulae 7,5,15\n"
	  "set palette model RGB rgbformulae 30,31,32\n"
	  "set palette model RGB rgbformulae 21,22,23\n"
	  "set palette model RGB rgbformulae 23,28,3\n"

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
	  "set cbrange [0.00001:1]\n"
	  //	  "set xlabel \"x\"\n"
	  "set ylabel \"y\"\n"
	  "set cblabel \"\"\n"
	  "set log cb\n"
	  "set title \"rho\"\n"
#ifdef MINKOWSKI
	  "splot \"%s\" u (($1)):(($2)):14 ti \"\" w l ls 1\n"
#else 
	  "splot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):14 ti \"\" w l ls 1\n"
#endif	  
	  "set lmargin at screen 0.55\n"
	  "set rmargin at screen 0.9\n"
	  "set bmargin at screen .45\n"
	  "set tmargin at screen .92\n"
	  "set cblabel \"\"\n"
	  "set ylabel \"\"\n"
	  "set title \"radial velocity\"\n"
	  "set cbrange [0.10:.5]\n"
#ifdef MINKOWSKI
	  "splot \"%s\" u (($1)):(($2)):16 ti \"\" w l ls 1\n"
#else
	  "splot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):(-$16) ti \"\" w l ls 1\n"
#endif
	  "unset pm3d\n"
	  "set lmargin at screen 0.08\n"
	  "set rmargin at screen 0.43\n"
	  "set bmargin at screen .10\n"
	  "set tmargin at screen .40\n"
	  "set title \"\"\n"
	  "set xlabel \"x\"\n"
	  "set yrange [0.001:1]\n"
	  //	  "set autoscale\n"
	  "set log y\n"
	  "plot \"%s\" u 1:(%f/($1*$1*(2./$1*(1.))**.5)) ti \"\" w l ls 11, \"%s\" u 1:14 ti \"\" w l ls 1 \n"

	  "set lmargin at screen 0.55\n"
	  "set rmargin at screen 0.92\n"
	  "set bmargin at screen .10\n"
	  "set tmargin at screen .40\n"
	  "set title \"\"\n"
	  "set yrange [0.09:1]\n"
	  "set autoscale\n"
	  "set log y\n"
	  "plot \"%s\" u 1:((2./$1)**.5*(1.-2./$1)*(1/(1 - 2/$1))**.5) ti \"\" w l ls 21,\"%s\" u 1:(-$16) ti \"\" w l ls 2 \n"
	  //"plot \"%s\" u 1:($15) ti \"\" w l ls 2, \"%s\" u 1:($15) ti \"\" w l ls 2 \n"
	  ,fname2,
#ifdef MINKOWSKI
	  get_xb(-NG,0),
	  get_xb(NX+NG,0),
	  get_xb(-NG,1),
	  get_xb(NY+NG,1),
#else
	  -.05*get_xb(NX+NG,0),
	  1.05*get_xb(NX+NG,0),
	  -.05*get_xb(NX+NG,0),
	  1.05*get_xb(NX+NG,0),
#endif
	  fname,fname,fname,PAR_D,fname,fname,fname); 
	    

/*
  fprintf(fgnu,"\n");
  fclose(fgnu);   
  
  int i=system("gnuplot plot.gp ");
  return 0;
*/
