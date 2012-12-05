//int
//convert_out2gif_1d(char *fname,char *fname2,int niter,ldouble t)
//{
//  FILE *fgnu=fopen("plot.gp","w");
//  char bufor[50];

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
	  "set term gif large size 800,600\n"
	  "set output \"%s\"\n"
	  "set size 1,1\n"
	  "set origin 0,0\n"
	  "set multiplot\n"
	  "set size .5,.5\n"

	  "set origin 0,.5\n"
	  "set xrange [%Lf:%Lf]\n"
	  //	  "set yrange [0.:1.1]\n"
	  "set xlabel \"x\"\n"
	  "unset log y\n"
	  "set ylabel \"rho\"\n"
	  "plot \"%s\" u 1:($14) w p ls 1 pt 7 ps .5 ti \"i=%d\", \"%s\" u 1:4 w l ls 2 ti \"anal\"\n"
	  //	  "set autoscale\n"
	  "set origin 0,0\n"
	  "unset log y\n"
	  //	  "set yrange [-120000:120000]\n"
	  "set ylabel \"velocity\"\n"
	  "plot \"%s\" u 1:($16/$14) w p ls 2 pt 7 ps .5  ti \"i=%d\"\n"
	  //"plot \"%s\" u 1:13 w p ls 2 pt 7 ps .5  ti \"i=%d\", \"%s\" u 1:14 w p ls 4 pt 7 ps .5  ti \"i=%d\"\n"

	  "set origin .5,.5\n"
	  //	  "set yrange [1e11:5e12]\n"
	  "unset log y\n"
	  "set ylabel \"unknown\"\n"
	  "plot \"%s\" u 1:5 w p ls 3 pt 7 ps .5 ti \"i=%d\"\n"

	  "set origin .5,0\n"
	  "unset log y\n"
	  //	  "set yrange [-2e6:2e6]\n"
	  //	  "set autoscale\n"
	  "set ylabel \"unknown\"\n"
	  "plot \"%s\" u 1:($9) w p ls 4 pt 7 ti \"i=%d\", \"%s\" u 1:9 w l ls 3 ti \"anal\"\n"
	  ,fname2,get_xb(-NG,0),get_xb(NX+NG,0),fname,niter,fname,fname,niter,fname,niter,fname,niter,fname);
 


/****************************************/
/****************************************/
/****************************************/

//  fprintf(fgnu,"\n");
//  fclose(fgnu);   
//  
//  int i=system("gnuplot plot.gp ");
//  return 0;
//}
