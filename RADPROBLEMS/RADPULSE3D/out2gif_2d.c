//int
//convert_out2gif_2d(char *fname,char *fname2,int niter,ldouble t)
//{
//  FILE *fgnu=fopen("plot.gp","w");
//  char bufor[50];
  fprintf(fgnu,
	  "set view map\n"
	  "set pm3d\n"
	  "unset key\n"
	  "set style line 1 lt 1 lw 2 lc 2\n"
	  "set style arrow 1 head nofilled size screen 0.002,35 ls 1\n"
	  "set palette model RGB rgbformulae 7,5,15\n"
	  "set palette model RGB rgbformulae 30,31,32\n"
	  "set palette model RGB rgbformulae 21,22,23\n"
	  "set palette model RGB rgbformulae 21,22,23\n"
	  "set palette model RGB rgbformulae 35,3,9\n"
	  "unset surface\n"
	  "set term gif large size 600,500\n"
	  "set output \"%s\"\n"
	  "set size 1,1\n"
	  "set origin 0,0\n"
	  "set log cb\n"
	  "set lmargin at screen 0.10\n"
	  "set rmargin at screen 0.8\n"
	  "set bmargin at screen .12\n"
	  "set tmargin at screen .95\n"
	  "set xrange [%Lf:%Lf]\n"
	  "set yrange [%Lf:%Lf]\n"
	  "set cbrange [1.e-37:1.e-33]\n"
	  "set xlabel \"x\"\n"
	  "set ylabel \"y\"\n"
	  "set cblabel \"\"\n"
	  "set title \"\" offset 0,-1\n"
	  "splot \"%s\" u 1:3:20 ti \"\"\n"
	  "unset pm3d\n"
	  "set format x \"\"\n"
	  "set format y \"\"\n" 
	  "set xlabel \"\"\n"
	  "set ylabel \"\"\n"
	  "plot \"%s\" u 1:3:($21/(($21*$21+$22*$22+$23*$23)**.5)/%Lf):($23/(($21*$21+$22*$22+$23*$23)**.5)/%Lf) every %d:%d w vectors arrowstyle 1 ti \"\"\n"
	  ,fname2,get_xb(0,0),get_xb(NX,0),get_xb(0,1),get_xb(NY,1),fname,fname,
	  30./(get_xb(NX,0)-get_xb(0,0)),30./(get_xb(NY,1)-get_xb(0,1)),(int)(NX/20),(int)(NZ/20));
	    

/*
  fprintf(fgnu,"\n");
  fclose(fgnu);   
  
  int i=system("gnuplot plot.gp ");
  return 0;
*/
