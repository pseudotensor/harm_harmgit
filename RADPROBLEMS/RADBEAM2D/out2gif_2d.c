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
	  "set term gif large size 1200,500\n"
	  "set output \"%s\"\n"
	  "set size 1,1\n"
	  "set origin 0,0\n"
	  "set multiplot\n"
	  "set lmargin at screen 0.05\n"
	  "set rmargin at screen 0.4\n"
	  "set bmargin at screen .12\n"
	  "set tmargin at screen .95\n"
	  "set size .5,1\n"
	  "set xrange [%Lf:%Lf]\n"
	  "set yrange [%Lf:%Lf]\n"
	  "set xlabel \"x\"\n"
	  "set ylabel \"y\"\n"
	  "set cblabel \"\"\n"
	  "set title \"rad. energy density \" offset 0,-1\n"
#ifdef MINKOWSKI
	  "splot \"%s\" u 1:3:20 ti \"\"\n"
#else
	  //	  "set autoscale\n"
	  "splot \"%s\" u (($1)*cos($3)):(($1)*sin($3)):20 ti \"\" w l \n"
#endif

	  //	  "set autoscale\n"
	  "set lmargin at screen 0.55\n"
	  "set rmargin at screen 0.9\n"
	  "set bmargin at screen .12\n"
	  "set tmargin at screen .95\n"
	  "unset log cb\n"
	  "set cblabel \"\"\n"
	  "set title \"Flux\" offset 0,-1\n"
	  //	  "set cbrange [0:1]\n"
#ifdef MINKOWSKI
	  "splot \"%s\" u 1:3:(($21*$21+$22*$22+$23*$23)**.5) ti \"\"\n"
#else
	  	  "splot \"%s\" u (($1)*cos($3)):(($1)*sin($3)):(($21*$21+$22*$22+$23*$23)**.5/1) ti \"\" w l \n"
	  //	  "splot \"%s\" u (($1)*cos($3)):(($1)*sin($3)):(($23/$21)) ti \"\" w l \n"
#endif
	  "unset pm3d\n"
	  "set isosam 10,10\n"
	  "set format x \"\"\n"
	  "set format y \"\"\n" 
	  "set xlabel \"\"\n"
	  "set ylabel \"\"\n"
	  "set title \"\" offset 0,-1\n"
#ifdef MINKOWSKI
	  "plot \"%s\" u 1:3:($21/(($21*$21+$22*$22+$23*$23)**.5)/%Lf):($23/(($21*$21+$22*$22+$23*$23)**.5)/%Lf) every %d:%d w vectors arrowstyle 1 ti \"\"\n"
	  ,fname2,get_xb(0,0),get_xb(NX,0),get_xb(-NG,2),get_xb(NZ,2),fname,fname,
	  fname,30./(get_xb(NX,0)-get_xb(0,0)),30./(get_xb(NZ,2)-get_xb(0,2)),(int)(NX/20),(int)(NZ/20));
#else
	  "plot \"%s\" u (($1)*cos($3)):(($1)*sin($3)):(($21*cos($3)-($23)*sin($3))/(($21*cos($3)*$21*cos($3)+$23*sin($3)*$23*sin($3)+1.e-30*$20)**.5)/%Lf):(($23*cos($3)+($21)*sin($3))/(($21*cos($3)*$21*cos($3)+$23*sin($3)*$23*sin($3)+1.e-30*$20)**.5)/%Lf) every %d:%d w vectors arrowstyle 1 ti \"\"\n"
	  ,fname2,
	    
#ifndef PLOTFULLPHI
	    get_xb( -NG,0)*cos(get_xb(NZ+NG,2)),
	    get_xb(NX+NG,0),
	    sin(get_xb(-NG,2))*get_xb(NX+NG,0),
	    sin(get_xb(NZ+NG,2))*get_xb(NX+NG,0)
#else
	    get_xb(NZ+NG,2)<Pi/2. ? 0.l : -get_xb(NX+NG,0)*cos(get_xb(NZ+NG,2)) ,
	    get_xb(NX+NG,0),
	    get_xb(NZ+NG,2)<Pi ? -0.l : -get_xb(NX+NG,0)*sin(get_xb(NZ+NG,2)>1.5*Pi ? Pi/2. : get_xb(NZ+NG,2)),
	    get_xb(NZ+NG,2)<Pi/2. ? get_xb(NX+NG,0)*cos(get_xb(NZ+NG,2)) : get_xb(NX+NG,0)
#endif
	    ,fname,fname,
	    fname,30./(get_xb(NX,0)-get_xb(NX/2,0)),30./(get_xb(NX,0)-get_xb(NX/2,0)),(int)(NX/10),(int)(NZ/10));
#endif
	    

/*
  fprintf(fgnu,"\n");
  fclose(fgnu);   
  
  int i=system("gnuplot plot.gp ");
  return 0;
*/
