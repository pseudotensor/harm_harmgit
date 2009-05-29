#include "defs.h"

// jon_interp_mnewt.c is same as mnewt.c except debug statements, such as those inside debugfail>? are removed since lots of functions called in debug (at end of mnewt.c).  Also removed DEBUGPOINT stuff and mpildsum0 commands
// nrutil.c, coord.c, lubksb.c, ludcmp.c same

// decs.h, defs.h, and global.h are mostly different with some things borrowed from HARM

// newt,brodyn have choice added for analytic jacobian
// rest are same




static void interp_init(void);
static void setup_zones(void);
static void interp_readcommandlineargs(int argc, char *argv[]);
static void readdata_preprocessdata(void);
static void input_header(void);
static void output_header(void);
static void output2file_postinterpolation(void);





int main(int argc, char *argv[])
{
  int doinginterpolation;

  // Initialize interpolation
  interp_init();

  // Read command line commands
  interp_readcommandlineargs(argc, argv);

  // not doing interpolation if input and output resolutions same with no gridtype change
  // doesn't matter if refining
  doinginterpolation=!(oN1==nN1 && oN2==nN2 && oN3==nN3 && (newgridtype==GRIDTYPENOCHANGE || oldgridtype==newgridtype));

  input_header();

  // read in coordinate parameters (after inputting header that may determine defcoord, etc.)
  set_coord_parms(defcoord); // set to default
  read_coord_parms(defcoord); // read if file has updated values

  // setup grid data information (in case header info different from commandline, use header info)
  setup_zones();

  // setup new grid parameters(before waited till after refinement)
  setup_newgrid();

  // now can output header for new grid
  output_header();

  // read data and preprocess data before interpolation
  readdata_preprocessdata();

  if(immediateoutput){
    // then done!
    // flush outputs
    fflush(stderr);
    fflush(stdout);
    return(0);
  }


  // setup grid data information (pre-refinement)
  setup_zones();


  if(doinginterpolation){
    // only refine or setup new grid if doing interpolation

    // refine data
    refine_data();

    // setup boundaries and zones for old grid (post refinement)
    setup_zones();
  }



  ///////////////////
  //
  // interpolation
  //
  ///////////////////
  if(doinginterpolation==0){
    // avoid interpolation if not necessary
    fprintf(stderr,"No spatial interpolation necessary since oN1=nN1=%d  oN2=nN2=%d  oN3=nN3=%d\n",oN1,oN2,oN3);
    copy_old2new();
  }
  else{
    compute_spatial_interpolation();
  }


  // output file after interpolation
  output2file_postinterpolation();


  // flush outputs
  fflush(stderr);
  fflush(stdout);
  return(0);
}







static void interp_init(void)
{

  //  NX1=N1;
  //NX2=N2;
  //NX3=N3;
  //NX1BND=N1BND;
  //NX2BND=N2BND;
  //NX3BND=N3BND;

 NUMEPSILONPOW23=pow(NUMEPSILON,2.0/3.0);


  // some dummy assignments to make menwt.c and coord.c work, even if don't use these quantities
  ncpux1=0;
  numprocs=1;
  mycpupos[1]=0;
  mycpupos[2]=0;
  mycpupos[3]=0;
  horizoni=horizoncpupos1=0;
  didstorepositiondata=0; // so don't really need dxdxpstore, idxdxpstore, Xstore, Vstore.  So don't need to (but one can) reorder N?M as in global.storage.h
  myid=0;
  debugfail=0;
  nstroke=0;
  failed=0;
  nstep=0;
  //  N1BND=N2BND=N3BND=0;

  // fire up random number generator
  ranc(1,0);


}







static void interp_readcommandlineargs(int argc, char *argv[])
{
  int i;

  /////////////////////////////
  //
  // get arguments
  //
  if(argc < 27) {
    for(i=0;i<argc;i++){
      fprintf(stderr,"argv[%d]=%s\n",i,argv[i]);
    }
    fprintf(stderr,"args (argc=%d should be 27+ (26+ user args)): DATATYPE,INTERPTYPE,READHEADER,WRITEHEADER,oN1,oN2,oN3,refinefactor,filter,sigma,oldgridtype,newgridtype,nN1,nN2,nN3,startxc,endxc,startyc,endyc,startzc,endzc,Rin,Rout,R0,hslope,defcoord\n",argc) ; 

    fprintf(stderr,"DATATYPE:\n"
	    "0=image (byte binary only, 1 column only)\n"
	    "1=data (text only, 1=scalar 1 column)\n"
	    "2,3,4,5=correspond to orthonormal vectors v^0,v^1,v^2,v^3 (inputting 4 columns of data u^0 u^1 u^2 u^3)\n"
	    "11=corresponds to \\detg T^x1_t[EM]/sin(\\theta) (inputting 7 columns of data: u^t v^1 v^2 v^3 B^1 B^2 B^3)\n"
	    "12=get lower component (inputting 4 columns of data: u^i)\n"
	    "100+x=corresponds to inputting x-number of 4-vectors and outputting all 4-vectors in orthonormal basis without any interpolation\n"
	    );
    fprintf(stderr,"INTERPTYPE: 0=nearest 1=bi-linear 2=planar 3=bicubic\n");
    fprintf(stderr,"READHEADER: 0=false 1=true\n");
    fprintf(stderr,"WRITEHEADER: 0=false 1=true\n");
    fprintf(stderr,"oN1: old N1 grid size\n");
    fprintf(stderr,"oN2: old N2 grid size\n");
    fprintf(stderr,"oN3: old N3 grid size\n");
    fprintf(stderr,"refinefactor: 1.0=no refinement, otherwise refines image before interpolation with this factor increase in size: standard is bicubic refinement\n");
    fprintf(stderr,"filter: 0=no filter #=filter image with surrounding # pixels per pixel with sigma width\n");
    fprintf(stderr,"sigma: only used if filter!=0, then sigma of gaussian filter, usually ~ filter value\n");
    fprintf(stderr,"oldgridtype (V in GRMHD code): 0=Cartesian  1=spherical polar 2=cylindrical 3=log(z) vs. log(R), 4=x'=sin\theta log(r) z'=cos\theta log(r)\n");
    fprintf(stderr,"newgridtype (output coord system): -1=no change (and rest same as above)\n");
    fprintf(stderr,"Note: Assume if oN3>1 and oldgridtype==1, then periodic in \\phi since otherwise extrapolate values and poor at low resolution.\n");
    
    fprintf(stderr,"nN1: new N1 grid size\n");
    fprintf(stderr,"nN2: new N2 grid size\n");
    fprintf(stderr,"nN3: new N3 grid size\n");

    fprintf(stderr,"startxc: inner interp x(x-Cart,R-Cyl)\n");
    fprintf(stderr,"endxc: outer x\n");
    fprintf(stderr,"startyc: inner interp y(z-Cart,z-Cyl)\n");
    fprintf(stderr,"endyc: outer y\n");
    fprintf(stderr,"startzc: inner interp z(y-Cart,y-Cyl)\n");
    fprintf(stderr,"endzc: outer z\n");

    // some basic grid parameters, but sometimes need specific coord.c file with its parameters
    fprintf(stderr,"Rin: Inner radial edge\n");
    fprintf(stderr,"Rout: Outer radial edge\n");
    fprintf(stderr,"R0: Radial inner-grid enhancement factor\n");
    fprintf(stderr,"hslope: theta grid refinement factor\n");
    fprintf(stderr,"defcoord: which coordinate system (see coord.c)\n");
    fprintf(stderr,"extrapolate: 0 = no, 1 = yes\n");
    fprintf(stderr,"defaultvaluetype: 0 = min if scalar 0 if vector, 1 = min, 2 = max, 3 = 0.0, 4 = 1E35 for v5d missingdata\n");
    fprintf(stderr,"gdumpfilepathname : only if vector type\n");

    fprintf(stderr,"e.g.\n");
    fprintf(stderr,"~/sm/iinterp 0 0 1 1 456 456 1  1 0 0  1 0 256 512 1  1.321 40 0 40 40 0.3 0 < im0p0s0l0000.r8 > ../iimages/iim0p0s0l0000.r8\n");

    exit(0) ;
  }
  for(i=0;i<argc;i++){
    fprintf(stderr,"argv[%d]=%s\n",i,argv[i]);
  }


  


  // get args
  i=1;
  sscanf(argv[i++],"%d",&DATATYPE) ; // 0,1, etc.

  // set number of columns of data
  if(DATATYPE>1){
    immediateoutput=0;// default is not immediate in-out

    // assume always want to transform vectors correctly
    if(DATATYPE>=2 && DATATYPE<=5){
      num4vectors=1; // default
      vectorcomponent=DATATYPE-2;
      outputvartype=1; // simple vector conversion to orthonormal basis
      fprintf(stderr,"Processing vector component %d\n",vectorcomponent);
    }
    else if(DATATYPE==11){
      num4vectors=1; // default
      vectorcomponent=1;
      outputvartype=2; // Compute EM Poynting term
      fprintf(stderr,"Computing poloidal EM Polynting angular flux density\n");
    }
    else if(DATATYPE==12){
      num4vectors=1; // default
      vectorcomponent=3; // B_\\phi
      outputvartype=3; // Compute B_\\phi [conserved current]
      fprintf(stderr,"Computing B_\\phi\n");
    }
    else if(DATATYPE>=101){
      num4vectors=DATATYPE-100;
      fprintf(stderr,"Transforming %d 4-vectors to orthonormal basis\n",num4vectors);
      // force input and output grid types to be the same
      outputvartype=1; // conversion to orthonormal basis
      immediateoutput=1;
      vectorcomponent=-1; // indicates to output all components
      // only valid for DATATYPE=1 type of data (i.e. not images)
    }
    else{
      dualfprintf(fail_file,"No such DATATYPE=%d\n",DATATYPE);
    }

    // finally set to normal data type from now on using vectorcomponent or outputvartype
    DATATYPE=1;
  }
  else{
    immediateoutput=0;// default is not immediate in-out
    outputvartype=0; // scalar
    vectorcomponent=0; // scalar
    fprintf(stderr,"Processing scalar\n");
  }



  sscanf(argv[i++],"%d",&INTERPTYPE) ; // 0,1,2,3

  // set totalbc in case required
  if(INTERPTYPE==0) totalbc=0;
  else if(INTERPTYPE==1 || INTERPTYPE==2) totalbc=1;
  else if(INTERPTYPE==3) totalbc=2;
  else totalbc=MAXBC; // maximum


  sscanf(argv[i++],"%d",&READHEADER) ; // 0 or 1
  sscanf(argv[i++],"%d",&WRITEHEADER) ; // 0 or 1
  sscanf(argv[i++],"%d",&oN1) ;
  sscanf(argv[i++],"%d",&oN2) ;
  sscanf(argv[i++],"%d",&oN3) ;

  // set sizes
  totalsize[1]=oN1;
  totalsize[2]=oN2;
  totalsize[3]=oN3;
  totalzones=totalsize[1]*totalsize[2]*totalsize[3];

  // set number of bc's per dimension
  numbc[0]=0;
  numbc[1]=totalbc*(oN1>1);
  numbc[2]=totalbc*(oN2>1);
  numbc[3]=totalbc*(oN3>1);


  sscanf(argv[i++],SCANARG,&refinefactor) ;// 1.0 then no refinement, just normal old image used
  sscanf(argv[i++],"%d",&filter) ;// 0=no filter #=filter given image within surrounding # pixels per pixel with sigma
  sscanf(argv[i++],SCANARG,&sigma) ;// only used if filter!=0, then sigma of gaussian filter, usually ~ filter value
  sscanf(argv[i++],"%d",&oldgridtype) ; // 0, 1, 2, 3, and 4 currently
  sscanf(argv[i++],"%d",&newgridtype) ; // -1, and above 0+ versions
  sscanf(argv[i++],"%d",&nN1) ; // arbitrary
  sscanf(argv[i++],"%d",&nN2) ; // arbitrary
  sscanf(argv[i++],"%d",&nN3) ; // arbitrary

  //(x=R-cyl, y=Z-cyl, z=Y-cyl)
  sscanf(argv[i++],SCANARG,&startxc) ; // arbitrary
  sscanf(argv[i++],SCANARG,&endxc) ; // arbitrary
  sscanf(argv[i++],SCANARG,&startyc) ; // arbitrary
  sscanf(argv[i++],SCANARG,&endyc) ; // arbitrary
  sscanf(argv[i++],SCANARG,&startzc) ; // arbitrary
  sscanf(argv[i++],SCANARG,&endzc) ; // arbitrary

  // often other coord.c dependent stuff needed
  sscanf(argv[i++],SCANARG,&Rin) ; // could use setRin()
  sscanf(argv[i++],SCANARG,&Rout) ;
  sscanf(argv[i++],SCANARG,&R0) ;
  sscanf(argv[i++],SCANARG,&hslope) ;
  sscanf(argv[i++],"%d",&defcoord) ;

  // conditionally read-in things (all in or out)
  getgdump=0;
  if(argc>=28){ // 29 or 30
    
    if(argc<29){
      fprintf(stderr,"argc=%d insufficient-pre\n",argc);
      exit(1);
    }
    sscanf(argv[i++],"%d",&EXTRAPOLATE);
    sscanf(argv[i++],"%d",&defaultvaluetype);

    // see if doing non-scalar output that needs metric, etc.
    if(outputvartype>0){
      if(argc!=30){
	fprintf(stderr,"argc=%d insufficient-vector type\n",argc);
	exit(1);
      }
      sscanf(argv[i++],"%s",&gdumpfilename[0]);
      getgdump=1;
    }// end if doing non-scalar interpolation/output

  }// end if conditionally read-in things exist to be read-in
  else{
    EXTRAPOLATE=1; // default to extrapolate
    defaultvaluetype=0; // assume typical default
  }


  fprintf(stderr,"done reading %d arguments\n",i-1); fflush(stderr);


  if(filter && oN3!=1){
    filter=0;
    fprintf(stderr,"Turned off filter since oN3=%d\n",oN3); fflush(stderr);
  }


  if(fabs(refinefactor-1.0)>0.1 && oN3!=1){
    refinefactor=1.0;
    fprintf(stderr,"Turned off refinement since oN3=%d\n",oN3); fflush(stderr);
  }

  if( (oN3>1) && (INTERPTYPE==3 || INTERPTYPE==2) ){
    fprintf(stderr,"PLANAR and BICUBIC interpolation not setup for oN3=%d>1 -- uses nearest for k",oN3);
  }


  if(READHEADER)  jonheader=1;
  else jonheader=0;



  // process case where not interpolating and just processing input and output per cell at once
  if(immediateoutput){
    // force this even if user messed up
    newgridtype=GRIDTYPENOCHANGE;
    // force output grid size to be same as input
    nN1=oN1;
    nN2=oN2;
    nN3=oN3;
  }

}





static void setup_zones(void)
{
  int j;
 
  ////////////////////////////
  //
  // setup boundaries and zones for old grid (pre refinement)
  //
  ///////////////////////////
    
  if(immediateoutput==0){
    totalsize[1]=oN1;
    totalsize[2]=oN2;// in case changed, setup for set_points
    totalsize[3]=oN3;
  }
  else{
    // then get totalsize from header and oN? from command line since can be different if only converting part of original dataset
  }

  fprintf(stderr,"ts1=%d ts2=%d ts3=%d :: oN1=%d oN2=%d oN3=%d\n",totalsize[1],totalsize[2],totalsize[3],oN1,oN2,oN3);
  fprintf(stderr,"sx1=%g sx2=%g sx3=%g\n",startx[1],startx[2],startx[3]);


  //  fprintf(stderr,"defcoord=%d\n",defcoord);
  set_points();

  fprintf(stderr,"after set_points(): sx1=%g sx2=%g sx3=%g\n",startx[1],startx[2],startx[3]);

  SLOOPA(j) dX[j]=dx[j];// just convert from coord.c code result

  startpos[1]=startpos[2]=startpos[3]=0;
  Xmax[1] = startx[1]+dX[1]*(FTYPE)oN1;
  Xmax[2] = startx[2]+dX[2]*(FTYPE)oN2;
  Xmax[3] = startx[3]+dX[3]*(FTYPE)oN3;
  
  // below 2 not really used (yet at least) -- needed for setRin() but that function isn't used
  a=spin;
  Rhor=rhor_calc(0);

  fprintf(stderr,"defcoord=%d Rin=%g R0=%g\n",defcoord,Rin,R0);

  fprintf(stderr,"startx,Xmax,dX: %g %g %g :: %g %g %g :: %g %g %g\n",startx[1],Xmax[1],dX[1],startx[2],Xmax[2],dX[2],startx[3],Xmax[3],dX[3]) ; fflush(stderr);

}




static void readdata_preprocessdata(void)
{
  int i,j,k;
  unsigned char tempuc;


  if(getgdump){
    gdumpin=fopen(gdumpfilename,"rt");
    if(gdumpin==NULL){
      fprintf(stderr,"No such gdump file %s\n",gdumpfilename);
      exit(1);
    }
    else{
      // skip first line assuming it's a header line
      while(fgetc(gdumpin)!='\n'); // go past end of line
    }
  }


  ///////////////////////////////////
  //
  // read data
  //
  if(DATATYPE==0){
    imagedata=0; // says treat as image and fit output between 0-255
    if(DOUBLEWORK) DATATYPE=1;

    /* make arrays for images */
    if(!DOUBLEWORK){
      oldimage0 = c3matrix(-numbc[1]+0,oN1-1+numbc[1],-numbc[2]+0,oN2-1+numbc[2],-numbc[3]+0,oN3-1+numbc[3]) ;
      newimage  = c3matrix(-numbc[1]+0,nN1-1+numbc[1],-numbc[2]+0,nN2-1+numbc[2],-numbc[3]+0,nN3-1+numbc[3]) ;
    }
    else{
      olddata0 = f3matrix(-numbc[1]+0,oN1-1+numbc[1],-numbc[2]+0,oN2-1+numbc[2],-numbc[3]+0,oN3-1+numbc[3]) ;   // olddata0[i][j][k]
      newdata  = f3matrix(-numbc[1]+0,nN1-1+numbc[1],-numbc[2]+0,nN2-1+numbc[2],-numbc[3]+0,nN3-1+numbc[3]) ;   // newdata[i][j][k]
    }
    /* read in old image */
    if(jonheader){
      // skip 4 lines
      for(i=1;i<=4;i++) while(fgetc(stdin)!='\n');
    }
    fprintf(stderr,"reading image\n"); fflush(stderr);
    
    // read order and write order same, so good image output
    totalmin=BIG;
    totalmax=-BIG;
    for(k=0;k<oN3;k++){
      for(j=0;j<oN2;j++){
	for(i=0;i<oN1;i++){
	  if(DOUBLEWORK){
	    fread(&tempuc, sizeof(unsigned char), 1, stdin) ;
	    olddata0[i][j][k]=(FTYPE)tempuc;
	    if(olddata0[i][j][k]>totalmax) totalmax=olddata0[i][j][k];
	    if(olddata0[i][j][k]<totalmin) totalmin=olddata0[i][j][k];
	  }
	  else{
	    fread(&oldimage0[i][j][k], sizeof(unsigned char), 1, stdin) ;
	    if(oldimage0[i][j][k]>totalmax) totalmax=oldimage0[i][j][k];
	    if(oldimage0[i][j][k]<totalmin) totalmin=oldimage0[i][j][k];
	  }
	}
      }
    }

    // filter or not
    if(filter){
      // filter not setup for periodic bc
      fprintf(stderr,"filter\n");
      gaussian_filter(filter,sigma,oN1,oN2,oN3,oldimage0,olddata0);      
    }

    if(PERIODICINPHI && oN3>1 && oldgridtype==GRIDTYPESPC){
      // then fill boundary cells for good interpolation rather than ad hoc extrapolation that leaves feature at \phi=0=2\pi boundary
      for(j=0;j<oN2;j++){
	for(i=0;i<oN1;i++){
	  if(DOUBLEWORK){
	    for(k=-numbc[3];k<0;k++) olddata0[i][j][k]=olddata0[i][j][k+oN3];
	    for(k=oN3;k<oN3+numbc[3];k++) olddata0[i][j][k]=olddata0[i][j][k-oN3];
	  }
	  else{
	    for(k=-numbc[3];k<0;k++) oldimage0[i][j][k]=oldimage0[i][j][k+oN3];
	    for(k=oN3;k<oN3+numbc[3];k++) oldimage0[i][j][k]=oldimage0[i][j][k-oN3];
	  }
	}
      }
    }// end if periodic


  }
  else if(DATATYPE==1){
    imagedata=1; // says treat as data

    if(immediateoutput==0){
      olddata0 = f3matrix(-numbc[1]+0,oN1-1+numbc[1],-numbc[2]+0,oN2-1+numbc[2],-numbc[3]+0,oN3-1+numbc[3]) ;   // olddata0[i][j][k]
      newdata  = f3matrix(-numbc[1]+0,nN1-1+numbc[1],-numbc[2]+0,nN2-1+numbc[2],-numbc[3]+0,nN3-1+numbc[3]) ;   // newdata[i][j][k]
    }


    if(immediateoutput==0){
      totalmin=BIG;
      totalmax=-BIG;
      // read it (Note the loop order!)
      LOOPOLDDATA{

	if(outputvartype==0){
	  fscanf(stdin,SCANARG,&olddata0[i][j][k]) ;
	}
	else{
	  compute_preprocess(gdumpin, &olddata0[i][j][k]);
	}
	if(olddata0[i][j][k]>totalmax) totalmax=olddata0[i][j][k];
	if(olddata0[i][j][k]<totalmin) totalmin=olddata0[i][j][k];
      
      }// end LOOPOLDDATA

      if(filter){
	// filter not setup for periodic bc
	fprintf(stderr,"filter\n");
	gaussian_filter(filter,sigma,oN1,oN2,oN3,oldimage0,olddata0);
      }

      if(PERIODICINPHI && oN3>1 && oldgridtype==GRIDTYPESPC){
	// then fill boundary cells for good interpolation rather than ad hoc extrapolation that leaves feature at \phi=0=2\pi boundary
	for(j=0;j<oN2;j++){
	  for(i=0;i<oN1;i++){
	    for(k=-numbc[3];k<0;k++) olddata0[i][j][k]=olddata0[i][j][k+oN3];
	    for(k=oN3;k<oN3+numbc[3];k++) olddata0[i][j][k]=olddata0[i][j][k-oN3];
	  }
	}
      }// end if PERIODIC
    }
    else{
      compute_preprocess(gdumpin, NULL);
    }

  }



  if(immediateoutput==0){
    ///////////////
    //
    // set default value of interpolation when nothing to interpolate if not extrapolating
    //
    ///////////////

    if(defaultvaluetype==0){
      if(outputvartype==0) defaultvalue=totalmin;
      else defaultvalue=0.0; // vector-like things otherwise around 0
    }
    else if(defaultvaluetype==1){
      defaultvalue=totalmin;
    }
    else if(defaultvaluetype==2){
      defaultvalue=totalmax;
    }
    else if(defaultvaluetype==3){
      defaultvalue=0.0;
    }
    else if(defaultvaluetype==4){
      defaultvalue=1E35; // for V5D missing data
      fprintf(stderr,"Using V5D missing data for defaultvalue=%g\n",defaultvalue);
    }
  }
}











static void input_header(void)
{
  FTYPE ftemp;


  /* read in old data (only designed for 1 column right now)*/
  if(READHEADER){
    // assumes header really has ALL this info (could tell user how many entries on header with wc and compare against desired.
    // GODMARK
    // If using gammie.m's interpsingle, must keep interpsingle macro's header output up-to-date
    fscanf(stdin, SCANHEADER,
	   &t,&totalsize[1],&totalsize[2],&totalsize[3],&startx[1],&startx[2],&startx[3],&dX[1],&dX[2],&dX[3],&readnstep,&gam,&spin,&R0,&Rin,&Rout,&hslope,&dt,&defcoord,&MBH,&QBH);



    realnstep=(long)readnstep;
    nstep=realnstep;
    if((totalsize[1]!=oN1)||(totalsize[2]!=oN2)||(totalsize[3]!=oN3)){
      fprintf(stderr,"expected %d x %d x %d and got %d x %d x %d resolution -- ok if totalsize sets grid and oN? sets data size in file itself\n",oN1,oN2,oN3,totalsize[1],totalsize[2],totalsize[3]);
    }
    while(fgetc(stdin)!='\n'); // go past end of line
  }

  // print header from file
  fprintf(stderr, PRINTSCANHEADER,
	  t,totalsize[1],totalsize[2],totalsize[3],startx[1],startx[2],startx[3],dX[1],dX[2],dX[3],readnstep,gam,spin,R0,Rin,Rout,hslope,dt,defcoord,MBH,QBH);


}


static void output_header(void)
{
  FTYPE ftemp;

 
  //////////////////////////
  //
  // output header
  //
  //////////////////////////
  fprintf(stderr,"header:\n");
  fprintf(stderr, "OLD: %22.16g :: %d %d %d :: %22.16g %22.16g %22.16g :: %22.16g %22.16g %22.16g :: %ld %22.16g %22.16g %22.16g %22.16g %22.16g %22.16g %22.16g %d %22.16g %22.16g\n",
	  t,oN1,oN2,oN3,startx[1],startx[2],startx[3],dX[1],dX[2],dX[3],realnstep,gam,spin,R0,Rin,Rout,hslope,dt,defcoord,MBH,QBH);
  fprintf(stderr, "NEW: %d %d %d :: %22.16g %22.16g %22.16g :: %22.16g %22.16g %22.16g\n",nN1,nN2,nN3,Xmax[1],Xmax[2],Xmax[3],fakedxc,fakedyc,fakedzc);
  fprintf(stderr, "OTHER: %22.16g %22.16g %22.16g %22.16g\n",fakeRin,dxc,dyc,dzc);
   

  if(WRITEHEADER){
    if(DATATYPE==1){
      // print out a header
      ftemp=0.0;
      fprintf(stdout, "%22.16g %d %d %d %22.16g %22.16g %22.16g %22.16g %22.16g %22.16g %ld %22.16g %22.16g %22.16g %22.16g %22.16g %22.16g %22.16g %d %22.16g %22.16g\n",
	      t, nN1, nN2, nN3, startxc, startyc, startzc, fakedxc,fakedyc,fakedzc,realnstep,gam,spin,ftemp,endxc,endyc,hslope,dt,defcoord,MBH,QBH);
    }
  }

}









static void output2file_postinterpolation(void)
{
  int i,j,k;
  unsigned char uctemp;
  FTYPE ftemp;

  // OUTPUT TO FILE
  fprintf(stderr,"Output to file\n"); fflush(stderr);

  // in principle could output in different order if wanted
  if(DATATYPE==0){
    for(k=0;k<nN3;k++) for(j=0;j<nN2;j++)      for(i=0;i<nN1;i++) {
      fwrite(&newimage[i][j][k], sizeof(unsigned char), 1, stdout) ;
    }
  }
  else if(DATATYPE==1){
    if(imagedata==0){
      for(k=0;k<nN3;k++) for(j=0;j<nN2;j++)      for(i=0;i<nN1;i++) {
	ftemp=newdata[i][j][k];
	//	if(ftemp<0.0) ftemp=0.0;
	//if(ftemp>255.0) ftemp=255.0;
	//uctemp=(unsigned char)ftemp;
	uctemp = FLOAT2IMAGE(ftemp);
	fwrite(&uctemp, sizeof(unsigned char), 1, stdout) ;
      }
    }
    else{
      if(sizeof(FTYPE)==sizeof(double)){
	for(k=0;k<nN3;k++) for(j=0;j<nN2;j++)      for(i=0;i<nN1;i++) {
	  //fprintf(stderr,"write: i=%d j=%d newdata=%22.16g\n",i,j,newdata[i][j][k]); fflush(stderr);
	  fprintf(stdout,"%22.16g\n",newdata[i][j][k]) ;
	}
      }
      else if(sizeof(FTYPE)==sizeof(float)){
	for(k=0;k<nN3;k++) for(j=0;j<nN2;j++)      for(i=0;i<nN1;i++) {
	  fprintf(stdout,"%15.7g\n",newdata[i][j][k]) ;
	}
      }
    }
  }
}




////////////////////////////////////////
//
// IMAGE WRITE FUNCTION
//
////////////////////////////////////////

void writeimage(char * name, unsigned char *** image,int nx, int ny, int nz)
{
  FILE * out;
  int i,j,k;

  if((out=fopen(name,"wb"))==NULL){
    fprintf(stderr,"Cannot open %s\n",name);
    exit(1);
  }

  for(k=0;k<nz;k++){
    for(j=0;j<ny;j++){
      for(i=0;i<nx;i++){
	fwrite(&image[i][j][k], sizeof(unsigned char), 1, out) ;
	//      fprintf(out, "%c",(unsigned char)((int)image[i][j][k]));
      }
    }
  }
  fclose(out);

}

