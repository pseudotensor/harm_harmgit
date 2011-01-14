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

static void old_usage(int argc, int basicargcnum);
static void usage(int argc, int basicargcnum);

static void interpret_commandlineresults(void);
void interpret_commandlineresults_subpart1(void);

static void old_parse_commandline(int argc, char *argv[]);
static void parse_commandline(int argc, char *argv[]);

static void defaultoptions(void);

static int post_coordsetup(void);

static void readelement(FILE *input, FTYPE *datain);




int main(int argc, char *argv[])
{
  int doinginterpolation;
  
  // Initialize interpolation
  interp_init();

  // Read and process command line commands
  interp_readcommandlineargs(argc, argv);

  input_header();

  // read in coordinate parameters (after inputting header that may determine defcoord, etc.)
  set_coord_parms(defcoord); // set to default
  read_coord_parms(defcoord); // read if file has updated values

  // do things once header read-in and coordinates setup
  doinginterpolation=post_coordsetup();

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
    fprintf(stderr,"No spatial interpolation necessary since oN0=nN0=%d oN1=nN1=%d  oN2=nN2=%d  oN3=nN3=%d\n",oN0,oN1,oN2,oN3);
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

  bytesize=sizeof(unsigned char);
  intsize=sizeof(int);
  longintsize=sizeof(long int);
  longlongintsize=sizeof(long long int);
  floatsize=sizeof(float);
  doublesize=sizeof(double);
  longdoublesize=sizeof(long double);


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




  // set default options
  defaultoptions();


  /////////////////////////////
  //
  // parse command-line arguments
  //
  /////////////////////////////
  parse_commandline(argc, argv); // will be switching to switches


  if(VERBOSITY>=2){
    // report arguments
    for(i=0;i<argc;i++){
      fprintf(stderr,"argv[%d]=%s\n",i,argv[i]);
    }
  }


  /////////////////////////////
  //
  // process results from command-line arguments into normal/expected parameters setup
  //
  /////////////////////////////
  interpret_commandlineresults();


}




// setup oldgrid zones
static void setup_zones(void)
{
  int j;
 
  ////////////////////////////
  //
  // setup boundaries and zones for old grid (pre refinement)
  //
  ///////////////////////////
    
  if(immediateoutput==0){
    totalsize[0]=oN0;
    totalsize[1]=oN1;
    totalsize[2]=oN2;// in case changed, setup for set_points
    totalsize[3]=oN3;
  }
  else{
    // then get totalsize from header and oN? from command line since can be different if only converting part of original dataset
  }

  fprintf(stderr,"ts0=%d ts1=%d ts2=%d ts3=%d :: oN0=%d oN1=%d oN2=%d oN3=%d\n",totalsize[0],totalsize[1],totalsize[2],totalsize[3],oN0,oN1,oN2,oN3);
  fprintf(stderr,"sx0=%g sx1=%g sx2=%g sx3=%g\n",startx[0],startx[1],startx[2],startx[3]);



  //////////////
  //
  // setup startx[] and dx[] and Diffx[] and endx[] and dV and dVF
  //
  //////////////
  set_points();


  /////////////
  //
  // check if user did something funny
  //
  ////////////
  DLOOPA(j){
    if(fabs(Diffx[j])<SMALL){
      fprintf(stderr,"User: You have startx and endx the same.  Did you interp after already reading an interp file and forget to load the original dump?\n");
      exit(1);
    }
  }


  fprintf(stderr,"after set_points(): sx0=%g sx1=%g sx2=%g sx3=%g\n",startx[0],startx[1],startx[2],startx[3]);

  DLOOPA(j) dX[j]=dx[j];// just convert from coord.c code result

  startpos[0]=startpos[1]=startpos[2]=startpos[3]=0;
  Xmax[0] = startx[0]+dX[0]*(FTYPE)oN0;
  Xmax[1] = startx[1]+dX[1]*(FTYPE)oN1;
  Xmax[2] = startx[2]+dX[2]*(FTYPE)oN2;
  Xmax[3] = startx[3]+dX[3]*(FTYPE)oN3;
  
  // below 2 not really used (yet at least) -- needed for setRin() but that function isn't used
  a=spin;
  Rhor=rhor_calc(0);

  fprintf(stderr,"defcoord=%d Rin=%g R0=%g dofull2pi=%d\n",defcoord,Rin,R0,dofull2pi);

  fprintf(stderr,"startx,Xmax,dX: %g %g %g :: %g %g %g :: %g %g %g :: %g %g %g\n",startx[0],Xmax[0],dX[0],startx[1],Xmax[1],dX[1],startx[2],Xmax[2],dX[2],startx[3],Xmax[3],dX[3]) ; fflush(stderr);

}




// read and process input data
static void readdata_preprocessdata(void)
{
  int h,i,j,k;
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


  ///////////////////////////////////
  //
  // DATATYPE==0
  //
  // only images with 1 column of byte size
  ///////////////////////////////////
  if(DATATYPE==0){
    imagedata=0; // says treat as image and fit output between 0-255
    if(DOUBLEWORK) DATATYPE=1;

    /* make arrays for images */
    if(!DOUBLEWORK){
      oldimage0 = c4matrix(-numbc[0]+0,oN0-1+numbc[0],-numbc[1]+0,oN1-1+numbc[1],-numbc[2]+0,oN2-1+numbc[2],-numbc[3]+0,oN3-1+numbc[3]) ;
      newimage  = c4matrix(-numbc[0]+0,nN0-1+numbc[0],-numbc[1]+0,nN1-1+numbc[1],-numbc[2]+0,nN2-1+numbc[2],-numbc[3]+0,nN3-1+numbc[3]) ;
    }
    else{
      olddata0 = f4matrix(-numbc[0]+0,oN0-1+numbc[0],-numbc[1]+0,oN1-1+numbc[1],-numbc[2]+0,oN2-1+numbc[2],-numbc[3]+0,oN3-1+numbc[3]) ;   // olddata0[h][i][j][k]
      newdata  = f4matrix(-numbc[0]+0,nN0-1+numbc[0],-numbc[1]+0,nN1-1+numbc[1],-numbc[2]+0,nN2-1+numbc[2],-numbc[3]+0,nN3-1+numbc[3]) ;   // newdata[h][i][j][k]
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
    for(h=0;h<oN0;h++){ // for files (reading-writing), time is slowest index
      for(k=0;k<oN3;k++){
	for(j=0;j<oN2;j++){
	  for(i=0;i<oN1;i++){
	    if(DOUBLEWORK){
	      fread(&tempuc, sizeof(unsigned char), 1, stdin) ;
	      olddata0[h][i][j][k]=(FTYPE)tempuc;
	      if(olddata0[h][i][j][k]>totalmax) totalmax=olddata0[h][i][j][k];
	      if(olddata0[h][i][j][k]<totalmin) totalmin=olddata0[h][i][j][k];
	    }
	    else{
	      fread(&oldimage0[h][i][j][k], sizeof(unsigned char), 1, stdin) ;
	      if(oldimage0[h][i][j][k]>totalmax) totalmax=oldimage0[h][i][j][k];
	      if(oldimage0[h][i][j][k]<totalmin) totalmin=oldimage0[h][i][j][k];
	    }
	  }
	}
      }
    }

    // filter or not
    if(filter){
      // filter not setup for periodic bc
      fprintf(stderr,"filter\n");
      gaussian_filter(filter,sigma,oN0,oN1,oN2,oN3,oldimage0,olddata0);      
    }

    /////////////
    //
    // set boundary conditions (as if scalars)
    //
    ///////////// 
    if(BOUNDARYEXTRAP==1){
      // lower and upper h
      for(i=0;i<oN1;i++){
	for(j=0;j<oN2;j++){
	  for(k=0;k<oN3;k++){
	    if(DOUBLEWORK){
	      for(h=-numbc[0];h<0;h++) olddata0[h][i][j][k]=olddata0[0][i][j][k];
	      for(h=oN0;h<oN0+numbc[0];h++) olddata0[h][i][j][k]=olddata0[oN0-1][i][j][k];
	    }
	    else{
	      for(h=-numbc[0];h<0;h++) oldimage0[h][i][j][k]=oldimage0[0][i][j][k];
	      for(h=oN0;h<oN0+numbc[0];h++) oldimage0[h][i][j][k]=oldimage0[oN0-1][i][j][k];
	    }
	  }
	}
      }
      // lower and upper i
      for(h=0;h<oN0;h++){
	for(j=0;j<oN2;j++){
	  for(k=0;k<oN3;k++){
	    if(DOUBLEWORK){
	      for(i=-numbc[1];i<0;i++) olddata0[h][i][j][k]=olddata0[h][0][j][k];
	      for(i=oN1;i<oN1+numbc[1];i++) olddata0[h][i][j][k]=olddata0[h][oN1-1][j][k];
	    }
	    else{
	      for(i=-numbc[1];i<0;i++) oldimage0[h][i][j][k]=oldimage0[h][0][j][k];
	      for(i=oN1;i<oN1+numbc[1];i++) oldimage0[h][i][j][k]=oldimage0[h][oN1-1][j][k];
	    }
	  }
	}
      }
      // lower and upper j
      for(h=0;h<oN0;h++){
	for(i=0;i<oN1;i++){
	  for(k=0;k<oN3;k++){
	    if(DOUBLEWORK){
	      for(j=-numbc[2];j<0;j++) olddata0[h][i][j][k]=olddata0[h][i][0][k];
	      for(j=oN2;j<oN2+numbc[2];j++) olddata0[h][i][j][k]=olddata0[h][i][oN2-1][k];
	    }
	    else{
	      for(j=-numbc[2];j<0;j++) oldimage0[h][i][j][k]=oldimage0[h][i][0][k];
	      for(j=oN2;j<oN2+numbc[2];j++) oldimage0[h][i][j][k]=oldimage0[h][i][oN2-1][k];
	    }
	  }
	}
      }
      // lower and upper k
      for(h=0;h<oN0;h++){
	for(j=0;j<oN2;j++){
	  for(i=0;i<oN1;i++){
	    if(DOUBLEWORK){
	      for(k=-numbc[3];k<0;k++) olddata0[h][i][j][k]=olddata0[h][i][j][0];
	      for(k=oN3;k<oN3+numbc[3];k++) olddata0[h][i][j][k]=olddata0[h][i][j][oN3-1];
	    }
	    else{
	      for(k=-numbc[3];k<0;k++) oldimage0[h][i][j][k]=oldimage0[h][i][j][0];
	      for(k=oN3;k<oN3+numbc[3];k++) oldimage0[h][i][j][k]=oldimage0[h][i][j][oN3-1];
	    }
	  }
	}
      }
    }


    if(PERIODICINPHI && oN3>1 && oldgridtype==GRIDTYPESPC){
      // then fill boundary cells for good interpolation rather than ad hoc extrapolation that leaves feature at \phi=0=2\pi boundary
      for(h=0;h<oN0;h++){
	for(j=0;j<oN2;j++){
	  for(i=0;i<oN1;i++){
	    if(DOUBLEWORK){
	      for(k=-numbc[3];k<0;k++) olddata0[h][i][j][k]=olddata0[h][i][j][k+oN3];
	      for(k=oN3;k<oN3+numbc[3];k++) olddata0[h][i][j][k]=olddata0[h][i][j][k-oN3];
	    }
	    else{
	      for(k=-numbc[3];k<0;k++) oldimage0[h][i][j][k]=oldimage0[h][i][j][k+oN3];
	      for(k=oN3;k<oN3+numbc[3];k++) oldimage0[h][i][j][k]=oldimage0[h][i][j][k-oN3];
	    }
	  }
	}
      }
    }// end if periodic
    
  }
  ///////////////////////////////////
  //
  // DATATYPE==1
  //
  ///////////////////////////////////
  else if(DATATYPE==1){
    imagedata=1; // says treat as data

    if(immediateoutput==0){
      olddata0 = f4matrix(-numbc[0]+0,oN0-1+numbc[0],-numbc[1]+0,oN1-1+numbc[1],-numbc[2]+0,oN2-1+numbc[2],-numbc[3]+0,oN3-1+numbc[3]) ;   // olddata0[h][i][j][k]
      newdata  = f4matrix(-numbc[0]+0,nN0-1+numbc[0],-numbc[1]+0,nN1-1+numbc[1],-numbc[2]+0,nN2-1+numbc[2],-numbc[3]+0,nN3-1+numbc[3]) ;   // newdata[h][i][j][k]
    }


    if(immediateoutput==0){
      totalmin=BIG;
      totalmax=-BIG;
      // read it (Note the loop order!) (see global.jon_interp.h)  time is slowest index for reading and writing files
      LOOPOLDDATA{

	if(outputvartype==0){
	  readelement(stdin,&olddata0[h][i][j][k]);
	}
	else{
	  compute_preprocess(outputvartype,gdumpin, &olddata0[h][i][j][k]);
	}
	if(olddata0[h][i][j][k]>totalmax) totalmax=olddata0[h][i][j][k];
	if(olddata0[h][i][j][k]<totalmin) totalmin=olddata0[h][i][j][k];
	
      }// end LOOPOLDDATA
      

      if(filter){
	// filter not setup for periodic bc
	fprintf(stderr,"filter\n");
	gaussian_filter(filter,sigma,oN0,oN1,oN2,oN3,oldimage0,olddata0);
      }

      /////////////
      //
      // set boundary conditions (as if scalars)
      //
      ///////////// 
      if(BOUNDARYEXTRAP==1){
	// lower and upper h
	for(i=0;i<oN1;i++){
	  for(j=0;j<oN2;j++){
	    for(k=0;k<oN3;k++){
	      for(h=-numbc[0];h<0;h++) olddata0[h][i][j][k]=olddata0[0][i][j][k];
	      for(h=oN0;h<oN0+numbc[0];h++) olddata0[h][i][j][k]=olddata0[oN0-1][i][j][k];
	    }
	  }
	}
	// lower and upper i
	for(h=0;h<oN0;h++){
	  for(j=0;j<oN2;j++){
	    for(k=0;k<oN3;k++){
	      for(i=-numbc[1];i<0;i++) olddata0[h][i][j][k]=olddata0[h][0][j][k];
	      for(i=oN1;i<oN1+numbc[1];i++) olddata0[h][i][j][k]=olddata0[h][oN1-1][j][k];
	    }
	  }
	}
	// lower and upper j
	for(h=0;h<oN0;h++){
	  for(i=0;i<oN1;i++){
	    for(k=0;k<oN3;k++){
	      for(j=-numbc[2];j<0;j++) olddata0[h][i][j][k]=olddata0[h][i][0][k];
	      for(j=oN2;j<oN2+numbc[2];j++) olddata0[h][i][j][k]=olddata0[h][i][oN2-1][k];
	    }
	  }
	}
	// lower and upper k
	for(h=0;h<oN0;h++){
	  for(j=0;j<oN2;j++){
	    for(i=0;i<oN1;i++){
	      for(k=-numbc[3];k<0;k++) olddata0[h][i][j][k]=olddata0[h][i][j][0];
	      for(k=oN3;k<oN3+numbc[3];k++) olddata0[h][i][j][k]=olddata0[h][i][j][oN3-1];
	    }
	  }
	}
      }

      // override and always have boundary cells if periodic in x3
      if(PERIODICINPHI && oN3>1 && oldgridtype==GRIDTYPESPC){
	// then fill boundary cells for good interpolation rather than ad hoc extrapolation that leaves feature at \phi=0=2\pi boundary
	for(h=0;h<oN0;h++){
	  for(j=0;j<oN2;j++){
	    for(i=0;i<oN1;i++){
	      for(k=-numbc[3];k<0;k++) olddata0[h][i][j][k]=olddata0[h][i][j][k+oN3];
	      for(k=oN3;k<oN3+numbc[3];k++) olddata0[h][i][j][k]=olddata0[h][i][j][k-oN3];
	    }
	  }
	}
      }// end if PERIODIC
    }
    else{
      compute_preprocess(outputvartype,gdumpin, NULL);
    }

  }



  if(immediateoutput==0){
    ///////////////
    //
    // set default value of interpolation when nothing to interpolate if not extrapolating
    //
    ///////////////

    if(defaultvaluetype==0){
      if(outputvartype==0 || (outputvartype==1||outputvartype==2) && vectorcomponent==0) defaultvalue=totalmin;
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


    if(VERBOSITY>=1){
      fprintf(stderr,"defaultvalue=%21.15g\n",defaultvalue);
    }

  }
}


// read single element from file in text or binary format and any element C type
void readelement(FILE *input, FTYPE *datain)
{
  
  if     (binaryinput==0 && strcmp(inFTYPE,"b")==0   ){ int dumi;             fscanf(stdin,"%d",&dumi) ;    *datain=(FTYPE)dumi;       }
  else if(binaryinput==0 && strcmp(inFTYPE,"i")==0   ){ int dumi;             fscanf(stdin,"%d",&dumi) ;    *datain=(FTYPE)dumi;       }
  else if(binaryinput==0 && strcmp(inFTYPE,"li")==0  ){ long int dumli;       fscanf(stdin,"%ld",&dumli);   *datain=(FTYPE)dumli;      }
  else if(binaryinput==0 && strcmp(inFTYPE,"lli")==0 ){ long long int dumlli; fscanf(stdin,"%lld",&dumlli); *datain=(FTYPE)dumlli;     }
  else if(binaryinput==0 && strcmp(inFTYPE,"f")==0   ){ float dumf;           fscanf(stdin,"%f",&dumf) ;    *datain=(FTYPE)dumf;       }
  else if(binaryinput==0 && strcmp(inFTYPE,"d")==0   ){ double dumd;          fscanf(stdin,"%lf",&dumd);    *datain=(FTYPE)dumd;       }
  else if(binaryinput==0 && strcmp(inFTYPE,"ld")==0  ){ long double dumld;    fscanf(stdin,"%Lf",&dumld);   *datain=(FTYPE)dumld;      }

  if     (binaryinput==1 && strcmp(inFTYPE,"b")==0   ){ int dumi;             fread(&dumi,bytesize,1,input);          *datain=(FTYPE)dumi;       }
  else if(binaryinput==1 && strcmp(inFTYPE,"i")==0   ){ int dumi;             fread(&dumi,intsize,1,input);           *datain=(FTYPE)dumi;       }
  else if(binaryinput==1 && strcmp(inFTYPE,"li")==0  ){ long int dumli;       fread(&dumli,longintsize,1,input);      *datain=(FTYPE)dumli;      }
  else if(binaryinput==1 && strcmp(inFTYPE,"lli")==0 ){ long long int dumlli; fread(&dumlli,longlongintsize,1,input); *datain=(FTYPE)dumlli;     }
  else if(binaryinput==1 && strcmp(inFTYPE,"f")==0   ){ float dumf;           fread(&dumf,floatsize,1,input);         *datain=(FTYPE)dumf;       }
  else if(binaryinput==1 && strcmp(inFTYPE,"d")==0   ){ double dumd;          fread(&dumd,doublesize,1,input);        *datain=(FTYPE)dumd;       }
  else if(binaryinput==1 && strcmp(inFTYPE,"ld")==0  ){ long double dumld;    fread(&dumld,longdoublesize,1,input);   *datain=(FTYPE)dumld;      }

}







// input header from dump file
// does not consider time component dimension
static void input_header(void)
{
  FTYPE ftemp;


  /* read in old data (only designed for 1 column right now)*/
  if(READHEADER){
    // assumes header really has ALL this info (could tell user how many entries on header with wc and compare against desired.
    // GODMARK
    // If using gammie.m's interpsingle, must keep interpsingle macro's header output up-to-date
    fscanf(stdin, SCANHEADER,
	   &tdump,&totalsize[1],&totalsize[2],&totalsize[3],&startx[1],&startx[2],&startx[3],&dX[1],&dX[2],&dX[3],&readnstep,&gam,&spin,&R0,&Rin,&Rout,&hslope,&dtdump,&defcoord,&MBH,&QBH,&is,&ie,&js,&je,&ks,&ke,&whichdump,&whichdumpversion,&numcolumns);


    // set other things not set by header, but not really used right now
    startx[0]=0;  // will be set by set_points() to be this way, but put here to user knows this has to be true
    //endx[0]=???


    realnstep=(long)readnstep;
    nstep=realnstep;
    if(oldparse){ // only a problem at this point if oldparse==1 and READHEADER==1 since old parse line needs to have correct oN?
      if((totalsize[1]!=oN1)||(totalsize[2]!=oN2)||(totalsize[3]!=oN3)){
	fprintf(stderr,"expected %d x %d x %d and got %d x %d x %d resolution -- ok if totalsize sets grid and oN? sets data size in file itself\n",oN1,oN2,oN3,totalsize[1],totalsize[2],totalsize[3]);
      }
    }
    while(fgetc(stdin)!='\n'); // go past end of line (so can add stuff to end of header, but won't be funneled to new interpolated file)
  }

  // print header from file
  fprintf(stderr,"PRINTSCANHEADER\n");
  fprintf(stderr, PRINTSCANHEADER,
	  tdump,totalsize[1],totalsize[2],totalsize[3],startx[1],startx[2],startx[3],dX[1],dX[2],dX[3],readnstep,gam,spin,R0,Rin,Rout,hslope,dtdump,defcoord,MBH,QBH,is,ie,js,je,ks,ke,whichdump,whichdumpversion,numcolumns);



}



// output dump header
// does not consider time component dimension
static void output_header(void)
{
  FTYPE ftemp;

 
  //////////////////////////
  //
  // output header
  //
  //////////////////////////
  fprintf(stderr,"header:\n");
  fprintf(stderr, "OLD: %22.16g :: %d %d %d :: %22.16g %22.16g %22.16g :: %22.16g %22.16g %22.16g :: %ld %22.16g %22.16g %22.16g %22.16g %22.16g %22.16g %22.16g %d %22.16g %22.16g %d %d %d %d %d %d %d %d %d\n",
	  tdump,oN1,oN2,oN3,startx[1],startx[2],startx[3],dX[1],dX[2],dX[3],realnstep,gam,spin,R0,Rin,Rout,hslope,dtdump,defcoord,MBH,QBH,is,ie,js,je,ks,ke,whichdump,whichdumpversion,numcolumns);
  fprintf(stderr, "NEW: %d %d %d :: %22.16g %22.16g %22.16g :: %22.16g %22.16g %22.16g\n",nN1,nN2,nN3,Xmax[1],Xmax[2],Xmax[3],fakedxc,fakedyc,fakedzc);
  fprintf(stderr, "OTHER: %22.16g %22.16g %22.16g %22.16g\n",fakeRin,dxc,dyc,dzc);
   

  if(WRITEHEADER){
    if(DATATYPE==1){
      // print out a header
      ftemp=0.0;
      fprintf(stdout, "%22.16g %d %d %d %22.16g %22.16g %22.16g %22.16g %22.16g %22.16g %ld %22.16g %22.16g %22.16g %22.16g %22.16g %22.16g %22.16g %d %22.16g %22.16g %d %d %d %d %d %d %d %d %d\n",
	      tdump, nN1, nN2, nN3, startxc, startyc, startzc, fakedxc,fakedyc,fakedzc,realnstep,gam,spin,ftemp,endxc,endyc,hslope,dtdump,defcoord,MBH,QBH,is,ie,js,je,ks,ke,whichdump,whichdumpversion,numcolumns);
    }
  }

}


// do things once header read-in and coordinates setup
int post_coordsetup(void)
{
  int doinginterpolation;

  if(READHEADER){
    // if read header, then have some things that user didn't have to set on command line if using new parse mode, so set them up
    oN1=totalsize[1];
    oN2=totalsize[2];
    oN3=totalsize[3];
    // also set is R0, Rin, Rout, hslope, defcoord, gam, spin, etc. (see input_header() ).
  }

  // not doing interpolation if input and output resolutions same with no gridtype change
  // doesn't matter if refining
  doinginterpolation=!(oN0==nN0 && oN1==nN1 && oN2==nN2 && oN3==nN3 && (newgridtype==GRIDTYPENOCHANGE || oldgridtype==newgridtype));


  return(doinginterpolation);

}





// output newgrid data to file
static void output2file_postinterpolation(void)
{
  int h,i,j,k;
  unsigned char uctemp;
  FTYPE ftemp;


  // OUTPUT TO FILE
  fprintf(stderr,"Output to file\n"); fflush(stderr);



  // in principle could output in different order if wanted
  if(DATATYPE==0){
    for(h=0;h<nN0;h++)  for(k=0;k<nN3;k++) for(j=0;j<nN2;j++)      for(i=0;i<nN1;i++) {
      fwrite(&newimage[h][i][j][k], sizeof(unsigned char), 1, stdout) ;
    }
  }
  else if(DATATYPE==1){
    if(imagedata==0){
      for(h=0;h<nN0;h++)  for(k=0;k<nN3;k++) for(j=0;j<nN2;j++)      for(i=0;i<nN1;i++) {
	ftemp=newdata[h][i][j][k];
	//	if(ftemp<0.0) ftemp=0.0;
	//if(ftemp>255.0) ftemp=255.0;
	//uctemp=(unsigned char)ftemp;
	uctemp = FLOAT2IMAGE(ftemp);
	fwrite(&uctemp, sizeof(unsigned char), 1, stdout) ;
      }
    }
    else{
      if(sizeof(FTYPE)==sizeof(double)){
	for(h=0;h<nN0;h++)  for(k=0;k<nN3;k++) for(j=0;j<nN2;j++)      for(i=0;i<nN1;i++) {
	  //fprintf(stderr,"write: i=%d j=%d newdata=%22.16g\n",i,j,newdata[h][i][j][k]); fflush(stderr);
	  fprintf(stdout,"%22.16g\n",newdata[h][i][j][k]) ;
	}
      }
      else if(sizeof(FTYPE)==sizeof(float)){
	for(h=0;h<nN0;h++) for(k=0;k<nN3;k++) for(j=0;j<nN2;j++)      for(i=0;i<nN1;i++) {
	  fprintf(stdout,"%15.7g\n",newdata[h][i][j][k]) ;
	}
      }
    }
  }
}




// IMAGE WRITE FUNCTION
void writeimage(char * name, unsigned char ****image, int nt, int nx, int ny, int nz)
{
  FILE * out;
  int h,i,j,k;

  if((out=fopen(name,"wb"))==NULL){
    fprintf(stderr,"Cannot open %s\n",name);
    exit(1);
  }

  for(h=0;h<nt;h++){
    for(k=0;k<nz;k++){
      for(j=0;j<ny;j++){
	for(i=0;i<nx;i++){
	  fwrite(&image[h][i][j][k], sizeof(unsigned char), 1, out) ;
	  //      fprintf(out, "%c",(unsigned char)((int)image[h][i][j][k]));
	}
      }
    }
  }
  fclose(out);

}


void gcov_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE *gcovinfunc, FTYPE *gcovpertinfunc)
{

  // filler function so can use metric_tools.c
}

void gcon_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE *gcov, FTYPE *gcon)
{

  // filler function so can use metric_tools.c
}

void get_geometry(int ii, int jj, int kk, int pp, struct of_geom *geom)
{
  // filler function so can use metric_tools.c

  // assume time coordinate (related to oN0 and nN0) do not affect anything in metric calculations (i.e. only coordinates are affected)
  geom->i=ii;
  geom->j=jj;
  geom->k=kk;
  geom->p=pp;

}


// filler function so can use metric_tools.c
void eomfunc_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE *EOMFUNCNAME)
{
  int pl,pliter;

  PLOOP(pliter,pl) EOMFUNCASSIGN(pl)=1.0;


}


void assign_eomfunc(struct of_geom *geom, FTYPE *EOMFUNCNAME)
{
  int pl,pliter;

  PLOOP(pliter,pl) geom->EOMFUNCMAC(pl)=geom->gdet;
}


// report usage of program to user and exit
void old_new(int argc, int basicargcnum)
{


}




// old report usage of program to user and exit
void old_usage(int argc, int basicargcnum)
{


  fprintf(stderr,"args (argc=%d should be %d+ (%d+ user args)): DATATYPE,INTERPTYPE, READHEADER,WRITEHEADER, oN0,oN1,oN2,oN3, refinefactor,filter,sigma, oldgridtype,newgridtype, nN0,nN1,nN2,nN3, starttc,endtc,startxc,endxc,startyc,endyc,startzc,endzc, Rin,Rout,R0,hslope,defcoord,dofull2pi, <starttdata,endtdata>,  tnrdegrees, <extrapolate,defaultvaluetype,gdumpfilepathname>\n",argc,basicargcnum,basicargcnum-1);

  fprintf(stderr,"DATATYPE:\n"
	  "0=image (input&output: byte binary only, 1 column only)\n"
	  "1=data (input&output: text only, 1=scalar 1 column)\n"
	  "2,3,4,5=correspond to output of orthonormal vectors v^0,v^1,v^2,v^3 (inputting all 4 columns of data u^0 u^1 u^2 u^3)\n"
	  "6,7,8,9=correspond to output of orthonormal vectors v^0,v^1,v^2,v^3 (inputting all 4 columns of data u_0 u_1 u_2 u_3)\n"
	  "11=corresponds to output of \\detg T^x1_t[EM]/sin(\\theta) (inputting all 7 columns of data: u^t v^1 v^2 v^3 B^1 B^2 B^3)\n"
	  "12=output lower component (inputting all 4 columns of data: u^i)\n"
	  "100+x=corresponds to inputting x-number of 4-vectors and outputting all 4-vectors in orthonormal basis without any interpolation\n"
	  );
  fprintf(stderr,"INTERPTYPE: 0=nearest 1=bi-linear 2=planar 3=bicubic\n");
  fprintf(stderr,"READHEADER: 0=false 1=true\n");
  fprintf(stderr,"WRITEHEADER: 0=false 1=true\n");

  fprintf(stderr,"oN0: old N0 grid size\n");
  fprintf(stderr,"oN1: old N1 grid size\n");
  fprintf(stderr,"oN2: old N2 grid size\n");
  fprintf(stderr,"oN3: old N3 grid size\n");

  fprintf(stderr,"refinefactor: 1.0=no refinement, otherwise refines image before interpolation with this factor increase in size: standard is bicubic refinement\n");
  fprintf(stderr,"filter: 0=no filter #=filter image with surrounding # pixels per pixel with sigma width\n");
  fprintf(stderr,"sigma: only used if filter!=0, then sigma of gaussian filter, usually ~ filter value\n");

  fprintf(stderr,"oldgridtype (V in GRMHD code): 0=Cartesian  1=spherical polar 2=cylindrical 3=log(z) vs. log(R), 4=x'=sin\\theta log(r) z'=cos\\theta log(r) 5=Cartesian w/ light travel time accounting\n");
  fprintf(stderr,"newgridtype (output coord system): -1=no change (and rest same as above)\n");
  fprintf(stderr,"Note: Assume if oN3>1 and oldgridtype==1, then periodic in \\phi since otherwise extrapolate values and poor at low resolution.\n");
    
  fprintf(stderr,"nN0: new N0 grid size\n");
  fprintf(stderr,"nN1: new N1 grid size\n");
  fprintf(stderr,"nN2: new N2 grid size\n");
  fprintf(stderr,"nN3: new N3 grid size\n");

  fprintf(stderr,"starttc: inner interp t(t-Cart,t-Cyl t=time)\n");
  fprintf(stderr,"endtc: outer t\n");
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
  fprintf(stderr,"dofull2pi: whether to do full 2pi or not (see coord.c)\n");

  // below are only inputted when oN0>1 || nN0>1
  fprintf(stderr,"starttdata: time 4D input dumps start (assume uniformly spaced in time, and corresponds to actual time when data exists, not FACE values, but CENT in terms of internal interpolation routines)\n");
  fprintf(stderr,"endtdata: time 4D input dumps end\n");
  fprintf(stderr,"tnrdegrees: angle [in degrees] between observer and z-axis of original grid.  0 degrees gives no transformation.  20degrees means data rotated 20degrees from z-axis towards x-axis around y-axis.\n");

  // below 1 are separately optional
  fprintf(stderr,"extrapolate: 0 = no, 1 = yes\n");
  fprintf(stderr,"defaultvaluetype: 0 = min if scalar 0 if vector, 1 = min, 2 = max, 3 = 0.0, 4 = 1E35 for v5d missingdata\n");
  // below is optional but requires above 2 to be read-in
  fprintf(stderr,"gdumpfilepathname : only if vector type (DATATYPE=(e.g.) 2,3,4,5,11,12,100...\n");

  fprintf(stderr,"e.g.\n");
  fprintf(stderr,"~/sm/iinterp 0 0 1 1 456 456 1  1 0 0  1 0 256 512 1  1.321 40 0 40 40 0.3 0 < im0p0s0l0000.r8 > ../iimages/iim0p0s0l0000.r8\n");
  fprintf(stderr,"~/bin/iinterp.rh39 1 1 1 1  12 256 128 32  1 0 0  1 5  32 32 32 32   500 3250  -300 300 -300 300 -950 950   1.1 1000 0 0.3  9 1  500 3250 10   0 1 < dumpalltimes.1.txt > observerdumpalltimes.1.true.txt\n");

  fprintf(stderr,"\nin SM: (see blandford.m and zakamska.m)\n");
  fprintf(stderr,"$program $DATATYPE $interptype $READHEADERDATA $WRITEHEADERDATA \\ "
	  "1 $nx $ny $nz $refinement 0 0  $oldgrid $igrid \\ \n"
	  "1 $iinx $iiny $iinz 0 0 $iixmin $iixmax $iiymin $iiymax $iizmin $iizmax \\ \n"
	  "$iRin $iRout $iR0 $ihslope $idefcoord $dofull2pi $EXTRAPOLATE $DEFAULTVALUETYPE < $1 > $2 \n");
  fprintf(stderr,"$program $1 $interptype $READHEADERDATA $WRITEHEADERDATA \\ \n"
	  "1 $nx $ny $nz $refinement 0 0  $oldgrid $igrid \\ \n"
	  "1 $iinx $iiny $iinz 0 0 $iixmin $iixmax $iiymin $iiymax $iizmin $iizmax \\ \n"
	  "$iRin $iRout $iR0 $ihslope $idefcoord $dofull2pi $EXTRAPOLATE $DEFAULTVALUETYPE dumps/gdump < $2 > $3 \n");
  
  exit(0) ;

}



// Parse command line options as a set of switches (just template so far)
void parse_commandline(int argc, char *argv[])
{
  int i;
  int checktdata=0;
  int usage=0;
  int goodarg=0;
  int argcheckstart=1;
  int setbinaryinput=0;


  if(argc>1){
    i=1;
    if (strcmp(argv[i],"-oldargs")==0) {
      // old type of parsing of command line
      old_parse_commandline(argc-1,argv+1); // add 1 to ignore oldargs as if doing old parsing method without any interior changes
      oldparse=1;
    }
  }



  // new type of parsing of command line
  if(oldparse==0){


    // check if usage should be printed, which means 
    if(argc==1){
      usage=1;
    }
    if(argc>1){
      for(i=1;i<argc;i++){
	if(strncmp(argv[i],"-usage",6)==0){ // usage check must come first so only 1 arg of -usage or even argc==1 will allow for usage print-out
	  usage=1;
	}
      }
    }

    if(usage){
      fprintf(stderr,"Usage (argc=%d): \n\n",argc);

      fprintf(stderr,"iinterp <args> < inputfile2stdin > outputfile2stdout\n");

      fprintf(stderr,"E.g. Normal Data Interpolation:\n iinterp -dtype 1 -itype 1 -head 1 1 -oN 1 128 64 128 -grids 1 0 -nN 1 128 128 1 -box 0 0 -40 40 -40 40 0 0 -dofull2pi 1 -defaultvaluetype 0 < infile > outfile \n\n");
      fprintf(stderr,"E.g. Full Diag:\n iinterp -binaryinput 1 -dtype 1001 -head 1 0 -dofull2pi 1 -gdump gdump < infile > outfile \n\n");

      fprintf(stderr,"<args> is any of the below:\n\n");
      fprintf(stderr,"-usage (or no args): Gives Usage Information and exit (don't actually do anything with rest of args)\n");
      if(argc==1) argcheckstart=0; // forces at least one instances in loop below.  "usage ||" will check against argv[0], but doesn't matter since always triggers conditional as true
    }


    // loop over all arguments, looking for space-separated switches.  For each found switch, grab next entry as the parameter
    for(i=argcheckstart;i<argc;i++){
      // assume won't find good argument by default
      goodarg=0;

      // don't use if else if -- just use if's since want "usage" to get triggered (i.e. common triggers and not always specific unique triggers)


      // checking for -usage to trigger goodarg==1
      if(argc>1 && strncmp(argv[i],"-usage",6)==0){
	goodarg++;
      }

      //////////////////
      //
      // normal args
      //
      // could also use atoi(), atof(), and direct assignment of string pointer.
      // can also use strcmp() in a sublevel (e.g. if (strncmp(argv[i+1],"gen",3)==0))
      //
      /////////////////
      if(usage || strncmp(argv[i],"-version",8)==0 || strncmp(argv[i],"-ver",4)==0 || strncmp(argv[i],"--version",9)==0 || strncmp(argv[i],"--ver",5)==0 || strncmp(argv[i],"-v",2)==0 || strncmp(argv[i],"--v",3)==0 ){
	if(usage==0){
	  goodarg++;
	  fprintf(stderr,"iinterp %s\n\n"
		 "Jon's interpolation code for HARM post-processing.\n",
		 VERSION);
	  fprintf(stderr,"         Home page: http://harm.unfuddle.com/projects/3 \n"
		 "Copyright (C) 2002-2010 Jonathan McKinney, et al.  This program is free\n"
		 "software; you can redistribute it and/or modify it under the terms of\n"
		 "the GNU General Public License, with some exceptions.\n\n"
		 "This program is distributed in the hope that it will be useful,\n"
		 "but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
		 "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n"
		 "For more details, see the file license.txt in code tree.\n");
	  exit(0);
	}
	else{
	  fprintf(stderr,"-version -ver --version --ver -v --v : Gives Version Information\n");
	}
      }  
      if(usage || strncmp(argv[i],"-compile",8)==0) {
	if(usage==0){
	  goodarg++;
	  fprintf(stderr,"To compile as part of HARM package, do: \n\n"
		 "make superclean ; make prepiinterp ; make iinterp \n");
	  exit(0);
	}
	else{
	  fprintf(stderr,"-compile : Gives Compilation Information\n");
	}
      }  
      if (usage || strcmp(argv[i],"-binaryinput")==0) {
	if(usage==0){
	  goodarg++;
	  if(i+1<argc) sscanf(argv[++i],"%d",&binaryinput) ;
	  setbinaryinput=1;
	}
	else{
	  fprintf(stderr,"-binaryinput <binaryinput>\n");
	  fprintf(stderr,"\t<binaryinput>: 0=text 1=binary (assumed little Endian or at least same Endian)\n");
	}
      }
      if (usage || strcmp(argv[i],"-binaryoutput")==0) {
	if(usage==0){
	  goodarg++;
	  if(i+1<argc) sscanf(argv[++i],"%d",&binaryoutput) ;
	}
	else{
	  fprintf(stderr,"-binaryoutput <binaryoutput>\n");
	  fprintf(stderr,"\t<binaryoutput>: 0=text 1=binary\n");
	}
      }
      if (usage || strcmp(argv[i],"-inFTYPE")==0) {
	if(usage==0){
	  goodarg++;
	  if(i+1<argc){
	    i++;
	    if(strncmp(argv[i],"byte",4)==0){
	      strcpy(inFTYPE,"b");
	    }
	    else if(strncmp(argv[i],"int",3)==0){
	      strcpy(inFTYPE,"i");
	    }
	    else if(strncmp(argv[i],"longint",7)==0){
	      strcpy(inFTYPE,"li");
	    }
	    else if(strncmp(argv[i],"longlongint",11)==0){
	      strcpy(inFTYPE,"lli");
	    }
	    else if(strncmp(argv[i],"float",5)==0){
	      strcpy(inFTYPE,"f");
	    }
	    else if(strncmp(argv[i],"double",6)==0){
	      strcpy(inFTYPE,"d");
	    }
	    else if(strncmp(argv[i],"longdouble",10)==0){
	      strcpy(inFTYPE,"ld");
	    }
	    else{
	      fprintf(stderr,"Unknown inFTYPE\n");
	      exit(1);
	    }
	  }
	}
	else{
	  fprintf(stderr,"-inFTYPE <inFTYPEstring>\n");
	  fprintf(stderr,"\t<inFTYPEstring>: byte, int, longint, longlongint, float, double, longdouble\n");
	}
      }
      if (usage || strcmp(argv[i],"-dtype")==0) {
	if(usage==0){
	  goodarg++;
	  if(i+1<argc) sscanf(argv[++i],"%d",&DATATYPE) ; // 0,1, etc.
	  // set specific defaults, but don't override user choice if already set
	  if(setbinaryinput==0){
	    if(DATATYPE==0) binaryinput=1;
	    else binaryinput=0;
	  }
	}
	else{
	  fprintf(stderr,"-dtype <DATATYPE>\n");
	  fprintf(stderr,"\t<DATATYPE> as one of:\n"
	  "\t0=image (input&output: byte binary default, 1 column only)\n"
	  "\t1=data (input&output: text default, 1=scalar 1 column)\n"
	  "\t2,3,4,5=correspond to output of orthonormal vectors v^0,v^1,v^2,v^3 (inputting all 4 columns of data u^0 u^1 u^2 u^3)\n"
	  "\t6,7,8,9=correspond to output of orthonormal vectors v^0,v^1,v^2,v^3 (inputting all 4 columns of data u_0 u_1 u_2 u_3)\n"
	  "\t11=corresponds to output of \\detg T^x1_t[EM]/sin(\\theta) (inputting all 7 columns of data: u^t v^1 v^2 v^3 B^1 B^2 B^3)\n"
	  "\t12=output lower component (inputting all 4 columns of data: u^i)\n"
	  "\t100+x=corresponds to inputting x-number of 4-vectors and outputting all 4-vectors in orthonormal basis without any interpolation\n"
	  "\t1001=Full Diag\n"
	  );
	}
      }
      if (usage || strcmp(argv[i],"-itype")==0) {
	if(usage==0){
	  goodarg++;
	  if(i+1<argc) sscanf(argv[++i],"%d",&INTERPTYPE) ; // 0,1,2,3
	}
	else{
	  fprintf(stderr,"-itype <INTERPTYPE>\n");
	  fprintf(stderr,"\t<INTERPTYPE>: 0=nearest 1=bi-linear 2=planar 3=bicubic\n");
	}
      }
      if (usage || strcmp(argv[i],"-head")==0) {
	if(usage==0){
	  goodarg++;
	  if(i+1<argc) sscanf(argv[++i],"%d",&READHEADER) ; // 0 or 1
	  if(i+1<argc) sscanf(argv[++i],"%d",&WRITEHEADER) ; // 0 or 1
	}
	else{
	  fprintf(stderr,"-head <DOREADHEADER> <DOWRITEHEADER>\n");
	  fprintf(stderr,"\t<DOREADHEADER> or <DOWRITEHEADER>: 0 or 1 for each\n");
	}
      }
      if (usage || strcmp(argv[i],"-oN")==0) {
	if(usage==0){
	  goodarg++;
	  if(i+1<argc) sscanf(argv[++i],"%d",&oN0) ;
	  if(i+1<argc) sscanf(argv[++i],"%d",&oN1) ;
	  if(i+1<argc) sscanf(argv[++i],"%d",&oN2) ;
	  if(i+1<argc) sscanf(argv[++i],"%d",&oN3) ;
	}
	else{
	  fprintf(stderr,"-oN <oN0> <oN1> <oN2> <oN3>\n");
	  fprintf(stderr,"\t<oN?> : old grid sizes for (e.g.) t,r,\\theta,\\phi \n");
	}
      }
      if (usage || strcmp(argv[i],"-refine")==0) {
	if(usage==0){
	  goodarg++;
	  if(i+1<argc) sscanf(argv[++i],SCANARG,&refinefactor) ;// 1.0 then no refinement, just normal old image used
	}
	else{
	  fprintf(stderr,"-refine <refinementfactor>\n");
	  fprintf(stderr,"\t<refinementfactor>=1.0=no refinement, otherwise refines image before interpolation with this factor increase in size: standard is bicubic refinement\n");
	}
      }
      if (usage || strcmp(argv[i],"-filter")==0) {
	if(usage==0){
	  goodarg++;
	  if(i+1<argc) sscanf(argv[++i],"%d",&filter) ;// 0=no filter #=filter given image within surrounding # pixels per pixel with sigma
	  if(filter!=0 && i+1<argc) sscanf(argv[++i],SCANARG,&sigma) ;// only used if filter!=0, then sigma of gaussian filter, usually ~ filter value
	}
	else{
	  fprintf(stderr,"-filter <filter> <sigma>\n");
	  fprintf(stderr,"\t<filter>: 0=no filter #=filter image with surrounding # pixels per pixel with sigma width\n");
	  fprintf(stderr,"\t<sigma>: only used if filter!=0, then sigma of gaussian filter, usually ~ filter value\n");
	}
      }
      if (usage || strcmp(argv[i],"-grids")==0) {
	if(usage==0){
	  goodarg++;
	  if(i+1<argc) sscanf(argv[++i],"%d",&oldgridtype) ; // 0, 1, 2, 3, 4, and 5 currently
	  if(i+1<argc) sscanf(argv[++i],"%d",&newgridtype) ; // -1, and above 0+ versions
	}
	else{
	  fprintf(stderr,"-grids <oldgridtype> <newgridtype>\n");
	  fprintf(stderr,"\t<oldgridtype>: (V in GRMHD code): 0=Cartesian  1=spherical polar 2=cylindrical 3=log(z) vs. log(R), 4=x'=sin\\theta log(r) z'=cos\\theta log(r) 5=Cartesian w/ light travel time accounting\n");
	  fprintf(stderr,"\t<newgridtype>: (output coord system): -1=no change (and rest same as above)\n");
	  fprintf(stderr,"\tNote: Assume if oN3>1 and oldgridtype==1, then periodic in \\phi since otherwise extrapolate values and poor at low resolution.\n");
	}
      }
      if (usage || strcmp(argv[i],"-nN")==0) {
	if(usage==0){
	  goodarg++;
	  if(i+1<argc) sscanf(argv[++i],"%d",&nN0) ;
	  if(i+1<argc) sscanf(argv[++i],"%d",&nN1) ;
	  if(i+1<argc) sscanf(argv[++i],"%d",&nN2) ;
	  if(i+1<argc) sscanf(argv[++i],"%d",&nN3) ;
	}
	else{
	  fprintf(stderr,"-nN <nN0> <nN1> <nN2> <nN3>\n");
	  fprintf(stderr,"\t<nN?> : new grid sizes for tc,xc,yc,zc (e.g. for Cart as for tprime,xprime,zprime,yprime as potentially rotated by tnrdegrees)  \n");
	}
      }
      if (usage || strcmp(argv[i],"-box")==0) {
	if(usage==0){
	  goodarg++;
	  //(t=time, x=R-cyl, y=Z-cyl, z=Y-cyl)
	  if(i+1<argc) sscanf(argv[++i],SCANARG,&starttc) ; // arbitrary
	  if(i+1<argc) sscanf(argv[++i],SCANARG,&endtc) ; // arbitrary
	  if(i+1<argc) sscanf(argv[++i],SCANARG,&startxc) ; // arbitrary
	  if(i+1<argc) sscanf(argv[++i],SCANARG,&endxc) ; // arbitrary
	  if(i+1<argc) sscanf(argv[++i],SCANARG,&startyc) ; // arbitrary
	  if(i+1<argc) sscanf(argv[++i],SCANARG,&endyc) ; // arbitrary
	  if(i+1<argc) sscanf(argv[++i],SCANARG,&startzc) ; // arbitrary
	  if(i+1<argc) sscanf(argv[++i],SCANARG,&endzc) ; // arbitrary
	}
	else{
	  fprintf(stderr,"-box <starttc> <endtc> <startxc> <endxc> <startyc> <endyc> <startzc> <endzc>\n");
	  fprintf(stderr,"\t<starttc>: inner interp t(t-Cart,t-Cyl t=time)\n");
	  fprintf(stderr,"\t<endtc>: outer t\n");
	  fprintf(stderr,"\t<startxc>: inner interp x(x-Cart,R-Cyl)\n");
	  fprintf(stderr,"\t<endxc>: outer x\n");
	  fprintf(stderr,"\t<startyc>: inner interp y(z-Cart,z-Cyl)\n");
	  fprintf(stderr,"\t<endyc>: outer y\n");
	  fprintf(stderr,"\t<startzc>: inner interp z(y-Cart,y-Cyl)\n");
	  fprintf(stderr,"\t<endzc>: outer z\n");
	}
      }
      if (usage || strcmp(argv[i],"-coord")==0) {
	if(usage==0){
	  goodarg++;
	  // often other coord.c dependent stuff needed
	  if(i+1<argc) sscanf(argv[++i],SCANARG,&Rin) ; // could use setRin()
	  if(i+1<argc) sscanf(argv[++i],SCANARG,&Rout) ;
	  if(i+1<argc) sscanf(argv[++i],SCANARG,&R0) ;
	  if(i+1<argc) sscanf(argv[++i],SCANARG,&hslope) ;
	}
	else{
	  fprintf(stderr,"-coord <Rin> <Rout> <R0> <hslope>\n");
	  // some basic grid parameters, but sometimes need specific coord.c file with its parameters
	  fprintf(stderr,"\t<Rin>: Inner radial edge\n");
	  fprintf(stderr,"\t<Rout>: Outer radial edge\n");
	  fprintf(stderr,"\t<R0>: Radial inner-grid enhancement factor\n");
	  fprintf(stderr,"\t<hslope>: theta grid refinement factor\n");
	  fprintf(stderr,"\tNote: -coord required if not reading header\n");
	}
      }
      if (usage || strcmp(argv[i],"-defcoord")==0) {
	if(usage==0){
	  goodarg++;
	  if(i+1<argc) sscanf(argv[++i],"%d",&defcoord) ;
	}
	else{
	  fprintf(stderr,"-defcoord <defcoord>\n");
	  fprintf(stderr,"\t<defcoord>: which coordinate system (see coord.c)\n");
	  fprintf(stderr,"\tNote: -defcoord required if not reading header\n");
	}
      }
      if (usage || strcmp(argv[i],"-dofull2pi")==0) {
	if(usage==0){
	  goodarg++;
	  if(i+1<argc) sscanf(argv[++i],"%d",&dofull2pi) ;
	}
	else{
	  fprintf(stderr,"-dofull2pi <dofull2pi>\n");
	  fprintf(stderr,"\t<dofull2pi>: whether to do full 2pi or not (see coord.c)\n");
	}
      }
      if (usage || strcmp(argv[i],"-tdata")==0) {
	if(usage==0){
	  goodarg++;
	  // below 2 for 4D data inputs, specifying start and end times for dumps used
	  if(i+2<argc) checktdata=1;
	  if(i+1<argc) sscanf(argv[++i],SCANARG,&starttdata) ;
	  if(i+1<argc) sscanf(argv[++i],SCANARG,&endtdata) ;
	}
	else{
	  fprintf(stderr,"-tdata <starttdata> <endtdata>\n");
	  // below are only inputted when oN0>1 || nN0>1
	  fprintf(stderr,"\t<starttdata>: time 4D input dumps start (assume uniformly spaced in time, and corresponds to actual time when data exists, not FACE values, but CENT in terms of internal interpolation routines)\n");
	  fprintf(stderr,"\t<endtdata>: time 4D input dumps end\n");
	  fprintf(stderr,"\tNote: -tdata required if oN0>1 || nN1>1\n");
	}
      }
      if (usage || strcmp(argv[i],"-tnrdeg")==0) {
	if(usage==0){
	  goodarg++;
	  // angle for CARTLIGHT or CART grid
	  if(i+1<argc) sscanf(argv[++i],SCANARG,&tnrdegrees) ;
	}
	else{
	  fprintf(stderr,"-tnrdeg <tnrdegrees>\n");
	  fprintf(stderr,"\t<tnrdegrees>: angle [in degrees] between observer and z-axis of original grid.  0 degrees gives no transformation.  20degrees means data rotated 20degrees from z-axis towards x-axis around y-axis.\n");
	  fprintf(stderr,"\tNote: -tnrdeg recommended if want rotation or using multiple time data (i.e. -tdata with oN0>1 or nN0>1)\n");
	}
      }
      if (usage || strcmp(argv[i],"-extrap")==0) {
	if(usage==0){
	  goodarg++;
	  sscanf(argv[++i],"%d",&EXTRAPOLATE);
	}
	else{
	  fprintf(stderr,"-extrap <extrapolate>\n");
	  // below 1 are separately optional
	  fprintf(stderr,"\t<extrapolate>: 0 = no, 1 = yes\n");
	}
      }
      if (usage || strcmp(argv[i],"-defaultvaluetype")==0) {
	if(usage==0){
	  goodarg++;
	  if(i+1<argc) sscanf(argv[++i],"%d",&defaultvaluetype);
	}
	else{
	  fprintf(stderr,"-defaultvaluetype <defaultvaluetype>\n");
	  fprintf(stderr,"\t<defaultvaluetype>: 0 = min if scalar 0 if vector, 1 = min, 2 = max, 3 = 0.0, 4 = 1E35 for v5d missingdata\n");
	}
      }
      if (usage || strcmp(argv[i],"-gdump")==0) {
	if(usage==0){
	  goodarg++;
	  if(i+1<argc) getgdump=1;
	  if(i+1<argc) sscanf(argv[++i],"%s",&gdumpfilename[0]);
	}
	else{
	  fprintf(stderr,"-gdump <gdumpfilepathname>\n");
	  // below is optional but requires above 2 to be read-in
	  fprintf(stderr,"\t<gdumpfilepathname> : only if vector type (DATATYPE=(e.g.) 2,3,4,5,11,12,100...\n");
	}
      }
      if (usage || strcmp(argv[i],"-verbose")==0) {
	if(usage==0){
	  goodarg++;
	  if(i+1<argc) sscanf(argv[++i],"%d",&VERBOSITY) ;
	}
	else{
	  fprintf(stderr,"-verbose <VERBOSITY>\n");
	  fprintf(stderr,"\t<VERBOSITY> : 0 = no extra info : 1 = some extra info : 2 = lots of extra info\n");
	}
      }
      if (usage || strcmp(argv[i],"-debug")==0){
	if(usage==0){
	  goodarg++;
	  if(i+1<argc) sscanf(argv[++i],"%d",&DEBUGINTERP) ;
	}
	else{
	  fprintf(stderr,"-debug <DEBUGINTERP>\n");
	  fprintf(stderr,"\t<DEBUGINTERP> : Level of debug for interpolation: 0 or 1 currently\n");
	}
      }
      if (usage || strcmp(argv[i],"-simpledebug")==0){
	if(usage==0){
	  goodarg++;
	  if(i+1<argc) sscanf(argv[++i],"%d",&SIMPLEDEBUGINTERP) ;
	}
	else{
	  fprintf(stderr,"-simpledebug <SIMPLEDEBUGINTERP>\n");
	  fprintf(stderr,"\t<SIMPLEDEBUGINTERP> : Level of simpledebug for interpolation: 0 or 1 currently\n");
	}
      }
      if (usage || argv[i][0]!='-' ){
	if(usage==0){
	  goodarg++;
	  //      v5dfile[gopointer] = argv[i];
	  //      gopointer++;
	}
	else{
	  fprintf(stderr,"\n\nWithout - in front, can code-up special options\n\n");
	}
      }


      if(argc>1 && goodarg==0){
	fprintf(stderr,"\n\nUnknown option: %s\n\n", argv[i] );
      }
	
    }// loop over command-line arguments
      



  }// end if new type of parsing of command line



  // exit if giving usage
  if(usage){
    exit(0);
  }


  ///////////////////////////
  //
  // If reach here, then actually got arguments instead of just showing usage
  //
  ///////////////////////////

  // some checks if actually got args
  if(nN0>1 || oN0>1 && checktdata==0){
    fprintf(stderr,"Have nN0>1 || oN0>1 but didn't set tdata, using default!\n");
  }


}




// Parse command line options as a string of input variables
void old_parse_commandline(int argc, char *argv[])
{
  int i;
  int basicargcnum=32+1;



  // see if enough arguments to do anything
  if(argc < basicargcnum){
    for(i=0;i<argc;i++){
      fprintf(stderr,"argv[%d]=%s\n",i,argv[i]);
    }
    old_usage(argc,basicargcnum); // report how to use program if not enough arguments
  }



  // get args
  i=1;
  sscanf(argv[i++],"%d",&DATATYPE) ; // 0,1, etc.
  if(DATATYPE==0) binaryinput=1; // assume image is binary
  else binaryinput=0; // assume text

  sscanf(argv[i++],"%d",&INTERPTYPE) ; // 0,1,2,3

  sscanf(argv[i++],"%d",&READHEADER) ; // 0 or 1
  sscanf(argv[i++],"%d",&WRITEHEADER) ; // 0 or 1

  sscanf(argv[i++],"%d",&oN0) ;
  sscanf(argv[i++],"%d",&oN1) ;
  sscanf(argv[i++],"%d",&oN2) ;
  sscanf(argv[i++],"%d",&oN3) ;


  sscanf(argv[i++],SCANARG,&refinefactor) ;// 1.0 then no refinement, just normal old image used
  sscanf(argv[i++],"%d",&filter) ;// 0=no filter #=filter given image within surrounding # pixels per pixel with sigma
  sscanf(argv[i++],SCANARG,&sigma) ;// only used if filter!=0, then sigma of gaussian filter, usually ~ filter value

  sscanf(argv[i++],"%d",&oldgridtype) ; // 0, 1, 2, 3, 4, and 5 currently
  sscanf(argv[i++],"%d",&newgridtype) ; // -1, and above 0+ versions

  sscanf(argv[i++],"%d",&nN0) ; // arbitrary
  sscanf(argv[i++],"%d",&nN1) ; // arbitrary
  sscanf(argv[i++],"%d",&nN2) ; // arbitrary
  sscanf(argv[i++],"%d",&nN3) ; // arbitrary

  //(t=time, x=R-cyl, y=Z-cyl, z=Y-cyl)
  sscanf(argv[i++],SCANARG,&starttc) ; // arbitrary
  sscanf(argv[i++],SCANARG,&endtc) ; // arbitrary
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
  sscanf(argv[i++],"%d",&dofull2pi) ;



  if(nN0>1 || oN0>1){
    // extend basicargcnum if needed
    basicargcnum+=2;
    // below 2 for 4D data inputs, specifying start and end times for dumps used
    sscanf(argv[i++],SCANARG,&starttdata) ;
    sscanf(argv[i++],SCANARG,&endtdata) ;
  }


  // angle for CARTLIGHT or CART grid
  sscanf(argv[i++],SCANARG,&tnrdegrees) ;



  // define (at least) outputvartype
  interpret_commandlineresults_subpart1();


  // conditionally read-in things (all in or out)
  getgdump=0;
  if(argc>basicargcnum){
    
    if(argc<basicargcnum+2){
      fprintf(stderr,"argc=%d insufficient-pre\n",argc);
      exit(1);
    }
    sscanf(argv[i++],"%d",&EXTRAPOLATE);
    sscanf(argv[i++],"%d",&defaultvaluetype);

    // see if doing non-scalar output that needs metric, etc.
    if(outputvartype>0){
      if(argc!=basicargcnum+2+1){
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

}



// Parse command line options as a string of input variables
void interpret_commandlineresults_subpart1(void)
{


  // set number of columns of data
  if(DATATYPE>1){
    immediateoutput=0;// default is not immediate in-out

    // assume always want to transform vectors correctly
    if(DATATYPE>=2 && DATATYPE<=5){
      num4vectors=1; // default
      vectorcomponent=DATATYPE-2;
      outputvartype=1; // simple vector conversion to orthonormal basis
      fprintf(stderr,"Processing vector component %d from u^\\mu\n",vectorcomponent);
    }
    else if(DATATYPE>=6 && DATATYPE<=9){
      num4vectors=1; // default
      vectorcomponent=DATATYPE-6;
      outputvartype=2; // simple vector conversion to orthonormal basis
      fprintf(stderr,"Processing vector component %d from u_\\mu\n",vectorcomponent);
    }
    else if(DATATYPE==11){
      num4vectors=1; // default
      vectorcomponent=1;
      outputvartype=11; // Compute EM Poynting term
      fprintf(stderr,"Computing poloidal EM Polynting angular flux density\n");
    }
    else if(DATATYPE==12){
      num4vectors=1; // default
      vectorcomponent=3; // B_\\phi
      outputvartype=12; // Compute B_\\phi [conserved current]
      fprintf(stderr,"Computing B_\\phi\n");
    }
    else if(DATATYPE>=101 && DATATYPE<1000){
      num4vectors=DATATYPE-100;
      fprintf(stderr,"Transforming %d contravariant 4-vectors to orthonormal basis\n",num4vectors);
      // force input and output grid types to be the same
      outputvartype=1; // conversion to orthonormal basis
      immediateoutput=1;
      vectorcomponent=-1; // indicates to output all components
      // only valid for DATATYPE=1 type of data (i.e. not images)
    }
    else if(DATATYPE>=1001 && DATATYPE<10000){
      fprintf(stderr,"Full Diagnostics for this data file (all except time derivative stuff)\n");
      // force input and output grid types to be the same
      outputvartype=1001; // Full Diag
      immediateoutput=0; // No interpolation, just process and output, but immediateoutput==0 since not processing per grid point but going over many grid points
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
    fprintf(stderr,"Processing scalar\n");
    immediateoutput=0;// default is not immediate in-out
    outputvartype=0; // scalar
    vectorcomponent=0; // scalar
  }
}



// Parse command line options as a string of input variables
void interpret_commandlineresults(void)
{

  // subpart needed for old parsing method checks of number of args
  interpret_commandlineresults_subpart1();

  ///////////////
  //
  // consistency checks
  //
  ///////////////

  // see if doing non-scalar output that needs metric, etc.
  if(outputvartype>0 && getgdump==0){
    fprintf(stderr,"Need gdump file with outputvartype>0\n");
    exit(1);
  }

  




  // set totalbc in case required
  if(INTERPTYPE==0) totalbc=0;
  else if(INTERPTYPE==1 || INTERPTYPE==2) totalbc=1;
  else if(INTERPTYPE==3) totalbc=2;
  else totalbc=MAXBC; // maximum


  // set sizes
  totalsize[0]=oN0;
  totalsize[1]=oN1;
  totalsize[2]=oN2;
  totalsize[3]=oN3;
  totalzones=totalsize[0]*totalsize[1]*totalsize[2]*totalsize[3];


  // set number of bc's per dimension
  numbc[0]=totalbc*(oN0>1);
  numbc[1]=totalbc*(oN1>1);
  numbc[2]=totalbc*(oN2>1);
  numbc[3]=totalbc*(oN3>1);






  ///////////
  //
  // set things not set by header or user that will be used (e.g. by set_points(), which sets dx[0])
  //
  ///////////


  // ensure new coordinate ignorable but no divisions by zero
  // if too narrow (or vanishing) grid width, then leave a bit resolved so no divisions by zero

  if(nN0==1){
    dtc=(endtc-starttc)/((FTYPE)nN0);
    if(fabs(dtc)<SMALL){
      dtc=NUMEPSILONPOW23;
      starttc-=dtc*0.5;
      endtc+=dtc*0.5;
    }    
  }

  if(nN1==1){
    dxc=(endxc-startxc)/((FTYPE)nN1);
    if(fabs(dxc)<SMALL){
      dxc=NUMEPSILONPOW23;
      startxc-=dxc*0.5;
      endxc+=dxc*0.5;
    }    
  }

  if(nN2==1){
    dyc=(endyc-startyc)/((FTYPE)nN2);
    if(fabs(dyc)<SMALL){
      dyc=NUMEPSILONPOW23;
      startyc-=dyc*0.5;
      endyc+=dyc*0.5;
    }    
  }

  if(nN3==1){
    dzc=(endzc-startzc)/((FTYPE)nN3);
    if(fabs(dzc)<SMALL){
      dzc=NUMEPSILONPOW23;
      startzc-=dzc*0.5;
      endzc+=dzc*0.5;
    }    
  }


  if(oN0>1){
    // dt crucial when using multiple time data.
    // Since set_points doesn't explicitly set startx[0] or dx[0] based upon any problem size information, need to set by hand.
    // set_points() will set startx[0]=0 and dx[0]=dt.
    // when oN0==1, this just keeps X[0] sane, but otherwise it's not used
    // assumes data is uniformly spaced in time (i.e. dumps were evenly spaced when merged)
    // then set_points() will set startx[0] but dx[0]=dt below, such that true lab time: tlab = starttc + dt*(hold+0.5); for hold=0..oN0-1
    // so note that if want exact correctness in time, need input starttdata and endtdata to be chosen so that data appear on CENT for given FACE versions of starttdata and endtdata
    dt=(endtdata-starttdata)/((FTYPE)oN0);
    fprintf(stderr,"Times given assumed user inputs where data sits (i.e. CENT in this code): starttdata=%21.15g endtdata=%21.15g dt=%21.15g\n",starttdata,endtdata,dt);
    // assume user inputs dump times corresopnding to location of data in time, while iinterp wants these to be FACE values where data located at CENT
    starttdata-=dt*0.5;// back up to FACE value
    endtdata+=dt*0.5;// bump up to FACE value
    dt=(endtdata-starttdata)/((FTYPE)oN0); // FACE type calculation
    fprintf(stderr,"Adjusted Times so inputted times (assumed colocated with time data) are actually FACE values that give a box with CENT values for the dump times: starttdata=%21.15g endtdata=%21.15g dt=%21.15g\n",starttdata,endtdata,dt);
    fprintf(stderr,"SUPERNOTE: User will still have to specify new-grid box size in time based upon FACE positions (i.e. not times wanted, but time box with values located at CENT)\n");
    //  dX[0]=dt; (will be set in set_points()
    // So X[0] = starttdata + (hold + startpos[0] + 0.5) dx[0]; and dx[0]=dX[0]=dt = (endtdata-starttdata)/oN0;
  }
  else{
    // need to make time coordinate ignorable in case user set the start and end times to be the same
    // avoids divisions by zero
    dt=1.0;
    starttdata=starttc;
    endtdata=endtc;
  }





  if(filter && (oN3!=1 || oN0!=1)){
    filter=0;
    fprintf(stderr,"Turned off filter since oN3=%d or  oN0=%d\n",oN3,oN0); fflush(stderr);
  }


  if(fabs(refinefactor-1.0)>0.1 && (oN3!=1 || oN0!=1) ){
    refinefactor=1.0;
    fprintf(stderr,"Turned off refinement since oN3=%d or oN0=%d\n",oN3,oN0); fflush(stderr);
  }

  if( (oN3>1 || oN0>1) && (INTERPTYPE==3 || INTERPTYPE==2) ){
    fprintf(stderr,"PLANAR and BICUBIC interpolation not setup for oN3=%d>1 or oN0=%d>1 -- uses nearest for k,h",oN3,oN0);
  }


  if(READHEADER)  jonheader=1;
  else jonheader=0;



  // process case where not interpolating and just processing input and output per cell at once
  if(immediateoutput){
    // force this even if user messed up
    newgridtype=GRIDTYPENOCHANGE;
    // force output grid size to be same as input
    nN0=oN0;
    nN1=oN1;
    nN2=oN2;
    nN3=oN3;
  }




}




// set default options
void defaultoptions(void)
{


  DATATYPE=1; binaryinput=0;

  INTERPTYPE=1;
  READHEADER=1;
  WRITEHEADER=1; binaryoutput=0;


  oN0=1;
  oN1=64;
  oN2=64;
  oN3=64;

  refinefactor=1.0;
  filter=0;
  sigma=1.0;

  oldgridtype=1;
  newgridtype=0;
  
  nN0=1;
  nN1=64;
  nN2=64;
  nN3=64;
  
  starttc=0;
  endtc=2000;
  startxc=-40;
  endxc=40;
  startyc=-40;
  endyc=40;
  startzc=-40;
  endzc=40;

  Rin=1.1;
  Rout=40;
  R0=0;
  hslope=0.3;

  defcoord=0;

  dofull2pi=1;

  starttdata=0;
  endtdata=2000;

  tnrdegrees=0;

  getgdump=0;
  strcpy(&gdumpfilename[0],"NOFILESET");

  EXTRAPOLATE=1; // default to extrapolate
  defaultvaluetype=0; // assume typical default

  immediateoutput=0;// default is not immediate in-out
  outputvartype=0; // scalar
  vectorcomponent=0; // scalar

  DEBUGINTERP=0;
  SIMPLEDEBUGINTERP=0;
  VERBOSITY=1; // default to moderate verbosity
  

}
