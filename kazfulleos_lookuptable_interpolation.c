

// interpolation stuff
// Should be only stuff that actually accesses global EOS arrays during simulation

// for function prototypes associated with function calls
#define INDEXPARAMETERSPROTOTYPES int kazii, int kazjj, int kazkk, int kazll, int kazmm, \
    int kaziio, int kazjjo, int kazkko, int kazllo, int kazmmo,			\
    int kazstartiii, int kazstartjjj, int kazstartkkk, int kazstartlll, int kazstartmmm, \
    int kazendiii, int kazendjjj, int kazendkkk, int kazendlll, int kazendmmm,	\
    FTYPE kazdi[], FTYPE kazdj[], FTYPE kazdk[], FTYPE kazdl[], FTYPE kazdm[]

// for function calls
#define INDEXPARAMETERSTOPASS kaziiwhichd[whichd],kazjjwhichd[whichd],kazkkwhichd[whichd],kazllwhichd[whichd],kazmmwhichd[whichd], \
    kaziiowhichd[whichd],kazjjowhichd[whichd],kazkkowhichd[whichd],kazllowhichd[whichd],kazmmowhichd[whichd], \
    kazstartiiiwhichd[whichd],kazstartjjjwhichd[whichd],kazstartkkkwhichd[whichd],kazstartlllwhichd[whichd],kazstartmmmwhichd[whichd], \
    kazendiiiwhichd[whichd],kazendjjjwhichd[whichd],kazendkkkwhichd[whichd],kazendlllwhichd[whichd],kazendmmmwhichd[whichd], \
    kazdiwhichd[whichd],kazdjwhichd[whichd],kazdkwhichd[whichd],kazdlwhichd[whichd],kazdmwhichd[whichd]

// for function calls inside main function call (e.g. _dumb() inside others)
#define INDEXPARAMETERSTOSUBPASS kazii,kazjj,kazkk,kazll,kazmm, \
    kaziio,kazjjo,kazkko,kazllo,kazmmo, \
    kazstartiii,kazstartjjj,kazstartkkk,kazstartlll,kazstartmmm, \
    kazendiii,kazendjjj,kazendkkk,kazendlll,kazendmmm, \
    kazdi,kazdj,kazdk,kazdl,kazdm




static int get_eos_fromlookup(int repeatedeos, int tabledimen, int degentable, int whichtable, int whichtablesubtype, int *iffun, int whichindep, FTYPE quant1, int *vartypearraylocal, FTYPE *indexarray, FTYPE *answers, int *badlookups);


static int get_eos_fromlookup_nearest_dumb(INDEXPARAMETERSPROTOTYPES, int repeatedeos, int tabledimen, int degentable, int whichtable, int whichtablesubtype, int *iffun, int whichindep, FTYPE quant1, int *vartypearraylocal, FTYPE *indexarray, FTYPE *answers, int *badlookups);
static int get_eos_fromlookup_nearest(INDEXPARAMETERSPROTOTYPES, int repeatedeos, int tabledimen, int degentable, int whichtable, int whichtablesubtype, int *iffun, int whichindep, FTYPE quant1, int *vartypearraylocal, FTYPE *indexarray, FTYPE *answers, int *badlookups);
static int get_eos_fromlookup_linear(INDEXPARAMETERSPROTOTYPES, int repeatedeos, int tabledimen, int degentable, int whichtable, int whichtablesubtype, int *iffun, int whichindep, FTYPE quant1, int *vartypearraylocal, FTYPE *indexarray, FTYPE *answers, int *badlookups);
static int get_eos_fromlookup_parabolic(INDEXPARAMETERSPROTOTYPES, int repeatedeos, int tabledimen, int degentable, int whichtable, int whichtablesubtype, int *iffun, int whichindep, FTYPE quant1, int *vartypearraylocal, FTYPE *indexarray, FTYPE *answers, int *badlookups);
static int get_eos_fromlookup_parabolicfull(INDEXPARAMETERSPROTOTYPES, int repeatedeos, int tabledimen, int degentable, int whichtable, int whichtablesubtype, int *iffun, int whichindep, FTYPE quant1, int *vartypearraylocal, FTYPE *indexarray, FTYPE *answers, int *badlookups);

static void get_arrays_eostable_direct(int whichd, int whichdegen, int whichtable, int whichtablesubtype, int *iffun, int mmm, int lll, int kkk, int jjj, int iii, double *values);
static void get_arrays_eostable_direct_temperature(int whichd, int whichtable, int mmm, int lll, int kkk, int jjj, int iii, double *values);









// General EOS interpolation wrapper
static int get_eos_fromlookup(int repeatedeos, int tabledimen, int degentable, int whichtable, int whichtablesubtype, int *iffun, int whichindep, FTYPE quant1, int *vartypearraylocal, FTYPE *indexarray, FTYPE *answers, int *badlookups)
{
  int failreturn;
  int whichd;

  ///////////
  //
  // definition consistent with numerical assignments of indecies of arrays
  //
  ///////////
  whichd = whichindep-1;

  
#if(0)
  // neartest dumb for normal table
  if(degentable==0){
    failreturn=get_eos_fromlookup_nearest_dumb(INDEXPARAMETERSTOPASS,repeatedeos,tabledimen, degentable, whichtable, whichtablesubtype, iffun, whichindep, quant1, vartypearraylocal, indexarray,answers,badlookups);
  }
  else{
    // linear
    failreturn=get_eos_fromlookup_linear(INDEXPARAMETERSTOPASS,repeatedeos,tabledimen,degentable, whichtable, whichtablesubtype, iffun, whichindep, quant1, vartypearraylocal, indexarray,answers,badlookups);
  }
#elif(0)
  // nearest dumb for all
  failreturn=get_eos_fromlookup_nearest_dumb(INDEXPARAMETERSTOPASS,repeatedeos,tabledimen,degentable, whichtable, whichtablesubtype, iffun, whichindep, quant1, vartypearraylocal, indexarray,answers,badlookups);
#elif(0)
  // nearest
  failreturn=get_eos_fromlookup_nearest(INDEXPARAMETERSTOPASS,repeatedeos,tabledimen,degentable, whichtable, whichtablesubtype, iffun, whichindep, quant1, vartypearraylocal, indexarray,answers,badlookups);

#elif(1)
  // linear
  failreturn=get_eos_fromlookup_linear(INDEXPARAMETERSTOPASS,repeatedeos,tabledimen,degentable, whichtable, whichtablesubtype, iffun, whichindep, quant1, vartypearraylocal, indexarray,answers,badlookups);

#elif(0)
  // parabolic for density and tri-linear otherwise
  failreturn=get_eos_fromlookup_parabolic(INDEXPARAMETERSTOPASS,repeatedeos,tabledimen,degentable, whichtable, whichtablesubtype, iffun, whichindep, quant1, vartypearraylocal, indexarray,answers,badlookups);
#elif(0)
  // parabolic for all quantities
  failreturn=get_eos_fromlookup_parabolicfull(INDEXPARAMETERSTOPASS,repeatedeos,tabledimen,degentable, whichtable, whichtablesubtype, iffun, whichindep, quant1, vartypearraylocal, indexarray,answers,badlookups);
#endif

  // return whether really bad failure
  return(failreturn);

}







// full parabolic interpolation
// Uses globals so can make them thread safe instead of using static's that are not
// Assumes H is not a dependent dimension as for whichdatatype==4
static int get_eos_fromlookup_parabolicfull(INDEXPARAMETERSPROTOTYPES, int repeatedeos, int tabledimen, int degentable, int whichtable, int whichtablesubtype, int *iffun, int whichindep, FTYPE quant1, int *vartypearraylocal, FTYPE *indexarray, FTYPE *answers, int *badlookups)
{
  FTYPE tempcheck;
  FTYPE tempanswers[MAXEOSPIPELINE];
  int shouldloginterp[MAXEOSPIPELINE];

  FTYPE totalf[MAXEOSPIPELINE][3][3][3][3]; // 3 values for parabolic interpolation
  FTYPE (*tfptr)[3][3][3][3];

  int includef[MAXEOSPIPELINE][3][3][3][3]; // 3 values for parabolic interpolation
  int (*includeptr)[3][3][3][3];

  int include2f[MAXEOSPIPELINE][3][3][3]; // 3 values for parabolic interpolation
  int (*include2ptr)[3][3][3];

  int include3f[MAXEOSPIPELINE][3][3]; // 3 values for parabolic interpolation
  int (*include3ptr)[3][3];

  int include4f[MAXEOSPIPELINE][3]; // 3 values for parabolic interpolation
  int (*include4ptr)[3];

  int iii,jjj,kkk,lll,mmm;
  FTYPE ieos,jeos,keos,leos,meos;

  int whichdegenfun;
  FTYPE xmx0,AA,BB;

  int numcols,coli;


  

  ///////////////
  //
  // get number of columns possible for this whichtablesubtype
  //
  ///////////////
  numcols = numcolintablesubtype[whichtablesubtype];




  ///////////////
  //
  // store index within table
  //
  ///////////////
  ieos=indexarray[1];
  jeos=indexarray[2];
  keos=indexarray[3];
  leos=indexarray[4];
  meos=0;


  ///////////////
  //
  // pointer shift
  //
  ///////////////
  tfptr=(FTYPE (*)[3][3][3][3]) (&(totalf[0][1][1][1][1])); // so tfptr[-1,0,1]
  includeptr=(int (*)[3][3][3][3]) (&(includef[0][1][1][1][1])); // so includeptr[-1,0,1]
  include2ptr=(int (*)[3][3][3]) (&(includef[0][1][1][1])); // so include2ptr[-1,0,1]
  include3ptr=(int (*)[3][3]) (&(includef[0][1][1])); // so include3ptr[-1,0,1]
  include4ptr=(int (*)[3]) (&(includef[0][1])); // so include4ptr[-1,0,1]


  ///////////
  //
  // definition consistent with numerical assignments of indecies of arrays
  //
  ///////////
  whichdegenfun = whichindep-1;


  //////////
  //
  // determine if should do log interpolation
  //
  //////////
  get_dologinterp_subtype_wrapper(degentable, whichtablesubtype, numcols, shouldloginterp);
  


  //////////
  //
  // Get loop ranges within grid for lookup
  //
  //////////
  if(repeatedeos==0){
    // assume jeos, keos, leos set correctly for given degentable and tabledimen (i.e. are 0 or integers in those cases so d's computed correctly)
    kazii=ROUND2INT(ieos); //ii=(int)ieos; // round instead when doing parabolic interpolation with 3 points
    // limit to within table
    if(kazii<1) kazii=1; // 0 is minimum
    if(kazii>tablesize[whichtable][vartypearraylocal[1]]-2) kazii=tablesize[whichtable][vartypearraylocal[1]]-2; // N-1 is maximum

    kazjj=ROUND2INT(jeos); //kazjj=(int)jeos; // round instead when doing parabolic interpolation with 3 points
    // limit to within table
    if(degentable==0){
      if(kazjj<1) kazjj=1; // 0 is minimum
      if(kazjj>tablesize[whichtable][whichindep]-2) kazjj=tablesize[whichtable][whichindep]-2; // N-1 is maximum
    }

    kazkk=ROUND2INT(keos); //kazkk=(int)keos; // round instead when doing parabolic interpolation with 3 points
    // limit to within table
    if(kazkk<1) kazkk=1; // 0 is minimum
    if(kazkk>tablesize[whichtable][vartypearraylocal[3]]-2) kazkk=tablesize[whichtable][vartypearraylocal[3]]-2; // N-1 is maximum

    kazll=ROUND2INT(leos); //kazll=(int)leos; // round instead when doing parabolic interpolation with 3 points
    // limit to within table
    if(kazll<1) kazll=1; // 0 is minimum
    if(kazll>tablesize[whichtable][vartypearraylocal[4]]-2) kazll=tablesize[whichtable][vartypearraylocal[4]]-2; // N-1 is maximum

    kazmm=0;

    // if only choosing 1 value, then choose rounded version
    kaziio=ROUND2INT(ieos);
    kazjjo=ROUND2INT(jeos);
    kazkko=ROUND2INT(keos);
    kazllo=ROUND2INT(leos);
    kazmmo=0;
    

    // set range of loops for different table types
    // tabledimen overrules table type (i.e. take section out of fuller table -- assumed to be 0 index)
    // GODMARK: these things couuld be stored as functions of whichtable/degentable/tabledimen
    if(tablesize[whichtable][vartypearraylocal[1]]==1) { kazstartiii=kazendiii=0; } // kazll will be 0
    else if(tablesize[whichtable][vartypearraylocal[1]]==2) { kazstartiii=0; kazendiii=1;} // kazll will be 0
    else if(tablesize[whichtable][vartypearraylocal[1]]>2) { kazstartiii=-1; kazendiii=1;}

    if(degentable==1 || tablesize[whichtable][whichindep]==1) {kazstartjjj=kazendjjj=0;}
    else if(tablesize[whichtable][whichindep]==2) { kazstartjjj=0; kazendjjj=1;}
    else if(tablesize[whichtable][whichindep]>2) { kazstartjjj=-1; kazendjjj=1;}

    if(tablesize[whichtable][vartypearraylocal[3]]==1) {kazstartkkk=kazendkkk=0;}
    else if(tablesize[whichtable][vartypearraylocal[3]]==2){ kazstartkkk=0; kazendkkk=1;}
    else if(tablesize[whichtable][vartypearraylocal[3]]>2){ kazstartkkk=-1; kazendkkk=1;}
     
    if(tablesize[whichtable][vartypearraylocal[4]]==1) {kazstartlll=kazendlll=0;}
    else if(tablesize[whichtable][vartypearraylocal[4]]==2){ kazstartlll=0; kazendlll=1;}
    else if(tablesize[whichtable][vartypearraylocal[4]]>2){ kazstartlll=-1; kazendlll=1;}

  }

  // overrides (must be placed outside if(repeatedeos==0) conditional so always done since kazendlll,kazendmmm, etc. are stored per whichd and not also per whichtablesubtype
  // only extra table is function of Ynu or H if whichdatatype==4
  if(WHICHDATATYPEGENERAL==4 && (whichtablesubtype!=SUBTYPEEXTRA)){kazll=kazmm=kazllo=kazmmo=meos=leos=kazstartlll=kazendlll=0;}




  //////////////////////
  //
  // Set that default return that no bad lookups
  //
  //////////////////////
  for(coli=0;coli<numcols;coli++){
    badlookups[coli]=0;
  }




#define LOOPKAZIJKL for(iii=kazstartiii;iii<=kazendiii;iii++) for(jjj=kazstartjjj;jjj<=kazendjjj;jjj++)for(kkk=kazstartkkk;kkk<=kazendkkk;kkk++)for(lll=kazstartlll;lll<=kazendlll;lll++)


  // Loop over nearby table values and determine bi-linearly interpolated value
  // 4-D means 2^4=16 positions
  // 2-D means 2^2=4 positions
  // get 3 values as function of density
  mmm=0;
  LOOPKAZIJKL{

    
#if(CHECKIFVALIDEOSDATA)
    if(degentable==0){
      // don't use values of table that have no inversion to temperature
      get_arrays_eostable_direct_temperature(whichdegenfun,whichtable,kazmm+mmm,kazll+lll,kazkk+kkk,kazjj+jjj,kazii+iii,&tempcheck);
    }
    if(degentable==1 || tempcheck>invalidtempcode) // Avoid invalid inversions if T>Tbad, but only deal with temperature if degentable==0
#else
    if(1)
#endif
    {

      // get value
      get_arrays_eostable_direct(whichdegenfun,degentable,whichtable,whichtablesubtype,iffun,kazmm+mmm,kazll+lll,kazkk+kkk,kazjj+jjj,kazii+iii,tempanswers);


      // non-pipelined section
      for(coli=0;coli<numcols;coli++){
	if(iffun[coli]){

	  if(shouldloginterp[coli] && tempanswers[coli]<=0.0){
	    // avoid point if should have been log but value is <=0.0
	    includeptr[coli][iii][jjj][kkk][lll]=0;
	  }
	  else{
	    // set to include this point in interpolation
	    includeptr[coli][iii][jjj][kkk][lll]=1;
	    // then include general offset
	    if(degentable==1) offsetquant2_general(whichdegenfun, quant1, tempanswers[coli], &tempanswers[coli]);
	    tfptr[coli][iii][jjj][kkk][lll]=tempanswers[coli];
	  }
	}// end iffun==1
      }//end over coli
    }// end if good temperature or doing degentable
    else{
      // no temperature found

      for(coli=0;coli<numcols;coli++){
	if(iffun[coli]){
	  // set to NOT include point in interpolation
	  includeptr[coli][iii][jjj][kkk][lll]=0;
	}
      }

    }// end else if no temperature lookup

  }// end loop over dimensions


  //////////////////////
  //
  // logify
  //
  //////////////////////
  for(coli=0;coli<numcols;coli++){
    if(iffun[coli]){
      if(shouldloginterp[coli]){
	LOOPKAZIJKL{
	  if(includeptr[coli][iii][jjj][kkk][lll]){
	    tfptr[coli][iii][jjj][kkk][lll] = log10(tfptr[coli][iii][jjj][kkk][lll]);
	  }
	}// end loop over dimensions
      }// end if logifying
    }// end iffun==1
  }// end over coli



#define KAZPARALOOP1 for(jjj=kazstartjjj;jjj<=kazendjjj;jjj++)for(kkk=kazstartkkk;kkk<=kazendkkk;kkk++)for(lll=kazstartlll;lll<=kazendlll;lll++)
#define KAZPARALOOP2 for(kkk=kazstartkkk;kkk<=kazendkkk;kkk++)for(lll=kazstartlll;lll<=kazendlll;lll++)
#define KAZPARALOOP3 for(lll=kazstartlll;lll<=kazendlll;lll++)
#define KAZPARALOOP4 if(1)

  //////////////
  //
  // perform interpolation over values
  //
  //////////////
  for(coli=0;coli<numcols;coli++){
    if(iffun[coli]){
  
      KAZPARALOOP1{


	// now use 3 data points to get density-parabolic distribution and value at ieos
	// have tfptr[coli][iii][jjj][kkk][lll] @ iii=-1,0,1 with i=0 meaning ii and offset being ROUND2INT(ieos)
	// Form parabolic answer
	include2ptr[coli][jjj][kkk][lll]=1; // default
	if(kazstartiii<=-1 && kazendiii>=1 && includeptr[coli][-1][jjj][kkk][lll] && includeptr[coli][0][jjj][kkk][lll]&& includeptr[coli][1][jjj][kkk][lll]){
	  xmx0 = (ieos-(FTYPE)kazii);
	  AA = 0.5*(tfptr[coli][1][jjj][kkk][lll]-tfptr[coli][-1][jjj][kkk][lll]);
	  BB = 0.5*(tfptr[coli][1][jjj][kkk][lll]+tfptr[coli][-1][jjj][kkk][lll]-2.0*tfptr[coli][0][jjj][kkk][lll]);
	  tfptr[coli][0][jjj][kkk][lll] = tfptr[coli][0][jjj][kkk][lll] + AA*xmx0 + BB*xmx0*xmx0;
	  // e.g. xmx0=0 -> tf[0]
	  // e.g. xmx0=-1 -> tf[0]-tf[1]/2+tf[-1]/2 + tf[1]/2+tf[-1]/2-tf[0] = tf[-1]
	  // e.g. xmx0=1 -> tf[0]+tf[1]/2-tf[-1]/2 + tf[1]/2+tf[-1]/2-tf[0] = tf[1]
	}
	else if(kazstartiii<=-1 && kazendiii>=0 && includeptr[coli][-1][jjj][kkk][lll] && includeptr[coli][0][jjj][kkk][lll]){
	  // reduce to linear interpolation
	  xmx0 = (ieos-(FTYPE)(kazii-1));
	  AA=(tfptr[coli][0][jjj][kkk][lll]-tfptr[coli][-1][jjj][kkk][lll]);
	  tfptr[coli][0][jjj][kkk][lll] = tfptr[coli][-1][jjj][kkk][lll] + AA*xmx0;
	  // e.g. xmx0=0 (ieos=kazii-1) -> tf[-1]
	  // e.g. xmx0=1 (ieos=kazii) -> tf[0]
	}
	else if(kazstartiii<=0 && kazendiii>=1 && includeptr[coli][0][jjj][kkk][lll] && includeptr[coli][1][jjj][kkk][lll]){
	  // reduce to linear interpolation
	  xmx0 = (ieos-(FTYPE)(kazii));
	  AA=(tfptr[coli][1][jjj][kkk][lll]-tfptr[coli][0][jjj][kkk][lll]);
	  tfptr[coli][0][jjj][kkk][lll] = tfptr[coli][0][jjj][kkk][lll] + AA*xmx0;
	}
	else if(kazstartiii<=0 && kazendiii>=0 && includeptr[coli][0][jjj][kkk][lll]){
	  // reduce to nearest neighbor
	  tfptr[coli][0][jjj][kkk][lll] = tfptr[coli][0][jjj][kkk][lll];
	}
	else if(kazstartiii<=-1 && kazendiii>=-1 && includeptr[coli][-1][jjj][kkk][lll]){
	  // reduce to nearest neighbor
	  tfptr[coli][0][jjj][kkk][lll] = tfptr[coli][-1][jjj][kkk][lll];
	}
	else if(kazstartiii<=1 && kazendiii>=1 && includeptr[coli][1][jjj][kkk][lll]){
	  // reduce to nearest neighbor
	  tfptr[coli][0][jjj][kkk][lll] = tfptr[coli][1][jjj][kkk][lll];
	}
	else{
	  include2ptr[coli][jjj][kkk][lll]=0;
	}

#if(PRODUCTION==0)
	if(!isfinite(tfptr[coli][0][jjj][kkk][lll])){
	  dualfprintf(fail_file,"1notfinite, :: %d %d %d\n",jjj,kkk,lll);
	}
#endif
    
      } // end over iii,jjj,kkk,lll

      // perform interpolation over values
      KAZPARALOOP2{

	include3ptr[coli][kkk][lll]=1; // default
	if(kazstartjjj<=-1 && kazendjjj>=1 && include2ptr[coli][-1][kkk][lll] && include2ptr[coli][0][kkk][lll] && include2ptr[coli][1][kkk][lll]){
	  xmx0 = (jeos-(FTYPE)kazjj);
	  AA = 0.5*(tfptr[coli][0][1][kkk][lll]-tfptr[coli][0][-1][kkk][lll]);
	  BB = 0.5*(tfptr[coli][0][1][kkk][lll]+tfptr[coli][0][-1][kkk][lll]-2.0*tfptr[coli][0][0][kkk][lll]);
	  tfptr[coli][0][0][kkk][lll] = tfptr[coli][0][0][kkk][lll] + AA*xmx0 + BB*xmx0*xmx0;
	}
	else if(kazstartjjj<=-1 && kazendjjj>=0 && include2ptr[coli][-1][kkk][lll] && include2ptr[coli][0][kkk][lll]){
	  // reduce to linear interpolation
	  xmx0 = (jeos-(FTYPE)(kazjj-1));
	  AA=(tfptr[coli][0][0][kkk][lll]-tfptr[coli][0][-1][kkk][lll]);
	  tfptr[coli][0][0][kkk][lll] = tfptr[coli][0][-1][kkk][lll] + AA*xmx0;
	}
	else if(kazstartjjj<=0 && kazendjjj>=1 && include2ptr[coli][0][kkk][lll] && include2ptr[coli][1][kkk][lll]){
	  // reduce to linear interpolation
	  xmx0 = (jeos-(FTYPE)(kazjj));
	  AA=(tfptr[coli][0][1][kkk][lll]-tfptr[coli][0][0][kkk][lll]);
	  tfptr[coli][0][0][kkk][lll] = tfptr[coli][0][0][kkk][lll] + AA*xmx0;
	}
	else if(kazstartjjj<=0 && kazendjjj>=0 && include2ptr[coli][0][kkk][lll]){
	  // reduce to nearest neighbor
	  tfptr[coli][0][0][kkk][lll] = tfptr[coli][0][0][kkk][lll];
	}
	else if(kazstartjjj<=-1 && kazendjjj>=-1 && include2ptr[coli][-1][kkk][lll]){
	  // reduce to nearest neighbor
	  tfptr[coli][0][0][kkk][lll] = tfptr[coli][0][-1][kkk][lll];
	}
	else if(kazstartjjj<=1 && kazendjjj>=1 && include2ptr[coli][1][kkk][lll]){
	  // reduce to nearest neighbor
	  tfptr[coli][0][0][kkk][lll] = tfptr[coli][0][1][kkk][lll];
	}
	else{
	  include3ptr[coli][kkk][lll]=0;
	}

#if(PRODUCTION==0)
	if(!isfinite(tfptr[coli][0][0][kkk][lll])){
	  dualfprintf(fail_file,"2notfinite, :: %d %d %d\n",0,kkk,lll);
	}
#endif


      } // end over kkk,lll

      // perform interpolation over values
      KAZPARALOOP3{

	include4ptr[coli][lll]=1; // default
	if(kazstartkkk<=-1 && kazendkkk>=1 && include3ptr[coli][-1][lll] && include3ptr[coli][0][lll] && include3ptr[coli][1][lll]){
	  xmx0 = (keos-(FTYPE)kazkk);
	  AA = 0.5*(tfptr[coli][0][0][1][lll]-tfptr[coli][0][0][-1][lll]);
	  BB = 0.5*(tfptr[coli][0][0][1][lll]+tfptr[coli][0][0][-1][lll]-2.0*tfptr[coli][0][0][0][lll]);
	  tfptr[coli][0][0][0][lll] = tfptr[coli][0][0][0][lll] + AA*xmx0 + BB*xmx0*xmx0;
	}
	else if(kazstartkkk<=-1 && kazendkkk>=0 && include3ptr[coli][-1][lll] && include3ptr[coli][0][lll]){
	  // reduce to linear interpolation
	  xmx0 = (keos-(FTYPE)(kazkk-1));
	  AA=(tfptr[coli][0][0][0][lll]-tfptr[coli][0][0][-1][lll]);
	  tfptr[coli][0][0][0][lll] = tfptr[coli][0][0][-1][lll] + AA*xmx0;
	}
	else if(kazstartkkk<=0 && kazendkkk>=1 && include3ptr[coli][0][lll] && include3ptr[coli][1][lll]){
	  // reduce to linear interpolation
	  xmx0 = (keos-(FTYPE)(kazkk));
	  AA=(tfptr[coli][0][0][1][lll]-tfptr[coli][0][0][0][lll]);
	  tfptr[coli][0][0][0][lll] = tfptr[coli][0][0][0][lll] + AA*xmx0;
	}
	else if(kazstartkkk<=0 && kazendkkk>=0 && include3ptr[coli][0][lll]){
	  // reduce to nearest neighbor
	  tfptr[coli][0][0][0][lll] = tfptr[coli][0][0][0][lll];
	}
	else if(kazstartkkk<=-1 && kazendkkk>=-1 && include3ptr[coli][-1][lll]){
	  // reduce to nearest neighbor
	  tfptr[coli][0][0][0][lll] = tfptr[coli][0][0][-1][lll];
	}
	else if(kazstartkkk<=1 && kazendkkk>=1 && include3ptr[coli][1][lll]){
	  // reduce to nearest neighbor
	  tfptr[coli][0][0][0][lll] = tfptr[coli][0][0][1][lll];
	}
	else{
	  include4ptr[coli][lll]=0;
	}

#if(PRODUCTION==0)
	if(!isfinite(tfptr[coli][0][0][0][lll])){
	  dualfprintf(fail_file,"3notfinite, :: %d %d %d\n",0,0,lll);
	}
#endif
    
      } // end over lll

      // final interpolation
      KAZPARALOOP4{
	if(kazstartlll<=-1 && kazendlll>=1 && include4ptr[coli][-1] && include4ptr[coli][0] && include4ptr[coli][1]){
	  xmx0 = (leos-(FTYPE)kazll);
	  AA = 0.5*(tfptr[coli][0][0][0][1]-tfptr[coli][0][0][0][-1]);
	  BB = 0.5*(tfptr[coli][0][0][0][1]+tfptr[coli][0][0][0][-1]-2.0*tfptr[coli][0][0][0][0]);
	  tfptr[coli][0][0][0][0] = tfptr[coli][0][0][0][0] + AA*xmx0 + BB*xmx0*xmx0;
	}
	else if(kazstartlll<=-1 && kazendlll>=0 && include4ptr[coli][-1] && include4ptr[coli][0]){
	  // reduce to linear interpolation
	  xmx0 = (leos-(FTYPE)(kazll-1));
	  AA=(tfptr[coli][0][0][0][0]-tfptr[coli][0][0][0][-1]);
	  tfptr[coli][0][0][0][0] = tfptr[coli][0][0][0][-1] + AA*xmx0;
	}
	else if(kazstartlll<=0 && kazendlll>=1 && include4ptr[coli][0] && include4ptr[coli][1]){
	  // reduce to linear interpolation
	  xmx0 = (leos-(FTYPE)(kazll));
	  AA=(tfptr[coli][0][0][0][1]-tfptr[coli][0][0][0][0]);
	  tfptr[coli][0][0][0][0] = tfptr[coli][0][0][0][0] + AA*xmx0;
	}
	else if(kazstartlll<=0 && kazendlll>=0 && include4ptr[coli][0]){
	  // reduce to nearest neighbor
	  tfptr[coli][0][0][0][0] = tfptr[coli][0][0][0][0];
	}
	else if(kazstartlll<=-1 && kazendlll>=-1 && include4ptr[coli][-1]){
	  // reduce to nearest neighbor
	  tfptr[coli][0][0][0][0] = tfptr[coli][0][0][0][-1];
	}
	else if(kazstartlll<=1 && kazendlll>=1 && include4ptr[coli][1]){
	  // reduce to nearest neighbor
	  tfptr[coli][0][0][0][0] = tfptr[coli][0][0][0][1];
	}
	else{

#if(0)
	  // DEBUG:
	  LOOPKAZIJKL{
	    dualfprintf(fail_file,"BOOB: %d : %d : %d %d %d %d : %d : %21.15g\n",degentable,shouldloginterp,iii,jjj,kkk,lll,includeptr[coli][iii][jjj][kkk][lll],tfptr[coli][iii][jjj][kkk][lll]);
	  }
	  dualfprintf(fail_file,"kaz=%d %d %d %d : %d %d : %d %d : %d %d : %d %d\n",kazii,kazjj,kazkk,kazll,kazstartiii,kazendiii,kazstartjjj,kazendjjj,kazstartkkk,kazendkkk,kazstartlll,kazendlll);
#endif


	  dualfprintf(fail_file,"No valid data points in table despite within table.\n");
	  badlookups[coli]=1;
	  // no valid data points, so return badlookup
	  // GODMARK: Should have caught this with the check if within table function, so could return failure
	  //	  return(1);

	}// end else if no valid points

#if(PRODUCTION==0)
	if(!isfinite(tfptr[coli][0][0][0][0])){
	  dualfprintf(fail_file,"4notfinite, :: %d %d %d\n",0,0,0);
	}
#endif

      }// end KAZPARALOOP4




      ///////////////////////
      //
      // final fully parabolic result
      //
      ///////////////////////
      answers[coli]=tfptr[coli][0][0][0][0];


      ////////////////////////////
      //
      // enforce no new extrema to over overshoots in sharp regions of table
      //
      ////////////////////////////
#if(1)
      FTYPE largest=-BIG;
      FTYPE smallest=BIG;
      LOOPKAZIJKL{
	if(includeptr[coli][iii][jjj][kkk][lll]){
	  if(tfptr[coli][iii][jjj][kkk][lll]>largest) largest=tfptr[coli][iii][jjj][kkk][lll];
	  if(tfptr[coli][iii][jjj][kkk][lll]<smallest) smallest=tfptr[coli][iii][jjj][kkk][lll];
	}
      }
      if(answers[coli]>largest) answers[coli]=largest;
      if(answers[coli]<smallest) answers[coli]=smallest;
#endif


      ///////////////////////////
      //
      // unlogify
      //
      //////////////////////////
      if(shouldloginterp[coli]){
	answers[coli]=pow(10.0,answers[coli]);
      }

      //////////////////////////
      //
      // invert offset
      //
      //////////////////////////
      if(degentable==1) offsetquant2_general_inverse(whichdegenfun, quant1, answers[coli], &answers[coli]);

    }// end iffun==1
  }// end over coli



  //////////////////////////////
  //
  // return that no failure
  //
  //////////////////////////////
  return(0);

}










// tri-linear + parabolic (density) interpolation
static int get_eos_fromlookup_parabolic(INDEXPARAMETERSPROTOTYPES, int repeatedeos, int tabledimen, int degentable, int whichtable, int whichtablesubtype, int *iffun, int whichindep, FTYPE quant1, int *vartypearraylocal, FTYPE *indexarray, FTYPE *answers, int *badlookups)
{
  FTYPE totaldist[MAXEOSPIPELINE][3];
  FTYPE (*tdist)[3];

  FTYPE totalanswers[MAXEOSPIPELINE][3]; // 3 values for parabolic interpolation
  FTYPE (*tfptr)[3];

  FTYPE tempanswers[MAXEOSPIPELINE];
  int shouldloginterp[MAXEOSPIPELINE];

  FTYPE dist;
  FTYPE tempcheck;
  int whichdegenfun;


  int iii,jjj,kkk,lll,mmm;
  FTYPE ieos,jeos,keos,leos,meos;

  FTYPE xmx0,AA,BB;
  int numcols,coli;
  


  // get number of columns possible for this whichtablesubtype
  numcols = numcolintablesubtype[whichtablesubtype];



  ///////////
  //
  // Store grid position within EOS table
  //
  ///////////
  ieos=indexarray[1];
  jeos=indexarray[2];
  keos=indexarray[3];
  leos=indexarray[4];
  meos=0;


  ///////////
  //
  // Pointer shift
  //
  ///////////
  tfptr=(FTYPE (*)[3]) (&(totalanswers[0][1])); // so tfptr[-1,0,1][]
  tdist=(FTYPE (*)[3]) (&(totaldist[0][1]));

  ///////////
  //
  // definition consistent with numerical assignments of indecies of arrays
  //
  ///////////
  whichdegenfun = whichindep-1;




  /////////////
  //
  // determine if should do log interpolation
  //
  ////////////
  get_dologinterp_subtype_wrapper(degentable, whichtablesubtype, numcols, shouldloginterp);





  /////////////
  //
  // Get loop ranges and other lookup parameters
  //
  ////////////
  if(repeatedeos==0){
    // assume jeos, keos, leos set correctly for given degentable and tabledimen (i.e. are 0 or integers in those cases so d's computed correctly)
    kazii=ROUND2INT(ieos); //kazii=(int)ieos; // round instead when doing parabolic interpolation with 3 points
    // limit to within table
    if(kazii<1) kazii=1; // 0 is minimum
    if(kazii>tablesize[whichtable][RHOEOS]-2) kazii=tablesize[whichtable][RHOEOS]-2; // N-1 is maximum

    // normal behavior for other independent variables
    kazjj=(int)jeos;
    kazkk=(int)keos;
    kazll=(int)leos;
    kazmm=0;


    kazdj[1]=jeos-(FTYPE)kazjj;
    kazdj[0]=1.0-kazdj[1];
    kazdk[1]=keos-(FTYPE)kazkk;
    kazdk[0]=1.0-kazdk[1];
    kazdl[1]=leos-(FTYPE)kazll;
    kazdl[0]=1.0-kazdl[1];
    kazdm[1]=meos-(FTYPE)kazmm;
    kazdm[0]=1.0-kazdm[1];
    
    // set range of loops for different table types
    // tabledimen overrules table type (i.e. take section out of fuller table -- assumed to be 0 index)
    // GODMARK: these things couuld be stored as functions of whichtable/degentable/tabledimen
    if(tablesize[whichtable][vartypearraylocal[1]]!=1) kazendiii=1;
    else kazendiii=0;

    if(degentable==0 && tablesize[whichtable][vartypearraylocal[2]]!=1) kazendjjj=1;
    else kazendjjj=0;

    if(tablesize[whichtable][vartypearraylocal[3]]!=1) kazendkkk=1;
    else kazendkkk=0;

    if(tablesize[whichtable][vartypearraylocal[4]]!=1) kazendlll=1;
    else kazendlll=0;

  }

  // overrides (must be placed outside if(repeatedeos==0) conditional so always done since kazendlll,kazendmmm, etc. are stored per whichd and not also per whichtablesubtype
  // only extra table is function of Ynu or H if whichdatatype==4
  if(WHICHDATATYPEGENERAL==4 && (whichtablesubtype!=SUBTYPEEXTRA)){ kazll=leos=kazendlll=kazdl[1]=kazmm=meos=kazendmmm=kazdm[1]=0; kazdl[0]=kazdm[0]=1.0;}


    

  ///////////////////////
  //
  // initialize cumulative result and distance for normalization
  //
  ///////////////////////
  for(coli=0;coli<numcols;coli++){
    if(iffun[coli]){
      for(iii=kazstartiii;iii<=kazendiii;iii++){
	tfptr[coli][iii]=0.0;
	tdist[coli][iii]=0.0;
      }
    }
  }



  ///////////////////////
  //
  // Get cumulative results
  //
  // Loop over nearby table values and determine bi-linearly interpolated value
  // 4-D means 2^4=16 positions
  // 2-D means 2^2=4 positions
  // get 3 values as function of density
  //
  ///////////////////////
  mmm=0;
  for(iii=kazstartiii;iii<=kazendiii;iii++) for(jjj=0;jjj<=kazendjjj;jjj++)for(kkk=0;kkk<=kazendkkk;kkk++)for(lll=0;lll<=kazendlll;lll++){


#if(CHECKIFVALIDEOSDATA)
    if(degentable==0){
      // don't use values of table that have no inversion to temperature
      get_arrays_eostable_direct_temperature(whichdegenfun,whichtable,kazmm+mmm,kazll+lll,kazkk+kkk,kazjj+jjj,kazii+iii,&tempcheck);
    }
    //    dualfprintf(fail_file,"temp=%21.15g\n",tempcheck);
    if(degentable==1 || tempcheck>invalidtempcode) // Avoid invalid inversions if T>Tbad, but only deal with temperature if degentable==0
#else
    if(1)
#endif
    {

      // get value
      get_arrays_eostable_direct(whichdegenfun,degentable,whichtable,whichtablesubtype,iffun,kazmm+mmm,kazll+lll,kazkk+kkk,kazjj+jjj,kazii+iii,tempanswers);

      // get distance
      dist = kazdj[jjj]*kazdk[kkk]*kazdl[lll];


      // non-pipelined section
      // first do offset and see if can use log for all points
      for(coli=0;coli<numcols;coli++){
	if(iffun[coli]){
	  if(degentable==1) offsetquant2_general(whichdegenfun, quant1, tempanswers[coli], &tempanswers[coli]);
  
	  if(shouldloginterp[coli]){
	    if(tempanswers[coli]<=0.0){
	      // then don't include since if should have been able to do log but couldn't, then remove point from interpolation
	    }
	    else{
	      // should and can do log
	      tfptr[coli][iii] += log10(tempanswers[coli])*dist; // interpolated answer that will replace answer[coli] in the end
	      tdist[coli][iii] += dist;
	    }
	  }// if should do log interp
	  else{
	    // shouldn't do log
	    tfptr[coli][iii] += tempanswers[coli]*dist;
	    tdist[coli][iii] += dist;
	  }// else if not log interp
	}// end if doing this coli
      }// end over coli      
    }// end if good temperature or doing degentable
    else{
      // no good temperature for this point
    }
  }// end loop over dimensions



  //////////////////////
  //
  // Set that default return that no bad lookups
  //
  //////////////////////
  for(coli=0;coli<numcols;coli++){
    badlookups[coli]=0;
  }




  ////////////////////////////
  //
  // finally normalize
  //
  ////////////////////////////
  for(iii=kazstartiii;iii<=kazendiii;iii++) {

    // see if any bad lookups
    int badany=0;
    for(coli=0;coli<numcols;coli++){
      if(iffun[coli]){
	// check if any bad lookups
	if(tdist[coli][iii]==0.0) badany++;
      }
    }
      
  
    // never found good temperature, so look up nearest_dumb()
    // GODMARK: Could instead reduce from parabolic to lower order since using nearest_dumb() is probably a bad idea
    if(badany){
      get_eos_fromlookup_nearest_dumb(INDEXPARAMETERSTOSUBPASS, repeatedeos, tabledimen, degentable, whichtable, whichtablesubtype, iffun, whichindep, quant1, vartypearraylocal, indexarray, tempanswers, badlookups); // this can overwrite default badlookups from 0->1

      // do offset and logify, but only if badlookup before but goodlookup from nearest_dumb()
      for(coli=0;coli<numcols;coli++){
	if(iffun[coli]){
	  if(tdist[coli][iii]==0.0 && badlookups[coli]==0){ // only do if fixing no-temperature lookup
	    if(degentable==1) offsetquant2_general(whichdegenfun, quant1, tempanswers[coli], &tempanswers[coli]);
	    if(shouldloginterp[coli]){
	      tfptr[coli][iii] = log10(tempanswers[coli]);
	    }
	    else{
	      tfptr[coli][iii] = tempanswers[coli];
	    }
	  }// end if badlookups==0
	}// end iffun==1
      }//end over coli
    }// end if any bad values that need to be looked-up with nearest_dumb()


    // for those cases when found good temperature, normalize based upon distance
    for(coli=0;coli<numcols;coli++){
      if(iffun[coli]){
	if(tdist[coli][iii]!=0.0){ // only do if found good temperature (and can't be badlookup)
	  tfptr[coli][iii] /= tdist[coli][iii];
	}
      } // end iffun==1
    }//end over coli

  } // end over iii


  //////////////////////////////
  //
  // now use 3 data points to get density-parabolic distribution and value at ieos
  // have tfptr[iii] @ iii=-1,0,1 with i=0 meaning ii and offset being ROUND2INT(ieos)
  // Form parabolic answer
  //
  //////////////////////////////
  for(coli=0;coli<numcols;coli++){
    if(iffun[coli]){
      if(badlookups[coli]==0){
	xmx0 = (ieos-(FTYPE)kazii);
	AA = 0.5*(tfptr[coli][1]-tfptr[coli][-1]);
	BB = 0.5*(tfptr[coli][1]+tfptr[coli][-1]-2.0*tfptr[coli][0]);
	answers[coli] = tfptr[coli][0] + AA*xmx0 + BB*xmx0*xmx0;
	
	if(shouldloginterp[coli]){
	  answers[coli]=pow(10.0,answers[coli]);
	}
	if(degentable==1) offsetquant2_general_inverse(whichdegenfun, quant1, answers[coli], &answers[coli]);
      }
    }
  }


  //////////////////////
  //
  // Return that no bad failure
  //
  //////////////////////
  return(0);

}




// full 5D linear interpolation
static int get_eos_fromlookup_linear(INDEXPARAMETERSPROTOTYPES, int repeatedeos, int tabledimen, int degentable, int whichtable, int whichtablesubtype, int *iffun, int whichindep, FTYPE quant1, int *vartypearraylocal, FTYPE *indexarray, FTYPE *answers, int *badlookups)
{
  FTYPE dist,totaldist[MAXEOSPIPELINE];
  FTYPE tempanswers[MAXEOSPIPELINE];
  FTYPE totalanswers[MAXEOSPIPELINE];
  FTYPE tempcheck;
  int iii,jjj,kkk,lll,mmm;
  int whichdegenfun;
  int shouldloginterp[MAXEOSPIPELINE];
  FTYPE ieos,jeos,keos,leos,meos;
  int qi;
  int numcols,coli;



  // get number of columns possible for this whichtablesubtype
  numcols = numcolintablesubtype[whichtablesubtype];


  ieos=indexarray[1];
  jeos=indexarray[2];
  keos=indexarray[3];
  leos=indexarray[4];
  meos=indexarray[5];
  



  // definition consistent with numerical assignments of indecies of arrays
  whichdegenfun = whichindep-1;




  // determine if should do log interpolation
  get_dologinterp_subtype_wrapper(degentable, whichtablesubtype, numcols, shouldloginterp);


  if(repeatedeos==0){
    // assume jeos, keos, leos set correctly for given degentable and tabledimen (i.e. are 0 or integers in those cases so d's computed correctly)
    kazii=(int)ieos;
    //if(degentable) kazjj=0; else kazjj=(int)jeos; // force to be 0
    kazjj=(int)jeos; // (int)(0.0) = 0, so ok as is
    kazkk=(int)keos;
    kazll=(int)leos;
    kazmm=(int)meos;
	
    kazdi[1]=ieos-(FTYPE)kazii;
    kazdi[0]=1.0-kazdi[1];
    kazdj[1]=jeos-(FTYPE)kazjj;
    kazdj[0]=1.0-kazdj[1];
    kazdk[1]=keos-(FTYPE)kazkk;
    kazdk[0]=1.0-kazdk[1];
    kazdl[1]=leos-(FTYPE)kazll;
    kazdl[0]=1.0-kazdl[1];
    kazdm[1]=meos-(FTYPE)kazmm;
    kazdm[0]=1.0-kazdm[1];


    // set range of loops for different table types
    // tabledimen overrules table type (i.e. take section out of fuller table -- assumed to be 0 index)
    // GODMARK: these things couuld be stored as functions of whichtable/degentable/tabledimen
    if(tablesize[whichtable][vartypearraylocal[1]]!=1) kazendiii=1;
    else kazendiii=0;

    if(degentable==0 && tablesize[whichtable][vartypearraylocal[2]]!=1) kazendjjj=1;
    else kazendjjj=0;

    if(tablesize[whichtable][vartypearraylocal[3]]!=1) kazendkkk=1;
    else kazendkkk=0;
    
    if(tablesize[whichtable][vartypearraylocal[4]]!=1) kazendlll=1;
    else kazendlll=0;

    if(tablesize[whichtable][vartypearraylocal[5]]!=1) kazendmmm=1;
    else kazendmmm=0;

    // avoid using data beyond table (i.e. reduce to nearest neighbor at edges of table as long as within table)
    if(kazii+kazendiii>=tablesize[whichtable][vartypearraylocal[1]]) kazendiii=kazii;
    if(kazjj+kazendjjj>=tablesize[whichtable][vartypearraylocal[2]]) kazendjjj=kazjj;
    if(kazkk+kazendkkk>=tablesize[whichtable][vartypearraylocal[3]]) kazendkkk=kazkk;
    if(kazll+kazendlll>=tablesize[whichtable][vartypearraylocal[4]]) kazendlll=kazll;
    if(kazmm+kazendmmm>=tablesize[whichtable][vartypearraylocal[5]]) kazendmmm=kazmm;

  }

  // overrides (must be placed outside if(repeatedeos==0) conditional so always done since kazendlll,kazendmmm, etc. are stored per whichd and not also per whichtablesubtype
  // only extra table is function of Ynu or H if whichdatatype==4
  if(WHICHDATATYPEGENERAL==4 && (whichtablesubtype!=SUBTYPEEXTRA)){leos=meos=kazll=kazmm=kazdl[1]=kazdm[1]=kazendlll=kazendmmm=0; kazdl[0]=kazdm[0]=1.0;}



    
  // Loop over nearby table values and determine bi-linearly interpolated value
  // 4-D means 2^4=16 positions
  // 2-D means 2^2=4 positions

  // initialize distance and value
  for(coli=0;coli<numcols;coli++){
    if(iffun[coli]){
      totaldist[coli]=0.0;
      totalanswers[coli]=0.0;
      tempanswers[coli]=0.0;
    }
  }


  // LOOP:
#define LOOPOVERKAZINDEX  for(mmm=0;mmm<=kazendmmm;mmm++)for(lll=0;lll<=kazendlll;lll++)for(kkk=0;kkk<=kazendkkk;kkk++)for(jjj=0;jjj<=kazendjjj;jjj++)for(iii=0;iii<=kazendiii;iii++)


  LOOPOVERKAZINDEX{


#if(CHECKIFVALIDEOSDATA)
    if(degentable==0){
      // don't use values of table that have no inversion to temperature
      get_arrays_eostable_direct_temperature(whichdegenfun,whichtable,kazmm+mmm,kazll+lll,kazkk+kkk,kazjj+jjj,kazii+iii,&tempcheck);
    }
    if(degentable==1 || tempcheck>invalidtempcode) // Avoid invalid inversions if T>Tbad, but only deal with temperature if degentable==0
#else
    if(1)
#endif
    {
      // get value
      get_arrays_eostable_direct(whichdegenfun,degentable,whichtable,whichtablesubtype,iffun,kazmm+mmm,kazll+lll,kazkk+kkk,kazjj+jjj,kazii+iii,tempanswers);

      // get distance to value
      dist = kazdi[iii]*kazdj[jjj]*kazdk[kkk]*kazdl[lll]*kazdl[mmm];

      // non-pipelined section:
      for(coli=0;coli<numcols;coli++){
	if(iffun[coli]){
	  if(degentable==1) offsetquant2_general(whichdegenfun, quant1, tempanswers[coli], &tempanswers[coli]);
  
	  if(shouldloginterp[coli]){
	    if(tempanswers[coli]<=0.0){
	      // then don't include since if should have been able to do log but couldn't, then remove point from interpolation
	    }
	    else{
	      // should and can do log
	      totalanswers[coli] += log10(tempanswers[coli])*dist; // interpolated answer that will replace answer[coli] in the end
	      totaldist[coli] += dist;
	    }
	  }// if should do log interp
	  else{
	    // shouldn't do log
	    totalanswers[coli] += tempanswers[coli]*dist;
	    totaldist[coli] += dist;
	  }// else if not log interp
	}// end if doing this coli
      }// end over coli      

    }// end if good temperature or doing degentable
    else{
      
      // out of bounds for temperature
      
    }// end if out of bounds for temperature
  }// end loop over dimensions



  ////////////////////////////
  //
  // finally normalize
  //
  ////////////////////////////
  // non-pipelined section:
  for(coli=0;coli<numcols;coli++){
    if(iffun[coli]){
      if(totaldist[coli]==0.0){

	// return that was bad lookup
	badlookups[coli]=1;

	if(debugfail>=2){
	  dualfprintf(fail_file,"No valid data points in table despite within table.\n");
	  dualfprintf(fail_file,"kaz=%d %d %d %d : %d %d : %d %d : %d %d : %d %d\n",kazii,kazjj,kazkk,kazll,kazstartiii,kazendiii,kazstartjjj,kazendjjj,kazstartkkk,kazendkkk,kazstartlll,kazendlll);
	}

      }
      else{
	// return that was good lookup
	badlookups[coli]=0;

	// get normalized answer
	answers[coli] = totalanswers[coli]/totaldist[coli];

	// invert log interpolation
	if(shouldloginterp[coli]) answers[coli]=pow(10.0,answers[coli]);
	
	// invert offset
	if(degentable==1) offsetquant2_general_inverse(whichdegenfun, quant1, answers[coli], &answers[coli]);
      }// end else if totaldist!=0.0
    }// end iffun==1
  }// end over coli



  ///////////////////
  //
  // return that no bad failure
  //
  ///////////////////
  return(0);

}










// nearest neighbor interpolation but caution to temperature defined or not
static int get_eos_fromlookup_nearest(INDEXPARAMETERSPROTOTYPES, int repeatedeos, int tabledimen, int degentable, int whichtable, int whichtablesubtype, int *iffun, int whichindep, FTYPE quant1, int *vartypearraylocal, FTYPE *indexarray, FTYPE *answers, int *badlookups)
{
  FTYPE totaldist;
  FTYPE tempcheck;
  int iii,jjj,kkk,lll,mmm;
  int whichdegenfun;
  FTYPE ieos,jeos,keos,leos,meos;
  int numcols,coli;


  // get number of columns possible for this whichtablesubtype
  numcols = numcolintablesubtype[whichtablesubtype];


  ieos=indexarray[1];
  jeos=indexarray[2];
  keos=indexarray[3];
  leos=indexarray[4];
  meos=indexarray[5];




  // definition consistent with numerical assignments of indecies of arrays
  whichdegenfun = whichindep-1;




  if(repeatedeos==0){
    // assume jeos, keos, leos set correctly for given degentable and tabledimen (i.e. are 0 or integers in those cases so d's computed correctly)
    // if only choosing 1 value, then choose rounded version
    kazii=ROUND2INT(ieos);
    kazjj=ROUND2INT(jeos);
    kazkk=ROUND2INT(keos);
    kazll=ROUND2INT(leos);
    kazmm=ROUND2INT(meos);
	

    if(kazii>=tablesize[whichtable][vartypearraylocal[1]]) kazii=tablesize[whichtable][vartypearraylocal[1]]-1;
    if(kazjj>=tablesize[whichtable][vartypearraylocal[2]]) kazjj=tablesize[whichtable][vartypearraylocal[2]]-1;
    if(kazkk>=tablesize[whichtable][vartypearraylocal[3]]) kazkk=tablesize[whichtable][vartypearraylocal[3]]-1;
    if(kazll>=tablesize[whichtable][vartypearraylocal[4]]) kazll=tablesize[whichtable][vartypearraylocal[4]]-1;
    if(kazmm>=tablesize[whichtable][vartypearraylocal[5]]) kazmm=tablesize[whichtable][vartypearraylocal[5]]-1;

    // set range of loops for different table types
    // tabledimen overrules table type (i.e. take section out of fuller table -- assumed to be 0 index)
    // GODMARK: these things couuld be stored as functions of whichtable/degentable/tabledimen
    if(tablesize[whichtable][vartypearraylocal[1]]!=1) kazendiii=1;
    else kazendiii=0;

    if(degentable==0 && tablesize[whichtable][vartypearraylocal[2]]!=1) kazendjjj=1;
    else kazendjjj=0;

    if(tablesize[whichtable][vartypearraylocal[3]]!=1) kazendkkk=1;
    else kazendkkk=0;
    
    if(tablesize[whichtable][vartypearraylocal[4]]!=1) kazendlll=1;
    else kazendlll=0;

    if(tablesize[whichtable][vartypearraylocal[5]]!=1) kazendmmm=1;
    else kazendmmm=0;

  }

  
  // overrides (must be placed outside if(repeatedeos==0) conditional so always done since kazendlll,kazendmmm, etc. are stored per whichd and not also per whichtablesubtype
  // only extra table is function of Ynu or H if whichdatatype==4
  if(WHICHDATATYPEGENERAL==4 && (whichtablesubtype!=SUBTYPEEXTRA)){leos=meos=kazll=kazmm=kazendlll=kazendmmm=0;}



  ////////////
  //
  // lookup value 
  // try ROUND2INT version first
  //
  ////////////
  int gotgood=0;
#if(CHECKIFVALIDEOSDATA)
  if(degentable==0){
    // don't use values of table that have no inversion to temperature
    get_arrays_eostable_direct_temperature(whichdegenfun,whichtable,kazmm,kazll,kazkk,kazjj,kazii,&tempcheck);
  }
  if(degentable==1 || tempcheck>invalidtempcode) // Avoid invalid inversions if T>Tbad, but only deal with temperature if degentable==0
#else
  if(1)
#endif
  {
    get_arrays_eostable_direct(whichdegenfun,degentable,whichtable,whichtablesubtype,iffun,kazmm,kazll,kazkk,kazjj,kazii,answers);
    gotgood=1;
  }


  ////////////
  //
  // if ROUND2INT looked-up in invalid part of EOS domain, try first-found neighbor
  //
  ////////////

  if(gotgood==0){
#define LOOPKAZIJKLM for(iii=0;iii<=kazendiii;iii++)for(jjj=0;jjj<=kazendjjj;jjj++)for(kkk=0;kkk<=kazendkkk;kkk++)for(lll=0;lll<=kazendlll;lll++)for(mmm=0;mmm<=kazendmmm;mmm++)

    // loop over to ensure obtained best nearest neighbor    
    LOOPKAZIJKLM{
    
#if(CHECKIFVALIDEOSDATA)
      if(degentable==0){
	// don't use values of table that have no inversion to temperature
	get_arrays_eostable_direct_temperature(whichdegenfun,whichtable,kazmm,kazll,kazkk,kazjj,kazii,&tempcheck);
      }
      if(degentable==1 || tempcheck>invalidtempcode) // Avoid invalid inversions if T>Tbad, but only deal with temperature if degentable==0
#else
      if(1)
#endif
      {
	get_arrays_eostable_direct(whichdegenfun,degentable,whichtable,whichtablesubtype,iffun,kazmm,kazll,kazkk,kazjj,kazii,answers);

	// terminate all 4 loops once have a single value
	iii=kazendiii;
	jjj=kazendjjj;
	kkk=kazendkkk;
	lll=kazendlll;
	mmm=kazendmmm;
	// set that got good
	gotgood=1;

      }// end if did find a good temperature or doing degentable
    }// end loop over dimensions
  }// else if ROUNT2INT failed so picking nearest neighbor  


  ////////////
  //
  // If still failed to find a valid EOS value, then report/return that bad lookup
  //
  ////////////
  if(gotgood==0){
    // return that was bad lookup
    for(coli=0;coli<numcols;coli++) badlookups[coli]=1;
  }
  else{
    // return that was not bad lookup
    for(coli=0;coli<numcols;coli++) badlookups[coli]=0;
  }

  /////////////
  //
  // return that no bad failure
  //
  /////////////
  return(0);


}


// nearest neighbor interpolation with no temperature check
// note that kazii-...-kazmm are globals in kazfulleos.defsglobalprivate.h
// They are globals because need to do repeated lookup that is multi-thread safe, so each thread as a global version of these
static int get_eos_fromlookup_nearest_dumb(INDEXPARAMETERSPROTOTYPES, int repeatedeos, int tabledimen, int degentable, int whichtable, int whichtablesubtype, int *iffun, int whichindep, FTYPE quant1, int *vartypearraylocal, FTYPE *indexarray, FTYPE *answers, int *badlookups)
{
  int whichdegenfun;
  int numcols,coli;


  // get number of columns possible for this whichtablesubtype
  numcols = numcolintablesubtype[whichtablesubtype];

  // definition consistent with numerical assignments of indeces of arrays
  whichdegenfun = whichindep-1;


  if(repeatedeos==0){

#if(0)
    kazii=(int)indexarray[1];
    kazjj=(int)indexarray[2];
    kazkk=(int)indexarray[3];
    kazll=(int)indexarray[4];
    kazmm=(int)indexarray[5];
#else
    // if only choosing 1 value, then choose rounded version
    kazii=ROUND2INT(indexarray[1]);
    kazjj=ROUND2INT(indexarray[2]);
    kazkk=ROUND2INT(indexarray[3]);
    kazll=ROUND2INT(indexarray[4]);
    kazmm=ROUND2INT(indexarray[5]);
#endif


    if(kazii>=tablesize[whichtable][vartypearraylocal[1]]) kazii=tablesize[whichtable][vartypearraylocal[1]]-1;
    if(kazjj>=tablesize[whichtable][vartypearraylocal[2]]) kazjj=tablesize[whichtable][vartypearraylocal[2]]-1;
    if(kazkk>=tablesize[whichtable][vartypearraylocal[3]]) kazkk=tablesize[whichtable][vartypearraylocal[3]]-1;
    if(kazll>=tablesize[whichtable][vartypearraylocal[4]]) kazll=tablesize[whichtable][vartypearraylocal[4]]-1;
    if(kazmm>=tablesize[whichtable][vartypearraylocal[5]]) kazmm=tablesize[whichtable][vartypearraylocal[5]]-1;
  }

  // overrides (must be placed outside if(repeatedeos==0) conditional so always done since kazendlll,kazendmmm, etc. are stored per whichd and not also per whichtablesubtype
  // only extra table is function of Ynu or H if whichdatatype==4
  if(WHICHDATATYPEGENERAL==4 && (whichtablesubtype!=SUBTYPEEXTRA)){kazll=kazmm=kazendlll=kazendmmm=0;}

  // get value
  get_arrays_eostable_direct(whichdegenfun,degentable,whichtable,whichtablesubtype,iffun,kazmm,kazll,kazkk,kazjj,kazii,answers);
  
  // return that was not bad lookup
  for(coli=0;coli<numcols;coli++) badlookups[coli]=0;
  
  // return that no bad failure
  return(0);

}





//#define thiseostable(name,fun,index5,index4,index3,index2,index1) EOSMAC(name,fun,index5,index4,index3,index2,index1)
//#define index2array(name,fun) name[fun][indexarray[5]][indexarray[4]][indexarray[3]][indexarray[2]][indexarray[1]]


// kinda like get_arrays_eostable() but directly grab multiple columns of data from a given table subtype
// simple and simplezoom table code is just copy-paste from full table with fulltable->simpletable and fulltable->simplezoomtable
static void get_arrays_eostable_direct(int whichd, int whichdegen, int whichtable, int whichtablesubtype, int *iffun, int mmm, int lll, int kkk, int jjj, int iii, double *values)
{
  int numcols,coli;
  //  int whichdlocal;


#if(PRODUCTION==0)
  // overrides based upon true memory allocated as only required
  if(WHICHDATATYPEGENERAL==4){
    if(whichtablesubtype!=SUBTYPEEXTRA){
      if(mmm!=0 || lll!=0){
	// Then non-extras have no Ynu or H dependence
	// should ensure loops don't use non-zero mmm or lll
	dualfprintf(fail_file,"Bad get_arrays_eostable_direct: mmm=%d lll=%d whichtablesubtype=%d\n",mmm,lll,whichtablesubtype);
	myexit(4633221);
      }
    }
  }
  if(whichtablesubtype==SUBTYPEDEGEN){
    if(jjj!=0){
      dualfprintf(fail_file,"Bad get_arrays_eostable_direct: degen: mmm=%d lll=%d whichtablesubtype=%d\n",mmm,lll,whichtablesubtype);
      myexit(4633222);
    }
  }
#endif



  // get number of columns possible for this whichtablesubtype
  numcols = numcolintablesubtype[whichtablesubtype];
  //  whichdlocal = whichdintablesubtype[whichtablesubtype]; // could ensure that whichd==whichdlocal for non-degen and non-temp




  if(0){
  }
#if(ALLOWFULLTABLE==1)
  else if(whichtable==FULLTABLE||whichtable==EXTRAFULLTABLE){
    if(whichdegen==ISNOTDEGENTABLE){
      // degentable not dealt with here
      if(whichtable==FULLTABLE){
	// assume whichtablesubtype==SUBTYPEEXTRA is not reached if WHICHDATATYPEGENERAL==4 (controlled by which_eostable())
	if(whichtablesubtype==SUBTYPESTANDARD){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltablestandard,0,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSSTANDARD);}
	else if(whichtablesubtype==SUBTYPEGUESS){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltableguess,0,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSGUESS);}
	else if(whichtablesubtype==SUBTYPEDISS){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltablediss,0,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSDISS);}
	else if(whichtablesubtype==SUBTYPEDP){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltabledp,0,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSDP);}
	else if(whichtablesubtype==SUBTYPESDEN){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltablesden,0,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSSDEN);}
	else if(whichtablesubtype==SUBTYPESSPEC){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltablesspec,0,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSSSPEC);}
	else if(whichtablesubtype==SUBTYPEPOFCHI){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltablepofchi,0,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSPOFCHI);}
	// temp table special like degen table such that only 1 returned temp associated with whichd value
	else if(whichtablesubtype==SUBTYPETEMP){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltabletemp,whichd,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSTEMP);} // whichd used
	else if(whichtablesubtype==SUBTYPEEXTRA){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltableextra,0,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSEXTRA);}
	// in principle can request degen data directly, so allow this
	else if(whichtablesubtype==SUBTYPEDEGEN){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltabledegen,whichd,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSDEGEN);} // whichd used
      }
      else if(whichtable==EXTRAFULLTABLE){
	// reach here if WHICHDATATYPEGENERAL==4 and want extras(controlled by which_eostable())
	if(whichtablesubtype==SUBTYPEEXTRA){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltableextra,0,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSEXTRA);}
	// in principle can request degen data directly, so allow this
	else if(whichtablesubtype==SUBTYPEDEGEN){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltableextradegen,whichd,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSDEGEN);} // whichd used
	else if(whichtablesubtype==SUBTYPETEMP){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltableextratemp,whichd,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSTEMP);}
      }
    }
    else{
      if(whichtable==FULLTABLE){
	// then ignore whichtablesubtype and iffun[] since must be degen lookup and need all quantities
	// get basic utotoffset
	// starts at 0 as coli=0 is start.  Later mapped to whichfun type
	values[0] = EOSMAC(eosfulltabledegen,whichd,mmm,lll,kkk,0,iii,EOSOFFSET);
	if(utotdegencut[whichtable]>=2){
	  // then also get other things
	  values[1] = EOSMAC(eosfulltabledegen,whichd,mmm,lll,kkk,0,iii,EOSIN);
	  values[2] = EOSMAC(eosfulltabledegen,whichd,mmm,lll,kkk,0,iii,EOSOUT);
	}
      }
      else if(whichtable==EXTRAFULLTABLE){
	// reach here if WHICHDATATYPEGENERAL==4 and want extras(controlled by which_eostable())
	values[0] = EOSMAC(eosfulltableextradegen,whichd,mmm,lll,kkk,0,iii,EOSOFFSET);
	if(utotdegencut[whichtable]>=2){
	  values[1] = EOSMAC(eosfulltableextradegen,whichd,mmm,lll,kkk,0,iii,EOSIN);
	  values[2] = EOSMAC(eosfulltableextradegen,whichd,mmm,lll,kkk,0,iii,EOSOUT);
	}
      }
    }// end else if degentable      
  }// end if fulltable
#endif
#if(ALLOWSIMPLETABLE==1)
  else if(whichtable==SIMPLETABLE||whichtable==EXTRASIMPLETABLE){
    if(whichdegen==ISNOTDEGENTABLE){
      // degentable not dealt with here
      if(whichtable==SIMPLETABLE){
	// assume whichtablesubtype==SUBTYPEEXTRA is not reached if WHICHDATATYPEGENERAL==4 (controlled by which_eostable())
	if(whichtablesubtype==SUBTYPESTANDARD){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletablestandard,0,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSSTANDARD);}
	else if(whichtablesubtype==SUBTYPEGUESS){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletableguess,0,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSGUESS);}
	else if(whichtablesubtype==SUBTYPEDISS){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletablediss,0,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSDISS);}
	else if(whichtablesubtype==SUBTYPEDP){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletabledp,0,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSDP);}
	else if(whichtablesubtype==SUBTYPESDEN){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletablesden,0,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSSDEN);}
	else if(whichtablesubtype==SUBTYPESSPEC){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletablesspec,0,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSSSPEC);}
	else if(whichtablesubtype==SUBTYPEPOFCHI){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletablepofchi,0,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSPOFCHI);}
	// temp table special like degen table such that only 1 returned temp associated with whichd value
	else if(whichtablesubtype==SUBTYPETEMP){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletabletemp,whichd,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSTEMP);} // whichd used
	else if(whichtablesubtype==SUBTYPEEXTRA){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletableextra,0,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSEXTRA);}
	// in principle can request degen data directly, so allow this
	else if(whichtablesubtype==SUBTYPEDEGEN){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletabledegen,whichd,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSDEGEN);} // whichd used
      }
      else if(whichtable==EXTRASIMPLETABLE){
	// reach here if WHICHDATATYPEGENERAL==4 and want extras(controlled by which_eostable())
	if(whichtablesubtype==SUBTYPEEXTRA){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletableextra,0,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSEXTRA);}
	// in principle can request degen data directly, so allow this
	else if(whichtablesubtype==SUBTYPEDEGEN){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletableextradegen,whichd,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSDEGEN);} // whichd used
	else if(whichtablesubtype==SUBTYPETEMP){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletableextratemp,whichd,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSTEMP);}
      }
    }
    else{
      if(whichtable==SIMPLETABLE){
	// then ignore whichtablesubtype and iffun[] since must be degen lookup and need all quantities
	// get basic utotoffset
	// starts at 0 as coli=0 is start.  Later mapped to whichfun type
	values[0] = EOSMAC(eossimpletabledegen,whichd,mmm,lll,kkk,0,iii,EOSOFFSET);
	if(utotdegencut[whichtable]>=2){
	  // then also get other things
	  values[1] = EOSMAC(eossimpletabledegen,whichd,mmm,lll,kkk,0,iii,EOSIN);
	  values[2] = EOSMAC(eossimpletabledegen,whichd,mmm,lll,kkk,0,iii,EOSOUT);
	}
      }
      else if(whichtable==EXTRASIMPLETABLE){
	// reach here if WHICHDATATYPEGENERAL==4 and want extras(controlled by which_eostable())
	values[0] = EOSMAC(eossimpletableextradegen,whichd,mmm,lll,kkk,0,iii,EOSOFFSET);
	if(utotdegencut[whichtable]>=2){
	  values[1] = EOSMAC(eossimpletableextradegen,whichd,mmm,lll,kkk,0,iii,EOSIN);
	  values[2] = EOSMAC(eossimpletableextradegen,whichd,mmm,lll,kkk,0,iii,EOSOUT);
	}
      }
    }// end else if degentable      
  }// end if simpletable
#endif
#if(ALLOWSIMPLEZOOMTABLE==1)
  else if(whichtable==SIMPLEZOOMTABLE||whichtable==EXTRASIMPLEZOOMTABLE){
    if(whichdegen==ISNOTDEGENTABLE){
      // degentable not dealt with here
      if(whichtable==SIMPLEZOOMTABLE){
	// assume whichtablesubtype==SUBTYPEEXTRA is not reached if WHICHDATATYPEGENERAL==4 (controlled by which_eostable())
	if(whichtablesubtype==SUBTYPESTANDARD){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtablestandard,0,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSSTANDARD);}
	else if(whichtablesubtype==SUBTYPEGUESS){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtableguess,0,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSGUESS);}
	else if(whichtablesubtype==SUBTYPEDISS){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtablediss,0,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSDISS);}
	else if(whichtablesubtype==SUBTYPEDP){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtabledp,0,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSDP);}
	else if(whichtablesubtype==SUBTYPESDEN){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtablesden,0,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSSDEN);}
	else if(whichtablesubtype==SUBTYPESSPEC){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtablesspec,0,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSSSPEC);}
	else if(whichtablesubtype==SUBTYPEPOFCHI){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtablepofchi,0,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSPOFCHI);}
	// temp table special like degen table such that only 1 returned temp associated with whichd value
	else if(whichtablesubtype==SUBTYPETEMP){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtabletemp,whichd,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSTEMP);} // whichd used
	else if(whichtablesubtype==SUBTYPEEXTRA){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtableextra,0,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSEXTRA);}
	// in principle can request degen data directly, so allow this
	else if(whichtablesubtype==SUBTYPEDEGEN){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtabledegen,whichd,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSDEGEN);} // whichd used
      }
      else if(whichtable==EXTRASIMPLEZOOMTABLE){
	// reach here if WHICHDATATYPEGENERAL==4 and want extras(controlled by which_eostable())
	if(whichtablesubtype==SUBTYPEEXTRA){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtableextra,0,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSEXTRA);}
	// in principle can request degen data directly, so allow this
	else if(whichtablesubtype==SUBTYPEDEGEN){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtableextradegen,whichd,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSDEGEN);} // whichd used
	else if(whichtablesubtype==SUBTYPETEMP){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtableextratemp,whichd,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSTEMP);}
      }
    }
    else{
      if(whichtable==SIMPLEZOOMTABLE){
	// then ignore whichtablesubtype and iffun[] since must be degen lookup and need all quantities
	// get basic utotoffset
	// starts at 0 as coli=0 is start.  Later mapped to whichfun type
	values[0] = EOSMAC(eossimplezoomtabledegen,whichd,mmm,lll,kkk,0,iii,EOSOFFSET);
	if(utotdegencut[whichtable]>=2){
	  // then also get other things
	  values[1] = EOSMAC(eossimplezoomtabledegen,whichd,mmm,lll,kkk,0,iii,EOSIN);
	  values[2] = EOSMAC(eossimplezoomtabledegen,whichd,mmm,lll,kkk,0,iii,EOSOUT);
	}
      }
      else if(whichtable==EXTRASIMPLEZOOMTABLE){
	// reach here if WHICHDATATYPEGENERAL==4 and want extras(controlled by which_eostable())
	values[0] = EOSMAC(eossimplezoomtableextradegen,whichd,mmm,lll,kkk,0,iii,EOSOFFSET);
	if(utotdegencut[whichtable]>=2){
	  values[1] = EOSMAC(eossimplezoomtableextradegen,whichd,mmm,lll,kkk,0,iii,EOSIN);
	  values[2] = EOSMAC(eossimplezoomtableextradegen,whichd,mmm,lll,kkk,0,iii,EOSOUT);
	}
      }
    }// end else if degentable      
  }// end if simplezoomtable
#endif



}






// like get_arrays_eostable_direct(), but only for temperature that ignores associated variable actually looking up
static void get_arrays_eostable_direct_temperature(int whichd, int whichtable, int mmm, int lll, int kkk, int jjj, int iii, double *values)
{
  int numcols,coli;


  // overrides based upon true memory allocated as only required
  // override instead of pre-condition because this function is called inside loop over other tables
  if(WHICHDATATYPEGENERAL==4){
    // Then temp has no Ynu or H dependence
    mmm=lll=0;
  }


  // get number of columns possible for this whichtablesubtype
  numcols = numcolintablesubtype[SUBTYPETEMP];

  if(0){
  }
#if(ALLOWFULLTABLE==1)
  else if(whichtable==FULLTABLE){
    for(coli=0;coli<numcols;coli++) values[coli]=EOSMAC(eosfulltabletemp,whichd,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSTEMP);
  }// end if fulltable
  else if(whichtable==EXTRAFULLTABLE){
    for(coli=0;coli<numcols;coli++) values[coli]=EOSMAC(eosfulltableextratemp,whichd,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSTEMP);
  }// end if extrafulltable
#endif
  // simple table is just copy-paste from full table with fulltable->simpletable
#if(ALLOWSIMPLETABLE==1)
  else if(whichtable==SIMPLETABLE){
    for(coli=0;coli<numcols;coli++) values[coli]=EOSMAC(eossimpletabletemp,whichd,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSTEMP);
  }// end if simpletable
  else if(whichtable==EXTRASIMPLETABLE){
    for(coli=0;coli<numcols;coli++) values[coli]=EOSMAC(eossimpletableextratemp,whichd,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSTEMP);
  }// end if extrasimpletable
#endif
  // simplezoom table is just copy-paste from full table with fulltable->simplezoomtable
#if(ALLOWSIMPLEZOOMTABLE==1)
  else if(whichtable==SIMPLEZOOMTABLE){
    for(coli=0;coli<numcols;coli++) values[coli]=EOSMAC(eossimplezoomtabletemp,whichd,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSTEMP);
  }// end if simplezoomtable
  else if(whichtable==EXTRASIMPLEZOOMTABLE){
    for(coli=0;coli<numcols;coli++) values[coli]=EOSMAC(eossimplezoomtableextratemp,whichd,mmm,lll,kkk,jjj,iii,coli+FIRSTEOSTEMP);
  }// end if extrasimplezoomtable
#endif


}














