

// interpolation stuff
// Should be only stuff that actually accesses global EOS arrays during simulation

// No longer store globally repated "kaz" indices since very simple computations and prone to error in repeated lookup


// for function prototypes associated with function calls
//#define INDEXPARAMETERSPROTOTYPES int kazii, int kazjj, int kazkk, int kazll, int kazmm, \
//    int kaziio, int kazjjo, int kazkko, int kazllo, int kazmmo, \
//    int kazstartiii, int kazstartjjj, int kazstartkkk, int kazstartlll, int kazstartmmm, \
//    int kazendiii, int kazendjjj, int kazendkkk, int kazendlll, int kazendmmm, \
//    FTYPEEOS kazdi[], FTYPEEOS kazdj[], FTYPEEOS kazdk[], FTYPEEOS kazdl[], FTYPEEOS kazdm[]

// for function calls
//#define INDEXPARAMETERSTOPASS kaziiwhichd[whichd],kazjjwhichd[whichd],kazkkwhichd[whichd],kazllwhichd[whichd],kazmmwhichd[whichd], \
//    kaziiowhichd[whichd],kazjjowhichd[whichd],kazkkowhichd[whichd],kazllowhichd[whichd],kazmmowhichd[whichd], \
//    kazstartiiiwhichd[whichd],kazstartjjjwhichd[whichd],kazstartkkkwhichd[whichd],kazstartlllwhichd[whichd],kazstartmmmwhichd[whichd], \
//    kazendiiiwhichd[whichd],kazendjjjwhichd[whichd],kazendkkkwhichd[whichd],kazendlllwhichd[whichd],kazendmmmwhichd[whichd], \
//    kazdiwhichd[whichd],kazdjwhichd[whichd],kazdkwhichd[whichd],kazdlwhichd[whichd],kazdmwhichd[whichd]

// for function calls inside main function call (e.g. _dumb() inside others)
//#define INDEXPARAMETERSTOSUBPASS kazii,kazjj,kazkk,kazll,kazmm, \
//    kaziio,kazjjo,kazkko,kazllo,kazmmo,    \
//    kazstartiii,kazstartjjj,kazstartkkk,kazstartlll,kazstartmmm, \
//    kazendiii,kazendjjj,kazendkkk,kazendlll,kazendmmm,  \
//    kazdi,kazdj,kazdk,kazdl,kazdm

#define INDEXPARAMETERSLOCALDEFS int kazii; int kazjj; int kazkk; int kazll; int kazmm; \
  int kaziio; int kazjjo; int kazkko; int kazllo; int kazmmo;           \
  int kazstartiii; int kazstartjjj; int kazstartkkk; int kazstartlll; int kazstartmmm; \
  int kazendiii; int kazendjjj; int kazendkkk; int kazendlll; int kazendmmm; \
  FTYPEEOS kazdi[2]; FTYPEEOS kazdj[2]; FTYPEEOS kazdk[2]; FTYPEEOS kazdl[2]; FTYPEEOS kazdm[2]



static int get_eos_fromlookup(int repeatedeos, int tabledimen, int degentable, int whichtable, int whichtablesubtype, int *iffun, int whichindep, FTYPE quant1, int *vartypearraylocal, FTYPEEOS *indexarrayvar, FTYPEEOS *answers, int *badlookups);


static int get_eos_fromlookup_nearest_dumb(int repeatedeos, int tabledimen, int degentable, int whichtable, int whichtablesubtype, int *iffun, int whichindep, FTYPE quant1, int *vartypearraylocal, FTYPEEOS *indexarrayvar, FTYPEEOS *answers, int *badlookups);
static int get_eos_fromlookup_nearest(int repeatedeos, int tabledimen, int degentable, int whichtable, int whichtablesubtype, int *iffun, int whichindep, FTYPE quant1, int *vartypearraylocal, FTYPEEOS *indexarrayvar, FTYPEEOS *answers, int *badlookups);
static int get_eos_fromlookup_linear(int repeatedeos, int tabledimen, int degentable, int whichtable, int whichtablesubtype, int *iffun, int whichindep, FTYPE quant1, int *vartypearraylocal, FTYPEEOS *indexarrayvar, FTYPEEOS *answers, int *badlookups);
static int get_eos_fromlookup_parabolicone(int repeatedeos, int tabledimen, int degentable, int whichtable, int whichtablesubtype, int *iffun, int whichindep, FTYPE quant1, int *vartypearraylocal, FTYPEEOS *indexarrayvar, FTYPEEOS *answers, int *badlookups);
static int get_eos_fromlookup_parabolicfull(int repeatedeos, int tabledimen, int degentable, int whichtable, int whichtablesubtype, int *iffun, int whichindep, FTYPE quant1, int *vartypearraylocal, FTYPEEOS *indexarrayvar, FTYPEEOS *answers, int *badlookups);

static void get_arrays_eostable_direct(int whichd, int whichdegen, int whichtable, int whichtablesubtype, int *iffun, int mmm, int lll, int kkk, int jjj, int iii, FTYPEEOS *values);
static void get_arrays_eostable_direct_temperature(int whichd, int whichtable, int mmm, int lll, int kkk, int jjj, int iii, FTYPEEOS *values);









// General EOS interpolation wrapper
static int get_eos_fromlookup(int repeatedeos, int tabledimen, int degentable, int whichtablevar, int whichtablesubtype, int *iffun, int whichindep, FTYPE quant1, int *vartypearraylocal, FTYPEEOS *indexarrayvar, FTYPEEOS *answers, int *badlookups)
{
  int failreturn;
  int whichd;

  ///////////
  //
  // definition consistent with numerical assignments of indices of arrays
  //
  ///////////
  whichd = whichindep-1;

  
#if(0)
  // neartest dumb for normal table
  if(degentable==0){
    failreturn=get_eos_fromlookup_nearest_dumb(repeatedeos,tabledimen, degentable, whichtablevar, whichtablesubtype, iffun, whichindep, quant1, vartypearraylocal, indexarrayvar,answers,badlookups);
  }
  else{
    // linear
    failreturn=get_eos_fromlookup_linear(repeatedeos,tabledimen,degentable, whichtablevar, whichtablesubtype, iffun, whichindep, quant1, vartypearraylocal, indexarrayvar,answers,badlookups);
  }
#elif(0)
  // nearest dumb for all
  failreturn=get_eos_fromlookup_nearest_dumb(repeatedeos,tabledimen,degentable, whichtablevar, whichtablesubtype, iffun, whichindep, quant1, vartypearraylocal, indexarrayvar,answers,badlookups);
#elif(0)
  // nearest
  failreturn=get_eos_fromlookup_nearest(repeatedeos,tabledimen,degentable, whichtablevar, whichtablesubtype, iffun, whichindep, quant1, vartypearraylocal, indexarrayvar,answers,badlookups);

#elif(0)

  // linear
  failreturn=get_eos_fromlookup_linear(repeatedeos,tabledimen,degentable, whichtablevar, whichtablesubtype, iffun, whichindep, quant1, vartypearraylocal, indexarrayvar,answers,badlookups);

#elif(1)

  // parabolic for one quantity, but linear for rest
  failreturn=get_eos_fromlookup_parabolicone(repeatedeos,tabledimen,degentable, whichtablevar, whichtablesubtype, iffun, whichindep, quant1, vartypearraylocal, indexarrayvar,answers,badlookups);

#elif(0)
  // parabolic for all quantities
  failreturn=get_eos_fromlookup_parabolicfull(repeatedeos,tabledimen,degentable, whichtablevar, whichtablesubtype, iffun, whichindep, quant1, vartypearraylocal, indexarrayvar,answers,badlookups);
#endif

  // return whether really bad failure
  return(failreturn);

}







// full parabolic interpolation
// Uses globals so can make them thread safe instead of using static's that are not
// Assumes H is not a dependent dimension as for whichdatatype==4
static int get_eos_fromlookup_parabolicfull(int repeatedeos, int tabledimen, int degentable, int whichtablevar, int whichtablesubtype, int *iffun, int whichindep, FTYPE quant1, int *vartypearraylocal, FTYPEEOS *indexarrayvar, FTYPEEOS *answers, int *badlookups)
{

  INDEXPARAMETERSLOCALDEFS;

  FTYPEEOS tempcheck;
  FTYPEEOS tempanswers[MAXEOSPIPELINE];
  int shouldloginterp[MAXEOSPIPELINE];
  int localbadlookups[MAXEOSPIPELINE];

  FTYPEEOS totalf[MAXEOSPIPELINE][3][3][3][3]; // 3 values for parabolic interpolation
  FTYPEEOS (*tfptr)[3][3][3][3];

  int includef[MAXEOSPIPELINE][3][3][3][3]; // 3 values for parabolic interpolation
  int (*includeptr)[3][3][3][3];

  int include2f[MAXEOSPIPELINE][3][3][3]; // 3 values for parabolic interpolation
  int (*include2ptr)[3][3][3];

  int include3f[MAXEOSPIPELINE][3][3]; // 3 values for parabolic interpolation
  int (*include3ptr)[3][3];

  int include4f[MAXEOSPIPELINE][3]; // 3 values for parabolic interpolation
  int (*include4ptr)[3];

  int iii,jjj,kkk,lll,mmm;
  FTYPEEOS ieos,jeos,keos,leos,meos;

  int whichdegenfun;
  FTYPEEOS xmx0,AA,BB;

  int numcols,coli;


  

  ///////////////
  //
  // get number of columns possible for this whichtablesubtype
  //
  ///////////////
  numcols = get_numcols(whichtablevar,whichtablesubtype);




  ///////////////
  //
  // store index within table
  //
  ///////////////
  ieos=indexarrayvar[1];
  jeos=indexarrayvar[2];
  keos=indexarrayvar[3];
  leos=indexarrayvar[4];
  meos=0;


  ///////////////
  //
  // pointer shift
  //
  ///////////////
  tfptr=(FTYPEEOS (*)[3][3][3][3]) (&(totalf[0][1][1][1][1])); // so tfptr[-1,0,1]
  includeptr=(int (*)[3][3][3][3]) (&(includef[0][1][1][1][1])); // so includeptr[-1,0,1]
  include2ptr=(int (*)[3][3][3]) (&(includef[0][1][1][1])); // so include2ptr[-1,0,1]
  include3ptr=(int (*)[3][3]) (&(includef[0][1][1])); // so include3ptr[-1,0,1]
  include4ptr=(int (*)[3]) (&(includef[0][1])); // so include4ptr[-1,0,1]


  ///////////
  //
  // definition consistent with numerical assignments of indices of arrays
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
  if(1||repeatedeos==0){
    // assume jeos, keos, leos set correctly for given degentable and tabledimen (i.e. are 0 or integers in those cases so d's computed correctly)
    kazii=ROUND2INT(ieos); //ii=(int)ieos; // round instead when doing parabolic interpolation with 3 points
    // limit to within table
    if(kazii<1) kazii=1; // 0 is minimum
    if(kazii>tablesize[whichtablevar][vartypearraylocal[1]]-2) kazii=tablesize[whichtablevar][vartypearraylocal[1]]-2; // N-1 is maximum

    kazjj=ROUND2INT(jeos); //kazjj=(int)jeos; // round instead when doing parabolic interpolation with 3 points
    // limit to within table
    if(degentable==0){
      if(kazjj<1) kazjj=1; // 0 is minimum
      if(kazjj>tablesize[whichtablevar][whichindep]-2) kazjj=tablesize[whichtablevar][whichindep]-2; // N-1 is maximum
    }

    kazkk=ROUND2INT(keos); //kazkk=(int)keos; // round instead when doing parabolic interpolation with 3 points
    // limit to within table
    if(kazkk<1) kazkk=1; // 0 is minimum
    if(kazkk>tablesize[whichtablevar][vartypearraylocal[3]]-2) kazkk=tablesize[whichtablevar][vartypearraylocal[3]]-2; // N-1 is maximum

    kazll=ROUND2INT(leos); //kazll=(int)leos; // round instead when doing parabolic interpolation with 3 points
    // limit to within table
    if(kazll<1) kazll=1; // 0 is minimum
    if(kazll>tablesize[whichtablevar][vartypearraylocal[4]]-2) kazll=tablesize[whichtablevar][vartypearraylocal[4]]-2; // N-1 is maximum

    kazmm=0;

    // if only choosing 1 value, then choose rounded version
    kaziio=ROUND2INT(ieos);
    kazjjo=ROUND2INT(jeos);
    kazkko=ROUND2INT(keos);
    kazllo=ROUND2INT(leos);
    kazmmo=0;
    

    // set range of loops for different table types
    // tabledimen overrules table type (i.e. take section out of fuller table -- assumed to be 0 index)
    // GODMARK: these things couuld be stored as functions of whichtablevar/degentable/tabledimen
    if(tablesize[whichtablevar][vartypearraylocal[1]]==1) { kazstartiii=kazendiii=0; } // kazll will be 0
    else if(tablesize[whichtablevar][vartypearraylocal[1]]==2) { kazstartiii=0; kazendiii=1;} // kazll will be 0
    else if(tablesize[whichtablevar][vartypearraylocal[1]]>2) { kazstartiii=-1; kazendiii=1;}

    if(degentable==1 || tablesize[whichtablevar][vartypearraylocal[2]]==1) {kazstartjjj=kazendjjj=0;}
    else if(tablesize[whichtablevar][vartypearraylocal[2]]==2) { kazstartjjj=0; kazendjjj=1;}
    else if(tablesize[whichtablevar][vartypearraylocal[2]]>2) { kazstartjjj=-1; kazendjjj=1;}

    if(tablesize[whichtablevar][vartypearraylocal[3]]==1) {kazstartkkk=kazendkkk=0;}
    else if(tablesize[whichtablevar][vartypearraylocal[3]]==2){ kazstartkkk=0; kazendkkk=1;}
    else if(tablesize[whichtablevar][vartypearraylocal[3]]>2){ kazstartkkk=-1; kazendkkk=1;}
     
    if(tablesize[whichtablevar][vartypearraylocal[4]]==1) {kazstartlll=kazendlll=0;}
    else if(tablesize[whichtablevar][vartypearraylocal[4]]==2){ kazstartlll=0; kazendlll=1;}
    else if(tablesize[whichtablevar][vartypearraylocal[4]]>2){ kazstartlll=-1; kazendlll=1;}

    if(tablesize[whichtablevar][vartypearraylocal[5]]==1) {kazstartmmm=kazendmmm=0;}
    else if(tablesize[whichtablevar][vartypearraylocal[5]]==2){ kazstartmmm=0; kazendmmm=1;}
    else if(tablesize[whichtablevar][vartypearraylocal[5]]>2){ kazstartmmm=-1; kazendmmm=1;}

    // avoid using data beyond table by shifting starting point and essentially extrapolating
    if(kazii+kazendiii>=tablesize[whichtablevar][vartypearraylocal[1]]){ kazii=tablesize[whichtablevar][vartypearraylocal[1]]-1-kazendiii; }
    if(kazjj+kazendjjj>=tablesize[whichtablevar][vartypearraylocal[2]]){ kazjj=tablesize[whichtablevar][vartypearraylocal[2]]-1-kazendjjj; }
    if(kazkk+kazendkkk>=tablesize[whichtablevar][vartypearraylocal[3]]){ kazkk=tablesize[whichtablevar][vartypearraylocal[3]]-1-kazendkkk; }
    if(kazll+kazendlll>=tablesize[whichtablevar][vartypearraylocal[4]]){ kazll=tablesize[whichtablevar][vartypearraylocal[4]]-1-kazendlll; }
    if(kazmm+kazendmmm>=tablesize[whichtablevar][vartypearraylocal[5]]){ kazmm=tablesize[whichtablevar][vartypearraylocal[5]]-1-kazendmmm; }


  }

  // overrides (must be placed outside if(repeatedeos==0) conditional so always done since kazendlll,kazendmmm, etc. are stored per whichd and not also per whichtablesubtype
  // only extra table is function of Ynu or H if whichdatatype==4
  if(WHICHDATATYPEGENERAL==4 && (whichtablesubtype!=SUBTYPEEXTRA)){kazll=kazmm=kazllo=kazmmo=leos=meos=kazstartlll=kazendlll=kazstartmmm=kazendmmm=0;}
  if(whichtablesubtype==SUBTYPEDEGEN){kazjj=kazjjo=jeos=kazstartjjj=kazendjjj=0;}




#define LOOPKAZIJKL for(lll=kazstartlll;lll<=kazendlll;lll++)for(kkk=kazstartkkk;kkk<=kazendkkk;kkk++)for(jjj=kazstartjjj;jjj<=kazendjjj;jjj++)for(iii=kazstartiii;iii<=kazendiii;iii++)


  // Loop over nearby table values and determine bi-linearly interpolated value
  // 4-D means 2^4=16 positions
  // 2-D means 2^2=4 positions
  // get 3 values as function of density
  mmm=0;
  LOOPKAZIJKL{

    
#if(CHECKIFVALIDEOSDATA)
    if(degentable==0){
      // don't use values of table that have no inversion to temperature
      get_arrays_eostable_direct_temperature(whichdegenfun,whichtablevar,kazmm+mmm,kazll+lll,kazkk+kkk,kazjj+jjj,kazii+iii,&tempcheck);
    }
    if(degentable==1 || tempcheck>invalidtempcode) // Avoid invalid inversions if T>Tbad, but only deal with temperature if degentable==0
#else
      if(1)
#endif
        {

          // get value
          get_arrays_eostable_direct(whichdegenfun,degentable,whichtablevar,whichtablesubtype,iffun,kazmm+mmm,kazll+lll,kazkk+kkk,kazjj+jjj,kazii+iii,tempanswers);


          // non-pipelined section
          for(coli=0;coli<numcols;coli++){
            if(iffun[coli]){

              if(shouldloginterp[coli] && (tempanswers[coli]<=0.0 && DOPRELOGIFY==0 ||tempanswers[coli]>=0.99999*OUTOFBOUNDSPRELOGIFY && DOPRELOGIFY==1) ){
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
#if(DOLOGINTERP && DOPRELOGIFY)
      // then above check for logification already sufficient
#else 
      if(shouldloginterp[coli]){
        LOOPKAZIJKL{
          if(includeptr[coli][iii][jjj][kkk][lll]){
            if(tfptr[coli][iii][jjj][kkk][lll]>0.0) tfptr[coli][iii][jjj][kkk][lll] = log10(tfptr[coli][iii][jjj][kkk][lll]);
            else includeptr[coli][iii][jjj][kkk][lll]=0; // then force to not include
          }
        }// end loop over dimensions
      }// end if logifying
#endif

    }// end iffun==1
  }// end over coli





#define KAZPARALOOP1 for(lll=kazstartlll;lll<=kazendlll;lll++)for(kkk=kazstartkkk;kkk<=kazendkkk;kkk++)for(jjj=kazstartjjj;jjj<=kazendjjj;jjj++)
#define KAZPARALOOP2 for(lll=kazstartlll;lll<=kazendlll;lll++)for(kkk=kazstartkkk;kkk<=kazendkkk;kkk++)
#define KAZPARALOOP3 for(lll=kazstartlll;lll<=kazendlll;lll++)
#define KAZPARALOOP4 if(1)

  //////////////
  //
  // perform interpolation over values
  //
  //////////////
  for(coli=0;coli<numcols;coli++){
    if(iffun[coli]){
      localbadlookups[coli]=0; // default is to succeed
  
      KAZPARALOOP1{


        // now use 3 data points to get density-parabolic distribution and value at ieos
        // have tfptr[coli][iii][jjj][kkk][lll] @ iii=-1,0,1 with i=0 meaning ii and offset being ROUND2INT(ieos)
        // Form parabolic answer
        include2ptr[coli][jjj][kkk][lll]=1; // default
        if(kazstartiii<=-1 && kazendiii>=1 && includeptr[coli][-1][jjj][kkk][lll] && includeptr[coli][0][jjj][kkk][lll]&& includeptr[coli][1][jjj][kkk][lll]){
          xmx0 = (ieos-(FTYPEEOS)kazii);
          AA = 0.5*(tfptr[coli][1][jjj][kkk][lll]-tfptr[coli][-1][jjj][kkk][lll]);
          BB = 0.5*(tfptr[coli][1][jjj][kkk][lll]+tfptr[coli][-1][jjj][kkk][lll]-2.0*tfptr[coli][0][jjj][kkk][lll]);
          tfptr[coli][0][jjj][kkk][lll] = tfptr[coli][0][jjj][kkk][lll] + AA*xmx0 + BB*xmx0*xmx0;
          // e.g. xmx0=0 -> tf[0]
          // e.g. xmx0=-1 -> tf[0]-tf[1]/2+tf[-1]/2 + tf[1]/2+tf[-1]/2-tf[0] = tf[-1]
          // e.g. xmx0=1 -> tf[0]+tf[1]/2-tf[-1]/2 + tf[1]/2+tf[-1]/2-tf[0] = tf[1]
        }
        else if(kazstartiii<=-1 && kazendiii>=0 && includeptr[coli][-1][jjj][kkk][lll] && includeptr[coli][0][jjj][kkk][lll]){
          // reduce to linear interpolation
          xmx0 = (ieos-(FTYPEEOS)(kazii-1));
          AA=(tfptr[coli][0][jjj][kkk][lll]-tfptr[coli][-1][jjj][kkk][lll]);
          tfptr[coli][0][jjj][kkk][lll] = tfptr[coli][-1][jjj][kkk][lll] + AA*xmx0;
          // e.g. xmx0=0 (ieos=kazii-1) -> tf[-1]
          // e.g. xmx0=1 (ieos=kazii) -> tf[0]
        }
        else if(kazstartiii<=0 && kazendiii>=1 && includeptr[coli][0][jjj][kkk][lll] && includeptr[coli][1][jjj][kkk][lll]){
          // reduce to linear interpolation
          xmx0 = (ieos-(FTYPEEOS)(kazii));
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
          xmx0 = (jeos-(FTYPEEOS)kazjj);
          AA = 0.5*(tfptr[coli][0][1][kkk][lll]-tfptr[coli][0][-1][kkk][lll]);
          BB = 0.5*(tfptr[coli][0][1][kkk][lll]+tfptr[coli][0][-1][kkk][lll]-2.0*tfptr[coli][0][0][kkk][lll]);
          tfptr[coli][0][0][kkk][lll] = tfptr[coli][0][0][kkk][lll] + AA*xmx0 + BB*xmx0*xmx0;
        }
        else if(kazstartjjj<=-1 && kazendjjj>=0 && include2ptr[coli][-1][kkk][lll] && include2ptr[coli][0][kkk][lll]){
          // reduce to linear interpolation
          xmx0 = (jeos-(FTYPEEOS)(kazjj-1));
          AA=(tfptr[coli][0][0][kkk][lll]-tfptr[coli][0][-1][kkk][lll]);
          tfptr[coli][0][0][kkk][lll] = tfptr[coli][0][-1][kkk][lll] + AA*xmx0;
        }
        else if(kazstartjjj<=0 && kazendjjj>=1 && include2ptr[coli][0][kkk][lll] && include2ptr[coli][1][kkk][lll]){
          // reduce to linear interpolation
          xmx0 = (jeos-(FTYPEEOS)(kazjj));
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
          xmx0 = (keos-(FTYPEEOS)kazkk);
          AA = 0.5*(tfptr[coli][0][0][1][lll]-tfptr[coli][0][0][-1][lll]);
          BB = 0.5*(tfptr[coli][0][0][1][lll]+tfptr[coli][0][0][-1][lll]-2.0*tfptr[coli][0][0][0][lll]);
          tfptr[coli][0][0][0][lll] = tfptr[coli][0][0][0][lll] + AA*xmx0 + BB*xmx0*xmx0;
        }
        else if(kazstartkkk<=-1 && kazendkkk>=0 && include3ptr[coli][-1][lll] && include3ptr[coli][0][lll]){
          // reduce to linear interpolation
          xmx0 = (keos-(FTYPEEOS)(kazkk-1));
          AA=(tfptr[coli][0][0][0][lll]-tfptr[coli][0][0][-1][lll]);
          tfptr[coli][0][0][0][lll] = tfptr[coli][0][0][-1][lll] + AA*xmx0;
        }
        else if(kazstartkkk<=0 && kazendkkk>=1 && include3ptr[coli][0][lll] && include3ptr[coli][1][lll]){
          // reduce to linear interpolation
          xmx0 = (keos-(FTYPEEOS)(kazkk));
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
          xmx0 = (leos-(FTYPEEOS)kazll);
          AA = 0.5*(tfptr[coli][0][0][0][1]-tfptr[coli][0][0][0][-1]);
          BB = 0.5*(tfptr[coli][0][0][0][1]+tfptr[coli][0][0][0][-1]-2.0*tfptr[coli][0][0][0][0]);
          tfptr[coli][0][0][0][0] = tfptr[coli][0][0][0][0] + AA*xmx0 + BB*xmx0*xmx0;
        }
        else if(kazstartlll<=-1 && kazendlll>=0 && include4ptr[coli][-1] && include4ptr[coli][0]){
          // reduce to linear interpolation
          xmx0 = (leos-(FTYPEEOS)(kazll-1));
          AA=(tfptr[coli][0][0][0][0]-tfptr[coli][0][0][0][-1]);
          tfptr[coli][0][0][0][0] = tfptr[coli][0][0][0][-1] + AA*xmx0;
        }
        else if(kazstartlll<=0 && kazendlll>=1 && include4ptr[coli][0] && include4ptr[coli][1]){
          // reduce to linear interpolation
          xmx0 = (leos-(FTYPEEOS)(kazll));
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
          // no valid data points, so return badlookup
          localbadlookups[coli]=1;
          // GODMARK: Should have caught this with the check if within table function, so could return failure
          //   return(1);

        }// end else if no valid points

#if(PRODUCTION==0)
        if(!isfinite(tfptr[coli][0][0][0][0])){
          dualfprintf(fail_file,"4notfinite, :: %d %d %d\n",0,0,0);
        }
#endif

      }// end KAZPARALOOP4




      if(localbadlookups[coli]==0){

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
        FTYPEEOS largest=-BIG;
        FTYPEEOS smallest=BIG;
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
        // unlogify (if either logified here or prelogified)
        //
        //////////////////////////
        if(shouldloginterp[coli]) answers[coli]=pow(10.0,answers[coli]);

        //////////////////////////
        //
        // invert offset
        //
        //////////////////////////
        if(degentable==1) offsetquant2_general_inverse(whichdegenfun, quant1, answers[coli], &answers[coli]);
        badlookups[coli]=0; // finally return that have good answer

      } // end if local was not bad lookup

    }// end iffun==1
  }// end over coli
  


  //////////////////////////////
  //
  // return that no failure
  //
  //////////////////////////////
  return(0);

}










// linear for all but one quantity
static int get_eos_fromlookup_parabolicone(int repeatedeos, int tabledimen, int degentable, int whichtablevar, int whichtablesubtype, int *iffun, int whichindep, FTYPE quant1, int *vartypearraylocal, FTYPEEOS *indexarrayvar, FTYPEEOS *answers, int *badlookups)
{

  INDEXPARAMETERSLOCALDEFS;

  FTYPEEOS totaldist[MAXEOSPIPELINE][3];
  FTYPEEOS (*tdist)[3];

  FTYPEEOS totalanswers[MAXEOSPIPELINE][3]; // 3 values for parabolic interpolation
  FTYPEEOS (*tfptr)[3];

  FTYPEEOS tempanswers[MAXEOSPIPELINE];
  int shouldloginterp[MAXEOSPIPELINE];

  int localiffun[MAXEOSPIPELINE]; // need because we do multiple lookups here
  int localbadlookups[MAXEOSPIPELINE];

  FTYPEEOS dist;
  FTYPEEOS tempcheck;
  int whichdegenfun;


  int iii,jjj,kkk,lll,mmm;
  FTYPEEOS ieos,jeos,keos,leos,meos;
  int genpi,genstart,genend;

  FTYPEEOS xmx0,AA,BB;
  int numcols,coli;

  int paravar;


  // get number of columns possible for this whichtablesubtype
  numcols = get_numcols(whichtablevar,whichtablesubtype);



  /////////////////////
  //
  // semi-user choice for which variable to use parabolic interpolation upon
  //
  // Note: If choose paravar==2, can't do so for whichdegen==1, so choose different paravar for that case
  //
  /////////////////////
  if(whichtablesubtype==SUBTYPEEXTRA){
    paravar=3; // Y_e=3
  }
  else{
    paravar=3; // Y_e=3 // also do Y_e for other tables
  }
  

  ///////////
  //
  // Store grid position within EOS table
  //
  ///////////
  ieos=indexarrayvar[1];
  jeos=indexarrayvar[2];
  keos=indexarrayvar[3];
  leos=indexarrayvar[4];
  meos=indexarrayvar[5];


  ///////////
  //
  // Pointer shift
  //
  ///////////
  tfptr=(FTYPEEOS (*)[3]) (&(totalanswers[0][1])); // so tfptr[-1,0,1][]
  tdist=(FTYPEEOS (*)[3]) (&(totaldist[0][1]));

  ///////////
  //
  // definition consistent with numerical assignments of indeces of arrays
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
  if(1||repeatedeos==0){
    // assume jeos, keos, leos set correctly for given degentable and tabledimen (i.e. are 0 or integers in those cases so d's computed correctly)

    ///////////////////////////
    //
    // Set nearest integer index to use
    //
    // Also set range of loops for different table types
    // tabledimen overrules table type (i.e. take section out of fuller table -- assumed to be 0 index)
    // GODMARK: these things couuld be stored as functions of whichtablevar/degentable/tabledimen
    //
    //////////////////////
    //
    ///////////////////////////

    // 1
    if(paravar==1){
      kazii=ROUND2INT(ieos); //kazii=(int)ieos; // round instead when doing parabolic interpolation with 3 points
      // limit to within table
      if(kazii<1) kazii=1; // 0 is minimum
      if(kazii>tablesize[whichtablevar][vartypearraylocal[1]]-2) kazii=tablesize[whichtablevar][vartypearraylocal[1]]-2; // N-1 is maximum

      // set range of loops for different table types
      // tabledimen overrules table type (i.e. take section out of fuller table -- assumed to be 0 index)
      // GODMARK: these things couuld be stored as functions of whichtablevar/degentable/tabledimen
      if(tablesize[whichtablevar][vartypearraylocal[1]]==1) { kazstartiii=kazendiii=0; } // kazll will be 0
      else if(tablesize[whichtablevar][vartypearraylocal[1]]==2) { kazstartiii=0; kazendiii=1;} // kazll will be 0
      else if(tablesize[whichtablevar][vartypearraylocal[1]]>2) { kazstartiii=-1; kazendiii=1;}
    }
    else{
      kazii=(int)ieos;

      kazstartiii=0;
      if(tablesize[whichtablevar][vartypearraylocal[1]]!=1) kazendiii=1;
      else kazendiii=0;
    }

    // 2
    if(paravar==2){
      kazjj=ROUND2INT(jeos);
      // limit to within table
      if(kazjj<1) kazjj=1; // 0 is minimum
      // was whichindep for vartypearraylocal[2]
      if(kazjj>tablesize[whichtablevar][vartypearraylocal[2]]-2) kazjj=tablesize[whichtablevar][vartypearraylocal[2]]-2; // N-1 is maximum

      if(degentable==1 || tablesize[whichtablevar][vartypearraylocal[2]]==1) {kazstartjjj=kazendjjj=0;}
      else if(tablesize[whichtablevar][vartypearraylocal[2]]==2) { kazstartjjj=0; kazendjjj=1;}
      else if(tablesize[whichtablevar][vartypearraylocal[2]]>2) { kazstartjjj=-1; kazendjjj=1;}
    }
    else{
      kazjj=(int)jeos;

      kazstartjjj=0;
      if(degentable==0 && tablesize[whichtablevar][vartypearraylocal[2]]!=1) kazendjjj=1;
      else kazendjjj=0;
    }

    // 3
    if(paravar==3){
      kazkk=ROUND2INT(keos);
      // limit to within table
      if(kazkk<1) kazkk=1; // 0 is minimum
      if(kazkk>tablesize[whichtablevar][vartypearraylocal[3]]-2) kazkk=tablesize[whichtablevar][vartypearraylocal[3]]-2; // N-1 is maximum

      if(tablesize[whichtablevar][vartypearraylocal[3]]==1) {kazstartkkk=kazendkkk=0;}
      else if(tablesize[whichtablevar][vartypearraylocal[3]]==2){ kazstartkkk=0; kazendkkk=1;}
      else if(tablesize[whichtablevar][vartypearraylocal[3]]>2){ kazstartkkk=-1; kazendkkk=1;}
    }
    else{
      kazkk=(int)keos;

      kazstartkkk=0;
      if(tablesize[whichtablevar][vartypearraylocal[3]]!=1) kazendkkk=1;
      else kazendkkk=0;
    }

    // 4
    if(paravar==4){
      kazll=ROUND2INT(leos);
      // limit to within table
      if(kazll<1) kazll=1; // 0 is minimum
      if(kazll>tablesize[whichtablevar][vartypearraylocal[4]]-2) kazll=tablesize[whichtablevar][vartypearraylocal[4]]-2; // N-1 is maximum

      if(tablesize[whichtablevar][vartypearraylocal[4]]==1) {kazstartlll=kazendlll=0;}
      else if(tablesize[whichtablevar][vartypearraylocal[4]]==2){ kazstartlll=0; kazendlll=1;}
      else if(tablesize[whichtablevar][vartypearraylocal[4]]>2){ kazstartlll=-1; kazendlll=1;}
    }
    else{
      kazll=(int)leos;

      kazstartlll=0;
      if(tablesize[whichtablevar][vartypearraylocal[4]]!=1) kazendlll=1;
      else kazendlll=0;
    }

    // 5
    if(paravar==5){
      kazmm=ROUND2INT(meos);
      // limit to within table
      if(kazmm<1) kazmm=1; // 0 is minimum
      if(kazmm>tablesize[whichtablevar][vartypearraylocal[5]]-2) kazmm=tablesize[whichtablevar][vartypearraylocal[5]]-2; // N-1 is maximum

      if(tablesize[whichtablevar][vartypearraylocal[5]]==1) {kazstartmmm=kazendmmm=0;}
      else if(tablesize[whichtablevar][vartypearraylocal[5]]==2){ kazstartmmm=0; kazendmmm=1;}
      else if(tablesize[whichtablevar][vartypearraylocal[5]]>2){ kazstartmmm=-1; kazendmmm=1;}
    }
    else{
      kazmm=(int)meos;

      kazstartmmm=0;
      if(tablesize[whichtablevar][vartypearraylocal[5]]!=1) kazendmmm=1;
      else kazendmmm=0;
    }


    

    /////////////////////
    //
    // avoid using data beyond table by shifting starting point and essentially extrapolating
    //
    /////////////////////
    if(kazii+kazendiii>=tablesize[whichtablevar][vartypearraylocal[1]]){ kazii=tablesize[whichtablevar][vartypearraylocal[1]]-1-kazendiii; }
    if(kazjj+kazendjjj>=tablesize[whichtablevar][vartypearraylocal[2]]){ kazjj=tablesize[whichtablevar][vartypearraylocal[2]]-1-kazendjjj; }
    if(kazkk+kazendkkk>=tablesize[whichtablevar][vartypearraylocal[3]]){ kazkk=tablesize[whichtablevar][vartypearraylocal[3]]-1-kazendkkk; }
    if(kazll+kazendlll>=tablesize[whichtablevar][vartypearraylocal[4]]){ kazll=tablesize[whichtablevar][vartypearraylocal[4]]-1-kazendlll; }
    if(kazmm+kazendmmm>=tablesize[whichtablevar][vartypearraylocal[5]]){ kazmm=tablesize[whichtablevar][vartypearraylocal[5]]-1-kazendmmm; }

    kazdi[1]=ieos-(FTYPEEOS)kazii;
    kazdi[0]=1.0-kazdi[1];
    kazdj[1]=jeos-(FTYPEEOS)kazjj;
    kazdj[0]=1.0-kazdj[1];
    kazdk[1]=keos-(FTYPEEOS)kazkk;
    kazdk[0]=1.0-kazdk[1];
    kazdl[1]=leos-(FTYPEEOS)kazll;
    kazdl[0]=1.0-kazdl[1];
    kazdm[1]=meos-(FTYPEEOS)kazmm;
    kazdm[0]=1.0-kazdm[1];
    
  }


  // overrides (must be placed outside if(repeatedeos==0) conditional so always done since kazendlll,kazendmmm, etc. are stored per whichd and not also per whichtablesubtype
  // only extra table is function of Ynu or H if whichdatatype==4
  if(WHICHDATATYPEGENERAL==4 && (whichtablesubtype!=SUBTYPEEXTRA)){kazll=kazmm=kazllo=kazmmo=leos=meos=kazstartlll=kazendlll=kazstartmmm=kazendmmm=kazdm[1]=kazdl[1]=0; kazdl[0]=kazdm[0]=1.0;}
  if(whichtablesubtype==SUBTYPEDEGEN){kazjj=kazjjo=jeos=kazstartjjj=kazendjjj=kazdj[1]=0; kazdj[0]=1.0; }



  // pick para-interpolation range for quantity to do para interpolation upon
  if(paravar==1){      genstart=kazstartiii; genend=kazendiii;  xmx0 = kazdi[1];}
  else if(paravar==2){ genstart=kazstartjjj; genend=kazendjjj;  xmx0 = kazdj[1];}
  else if(paravar==3){ genstart=kazstartkkk; genend=kazendkkk;  xmx0 = kazdk[1];}
  else if(paravar==4){ genstart=kazstartlll; genend=kazendlll;  xmx0 = kazdl[1];}
  else if(paravar==5){ genstart=kazstartmmm; genend=kazendmmm;  xmx0 = kazdm[1];}
    



  ///////////////////////
  //
  // initialize cumulative result and distance for normalization
  //
  ///////////////////////
  for(coli=0;coli<numcols;coli++){
    if(iffun[coli]){
      for(genpi=genstart;genpi<=genend;genpi++){
        tfptr[coli][genpi]=0.0;
        tdist[coli][genpi]=0.0;
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

  // m,l,k,j, then "i" so that "i" is fastest as best for memory accesses
#define PARAONEKAZLOOP for(mmm=kazstartmmm;mmm<=kazendmmm;mmm++)for(lll=kazstartlll;lll<=kazendlll;lll++)for(kkk=kazstartkkk;kkk<=kazendkkk;kkk++)for(jjj=kazstartjjj;jjj<=kazendjjj;jjj++)for(iii=kazstartiii;iii<=kazendiii;iii++)



  PARAONEKAZLOOP{

    // pick which variable to do parabolic interpolation on
    if(paravar==1) genpi=iii;
    else if(paravar==2) genpi=jjj;
    else if(paravar==3) genpi=kkk;
    else if(paravar==4) genpi=lll;
    else if(paravar==5) genpi=mmm;



#if(CHECKIFVALIDEOSDATA)
    if(degentable==0){
      // don't use values of table that have no inversion to temperature
      get_arrays_eostable_direct_temperature(whichdegenfun,whichtablevar,kazmm+mmm,kazll+lll,kazkk+kkk,kazjj+jjj,kazii+iii,&tempcheck);
    }
    //    dualfprintf(fail_file,"temp=%21.15g\n",tempcheck);
    if(degentable==1 || tempcheck>invalidtempcode) // Avoid invalid inversions if T>Tbad, but only deal with temperature if degentable==0
#else
      if(1)
#endif
        {

          // get value
          get_arrays_eostable_direct(whichdegenfun,degentable,whichtablevar,whichtablesubtype,iffun,kazmm+mmm,kazll+lll,kazkk+kkk,kazjj+jjj,kazii+iii,tempanswers);

          // get distance
          if(paravar==1)      dist = kazdj[jjj]*kazdk[kkk]*kazdl[lll]*kazdm[mmm]; // dist for kazi
          else if(paravar==2) dist = kazdi[iii]*kazdk[kkk]*kazdl[lll]*kazdm[mmm]; // dist for kazj
          else if(paravar==3) dist = kazdi[iii]*kazdj[jjj]*kazdl[lll]*kazdm[mmm]; // dist for kazk
          else if(paravar==4) dist = kazdi[iii]*kazdj[jjj]*kazdk[kkk]*kazdm[mmm]; // dist for kazl
          else if(paravar==5) dist = kazdi[iii]*kazdj[jjj]*kazdk[kkk]*kazdl[lll]; // dist for kazm



          // non-pipelined section
          // first do offset and see if can use log for all points
          for(coli=0;coli<numcols;coli++){
            if(iffun[coli]){
              if(degentable==1) offsetquant2_general(whichdegenfun, quant1, tempanswers[coli], &tempanswers[coli]);
  
#if(DOLOGINTERP && DOPRELOGIFY)
              // only add if non-log or if should be log and not out of bounds
              if(shouldloginterp[coli]==0 || shouldloginterp[coli] && tempanswers[coli]<0.99999*OUTOFBOUNDSPRELOGIFY){
                tfptr[coli][genpi] += tempanswers[coli]*dist;
                tdist[coli][genpi] += dist;
              }
#else 
              if(shouldloginterp[coli]){
                if(tempanswers[coli]<=0.0){
                  // then don't include since if should have been able to do log but couldn't, then remove point from interpolation
                }
                else{
                  // should and can do log
                  tfptr[coli][genpi] += log10(tempanswers[coli])*dist; // interpolated answer that will replace answer[coli] in the end
                  tdist[coli][genpi] += dist;
                }
              }// if should do log interp
              else{
                // shouldn't do log
                tfptr[coli][genpi] += tempanswers[coli]*dist;
                tdist[coli][genpi] += dist;
              }// else if not log interp
#endif

            }// end if doing this coli
          }// end over coli      
        }// end if good temperature or doing degentable
      else{
        // no good temperature for this point
      }
  }// end loop over dimensions






  ////////////////////////////
  //
  // finally normalize
  //
  ////////////////////////////
  for(genpi=genstart;genpi<=genend;genpi++) {


    // GODMARK: Probably want to try reverting to linear instead of using _dumb() value

    // see if any bad lookups
    int badany=0;
    for(coli=0;coli<numcols;coli++){
      if(iffun[coli]){
        // check if any bad lookups
        if(tdist[coli][genpi]==0.0){
          badany++;
          localiffun[coli]=1;
          localbadlookups[coli]=1;
        }
        else{
          localiffun[coli]=0;
          localbadlookups[coli]=0; // already good
        }
      }
    }
      
  
    // never found good temperature, so look up nearest_dumb()
    // GODMARK: Could instead reduce from parabolic to lower order since using nearest_dumb() is probably a bad idea
    if(badany){
      get_eos_fromlookup_nearest_dumb(repeatedeos, tabledimen, degentable, whichtablevar, whichtablesubtype, localiffun, whichindep, quant1, vartypearraylocal, indexarrayvar, tempanswers, localbadlookups); // this can overwrite default localbadlookups from 0->1

      // do offset and logify, but only if badlookup before but goodlookup from nearest_dumb()
      for(coli=0;coli<numcols;coli++){
        if(tdist[coli][genpi]==0.0 && localiffun[coli] && localbadlookups[coli]==0){
          // only do if fixing no-temperature lookup
          if(degentable==1) offsetquant2_general(whichdegenfun, quant1, tempanswers[coli], &tempanswers[coli]);
   
#if(DOLOGINTERP && DOPRELOGIFY)
          // only add if non-log or if should be log and not out of bounds
          if(shouldloginterp[coli]==0 || shouldloginterp[coli] && tempanswers[coli]<0.99999*OUTOFBOUNDSPRELOGIFY) tfptr[coli][genpi] = tempanswers[coli];
#else 
          if(shouldloginterp[coli]) tfptr[coli][genpi] = log10(tempanswers[coli]);
          else tfptr[coli][genpi] = tempanswers[coli];
#endif
   
        }// end if localbadlookups==0
      }//end over coli
    }// end if any bad values that need to be looked-up with nearest_dumb()


    // for those cases when found good temperature, normalize based upon distance
    for(coli=0;coli<numcols;coli++){
      if(iffun[coli]){
        if(tdist[coli][genpi]!=0.0){ // only do if found good temperature (and can't be badlookup)
          //   tfptr[coli][genpi] /= tdist[coli][genpi];
          // already normalized properly (old wrong normalization)
        }
      } // end iffun==1
    }//end over coli

  } // end over genpi




  //////////////////////////////
  //
  // now use 3 data points to get parabolic distribution and value
  // have tfptr[genpi] @ genpi=-1,0,1 with i=0 meaning ii and offset being ROUND2INT(ieos)
  // Form parabolic answer
  //
  //////////////////////////////
  FTYPE minanswer,maxanswer;
  for(coli=0;coli<numcols;coli++){
    if(iffun[coli]){
      if(localbadlookups[coli]==0){

        AA = 0.5*(tfptr[coli][1]-tfptr[coli][-1]);
        BB = 0.5*(tfptr[coli][1]+tfptr[coli][-1]-2.0*tfptr[coli][0]);
        answers[coli] = tfptr[coli][0] + AA*xmx0 + BB*xmx0*xmx0;

        // enforce monotonicity to avoid exaggerated answers near sharp regions
        maxanswer=MAX(MAX(tfptr[coli][1],tfptr[coli][-1]),tfptr[coli][0]);
        minanswer=MIN(MIN(tfptr[coli][1],tfptr[coli][-1]),tfptr[coli][0]);
        if(answers[coli]>maxanswer) answers[coli]=maxanswer;
        if(answers[coli]<minanswer) answers[coli]=minanswer;
 
        // either prelogify or not, we have to unlogify if should have logified
        if(shouldloginterp[coli]) answers[coli]=pow(10.0,answers[coli]);
        if(degentable==1) offsetquant2_general_inverse(whichdegenfun, quant1, answers[coli], &answers[coli]);
        badlookups[coli]=0; // finally return that was good lookup

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
static int get_eos_fromlookup_linear(int repeatedeos, int tabledimen, int degentable, int whichtablevar, int whichtablesubtype, int *iffun, int whichindep, FTYPE quant1, int *vartypearraylocal, FTYPEEOS *indexarrayvar, FTYPEEOS *answers, int *badlookups)
{

  INDEXPARAMETERSLOCALDEFS;

  FTYPEEOS dist,totaldist[MAXEOSPIPELINE];
  FTYPEEOS tempanswers[MAXEOSPIPELINE];
  FTYPEEOS totalanswers[MAXEOSPIPELINE];
  FTYPEEOS tempcheck;
  int iii,jjj,kkk,lll,mmm;
  int whichdegenfun;
  int shouldloginterp[MAXEOSPIPELINE];
  FTYPEEOS ieos,jeos,keos,leos,meos;
  int qi;
  int numcols,coli;



  // get number of columns possible for this whichtablesubtype
  numcols = get_numcols(whichtablevar,whichtablesubtype);


  ieos=indexarrayvar[1];
  jeos=indexarrayvar[2];
  keos=indexarrayvar[3];
  leos=indexarrayvar[4];
  meos=indexarrayvar[5];
  



  // definition consistent with numerical assignments of indices of arrays
  whichdegenfun = whichindep-1;




  // determine if should do log interpolation
  get_dologinterp_subtype_wrapper(degentable, whichtablesubtype, numcols, shouldloginterp);


  if(1||repeatedeos==0){
    // assume jeos, keos, leos set correctly for given degentable and tabledimen (i.e. are 0 or integers in those cases so d's computed correctly)
    kazii=(int)ieos;
    //if(degentable) kazjj=0; else kazjj=(int)jeos; // force to be 0
    kazjj=(int)jeos; // (int)(0.0) = 0, so ok as is
    kazkk=(int)keos;
    kazll=(int)leos;
    kazmm=(int)meos;

    // set range of loops for different table types
    // tabledimen overrules table type (i.e. take section out of fuller table -- assumed to be 0 index)
    // GODMARK: these things couuld be stored as functions of whichtablevar/degentable/tabledimen
    if(tablesize[whichtablevar][vartypearraylocal[1]]!=1) kazendiii=1;
    else kazendiii=0;

    if(degentable==0 && tablesize[whichtablevar][vartypearraylocal[2]]!=1) kazendjjj=1;
    else kazendjjj=0;

    if(tablesize[whichtablevar][vartypearraylocal[3]]!=1) kazendkkk=1;
    else kazendkkk=0;
    
    if(tablesize[whichtablevar][vartypearraylocal[4]]!=1) kazendlll=1;
    else kazendlll=0;

    if(tablesize[whichtablevar][vartypearraylocal[5]]!=1) kazendmmm=1;
    else kazendmmm=0;


    // avoid using data beyond table by shifting starting point and essentially extrapolating
    if(kazii+kazendiii>=tablesize[whichtablevar][vartypearraylocal[1]]){ kazii=tablesize[whichtablevar][vartypearraylocal[1]]-1-kazendiii; }
    if(kazjj+kazendjjj>=tablesize[whichtablevar][vartypearraylocal[2]]){ kazjj=tablesize[whichtablevar][vartypearraylocal[2]]-1-kazendjjj; }
    if(kazkk+kazendkkk>=tablesize[whichtablevar][vartypearraylocal[3]]){ kazkk=tablesize[whichtablevar][vartypearraylocal[3]]-1-kazendkkk; }
    if(kazll+kazendlll>=tablesize[whichtablevar][vartypearraylocal[4]]){ kazll=tablesize[whichtablevar][vartypearraylocal[4]]-1-kazendlll; }
    if(kazmm+kazendmmm>=tablesize[whichtablevar][vartypearraylocal[5]]){ kazmm=tablesize[whichtablevar][vartypearraylocal[5]]-1-kazendmmm; }


 
    kazdi[1]=ieos-(FTYPEEOS)kazii;
    kazdi[0]=1.0-kazdi[1];
    kazdj[1]=jeos-(FTYPEEOS)kazjj;
    kazdj[0]=1.0-kazdj[1];
    kazdk[1]=keos-(FTYPEEOS)kazkk;
    kazdk[0]=1.0-kazdk[1];
    kazdl[1]=leos-(FTYPEEOS)kazll;
    kazdl[0]=1.0-kazdl[1];
    kazdm[1]=meos-(FTYPEEOS)kazmm;
    kazdm[0]=1.0-kazdm[1];
    
  }

  // overrides (must be placed outside if(repeatedeos==0) conditional so always done since kazendlll,kazendmmm, etc. are stored per whichd and not also per whichtablesubtype
  // only extra table is function of Ynu or H if whichdatatype==4
  if(WHICHDATATYPEGENERAL==4 && (whichtablesubtype!=SUBTYPEEXTRA)){leos=meos=kazll=kazmm=kazdl[1]=kazdm[1]=kazendlll=kazendmmm=0; kazdl[0]=kazdm[0]=1.0;}
  if(whichtablesubtype==SUBTYPEDEGEN){jeos=kazjj=kazdj[1]=kazendjjj=0; kazdj[0]=1.0;}


    



    
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



    //      // DEBUG:
    //      for(coli=0;coli<numcols;coli++){
    // if(iffun[coli]){
    //         if(whichtablesubtype==SUBTYPEDEGEN || whichtablesubtype==SUBTYPESTANDARD && coli==0) dualfprintf(fail_file,"coli=%d :: ieos=%21.15g jeos=%21.15g keos=%21.15g leos=%21.15g meos=%21.15g :: kazii=%d kazjj=%d kazkk=%d kazll=%d kazmm=%d : kazdi[%d]=%21.15g kazdj[%d]=%21.15g kazdk[%d]=%21.15g kazdl[%d]=%21.15g kazdm[%d]=%21.15g\n",coli,ieos,jeos,keos,leos,meos,kazii,kazjj,kazkk,kazll,kazmm,iii,kazdi[iii],jjj,kazdj[jjj],kkk,kazdk[kkk],lll,kazdl[lll],mmm,kazdm[mmm]);
    // }
    //      }




#if(CHECKIFVALIDEOSDATA)
    if(degentable==0){
      // don't use values of table that have no inversion to temperature
      get_arrays_eostable_direct_temperature(whichdegenfun,whichtablevar,kazmm+mmm,kazll+lll,kazkk+kkk,kazjj+jjj,kazii+iii,&tempcheck);
    }
    if(degentable==1 || tempcheck>invalidtempcode) // Avoid invalid inversions if T>Tbad, but only deal with temperature if degentable==0
#else
      if(1)
#endif
        {
          // get value
          get_arrays_eostable_direct(whichdegenfun,degentable,whichtablevar,whichtablesubtype,iffun,kazmm+mmm,kazll+lll,kazkk+kkk,kazjj+jjj,kazii+iii,tempanswers);


          // get distance to value
          dist = kazdi[iii]*kazdj[jjj]*kazdk[kkk]*kazdl[lll]*kazdm[mmm];



          //      // DEBUG:
          //      for(coli=0;coli<numcols;coli++){
          // if(iffun[coli]){
          //         if(whichtablesubtype==SUBTYPEDEGEN || whichtablesubtype==SUBTYPESTANDARD && coli==0) dualfprintf(fail_file,"coli=%d :: dist=%21.15g : answer=%21.15g tempcheck=%21.15g\n",coli,dist,tempanswers[coli],tempcheck);
          //    }
          //      }



          // non-pipelined section:
          for(coli=0;coli<numcols;coli++){

            if(iffun[coli]){

              if(degentable==1) offsetquant2_general(whichdegenfun, quant1, tempanswers[coli], &tempanswers[coli]);

#if(DOLOGINTERP && DOPRELOGIFY)
              // only add if non-log or if should be log and not out of bounds
              if(shouldloginterp[coli]==0 || shouldloginterp[coli] && tempanswers[coli]<0.99999*OUTOFBOUNDSPRELOGIFY){
                totalanswers[coli] += tempanswers[coli]*dist;
                totaldist[coli] += dist;
              }
#else 
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
#endif

            }// end if doing this coli
          }// end over coli      

        }// end if good temperature or doing degentable
      else{
      
        // out of bounds for temperature
        dualfprintf(fail_file,"Out of bounds based upon temperature: whichtablesubtype=%d tempcheck=%21.15g\n",whichtablesubtype,tempcheck);
      
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

        if(debugfail>=2){
          dualfprintf(fail_file,"No valid data points in table despite within table.\n");
          dualfprintf(fail_file,"kaz=%d %d %d %d : %d %d : %d %d : %d %d : %d %d : %d %d\n",kazii,kazjj,kazkk,kazll,kazstartiii,kazendiii,kazstartjjj,kazendjjj,kazstartkkk,kazendkkk,kazstartlll,kazendlll,kazstartmmm,kazendmmm);
          dualfprintf(fail_file,"%d %d %d %d %d X %d %21.15g X X X X\n",repeatedeos, tabledimen, degentable, whichtablevar, whichtablesubtype, whichindep, quant1);

          dualfprintf(fail_file,"coli=%d tempanswers=%21.15g\n",coli,tempanswers[coli]);

        }

      }
      else{

        // get normalized answer
        // answers[coli] = totalanswers[coli]/totaldist[coli];
        // true normalization for multid-linear interpolation
        answers[coli] = totalanswers[coli];

        // invert log interpolation
        // either prelogify or not, we have to unlogify if should have logified
        if(shouldloginterp[coli]) answers[coli]=pow(10.0,answers[coli]);
 
        // invert offset
        if(degentable==1) offsetquant2_general_inverse(whichdegenfun, quant1, answers[coli], &answers[coli]);
 
        // finally return that was good lookup
        badlookups[coli]=0;

      }// end else if totaldist!=0.0
    }// end iffun==1
  }// end over coli


  // DEBUG:
  //  if(whichtablesubtype==SUBTYPEEXTRA){
  //    dualfprintf(fail_file,"kaz=%d %d %d %d %d : %d %d : %d %d : %d %d : %d %d : %d %d\n",kazii,kazjj,kazkk,kazll,kazmm,kazstartiii,kazendiii,kazstartjjj,kazendjjj,kazstartkkk,kazendkkk,kazstartlll,kazendlll,kazstartmmm,kazendmmm);
  //    dualfprintf(fail_file,"%d %d %d %d %d X %d %21.15g X X X X\n",repeatedeos, tabledimen, degentable, whichtablevar, whichtablesubtype, whichindep, quant1);
  //  }


  ///////////////////
  //
  // return that no bad failure
  //
  ///////////////////
  return(0);

}










// nearest neighbor interpolation but caution to temperature defined or not
static int get_eos_fromlookup_nearest(int repeatedeos, int tabledimen, int degentable, int whichtablevar, int whichtablesubtype, int *iffun, int whichindep, FTYPE quant1, int *vartypearraylocal, FTYPEEOS *indexarrayvar, FTYPEEOS *answers, int *badlookups)
{

  INDEXPARAMETERSLOCALDEFS;

  FTYPEEOS totaldist;
  FTYPEEOS tempcheck;
  int whichdegenfun;
  int iii,jjj,kkk,lll,mmm;
  FTYPEEOS ieos,jeos,keos,leos,meos;
  int shouldloginterp[MAXEOSPIPELINE];
  int localbadlookups[MAXEOSPIPELINE];
  int numcols,coli;
  int gotgood;

  // get number of columns possible for this whichtablesubtype
  numcols = get_numcols(whichtablevar,whichtablesubtype);


  ieos=indexarrayvar[1];
  jeos=indexarrayvar[2];
  keos=indexarrayvar[3];
  leos=indexarrayvar[4];
  meos=indexarrayvar[5];




  // definition consistent with numerical assignments of indices of arrays
  whichdegenfun = whichindep-1;


  // determine if should do log interpolation
  get_dologinterp_subtype_wrapper(degentable, whichtablesubtype, numcols, shouldloginterp);



  if(1||repeatedeos==0){
    // assume jeos, keos, leos set correctly for given degentable and tabledimen (i.e. are 0 or integers in those cases so d's computed correctly)
    // if only choosing 1 value, then choose rounded version
    kazii=ROUND2INT(ieos);
    kazjj=ROUND2INT(jeos);
    kazkk=ROUND2INT(keos);
    kazll=ROUND2INT(leos);
    kazmm=ROUND2INT(meos);
 

  }

  
  // overrides (must be placed outside if(repeatedeos==0) conditional so always done since kazendlll,kazendmmm, etc. are stored per whichd and not also per whichtablesubtype
  // only extra table is function of Ynu or H if whichdatatype==4
  if(WHICHDATATYPEGENERAL==4 && (whichtablesubtype!=SUBTYPEEXTRA)){leos=meos=kazll=kazmm=kazendlll=kazendmmm=0;}
  if(whichtablesubtype==SUBTYPEDEGEN){jeos=kazjj=kazendjjj=0;}



  ////////////
  //
  // lookup value 
  // try ROUND2INT version first
  //
  ////////////
  // default is that bad lookup
  for(coli=0;coli<numcols;coli++) localbadlookups[coli]=1;

#if(CHECKIFVALIDEOSDATA)
  if(degentable==0){
    // don't use values of table that have no inversion to temperature
    get_arrays_eostable_direct_temperature(whichdegenfun,whichtablevar,kazmm,kazll,kazkk,kazjj,kazii,&tempcheck);
  }
  if(degentable==1 || tempcheck>invalidtempcode) // Avoid invalid inversions if T>Tbad, but only deal with temperature if degentable==0
#else
    if(1)
#endif
      {
        get_arrays_eostable_direct(whichdegenfun,degentable,whichtablevar,whichtablesubtype,iffun,kazmm,kazll,kazkk,kazjj,kazii,answers);


        for(coli=0;coli<numcols;coli++) {
          if(iffun[coli]){
#if(DOLOGINTERP && DOPRELOGIFY)
            // only add if non-log or if should be log and not out of bounds
            if(shouldloginterp[coli]==0 || shouldloginterp[coli] && answers[coli]<0.99999*OUTOFBOUNDSPRELOGIFY) localbadlookups[coli]=0;
#else
            if(shouldloginterp[coli]==0 || shouldloginterp[coli] && answers[coli]>0.0) localbadlookups[coli]=0;
#endif
          }// end iffun==1
        }// end over coli
      }


  ////////////////////////
  // check if all good
  ////////////////////////
  gotgood=1; for(coli=0;coli<numcols;coli++) if(iffun[coli]) if(localbadlookups[coli]==1) gotgood=0;



  ////////////
  //
  // if ROUND2INT looked-up in invalid part of EOS domain, try first-found neighbor
  //
  ////////////

  if(gotgood==0){


    if(1||repeatedeos==0){
      
      kazii=(int)(ieos);
      kazjj=(int)(jeos);
      kazkk=(int)(keos);
      kazll=(int)(leos);
      kazmm=(int)(meos);
 

      // set range of loops for different table types
      // tabledimen overrules table type (i.e. take section out of fuller table -- assumed to be 0 index)
      // GODMARK: these things couuld be stored as functions of whichtablevar/degentable/tabledimen
      if(tablesize[whichtablevar][vartypearraylocal[1]]!=1) kazendiii=1;
      else kazendiii=0;

      if(degentable==0 && tablesize[whichtablevar][vartypearraylocal[2]]!=1) kazendjjj=1;
      else kazendjjj=0;

      if(tablesize[whichtablevar][vartypearraylocal[3]]!=1) kazendkkk=1;
      else kazendkkk=0;
    
      if(tablesize[whichtablevar][vartypearraylocal[4]]!=1) kazendlll=1;
      else kazendlll=0;

      if(tablesize[whichtablevar][vartypearraylocal[5]]!=1) kazendmmm=1;
      else kazendmmm=0;

      // avoid using data beyond table by shifting starting point and essentially extrapolating
      if(kazii+kazendiii>=tablesize[whichtablevar][vartypearraylocal[1]]){ kazii=tablesize[whichtablevar][vartypearraylocal[1]]-1-kazendiii; }
      if(kazjj+kazendjjj>=tablesize[whichtablevar][vartypearraylocal[2]]){ kazjj=tablesize[whichtablevar][vartypearraylocal[2]]-1-kazendjjj; }
      if(kazkk+kazendkkk>=tablesize[whichtablevar][vartypearraylocal[3]]){ kazkk=tablesize[whichtablevar][vartypearraylocal[3]]-1-kazendkkk; }
      if(kazll+kazendlll>=tablesize[whichtablevar][vartypearraylocal[4]]){ kazll=tablesize[whichtablevar][vartypearraylocal[4]]-1-kazendlll; }
      if(kazmm+kazendmmm>=tablesize[whichtablevar][vartypearraylocal[5]]){ kazmm=tablesize[whichtablevar][vartypearraylocal[5]]-1-kazendmmm; }

    }


#define LOOPKAZIJKLM for(mmm=0;mmm<=kazendmmm;mmm++)for(lll=0;lll<=kazendlll;lll++)for(kkk=0;kkk<=kazendkkk;kkk++)for(jjj=0;jjj<=kazendjjj;jjj++)for(iii=0;iii<=kazendiii;iii++)

    // loop over to ensure obtained best nearest neighbor    
    LOOPKAZIJKLM{
    
#if(CHECKIFVALIDEOSDATA)
      if(degentable==0){
        // don't use values of table that have no inversion to temperature
        get_arrays_eostable_direct_temperature(whichdegenfun,whichtablevar,kazmm,kazll,kazkk,kazjj,kazii,&tempcheck);
      }
      if(degentable==1 || tempcheck>invalidtempcode) // Avoid invalid inversions if T>Tbad, but only deal with temperature if degentable==0
#else
        if(1)
#endif
          {
            get_arrays_eostable_direct(whichdegenfun,degentable,whichtablevar,whichtablesubtype,iffun,kazmm,kazll,kazkk,kazjj,kazii,answers);

            for(coli=0;coli<numcols;coli++){
              if(iffun[coli]){
#if(DOLOGINTERP && DOPRELOGIFY)
                // only add if non-log or if should be log and not out of bounds
                if(shouldloginterp[coli]==0 || shouldloginterp[coli] && answers[coli]<0.99999*OUTOFBOUNDSPRELOGIFY) localbadlookups[coli]=0;
#else
                if(shouldloginterp[coli]==0 || shouldloginterp[coli] && answers[coli]>0.0) localbadlookups[coli]=0;
#endif
              }// end iffun==1
            }// end over coli

            ////////////////////////
            // see if should keep going and trying to get all or some values
            // check if all good
            ////////////////////////
            gotgood=1; for(coli=0;coli<numcols;coli++) if(iffun[coli]) if(localbadlookups[coli]==1) gotgood=0;

            if(gotgood){
              // terminate all 4 loops once have all good values
              iii=kazendiii;
              jjj=kazendjjj;
              kkk=kazendkkk;
              lll=kazendlll;
              mmm=kazendmmm;
            }

          }// end if did find a good temperature or doing degentable
    }// end loop over dimensions
  }// else if ROUNT2INT failed so picking nearest neighbor  




#if(DOLOGINTERP && DOPRELOGIFY)
  // then still need to unlogify even if didn't logify here since no interpolation, but prelogify did occur still.
  for(coli=0;coli<numcols;coli++){
    if(shouldloginterp[coli] && iffun[coli] && localbadlookups[coli]==0){
      answers[coli]=pow(10.0,answers[coli]);
      // finall return that good lookup
      badlookups[coli]=0;
    }
  }
#endif


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
static int get_eos_fromlookup_nearest_dumb(int repeatedeos, int tabledimen, int degentable, int whichtablevar, int whichtablesubtype, int *iffun, int whichindep, FTYPE quant1, int *vartypearraylocal, FTYPEEOS *indexarrayvar, FTYPEEOS *answers, int *badlookups)
{

  INDEXPARAMETERSLOCALDEFS;

  int whichdegenfun;
  int shouldloginterp[MAXEOSPIPELINE];
  int localbadlookups[MAXEOSPIPELINE];
  int numcols,coli;


  // get number of columns possible for this whichtablesubtype
  numcols = get_numcols(whichtablevar,whichtablesubtype);

  // definition consistent with numerical assignments of indeces of arrays
  whichdegenfun = whichindep-1;


  // determine if should do log interpolation
  get_dologinterp_subtype_wrapper(degentable, whichtablesubtype, numcols, shouldloginterp);


  if(1||repeatedeos==0){

#if(0)
    kazii=(int)indexarrayvar[1];
    kazjj=(int)indexarrayvar[2];
    kazkk=(int)indexarrayvar[3];
    kazll=(int)indexarrayvar[4];
    kazmm=(int)indexarrayvar[5];
#else
    // if only choosing 1 value, then choose rounded version
    kazii=ROUND2INT(indexarrayvar[1]);
    kazjj=ROUND2INT(indexarrayvar[2]);
    kazkk=ROUND2INT(indexarrayvar[3]);
    kazll=ROUND2INT(indexarrayvar[4]);
    kazmm=ROUND2INT(indexarrayvar[5]);
#endif


    if(kazii>=tablesize[whichtablevar][vartypearraylocal[1]]) kazii=tablesize[whichtablevar][vartypearraylocal[1]]-1;
    if(kazjj>=tablesize[whichtablevar][vartypearraylocal[2]]) kazjj=tablesize[whichtablevar][vartypearraylocal[2]]-1;
    if(kazkk>=tablesize[whichtablevar][vartypearraylocal[3]]) kazkk=tablesize[whichtablevar][vartypearraylocal[3]]-1;
    if(kazll>=tablesize[whichtablevar][vartypearraylocal[4]]) kazll=tablesize[whichtablevar][vartypearraylocal[4]]-1;
    if(kazmm>=tablesize[whichtablevar][vartypearraylocal[5]]) kazmm=tablesize[whichtablevar][vartypearraylocal[5]]-1;
  }

  // overrides (must be placed outside if(repeatedeos==0) conditional so always done since kazendlll,kazendmmm, etc. are stored per whichd and not also per whichtablesubtype
  // only extra table is function of Ynu or H if whichdatatype==4
  if(WHICHDATATYPEGENERAL==4 && (whichtablesubtype!=SUBTYPEEXTRA)){kazll=kazmm=kazendlll=kazendmmm=0;}
  if(whichtablesubtype==SUBTYPEDEGEN){kazjj=kazendjjj=0;}


  // get value
  get_arrays_eostable_direct(whichdegenfun,degentable,whichtablevar,whichtablesubtype,iffun,kazmm,kazll,kazkk,kazjj,kazii,answers);


  for(coli=0;coli<numcols;coli++) localbadlookups[coli]=1; // default to bad

  // return whether was bad lookup
#if(DOLOGINTERP && DOPRELOGIFY)
  // only add if non-log or if should be log and not out of bounds
  for(coli=0;coli<numcols;coli++) if(shouldloginterp[coli]==0 || shouldloginterp[coli] && answers[coli]<0.99999*OUTOFBOUNDSPRELOGIFY) localbadlookups[coli]=0;
#else
  for(coli=0;coli<numcols;coli++) if(shouldloginterp[coli]==0 || shouldloginterp[coli] && answers[coli]>0.0) localbadlookups[coli]=0;
#endif


#if(DOLOGINTERP && DOPRELOGIFY)
  // then still need to unlogify even if didn't logify here since no interpolation, but prelogify did occur still.
  for(coli=0;coli<numcols;coli++) if(shouldloginterp[coli] && localbadlookups[coli]==0)  answers[coli]=pow(10.0,answers[coli]);
#endif

  // return that was good answer
  for(coli=0;coli<numcols;coli++) if(localbadlookups[coli]==0) badlookups[coli]=0;



  // DEBUG:
  //  if(whichtablesubtype==SUBTYPESTANDARD||whichtablesubtype==SUBTYPEDEGEN){
  //    for(coli=0;coli<numcols;coli++) if(iffun[coli]) dualfprintf(fail_file,"kazii=%d kazjj=%d kazkk=%d kazll=%d kazmm=%d whichtablesubtype=%d answers[coli=%d]=%21.15g bad=%d\n",kazii,kazjj,kazkk,kazll,kazmm,whichtablesubtype,coli,answers[coli],badlookups[coli]);
  //    dualfprintf(fail_file,"%d %d %d %d %d %d %21.15g\n",repeatedeos, tabledimen, degentable, whichtablevar, whichtablesubtype, whichindep, quant1);
  //  }
  
  
  // return that no bad failure
  return(0);

}





//#define thiseostable(name,fun,index5,index4,index3,index2,index1) EOSMAC(name,fun,index5,index4,index3,index2,index1)
//#define index2array(name,fun) name[fun][indexarrayvar[5]][indexarrayvar[4]][indexarrayvar[3]][indexarrayvar[2]][indexarrayvar[1]]


// kinda like get_arrays_eostable() but directly grab multiple columns of data from a given table subtype
// simple and simplezoom table code is just copy-paste from full table with fulltable->simpletable and fulltable->simplezoomtable
static void get_arrays_eostable_direct(int whichd, int whichdegen, int whichtablevar, int whichtablesubtype, int *iffun, int mmm, int lll, int kkk, int jjj, int iii, FTYPEEOS *values)
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
  numcols = get_numcols(whichtablevar,whichtablesubtype);
  //  whichdlocal = whichdintablesubtype[whichtablesubtype]; // could ensure that whichd==whichdlocal for non-degen and non-temp




  if(0){
  }
#if(ALLOWFULLTABLE==1)
  else if(whichtablevar==FULLTABLE||whichtablevar==FULLTABLEEXTRA){
    if(whichdegen==ISNOTDEGENTABLE){
      // degentable not dealt with here
      if(whichtablevar==FULLTABLE){
        // assume whichtablesubtype==SUBTYPEEXTRA is not reached if WHICHDATATYPEGENERAL==4 (controlled by which_eostable())
        if(whichtablesubtype==SUBTYPESTANDARD){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltablestandard,0,mmm,lll,kkk,jjj,iii,coli);}
        else if(whichtablesubtype==SUBTYPEGUESS){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltableguess,0,mmm,lll,kkk,jjj,iii,coli);}
        else if(whichtablesubtype==SUBTYPEDISS){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltablediss,0,mmm,lll,kkk,jjj,iii,coli);}
        else if(whichtablesubtype==SUBTYPEDP){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltabledp,0,mmm,lll,kkk,jjj,iii,coli);}
        else if(whichtablesubtype==SUBTYPESDEN){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltablesden,0,mmm,lll,kkk,jjj,iii,coli);}
        else if(whichtablesubtype==SUBTYPESSPEC){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltablesspec,0,mmm,lll,kkk,jjj,iii,coli);}
        else if(whichtablesubtype==SUBTYPEPOFCHI){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltablepofchi,0,mmm,lll,kkk,jjj,iii,coli);}
        // temp table special like degen table such that only 1 returned temp associated with whichd value
        else if(whichtablesubtype==SUBTYPETEMP){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltabletemp,whichd,mmm,lll,kkk,jjj,iii,coli);} // whichd used
        else if(whichtablesubtype==SUBTYPEEXTRA){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltableextra,0,mmm,lll,kkk,jjj,iii,coli);}
        // in principle can request degen data directly, so allow this
        else if(whichtablesubtype==SUBTYPEDEGEN){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltabledegen,whichd,mmm,lll,kkk,jjj,iii,coli);} // whichd used
      }
      else if(whichtablevar==FULLTABLEEXTRA){
        // reach here if WHICHDATATYPEGENERAL==4 and want extras(controlled by which_eostable())
        if(whichtablesubtype==SUBTYPEEXTRA){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltableextra,0,mmm,lll,kkk,jjj,iii,coli);}
        // in principle can request degen data directly, so allow this
        else if(whichtablesubtype==SUBTYPEDEGEN){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltableextradegen,whichd,mmm,lll,kkk,jjj,iii,coli);} // whichd used
        else if(whichtablesubtype==SUBTYPETEMP){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltableextratemp,whichd,mmm,lll,kkk,jjj,iii,coli);}
      }
    }
    else{
      // then ignore whichtablesubtype
      if(whichtablevar==FULLTABLE){
        for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltabledegen,whichd,mmm,lll,kkk,jjj,iii,coli); // whichd used
      }
      else if(whichtablevar==FULLTABLEEXTRA){
        for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eosfulltableextradegen,whichd,mmm,lll,kkk,jjj,iii,coli); // whichd used
      }
    }// end else if degentable      
  }// end if fulltable
#endif
#if(ALLOWSIMPLETABLE==1)
  else if(whichtablevar==SIMPLETABLE||whichtablevar==SIMPLETABLEEXTRA){
    if(whichdegen==ISNOTDEGENTABLE){
      // degentable not dealt with here
      if(whichtablevar==SIMPLETABLE){
        // assume whichtablesubtype==SUBTYPEEXTRA is not reached if WHICHDATATYPEGENERAL==4 (controlled by which_eostable())
        if(whichtablesubtype==SUBTYPESTANDARD){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletablestandard,0,mmm,lll,kkk,jjj,iii,coli);}
        else if(whichtablesubtype==SUBTYPEGUESS){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletableguess,0,mmm,lll,kkk,jjj,iii,coli);}
        else if(whichtablesubtype==SUBTYPEDISS){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletablediss,0,mmm,lll,kkk,jjj,iii,coli);}
        else if(whichtablesubtype==SUBTYPEDP){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletabledp,0,mmm,lll,kkk,jjj,iii,coli);}
        else if(whichtablesubtype==SUBTYPESDEN){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletablesden,0,mmm,lll,kkk,jjj,iii,coli);}
        else if(whichtablesubtype==SUBTYPESSPEC){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletablesspec,0,mmm,lll,kkk,jjj,iii,coli);}
        else if(whichtablesubtype==SUBTYPEPOFCHI){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletablepofchi,0,mmm,lll,kkk,jjj,iii,coli);}
        // temp table special like degen table such that only 1 returned temp associated with whichd value
        else if(whichtablesubtype==SUBTYPETEMP){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletabletemp,whichd,mmm,lll,kkk,jjj,iii,coli);} // whichd used
        else if(whichtablesubtype==SUBTYPEEXTRA){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletableextra,0,mmm,lll,kkk,jjj,iii,coli);}
        // in principle can request degen data directly, so allow this
        else if(whichtablesubtype==SUBTYPEDEGEN){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletabledegen,whichd,mmm,lll,kkk,jjj,iii,coli);} // whichd used
      }
      else if(whichtablevar==SIMPLETABLEEXTRA){
        // reach here if WHICHDATATYPEGENERAL==4 and want extras(controlled by which_eostable())
        if(whichtablesubtype==SUBTYPEEXTRA){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletableextra,0,mmm,lll,kkk,jjj,iii,coli);}
        // in principle can request degen data directly, so allow this
        else if(whichtablesubtype==SUBTYPEDEGEN){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletableextradegen,whichd,mmm,lll,kkk,jjj,iii,coli);} // whichd used
        else if(whichtablesubtype==SUBTYPETEMP){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletableextratemp,whichd,mmm,lll,kkk,jjj,iii,coli);}
      }
    }
    else{
      // then ignore whichtablesubtype
      if(whichtablevar==SIMPLETABLE){
        for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletabledegen,whichd,mmm,lll,kkk,jjj,iii,coli); // whichd used
      }
      else if(whichtablevar==SIMPLETABLEEXTRA){
        for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimpletableextradegen,whichd,mmm,lll,kkk,jjj,iii,coli); // whichd used
      }
    }// end else if degentable      
  }// end if simpletable
#endif
#if(ALLOWSIMPLEZOOMTABLE==1)
  else if(whichtablevar==SIMPLEZOOMTABLE||whichtablevar==SIMPLEZOOMEXTRATABLE){
    if(whichdegen==ISNOTDEGENTABLE){
      // degentable not dealt with here
      if(whichtablevar==SIMPLEZOOMTABLE){
        // assume whichtablesubtype==SUBTYPEEXTRA is not reached if WHICHDATATYPEGENERAL==4 (controlled by which_eostable())
        if(whichtablesubtype==SUBTYPESTANDARD){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtablestandard,0,mmm,lll,kkk,jjj,iii,coli);}
        else if(whichtablesubtype==SUBTYPEGUESS){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtableguess,0,mmm,lll,kkk,jjj,iii,coli);}
        else if(whichtablesubtype==SUBTYPEDISS){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtablediss,0,mmm,lll,kkk,jjj,iii,coli);}
        else if(whichtablesubtype==SUBTYPEDP){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtabledp,0,mmm,lll,kkk,jjj,iii,coli);}
        else if(whichtablesubtype==SUBTYPESDEN){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtablesden,0,mmm,lll,kkk,jjj,iii,coli);}
        else if(whichtablesubtype==SUBTYPESSPEC){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtablesspec,0,mmm,lll,kkk,jjj,iii,coli);}
        else if(whichtablesubtype==SUBTYPEPOFCHI){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtablepofchi,0,mmm,lll,kkk,jjj,iii,coli);}
        // temp table special like degen table such that only 1 returned temp associated with whichd value
        else if(whichtablesubtype==SUBTYPETEMP){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtabletemp,whichd,mmm,lll,kkk,jjj,iii,coli);} // whichd used
        else if(whichtablesubtype==SUBTYPEEXTRA){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtableextra,0,mmm,lll,kkk,jjj,iii,coli);}
        // in principle can request degen data directly, so allow this
        else if(whichtablesubtype==SUBTYPEDEGEN){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtabledegen,whichd,mmm,lll,kkk,jjj,iii,coli);} // whichd used
      }
      else if(whichtablevar==SIMPLEZOOMEXTRATABLE){
        // reach here if WHICHDATATYPEGENERAL==4 and want extras(controlled by which_eostable())
        if(whichtablesubtype==SUBTYPEEXTRA){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtableextra,0,mmm,lll,kkk,jjj,iii,coli);}
        // in principle can request degen data directly, so allow this
        else if(whichtablesubtype==SUBTYPEDEGEN){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtableextradegen,whichd,mmm,lll,kkk,jjj,iii,coli);} // whichd used
        else if(whichtablesubtype==SUBTYPETEMP){ for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtableextratemp,whichd,mmm,lll,kkk,jjj,iii,coli);}
      }
    }
    else{
      // then ignore whichtablesubtype
      if(whichtablevar==SIMPLEZOOMTABLE){
        for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtabledegen,whichd,mmm,lll,kkk,jjj,iii,coli); // whichd used
      }
      else if(whichtablevar==SIMPLEZOOMEXTRATABLE){
        for(coli=0;coli<numcols;coli++) if(iffun[coli]) values[coli]=EOSMAC(eossimplezoomtableextradegen,whichd,mmm,lll,kkk,jjj,iii,coli); // whichd used
      }
    }// end else if degentable      
  }// end if simplezoomtable
#endif



}






// like get_arrays_eostable_direct(), but only for temperature that ignores associated variable actually looking up
// memory accessed using "coli" not "whichfun" which keeps code in kazfulleos_set_arrays.c simpler (before I included FIRSTEOS??? offsets)
static void get_arrays_eostable_direct_temperature(int whichd, int whichtablevar, int mmm, int lll, int kkk, int jjj, int iii, FTYPEEOS *values)
{
  int numcols,coli;


  // overrides based upon true memory allocated as only required
  // override instead of pre-condition because this function is called inside loop over other tables
  if(WHICHDATATYPEGENERAL==4){
    // Then temp has no Ynu or H dependence
    mmm=lll=0;
  }


  // get number of columns possible for this whichtablesubtype
  numcols = get_numcols(whichtablevar,SUBTYPETEMP);


  if(0){
  }
#if(ALLOWFULLTABLE==1)
  else if(whichtablevar==FULLTABLE){
    for(coli=0;coli<numcols;coli++) values[coli]=EOSMAC(eosfulltabletemp,whichd,mmm,lll,kkk,jjj,iii,coli);
  }// end if fulltable
  else if(whichtablevar==FULLTABLEEXTRA){
    for(coli=0;coli<numcols;coli++) values[coli]=EOSMAC(eosfulltableextratemp,whichd,mmm,lll,kkk,jjj,iii,coli);
  }// end if fulltableextra
#endif
  // simple table is just copy-paste from full table with fulltable->simpletable
#if(ALLOWSIMPLETABLE==1)
  else if(whichtablevar==SIMPLETABLE){
    for(coli=0;coli<numcols;coli++) values[coli]=EOSMAC(eossimpletabletemp,whichd,mmm,lll,kkk,jjj,iii,coli);
  }// end if simpletable
  else if(whichtablevar==SIMPLETABLEEXTRA){
    for(coli=0;coli<numcols;coli++) values[coli]=EOSMAC(eossimpletableextratemp,whichd,mmm,lll,kkk,jjj,iii,coli);
  }// end if simpletableextra
#endif
  // simplezoom table is just copy-paste from full table with fulltable->simplezoomtable
#if(ALLOWSIMPLEZOOMTABLE==1)
  else if(whichtablevar==SIMPLEZOOMTABLE){
    for(coli=0;coli<numcols;coli++) values[coli]=EOSMAC(eossimplezoomtabletemp,whichd,mmm,lll,kkk,jjj,iii,coli);
  }// end if simplezoomtable
  else if(whichtablevar==SIMPLEZOOMEXTRATABLE){
    for(coli=0;coli<numcols;coli++) values[coli]=EOSMAC(eossimplezoomtableextratemp,whichd,mmm,lll,kkk,jjj,iii,coli);
  }// end if simplezoomextratable
#endif


}














