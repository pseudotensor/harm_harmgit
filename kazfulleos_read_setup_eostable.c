/////////////////////////
//
// initialization and read-in of tables
// accesses global EOS arrays but only once for entire simulation
//
/////////////////////////

static int isextratable(int i);
static void assigntablecolumnnumbers(int i);
static void setvartypeinfo(void);
static void settablesubtypeinfo(void);
static void settablesizeexpected(int tablesizeexpected[][NUMEOSINDEPS]);


static void set_arrays_eostable(int whichdegen, int whichtable, int mmm, int lll, int kkk, int jjj, int iii, int incol, double value);
static void get_arrays_eostable(int whichdegen, int whichtable, int mmm, int lll, int kkk, int jjj, int iii, int incol, double *value);

static int get_dologinterp_subtype(int whichtablesubtype, int coli);
static int get_dologinterp_subtype_wrapper(int degentable, int whichtablesubtype, int numcols, int *shouldloginterp);

static void bcast_kazeos(void);




// sets up number of quantities to read-in from file
// numeosdegenquantitiesinfile[i] used for reading in degen data that will be stored
// NUMEOSDEGENQUANTITIESMEM1 : UofUdiffout, PofPdiffout, CHIofCHIdiffout, SofSdiffout read-in to be checked, not stored
// extralimits[i][0/1] used for interpolation
//
// Sets following globals:
//    numeosquantitiesinfile[i]
//    numeosdegenquantitiesinfile[i]
//    numallquantitiesinfile[i]
//    numalldegenquantitiesinfile[i]
//    extralimits[i][0,1]
//
//    numeosquantitiesinstandardlist[i]
//    numeosdegenquantitiesinstandardlist[i]
//
//

// before storing, reorder to be like standard "in" format for EOS quantities:
// Should be ONLY part of code that refers to the "inextra" versions of input macro labels except for use of "NUMEOSQUANTITIESBASEinextra" for total number of such quantities to read-in from file
void reorder_extratype(int isextratype, int numeosquantitiesinfile, FTYPE *values)
{
  int i,j;


  if(isextratype){
    // then translate "inextra" positions to "in" positions so rest of code can be the same
    FTYPE oldvalues[MAXEOSPIPELINE];

    for(i=0;i<numeosquantitiesinfile;i++){
      oldvalues[i]=values[i];
    }
    
    // copy over temperatures
    for(i=TEMPUinextra;i<=TEMPSinextra;i++){
      j = i-TEMPUinextra + TEMPUin;
      values[j] = oldvalues[i];
    }

    // copy over extras
    for(i=FIRSTEXTRAinextra;i<=LASTEXTRAinextra;i++){
      // this copies up to all possible extras, which is excessive for older whichdatatype's, but still within sizes of arrays
      j = i-FIRSTEXTRAinextra + FIRSTEXTRAin;
      values[j] = oldvalues[i];
    }    

  }
  else{
    // then already in canonical order
  }

  // now values is in canonical order so only have to refer to "in" instead of "inextra"

}


static int isextratable(int i)
{
  int extratype;

  if(i==FULLTABLE || i==SIMPLETABLE || i==SIMPLEZOOMTABLE) extratype=0;
  else if(i==EXTRAFULLTABLE || i==EXTRASIMPLETABLE || i==EXTRASIMPLEZOOMTABLE) extratype=1;

  return(extratype);
}


static void assigntablecolumnnumbers(int i)
{
  int extratype;

  // get whether extra type or not
  extratype=isextratable(i);


  if( ((i==FULLTABLE || i==EXTRAFULLTABLE) && ALLOWFULLTABLE) || ((i==SIMPLETABLE || i==EXTRASIMPLETABLE) && ALLOWSIMPLETABLE) || ((i==SIMPLEZOOMTABLE || i==EXTRASIMPLEZOOMTABLE) && ALLOWSIMPLEZOOMTABLE) ){

    if(extratype==0 && WHICHDATATYPEGENERAL!=4){
      // THEN FULL SET OF QUANTITIES IN TABLE

      // stored-normal quantities:
      numeosquantitiesinfile[i]=NUMEOSQUANTITIESBASEin;
      // can't split with other data types
      if(WHICHDATATYPEGENERAL==4) numeosquantitiesinfile[i]+=NUMEXTRAEOSQUANTITIESTYPE4;
      else if(WHICHDATATYPEGENERAL==3) numeosquantitiesinfile[i]+=NUMEXTRAEOSQUANTITIESTYPE3;
      else if(WHICHDATATYPEGENERAL==2) numeosquantitiesinfile[i]+=NUMEXTRAEOSQUANTITIESTYPE2;
      else if(WHICHDATATYPEGENERAL==1) numeosquantitiesinfile[i]+=NUMEXTRAEOSQUANTITIESTYPE1;

      // stored degen quantities:
      numeosdegenquantitiesinfile[i]=NUMEOSDEGENQUANTITIESin;

      // non-stored + stored:
      numallquantitiesinfile[i]=NUMEOSQUANTITIESNOTSTOREDin + numeosquantitiesinfile[i];
      // non-stored + stored:
      numalldegenquantitiesinfile[i]=NUMEOSDEGENQUANTITIESNOTSTOREDin + numeosdegenquantitiesinfile[i];
      if(WHICHDATATYPEGENERAL==1){
	extralimits[i][0]=EXTRA1; extralimits[i][1]=DATATYPE1_EXTRAFINAL;
      }
      else if(WHICHDATATYPEGENERAL==2){
	extralimits[i][0]=EXTRA2; extralimits[i][1]=DATATYPE2_EXTRAFINAL;
      }
      else if(WHICHDATATYPEGENERAL==3){
	extralimits[i][0]=EXTRA3; extralimits[i][1]=DATATYPE3_EXTRAFINAL;
      }
      else if(WHICHDATATYPEGENERAL==4){
	extralimits[i][0]=EXTRA4; extralimits[i][1]=DATATYPE4_EXTRAFINAL;
      }
    }
    else if(extratype==0 && WHICHDATATYPEGENERAL==4){
      // NO EXTRAS IN TABLE THEN

      numeosquantitiesinfile[i]=NUMEOSQUANTITIESBASEin;
      // always same number
      numeosdegenquantitiesinfile[i]=NUMEOSDEGENQUANTITIESin;
      numallquantitiesinfile[i]=NUMEOSQUANTITIESNOTSTOREDin + numeosquantitiesinfile[i];
      numalldegenquantitiesinfile[i]=NUMEOSDEGENQUANTITIESNOTSTOREDin + numeosdegenquantitiesinfile[i];
      extralimits[i][0]=0; extralimits[i][1]=-1; // so never touches extras
    }
    else if(extratype==1 && WHICHDATATYPEGENERAL!=4){
      // THEN NOTHING IN EXTRA TABLE SINCE ALL IN FULL TABLE

      numeosquantitiesinfile[i]=0;
      numeosdegenquantitiesinfile[i]=0;
      numallquantitiesinfile[i]=0;
      numalldegenquantitiesinfile[i]=0;
      extralimits[i][0]=0; extralimits[i][1]=-1; // so never touches extras
    }
    else if(extratype==1 && WHICHDATATYPEGENERAL==4){
      // THEN ONLY EXTRAS (with temperature as required to validate lookup)
      numeosquantitiesinfile[i]=NUMEOSQUANTITIESBASEinextra; // still must have temperature and that is present in same file as extras
      // then all in extra table
      if(WHICHDATATYPEGENERAL==4) numeosquantitiesinfile[i]+=NUMEXTRAEOSQUANTITIESTYPE4;
      else if(WHICHDATATYPEGENERAL==3) numeosquantitiesinfile[i]+=NUMEXTRAEOSQUANTITIESTYPE3;
      else if(WHICHDATATYPEGENERAL==2) numeosquantitiesinfile[i]+=NUMEXTRAEOSQUANTITIESTYPE2;
      else if(WHICHDATATYPEGENERAL==1) numeosquantitiesinfile[i]+=NUMEXTRAEOSQUANTITIESTYPE1;

      numeosdegenquantitiesinfile[i]=NUMEOSDEGENQUANTITIESin;
      numallquantitiesinfile[i]=NUMEOSQUANTITIESNOTSTOREDin + numeosquantitiesinfile[i];
      numalldegenquantitiesinfile[i]=NUMEOSDEGENQUANTITIESNOTSTOREDin + numeosdegenquantitiesinfile[i];
      if(WHICHDATATYPEGENERAL==1){
	extralimits[i][0]=EXTRA1; extralimits[i][1]=DATATYPE1_EXTRAFINAL;
      }
      else if(WHICHDATATYPEGENERAL==2){
	extralimits[i][0]=EXTRA2; extralimits[i][1]=DATATYPE2_EXTRAFINAL;
      }
      else if(WHICHDATATYPEGENERAL==3){
	extralimits[i][0]=EXTRA3; extralimits[i][1]=DATATYPE3_EXTRAFINAL;
      }
      else if(WHICHDATATYPEGENERAL==4){
	extralimits[i][0]=EXTRA4; extralimits[i][1]=DATATYPE4_EXTRAFINAL;
      }
    }
  }
  else{
    // then nothing in table
    numeosquantitiesinfile[i]=0;
    numeosdegenquantitiesinfile[i]=0;
    numallquantitiesinfile[i]=0;
    numalldegenquantitiesinfile[i]=0;
    extralimits[i][0]=0; extralimits[i][1]=-1; // so never touches extras
  }


  ////////////////////
  // standard canonical list

  // stored-normal quantities:
  numeosquantitiesinstandardlist[i]=NUMEOSQUANTITIESBASEin;
  // can't split with other data types
  if(WHICHDATATYPEGENERAL==4) numeosquantitiesinstandardlist[i]+=NUMEXTRAEOSQUANTITIESTYPE4;
  else if(WHICHDATATYPEGENERAL==3) numeosquantitiesinstandardlist[i]+=NUMEXTRAEOSQUANTITIESTYPE3;
  else if(WHICHDATATYPEGENERAL==2) numeosquantitiesinstandardlist[i]+=NUMEXTRAEOSQUANTITIESTYPE2;
  else if(WHICHDATATYPEGENERAL==1) numeosquantitiesinstandardlist[i]+=NUMEXTRAEOSQUANTITIESTYPE1;
  
  // stored degen quantities:
  numeosdegenquantitiesinstandardlist[i]=NUMEOSDEGENQUANTITIESin;


}





// setup expected read-in sizes for tables
// note that each primary table has a degen and temperature partner separated in memory (although degen is separate file while temperature is in file)
static void settablesizeexpected(int tablesizeexpected[][NUMEOSINDEPS])
{
  int i;
  int templikeiter;

  ///////////
  //
  // number of things corresponds to read-in number of quantities, not final table sizes
  //
  ///////////

  // FULLTABLE
  i=0;
  tablesizeexpected[FULLTABLE][i]=EOSFULLN1; i++;
  for(templikeiter=FIRSTTKLIKE;templikeiter<=LASTTKLIKE;templikeiter++){
    tablesizeexpected[FULLTABLE][i]=EOSFULLN2; i++;
  }
  tablesizeexpected[FULLTABLE][i]=EOSFULLN3; i++;
  tablesizeexpected[FULLTABLE][i]=EOSFULLN4; i++;
  tablesizeexpected[FULLTABLE][i]=EOSFULLN5; i++;
  if(i!=NUMEOSINDEPS){
    dualfprintf(fail_file,"tablesizeexpected(full) not setup for that many indepdimens: i=%d\n",i);
    myexit(3206882);
  }

  // EXTRAFULLTABLE
  i=0;
  tablesizeexpected[EXTRAFULLTABLE][i]=EOSEXTRAFULLN1; i++;
  for(templikeiter=FIRSTTKLIKE;templikeiter<=LASTTKLIKE;templikeiter++){
    tablesizeexpected[EXTRAFULLTABLE][i]=EOSEXTRAFULLN2; i++;
  }
  tablesizeexpected[EXTRAFULLTABLE][i]=EOSEXTRAFULLN3; i++;
  tablesizeexpected[EXTRAFULLTABLE][i]=EOSEXTRAFULLN4; i++;
  tablesizeexpected[EXTRAFULLTABLE][i]=EOSEXTRAFULLN5; i++;
  if(i!=NUMEOSINDEPS){
    dualfprintf(fail_file,"tablesizeexpected(full) not setup for that many indepdimens: i=%d\n",i);
    myexit(3206882);
  }

  // SIMPLETABLE
  i=0;
  tablesizeexpected[SIMPLETABLE][i]=EOSSIMPLEN1; i++;
  for(templikeiter=FIRSTTKLIKE;templikeiter<=LASTTKLIKE;templikeiter++){
    tablesizeexpected[SIMPLETABLE][i]=EOSSIMPLEN2; i++;
  }
  tablesizeexpected[SIMPLETABLE][i]=EOSSIMPLEN3; i++;
  tablesizeexpected[SIMPLETABLE][i]=EOSSIMPLEN4; i++;
  tablesizeexpected[SIMPLETABLE][i]=EOSSIMPLEN5; i++;
  if(i!=NUMEOSINDEPS){
    dualfprintf(fail_file,"tablesizeexpected(simple) not setup for that many indepdimens: i=%d\n",i);
    myexit(3206883);
  }
  // EXTRASIMPLETABLE
  i=0;
  tablesizeexpected[EXTRASIMPLETABLE][i]=EOSEXTRASIMPLEN1; i++;
  for(templikeiter=FIRSTTKLIKE;templikeiter<=LASTTKLIKE;templikeiter++){
    tablesizeexpected[EXTRASIMPLETABLE][i]=EOSEXTRASIMPLEN2; i++;
  }
  tablesizeexpected[EXTRASIMPLETABLE][i]=EOSEXTRASIMPLEN3; i++;
  tablesizeexpected[EXTRASIMPLETABLE][i]=EOSEXTRASIMPLEN4; i++;
  tablesizeexpected[EXTRASIMPLETABLE][i]=EOSEXTRASIMPLEN5; i++;
  if(i!=NUMEOSINDEPS){
    dualfprintf(fail_file,"tablesizeexpected(simple) not setup for that many indepdimens: i=%d\n",i);
    myexit(3206883);
  }

  // SIMPLEZOOMTABLE
  i=0;
  tablesizeexpected[SIMPLEZOOMTABLE][i]=EOSSIMPLEZOOMN1; i++;
  for(templikeiter=FIRSTTKLIKE;templikeiter<=LASTTKLIKE;templikeiter++){
    tablesizeexpected[SIMPLEZOOMTABLE][i]=EOSSIMPLEZOOMN2; i++;
  }
  tablesizeexpected[SIMPLEZOOMTABLE][i]=EOSSIMPLEZOOMN3; i++;
  tablesizeexpected[SIMPLEZOOMTABLE][i]=EOSSIMPLEZOOMN4; i++;
  tablesizeexpected[SIMPLEZOOMTABLE][i]=EOSSIMPLEZOOMN5; i++;
  if(i!=NUMEOSINDEPS){
    dualfprintf(fail_file,"tablesizeexpected(simplezoom) not setup for that many indepdimens: i=%d\n",i);
    myexit(3206884);
  }

  // EXTRASIMPLEZOOMTABLE
  i=0;
  tablesizeexpected[EXTRASIMPLEZOOMTABLE][i]=EOSEXTRASIMPLEZOOMN1; i++;
  for(templikeiter=FIRSTTKLIKE;templikeiter<=LASTTKLIKE;templikeiter++){
    tablesizeexpected[EXTRASIMPLEZOOMTABLE][i]=EOSEXTRASIMPLEZOOMN2; i++;
  }
  tablesizeexpected[EXTRASIMPLEZOOMTABLE][i]=EOSEXTRASIMPLEZOOMN3; i++;
  tablesizeexpected[EXTRASIMPLEZOOMTABLE][i]=EOSEXTRASIMPLEZOOMN4; i++;
  tablesizeexpected[EXTRASIMPLEZOOMTABLE][i]=EOSEXTRASIMPLEZOOMN5; i++;
  if(i!=NUMEOSINDEPS){
    dualfprintf(fail_file,"tablesizeexpected(simplezoom) not setup for that many indepdimens: i=%d\n",i);
    myexit(3206884);
  }




}





static void settablesubtypeinfo(void)
{
  int iteri;


  ///////////////////////////////////////////////
  //
  // Set table subtype sizes [each subtype has same independent variable for simplicity of the lookup table -- rarely need info from multiple subtables anyways]
  //
  ///////////////////////////////////////////////

  numcolintablesubtype[SUBTYPEDEGEN]=NUMEOSDEGENQUANTITIESMEM; // more complicated table, dealt with in special way, but can still access table quantities directly per whichd
  numcolintablesubtype[SUBTYPESTANDARD]=NUMEOSSTANDARDQUANTITIESMEM;
  numcolintablesubtype[SUBTYPEGUESS]=NUMEOSGUESSQUANTITIESMEM;
  numcolintablesubtype[SUBTYPEDISS]=NUMEOSDISSQUANTITIESMEM;
  numcolintablesubtype[SUBTYPEDP]=NUMEOSDPQUANTITIESMEM;
  numcolintablesubtype[SUBTYPESDEN]=NUMEOSSDENQUANTITIESMEM;
  numcolintablesubtype[SUBTYPESSPEC]=NUMEOSSSPECQUANTITIESMEM;
  numcolintablesubtype[SUBTYPEPOFCHI]=NUMEOSPOFCHIQUANTITIESMEM;
  numcolintablesubtype[SUBTYPETEMP]=NUMEOSTEMPQUANTITIESMEM; // table accessed in special way based upon "whichd" during lookup

  // below is one of UTOTDIFF,PTOTDIFF,CHIDIFF,STOTDIFF (there are NUMEOSDEGENQUANTITIESMEM1 of these)
  whichdintablesubtype[SUBTYPEDEGEN]=NOSUCHDIFF; // more complicated table, dealt with in special way
  whichdintablesubtype[SUBTYPESTANDARD]=UTOTDIFF;
  whichdintablesubtype[SUBTYPEGUESS]=PTOTDIFF;
  whichdintablesubtype[SUBTYPEDISS]=STOTDIFF;
  whichdintablesubtype[SUBTYPEDP]=UTOTDIFF;
  whichdintablesubtype[SUBTYPESDEN]=UTOTDIFF;
  whichdintablesubtype[SUBTYPESSPEC]=CHIDIFF;
  whichdintablesubtype[SUBTYPEPOFCHI]=CHIDIFF;
  whichdintablesubtype[SUBTYPETEMP]=NOSUCHDIFF; // table accessed in special way based upon "whichd" during lookup

  // Get "firsteos" for setting whichfun=coli+firsteos where coli starts from 0 for all subtables
  firsteosintablesubtype[SUBTYPEDEGEN]=FIRSTEOSDEGEN;
  firsteosintablesubtype[SUBTYPESTANDARD]=FIRSTEOSSTANDARD;
  firsteosintablesubtype[SUBTYPEGUESS]=FIRSTEOSGUESS;
  firsteosintablesubtype[SUBTYPEDISS]=FIRSTEOSDISS;
  firsteosintablesubtype[SUBTYPEDP]=FIRSTEOSDP;
  firsteosintablesubtype[SUBTYPESDEN]=FIRSTEOSSDEN;
  firsteosintablesubtype[SUBTYPESSPEC]=FIRSTEOSSSPEC;
  firsteosintablesubtype[SUBTYPEPOFCHI]=FIRSTEOSPOFCHI;
  firsteosintablesubtype[SUBTYPETEMP]=FIRSTEOSTEMP;

  
  // mapping from general "whichfun" list to "whichtablesubtype" list
  for(iteri=FIRSTEOSDEGEN;iteri<=LASTEOSDEGEN;iteri++) whichtablesubtypeinquantity[iteri] = SUBTYPEDEGEN;
  for(iteri=FIRSTEOSSTANDARD;iteri<=LASTEOSSTANDARD;iteri++) whichtablesubtypeinquantity[iteri] = SUBTYPESTANDARD;
  for(iteri=FIRSTEOSGUESS;iteri<=LASTEOSGUESS;iteri++) whichtablesubtypeinquantity[iteri] = SUBTYPEGUESS;
  for(iteri=FIRSTEOSDISS;iteri<=LASTEOSDISS;iteri++) whichtablesubtypeinquantity[iteri] = SUBTYPEDISS;
  for(iteri=FIRSTEOSDP;iteri<=LASTEOSDP;iteri++) whichtablesubtypeinquantity[iteri] = SUBTYPEDP;
  for(iteri=FIRSTEOSSDEN;iteri<=LASTEOSSDEN;iteri++) whichtablesubtypeinquantity[iteri] = SUBTYPESDEN;
  for(iteri=FIRSTEOSSSPEC;iteri<=LASTEOSSSPEC;iteri++) whichtablesubtypeinquantity[iteri] = SUBTYPESSPEC;
  for(iteri=FIRSTEOSPOFCHI;iteri<=LASTEOSPOFCHI;iteri++) whichtablesubtypeinquantity[iteri] = SUBTYPEPOFCHI;
  for(iteri=FIRSTEOSTEMP;iteri<=LASTEOSTEMP;iteri++) whichtablesubtypeinquantity[iteri] = SUBTYPETEMP;
  for(iteri=FIRSTEOSEXTRA;iteri<=LASTEOSEXTRA;iteri++) whichtablesubtypeinquantity[iteri] = SUBTYPEEXTRA;

  // mapping from general "whichfun" list to "whichcol" for a given "whichtablesubtype"
  for(iteri=FIRSTEOSDEGEN;iteri<=LASTEOSDEGEN;iteri++) whichcolinquantity[iteri] = iteri-FIRSTEOSDEGEN;
  for(iteri=FIRSTEOSSTANDARD;iteri<=LASTEOSSTANDARD;iteri++) whichcolinquantity[iteri] = iteri-FIRSTEOSSTANDARD;
  for(iteri=FIRSTEOSGUESS;iteri<=LASTEOSGUESS;iteri++) whichcolinquantity[iteri] = iteri-FIRSTEOSGUESS;
  for(iteri=FIRSTEOSDISS;iteri<=LASTEOSDISS;iteri++) whichcolinquantity[iteri] = iteri-FIRSTEOSDISS;
  for(iteri=FIRSTEOSDP;iteri<=LASTEOSDP;iteri++) whichcolinquantity[iteri] = iteri-FIRSTEOSDP;
  for(iteri=FIRSTEOSSDEN;iteri<=LASTEOSSDEN;iteri++) whichcolinquantity[iteri] = iteri-FIRSTEOSSDEN;
  for(iteri=FIRSTEOSSSPEC;iteri<=LASTEOSSSPEC;iteri++) whichcolinquantity[iteri] = iteri-FIRSTEOSSSPEC;
  for(iteri=FIRSTEOSPOFCHI;iteri<=LASTEOSPOFCHI;iteri++) whichcolinquantity[iteri] = iteri-FIRSTEOSPOFCHI;
  for(iteri=FIRSTEOSTEMP;iteri<=LASTEOSTEMP;iteri++) whichcolinquantity[iteri] = iteri-FIRSTEOSTEMP;
  for(iteri=FIRSTEOSEXTRA;iteri<=LASTEOSEXTRA;iteri++) whichcolinquantity[iteri] = iteri-FIRSTEOSEXTRA;

  //  for(iteri=FIRSTEOSQUANTITY;iteri<=LASTEOSQUANTITY;iteri++){ // NUMEOSQUANTITIES of them
  //    get_dologinterp_subtype(s, int tabledimen, int degentable, int whichtable, int whichtablesubtype, int coli, int whichindep)
  //  }

  int itersubtype,coli;
  for(itersubtype=0;itersubtype<NUMTABLESUBTYPES;itersubtype++){
    for(coli=0;coli<numcolintablesubtype[itersubtype];coli++){
      dologinterp_sub_coli[itersubtype][coli]=get_dologinterp_subtype(itersubtype, coli);
    }
  }




}








// note we define over full LASTINDEPDIMEN not just up to LASTINDEPDIMENUSED, so anything that tries to access mapping can still do it, even for non-existent dimension
static void setvartypeinfo(void)
{
  int i;

  ///////////
  //
  // setup what q1-q5 mean associated with the 5 dimensions of the arrays related to NUMINDEPDIMENS and not NUMEOSINDEPS
  //
  ///////////
  i=1;
  vartypearray[i]=RHOEOS;  i++;
  vartypearray[i]=UEOS;    i++; // UEOS used for reading table, but later used to change vartypearraylocal for each whichindep or whichd
  vartypearray[i]=YEEOS;   i++;
  vartypearray[i]=YNUEOS;  i++;
  vartypearray[i]=HEOS;    i++;
  if(i!=NUMINDEPDIMENS+1){
    dualfprintf(fail_file,"vartypearray not setup for that many indepdimens: i=%d\n",i);
    myexit(3206885);
  }


  ////////////
  //
  // for accessing EOSextra, setup what q1-q5 mean associated with the 5 dimensions of the arrays related to NUMINDEPDIMENS and not NUMEOSINDEPS
  // note that vartypeeosextraarray[0,1] should never be used since would map to EOSextraglobal[-2,-1] that doesn't exist
  //
  ////////////
  i=1;
  vartypeeosextraarray[i]=RHOGLOBAL;       i++;
  vartypeeosextraarray[i]=UGLOBAL;         i++;
  vartypeeosextraarray[i]=TDYNORYEGLOBAL;  i++;
  vartypeeosextraarray[i]=YNU0GLOBAL;      i++;
  vartypeeosextraarray[i]=HGLOBAL;         i++;
  if(i!=NUMINDEPDIMENS+1){
    dualfprintf(fail_file,"vartypeeosextraarray not setup for that many indepdimens: i=%d\n",i);
    myexit(9213825);
  }


  ///////////
  //
  // for accessing EOSextra's height's
  //
  ///////////
  i=1;
  vartypeheightarray[i]=HGLOBAL;   i++;
  vartypeheightarray[i]=H2GLOBAL;  i++;
  vartypeheightarray[i]=H3GLOBAL;  i++;
  vartypeheightarray[i]=H4GLOBAL;  i++;
  if(i!=NUMHDIRECTIONS+1){
    dualfprintf(fail_file,"vartypeheightarray not setup for that many heights: i=%d\n",i);
    myexit(2982352);
  }


  ///////////
  //
  // setup what degen type q1-q4 mean associated with the 5 dimensions of the arrays related to NUMEOSDEGENINDEPS (similar to NUMINDEPDIMENS since all difference between NUMINDEPDIMENS and not NUMEOSINDEPS are temperature-like quantities)
  //
  ///////////
  i=1;
  vardegentypearray[i]=RHOEOSDEGEN;  varnormalcompare2degentypearray[i] = RHOEOS;  i++;
  vardegentypearray[i]=YEEOSDEGEN;   varnormalcompare2degentypearray[i] = YEEOS;   i++;
  vardegentypearray[i]=YNUEOSDEGEN;  varnormalcompare2degentypearray[i] = YNUEOS;  i++;
  vardegentypearray[i]=HEOSDEGEN;    varnormalcompare2degentypearray[i] = HEOS;    i++;
  if(i!=NUMEOSDEGENINDEPS+1){
    dualfprintf(fail_file,"vardegentypearray and varnormalcompare2degentypearray not setup for that many indepdimens: i=%d\n",i);
    myexit(3206886);
  }



}











////////////////////////
//
// reads in table for EOS and sets up the table parameters
//
////////////////////////
void read_setup_eostable(void)
{
  FILE *inhead;
  FILE *intable;
  FILE *indegentable;
  int ii,jj;
  int m,n,o,p,q; // 5 dimension labels
  int iii,jjj,kkk,lll,mmm; // for 5 dimensions of table
  int ppp,qqq;
  int totalindex[NUMEOSINDEPS];
  FTYPE indep[NUMEOSINDEPS];
  FTYPE indepplusdegen[NUMEOSDEGENINDEPS];
  FTYPE indepdegen[NUMEOSDEGENINDEPS];
  FTYPE lstepdep,diff,lindep,lindeptry;
  char headername[NUMTBLS][MAXFILENAME];
  char tablename[NUMTBLS][MAXFILENAME];
  char degentablename[NUMTBLS][MAXFILENAME];
  int tableiter;
  int tablesizeexpected[NUMTBLS][NUMEOSINDEPS];
  double tabletemp[NUMEOSQUANTITIESMEM]; // size of read-in and stored columns
  double degentabletemp[NUMEOSDEGENQUANTITIESin]; // size of read-in and stored columns
  int sizematch; // 0 = table sizes don't match  1 = they match
  double valuetemp[MAXNUMEOSQUANTITIESin],degenvaluetemp[NUMEOSDEGENQUANTITIESin];
  FTYPE errordegen,errordegen2;
  int i;
  int testnumfscanfquantities,testnumfscanfquantitiesdegen;
  int templikeiter;
  long long int numberofdegenerrors[4],totalnumberofdegen; // 0=TABLETOL 1=0.01 2=0.1 3=TABLETOLTRUNCATION




  trifprintf("Setting up Kaz EOS table\n");



  if(ALLOWKAZEOS==0){
    dualfprintf(fail_file,"Need to have ALLOWKAZEOS==1 if using Kaz EOS\n");
    myexit(93458635);
  }



  ///////////////////////////////////////////////
  //
  // Set which table type is primary
  //
  ///////////////////////////////////////////////
  if(ALLOWFULLTABLE) primarytable=FULLTABLE;
  else if(ALLOWSIMPLETABLE) primarytable=SIMPLETABLE;
  else if(ALLOWSIMPLEZOOMTABLE) primarytable=SIMPLEZOOMTABLE;



  ///////////////////////////////////////////////
  //
  // Set table subtype information
  //
  ///////////////////////////////////////////////
  settablesubtypeinfo();


  ///////////////////////////////////////////////
  //
  // Set expected number of columns of dataf or each table
  //
  ///////////////////////////////////////////////
  for(i=0;i<NUMTBLS;i++) assigntablecolumnnumbers(i);

  
  ///////////////////////////////////////////////
  //
  // Set expected sizes of tables
  //
  ///////////////////////////////////////////////
  settablesizeexpected(tablesizeexpected);
  

  ///////////////////////////////////////////////
  //
  // set vartype translation arrays
  //
  ///////////////////////////////////////////////
  setvartypeinfo();


  ///////////////
  //
  // set pointer shifts for EOS tables and global spatial array
  //
  ///////////////
  kazfulleos_set_arrays();


  







  ///////////////////////////////////////////////
  //
  // Get tables from files
  //
  //////////////////////////////////////////////

  // do only over first CPU (myid==0) and then send data to rest of CPUs
  if(myid==0){



    ///////////////////////////////////////////////
    //
    // Set file names
    //
    ///////////////////////////////////////////////

    strcpy(headername[FULLTABLE],EOSFULLHEADNAME);
    strcpy(tablename[FULLTABLE],EOSFULLTABLENAME);
    strcpy(degentablename[FULLTABLE],EOSFULLTABLEDEGENNAME);

    strcpy(headername[EXTRAFULLTABLE],EOSEXTRAFULLHEADNAME);
    strcpy(tablename[EXTRAFULLTABLE],EOSEXTRAFULLTABLENAME);
    strcpy(degentablename[EXTRAFULLTABLE],EOSEXTRAFULLTABLEDEGENNAME);

    strcpy(headername[SIMPLETABLE],EOSSIMPLEHEADNAME);
    strcpy(tablename[SIMPLETABLE],EOSSIMPLETABLENAME);
    strcpy(degentablename[SIMPLETABLE],EOSSIMPLETABLEDEGENNAME);

    strcpy(headername[EXTRASIMPLETABLE],EOSEXTRASIMPLEHEADNAME);
    strcpy(tablename[EXTRASIMPLETABLE],EOSEXTRASIMPLETABLENAME);
    strcpy(degentablename[EXTRASIMPLETABLE],EOSEXTRASIMPLETABLEDEGENNAME);

    strcpy(headername[SIMPLEZOOMTABLE],EOSSIMPLEZOOMHEADNAME);
    strcpy(tablename[SIMPLEZOOMTABLE],EOSSIMPLEZOOMTABLENAME);
    strcpy(degentablename[SIMPLEZOOMTABLE],EOSSIMPLEZOOMTABLEDEGENNAME);

    strcpy(headername[EXTRASIMPLEZOOMTABLE],EOSEXTRASIMPLEZOOMHEADNAME);
    strcpy(tablename[EXTRASIMPLEZOOMTABLE],EOSEXTRASIMPLEZOOMTABLENAME);
    strcpy(degentablename[EXTRASIMPLEZOOMTABLE],EOSEXTRASIMPLEZOOMTABLEDEGENNAME);
 




    ///////////////
    //
    // loop over normal tables (degen table read-in with its normal version)
    //
    //////////////

    for(tableiter=0;tableiter<NUMTBLS;tableiter++){


      // avoid tables not needed since should not expect user has copied or created such tables if not using them
      if(tableiter==FULLTABLE && ALLOWFULLTABLE==0) continue;
      if(tableiter==EXTRAFULLTABLE && (ALLOWFULLTABLE==0 || WHICHDATATYPEGENERAL!=4) ) continue;
      if(tableiter==SIMPLETABLE && ALLOWSIMPLETABLE==0) continue;
      if(tableiter==EXTRASIMPLETABLE && (ALLOWSIMPLETABLE==0 || WHICHDATATYPEGENERAL!=4) ) continue;
      if(tableiter==SIMPLEZOOMTABLE && ALLOWSIMPLEZOOMTABLE==0) continue;
      if(tableiter==EXTRASIMPLEZOOMTABLE && (ALLOWSIMPLEZOOMTABLE==0 || WHICHDATATYPEGENERAL!=4) ) continue;



      ///////////////////////////////////
      //
      // first read-in the header file (only normal table has header file since degen table is same except one dimension (u/p/chi))
      //
      //////////////////////////////////
      
      // open header file
      if( (inhead = fopen(headername[tableiter],"rb"))==NULL){
	dualfprintf(fail_file,"No such file %s\n",headername[tableiter]);
	myexit(16622);
      }
      
      // get number of output colums and check that code and data file agree
      // also get whichmethod using
      fscanf(inhead,"%d %d %d",&whichrnpmethod[tableiter],&whichynumethod[tableiter],&whichhcmmethod[tableiter]);
      fscanf(inhead,"%d %d %d %d",&whichdatatype[tableiter],&utotdegencut[tableiter],&numc[tableiter],&numextras[tableiter]);




      ////////////////////////////
      //
      // perform some checks
      //
      ////////////////////////////

      if(whichrnpmethod[tableiter]==0){
	fprintf(fail_file,"This method is not setup\n");
	myexit(3966738);
      }

      if(numextras[tableiter]>MAXNUMEXTRAS){
	dualfprintf(fail_file,"Increase MAXNUMEXTRAS to %d\n",numextras[tableiter]);
	myexit(626236);
      }

      // check that number of extras is correct
      if(isextratable(tableiter) && WHICHDATATYPEGENERAL==4 || WHICHDATATYPEGENERAL!=4){
	// only check if extras should be there
	
	// now check how many extras should be there
	if(
	   (WHICHDATATYPEGENERAL==1 && numextras[tableiter]!=NUMEXTRAEOSQUANTITIESTYPE1) ||
	   (WHICHDATATYPEGENERAL==2 && numextras[tableiter]!=NUMEXTRAEOSQUANTITIESTYPE2) ||
	   (WHICHDATATYPEGENERAL==3 && numextras[tableiter]!=NUMEXTRAEOSQUANTITIESTYPE3) ||
	   (WHICHDATATYPEGENERAL==4 && numextras[tableiter]!=NUMEXTRAEOSQUANTITIESTYPE4)
	   ){
	  dualfprintf(fail_file,"Wrong number of extras.  header says: %d\n",numextras[tableiter]);
	  myexit(8923742);
	}
	// else correct
      }
      else{
	// then should be 0 extras
	if(numextras[tableiter]!=0){
	  dualfprintf(fail_file,"Wrong number of extras.  header says: %d but should be 0\n",numextras[tableiter]);
	  myexit(8923743);
	}
      }



      if(whichdatatype[tableiter]>MAXNUMDATATYPES){ // no longer relevant code
	dualfprintf(fail_file,"Increase MAXNUMDATATYPES to %d\n",whichdatatype[tableiter]);
	myexit(626237);
      }

      // new check since optimized everything for only 1 datatype per run
      if(whichdatatype[tableiter]!=WHICHDATATYPEGENERAL){
	dualfprintf(fail_file,"whichdatatype=%d while WHICHDATATYPEGENERAL=%d\n",whichdatatype[tableiter],WHICHDATATYPEGENERAL);
	myexit(239752266);
      }



      //////////////
      //
      // check number of quantities in normal (non-degen table)
      //
      //////////////
      if(numc[tableiter]!=numallquantitiesinfile[tableiter]){
	dualfprintf(fail_file,"numcolumns=%d but code expected %d\n",numc[tableiter],numallquantitiesinfile[tableiter]);
	dualfprintf(fail_file,"tableiter=%d whichdatatype=%d :: %d %d %d %d\n",tableiter,whichdatatype[tableiter],NUMINDEPDIMENS,NUMEOSINDEPS,NUMEOSDEGENQUANTITIESMEM1,numeosquantitiesinfile[tableiter]);
	myexit(16623);
      }
      


      //////////////
      //
      // Read-in table header
      //
      //////////////
      for(ii=0;ii<NUMEOSINDEPS;ii++){
	fscanf(inhead,"%d",&tablesize[tableiter][ii]);
	// second [4] : 0 = lower log_base limit, 1 = upper log_base limit, 2=step, 3 = divisor of grid position 4=base of log, 5 = linear value of offset for log_base stepping so can control how resolved
 	fscanf(inhead,HEADERONEIN,&inputtablelimits[tableiter][ii][0]); // start in log
	fscanf(inhead,HEADERONEIN,&inputtablelimits[tableiter][ii][1]); // end in log
	fscanf(inhead,HEADERONEIN,&inputtablelimits[tableiter][ii][2]); // step in log
	fscanf(inhead,HEADERONEIN,&inputtablelimits[tableiter][ii][4]); // base of log offset
	fscanf(inhead,HEADERONEIN,&inputtablelimits[tableiter][ii][5]); // linear offset

	if(inputtablelimits[tableiter][ii][5]>=0.0 && ii==YNUEOS){
	  dualfprintf(fail_file,"SUPERWARNING: Note that table setup so Ynu<=0 not allowed since log-indexing Ynu as independent variable.\n");
	}


	// convert base 10 to actual base with actual offset so can be easily used for mapping
	// x_in:=
	tablelimits[tableiter][ii][0] = log10(pow(10.0,inputtablelimits[tableiter][ii][0])-inputtablelimits[tableiter][ii][5])/log10(inputtablelimits[tableiter][ii][4]);
	// x_out:=
	tablelimits[tableiter][ii][1] = log10(pow(10.0,inputtablelimits[tableiter][ii][1])-inputtablelimits[tableiter][ii][5])/log10(inputtablelimits[tableiter][ii][4]);
	
	tablelimits[tableiter][ii][4]=inputtablelimits[tableiter][ii][4]; // just copy base, assumed to always be log_base() where the base itself is in base 10 format
	tablelimits[tableiter][ii][5]=inputtablelimits[tableiter][ii][5]; // shouldn't need to use linear version of offset, but just copy for now
	


	if(fabs(inputtablelimits[tableiter][ii][2])<SMALL || tablesize[tableiter][ii]==1 ){
	  // then assume that direction has no dimensionality
	  inputtablelimits[tableiter][ii][2] = 1.0;

	  tablelimits[tableiter][ii][2] = 0.0; // forced
      
	  tablelimits[tableiter][ii][3] = 0.0; // forced so index=0 always

	}
	else{

	  // compute log_base step to be used in formula:
	  // r = r_0 + base^(x_in+dx*i) for i=0..N-1 where x = x_in+dx*i
	  // dx = (x_out-x_in)/(N-1)
	  // dx=:
	  tablelimits[tableiter][ii][2] = (tablelimits[tableiter][ii][1]-tablelimits[tableiter][ii][0])/((FTYPE)tablesize[tableiter][ii]-1.0);
      
	  // the below definition is consistent with Kaz's code, matlab eos_extract.m and elsewhere in this code
	  //	tablelimits[tableiter][ii][3] = ((FTYPE)tablesize[tableiter][ii]-1.0)/(tablelimits[tableiter][ii][1] - tablelimits[tableiter][ii][0]);


	  // i = (log_b (r-r_0) - x_in)/dx
	  // 1/dx:=
	  tablelimits[tableiter][ii][3] = 1.0/tablelimits[tableiter][ii][2];

	}
	
	// below used to truncate(limit) input values so lookup doesn't have as many conditionals
	// lineartablelimits does NOT include degen offset, so these are combined with offset when used since offset is function of rhob,hcm,tdynorye
	if((utotdegencut[tableiter]>DEGENCUTLASTOLDVERSION)&&(ii>=FIRSTTKLIKE && ii<=LASTTKLIKE)){
	  // no log-conversion for limits
	  lineartablelimits[tableiter][ii][0]=inputtablelimits[tableiter][ii][0];
	  lineartablelimits[tableiter][ii][1]=inputtablelimits[tableiter][ii][1];
	  lineartablelimits[tableiter][ii][2]=1.0/log10(inputtablelimits[tableiter][ii][4]); // 1.0 divided by "normal log(base)" (used for mapping)
	  lineartablelimits[tableiter][ii][3]=inputtablelimits[tableiter][ii][5]; // linear offset
	}
	else{
	  lineartablelimits[tableiter][ii][0]=pow(10.0,inputtablelimits[tableiter][ii][0]); // linear inner
	  lineartablelimits[tableiter][ii][1]=pow(10.0,inputtablelimits[tableiter][ii][1]); // linear outer
	  lineartablelimits[tableiter][ii][2]=1.0/log10(inputtablelimits[tableiter][ii][4]); // 1.0 divided by "normal log(base)" (used for mapping)
	  lineartablelimits[tableiter][ii][3]=inputtablelimits[tableiter][ii][5]; // linear offset
	}
	
      }// done over NUMEOSINDEPS





      ////////////////////////
      //
      // Also read in lsoffset and fakelsoffset
      fscanf(inhead,HEADERONEIN,&lsoffset);
      fscanf(inhead,HEADERONEIN,&fakelsoffset);

      // GODMARK: should really read-in from table, although expected tabular value is 9.14Mev and this gives 7.108 for smooth connection for initial conditions
      // GODMARK: The below is for Shen EOS only
      FAKE2IDEALNUCLEAROFFSET = (-7.57E-3);

      // Exactly 9.14MeV/baryon as in helm/jon_lsbox.f
      // This offset was subtracted before tabulating in log-log.  After lookup, this energy/baryon needs to be added back in the correc tunits.
      // ergPmev=1.60218E-6
      //  TRUENUCLEAROFFSET=(9.14*ergPmev/(mb*C*C));// cgs/cgs = dimensionless quantity
      // above gives 9.7346E-3
      //  TRUENUCLEAROFFSET /= (1.0); // need to get code version of energy/baryon however used to add to internal energy
      // assume degen offset accounts for offset, that cannot be put into EOS itself!
      // Now read-in unphysical fakelsoffset used to avoid negative energies is required to reoffset internal energy back so only physical "lsoffset" is included
      TRUENUCLEAROFFSET=fakelsoffset*ergPmev/(mb*C*C);

      // degeneracy global offset:
      // ergPmev=1.60218E-6
      //  DEGENNUCLEAROFFSET=(50.0*ergPmev/(mb*C*C));// cgs/cgs = dimensionless quantity
      DEGENNUCLEAROFFSET=0.0; // avoid this approach for now -- try putting in offset into original HELM table first

      //
      ////////////////////////


      /////////////////
      //
      // Temperature table limits in header not read in, just for user use
      //
      /////////////////
      





      //////////////////////
      //
      // check size of tables
      //
      //////////////////////
      sizematch=1; // assume match
      for(ii=0;ii<NUMEOSINDEPS;ii++){
	sizematch *= tablesizeexpected[tableiter][ii]==tablesize[tableiter][ii];
      }

      if(!sizematch){
	dualfprintf(fail_file,"Size of table does not match (tableiter=%d)\n",tableiter);
	for(ii=0;ii<NUMEOSINDEPS;ii++){
	  dualfprintf(fail_file,"read-in tablesize[%d]=%d\n",ii,tablesize[tableiter][ii]);
	}
	dualfprintf(fail_file,"EOSN's = ");
	for(ii=0;ii<NUMEOSINDEPS;ii++){
	  dualfprintf(fail_file," %d",tablesizeexpected[tableiter][ii]);
	}
	dualfprintf(fail_file,"\n");
	myexit(16624);
      }




      ////////////
      //
      // done with header, so close file
      //
      ////////////
      fclose(inhead);
      


    



      ////////////////////////////
      //
      // at this point header is consistent with expectations, so read in the table
      //
      // now read-in table as written
      //
      // Each table has its own degen table, so process both at the same time
      //
      //////////////////////////////////
      
      // open eostable file
      if( (intable = fopen(tablename[tableiter],"rb"))==NULL){
	dualfprintf(fail_file,"No such file %s\n",tablename[tableiter]);
	myexit(16625);
      }


      // open eosdegentable file
      if( (indegentable = fopen(degentablename[tableiter],"rb"))==NULL){
	dualfprintf(fail_file,"No such file %s\n",degentablename[tableiter]);
	myexit(2382662);
      }


      // initialize table degen errors:
      numberofdegenerrors[0]=numberofdegenerrors[1]=numberofdegenerrors[2]=numberofdegenerrors[3]=0;
      totalnumberofdegen=0;


      // file has rhob as fastest index
      // here assumes tablesize of UEOS, PEOS, CHIEOS, and SEOS are same
      // assume jjj=0 to start since degen checks below depend on that
      // notice that vartypearray has same size as dimension of arrays and loops.
      for(mmm=0;mmm<tablesize[tableiter][vartypearray[5]];mmm++)for(lll=0;lll<tablesize[tableiter][vartypearray[4]];lll++)for(kkk=0;kkk<tablesize[tableiter][vartypearray[3]];kkk++)for(jjj=0;jjj<tablesize[tableiter][vartypearray[2]];jjj++)for(iii=0;iii<tablesize[tableiter][vartypearray[1]];iii++){
	testnumfscanfquantities=0;
	testnumfscanfquantitiesdegen=0;

	// first, get positions to make sure everything is consistent
	fscanf(intable,"%d %d %d %d %d",&m,&n,&o,&p,&q); // indexing related to NUMINDEPDIMENS not NUMEOSINDEPS
	testnumfscanfquantities += NUMINDEPDIMENS;

	if(m!=iii || n!=jjj || o!=kkk || p!=lll || q!=mmm){
	  dualfprintf(fail_file,"Read-in table (%d) indicies inconsistent with expected indicies: m=%d iii=%d n=%d jjj=%d o=%d kkk=%d p=%d lll=%d q=%d mmm=%d\n",tableiter,m,iii,n,jjj,o,kkk,p,lll,q,mmm);
	  dualfprintf(fail_file,"whichrnpmethod=%d whichynumethod=%d whichhcmmethod=%d\n",whichrnpmethod[tableiter],whichynumethod[tableiter],whichhcmmethod[tableiter]);
	  dualfprintf(fail_file,"whichdatatype=%d\n",whichdatatype[tableiter]);
	  for(jj=0;jj<NUMEOSINDEPS;jj++){
	    dualfprintf(fail_file,"tablesize[%d][%d]=%d\n",tableiter,jj,tablesize[tableiter][jj]);
	  }
	  dualfprintf(fail_file,"Total number of EOS quantities=%d (Ensure file has correct number of columns!)\n",numallquantitiesinfile[tableiter]);
	  
	  myexit(16626);
	}

	if(jjj==0){ // degen table
	  n=jjj;
	  fscanf(indegentable,"%d %d %d %d",&m,&o,&p,&q); // related to NUMEOSDEGENINDEPS
	  testnumfscanfquantitiesdegen += NUMEOSDEGENINDEPS;
	  if(m!=iii || n!=jjj || o!=kkk || p!=lll || q!=mmm){
	    dualfprintf(fail_file,"Read-in degentable indicies inconsistent with expected indicies: m=%d iii=%d n=%d jjj=%d o=%d kkk=%d p=%d lll=%d q=%d mmm=%d\n",tableiter,m,iii,n,jjj,o,kkk,p,lll,q,mmm);
	    myexit(166265);
	  }
	}




	/////////////////////
	//
	// second, read in the independent variable values and compare with expected value
	//
	/////////////////////
	for(ii=0;ii<NUMEOSINDEPS;ii++){
	  fscanf(intable,"%lf",&indep[ii]); // rhob, Udiff, Pdiff, CHIdiff, Sdiff, tdynorye, tdynorynu, H   for given grid value
	  testnumfscanfquantities += 1;
	}
	// read-in UofUdiff, PofPdiff, CHIofCHIdiff, SofSdiff -- used to check degen offset calculation in HARM
	for(ii=0;ii<NUMEOSDEGENQUANTITIESMEM1;ii++){
	  fscanf(intable,"%lf",&indepplusdegen[ii]);
	  testnumfscanfquantities += 1;
	}

	// true independent dimensions associated with free index
	totalindex[vartypearray[1]]   = iii;
	for(templikeiter=FIRSTTKLIKE;templikeiter<=LASTTKLIKE;templikeiter++){
	  totalindex[templikeiter] = jjj; // notice that UEOS,PEOS,CHIEOS,SEOS are actually associated with same index
	}
	totalindex[vartypearray[3]]   = kkk;
	totalindex[vartypearray[4]]   = lll;
	totalindex[vartypearray[5]]   = mmm;
	


	if(jjj==0){// inside degen table
	  for(ii=0;ii<NUMEOSDEGENINDEPS;ii++){
	    fscanf(indegentable,"%lf",&indepdegen[ii]); // rho, tdynorye, tdynorynu, H
	    testnumfscanfquantitiesdegen += 1;
	    // check consistency between normal and degen tablef or independent variables (assumes jjj!=0 in normal table is same for these quantities as jjj==0)
	  
	    if(fabs(indepdegen[vardegentypearray[ii]]-indep[varnormalcompare2degentypearray[ii]])>TABLETOL){
	      dualfprintf(fail_file,"degen table not consistent with normal table (tableiter=%d ii=%d iii=%d jjj=%d kkk=%d lll=%d mmm=%d) for %d %d: %21.15g %21.15g\n",tableiter,ii,iii,jjj,kkk,lll,mmm,vardegentypearray[ii],vartypearray[ii],indepdegen[vardegentypearray[ii]],indep[varnormalcompare2degentypearray[ii]]);
	    }
	  }
	}


	////////////
	//
	// check that read-in quantities agrees with expected number so far:
	//
	////////////
 	if(testnumfscanfquantities!=NUMEOSQUANTITIESNOTSTOREDin){
	  dualfprintf(fail_file,"NOTSTORED number of scanned EOS quantities=%d doesn't agree with expected amount=%d\n",testnumfscanfquantities,NUMEOSQUANTITIESNOTSTOREDin);
	  myexit(1873626);
	}

	if(jjj==0){
	  if(testnumfscanfquantitiesdegen!=NUMEOSDEGENQUANTITIESNOTSTOREDin){
	    dualfprintf(fail_file,"NOTSTORED number of scanned degen EOS quantities=%d doesn't agree with expected amount=%d\n",testnumfscanfquantitiesdegen,NUMEOSDEGENQUANTITIESNOTSTOREDin);
	    myexit(1873627);
	  }
	}



	/////////////////
	//
	// third, check gridding of independent variables in normal table
	//
	/////////////////
	for(ii=0;ii<NUMEOSINDEPS;ii++){


	  if(!isextratable(tableiter) && WHICHDATATYPEGENERAL==4){
	    // then do not check consistency of table for Ynu or H independents since redudant information was removed so no such dimensions even though range is the same as before.
	    if(ii==YNUEOS || ii==HEOS) continue; // skip till iterate ii again
	    // actually, this is already avoided by if(tablesize>1) check below
	  }



	  if(tablesize[tableiter][ii]>1){ // only check actual table and assume degen table consistent

	    // get step (consistent with how step is computed in Kaz's code and in matlab script eos_extract.m)
	    // really only has to be consistent with eos_extract.m
	    //	    lstepdep = (-tablelimits[tableiter][ii][0])/((FTYPE)tablesize[tableiter][ii]-1.0);
	    lstepdep = inputtablelimits[tableiter][ii][2];
	    // compare step sizes to read-in step sizes
	    diff = fabs(lstepdep - tablelimits[tableiter][ii][2])/(fabs(lstepdep)+fabs(tablelimits[tableiter][ii][2]));
	    if(diff > TABLETOL){
	      dualfprintf(fail_file,"Grid step size is incorrect: mmm=%d lll=%d kkk=%d jjj=%d iii=%d :: ii=%d readin-value=%21.15g lstepdep=%21.15g\n",mmm,lll,kkk,jjj,iii,ii,tablelimits[tableiter][ii][2],lstepdep);
	      dualfprintf(fail_file,"tablelimits[%d][%d][0]=%21.15g tablelimits[%d][%d][1]=%21.15g\n",tableiter,ii,tablelimits[tableiter][ii][0],tableiter,ii,tablelimits[tableiter][ii][1]);
	      myexit(16627);
	    }


	    
	    if((utotdegencut[tableiter]>DEGENCUTLASTOLDVERSION)&&(ii>=FIRSTTKLIKE && ii<=LASTTKLIKE)){
	      // grid is just 0..1
	      // then indep already logified and contains all offsets (i.e. lineartablelimits[tableiter][ii][3]=0 required and lineartablelimits[tableiter][ii][2]=1/log10(10.0)=1 required for now)
	      lindep = indep[ii];
	    }
	    else{
	      // grid used is:
	      //
	      // x = log(r-r_0)/log(base) such that r = r_0 + base^x
	      //
	      //
	      // get read-in value of independent variable
	      // get x (here x is now such things as rhob, utotoffset, ptotoffset, chioffset, hcm, tdynorye)
	      lindep=log10(indep[ii]-lineartablelimits[tableiter][ii][3])*lineartablelimits[tableiter][ii][2];
	    }

  
	    // get computed value of independent variable (used later for lookup, so verifies lookup method)
	    // get x using lookup from tablular index (didn't need to change)
	    lindeptry=tablelimits[tableiter][ii][0] + totalindex[ii]*lstepdep;
	    // compare to be sure same
	    //	    diff = fabs(lindep-lindeptry)/(fabs(lindep)+fabs(lindeptry));
	    // normalize below by range of limits instead of local values since otherwise if lindeptry or lindep is near 0 then relative error erroneously is large with old diff above
	    diff = fabs(lindep-lindeptry)/(tablelimits[tableiter][ii][1]-tablelimits[tableiter][ii][0]);
	    if(diff>TABLETOL){
	      dualfprintf(fail_file,"Grid position data is incorrect: mmm=%d lll=%d kkk=%d jjj=%d iii=%d :: ii=%d readin-lindep=%21.15g lindeptry=%21.15g diff=%21.15g\n",mmm,lll,kkk,jjj,iii,ii,lindep,lindeptry,diff);
	      dualfprintf(fail_file,"tablelimits[%d][%d][0]=%21.15g totalindex[%d]=%d lstepdep=%21.15g\n",tableiter,ii,tablelimits[tableiter][ii][0],ii,totalindex[ii],lstepdep);
	      dualfprintf(fail_file,"indep=%21.15g lindep=%21.15g [3]=%21.15g [2]=%21.15g\n",indep[ii],lindep,lineartablelimits[tableiter][ii][3],lineartablelimits[tableiter][ii][2]);
	      dualfprintf(fail_file,"Trigger: %d\n",(utotdegencut[tableiter]>DEGENCUTLASTOLDVERSION)&&(ii>=FIRSTTKLIKE && ii<=LASTTKLIKE));
	      myexit(16628);
	    }
	  }
	}




	////////////////////////
	//
	// fourth, since table gridding is consistent, now read in columns of actual dependent variable values
	//
	////////////////////////


	// normal table
	for(ppp=0;ppp<numeosquantitiesinfile[tableiter];ppp++){ // look over only to-be stored quantities
	  fscanf(intable,DOUBLEINPUT,&valuetemp[ppp]); // double values
	  testnumfscanfquantities += 1;

	  // before storing, reorder to be like standard "in" format for EOS quantities:
	  reorder_extratype(isextratable(i),numeosquantitiesinfile[tableiter],valuetemp);
	  // not in canonical list order/position
	}

	for(ppp=0;ppp<numeosquantitiesinstandardlist[tableiter];ppp++){ // look over only to-be stored quantities
	  set_arrays_eostable(ISNOTDEGENTABLE,tableiter,mmm,lll,kkk,jjj,iii,ppp,valuetemp[ppp]);
	}



	// degen table
	if(jjj==0){
	  for(ppp=0;ppp<numeosdegenquantitiesinfile[tableiter];ppp++){
	    fscanf(indegentable,DOUBLEINPUT,&degenvaluetemp[ppp]); // double values

	    if(degenvaluetemp[ppp]<=0.0 && DOLOGINTERP){
	      dualfprintf(fail_file,"Degenerate table contains non-positive offsets: degenvaluetemp[%d]=%21.15g (tableiter=%d mmm=%d lll=%d kkk=%d iii=%d)\n",ppp,degenvaluetemp[ppp],tableiter,mmm,lll,kkk,iii);
	      dualfprintf(fail_file,"Should use linear interpolation for degenerate table offsets and remove this check.  But note that degenerate offsets are well-described by linear in log-log\n");
	      dualfprintf(fail_file,"If table is too small, then eos_extract.m used to create the degenerate offsets could result in negative degenerate offsets.  Again, either use linear interpolation or use a larger table.\n");
	      myexit(248973463);
	    }

	    testnumfscanfquantitiesdegen += 1;
	  }

	  for(ppp=0;ppp<numeosdegenquantitiesinstandardlist[tableiter];ppp++){
	    // jjj==0 already
	    set_arrays_eostable(ISDEGENTABLE,tableiter,mmm,lll,kkk,jjj,iii,ppp,degenvaluetemp[ppp]);
	  }
	}
	else{
	  // this is just for checks, store back into simple array the degen table
	  for(ppp=0;ppp<numeosdegenquantitiesinstandardlist[tableiter];ppp++){
	    get_arrays_eostable(ISDEGENTABLE,tableiter,mmm,lll,kkk,0,iii,ppp,&degenvaluetemp[ppp]); // must choose jjj=0 -- actually, get_arrays_eostable() handles "jjj" issue
	  }
	}


	// check that read-in quantities agrees with expected number so far:
 	if(testnumfscanfquantities!=numallquantitiesinfile[tableiter]){
	  dualfprintf(fail_file,"Total number of scanned EOS quantities=%d doesn't agree with expected amount=%d (whichdatatype=%d tableiter=%d)\n",testnumfscanfquantities,numallquantitiesinfile[tableiter],whichdatatype[tableiter],tableiter);
	  myexit(2873626);
	}

	if(jjj==0){
	  if(testnumfscanfquantitiesdegen!=numalldegenquantitiesinfile[tableiter]){
	    dualfprintf(fail_file,"Total number of scanned degen EOS quantities=%d doesn't agree with expected amount=%d\n",testnumfscanfquantitiesdegen,numalldegenquantitiesinfile[tableiter]);
	    myexit(2873627);
	  }
	}



	//////////////////////////////
	//
	// fifth, check degen table against normal table for proper use of independent variable offset
	//
	// Note that indepplusdegen was independently interpolated quantity instead of being constructed after the fact.  So this should only result in agreement to *truncation* error in the Matlab interpolation and NOT machine error.  This is the point of this check -- not to check that table is formatted correctly, but to test actual accuracy of interpolations.  That is, indep[] is from original HELM table.  So this check tests accuracy of Matlab interpolations from T->U,P,CHI,S space *and* the subtraction of degeneracy cut *and* interpolation of original u,p,chi,s from HELM table to HARM table.
	//
	//////////////////////////////
	// check for all jjj

	// Check that use of degen on grid is consistent with Matlab eos_extract.m script:
	// U = U0 + (Uin-U0)*pow( (Uout-U0)/(Uin-U0) , i/N ) for N points with "i" as the index, where i=jjj and N=EOSFULLN2, etc.
	//
	// U = indepplusdegen[UTOTDIFF,PTOTDIFF,CHIDIFF,STOTDIFF] : UofUdiff, PofPdiff, CHIofCHIdiff, SofSdiff: true value on the grid with respective independent variable
	// U0 = EOSOFFSET
	// UIN = EOSIN
	// UOUT = EOSOUT
	// indep = "i/N" = lutotdiff = log10( (utot - utotoffset)./(utotin - utotoffset))./log10( (utotout - utotoffset)./(utotin - utotoffset));
	//
	// now ensure that indep[normal table UEOS,PEOS,CHIEOS,SEOS] + eosdegentable[UTOTDIFF,PTOTDIFF,CHIDIFF,STOTDIFF] = indepplusdegen[UTOTDIFF,PTOTDIFF,CHIDIFF,STOTDIFF] for degen'ed quantities
	// first and last are just read-in, while eosdegentable should be for all jjj, so use stored array for each tableiter
	for(ii=0;ii<NUMEOSDEGENQUANTITIESMEM1;ii++){
	  // assume hit jjj==0 first time in loop

	  FTYPE U0,UIN,UOUT,NN,testindepplusdegen;
	  if(utotdegencut[tableiter]>DEGENCUTLASTOLDVERSION){ // all temperature-like quantities here
	    // access read-in quantities
	    U0=degenvaluetemp[ii];
	    UIN=degenvaluetemp[ii+NUMEOSDEGENQUANTITIESMEM1];
	    UOUT=degenvaluetemp[ii+2*NUMEOSDEGENQUANTITIESMEM1];
	    NN = tablesizeexpected[tableiter][UEOS]; // UEOS,PEOS,CHIEOS,SEOS all same dimension
	    // already checked that indep is consistent, so can use indep[RHOEOS+ii] or (FTYPE)jjj/NN below for "i/N"
	    testindepplusdegen = U0 + (UIN-U0)*pow( (UOUT-U0)/(UIN-U0), (FTYPE)jjj/NN);
	  }
	  else{
	    // u = utotdiff + utotoffset 
	    U0=degenvaluetemp[ii];
	    testindepplusdegen = U0 + indep[RHOEOS+ii];
	  }

	  errordegen=fabs(testindepplusdegen-indepplusdegen[ii])/(fabs(testindepplusdegen)+(indepplusdegen[ii])+SMALL);

	  if(errordegen>=TABLETOL) numberofdegenerrors[0]++;
	  if(errordegen>=0.01) numberofdegenerrors[1]++;
	  if(errordegen>=0.1) numberofdegenerrors[2]++;
	  if(errordegen>=TABLETOLTRUNCATION) numberofdegenerrors[3]++;
	  totalnumberofdegen++;

	  if(errordegen>=TABLETOLTRUNCATION){
	    dualfprintf(fail_file,"Degen not correct: iii=%d jjj=%d kkk=%d lll=%d mmm=%d:: ii=%d :: error=%21.15g :: U0=%21.15g  UIN=%21.15g UOUT=%21.15g\n",iii,jjj,kkk,lll,mmm,ii,errordegen,U0,UIN,UOUT);
	    dualfprintf(fail_file,"file has indepplusdegen=%21.15g while testindepplusdegen=%21.15g\n",indepplusdegen[ii],testindepplusdegen);
	  }

	}






	// continue onto next row
      }// end loop over all rows
      

      ////////////////
      //
      // done reading in table from file
      //
      ////////////////
      fclose(intable);
      fclose(indegentable);
      
      
      // report read-in
      trifprintf("Done reading in EOS table: tableiter=%d of %d\n",tableiter,NUMTBLS-1);
      trifprintf("Number of degen errors: %lld %lld %lld %lld out of total=%lld\n",numberofdegenerrors[0],numberofdegenerrors[1],numberofdegenerrors[2],numberofdegenerrors[3],totalnumberofdegen);








      //////////////////////////////////////////
      //
      // convert quantities to code units
      //
      //////////////////////////////////////////
    
      for(jj=0;jj<TBLLINEARITEMS;jj++){
	if(jj!=2){ // log of base has no units conversion, but rest do
	  if(rho0unittype==0) lineartablelimits[tableiter][RHOEOS][jj]/=rhounit;
	  else lineartablelimits[tableiter][RHOEOS][jj]/=rhomassunit;

	  if(utotdegencut[tableiter]>DEGENCUTLASTOLDVERSION){ // all temperature-like quantities here
	    // For utotdegencut==2 or 3 the lineartablelimits are dimensionless as log(index) that is really meaningless
	    // essentially, for utotdegencut==2 or 3, the units are emedded in storing U0, UIN, and UOUT like quantities
	  }
	  else{
	    lineartablelimits[tableiter][UEOS][jj]/=Pressureunit;
	    lineartablelimits[tableiter][PEOS][jj]/=Pressureunit;
	    lineartablelimits[tableiter][CHIEOS][jj]/=Pressureunit;
	    lineartablelimits[tableiter][SEOS][jj]/=(1.0/pow(Lunit,3.0)); // Sden = 1/cc
	  }
	  if(whichrnpmethod[tableiter]==0) lineartablelimits[tableiter][YEEOS][jj]/=Tunit; // otherwise no conversion for dimensionless Y_e
	  if(whichynumethod[tableiter]==0) lineartablelimits[tableiter][YNUEOS][jj]/=Tunit; // otherwise no conversion for dimensionless Y_\nu
	  lineartablelimits[tableiter][HEOS][jj]/=Lunit;
	}
      }
      
      // recompute tablelimits (UPDOWN=[0,1])
      for(jj=0;jj<NUMEOSINDEPS;jj++){
	for(ii=0;ii<UPDOWN;ii++){
	  // log_base(Rin-R0), log_base(Rout-R0), which gives same result as if shifted x_in and x_out by -log_base(r_units)
	  // log(Rin-R0)/log(base) for both Rout and Rin
	  if((utotdegencut[tableiter]>DEGENCUTLASTOLDVERSION)&&(jj>=FIRSTTKLIKE && jj<=LASTTKLIKE)){
	    // for temperature-like jj's, no need to recompute since already lutotdiff="i/N" like quantities
	    tablelimits[tableiter][jj][ii]=lineartablelimits[tableiter][jj][ii];
	  }
	  else{
	    tablelimits[tableiter][jj][ii]=log10(lineartablelimits[tableiter][jj][ii]-lineartablelimits[tableiter][jj][3])*lineartablelimits[tableiter][jj][2];
	  }
	}
      }
      
      // convert log10 step to code units
      // unit division results in constant uniform offset in log-space
      // below ??[2] not ever used
      
      // for new generalized log grid, transformation is that \tilde{r} = r/r_units such that:
      // \tilde{x_in} = x_in - log_b(r_units)
      // new formula is just:
      // 
      // \tilde{r} = \tilde{r_0} + b^{\tilde{x_in} + i dx}
      // that is, dx and b are unchanged, there was only an offset of the starting position
      // shift x_in and x_out the same way is consistent with above recomputation of x_in and x_out

      // dx stays the same now [2]
      //      if(rho0unittype==0) tablelimits[tableiter][RHOEOS][2]=log10(pow(10.0,tablelimits[tableiter][RHOEOS][2])/rhounit);
      //      else tablelimits[tableiter][RHOEOS][2]=log10(pow(10.0,tablelimits[tableiter][RHOEOS][2])/rhomassunit);
      //      tablelimits[tableiter][UEOS][2]=log10(pow(10.0,tablelimits[tableiter][UEOS][2])/Pressureunit);
      //      tablelimits[tableiter][YEEOS][2]=log10(pow(10.0,tablelimits[tableiter][YEEOS][2])/Tunit);
      //      tablelimits[tableiter][YNUEOS][2]=log10(pow(10.0,tablelimits[tableiter][YNUEOS][2])/Tunit);
      //      tablelimits[tableiter][HEOS][2]=log10(pow(10.0,tablelimits[tableiter][HEOS][2])/Lunit);
      
      // just recompute the below
      for(ii=0;ii<NUMEOSINDEPS;ii++){
	tablelimits[tableiter][jj][3] = ((FTYPE)tablesize[tableiter][ii]-1.0)/(tablelimits[tableiter][ii][1] - tablelimits[tableiter][ii][0]);
      }

      // assume base stays same [4]
      // shouldn't need log10 version of offset, so assume not changed or needed anymore

      // now 0-3 (and 4-5) of tablelimits set, lineartablelimits set, tablesize is same
      
      
      // output to logfull_file so have information
      trifprintf("Table information in code units:\n");
      for(jj=0;jj<NUMEOSINDEPS;jj++){
	trifprintf("tablesize[%d][%d]=%d\n",tableiter,jj,tablesize[tableiter][jj]);
	for(ii=0;ii<TBLLINEARITEMS;ii++) trifprintf("lineartablelimits[%d][%d][%d]=%21.15g\n",tableiter,jj,ii,lineartablelimits[tableiter][jj][ii]);
	for(ii=0;ii<TBLITEMS;ii++) trifprintf("tablelimits[%d][%d][%d]=%21.15g\n",tableiter,jj,ii,tablelimits[tableiter][jj][ii]);
      }
      
      // report read-in
      trifprintf("Done converting units for header information: tableiter=%d of %d\n",tableiter,NUMTBLS-1);
      
          


      
      ///////////////////////////////
      //
      // convert eostable units
      //
      ///////////////////////////////

#define EOSCONVERTLOOP for(mmm=0;mmm<tablesize[tableiter][vartypearray[5]];mmm++)for(lll=0;lll<tablesize[tableiter][vartypearray[4]];lll++)for(kkk=0;kkk<tablesize[tableiter][vartypearray[3]];kkk++)for(jjj=0;jjj<tablesize[tableiter][vartypearray[2]];jjj++)for(iii=0;iii<tablesize[tableiter][vartypearray[1]];iii++)

      // UEOS gives same size as PEOS, CHIEOS, and SEOS
      EOSCONVERTLOOP{

	////////////////
	//
	// temp set
	//
	///////////////
	
	// normal table
	for(jj=0;jj<numeosquantitiesinstandardlist[tableiter];jj++) get_arrays_eostable(ISNOTDEGENTABLE,tableiter,mmm,lll,kkk,jjj,iii,jj,&tabletemp[jj]);
	// degen table
	if(jjj==0) for(jj=0;jj<numeosdegenquantitiesinstandardlist[tableiter];jj++) get_arrays_eostable(ISDEGENTABLE,tableiter,mmm,lll,kkk,jjj,iii,jj,&degentabletemp[jj]);



	// note that we process all entries for "in" in tabletemp since doesn't matter if certain table has no such quantity since will ignore irrelevant converted value later.  Just keeps code simple and avoiding conditionals

	tabletemp[PofRHOUin]/=Pressureunit;
	tabletemp[UofRHOPin]/=Pressureunit;
	// dPdRHO0ofRHOU dimensionless
	// dPdUofRHOU dimensionless

	tabletemp[UofRHOSin]/=Pressureunit;

	tabletemp[CS2ofRHOUin]/=(Vunit*Vunit);


	// entropy density (erg/K/cc)
	// kb doesn't convert units, but changes kb T to erg
	// presumes entropy is used with energy as in first-law: dQ = (kbT)dS where kbT is in ergs
	// previously had units of entropy density as erg/K/cc, but now units read-in are just 1/cc to avoid use of pointless kb
	tabletemp[SofRHOUin]/=(1.0/pow(Lunit,3.0));
	// Note that often people will plot "entropy per baryon" that is "SofRHOU"/(\rho/m_b) that is dimensionless entropy per baryon
	// From HARM quantities, convert back to cgs and then compute above
	// Note that in HARM, we *define* $Sden = Ss/(\rho_0 c^2)$ with the extra $m_b c^2$ baggage, so units should account for that.

	// below is (1/cc) / (erg/cc) \propto 1/erg since we *input* rho as rho c^2
	tabletemp[DSDRHOofRHOUin]/=(1.0/energyunit); 
	tabletemp[DSDUofRHOUin]/=(1.0/energyunit);

	// SSofRHOCHI is nearly dimensionless (really has units of 1/(mb c^2) since obtained by division of Sden by (\rho_0 c^2) instead of n_b)
	tabletemp[SSofRHOCHIin]/=(1.0/energyunit);

	// DSSDRHOofRHOCHI, DSSDCHIofRHOCHI have units of Ss/(rho0 c^2) since Ss has units of 1/(mb c^2) and input rho with units of rho*c^2
	tabletemp[DSSDRHOofRHOCHIin]/=(1.0/(energyunit*Pressureunit));
	tabletemp[DSSDCHIofRHOCHIin]/=(1.0/(energyunit*Pressureunit));

	
	tabletemp[PofRHOCHIin]/=Pressureunit;
	// IDRHO0DP is dimensionless
	// IDCHIDP is dimensionless


	// *always* have temperature in table
	// TEMP used for diagnostics, not used otherwise
	tabletemp[TEMPUin]/=Tempunit;
	tabletemp[TEMPPin]/=Tempunit;
	tabletemp[TEMPCHIin]/=Tempunit;
	tabletemp[TEMPSin]/=Tempunit;
	  
	


	////////////////////////////////
	//
	// deal with extra quantities
	if(whichdatatype[tableiter]==1){
	  // Qm is in erg/s/cm^3 (Qvol, not Qsurf)
	  // this is divided by H when used as a volume rate
	  tabletemp[EXTRA1in]/=(edotunit/(Lunit*Lunit*Lunit));
	}
	else if(whichdatatype[tableiter]==2){
	  // \tau/H
	  tabletemp[EXTRA1in]/=(1.0/Lunit);
	  tabletemp[EXTRA2in]/=(1.0/Lunit);
	  tabletemp[EXTRA3in]/=(1.0/Lunit);
	  tabletemp[EXTRA4in]/=(1.0/Lunit);
	  tabletemp[EXTRA5in]/=(1.0/Lunit);
	  tabletemp[EXTRA6in]/=(1.0/Lunit);
	  tabletemp[EXTRA7in]/=(1.0/Lunit);
	  tabletemp[EXTRA8in]/=(1.0/Lunit);
	  tabletemp[EXTRA9in]/=(1.0/Lunit);
	  tabletemp[EXTRA10in]/=(1.0/Lunit);
	  tabletemp[EXTRA11in]/=(1.0/Lunit);
	  tabletemp[EXTRA12in]/=(1.0/Lunit);
	  //	  \Gamma = 1/s
	  tabletemp[EXTRA13in]/=(1.0/Tunit);
	  tabletemp[EXTRA14in]/=(1.0/Tunit);
	  tabletemp[EXTRA15in]/=(1.0/Tunit);
	  tabletemp[EXTRA16in]/=(1.0/Tunit);
	}
	else if(whichdatatype[tableiter]==3){
	  // Qphoton=erg/s/cm^3
	  tabletemp[EXTRA1in]/=(edotunit/(Lunit*Lunit*Lunit));

	  // Qm=erg/s/cm^3
	  tabletemp[EXTRA2in]/=(edotunit/(Lunit*Lunit*Lunit));
	  
	  // graddotrhouyl=\rho/sec = m_b/s/cm^3
	  // \nabla_\mu (\rho_0 u^\mu Y_e) =  (m_b/H) (\dot{N}_{\bar{\nu}_e} - \dot{N}_{\nu_e}) 
	  if(rho0unittype==0) tabletemp[EXTRA3in]/=(rhounit/Tunit);
	  else tabletemp[EXTRA3in]/=(rhomassunit/Tunit);
	  
	  // Tthermaltot
	  tabletemp[EXTRA4in]/=Tunit;
	  tabletemp[EXTRA5in]/=Tunit;

	  // lambdatot = mean free path such that H = [\int (dr/\lambda)] / [1/\lambda]
	  // and \tau = \int dr/\lambda
	  tabletemp[EXTRA6in]/=Lunit;
	  tabletemp[EXTRA7in]/=Lunit;

	  // Enuglobal : erg
	  tabletemp[EXTRA8in]/=energyunit;
	  tabletemp[EXTRA9in]/=energyunit;
	  tabletemp[EXTRA10in]/=energyunit;

	  // no conversion for Ynuthermal that's dimensionless EXTRA11
	}
	else if(whichdatatype[tableiter]==4){
	  // \tau/H
	  tabletemp[EXTRA1in]/=(1.0/Lunit);
	  tabletemp[EXTRA2in]/=(1.0/Lunit);
	  tabletemp[EXTRA3in]/=(1.0/Lunit);
	  tabletemp[EXTRA4in]/=(1.0/Lunit);
	  tabletemp[EXTRA5in]/=(1.0/Lunit);
	  tabletemp[EXTRA6in]/=(1.0/Lunit);
	  tabletemp[EXTRA7in]/=(1.0/Lunit);
	  tabletemp[EXTRA8in]/=(1.0/Lunit);
	  tabletemp[EXTRA9in]/=(1.0/Lunit);
	  tabletemp[EXTRA10in]/=(1.0/Lunit);
	  tabletemp[EXTRA11in]/=(1.0/Lunit);
	  tabletemp[EXTRA12in]/=(1.0/Lunit);
	  // energy densities
	  tabletemp[EXTRA13in]/=(Pressureunit);
	  tabletemp[EXTRA14in]/=(Pressureunit);
	  tabletemp[EXTRA15in]/=(Pressureunit);
	  // number densities
	  tabletemp[EXTRA16in]/=(1.0/pow(Lunit,3.0));
	  tabletemp[EXTRA17in]/=(1.0/pow(Lunit,3.0));
	  tabletemp[EXTRA18in]/=(1.0/pow(Lunit,3.0));
	  // mean free paths (length)
	  tabletemp[EXTRA19in]/=(Lunit);
	  tabletemp[EXTRA20in]/=(Lunit);
	  // \tau/H for photons
	  tabletemp[EXTRA21in]/=(1.0/Lunit);
	  tabletemp[EXTRA22in]/=(1.0/Lunit);

	  // DEBUG: log(fun) \propto A+B*log(rho) -> fun = 10^(A+B*log(rho)) = Q*rho^B
	  // rho = exp(A+B*iii)
	  //	  FTYPE fakerho = exp(0.0+1.0*iii);
	  // tabletemp[EXTRA22in]=1.0*pow(fakerho,3.0);
	  //dualfprintf(fail_file,"GOD: iii=%d %21.15g\n",iii,tabletemp[EXTRA22in]);

	  // thermal number densities
	  tabletemp[EXTRA23in]/=(1.0/pow(Lunit,3.0));
	  tabletemp[EXTRA24in]/=(1.0/pow(Lunit,3.0));
	}

	


	////////////////////
	//
	// convert units for degen table
	//
	////////////////////
	if(jjj==0){

	  degentabletemp[UTOTOFFSETin]/=Pressureunit; // pressure units
	  degentabletemp[PTOTOFFSETin]/=Pressureunit; // pressure units
	  degentabletemp[CHIOFFSETin]/=Pressureunit; // pressure units
	  degentabletemp[STOTOFFSETin]/=(1.0/pow(Lunit,3.0)); // 1/cc

	  if(utotdegencut[tableiter]>DEGENCUTLASTOLDVERSION){
	    degentabletemp[UTOTINin]/=Pressureunit; // pressure units
	    degentabletemp[PTOTINin]/=Pressureunit; // pressure units
	    degentabletemp[CHIINin]/=Pressureunit; // pressure units
	    degentabletemp[STOTINin]/=(1.0/pow(Lunit,3.0)); // 1/cc

	    degentabletemp[UTOTOUTin]/=Pressureunit; // pressure units
	    degentabletemp[PTOTOUTin]/=Pressureunit; // pressure units
	    degentabletemp[CHIOUTin]/=Pressureunit; // pressure units
	    degentabletemp[STOTOUTin]/=(1.0/pow(Lunit,3.0)); // 1/cc
	  }

	  if(numeosdegenquantitiesinstandardlist[tableiter]!=NUMEOSDEGENQUANTITIESin){
	    dualfprintf(fail_file,"Degen table: expected %d but got %d\n",NUMEOSDEGENQUANTITIESin,numeosdegenquantitiesinstandardlist[tableiter]);
	    myexit(13905826);
	  }
	}
	

	////////////////
	//
	// temp restore back into arrays
	//
	///////////////

	// normal table
	for(jj=0;jj<numeosquantitiesinstandardlist[tableiter];jj++) set_arrays_eostable(ISNOTDEGENTABLE,tableiter,mmm,lll,kkk,jjj,iii,jj,tabletemp[jj]);
	// degen table
	if(jjj==0) for(jj=0;jj<numeosdegenquantitiesinstandardlist[tableiter];jj++) set_arrays_eostable(ISDEGENTABLE,tableiter,mmm,lll,kkk,jjj,iii,jj,degentabletemp[jj]);



      } // end loop over table

            
      // report conversion
      trifprintf("Done converting units of EOS table: tableiter=%d of %d\n",tableiter,NUMTBLS-1);

    }// end loop over tables


    // compute code version of invalid temperature
    invalidtempcode=INVALIDTEMP/Tempunit;


    trifprintf("Sden_convertfactor=%21.15g\n",(1.0/pow(Lunit,3.0)));

  } // end if myid==0
  


  ///////////
  //
  // tell code that did setup Kaz EOS
  //
  ////////////
  didsetupkazeos=1;


  ////////////////
  //
  // broadcast global variables to other CPUs in case only created by myid==0
  //
  ///////////////
  bcast_kazeos();


  ////////////////
  //
  // Report that done with setting up EOS
  //
  ///////////////

  trifprintf("Done with reading in EOS tables\n");
  dualfprintf(log_file,"proc: %d: Done reading in EOS table\n",myid);
  
  
  
  
}





static void bcast_kazeos(void)
{





  //////////////////////////////////////
  //
  // send data to all CPUs
  // so all CPUs have parameters and full eostable data
  //
  // Anything in kazfulleos.defs.h and kazfulleos.eostablesdefs.h should be Bcast()'ed here!
  //
  ///////////////////////////////////////
#if(USEMPI)

  /////////////////////////////
  //
  // kazfulleos.defs.h globals:
  //
  // regexp from kazfulleos.defs.h:
  // static \([_a-zA-Z]+\) \([_a-zA-Z0-9]+\)\[\(.*\)\]; -> MPI_Bcast(&\2[0],\3,MPI_\1,MPIid[0],MPI_COMM_GRMHD);
  // THEN:
  // 1) raise (e.g.) MPI_int -> MPI_INT
  // 2) replace ]['s with * in data size and add [0]'s to &array as required (same number of ]['s replaced)
  //
  /////////////////////////////

  MPI_Bcast(&numeosquantitiesinstandardlist[0],NUMTBLS,MPI_INT,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&numeosdegenquantitiesinstandardlist[0],NUMTBLS,MPI_INT,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&numeosquantitiesinfile[0],NUMTBLS,MPI_INT,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&numeosdegenquantitiesinfile[0],NUMTBLS,MPI_INT,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&numallquantitiesinfile[0],NUMTBLS,MPI_INT,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&numalldegenquantitiesinfile[0],NUMTBLS,MPI_INT,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&extralimits[0][0],NUMTBLS*2,MPI_INT,MPIid[0],MPI_COMM_GRMHD);

  MPI_Bcast(&inputtablelimits[0][0][0],NUMTBLS*NUMEOSINDEPS*TBLITEMS,MPI_FTYPE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&tablelimits[0][0][0],NUMTBLS*NUMEOSINDEPS*TBLITEMS,MPI_FTYPE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&lineartablelimits[0][0][0],NUMTBLS*NUMEOSINDEPS*TBLLINEARITEMS,MPI_FTYPE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&tablesize[0][0],NUMTBLS*NUMEOSINDEPS,MPI_INT,MPIid[0],MPI_COMM_GRMHD);

  MPI_Bcast(&vartypeeosextraarray[0],NUMINDEPDIMENSMEM,MPI_INT,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&vartypearray[0],NUMINDEPDIMENSMEM,MPI_INT,MPIid[0],MPI_COMM_GRMHD);

  MPI_Bcast(&vartypeheightarray[0],NUMHDIRECTIONS+1,MPI_INT,MPIid[0],MPI_COMM_GRMHD);

  MPI_Bcast(&vardegentypearray[0],NUMEOSDEGENINDEPS+1,MPI_INT,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&varnormalcompare2degentypearray[0],NUMEOSDEGENINDEPS+1,MPI_INT,MPIid[0],MPI_COMM_GRMHD);


  MPI_Bcast(&numcolintablesubtype[0],NUMTABLESUBTYPES,MPI_INT,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&whichdintablesubtype[0],NUMTABLESUBTYPES,MPI_INT,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&firsteosintablesubtype[0],NUMTABLESUBTYPES,MPI_INT,MPIid[0],MPI_COMM_GRMHD);

  MPI_Bcast(&whichtablesubtypeinquantity[0],NUMEOSQUANTITIESMEM,MPI_INT,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&whichcolinquantity[0],NUMEOSQUANTITIESMEM,MPI_INT,MPIid[0],MPI_COMM_GRMHD);

  MPI_Bcast(&dologinterp_sub_coli[0][0],NUMTABLESUBTYPES*MAXEOSPIPELINE,MPI_INT,MPIid[0],MPI_COMM_GRMHD);

  MPI_Bcast(&invalidtempcode,1,MPI_FTYPE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&invalidlogtempcode,1,MPI_FTYPE,MPIid[0],MPI_COMM_GRMHD);

  MPI_Bcast(&whichrnpmethod[0],NUMTBLS,MPI_INT,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&whichynumethod[0],NUMTBLS,MPI_INT,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&whichhcmmethod[0],NUMTBLS,MPI_INT,MPIid[0],MPI_COMM_GRMHD);

  MPI_Bcast(&whichdatatype[0],NUMTBLS,MPI_INT,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&utotdegencut[0],NUMTBLS,MPI_INT,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&numc[0],NUMTBLS,MPI_INT,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&numextras[0],NUMTBLS,MPI_INT,MPIid[0],MPI_COMM_GRMHD);

  MPI_Bcast(&primarytable,1,MPI_INT,MPIid[0],MPI_COMM_GRMHD);


  MPI_Bcast(&FAKE2IDEALNUCLEAROFFSET,1,MPI_FTYPE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&TRUENUCLEAROFFSET,1,MPI_FTYPE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&DEGENNUCLEAROFFSET,1,MPI_FTYPE,MPIid[0],MPI_COMM_GRMHD);

  MPI_Bcast(&lsoffset,1,MPI_FTYPE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&fakelsoffset,1,MPI_FTYPE,MPIid[0],MPI_COMM_GRMHD);

  MPI_Bcast(&didsetupkazeos,1,MPI_INT,MPIid[0],MPI_COMM_GRMHD);




  /////////////////////////////
  //
  // kazfulleos.eostablesdefs.h globals:
  //
  /////////////////////////////


  // Bcast table data
  // eostable is double, so MPI_DOUBLE
  // to generate below, take part of kazfull.eostablesdefs.h within WHICHEOS==KAZFULL and regexp:
  // double BASEEOSMAC(\([a-zA-Z0-9]+\),\([a-zA-Z0-9]+\),\([a-zA-Z0-9]+\),\([a-zA-Z0-9]+\),\([a-zA-Z0-9]+\),\([a-zA-Z0-9]+\),\([a-zA-Z0-9]+\),\([a-zA-Z0-9]+\)); ->
  // MPI_Bcase(&(BASEEOSMAC(\1,0,0,0,0,0,0,0)),\2*\3*\4*\5*\6*\7*\8,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  // NOTE!  We are using BASEEOSMAC() so can just use BASEEOSMAC(name,0,0,0,0,0,0,0) instead of having to add shift

#if(WHICHEOS==KAZFULL)

  // full
  MPI_Bcase(&(BASEEOSMAC(eosfulltabledegen,0,0,0,0,0,0,0)),NUMEOSDEGENQUANTITIESMEM1*EOSFULLDEGENN5*EOSFULLDEGENN4*EOSFULLDEGENN3*EOSFULLDEGENN2*EOSFULLDEGENN1*NUMEOSDEGENQUANTITIESMEM2,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcase(&(BASEEOSMAC(eosfulltablestandard,0,0,0,0,0,0,0)),1*EOSFULLN5*EOSFULLN4*EOSFULLN3*EOSFULLN2*EOSFULLN1*NUMEOSSTANDARDQUANTITIESMEM,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcase(&(BASEEOSMAC(eosfulltableguess,0,0,0,0,0,0,0)),1*EOSFULLN5*EOSFULLN4*EOSFULLN3*EOSFULLN2*EOSFULLN1*NUMEOSGUESSQUANTITIESMEM,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcase(&(BASEEOSMAC(eosfulltablediss,0,0,0,0,0,0,0)),1*EOSFULLN5*EOSFULLN4*EOSFULLN3*EOSFULLN2*EOSFULLN1*NUMEOSDISSQUANTITIESMEM,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcase(&(BASEEOSMAC(eosfulltabledp,0,0,0,0,0,0,0)),1*EOSFULLN5*EOSFULLN4*EOSFULLN3*EOSFULLN2*EOSFULLN1*NUMEOSDPQUANTITIESMEM,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcase(&(BASEEOSMAC(eosfulltablesden,0,0,0,0,0,0,0)),1*EOSFULLN5*EOSFULLN4*EOSFULLN3*EOSFULLN2*EOSFULLN1*NUMEOSSDENQUANTITIESMEM,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcase(&(BASEEOSMAC(eosfulltablesspec,0,0,0,0,0,0,0)),1*EOSFULLN5*EOSFULLN4*EOSFULLN3*EOSFULLN2*EOSFULLN1*NUMEOSSSPECQUANTITIESMEM,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcase(&(BASEEOSMAC(eosfulltablepofchi,0,0,0,0,0,0,0)),1*EOSFULLN5*EOSFULLN4*EOSFULLN3*EOSFULLN2*EOSFULLN1*NUMEOSPOFCHIQUANTITIESMEM,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcase(&(BASEEOSMAC(eosfulltabletemp,0,0,0,0,0,0,0)),NUMEOSTEMPQUANTITIESMEM1*EOSFULLN5*EOSFULLN4*EOSFULLN3*EOSFULLN2*EOSFULLN1*NUMEOSTEMPQUANTITIESMEM2,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);

  // different size for extra table
  MPI_Bcase(&(BASEEOSMAC(eosfulltableextra,0,0,0,0,0,0,0)),1*EOSEXTRAFULLN5*EOSEXTRAFULLN4*EOSEXTRAFULLN3*EOSEXTRAFULLN2*EOSEXTRAFULLN1*NUMEOSEXTRAQUANTITIESMEM,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
#if(WHICHDATATYPEGENERAL==4)
  MPI_Bcase(&(BASEEOSMAC(eosfulltableextradegen,0,0,0,0,0,0,0)),NUMEOSDEGENQUANTITIESMEM1*EOSEXTRAFULLDEGENN5*EOSEXTRAFULLDEGENN4*EOSEXTRAFULLDEGENN3*EOSEXTRAFULLDEGENN2*EOSEXTRAFULLDEGENN1*NUMEOSDEGENQUANTITIESMEM2,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcase(&(BASEEOSMAC(eosfulltableextratemp,0,0,0,0,0,0,0)),NUMEOSTEMPQUANTITIESMEM1*EOSEXTRAFULLTEMPN5*EOSEXTRAFULLTEMPN4*EOSEXTRAFULLTEMPN3*EOSEXTRAFULLTEMPN2*EOSEXTRAFULLTEMPN1*NUMEOSTEMPQUANTITIESMEM2,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
#endif

  // simple
  MPI_Bcase(&(BASEEOSMAC(eossimpletabledegen,0,0,0,0,0,0,0)),NUMEOSDEGENQUANTITIESMEM1*EOSSIMPLEDEGENN5*EOSSIMPLEDEGENN4*EOSSIMPLEDEGENN3*EOSSIMPLEDEGENN2*EOSSIMPLEDEGENN1*NUMEOSDEGENQUANTITIESMEM2,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcase(&(BASEEOSMAC(eossimpletablestandard,0,0,0,0,0,0,0)),1*EOSSIMPLEN5*EOSSIMPLEN4*EOSSIMPLEN3*EOSSIMPLEN2*EOSSIMPLEN1*NUMEOSSTANDARDQUANTITIESMEM,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcase(&(BASEEOSMAC(eossimpletableguess,0,0,0,0,0,0,0)),1*EOSSIMPLEN5*EOSSIMPLEN4*EOSSIMPLEN3*EOSSIMPLEN2*EOSSIMPLEN1*NUMEOSGUESSQUANTITIESMEM,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcase(&(BASEEOSMAC(eossimpletablediss,0,0,0,0,0,0,0)),1*EOSSIMPLEN5*EOSSIMPLEN4*EOSSIMPLEN3*EOSSIMPLEN2*EOSSIMPLEN1*NUMEOSDISSQUANTITIESMEM,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcase(&(BASEEOSMAC(eossimpletabledp,0,0,0,0,0,0,0)),1*EOSSIMPLEN5*EOSSIMPLEN4*EOSSIMPLEN3*EOSSIMPLEN2*EOSSIMPLEN1*NUMEOSDPQUANTITIESMEM,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcase(&(BASEEOSMAC(eossimpletablesden,0,0,0,0,0,0,0)),1*EOSSIMPLEN5*EOSSIMPLEN4*EOSSIMPLEN3*EOSSIMPLEN2*EOSSIMPLEN1*NUMEOSSDENQUANTITIESMEM,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcase(&(BASEEOSMAC(eossimpletablesspec,0,0,0,0,0,0,0)),1*EOSSIMPLEN5*EOSSIMPLEN4*EOSSIMPLEN3*EOSSIMPLEN2*EOSSIMPLEN1*NUMEOSSSPECQUANTITIESMEM,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcase(&(BASEEOSMAC(eossimpletablepofchi,0,0,0,0,0,0,0)),1*EOSSIMPLEN5*EOSSIMPLEN4*EOSSIMPLEN3*EOSSIMPLEN2*EOSSIMPLEN1*NUMEOSPOFCHIQUANTITIESMEM,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcase(&(BASEEOSMAC(eossimpletabletemp,0,0,0,0,0,0,0)),NUMEOSTEMPQUANTITIESMEM1*EOSSIMPLEN5*EOSSIMPLEN4*EOSSIMPLEN3*EOSSIMPLEN2*EOSSIMPLEN1*NUMEOSTEMPQUANTITIESMEM2,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);

  // different size for extra table
  MPI_Bcase(&(BASEEOSMAC(eossimpletableextra,0,0,0,0,0,0,0)),1*EOSEXTRASIMPLEN5*EOSEXTRASIMPLEN4*EOSEXTRASIMPLEN3*EOSEXTRASIMPLEN2*EOSEXTRASIMPLEN1*NUMEOSEXTRAQUANTITIESMEM,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
#if(WHICHDATATYPEGENERAL==4)
  MPI_Bcase(&(BASEEOSMAC(eossimpletableextradegen,0,0,0,0,0,0,0)),NUMEOSDEGENQUANTITIESMEM1*EOSEXTRASIMPLEDEGENN5*EOSEXTRASIMPLEDEGENN4*EOSEXTRASIMPLEDEGENN3*EOSEXTRASIMPLEDEGENN2*EOSEXTRASIMPLEDEGENN1*NUMEOSDEGENQUANTITIESMEM2,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcase(&(BASEEOSMAC(eossimpletableextratemp,0,0,0,0,0,0,0)),NUMEOSTEMPQUANTITIESMEM1*EOSEXTRASIMPLETEMPN5*EOSEXTRASIMPLETEMPN4*EOSEXTRASIMPLETEMPN3*EOSEXTRASIMPLETEMPN2*EOSEXTRASIMPLETEMPN1*NUMEOSTEMPQUANTITIESMEM2,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
#endif



  // simple zoom
  MPI_Bcase(&(BASEEOSMAC(eossimplezoomtabledegen,0,0,0,0,0,0,0)),NUMEOSDEGENQUANTITIESMEM1*EOSSIMPLEZOOMDEGENN5*EOSSIMPLEZOOMDEGENN4*EOSSIMPLEZOOMDEGENN3*EOSSIMPLEZOOMDEGENN2*EOSSIMPLEZOOMDEGENN1*NUMEOSDEGENQUANTITIESMEM2,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcase(&(BASEEOSMAC(eossimplezoomtablestandard,0,0,0,0,0,0,0)),1*EOSSIMPLEZOOMN5*EOSSIMPLEZOOMN4*EOSSIMPLEZOOMN3*EOSSIMPLEZOOMN2*EOSSIMPLEZOOMN1*NUMEOSSTANDARDQUANTITIESMEM,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcase(&(BASEEOSMAC(eossimplezoomtableguess,0,0,0,0,0,0,0)),1*EOSSIMPLEN5*EOSSIMPLEN4*EOSSIMPLEN3*EOSSIMPLEN2*EOSSIMPLEN1*NUMEOSGUESSQUANTITIESMEM,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcase(&(BASEEOSMAC(eossimplezoomtablediss,0,0,0,0,0,0,0)),1*EOSSIMPLEN5*EOSSIMPLEN4*EOSSIMPLEN3*EOSSIMPLEN2*EOSSIMPLEN1*NUMEOSDISSQUANTITIESMEM,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcase(&(BASEEOSMAC(eossimplezoomtabledp,0,0,0,0,0,0,0)),1*EOSSIMPLEZOOMN5*EOSSIMPLEZOOMN4*EOSSIMPLEZOOMN3*EOSSIMPLEZOOMN2*EOSSIMPLEZOOMN1*NUMEOSDPQUANTITIESMEM,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcase(&(BASEEOSMAC(eossimplezoomtablesden,0,0,0,0,0,0,0)),1*EOSSIMPLEZOOMN5*EOSSIMPLEZOOMN4*EOSSIMPLEZOOMN3*EOSSIMPLEZOOMN2*EOSSIMPLEZOOMN1*NUMEOSSDENQUANTITIESMEM,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcase(&(BASEEOSMAC(eossimplezoomtablesspec,0,0,0,0,0,0,0)),1*EOSSIMPLEZOOMN5*EOSSIMPLEZOOMN4*EOSSIMPLEZOOMN3*EOSSIMPLEZOOMN2*EOSSIMPLEZOOMN1*NUMEOSSSPECQUANTITIESMEM,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcase(&(BASEEOSMAC(eossimplezoomtablepofchi,0,0,0,0,0,0,0)),1*EOSSIMPLEZOOMN5*EOSSIMPLEZOOMN4*EOSSIMPLEZOOMN3*EOSSIMPLEZOOMN2*EOSSIMPLEZOOMN1*NUMEOSPOFCHIQUANTITIESMEM,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcase(&(BASEEOSMAC(eossimplezoomtabletemp,0,0,0,0,0,0,0)),NUMEOSTEMPQUANTITIESMEM1*EOSSIMPLEZOOMN5*EOSSIMPLEZOOMN4*EOSSIMPLEZOOMN3*EOSSIMPLEZOOMN2*EOSSIMPLEZOOMN1*NUMEOSTEMPQUANTITIESMEM2,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
 
  // different size for extra table
  MPI_Bcase(&(BASEEOSMAC(eossimplezoomtableextra,0,0,0,0,0,0,0)),1*EOSEXTRASIMPLEZOOMN5*EOSEXTRASIMPLEZOOMN4*EOSEXTRASIMPLEZOOMN3*EOSEXTRASIMPLEZOOMN2*EOSEXTRASIMPLEZOOMN1*NUMEOSEXTRAQUANTITIESMEM,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
#if(WHICHDATATYPEGENERAL==4)
  MPI_Bcase(&(BASEEOSMAC(eossimplezoomtableextradegen,0,0,0,0,0,0,0)),NUMEOSDEGENQUANTITIESMEM1*EOSEXTRASIMPLEZOOMDEGENN5*EOSEXTRASIMPLEZOOMDEGENN4*EOSEXTRASIMPLEZOOMDEGENN3*EOSEXTRASIMPLEZOOMDEGENN2*EOSEXTRASIMPLEZOOMDEGENN1*NUMEOSDEGENQUANTITIESMEM2,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcase(&(BASEEOSMAC(eossimplezoomtableextratemp,0,0,0,0,0,0,0)),NUMEOSTEMPQUANTITIESMEM1*EOSEXTRASIMPLEZOOMTEMPN5*EOSEXTRASIMPLEZOOMTEMPN4*EOSEXTRASIMPLEZOOMTEMPN3*EOSEXTRASIMPLEZOOMTEMPN2*EOSEXTRASIMPLEZOOMTEMPN1*NUMEOSTEMPQUANTITIESMEM2,MPI_DOUBLE,MPIid[0],MPI_COMM_GRMHD);
#endif

#endif // end if(WHICHEOS==KAZFULL)



#endif // end if USEMPI
  




}






// translate input EOS table columns into HARM EOS arrays
// input jjj should be 0 if inputting degen table data
// only used by read_setup_eostable()
// for simple and simplezoom, just use full block and replace "fulltable" -> "simpletable"
// note that name of array for extra stuff is "tableextra" and "tableextradegen", so the "extra" is at end of name rather than as "eosextra..." as in macro names
// no longer refer to macro label "inextra" since put into canonical order and position before calling set_ or get_ functions
static void set_arrays_eostable(int whichdegen, int whichtable, int mmm, int lll, int kkk, int jjj, int iii, int incol, double value)
{


  // overrides:
  // then modify indices since different sub tables have different dimensions.  This will cause repeated-reads to put into same locations in memory.
  if(WHICHDATATYPEGENERAL==4){
    if(
       ((whichtable==FULLTABLE || whichtable==SIMPLETABLE || whichtable==SIMPLEZOOMTABLE) && (incol<FIRSTEXTRAin && incol>LASTEXTRAin)) || // for all non-extras in normal table
       ((whichtable==EXTRAFULLTABLE || whichtable==EXTRASIMPLETABLE || whichtable==EXTRASIMPLEZOOMTABLE) && (incol<FIRSTEXTRAin && incol>LASTEXTRAin)) || // for temperature in extra table
       (whichdegen==1)
       ){
      // then table not storing mmm or lll, so set to zero
      mmm=lll=0;
    }
  }
  if(whichdegen){
    jjj=0; // in case not already set
  }


  if(0){
  }
#if(ALLOWFULLTABLE==1)
  else if(whichtable==FULLTABLE || whichtable==EXTRAFULLTABLE){
    if(whichdegen==ISNOTDEGENTABLE){

      if(whichtable==FULLTABLE){
	if(incol==PofRHOUin) EOSMAC(eosfulltablestandard,0,mmm,lll,kkk,jjj,iii,PofRHOU)=value;
	if(incol==CS2ofRHOUin) EOSMAC(eosfulltablestandard,0,mmm,lll,kkk,jjj,iii,CS2ofRHOU)=value;

	if(incol==UofRHOPin) EOSMAC(eosfulltableguess,0,mmm,lll,kkk,jjj,iii,UofRHOP)=value;

	if(incol==UofRHOSin) EOSMAC(eosfulltablediss,0,mmm,lll,kkk,jjj,iii,UofRHOS)=value;

	if(incol==DPDRHOofRHOUin) EOSMAC(eosfulltabledp,0,mmm,lll,kkk,jjj,iii,DPDRHOofRHOU)=value;
	if(incol==DPDUofRHOUin) EOSMAC(eosfulltabledp,0,mmm,lll,kkk,jjj,iii,DPDUofRHOU)=value;

	if(incol==SofRHOUin) EOSMAC(eosfulltablesden,0,mmm,lll,kkk,jjj,iii,SofRHOU)=value;
	if(incol==DSDRHOofRHOUin) EOSMAC(eosfulltablesden,0,mmm,lll,kkk,jjj,iii,DSDRHOofRHOU)=value;
	if(incol==DSDUofRHOUin) EOSMAC(eosfulltablesden,0,mmm,lll,kkk,jjj,iii,DSDUofRHOU)=value;

	if(incol==SSofRHOCHIin) EOSMAC(eosfulltablesspec,0,mmm,lll,kkk,jjj,iii,SSofRHOCHI)=value;
	if(incol==DSSDRHOofRHOCHIin) EOSMAC(eosfulltablesspec,0,mmm,lll,kkk,jjj,iii,DSSDRHOofRHOCHI)=value;
	if(incol==DSSDCHIofRHOCHIin) EOSMAC(eosfulltablesspec,0,mmm,lll,kkk,jjj,iii,DSSDCHIofRHOCHI)=value;

	if(incol==PofRHOCHIin) EOSMAC(eosfulltablepofchi,0,mmm,lll,kkk,jjj,iii,PofRHOCHI)=value;
	if(incol==IDRHO0DPin) EOSMAC(eosfulltablepofchi,0,mmm,lll,kkk,jjj,iii,IDRHO0DP)=value;
	if(incol==IDCHIDPin) EOSMAC(eosfulltablepofchi,0,mmm,lll,kkk,jjj,iii,IDCHIDP)=value;

	if(incol==TEMPUin) EOSMAC(eosfulltabletemp,UTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN)=value;
	if(incol==TEMPPin) EOSMAC(eosfulltabletemp,PTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN)=value;
	if(incol==TEMPCHIin) EOSMAC(eosfulltabletemp,CHIDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN)=value;
	if(incol==TEMPSin) EOSMAC(eosfulltabletemp,STOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN)=value;

	// note that if WHICHDATATYPEGENERAL==4, below is just not used
	if(incol>=FIRSTEXTRAin && incol<=LASTEXTRAin) EOSMAC(eosfulltableextra,0,mmm,lll,kkk,jjj,iii,incol-FIRSTEXTRAin+FIRSTEOSEXTRA)=value; // assumes extra's are ordered in sequence
      }
      else{
	// incol is different then since extras start with 4 temperature quantities and then rest of normal extras
	if(incol>=FIRSTEXTRAin && incol<=LASTEXTRAin) EOSMAC(eosfulltableextra,0,mmm,lll,kkk,jjj,iii,incol-FIRSTEXTRAin+FIRSTEOSEXTRA)=value; // assumes extra's are ordered in sequence
	if(incol==TEMPUin) EOSMAC(eosfulltableextratemp,UTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN)=value;
	if(incol==TEMPPin) EOSMAC(eosfulltableextratemp,PTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN)=value;
	if(incol==TEMPCHIin) EOSMAC(eosfulltableextratemp,CHIDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN)=value;
	if(incol==TEMPSin) EOSMAC(eosfulltableextratemp,STOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN)=value;
      }
    }
    else{
      // degen tables the same but with different name
      if(whichtable==FULLTABLE){
	if(utotdegencut[whichtable]<=DEGENCUTLASTOLDVERSION){
	  if(incol==UTOTOFFSETin) EOSMAC(eosfulltabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==PTOTOFFSETin) EOSMAC(eosfulltabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==CHIOFFSETin) EOSMAC(eosfulltabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==STOTOFFSETin) EOSMAC(eosfulltabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	}
	else{ // utotdegencut[whichtable]>DEGENCUTLASTOLDVERSION
	  if(incol==UTOTOFFSETin) EOSMAC(eosfulltabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==PTOTOFFSETin) EOSMAC(eosfulltabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==CHIOFFSETin) EOSMAC(eosfulltabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==STOTOFFSETin) EOSMAC(eosfulltabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;

	  if(incol==UTOTINin) EOSMAC(eosfulltabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN)=value;
	  if(incol==PTOTINin) EOSMAC(eosfulltabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN)=value;
	  if(incol==CHIINin) EOSMAC(eosfulltabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSIN)=value;
	  if(incol==STOTINin) EOSMAC(eosfulltabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN)=value;

	  if(incol==UTOTOUTin) EOSMAC(eosfulltabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT)=value;
	  if(incol==PTOTOUTin) EOSMAC(eosfulltabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT)=value;
	  if(incol==CHIOUTin) EOSMAC(eosfulltabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOUT)=value;
	  if(incol==STOTOUTin) EOSMAC(eosfulltabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT)=value;
	}
      }
      else{
	if(utotdegencut[whichtable]<=DEGENCUTLASTOLDVERSION){
	  if(incol==UTOTOFFSETin) EOSMAC(eosfulltableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==PTOTOFFSETin) EOSMAC(eosfulltableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==CHIOFFSETin) EOSMAC(eosfulltableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==STOTOFFSETin) EOSMAC(eosfulltableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	}
	else{ // utotdegencut[whichtable]>DEGENCUTLASTOLDVERSION
	  if(incol==UTOTOFFSETin) EOSMAC(eosfulltableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==PTOTOFFSETin) EOSMAC(eosfulltableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==CHIOFFSETin) EOSMAC(eosfulltableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==STOTOFFSETin) EOSMAC(eosfulltableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;

	  if(incol==UTOTINin) EOSMAC(eosfulltableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN)=value;
	  if(incol==PTOTINin) EOSMAC(eosfulltableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN)=value;
	  if(incol==CHIINin) EOSMAC(eosfulltableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSIN)=value;
	  if(incol==STOTINin) EOSMAC(eosfulltableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN)=value;

	  if(incol==UTOTOUTin) EOSMAC(eosfulltableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT)=value;
	  if(incol==PTOTOUTin) EOSMAC(eosfulltableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT)=value;
	  if(incol==CHIOUTin) EOSMAC(eosfulltableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOUT)=value;
	  if(incol==STOTOUTin) EOSMAC(eosfulltableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT)=value;
	}
      }
    }
  }
#endif
#if(ALLOWSIMPLETABLE==1)
  else if(whichtable==SIMPLETABLE || whichtable==EXTRASIMPLETABLE){
    if(whichdegen==ISNOTDEGENTABLE){

      if(whichtable==SIMPLETABLE){
	if(incol==PofRHOUin) EOSMAC(eossimpletablestandard,0,mmm,lll,kkk,jjj,iii,PofRHOU)=value;
	if(incol==CS2ofRHOUin) EOSMAC(eossimpletablestandard,0,mmm,lll,kkk,jjj,iii,CS2ofRHOU)=value;

	if(incol==UofRHOPin) EOSMAC(eossimpletableguess,0,mmm,lll,kkk,jjj,iii,UofRHOP)=value;

	if(incol==UofRHOSin) EOSMAC(eossimpletablediss,0,mmm,lll,kkk,jjj,iii,UofRHOS)=value;

	if(incol==DPDRHOofRHOUin) EOSMAC(eossimpletabledp,0,mmm,lll,kkk,jjj,iii,DPDRHOofRHOU)=value;
	if(incol==DPDUofRHOUin) EOSMAC(eossimpletabledp,0,mmm,lll,kkk,jjj,iii,DPDUofRHOU)=value;

	if(incol==SofRHOUin) EOSMAC(eossimpletablesden,0,mmm,lll,kkk,jjj,iii,SofRHOU)=value;
	if(incol==DSDRHOofRHOUin) EOSMAC(eossimpletablesden,0,mmm,lll,kkk,jjj,iii,DSDRHOofRHOU)=value;
	if(incol==DSDUofRHOUin) EOSMAC(eossimpletablesden,0,mmm,lll,kkk,jjj,iii,DSDUofRHOU)=value;

	if(incol==SSofRHOCHIin) EOSMAC(eossimpletablesspec,0,mmm,lll,kkk,jjj,iii,SSofRHOCHI)=value;
	if(incol==DSSDRHOofRHOCHIin) EOSMAC(eossimpletablesspec,0,mmm,lll,kkk,jjj,iii,DSSDRHOofRHOCHI)=value;
	if(incol==DSSDCHIofRHOCHIin) EOSMAC(eossimpletablesspec,0,mmm,lll,kkk,jjj,iii,DSSDCHIofRHOCHI)=value;

	if(incol==PofRHOCHIin) EOSMAC(eossimpletablepofchi,0,mmm,lll,kkk,jjj,iii,PofRHOCHI)=value;
	if(incol==IDRHO0DPin) EOSMAC(eossimpletablepofchi,0,mmm,lll,kkk,jjj,iii,IDRHO0DP)=value;
	if(incol==IDCHIDPin) EOSMAC(eossimpletablepofchi,0,mmm,lll,kkk,jjj,iii,IDCHIDP)=value;

	if(incol==TEMPUin) EOSMAC(eossimpletabletemp,UTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN)=value;
	if(incol==TEMPPin) EOSMAC(eossimpletabletemp,PTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN)=value;
	if(incol==TEMPCHIin) EOSMAC(eossimpletabletemp,CHIDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN)=value;
	if(incol==TEMPSin) EOSMAC(eossimpletabletemp,STOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN)=value;

	// note that if WHICHDATATYPEGENERAL==4, below is just not used
	if(incol>=FIRSTEXTRAin && incol<=LASTEXTRAin) EOSMAC(eossimpletableextra,0,mmm,lll,kkk,jjj,iii,incol-FIRSTEXTRAin+FIRSTEOSEXTRA)=value; // assumes extra's are ordered in sequence
      }
      else{
	// incol is different then since extras start with 4 temperature quantities and then rest of normal extras
	if(incol>=FIRSTEXTRAin && incol<=LASTEXTRAin) EOSMAC(eossimpletableextra,0,mmm,lll,kkk,jjj,iii,incol-FIRSTEXTRAin+FIRSTEOSEXTRA)=value; // assumes extra's are ordered in sequence
	if(incol==TEMPUin) EOSMAC(eossimpletableextratemp,UTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN)=value;
	if(incol==TEMPPin) EOSMAC(eossimpletableextratemp,PTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN)=value;
	if(incol==TEMPCHIin) EOSMAC(eossimpletableextratemp,CHIDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN)=value;
	if(incol==TEMPSin) EOSMAC(eossimpletableextratemp,STOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN)=value;
      }
    }
    else{
      // degen tables the same but with different name
      if(whichtable==SIMPLETABLE){
	if(utotdegencut[whichtable]<=DEGENCUTLASTOLDVERSION){
	  if(incol==UTOTOFFSETin) EOSMAC(eossimpletabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==PTOTOFFSETin) EOSMAC(eossimpletabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==CHIOFFSETin) EOSMAC(eossimpletabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==STOTOFFSETin) EOSMAC(eossimpletabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	}
	else{ // utotdegencut[whichtable]>DEGENCUTLASTOLDVERSION
	  if(incol==UTOTOFFSETin) EOSMAC(eossimpletabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==PTOTOFFSETin) EOSMAC(eossimpletabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==CHIOFFSETin) EOSMAC(eossimpletabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==STOTOFFSETin) EOSMAC(eossimpletabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;

	  if(incol==UTOTINin) EOSMAC(eossimpletabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN)=value;
	  if(incol==PTOTINin) EOSMAC(eossimpletabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN)=value;
	  if(incol==CHIINin) EOSMAC(eossimpletabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSIN)=value;
	  if(incol==STOTINin) EOSMAC(eossimpletabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN)=value;

	  if(incol==UTOTOUTin) EOSMAC(eossimpletabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT)=value;
	  if(incol==PTOTOUTin) EOSMAC(eossimpletabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT)=value;
	  if(incol==CHIOUTin) EOSMAC(eossimpletabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOUT)=value;
	  if(incol==STOTOUTin) EOSMAC(eossimpletabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT)=value;
	}
      }
      else{
	if(utotdegencut[whichtable]<=DEGENCUTLASTOLDVERSION){
	  if(incol==UTOTOFFSETin) EOSMAC(eossimpletableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==PTOTOFFSETin) EOSMAC(eossimpletableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==CHIOFFSETin) EOSMAC(eossimpletableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==STOTOFFSETin) EOSMAC(eossimpletableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	}
	else{ // utotdegencut[whichtable]>DEGENCUTLASTOLDVERSION
	  if(incol==UTOTOFFSETin) EOSMAC(eossimpletableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==PTOTOFFSETin) EOSMAC(eossimpletableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==CHIOFFSETin) EOSMAC(eossimpletableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==STOTOFFSETin) EOSMAC(eossimpletableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;

	  if(incol==UTOTINin) EOSMAC(eossimpletableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN)=value;
	  if(incol==PTOTINin) EOSMAC(eossimpletableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN)=value;
	  if(incol==CHIINin) EOSMAC(eossimpletableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSIN)=value;
	  if(incol==STOTINin) EOSMAC(eossimpletableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN)=value;

	  if(incol==UTOTOUTin) EOSMAC(eossimpletableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT)=value;
	  if(incol==PTOTOUTin) EOSMAC(eossimpletableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT)=value;
	  if(incol==CHIOUTin) EOSMAC(eossimpletableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOUT)=value;
	  if(incol==STOTOUTin) EOSMAC(eossimpletableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT)=value;
	}
      }
    }
  }
#endif
#if(ALLOWSIMPLEZOOMTABLE==1)
  else if(whichtable==SIMPLEZOOMTABLE || whichtable==EXTRASIMPLEZOOMTABLE){
    if(whichdegen==ISNOTDEGENTABLE){

      if(whichtable==SIMPLEZOOMTABLE){
	if(incol==PofRHOUin) EOSMAC(eossimplezoomtablestandard,0,mmm,lll,kkk,jjj,iii,PofRHOU)=value;
	if(incol==CS2ofRHOUin) EOSMAC(eossimplezoomtablestandard,0,mmm,lll,kkk,jjj,iii,CS2ofRHOU)=value;

	if(incol==UofRHOPin) EOSMAC(eossimplezoomtableguess,0,mmm,lll,kkk,jjj,iii,UofRHOP)=value;

	if(incol==UofRHOSin) EOSMAC(eossimplezoomtablediss,0,mmm,lll,kkk,jjj,iii,UofRHOS)=value;

	if(incol==DPDRHOofRHOUin) EOSMAC(eossimplezoomtabledp,0,mmm,lll,kkk,jjj,iii,DPDRHOofRHOU)=value;
	if(incol==DPDUofRHOUin) EOSMAC(eossimplezoomtabledp,0,mmm,lll,kkk,jjj,iii,DPDUofRHOU)=value;

	if(incol==SofRHOUin) EOSMAC(eossimplezoomtablesden,0,mmm,lll,kkk,jjj,iii,SofRHOU)=value;
	if(incol==DSDRHOofRHOUin) EOSMAC(eossimplezoomtablesden,0,mmm,lll,kkk,jjj,iii,DSDRHOofRHOU)=value;
	if(incol==DSDUofRHOUin) EOSMAC(eossimplezoomtablesden,0,mmm,lll,kkk,jjj,iii,DSDUofRHOU)=value;

	if(incol==SSofRHOCHIin) EOSMAC(eossimplezoomtablesspec,0,mmm,lll,kkk,jjj,iii,SSofRHOCHI)=value;
	if(incol==DSSDRHOofRHOCHIin) EOSMAC(eossimplezoomtablesspec,0,mmm,lll,kkk,jjj,iii,DSSDRHOofRHOCHI)=value;
	if(incol==DSSDCHIofRHOCHIin) EOSMAC(eossimplezoomtablesspec,0,mmm,lll,kkk,jjj,iii,DSSDCHIofRHOCHI)=value;

	if(incol==PofRHOCHIin) EOSMAC(eossimplezoomtablepofchi,0,mmm,lll,kkk,jjj,iii,PofRHOCHI)=value;
	if(incol==IDRHO0DPin) EOSMAC(eossimplezoomtablepofchi,0,mmm,lll,kkk,jjj,iii,IDRHO0DP)=value;
	if(incol==IDCHIDPin) EOSMAC(eossimplezoomtablepofchi,0,mmm,lll,kkk,jjj,iii,IDCHIDP)=value;

	if(incol==TEMPUin) EOSMAC(eossimplezoomtabletemp,UTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN)=value;
	if(incol==TEMPPin) EOSMAC(eossimplezoomtabletemp,PTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN)=value;
	if(incol==TEMPCHIin) EOSMAC(eossimplezoomtabletemp,CHIDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN)=value;
	if(incol==TEMPSin) EOSMAC(eossimplezoomtabletemp,STOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN)=value;

	// note that if WHICHDATATYPEGENERAL==4, below is just not used
	if(incol>=FIRSTEXTRAin && incol<=LASTEXTRAin) EOSMAC(eossimplezoomtableextra,0,mmm,lll,kkk,jjj,iii,incol-FIRSTEXTRAin+FIRSTEOSEXTRA)=value; // assumes extra's are ordered in sequence
      }
      else{
	// incol is different then since extras start with 4 temperature quantities and then rest of normal extras
	if(incol>=FIRSTEXTRAin && incol<=LASTEXTRAin) EOSMAC(eossimplezoomtableextra,0,mmm,lll,kkk,jjj,iii,incol-FIRSTEXTRAin+FIRSTEOSEXTRA)=value; // assumes extra's are ordered in sequence
	if(incol==TEMPUin) EOSMAC(eossimplezoomtableextratemp,UTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN)=value;
	if(incol==TEMPPin) EOSMAC(eossimplezoomtableextratemp,PTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN)=value;
	if(incol==TEMPCHIin) EOSMAC(eossimplezoomtableextratemp,CHIDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN)=value;
	if(incol==TEMPSin) EOSMAC(eossimplezoomtableextratemp,STOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN)=value;
      }
    }
    else{
      // degen tables the same but with different name
      if(whichtable==SIMPLEZOOMTABLE){
	if(utotdegencut[whichtable]<=DEGENCUTLASTOLDVERSION){
	  if(incol==UTOTOFFSETin) EOSMAC(eossimplezoomtabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==PTOTOFFSETin) EOSMAC(eossimplezoomtabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==CHIOFFSETin) EOSMAC(eossimplezoomtabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==STOTOFFSETin) EOSMAC(eossimplezoomtabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	}
	else{ // utotdegencut[whichtable]>DEGENCUTLASTOLDVERSION
	  if(incol==UTOTOFFSETin) EOSMAC(eossimplezoomtabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==PTOTOFFSETin) EOSMAC(eossimplezoomtabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==CHIOFFSETin) EOSMAC(eossimplezoomtabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==STOTOFFSETin) EOSMAC(eossimplezoomtabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;

	  if(incol==UTOTINin) EOSMAC(eossimplezoomtabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN)=value;
	  if(incol==PTOTINin) EOSMAC(eossimplezoomtabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN)=value;
	  if(incol==CHIINin) EOSMAC(eossimplezoomtabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSIN)=value;
	  if(incol==STOTINin) EOSMAC(eossimplezoomtabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN)=value;

	  if(incol==UTOTOUTin) EOSMAC(eossimplezoomtabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT)=value;
	  if(incol==PTOTOUTin) EOSMAC(eossimplezoomtabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT)=value;
	  if(incol==CHIOUTin) EOSMAC(eossimplezoomtabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOUT)=value;
	  if(incol==STOTOUTin) EOSMAC(eossimplezoomtabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT)=value;
	}
      }
      else{
	if(utotdegencut[whichtable]<=DEGENCUTLASTOLDVERSION){
	  if(incol==UTOTOFFSETin) EOSMAC(eossimplezoomtableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==PTOTOFFSETin) EOSMAC(eossimplezoomtableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==CHIOFFSETin) EOSMAC(eossimplezoomtableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==STOTOFFSETin) EOSMAC(eossimplezoomtableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	}
	else{ // utotdegencut[whichtable]>DEGENCUTLASTOLDVERSION
	  if(incol==UTOTOFFSETin) EOSMAC(eossimplezoomtableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==PTOTOFFSETin) EOSMAC(eossimplezoomtableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==CHIOFFSETin) EOSMAC(eossimplezoomtableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;
	  if(incol==STOTOFFSETin) EOSMAC(eossimplezoomtableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET)=value;

	  if(incol==UTOTINin) EOSMAC(eossimplezoomtableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN)=value;
	  if(incol==PTOTINin) EOSMAC(eossimplezoomtableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN)=value;
	  if(incol==CHIINin) EOSMAC(eossimplezoomtableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSIN)=value;
	  if(incol==STOTINin) EOSMAC(eossimplezoomtableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN)=value;

	  if(incol==UTOTOUTin) EOSMAC(eossimplezoomtableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT)=value;
	  if(incol==PTOTOUTin) EOSMAC(eossimplezoomtableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT)=value;
	  if(incol==CHIOUTin) EOSMAC(eossimplezoomtableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOUT)=value;
	  if(incol==STOTOUTin) EOSMAC(eossimplezoomtableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT)=value;
	}
      }
    }
  }
#endif

}


















// translate HARM EOS table columns into input type columns (i.e. eosarray -> value)
// note that "incol" is *input* column number, not HARM EOS numbers, so this function should be just like seteostable() except LHS and RHS are flipped,
// so regexp it:
// EOSMAC(\(.*\))=value; -> *value=EOSMAC(\1);
// and last function parameter goes from "double value" -> "double *value"
//
// input jjj should be 0 if inputting degen table data
// only used by read_setup_eostable()
// no longer refer to macro label "inextra" since put into canonical order and position before calling set_ or get_ functions
static void get_arrays_eostable(int whichdegen, int whichtable, int mmm, int lll, int kkk, int jjj, int iii, int incol, double *value)
{

  
  // overrides:
  // then modify indices since different sub tables have different dimensions.  This will cause repeated-reads to put into same locations in memory.
  if(WHICHDATATYPEGENERAL==4){
    if(
       ((whichtable==FULLTABLE || whichtable==SIMPLETABLE || whichtable==SIMPLEZOOMTABLE) && (incol<FIRSTEXTRAin && incol>LASTEXTRAin)) || // for all non-extras in normal table
       ((whichtable==EXTRAFULLTABLE || whichtable==EXTRASIMPLETABLE || whichtable==EXTRASIMPLEZOOMTABLE) && (incol<FIRSTEXTRAin && incol>LASTEXTRAin)) || // for temperature in extra table
       (whichdegen==1)
       ){
      // then table not storing mmm or lll, so set to zero
      mmm=lll=0;
    }
  }
  if(whichdegen){
    jjj=0; // in case not already set
  }


  if(0){
  }
#if(ALLOWFULLTABLE==1)
  else if(whichtable==FULLTABLE || whichtable==EXTRAFULLTABLE){
    if(whichdegen==ISNOTDEGENTABLE){

      if(whichtable==FULLTABLE){
	if(incol==PofRHOUin) *value=EOSMAC(eosfulltablestandard,0,mmm,lll,kkk,jjj,iii,PofRHOU);
	if(incol==CS2ofRHOUin) *value=EOSMAC(eosfulltablestandard,0,mmm,lll,kkk,jjj,iii,CS2ofRHOU);

	if(incol==UofRHOPin) *value=EOSMAC(eosfulltableguess,0,mmm,lll,kkk,jjj,iii,UofRHOP);

	if(incol==UofRHOSin) *value=EOSMAC(eosfulltablediss,0,mmm,lll,kkk,jjj,iii,UofRHOS);

	if(incol==DPDRHOofRHOUin) *value=EOSMAC(eosfulltabledp,0,mmm,lll,kkk,jjj,iii,DPDRHOofRHOU);
	if(incol==DPDUofRHOUin) *value=EOSMAC(eosfulltabledp,0,mmm,lll,kkk,jjj,iii,DPDUofRHOU);

	if(incol==SofRHOUin) *value=EOSMAC(eosfulltablesden,0,mmm,lll,kkk,jjj,iii,SofRHOU);
	if(incol==DSDRHOofRHOUin) *value=EOSMAC(eosfulltablesden,0,mmm,lll,kkk,jjj,iii,DSDRHOofRHOU);
	if(incol==DSDUofRHOUin) *value=EOSMAC(eosfulltablesden,0,mmm,lll,kkk,jjj,iii,DSDUofRHOU);

	if(incol==SSofRHOCHIin) *value=EOSMAC(eosfulltablesspec,0,mmm,lll,kkk,jjj,iii,SSofRHOCHI);
	if(incol==DSSDRHOofRHOCHIin) *value=EOSMAC(eosfulltablesspec,0,mmm,lll,kkk,jjj,iii,DSSDRHOofRHOCHI);
	if(incol==DSSDCHIofRHOCHIin) *value=EOSMAC(eosfulltablesspec,0,mmm,lll,kkk,jjj,iii,DSSDCHIofRHOCHI);

	if(incol==PofRHOCHIin) *value=EOSMAC(eosfulltablepofchi,0,mmm,lll,kkk,jjj,iii,PofRHOCHI);
	if(incol==IDRHO0DPin) *value=EOSMAC(eosfulltablepofchi,0,mmm,lll,kkk,jjj,iii,IDRHO0DP);
	if(incol==IDCHIDPin) *value=EOSMAC(eosfulltablepofchi,0,mmm,lll,kkk,jjj,iii,IDCHIDP);

	if(incol==TEMPUin) *value=EOSMAC(eosfulltabletemp,UTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN);
	if(incol==TEMPPin) *value=EOSMAC(eosfulltabletemp,PTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN);
	if(incol==TEMPCHIin) *value=EOSMAC(eosfulltabletemp,CHIDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN);
	if(incol==TEMPSin) *value=EOSMAC(eosfulltabletemp,STOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN);

	// note that if WHICHDATATYPEGENERAL==4, below is just not used
	if(incol>=FIRSTEXTRAin && incol<=LASTEXTRAin) *value=EOSMAC(eosfulltableextra,0,mmm,lll,kkk,jjj,iii,incol-FIRSTEXTRAin+FIRSTEOSEXTRA); // assumes extra's are ordered in sequence
      }
      else{
	// incol is different then since extras start with 4 temperature quantities and then rest of normal extras
	if(incol>=FIRSTEXTRAin && incol<=LASTEXTRAin) *value=EOSMAC(eosfulltableextra,0,mmm,lll,kkk,jjj,iii,incol-FIRSTEXTRAin+FIRSTEOSEXTRA); // assumes extra's are ordered in sequence
	if(incol==TEMPUin) *value=EOSMAC(eosfulltableextratemp,UTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN);
	if(incol==TEMPPin) *value=EOSMAC(eosfulltableextratemp,PTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN);
	if(incol==TEMPCHIin) *value=EOSMAC(eosfulltableextratemp,CHIDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN);
	if(incol==TEMPSin) *value=EOSMAC(eosfulltableextratemp,STOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN);
      }
    }
    else{
      // degen tables the same but with different name
      if(whichtable==FULLTABLE){
	if(utotdegencut[whichtable]<=DEGENCUTLASTOLDVERSION){
	  if(incol==UTOTOFFSETin) *value=EOSMAC(eosfulltabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==PTOTOFFSETin) *value=EOSMAC(eosfulltabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==CHIOFFSETin) *value=EOSMAC(eosfulltabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==STOTOFFSETin) *value=EOSMAC(eosfulltabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	}
	else{ // utotdegencut[whichtable]>DEGENCUTLASTOLDVERSION
	  if(incol==UTOTOFFSETin) *value=EOSMAC(eosfulltabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==PTOTOFFSETin) *value=EOSMAC(eosfulltabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==CHIOFFSETin) *value=EOSMAC(eosfulltabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==STOTOFFSETin) *value=EOSMAC(eosfulltabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);

	  if(incol==UTOTINin) *value=EOSMAC(eosfulltabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN);
	  if(incol==PTOTINin) *value=EOSMAC(eosfulltabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN);
	  if(incol==CHIINin) *value=EOSMAC(eosfulltabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSIN);
	  if(incol==STOTINin) *value=EOSMAC(eosfulltabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN);

	  if(incol==UTOTOUTin) *value=EOSMAC(eosfulltabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT);
	  if(incol==PTOTOUTin) *value=EOSMAC(eosfulltabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT);
	  if(incol==CHIOUTin) *value=EOSMAC(eosfulltabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOUT);
	  if(incol==STOTOUTin) *value=EOSMAC(eosfulltabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT);
	}
      }
      else{
	if(utotdegencut[whichtable]<=DEGENCUTLASTOLDVERSION){
	  if(incol==UTOTOFFSETin) *value=EOSMAC(eosfulltableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==PTOTOFFSETin) *value=EOSMAC(eosfulltableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==CHIOFFSETin) *value=EOSMAC(eosfulltableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==STOTOFFSETin) *value=EOSMAC(eosfulltableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	}
	else{ // utotdegencut[whichtable]>DEGENCUTLASTOLDVERSION
	  if(incol==UTOTOFFSETin) *value=EOSMAC(eosfulltableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==PTOTOFFSETin) *value=EOSMAC(eosfulltableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==CHIOFFSETin) *value=EOSMAC(eosfulltableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==STOTOFFSETin) *value=EOSMAC(eosfulltableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);

	  if(incol==UTOTINin) *value=EOSMAC(eosfulltableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN);
	  if(incol==PTOTINin) *value=EOSMAC(eosfulltableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN);
	  if(incol==CHIINin) *value=EOSMAC(eosfulltableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSIN);
	  if(incol==STOTINin) *value=EOSMAC(eosfulltableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN);

	  if(incol==UTOTOUTin) *value=EOSMAC(eosfulltableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT);
	  if(incol==PTOTOUTin) *value=EOSMAC(eosfulltableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT);
	  if(incol==CHIOUTin) *value=EOSMAC(eosfulltableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOUT);
	  if(incol==STOTOUTin) *value=EOSMAC(eosfulltableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT);
	}
      }
    }
  }
#endif
#if(ALLOWSIMPLETABLE==1)
  else if(whichtable==SIMPLETABLE || whichtable==EXTRASIMPLETABLE){
    if(whichdegen==ISNOTDEGENTABLE){

      if(whichtable==SIMPLETABLE){
	if(incol==PofRHOUin) *value=EOSMAC(eossimpletablestandard,0,mmm,lll,kkk,jjj,iii,PofRHOU);
	if(incol==CS2ofRHOUin) *value=EOSMAC(eossimpletablestandard,0,mmm,lll,kkk,jjj,iii,CS2ofRHOU);

	if(incol==UofRHOPin) *value=EOSMAC(eossimpletableguess,0,mmm,lll,kkk,jjj,iii,UofRHOP);

	if(incol==UofRHOSin) *value=EOSMAC(eossimpletablediss,0,mmm,lll,kkk,jjj,iii,UofRHOS);

	if(incol==DPDRHOofRHOUin) *value=EOSMAC(eossimpletabledp,0,mmm,lll,kkk,jjj,iii,DPDRHOofRHOU);
	if(incol==DPDUofRHOUin) *value=EOSMAC(eossimpletabledp,0,mmm,lll,kkk,jjj,iii,DPDUofRHOU);

	if(incol==SofRHOUin) *value=EOSMAC(eossimpletablesden,0,mmm,lll,kkk,jjj,iii,SofRHOU);
	if(incol==DSDRHOofRHOUin) *value=EOSMAC(eossimpletablesden,0,mmm,lll,kkk,jjj,iii,DSDRHOofRHOU);
	if(incol==DSDUofRHOUin) *value=EOSMAC(eossimpletablesden,0,mmm,lll,kkk,jjj,iii,DSDUofRHOU);

	if(incol==SSofRHOCHIin) *value=EOSMAC(eossimpletablesspec,0,mmm,lll,kkk,jjj,iii,SSofRHOCHI);
	if(incol==DSSDRHOofRHOCHIin) *value=EOSMAC(eossimpletablesspec,0,mmm,lll,kkk,jjj,iii,DSSDRHOofRHOCHI);
	if(incol==DSSDCHIofRHOCHIin) *value=EOSMAC(eossimpletablesspec,0,mmm,lll,kkk,jjj,iii,DSSDCHIofRHOCHI);

	if(incol==PofRHOCHIin) *value=EOSMAC(eossimpletablepofchi,0,mmm,lll,kkk,jjj,iii,PofRHOCHI);
	if(incol==IDRHO0DPin) *value=EOSMAC(eossimpletablepofchi,0,mmm,lll,kkk,jjj,iii,IDRHO0DP);
	if(incol==IDCHIDPin) *value=EOSMAC(eossimpletablepofchi,0,mmm,lll,kkk,jjj,iii,IDCHIDP);

	if(incol==TEMPUin) *value=EOSMAC(eossimpletabletemp,UTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN);
	if(incol==TEMPPin) *value=EOSMAC(eossimpletabletemp,PTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN);
	if(incol==TEMPCHIin) *value=EOSMAC(eossimpletabletemp,CHIDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN);
	if(incol==TEMPSin) *value=EOSMAC(eossimpletabletemp,STOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN);

	// note that if WHICHDATATYPEGENERAL==4, below is just not used
	if(incol>=FIRSTEXTRAin && incol<=LASTEXTRAin) *value=EOSMAC(eossimpletableextra,0,mmm,lll,kkk,jjj,iii,incol-FIRSTEXTRAin+FIRSTEOSEXTRA); // assumes extra's are ordered in sequence
      }
      else{
	// incol is different then since extras start with 4 temperature quantities and then rest of normal extras
	if(incol>=FIRSTEXTRAin && incol<=LASTEXTRAin) *value=EOSMAC(eossimpletableextra,0,mmm,lll,kkk,jjj,iii,incol-FIRSTEXTRAin+FIRSTEOSEXTRA); // assumes extra's are ordered in sequence
	if(incol==TEMPUin) *value=EOSMAC(eossimpletableextratemp,UTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN);
	if(incol==TEMPPin) *value=EOSMAC(eossimpletableextratemp,PTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN);
	if(incol==TEMPCHIin) *value=EOSMAC(eossimpletableextratemp,CHIDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN);
	if(incol==TEMPSin) *value=EOSMAC(eossimpletableextratemp,STOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN);
      }
    }
    else{
      // degen tables the same but with different name
      if(whichtable==SIMPLETABLE){
	if(utotdegencut[whichtable]<=DEGENCUTLASTOLDVERSION){
	  if(incol==UTOTOFFSETin) *value=EOSMAC(eossimpletabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==PTOTOFFSETin) *value=EOSMAC(eossimpletabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==CHIOFFSETin) *value=EOSMAC(eossimpletabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==STOTOFFSETin) *value=EOSMAC(eossimpletabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	}
	else{ // utotdegencut[whichtable]>DEGENCUTLASTOLDVERSION
	  if(incol==UTOTOFFSETin) *value=EOSMAC(eossimpletabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==PTOTOFFSETin) *value=EOSMAC(eossimpletabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==CHIOFFSETin) *value=EOSMAC(eossimpletabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==STOTOFFSETin) *value=EOSMAC(eossimpletabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);

	  if(incol==UTOTINin) *value=EOSMAC(eossimpletabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN);
	  if(incol==PTOTINin) *value=EOSMAC(eossimpletabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN);
	  if(incol==CHIINin) *value=EOSMAC(eossimpletabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSIN);
	  if(incol==STOTINin) *value=EOSMAC(eossimpletabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN);

	  if(incol==UTOTOUTin) *value=EOSMAC(eossimpletabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT);
	  if(incol==PTOTOUTin) *value=EOSMAC(eossimpletabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT);
	  if(incol==CHIOUTin) *value=EOSMAC(eossimpletabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOUT);
	  if(incol==STOTOUTin) *value=EOSMAC(eossimpletabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT);
	}
      }
      else{
	if(utotdegencut[whichtable]<=DEGENCUTLASTOLDVERSION){
	  if(incol==UTOTOFFSETin) *value=EOSMAC(eossimpletableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==PTOTOFFSETin) *value=EOSMAC(eossimpletableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==CHIOFFSETin) *value=EOSMAC(eossimpletableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==STOTOFFSETin) *value=EOSMAC(eossimpletableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	}
	else{ // utotdegencut[whichtable]>DEGENCUTLASTOLDVERSION
	  if(incol==UTOTOFFSETin) *value=EOSMAC(eossimpletableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==PTOTOFFSETin) *value=EOSMAC(eossimpletableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==CHIOFFSETin) *value=EOSMAC(eossimpletableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==STOTOFFSETin) *value=EOSMAC(eossimpletableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);

	  if(incol==UTOTINin) *value=EOSMAC(eossimpletableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN);
	  if(incol==PTOTINin) *value=EOSMAC(eossimpletableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN);
	  if(incol==CHIINin) *value=EOSMAC(eossimpletableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSIN);
	  if(incol==STOTINin) *value=EOSMAC(eossimpletableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN);

	  if(incol==UTOTOUTin) *value=EOSMAC(eossimpletableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT);
	  if(incol==PTOTOUTin) *value=EOSMAC(eossimpletableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT);
	  if(incol==CHIOUTin) *value=EOSMAC(eossimpletableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOUT);
	  if(incol==STOTOUTin) *value=EOSMAC(eossimpletableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT);
	}
      }
    }
  }
#endif
#if(ALLOWSIMPLEZOOMTABLE==1)
  else if(whichtable==SIMPLEZOOMTABLE || whichtable==EXTRASIMPLEZOOMTABLE){
    if(whichdegen==ISNOTDEGENTABLE){

      if(whichtable==SIMPLEZOOMTABLE){
	if(incol==PofRHOUin) *value=EOSMAC(eossimplezoomtablestandard,0,mmm,lll,kkk,jjj,iii,PofRHOU);
	if(incol==CS2ofRHOUin) *value=EOSMAC(eossimplezoomtablestandard,0,mmm,lll,kkk,jjj,iii,CS2ofRHOU);

	if(incol==UofRHOPin) *value=EOSMAC(eossimplezoomtableguess,0,mmm,lll,kkk,jjj,iii,UofRHOP);

	if(incol==UofRHOSin) *value=EOSMAC(eossimplezoomtablediss,0,mmm,lll,kkk,jjj,iii,UofRHOS);

	if(incol==DPDRHOofRHOUin) *value=EOSMAC(eossimplezoomtabledp,0,mmm,lll,kkk,jjj,iii,DPDRHOofRHOU);
	if(incol==DPDUofRHOUin) *value=EOSMAC(eossimplezoomtabledp,0,mmm,lll,kkk,jjj,iii,DPDUofRHOU);

	if(incol==SofRHOUin) *value=EOSMAC(eossimplezoomtablesden,0,mmm,lll,kkk,jjj,iii,SofRHOU);
	if(incol==DSDRHOofRHOUin) *value=EOSMAC(eossimplezoomtablesden,0,mmm,lll,kkk,jjj,iii,DSDRHOofRHOU);
	if(incol==DSDUofRHOUin) *value=EOSMAC(eossimplezoomtablesden,0,mmm,lll,kkk,jjj,iii,DSDUofRHOU);

	if(incol==SSofRHOCHIin) *value=EOSMAC(eossimplezoomtablesspec,0,mmm,lll,kkk,jjj,iii,SSofRHOCHI);
	if(incol==DSSDRHOofRHOCHIin) *value=EOSMAC(eossimplezoomtablesspec,0,mmm,lll,kkk,jjj,iii,DSSDRHOofRHOCHI);
	if(incol==DSSDCHIofRHOCHIin) *value=EOSMAC(eossimplezoomtablesspec,0,mmm,lll,kkk,jjj,iii,DSSDCHIofRHOCHI);

	if(incol==PofRHOCHIin) *value=EOSMAC(eossimplezoomtablepofchi,0,mmm,lll,kkk,jjj,iii,PofRHOCHI);
	if(incol==IDRHO0DPin) *value=EOSMAC(eossimplezoomtablepofchi,0,mmm,lll,kkk,jjj,iii,IDRHO0DP);
	if(incol==IDCHIDPin) *value=EOSMAC(eossimplezoomtablepofchi,0,mmm,lll,kkk,jjj,iii,IDCHIDP);

	if(incol==TEMPUin) *value=EOSMAC(eossimplezoomtabletemp,UTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN);
	if(incol==TEMPPin) *value=EOSMAC(eossimplezoomtabletemp,PTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN);
	if(incol==TEMPCHIin) *value=EOSMAC(eossimplezoomtabletemp,CHIDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN);
	if(incol==TEMPSin) *value=EOSMAC(eossimplezoomtabletemp,STOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN);

	// note that if WHICHDATATYPEGENERAL==4, below is just not used
	if(incol>=FIRSTEXTRAin && incol<=LASTEXTRAin) *value=EOSMAC(eossimplezoomtableextra,0,mmm,lll,kkk,jjj,iii,incol-FIRSTEXTRAin+FIRSTEOSEXTRA); // assumes extra's are ordered in sequence
      }
      else{
	// incol is different then since extras start with 4 temperature quantities and then rest of normal extras
	if(incol>=FIRSTEXTRAin && incol<=LASTEXTRAin) *value=EOSMAC(eossimplezoomtableextra,0,mmm,lll,kkk,jjj,iii,incol-FIRSTEXTRAin+FIRSTEOSEXTRA); // assumes extra's are ordered in sequence
	if(incol==TEMPUin) *value=EOSMAC(eossimplezoomtableextratemp,UTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN);
	if(incol==TEMPPin) *value=EOSMAC(eossimplezoomtableextratemp,PTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN);
	if(incol==TEMPCHIin) *value=EOSMAC(eossimplezoomtableextratemp,CHIDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN);
	if(incol==TEMPSin) *value=EOSMAC(eossimplezoomtableextratemp,STOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN);
      }
    }
    else{
      // degen tables the same but with different name
      if(whichtable==SIMPLEZOOMTABLE){
	if(utotdegencut[whichtable]<=DEGENCUTLASTOLDVERSION){
	  if(incol==UTOTOFFSETin) *value=EOSMAC(eossimplezoomtabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==PTOTOFFSETin) *value=EOSMAC(eossimplezoomtabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==CHIOFFSETin) *value=EOSMAC(eossimplezoomtabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==STOTOFFSETin) *value=EOSMAC(eossimplezoomtabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	}
	else{ // utotdegencut[whichtable]>DEGENCUTLASTOLDVERSION
	  if(incol==UTOTOFFSETin) *value=EOSMAC(eossimplezoomtabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==PTOTOFFSETin) *value=EOSMAC(eossimplezoomtabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==CHIOFFSETin) *value=EOSMAC(eossimplezoomtabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==STOTOFFSETin) *value=EOSMAC(eossimplezoomtabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);

	  if(incol==UTOTINin) *value=EOSMAC(eossimplezoomtabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN);
	  if(incol==PTOTINin) *value=EOSMAC(eossimplezoomtabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN);
	  if(incol==CHIINin) *value=EOSMAC(eossimplezoomtabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSIN);
	  if(incol==STOTINin) *value=EOSMAC(eossimplezoomtabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN);

	  if(incol==UTOTOUTin) *value=EOSMAC(eossimplezoomtabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT);
	  if(incol==PTOTOUTin) *value=EOSMAC(eossimplezoomtabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT);
	  if(incol==CHIOUTin) *value=EOSMAC(eossimplezoomtabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOUT);
	  if(incol==STOTOUTin) *value=EOSMAC(eossimplezoomtabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT);
	}
      }
      else{
	if(utotdegencut[whichtable]<=DEGENCUTLASTOLDVERSION){
	  if(incol==UTOTOFFSETin) *value=EOSMAC(eossimplezoomtableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==PTOTOFFSETin) *value=EOSMAC(eossimplezoomtableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==CHIOFFSETin) *value=EOSMAC(eossimplezoomtableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==STOTOFFSETin) *value=EOSMAC(eossimplezoomtableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	}
	else{ // utotdegencut[whichtable]>DEGENCUTLASTOLDVERSION
	  if(incol==UTOTOFFSETin) *value=EOSMAC(eossimplezoomtableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==PTOTOFFSETin) *value=EOSMAC(eossimplezoomtableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==CHIOFFSETin) *value=EOSMAC(eossimplezoomtableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);
	  if(incol==STOTOFFSETin) *value=EOSMAC(eossimplezoomtableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET);

	  if(incol==UTOTINin) *value=EOSMAC(eossimplezoomtableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN);
	  if(incol==PTOTINin) *value=EOSMAC(eossimplezoomtableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN);
	  if(incol==CHIINin) *value=EOSMAC(eossimplezoomtableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSIN);
	  if(incol==STOTINin) *value=EOSMAC(eossimplezoomtableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN);

	  if(incol==UTOTOUTin) *value=EOSMAC(eossimplezoomtableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT);
	  if(incol==PTOTOUTin) *value=EOSMAC(eossimplezoomtableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT);
	  if(incol==CHIOUTin) *value=EOSMAC(eossimplezoomtableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOUT);
	  if(incol==STOTOUTin) *value=EOSMAC(eossimplezoomtableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT);
	}
      }
    }
  }
#endif

}






// determine whether should do log interpolation based upon whichtablesubtype and coli
static int get_dologinterp_subtype(int whichtablesubtype, int coli)
{
  


#if(DOLOGINTERP)

  // translate from coli to whichfun numbering system
  int whichfun = coli+firsteosintablesubtype[whichtablesubtype];

  if(whichtablesubtype==SUBTYPEDEGEN){
    // all degens are log interpolated.
    if(whichfun>=EOSOFFSET && whichfun<=EOSOUT) return(1);
  }
  else if(whichtablesubtype==SUBTYPESTANDARD){
    if(whichfun==PofRHOU) return(1);
    else if(whichfun==CS2ofRHOU) return(1);
  }
  else if(whichtablesubtype==SUBTYPEGUESS){
    if(whichfun==UofRHOP) return(1);
  }
  else if(whichtablesubtype==SUBTYPEDISS){
    if(whichfun==UofRHOS) return(1);
  }
  else if(whichtablesubtype==SUBTYPEDP){
    if(whichfun==DPDRHOofRHOU) return(0);
    else if(whichfun==DPDUofRHOU) return(0);
  }
  else if(whichtablesubtype==SUBTYPESDEN){
    if(whichfun==SofRHOU) return(1);
    else if(whichfun==DSDRHOofRHOU) return(0);
    else if(whichfun==DSDUofRHOU) return(0);
  }
  else if(whichtablesubtype==SUBTYPESSPEC){
    if(whichfun==SSofRHOCHI) return(1);
    else if(whichfun==DSSDRHOofRHOCHI) return(0);
    else if(whichfun==DSSDCHIofRHOCHI) return(0);
  }
  else if(whichtablesubtype==SUBTYPEPOFCHI){
    if(whichfun==PofRHOCHI) return(1);
    else if(whichfun==IDRHO0DP) return(0);
    else if(whichfun==IDCHIDP) return(0);
  }
  else if(whichtablesubtype==SUBTYPETEMP){
    // all log interpolated
    if(whichfun==TEMPGEN) return(1);
  }
  else if(whichtablesubtype==SUBTYPEEXTRA){
    // all log interpolated
    if(whichfun>=FIRSTEOSEXTRA && whichfun<=LASTEOSEXTRA) return(1);
  }


  // if here, then no matching case found
  dualfprintf(fail_file,"Undefined whichtablesubtype=%d whichfun=%d in get_dologinterp_subtype()\n",whichtablesubtype, whichfun);
  myexit(62662);
  return(-1);

#else
  // not doing log interpolation
  return(0);
#endif


}



// determine whether should do log interpolation based upon whichtablesubtype and coli
static int get_dologinterp_subtype_wrapper(int degentable, int whichtablesubtype, int numcols, int *shouldloginterp)
{
  int logwhichtablesubtype,logcoli;
  int coli;

  logwhichtablesubtype = (degentable==ISNOTDEGENTABLE ? whichtablesubtype : SUBTYPEDEGEN);
  for(coli=0;coli<numcols;coli++){
    logcoli = (degentable==ISNOTDEGENTABLE ? coli : FIRSTEOSDEGEN+coli);
    shouldloginterp[coli]=get_dologinterp_subtype(logwhichtablesubtype, logcoli);
  }

  return(0);

}
