/////////////////////////
//
// initialization and read-in of tables
// accesses global EOS arrays but only once for entire simulation
//
/////////////////////////

static void reorder_extratype(int isextratype, int numeosquantitiesinfilevar, int numeosquantitiesinstandardlist, FTYPEEOS *values);
static int get_numcols(int whichtable, int whichtablesubtype);
static int isextratable(int i);
static void assigntablecolumnnumbers(int i);
static void setvartypeinfo(void);
static void settablesubtypeinfo(void);
static void settablesizeexpected(int tablesizeexpected[][NUMEOSINDEPS]);


static void set_arrays_eostable(int whichdegen, int whichtable, int mmm, int lll, int kkk, int jjj, int iii, int incol, FTYPEEOS value);
static void get_arrays_eostable(int whichdegen, int whichtable, int mmm, int lll, int kkk, int jjj, int iii, int incol, FTYPEEOS *value);

static int get_dologinterp_in(int degentable, int whichfun);
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
static void reorder_extratype(int isextratype, int numeosquantitiesinfilevar, int numeosquantitiesinstandardlistvar, FTYPEEOS *values)
{
  int i,j;


  if(isextratype){
    // then translate "inextra" positions to "in" positions so rest of code can be the same
    FTYPEEOS oldvalues[MAXEOSPIPELINE];

    for(i=0;i<numeosquantitiesinfilevar;i++){
      oldvalues[i]=values[i];
    }

    // initialize new values to zero just to see in output if checking what things are not set
    for(j=0;j<numeosquantitiesinstandardlistvar;j++){
      values[j]=0;
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
  else if(i==FULLTABLEEXTRA || i==SIMPLETABLEEXTRA || i==EXTRASIMPLEZOOMTABLE) extratype=1;

  return(extratype);
}


static void assigntablecolumnnumbers(int i)
{
  int extratype;

  // get whether extra type or not
  extratype=isextratable(i);


  if( ((i==FULLTABLE || i==FULLTABLEEXTRA) && ALLOWFULLTABLE) || ((i==SIMPLETABLE || i==SIMPLETABLEEXTRA) && ALLOWSIMPLETABLE) || ((i==SIMPLEZOOMTABLE || i==EXTRASIMPLEZOOMTABLE) && ALLOWSIMPLEZOOMTABLE) ){

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

  if(ALLOWFULLTABLE){
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

    // FULLTABLEEXTRA
    i=0;
    tablesizeexpected[FULLTABLEEXTRA][i]=EOSFULLEXTRAN1; i++;
    for(templikeiter=FIRSTTKLIKE;templikeiter<=LASTTKLIKE;templikeiter++){
      tablesizeexpected[FULLTABLEEXTRA][i]=EOSFULLEXTRAN2; i++;
    }
    tablesizeexpected[FULLTABLEEXTRA][i]=EOSFULLEXTRAN3; i++;
    tablesizeexpected[FULLTABLEEXTRA][i]=EOSFULLEXTRAN4; i++;
    tablesizeexpected[FULLTABLEEXTRA][i]=EOSFULLEXTRAN5; i++;
    if(i!=NUMEOSINDEPS){
      dualfprintf(fail_file,"tablesizeexpected(full) not setup for that many indepdimens: i=%d\n",i);
      myexit(3206882);
    }
  }



  if(ALLOWSIMPLETABLE){
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
    // SIMPLETABLEEXTRA
    i=0;
    tablesizeexpected[SIMPLETABLEEXTRA][i]=EOSSIMPLEEXTRAN1; i++;
    for(templikeiter=FIRSTTKLIKE;templikeiter<=LASTTKLIKE;templikeiter++){
      tablesizeexpected[SIMPLETABLEEXTRA][i]=EOSSIMPLEEXTRAN2; i++;
    }
    tablesizeexpected[SIMPLETABLEEXTRA][i]=EOSSIMPLEEXTRAN3; i++;
    tablesizeexpected[SIMPLETABLEEXTRA][i]=EOSSIMPLEEXTRAN4; i++;
    tablesizeexpected[SIMPLETABLEEXTRA][i]=EOSSIMPLEEXTRAN5; i++;
    if(i!=NUMEOSINDEPS){
      dualfprintf(fail_file,"tablesizeexpected(simple) not setup for that many indepdimens: i=%d\n",i);
      myexit(3206883);
    }
  }


  if(ALLOWSIMPLEZOOMTABLE){

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
    tablesizeexpected[EXTRASIMPLEZOOMTABLE][i]=EOSSIMPLEZOOMEXTRAN1; i++;
    for(templikeiter=FIRSTTKLIKE;templikeiter<=LASTTKLIKE;templikeiter++){
      tablesizeexpected[EXTRASIMPLEZOOMTABLE][i]=EOSSIMPLEZOOMEXTRAN2; i++;
    }
    tablesizeexpected[EXTRASIMPLEZOOMTABLE][i]=EOSSIMPLEZOOMEXTRAN3; i++;
    tablesizeexpected[EXTRASIMPLEZOOMTABLE][i]=EOSSIMPLEZOOMEXTRAN4; i++;
    tablesizeexpected[EXTRASIMPLEZOOMTABLE][i]=EOSSIMPLEZOOMEXTRAN5; i++;
    if(i!=NUMEOSINDEPS){
      dualfprintf(fail_file,"tablesizeexpected(simplezoom) not setup for that many indepdimens: i=%d\n",i);
      myexit(3206884);
    }

  }



}





// get true number of columns in subtype (exception is degen case that depends upon read-in utotdegencut)
static int get_numcols(int whichtablevar, int whichtablesubtype)
{

  if(whichtablesubtype!=SUBTYPEDEGEN){
    // doesn't depend upon whichtablevar
    return(numcolintablesubtype[whichtablesubtype]); // ok use of numcolintablesubtype
  }
  else{
    // depends upon whichtablevar in general
    if(utotdegencut[whichtablevar]>DEGENCUTLASTOLDVERSION){
      return(NUMEOSDEGENQUANTITIESMEMNEW);
    }
    else{
      return(NUMEOSDEGENQUANTITIESMEMOLD);
    }
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
  // ok use of numcolintablesubtype

  numcolintablesubtype[SUBTYPEDEGEN]=MAX(NUMEOSDEGENQUANTITIESMEMNEW,NUMEOSDEGENQUANTITIESMEMOLD); // more complicated table, dealt with in special way, but can still access table quantities directly per whichd.  Choose largest range possible in case setup (e.g.) loginterp or whatnot
  numcolintablesubtype[SUBTYPESTANDARD]=NUMEOSSTANDARDQUANTITIESMEM;
  numcolintablesubtype[SUBTYPEGUESS]=NUMEOSGUESSQUANTITIESMEM;
  numcolintablesubtype[SUBTYPEDISS]=NUMEOSDISSQUANTITIESMEM;
  numcolintablesubtype[SUBTYPEDP]=NUMEOSDPQUANTITIESMEM;
  numcolintablesubtype[SUBTYPESDEN]=NUMEOSSDENQUANTITIESMEM;
  numcolintablesubtype[SUBTYPESSPEC]=NUMEOSSSPECQUANTITIESMEM;
  numcolintablesubtype[SUBTYPEPOFCHI]=NUMEOSPOFCHIQUANTITIESMEM;
  numcolintablesubtype[SUBTYPETEMP]=NUMEOSTEMPQUANTITIESMEM; // table accessed in special way based upon "whichd" during lookup
  if(WHICHDATATYPEGENERAL==4) numcolintablesubtype[SUBTYPEEXTRA]=NUMEXTRAEOSQUANTITIESTYPE4;
  else if(WHICHDATATYPEGENERAL==3) numcolintablesubtype[SUBTYPEEXTRA]=NUMEXTRAEOSQUANTITIESTYPE3;
  else if(WHICHDATATYPEGENERAL==2) numcolintablesubtype[SUBTYPEEXTRA]=NUMEXTRAEOSQUANTITIESTYPE2;
  else if(WHICHDATATYPEGENERAL==1) numcolintablesubtype[SUBTYPEEXTRA]=NUMEXTRAEOSQUANTITIESTYPE1;
  

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
  whichdintablesubtype[SUBTYPEEXTRA]=UTOTDIFF; // all extras are same whichd=UTOTDIFF

  isextraintablesubtype[SUBTYPEDEGEN]=ISNOTEXTRATABLETYPE;
  isextraintablesubtype[SUBTYPESTANDARD]=ISNOTEXTRATABLETYPE;
  isextraintablesubtype[SUBTYPEGUESS]=ISNOTEXTRATABLETYPE;
  isextraintablesubtype[SUBTYPEDISS]=ISNOTEXTRATABLETYPE;
  isextraintablesubtype[SUBTYPEDP]=ISNOTEXTRATABLETYPE;
  isextraintablesubtype[SUBTYPESDEN]=ISNOTEXTRATABLETYPE;
  isextraintablesubtype[SUBTYPESSPEC]=ISNOTEXTRATABLETYPE;
  isextraintablesubtype[SUBTYPEPOFCHI]=ISNOTEXTRATABLETYPE;
  isextraintablesubtype[SUBTYPETEMP]=ISNOTEXTRATABLETYPE;
  isextraintablesubtype[SUBTYPEEXTRA]=ISEXTRATABLETYPE;


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
  firsteosintablesubtype[SUBTYPEEXTRA]=FIRSTEOSEXTRA;

  
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
    for(coli=0;coli<numcolintablesubtype[itersubtype];coli++){ // ok use of numcolintablesubtype
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
  FTYPEEOS indep[NUMEOSINDEPS];
  FTYPEEOS indepplusdegen[NUMEOSDEGENINDEPS];
  FTYPEEOS indepdegen[NUMEOSDEGENINDEPS];
  FTYPEEOS lstepdep,diff,lindep,lindeptry;
  char headername[NUMTBLS][MAXFILENAME];
  char tablename[NUMTBLS][MAXFILENAME];
  char degentablename[NUMTBLS][MAXFILENAME];
  int tableiter;
  int tablesizeexpected[NUMTBLS][NUMEOSINDEPS];
  FTYPEEOS tabletemp[MAXNUMEOSQUANTITIESin]; // size of read-in and stored columns
  FTYPEEOS degentabletemp[NUMEOSDEGENQUANTITIESin]; // size of read-in and stored columns
  int sizematch; // 0 = table sizes don't match  1 = they match
  FTYPEEOS valuetemp[MAXNUMEOSQUANTITIESin],degenvaluetemp[NUMEOSDEGENQUANTITIESin]; // FTYPEEOS since used to process tables
  FTYPEEOS errordegen,errordegen2;
  int i;
  int testnumfscanfquantities,testnumfscanfquantitiesdegen;
  int templikeiter;
  long long int numberofdegenerrors[4],totalnumberofdegen; // 0=TABLETOL 1=0.01 2=0.1 3=TABLETOLTRUNCATION
  FTYPEEOS indexcheck;



  trifprintf("Setting up Kaz EOS table\n");


  trifprintf("SUPER TODO: derivatives need to be computed accounting for LS offsets, which has to be done in Matlab.  For now, derivatives used only for Newton's method, and works fine as is.\n");
  trifprintf("SUPER TODO: Might want to only go to 1E-10 or so for inversion if doing this EOS.\n");
  trifprintf("SUPER TODO: Make H calculation with MPI.\n");
  trifprintf("SUPER TODO: Consider larger radial domain to get H calculation to be consistent with stellar model so Ynu0 is computed better.\n");
  trifprintf("SUPER TODO: Check and turn on source code that's been turned off so far.\n");



  if(ALLOWKAZEOS==0){
    dualfprintf(fail_file,"Need to have ALLOWKAZEOS==1 if using Kaz EOS\n");
    myexit(93458635);
  }


  // special initialization of Kaz EOS lookup table repeatedeos storage
  initeos_kazfulleos();




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


    
    
    trifprintf("Sden_convertfactor=%21.15g\n",(1.0/pow(Lunit,3.0)));
    

    ///////////////////////////////////////////////
    //
    // Set file names
    //
    ///////////////////////////////////////////////

    if(ALLOWFULLTABLE){
      strcpy(headername[FULLTABLE],EOSFULLHEADNAME);
      strcpy(tablename[FULLTABLE],EOSFULLTABLENAME);
      strcpy(degentablename[FULLTABLE],EOSFULLTABLEDEGENNAME);
      
      strcpy(headername[FULLTABLEEXTRA],EOSFULLEXTRAHEADNAME);
      strcpy(tablename[FULLTABLEEXTRA],EOSFULLTABLEEXTRANAME);
      strcpy(degentablename[FULLTABLEEXTRA],EOSFULLTABLEEXTRADEGENNAME);
    }

    if(ALLOWSIMPLETABLE){
      strcpy(headername[SIMPLETABLE],EOSSIMPLEHEADNAME);
      strcpy(tablename[SIMPLETABLE],EOSSIMPLETABLENAME);
      strcpy(degentablename[SIMPLETABLE],EOSSIMPLETABLEDEGENNAME);
      
      strcpy(headername[SIMPLETABLEEXTRA],EOSSIMPLEEXTRAHEADNAME);
      strcpy(tablename[SIMPLETABLEEXTRA],EOSSIMPLETABLEEXTRANAME);
      strcpy(degentablename[SIMPLETABLEEXTRA],EOSSIMPLETABLEEXTRADEGENNAME);
    }

    if(ALLOWSIMPLEZOOMTABLE){
      strcpy(headername[SIMPLEZOOMTABLE],EOSSIMPLEZOOMHEADNAME);
      strcpy(tablename[SIMPLEZOOMTABLE],EOSSIMPLEZOOMTABLENAME);
      strcpy(degentablename[SIMPLEZOOMTABLE],EOSSIMPLEZOOMTABLEDEGENNAME);
      
      strcpy(headername[EXTRASIMPLEZOOMTABLE],EOSSIMPLEZOOMEXTRAHEADNAME);
      strcpy(tablename[EXTRASIMPLEZOOMTABLE],EOSSIMPLEZOOMEXTRATABLENAME);
      strcpy(degentablename[EXTRASIMPLEZOOMTABLE],EOSSIMPLEZOOMEXTRATABLEDEGENNAME);
    }




    ///////////////
    //
    // loop over normal tables (degen table read-in with its normal version)
    //
    //////////////

    for(tableiter=0;tableiter<NUMTBLS;tableiter++){


      // avoid tables not needed since should not expect user has copied or created such tables if not using them
      if(tableiter==FULLTABLE && ALLOWFULLTABLE==0) continue;
      if(tableiter==FULLTABLEEXTRA && (ALLOWFULLTABLE==0 || WHICHDATATYPEGENERAL!=4) ) continue;
      if(tableiter==SIMPLETABLE && ALLOWSIMPLETABLE==0) continue;
      if(tableiter==SIMPLETABLEEXTRA && (ALLOWSIMPLETABLE==0 || WHICHDATATYPEGENERAL!=4) ) continue;
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
      fscanf(inhead,"%d %d %d %d",&whichrnpmethod[tableiter],&whichynumethod[tableiter],&whichhcmmethod[tableiter],&whichyelooptype[tableiter]);
      fscanf(inhead,"%d %d %d %d",&whichdatatype[tableiter],&utotdegencut[tableiter],&numc[tableiter],&numextras[tableiter]);




      ////////////////////////////
      //
      // perform some checks
      //
      ////////////////////////////

      if(whichrnpmethod[tableiter]==0){
        dualfprintf(fail_file,"This method is not setup\n");
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
        fscanf(inhead,EOSHEADERONEIN,&inputtablelimits[tableiter][ii][0]); // start in log
        fscanf(inhead,EOSHEADERONEIN,&inputtablelimits[tableiter][ii][1]); // end in log
        fscanf(inhead,EOSHEADERONEIN,&inputtablelimits[tableiter][ii][2]); // step in log
        fscanf(inhead,EOSHEADERONEIN,&inputtablelimits[tableiter][ii][4]); // base of log offset
        fscanf(inhead,EOSHEADERONEIN,&inputtablelimits[tableiter][ii][5]); // linear offset

        if(inputtablelimits[tableiter][ii][5]>=0.0 && ii==YNUEOS){
          dualfprintf(fail_file,"SUPERWARNING: Note that table setup so Ynu (or Ynu0)<=0 not allowed since log-indexing Ynu (or Ynu0) as independent variable.\n");
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
          tablelimits[tableiter][ii][2] = (tablelimits[tableiter][ii][1]-tablelimits[tableiter][ii][0])/((FTYPEEOS)tablesize[tableiter][ii]-1.0);
      
          // the below definition is consistent with Kaz's code, matlab eos_extract.m and elsewhere in this code
          // tablelimits[tableiter][ii][3] = ((FTYPEEOS)tablesize[tableiter][ii]-1.0)/(tablelimits[tableiter][ii][1] - tablelimits[tableiter][ii][0]);


          // i = (log_b (r-r_0) - x_in)/dx
          // 1/dx:=
          tablelimits[tableiter][ii][3] = 1.0/tablelimits[tableiter][ii][2];
          // but avoid nan or inf from division by zero if no such dimension
          if(!isfinite(tablelimits[tableiter][ii][3])) tablelimits[tableiter][ii][3]=0.0;

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
      fscanf(inhead,EOSHEADERONEIN,&lsoffset[tableiter]);
      fscanf(inhead,EOSHEADERONEIN,&fakelsoffset[tableiter]);
      fscanf(inhead,EOSHEADERONEIN,&fakeentropylsoffset[tableiter]);



      // GODMARK: should really read-in from table, although expected tabular value is 9.14Mev and this gives 7.108 for smooth connection for initial conditions
      // GODMARK: The below is for Shen EOS only
      FAKE2IDEALNUCLEAROFFSET[tableiter] = (-7.57E-3);




      /////////////
      //
      // Convert "LS" quantities from LS format to KAZ to Matlab format
      //
      /////////////



      // Exactly 9.14MeV/baryon as in helm/jon_lsbox.f
      // This offset was subtracted before tabulating in log-log.  After lookup, this energy/baryon needs to be added back in the correct units.
      // ergPmev=1.60218E-6
      //  TRUENUCLEAROFFSET=(9.14*ergPmev/(mb*C*C));// cgs/cgs = dimensionless quantity
      // above gives 9.7346E-3
      //  TRUENUCLEAROFFSET /= (1.0); // need to get code version of energy/baryon however used to add to internal energy
      // assume degen offset accounts for offset, that cannot be put into EOS itself!
      // Now read-in unphysical fakelsoffset used to avoid negative energies is required to reoffset internal energy back so only physical "lsoffset" is included
      //      TRUENUCLEAROFFSET[tableiter]=fakelsoffset[tableiter]*ergPmev/energyunit;
     
      // i.e. (u/(\rho_0 c^2)) = (u/c^2)/(\rho_0) = (u/c^2)/(m_b n_b) = (u/n_b) 1/(mb c^2), so always divide by $m_b c^2$ in cgs units if original in cgs units
      TRUENUCLEAROFFSET[tableiter]=(fakelsoffset[tableiter]*ergPmev)/(mb*C*C);

      // fakeentropylsoffset is in dimensionless entropy per baryon=Ss where Sden=1/cc
      // Ssjon = SJden/(rho0 c^2) = Ssdimenless nb/(rho0 c^2) = Ssdimenless/(mb c^2)
      TRUEENTROPYNUCLEAROFFSET[tableiter]=fakeentropylsoffset[tableiter]/(mb*C*C);

      // degeneracy global offset:
      // ergPmev=1.60218E-6
      //  DEGENNUCLEAROFFSET=(50.0*ergPmev/(mb*C*C));// cgs/cgs = dimensionless quantity
      DEGENNUCLEAROFFSET[tableiter]=0.0; // avoid this approach for now -- try putting in offset into original HELM table first



      /////////////
      //
      // Convert to code units
      //
      /////////////
      TRUENUCLEAROFFSET[tableiter]*=1.0; // already dimensionless

      // nearly dimensionless.  Treat like HARM meaning of specific entropy = Ssjon = Sden/(rho0 c^2)
      // So treat like SSofCHIin
      // Then will be like SSofCHIin = Sden/(rho0*c^2)
      TRUEENTROPYNUCLEAROFFSET[tableiter]/=(1.0/energyunit);



      // Table has Ss = Sden/(rho c^2) always.
      // if rho0unittype==0, HARM has specific entropy in Sden/(rho c^2) = Ss/(mb c^2), so just divide by (mb c^2)
      // if rho0unittype==1, HARM has specific entropy in Sden/(rho) = Ss/(mb), so just divide by (mb)
      // fakeentropylsoffset is in dimensionless entropy per baryon.  Need to convert to HARM Ss
      // Ss[harm] = Sden/(rho_0 c^2) = Sden/(mb c^2 nb) = Ss[dimenless]/(mb c^2)
      // Ss[harm] = Sden/rho_0 = Sden/(mb nb) = Ss[dimenless]/mb


      //
      ////////////////////////




      ////////////////////////
      //
      // Also read in ye grid parameters
      fscanf(inhead,EOSHEADERONEIN,&eosyegrid1[tableiter]);
      fscanf(inhead,EOSHEADERONEIN,&eosyegrid2[tableiter]);
      fscanf(inhead,EOSHEADERONEIN,&eosxgrid1[tableiter]);
      fscanf(inhead,EOSHEADERONEIN,&eosxgrid2[tableiter]);


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
#define KAZDATAREADLOOP for(mmm=0;mmm<tablesize[tableiter][vartypearray[5]];mmm++)for(lll=0;lll<tablesize[tableiter][vartypearray[4]];lll++)for(kkk=0;kkk<tablesize[tableiter][vartypearray[3]];kkk++)for(jjj=0;jjj<tablesize[tableiter][vartypearray[2]];jjj++)for(iii=0;iii<tablesize[tableiter][vartypearray[1]];iii++)


      KAZDATAREADLOOP{
        testnumfscanfquantities=0;
        testnumfscanfquantitiesdegen=0;

        // first, get positions to make sure everything is consistent
        fscanf(intable,"%d %d %d %d %d",&m,&n,&o,&p,&q); // indexing related to NUMINDEPDIMENS not NUMEOSINDEPS
        testnumfscanfquantities += NUMINDEPDIMENS;

        if(m!=iii || n!=jjj || o!=kkk || p!=lll || q!=mmm){
          dualfprintf(fail_file,"Read-in table (%d) indicies inconsistent with expected indicies: m=%d iii=%d n=%d jjj=%d o=%d kkk=%d p=%d lll=%d q=%d mmm=%d\n",tableiter,m,iii,n,jjj,o,kkk,p,lll,q,mmm);
          dualfprintf(fail_file,"whichrnpmethod=%d whichynumethod=%d whichhcmmethod=%d whichyelooptype=%d\n",whichrnpmethod[tableiter],whichynumethod[tableiter],whichhcmmethod[tableiter],whichyelooptype[tableiter]);
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
          fscanf(intable,EOSHEADERONEIN,&indep[ii]); // rhob, Udiff, Pdiff, CHIdiff, Sdiff, tdynorye, tdynorynu, H   for given grid value
          testnumfscanfquantities += 1;
        }
        // read-in UofUdiff, PofPdiff, CHIofCHIdiff, SofSdiff -- used to check degen offset calculation in HARM
        for(ii=0;ii<NUMEOSDEGENQUANTITIESMEM1;ii++){
          fscanf(intable,EOSHEADERONEIN,&indepplusdegen[ii]);
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
            fscanf(indegentable,EOSHEADERONEIN,&indepdegen[ii]); // rho, tdynorye, tdynorynu, H
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



            if(ii==YEEOS && whichyelooptype[tableiter]==YELOOPTYPESPECIAL){
              // then special gridding
              lookup_yespecial(tableiter,vartypearray[YEINDEP],indep[ii],&indexcheck);
              diff=(fabs(indexcheck-(FTYPEEOS)kkk))/(fabs(indexcheck)+fabs((FTYPEEOS)kkk));
              if(diff>TABLETOLYEINDEX){
                dualfprintf(fail_file,"Special Ye grid position not correct: diff=%21.15g kkk=%d indexcheck=%21.15g\n",diff,kkk,indexcheck);
                myexit(978351235);
              }
            }
            else{
     
              // get step (consistent with how step is computed in Kaz's code and in matlab script eos_extract.m)
              // really only has to be consistent with eos_extract.m
              //     lstepdep = (-tablelimits[tableiter][ii][0])/((FTYPEEOS)tablesize[tableiter][ii]-1.0);
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
              //     diff = fabs(lindep-lindeptry)/(fabs(lindep)+fabs(lindeptry));
              // normalize below by range of limits instead of local values since otherwise if lindeptry or lindep is near 0 then relative error erroneously is large with old diff above
              diff = fabs(lindep-lindeptry)/(tablelimits[tableiter][ii][1]-tablelimits[tableiter][ii][0]);
              if(diff>TABLETOL){
                dualfprintf(fail_file,"Grid position data is incorrect: mmm=%d lll=%d kkk=%d jjj=%d iii=%d :: ii=%d readin-lindep=%21.15g lindeptry=%21.15g diff=%21.15g\n",mmm,lll,kkk,jjj,iii,ii,lindep,lindeptry,diff);
                dualfprintf(fail_file,"tablelimits[%d][%d][0]=%21.15g totalindex[%d]=%d lstepdep=%21.15g\n",tableiter,ii,tablelimits[tableiter][ii][0],ii,totalindex[ii],lstepdep);
                dualfprintf(fail_file,"indep=%21.15g lindep=%21.15g [3]=%21.15g [2]=%21.15g\n",indep[ii],lindep,lineartablelimits[tableiter][ii][3],lineartablelimits[tableiter][ii][2]);
                dualfprintf(fail_file,"Trigger: %d\n",(utotdegencut[tableiter]>DEGENCUTLASTOLDVERSION)&&(ii>=FIRSTTKLIKE && ii<=LASTTKLIKE));
                myexit(16628);
              }

            }// end else if not yeindep with special loop type
          }
        }




        ////////////////////////
        //
        // fourth, since table gridding is consistent, now read in columns of actual dependent variable values
        //
        ////////////////////////


        // normal table
        for(ppp=0;ppp<numeosquantitiesinfile[tableiter];ppp++){ // look over only to-be stored quantities
          fscanf(intable,EOSHEADERONEIN,&valuetemp[ppp]); // FTYPEEOS values
          testnumfscanfquantities += 1;
        }


        /// DEBUG:
        // if(isextratable(tableiter)){
        //   if(mmm==0 && lll==3 && kkk==46 && jjj==11 && iii==31){
        //     for(ppp=0;ppp<numeosquantitiesinfile[tableiter];ppp++){
        //       dualfprintf(fail_file,"in file checking valuetemp[%d]=%21.15g\n",ppp,valuetemp[ppp]);
        //     }
        //   }
        // }


        // before storing, reorder to be like standard "in" format for EOS quantities:
        reorder_extratype(isextratable(tableiter),numeosquantitiesinfile[tableiter],numeosquantitiesinstandardlist[tableiter],valuetemp);
        // not in canonical list order/position

        for(ppp=0;ppp<numeosquantitiesinstandardlist[tableiter];ppp++){ // look over only to-be stored quantities
          set_arrays_eostable(ISNOTDEGENTABLE,tableiter,mmm,lll,kkk,jjj,iii,ppp,valuetemp[ppp]);
        }


        // DEBUG:
        // if(isextratable(tableiter)){
        //   if(mmm==0 && lll==3 && kkk==46 && jjj==11 && iii==31){
        //     for(ppp=0;ppp<numeosquantitiesinstandardlist[tableiter];ppp++){ // look over only to-be stored quantities
        //       dualfprintf(fail_file,"checking valuetemp[%d]=%21.15g\n",ppp,valuetemp[ppp]);
        //     }
        //   }
        // }



        // degen table
        if(jjj==0){
          for(ppp=0;ppp<numeosdegenquantitiesinfile[tableiter];ppp++){
            fscanf(indegentable,EOSHEADERONEIN,&degenvaluetemp[ppp]); // FTYPEEOS values

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

          FTYPEEOS U0,UIN,UOUT,NN,testindepplusdegen;
          if(utotdegencut[tableiter]>DEGENCUTLASTOLDVERSION){ // all temperature-like quantities here
            // access read-in quantities
            U0=degenvaluetemp[ii];
            UIN=degenvaluetemp[ii+NUMEOSDEGENQUANTITIESMEM1];
            UOUT=degenvaluetemp[ii+2*NUMEOSDEGENQUANTITIESMEM1];
            NN = tablesizeexpected[tableiter][UEOS]; // UEOS,PEOS,CHIEOS,SEOS all same dimension
            // already checked that indep is consistent, so can use indep[RHOEOS+ii] or (FTYPEEOS)jjj/NN below for "i/N"
            testindepplusdegen = U0 + (UIN-U0)*pow( (UOUT-U0)/(UIN-U0), (FTYPEEOS)jjj/NN);
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

#if(PRODUCTION==0&&0) /// DEATHMARK -- should probably see if can see why not accurate for extra table -- GODMARK SUPERGODMARK
          if(errordegen>=TABLETOLTRUNCATION){
            dualfprintf(fail_file,"Degen not correct: iii=%d jjj=%d kkk=%d lll=%d mmm=%d:: ii=%d :: error=%21.15g :: U0=%21.15g  UIN=%21.15g UOUT=%21.15g\n",iii,jjj,kkk,lll,mmm,ii,errordegen,U0,UIN,UOUT);
            dualfprintf(fail_file,"file has indepplusdegen=%21.15g while testindepplusdegen=%21.15g\n",indepplusdegen[ii],testindepplusdegen);
          }
#endif

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
   

      if(TBLLINEARITEMS>TBLITEMS){
        dualfprintf(fail_file,"Conversion code assumes TBLLINEARITEMS=%d <= TBLITEMS=%d\n",TBLLINEARITEMS,TBLITEMS);
        myexit(19825625);
      }
 
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
      for(jj=0;jj<NUMEOSINDEPS;jj++){
        tablelimits[tableiter][jj][3] = ((FTYPEEOS)tablesize[tableiter][jj]-1.0)/(tablelimits[tableiter][jj][1] - tablelimits[tableiter][jj][0]);
        if(!isfinite(tablelimits[tableiter][jj][3])) tablelimits[tableiter][jj][3]=0.0; // if no such dimension, then reduce step to 0 instead of infinity
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

      // UEOS gives same size as PEOS, CHIEOS, and SEOS
      KAZDATAREADLOOP{

        ////////////////
        //
        // temp set
        //
        ///////////////
 
        // normal table
        for(jj=0;jj<numeosquantitiesinstandardlist[tableiter];jj++) get_arrays_eostable(ISNOTDEGENTABLE,tableiter,mmm,lll,kkk,jjj,iii,jj,&tabletemp[jj]);
        // degen table
        if(jjj==0) for(jj=0;jj<numeosdegenquantitiesinstandardlist[tableiter];jj++) get_arrays_eostable(ISDEGENTABLE,tableiter,mmm,lll,kkk,jjj,iii,jj,&degentabletemp[jj]);



        if(
           (WHICHDATATYPEGENERAL!=4) // then no "extra" tables to check for -- all in table
           || (WHICHDATATYPEGENERAL==4 && isextratable(tableiter)==0 && mmm==0 && lll==0)
           ){

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






          ///////////////////////
          //
          // Check values of some quantities
          // Like derivatives, that can be wrong in Newton's method, but not very wrong
          // Compare with constraints placed in eos_extract.m
          // only check if temperature not out of bounds.  At this point, temperature still not converted to code units
          //
          ///////////////////////
      
      
          if(tabletemp[TEMPUin]>1 && (tabletemp[DPDUofRHOUin]<0.01 || tabletemp[DPDUofRHOUin]>2.0)){
            dualfprintf(fail_file,"Probably bad dPdU|rho = %21.15g :: tableiter=%d mmm=%d lll=%d kkk=%d jjj=%d iii=%d\n",tabletemp[DPDUofRHOUin],tableiter,mmm,lll,kkk,jjj,iii,jj);
          }

          if(tabletemp[TEMPCHIin]>1 && (tabletemp[IDCHIDPin]<0.01 || tabletemp[IDCHIDPin]>1.0)){
            dualfprintf(fail_file,"Probably bad dPd\\chi|rho = %21.15g :: tableiter=%d mmm=%d lll=%d kkk=%d jjj=%d iii=%d\n",tabletemp[IDCHIDPin],tableiter,mmm,lll,kkk,jjj,iii,jj);
          }

          // unsure how to check entropy.  Just make sure not very small?
          if(tabletemp[TEMPUin]>1 && (tabletemp[DSDUofRHOUin]<SMALL )){
            dualfprintf(fail_file,"Probably bad dSdendU|rho = %21.15g :: tableiter=%d mmm=%d lll=%d kkk=%d jjj=%d iii=%d\n",tabletemp[DSDUofRHOUin],tableiter,mmm,lll,kkk,jjj,iii,jj);
          }

          if(tabletemp[TEMPCHIin]>1 && (tabletemp[DSSDCHIofRHOCHIin]<SMALL )){
            dualfprintf(fail_file,"Probably bad dSSdCHI|rho = %21.15g :: tableiter=%d mmm=%d lll=%d kkk=%d jjj=%d iii=%d\n",tabletemp[DSSDCHIofRHOCHIin],tableiter,mmm,lll,kkk,jjj,iii,jj);
          }

          // below is (1/cc) / (erg/cc) \propto 1/erg since we *input* rho as rho c^2
          //tabletemp[DSDRHOofRHOUin]/=(1.0/energyunit); 
          //      /=(1.0/energyunit);
      
          // DSSDRHOofRHOCHI, DSSDCHIofRHOCHI have units of Ss/(rho0 c^2) since Ss has units of 1/(mb c^2) and input rho with units of rho*c^2
          //      tabletemp[DSSDRHOofRHOCHIin]/=(1.0/(energyunit*Pressureunit));
          //      tabletemp[DSSDCHIofRHOCHIin]/=(1.0/(energyunit*Pressureunit));




        }


        if(
           (WHICHDATATYPEGENERAL!=4) // then no "extra" tables to check for -- all in table
           || (WHICHDATATYPEGENERAL==4 && mmm==0 && lll==0)
           ){

          // *always* have temperature in table
          // TEMP used for table validity check
          tabletemp[TEMPUin]/=Tempunit;
          tabletemp[TEMPPin]/=Tempunit;
          tabletemp[TEMPCHIin]/=Tempunit;
          tabletemp[TEMPSin]/=Tempunit;
   
        } 



        if(
           (WHICHDATATYPEGENERAL!=4) // then no "extra" tables to check for -- all in table
           || (WHICHDATATYPEGENERAL==4 && isextratable(tableiter)==1 && mmm==0)
           ){

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
            //   \Gamma = 1/s
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
            // no conversion for Ynuthermal0 that's dimensionless EXTRA12 GODMARK: If used, should add!
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
            //   FTYPEEOS fakerho = exp(0.0+1.0*iii);
            // tabletemp[EXTRA22in]=1.0*pow(fakerho,3.0);
            //dualfprintf(fail_file,"GOD: iii=%d %21.15g\n",iii,tabletemp[EXTRA22in]);

            // thermal number densities
            tabletemp[EXTRA23in]/=(1.0/pow(Lunit,3.0));
            tabletemp[EXTRA24in]/=(1.0/pow(Lunit,3.0));
          }

 
        }// end if to convert extras


        ////////////////////
        //
        // convert units for degen table
        //
        ////////////////////
        if(
           (WHICHDATATYPEGENERAL!=4 && jjj==0)
           || (WHICHDATATYPEGENERAL==4 && jjj==0 && mmm==0 && lll==0) // degen table then never function of T-like, Ynu, or H
           ){

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

        }// end if to convert degens


        ////////////////
        //
        // temp restore back into arrays (can go over all quantities even if don't exist for a given table since taken care of inside.  However, can't convert existing table values multiple times, so why avoided certain indices above.)
        //
        ///////////////

        // normal table
        for(jj=0;jj<numeosquantitiesinstandardlist[tableiter];jj++) set_arrays_eostable(ISNOTDEGENTABLE,tableiter,mmm,lll,kkk,jjj,iii,jj,tabletemp[jj]);
        // degen table
        if(jjj==0) for(jj=0;jj<numeosdegenquantitiesinstandardlist[tableiter];jj++) set_arrays_eostable(ISDEGENTABLE,tableiter,mmm,lll,kkk,jjj,iii,jj,degentabletemp[jj]);



      } // end loop over table



            
      // report conversion
      trifprintf("Done converting units of EOS table: tableiter=%d of %d\n",tableiter,NUMTBLS-1);



      // compute pre-logify code version of invalid temperature
      invalidtempcode=INVALIDTEMP/Tempunit;



      

     
      ///////////////////////////////
      //
      // Pre-logify
      //
      ///////////////////////////////
      if(DOLOGINTERP && DOPRELOGIFY){

        long long int goodlog[MAXEOSPIPELINE],badlog[MAXEOSPIPELINE],totalall[MAXEOSPIPELINE];
        long long int gooddegenlog[MAXEOSPIPELINE],baddegenlog[MAXEOSPIPELINE],totaldegenall[MAXEOSPIPELINE];
        FTYPEEOS lowestpoint;

        for(jj=0;jj<numeosquantitiesinstandardlist[tableiter];jj++){ goodlog[jj]=badlog[jj]=totalall[jj]=0; }
        for(jj=0;jj<numeosdegenquantitiesinstandardlist[tableiter];jj++){ gooddegenlog[jj]=baddegenlog[jj]=totaldegenall[jj]=0; }


        //////////
        // LOOP
        //////////
        KAZDATAREADLOOP{
   
          // normal table
          for(jj=0;jj<numeosquantitiesinstandardlist[tableiter];jj++){
            // first check if quantity even really exists (i.e. was read-in and really stored)
            // can't multi-logify (i.e. temperature does not depend upon Ynu,H unlike extras that depend upon Ynu for ==4)
            if(
               (WHICHDATATYPEGENERAL!=4) // then no "extra" tables to check for -- all in table
               || (WHICHDATATYPEGENERAL==4 && isextratable(tableiter)==0 && (jj>=FIRSTEOSQUANTITIESBASEin && jj<=LASTEOSQUANTITIESBASEin) && (mmm==0 && lll==0)) // can't multi-logify!
               || (WHICHDATATYPEGENERAL==4 && isextratable(tableiter)==1 && (jj>=EXTRA1in && jj<=EXTRA24in && mmm==0 || jj>=TEMPUin && jj<=TEMPSin && mmm==0 && lll==0)) // no H-dep if ==4 even for extra table
               ){
              totalall[jj]++;
              get_arrays_eostable(ISNOTDEGENTABLE,tableiter,mmm,lll,kkk,jjj,iii,jj,&tabletemp[jj]);
              if(get_dologinterp_in(ISNOTDEGENTABLE, jj)){

                if(jj>=TEMPUin && jj<=TEMPSin){

                  lowestpoint=invalidtempcode; // after units conversion
                  if(tabletemp[jj]>lowestpoint){
                    tabletemp[jj]=log10(tabletemp[jj]);
                    goodlog[jj]++;
                  }
                  else{
                    tabletemp[jj]=log10(invalidtempcode); // indicates invalidate temperature -- will end up much much smaller than any normal temperature
                    badlog[jj]++;
                  }// end else if bad log
                }
                else{
                  lowestpoint=+0.0;
                  if(tabletemp[jj]>lowestpoint){
                    tabletemp[jj]=log10(tabletemp[jj]);
                    goodlog[jj]++;
                  }
                  else{
                    tabletemp[jj]=OUTOFBOUNDSPRELOGIFY;
                    badlog[jj]++;
                  }// end else if bad log
                }
              }// end if need to do log
              set_arrays_eostable(ISNOTDEGENTABLE,tableiter,mmm,lll,kkk,jjj,iii,jj,tabletemp[jj]);
            }
          }
   
          // degen table
          if(jjj==0 && mmm==0 && lll==0){ // can't multi-logify
            for(jj=0;jj<numeosdegenquantitiesinstandardlist[tableiter];jj++){
              // all tables have all degen quantities

              totaldegenall[jj]++;
              get_arrays_eostable(ISDEGENTABLE,tableiter,mmm,lll,kkk,jjj,iii,jj,&degentabletemp[jj]);
              if(get_dologinterp_in(ISDEGENTABLE, jj)){
                lowestpoint=+0.0;
                if(degentabletemp[jj]>lowestpoint){
                  degentabletemp[jj]=log10(degentabletemp[jj]);
                  gooddegenlog[jj]++;
                }
                else{
                  degentabletemp[jj]=OUTOFBOUNDSPRELOGIFY;
                  baddegenlog[jj]++;
                }
              }
              set_arrays_eostable(ISDEGENTABLE,tableiter,mmm,lll,kkk,jjj,iii,jj,degentabletemp[jj]);
            }// end over jj
          }// end if jjj==0
   
   
        }// end over prelogify grid loop


        // compute code version of invalid temperature
        invalidtempcode=log10(invalidtempcode); // unlike other quantities, check if temperature is *below* certain limit.  Just keep to log of original limit.


        trifprintf("Done with pre-logify\n");
        for(jj=0;jj<numeosquantitiesinstandardlist[tableiter];jj++) trifprintf("non-degen[%d]: bad=%lld good=%lld total=%lld\n",jj,badlog[jj],goodlog[jj],totalall[jj]);
        for(jj=0;jj<numeosdegenquantitiesinstandardlist[tableiter];jj++) trifprintf("degen[%d]: bad=%lld good=%lld total=%lld\n",jj,baddegenlog[jj],gooddegenlog[jj],totaldegenall[jj]);

 
      }// end if prelogifying










      trifprintf("tableiter=%d invalidtempcode=%21.15g Tempunit=%21.15g\n",tableiter,invalidtempcode,Tempunit);


    }// end loop over tables



  } // end if myid==0
  



  ///////////
  //
  // Set some global variables used by rest of non-EOS code
  //
  ////////////
  TRUENUCLEAROFFSETPRIMARY=TRUENUCLEAROFFSET[primarytable];


  ///////////////////////////////////////////////
  //
  // set maximum rho in table
  //
  ///////////////////////////////////////////////
  
  rhoupperlimit=-BIG;
  for(tableiter=0;tableiter<NUMTBLS;tableiter++){
    if(tableiter==FULLTABLE && ALLOWFULLTABLE==1 || tableiter==SIMPLETABLE && ALLOWSIMPLETABLE==1 || tableiter==SIMPLEZOOMTABLE && ALLOWSIMPLEZOOMTABLE==1){
      rhoupperlimit = MAX(rhoupperlimit,lineartablelimits[tableiter][RHOEOS][1]);
    }
  }



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

  MPI_Bcast(&inputtablelimits[0][0][0],NUMTBLS*NUMEOSINDEPS*TBLITEMS,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&tablelimits[0][0][0],NUMTBLS*NUMEOSINDEPS*TBLITEMS,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&lineartablelimits[0][0][0],NUMTBLS*NUMEOSINDEPS*TBLLINEARITEMS,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&tablesize[0][0],NUMTBLS*NUMEOSINDEPS,MPI_INT,MPIid[0],MPI_COMM_GRMHD);

  MPI_Bcast(&vartypeeosextraarray[0],NUMINDEPDIMENSMEM,MPI_INT,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&vartypearray[0],NUMINDEPDIMENSMEM,MPI_INT,MPIid[0],MPI_COMM_GRMHD);

  MPI_Bcast(&vartypeheightarray[0],NUMHDIRECTIONS+1,MPI_INT,MPIid[0],MPI_COMM_GRMHD);

  MPI_Bcast(&vardegentypearray[0],NUMEOSDEGENINDEPS+1,MPI_INT,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&varnormalcompare2degentypearray[0],NUMEOSDEGENINDEPS+1,MPI_INT,MPIid[0],MPI_COMM_GRMHD);


  MPI_Bcast(&numcolintablesubtype[0],NUMTABLESUBTYPES,MPI_INT,MPIid[0],MPI_COMM_GRMHD); // ok use of numcolintablesubtype
  MPI_Bcast(&whichdintablesubtype[0],NUMTABLESUBTYPES,MPI_INT,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&firsteosintablesubtype[0],NUMTABLESUBTYPES,MPI_INT,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&isextraintablesubtype[0],NUMTABLESUBTYPES,MPI_INT,MPIid[0],MPI_COMM_GRMHD);

  MPI_Bcast(&whichtablesubtypeinquantity[0],NUMEOSQUANTITIESMEM,MPI_INT,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&whichcolinquantity[0],NUMEOSQUANTITIESMEM,MPI_INT,MPIid[0],MPI_COMM_GRMHD);

  MPI_Bcast(&dologinterp_sub_coli[0][0],NUMTABLESUBTYPES*MAXEOSPIPELINE,MPI_INT,MPIid[0],MPI_COMM_GRMHD);

  MPI_Bcast(&invalidtempcode,1,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&invalidlogtempcode,1,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);

  MPI_Bcast(&whichrnpmethod[0],NUMTBLS,MPI_INT,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&whichynumethod[0],NUMTBLS,MPI_INT,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&whichhcmmethod[0],NUMTBLS,MPI_INT,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&whichyelooptype[0],NUMTBLS,MPI_INT,MPIid[0],MPI_COMM_GRMHD);

  MPI_Bcast(&whichdatatype[0],NUMTBLS,MPI_INT,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&utotdegencut[0],NUMTBLS,MPI_INT,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&numc[0],NUMTBLS,MPI_INT,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&numextras[0],NUMTBLS,MPI_INT,MPIid[0],MPI_COMM_GRMHD);

  MPI_Bcast(&primarytable,1,MPI_INT,MPIid[0],MPI_COMM_GRMHD);


  MPI_Bcast(&FAKE2IDEALNUCLEAROFFSET,NUMTBLS,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&TRUENUCLEAROFFSET,NUMTBLS,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&TRUEENTROPYNUCLEAROFFSET,NUMTBLS,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&TRUENUCLEAROFFSETPRIMARY,1,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&DEGENNUCLEAROFFSET,NUMTBLS,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);

  MPI_Bcast(&lsoffset,NUMTBLS,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&fakelsoffset,NUMTBLS,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&fakeentropylsoffset,NUMTBLS,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);

  MPI_Bcast(&eosyegrid1,NUMTBLS,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&eosyegrid2,NUMTBLS,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&eosxgrid1,NUMTBLS,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&eosxgrid2,NUMTBLS,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);

  MPI_Bcast(&rhoupperlimit,1,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);

  MPI_Bcast(&didsetupkazeos,1,MPI_INT,MPIid[0],MPI_COMM_GRMHD);




  /////////////////////////////
  //
  // kazfulleos.eostablesdefs.h globals:
  //
  /////////////////////////////


  // Bcast table data
  // eostable is FTYPEEOS, so MPI_FTYPEEOS
  // to generate below, take part of kazfull.eostablesdefs.h within WHICHEOS==KAZFULL and regexp:
  // FTYPEEOS BASEEOSMAC(\([a-zA-Z0-9]+\),\([a-zA-Z0-9]+\),\([a-zA-Z0-9]+\),\([a-zA-Z0-9]+\),\([a-zA-Z0-9]+\),\([a-zA-Z0-9]+\),\([a-zA-Z0-9]+\),\([a-zA-Z0-9]+\)); ->
  // MPI_Bcast(&(BASEEOSMAC(\1,0,0,0,0,0,0,0)),\2*\3*\4*\5*\6*\7*\8,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  // NOTE!  We are using BASEEOSMAC() so can just use BASEEOSMAC(name,0,0,0,0,0,0,0) instead of having to add shift

#if(WHICHEOS==KAZFULL)

  // full
  MPI_Bcast(&(BASEEOSMAC(eosfulltabledegen,0,0,0,0,0,0,0)),NUMEOSDEGENQUANTITIESMEM1*EOSFULLDEGENN5*EOSFULLDEGENN4*EOSFULLDEGENN3*EOSFULLDEGENN2*EOSFULLDEGENN1*NUMEOSDEGENQUANTITIESMEM2,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&(BASEEOSMAC(eosfulltablestandard,0,0,0,0,0,0,0)),1*EOSFULLN5*EOSFULLN4*EOSFULLN3*EOSFULLN2*EOSFULLN1*NUMEOSSTANDARDQUANTITIESMEM,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&(BASEEOSMAC(eosfulltableguess,0,0,0,0,0,0,0)),1*EOSFULLN5*EOSFULLN4*EOSFULLN3*EOSFULLN2*EOSFULLN1*NUMEOSGUESSQUANTITIESMEM,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&(BASEEOSMAC(eosfulltablediss,0,0,0,0,0,0,0)),1*EOSFULLN5*EOSFULLN4*EOSFULLN3*EOSFULLN2*EOSFULLN1*NUMEOSDISSQUANTITIESMEM,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&(BASEEOSMAC(eosfulltabledp,0,0,0,0,0,0,0)),1*EOSFULLN5*EOSFULLN4*EOSFULLN3*EOSFULLN2*EOSFULLN1*NUMEOSDPQUANTITIESMEM,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&(BASEEOSMAC(eosfulltablesden,0,0,0,0,0,0,0)),1*EOSFULLN5*EOSFULLN4*EOSFULLN3*EOSFULLN2*EOSFULLN1*NUMEOSSDENQUANTITIESMEM,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&(BASEEOSMAC(eosfulltablesspec,0,0,0,0,0,0,0)),1*EOSFULLN5*EOSFULLN4*EOSFULLN3*EOSFULLN2*EOSFULLN1*NUMEOSSSPECQUANTITIESMEM,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&(BASEEOSMAC(eosfulltablepofchi,0,0,0,0,0,0,0)),1*EOSFULLN5*EOSFULLN4*EOSFULLN3*EOSFULLN2*EOSFULLN1*NUMEOSPOFCHIQUANTITIESMEM,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&(BASEEOSMAC(eosfulltabletemp,0,0,0,0,0,0,0)),NUMEOSTEMPQUANTITIESMEM1*EOSFULLN5*EOSFULLN4*EOSFULLN3*EOSFULLN2*EOSFULLN1*NUMEOSTEMPQUANTITIESMEM2,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);

  // different size for extra table
  MPI_Bcast(&(BASEEOSMAC(eosfulltableextra,0,0,0,0,0,0,0)),1*EOSFULLEXTRAN5*EOSFULLEXTRAN4*EOSFULLEXTRAN3*EOSFULLEXTRAN2*EOSFULLEXTRAN1*NUMEOSEXTRAQUANTITIESMEM,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
#if(WHICHDATATYPEGENERAL==4)
  MPI_Bcast(&(BASEEOSMAC(eosfulltableextradegen,0,0,0,0,0,0,0)),NUMEOSDEGENQUANTITIESMEM1*EOSFULLEXTRADEGENN5*EOSFULLEXTRADEGENN4*EOSFULLEXTRADEGENN3*EOSFULLEXTRADEGENN2*EOSFULLEXTRADEGENN1*NUMEOSDEGENQUANTITIESMEM2,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&(BASEEOSMAC(eosfulltableextratemp,0,0,0,0,0,0,0)),NUMEOSTEMPQUANTITIESMEM1*EOSFULLEXTRATEMPN5*EOSFULLEXTRATEMPN4*EOSFULLEXTRATEMPN3*EOSFULLEXTRATEMPN2*EOSFULLEXTRATEMPN1*NUMEOSTEMPQUANTITIESMEM2,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
#endif

  // simple
  MPI_Bcast(&(BASEEOSMAC(eossimpletabledegen,0,0,0,0,0,0,0)),NUMEOSDEGENQUANTITIESMEM1*EOSSIMPLEDEGENN5*EOSSIMPLEDEGENN4*EOSSIMPLEDEGENN3*EOSSIMPLEDEGENN2*EOSSIMPLEDEGENN1*NUMEOSDEGENQUANTITIESMEM2,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&(BASEEOSMAC(eossimpletablestandard,0,0,0,0,0,0,0)),1*EOSSIMPLEN5*EOSSIMPLEN4*EOSSIMPLEN3*EOSSIMPLEN2*EOSSIMPLEN1*NUMEOSSTANDARDQUANTITIESMEM,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&(BASEEOSMAC(eossimpletableguess,0,0,0,0,0,0,0)),1*EOSSIMPLEN5*EOSSIMPLEN4*EOSSIMPLEN3*EOSSIMPLEN2*EOSSIMPLEN1*NUMEOSGUESSQUANTITIESMEM,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&(BASEEOSMAC(eossimpletablediss,0,0,0,0,0,0,0)),1*EOSSIMPLEN5*EOSSIMPLEN4*EOSSIMPLEN3*EOSSIMPLEN2*EOSSIMPLEN1*NUMEOSDISSQUANTITIESMEM,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&(BASEEOSMAC(eossimpletabledp,0,0,0,0,0,0,0)),1*EOSSIMPLEN5*EOSSIMPLEN4*EOSSIMPLEN3*EOSSIMPLEN2*EOSSIMPLEN1*NUMEOSDPQUANTITIESMEM,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&(BASEEOSMAC(eossimpletablesden,0,0,0,0,0,0,0)),1*EOSSIMPLEN5*EOSSIMPLEN4*EOSSIMPLEN3*EOSSIMPLEN2*EOSSIMPLEN1*NUMEOSSDENQUANTITIESMEM,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&(BASEEOSMAC(eossimpletablesspec,0,0,0,0,0,0,0)),1*EOSSIMPLEN5*EOSSIMPLEN4*EOSSIMPLEN3*EOSSIMPLEN2*EOSSIMPLEN1*NUMEOSSSPECQUANTITIESMEM,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&(BASEEOSMAC(eossimpletablepofchi,0,0,0,0,0,0,0)),1*EOSSIMPLEN5*EOSSIMPLEN4*EOSSIMPLEN3*EOSSIMPLEN2*EOSSIMPLEN1*NUMEOSPOFCHIQUANTITIESMEM,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&(BASEEOSMAC(eossimpletabletemp,0,0,0,0,0,0,0)),NUMEOSTEMPQUANTITIESMEM1*EOSSIMPLEN5*EOSSIMPLEN4*EOSSIMPLEN3*EOSSIMPLEN2*EOSSIMPLEN1*NUMEOSTEMPQUANTITIESMEM2,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);

  // different size for extra table
  MPI_Bcast(&(BASEEOSMAC(eossimpletableextra,0,0,0,0,0,0,0)),1*EOSSIMPLEEXTRAN5*EOSSIMPLEEXTRAN4*EOSSIMPLEEXTRAN3*EOSSIMPLEEXTRAN2*EOSSIMPLEEXTRAN1*NUMEOSEXTRAQUANTITIESMEM,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
#if(WHICHDATATYPEGENERAL==4)
  MPI_Bcast(&(BASEEOSMAC(eossimpletableextradegen,0,0,0,0,0,0,0)),NUMEOSDEGENQUANTITIESMEM1*EOSSIMPLEEXTRADEGENN5*EOSSIMPLEEXTRADEGENN4*EOSSIMPLEEXTRADEGENN3*EOSSIMPLEEXTRADEGENN2*EOSSIMPLEEXTRADEGENN1*NUMEOSDEGENQUANTITIESMEM2,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&(BASEEOSMAC(eossimpletableextratemp,0,0,0,0,0,0,0)),NUMEOSTEMPQUANTITIESMEM1*EOSSIMPLEEXTRATEMPN5*EOSSIMPLEEXTRATEMPN4*EOSSIMPLEEXTRATEMPN3*EOSSIMPLEEXTRATEMPN2*EOSSIMPLEEXTRATEMPN1*NUMEOSTEMPQUANTITIESMEM2,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
#endif



  // simple zoom
  MPI_Bcast(&(BASEEOSMAC(eossimplezoomtabledegen,0,0,0,0,0,0,0)),NUMEOSDEGENQUANTITIESMEM1*EOSSIMPLEZOOMDEGENN5*EOSSIMPLEZOOMDEGENN4*EOSSIMPLEZOOMDEGENN3*EOSSIMPLEZOOMDEGENN2*EOSSIMPLEZOOMDEGENN1*NUMEOSDEGENQUANTITIESMEM2,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&(BASEEOSMAC(eossimplezoomtablestandard,0,0,0,0,0,0,0)),1*EOSSIMPLEZOOMN5*EOSSIMPLEZOOMN4*EOSSIMPLEZOOMN3*EOSSIMPLEZOOMN2*EOSSIMPLEZOOMN1*NUMEOSSTANDARDQUANTITIESMEM,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&(BASEEOSMAC(eossimplezoomtableguess,0,0,0,0,0,0,0)),1*EOSSIMPLEN5*EOSSIMPLEN4*EOSSIMPLEN3*EOSSIMPLEN2*EOSSIMPLEN1*NUMEOSGUESSQUANTITIESMEM,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&(BASEEOSMAC(eossimplezoomtablediss,0,0,0,0,0,0,0)),1*EOSSIMPLEN5*EOSSIMPLEN4*EOSSIMPLEN3*EOSSIMPLEN2*EOSSIMPLEN1*NUMEOSDISSQUANTITIESMEM,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&(BASEEOSMAC(eossimplezoomtabledp,0,0,0,0,0,0,0)),1*EOSSIMPLEZOOMN5*EOSSIMPLEZOOMN4*EOSSIMPLEZOOMN3*EOSSIMPLEZOOMN2*EOSSIMPLEZOOMN1*NUMEOSDPQUANTITIESMEM,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&(BASEEOSMAC(eossimplezoomtablesden,0,0,0,0,0,0,0)),1*EOSSIMPLEZOOMN5*EOSSIMPLEZOOMN4*EOSSIMPLEZOOMN3*EOSSIMPLEZOOMN2*EOSSIMPLEZOOMN1*NUMEOSSDENQUANTITIESMEM,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&(BASEEOSMAC(eossimplezoomtablesspec,0,0,0,0,0,0,0)),1*EOSSIMPLEZOOMN5*EOSSIMPLEZOOMN4*EOSSIMPLEZOOMN3*EOSSIMPLEZOOMN2*EOSSIMPLEZOOMN1*NUMEOSSSPECQUANTITIESMEM,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&(BASEEOSMAC(eossimplezoomtablepofchi,0,0,0,0,0,0,0)),1*EOSSIMPLEZOOMN5*EOSSIMPLEZOOMN4*EOSSIMPLEZOOMN3*EOSSIMPLEZOOMN2*EOSSIMPLEZOOMN1*NUMEOSPOFCHIQUANTITIESMEM,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&(BASEEOSMAC(eossimplezoomtabletemp,0,0,0,0,0,0,0)),NUMEOSTEMPQUANTITIESMEM1*EOSSIMPLEZOOMN5*EOSSIMPLEZOOMN4*EOSSIMPLEZOOMN3*EOSSIMPLEZOOMN2*EOSSIMPLEZOOMN1*NUMEOSTEMPQUANTITIESMEM2,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
 
  // different size for extra table
  MPI_Bcast(&(BASEEOSMAC(eossimplezoomtableextra,0,0,0,0,0,0,0)),1*EOSSIMPLEZOOMEXTRAN5*EOSSIMPLEZOOMEXTRAN4*EOSSIMPLEZOOMEXTRAN3*EOSSIMPLEZOOMEXTRAN2*EOSSIMPLEZOOMEXTRAN1*NUMEOSEXTRAQUANTITIESMEM,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
#if(WHICHDATATYPEGENERAL==4)
  MPI_Bcast(&(BASEEOSMAC(eossimplezoomtableextradegen,0,0,0,0,0,0,0)),NUMEOSDEGENQUANTITIESMEM1*EOSSIMPLEZOOMEXTRADEGENN5*EOSSIMPLEZOOMEXTRADEGENN4*EOSSIMPLEZOOMEXTRADEGENN3*EOSSIMPLEZOOMEXTRADEGENN2*EOSSIMPLEZOOMEXTRADEGENN1*NUMEOSDEGENQUANTITIESMEM2,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
  MPI_Bcast(&(BASEEOSMAC(eossimplezoomtableextratemp,0,0,0,0,0,0,0)),NUMEOSTEMPQUANTITIESMEM1*EOSSIMPLEZOOMEXTRATEMPN5*EOSSIMPLEZOOMEXTRATEMPN4*EOSSIMPLEZOOMEXTRATEMPN3*EOSSIMPLEZOOMEXTRATEMPN2*EOSSIMPLEZOOMEXTRATEMPN1*NUMEOSTEMPQUANTITIESMEM2,MPI_FTYPEEOS,MPIid[0],MPI_COMM_GRMHD);
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
// Note that arrays are accessed relative to 0 correspondingn to "coli" and not "whichfun"
static void set_arrays_eostable(int whichdegen, int whichtablevar, int mmm, int lll, int kkk, int jjj, int iii, int incol, FTYPEEOS value)
{


  // overrides:
  // then modify indices since different sub tables have different dimensions.  This will cause repeated-reads to put into same locations in memory.
  if(WHICHDATATYPEGENERAL==4){
    if(
       ((whichtablevar==FULLTABLE || whichtablevar==SIMPLETABLE || whichtablevar==SIMPLEZOOMTABLE) && (incol<FIRSTEXTRAin || incol>LASTEXTRAin)) || // for all non-extras in normal table
       ((whichtablevar==FULLTABLEEXTRA || whichtablevar==SIMPLETABLEEXTRA || whichtablevar==EXTRASIMPLEZOOMTABLE) && (incol<FIRSTEXTRAin || incol>LASTEXTRAin)) || // for temperature in extra table
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
  else if(whichtablevar==FULLTABLE || whichtablevar==FULLTABLEEXTRA){
    if(whichdegen==ISNOTDEGENTABLE){

      if(whichtablevar==FULLTABLE){
        if(incol==PofRHOUin) EOSMAC(eosfulltablestandard,0,mmm,lll,kkk,jjj,iii,PofRHOU-FIRSTEOSSTANDARD)=value;
        if(incol==CS2ofRHOUin) EOSMAC(eosfulltablestandard,0,mmm,lll,kkk,jjj,iii,CS2ofRHOU-FIRSTEOSSTANDARD)=value;

        if(incol==UofRHOPin) EOSMAC(eosfulltableguess,0,mmm,lll,kkk,jjj,iii,UofRHOP-FIRSTEOSGUESS)=value;

        if(incol==UofRHOSin) EOSMAC(eosfulltablediss,0,mmm,lll,kkk,jjj,iii,UofRHOS-FIRSTEOSDISS)=value;

        if(incol==DPDRHOofRHOUin) EOSMAC(eosfulltabledp,0,mmm,lll,kkk,jjj,iii,DPDRHOofRHOU-FIRSTEOSDP)=value;
        if(incol==DPDUofRHOUin) EOSMAC(eosfulltabledp,0,mmm,lll,kkk,jjj,iii,DPDUofRHOU-FIRSTEOSDP)=value;

        if(incol==SofRHOUin) EOSMAC(eosfulltablesden,0,mmm,lll,kkk,jjj,iii,SofRHOU-FIRSTEOSSDEN)=value;
        if(incol==DSDRHOofRHOUin) EOSMAC(eosfulltablesden,0,mmm,lll,kkk,jjj,iii,DSDRHOofRHOU-FIRSTEOSSDEN)=value;
        if(incol==DSDUofRHOUin) EOSMAC(eosfulltablesden,0,mmm,lll,kkk,jjj,iii,DSDUofRHOU-FIRSTEOSSDEN)=value;

        if(incol==SSofRHOCHIin) EOSMAC(eosfulltablesspec,0,mmm,lll,kkk,jjj,iii,SSofRHOCHI-FIRSTEOSSSPEC)=value;
        if(incol==DSSDRHOofRHOCHIin) EOSMAC(eosfulltablesspec,0,mmm,lll,kkk,jjj,iii,DSSDRHOofRHOCHI-FIRSTEOSSSPEC)=value;
        if(incol==DSSDCHIofRHOCHIin) EOSMAC(eosfulltablesspec,0,mmm,lll,kkk,jjj,iii,DSSDCHIofRHOCHI-FIRSTEOSSSPEC)=value;

        if(incol==PofRHOCHIin) EOSMAC(eosfulltablepofchi,0,mmm,lll,kkk,jjj,iii,PofRHOCHI-FIRSTEOSPOFCHI)=value;
        if(incol==IDRHO0DPin) EOSMAC(eosfulltablepofchi,0,mmm,lll,kkk,jjj,iii,IDRHO0DP-FIRSTEOSPOFCHI)=value;
        if(incol==IDCHIDPin) EOSMAC(eosfulltablepofchi,0,mmm,lll,kkk,jjj,iii,IDCHIDP-FIRSTEOSPOFCHI)=value;

        if(incol==TEMPUin) EOSMAC(eosfulltabletemp,UTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP)=value;
        if(incol==TEMPPin) EOSMAC(eosfulltabletemp,PTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP)=value;
        if(incol==TEMPCHIin) EOSMAC(eosfulltabletemp,CHIDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP)=value;
        if(incol==TEMPSin) EOSMAC(eosfulltabletemp,STOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP)=value;

        // note that if WHICHDATATYPEGENERAL==4, below is just not used
        if(incol>=FIRSTEXTRAin && incol<=LASTEXTRAin) EOSMAC(eosfulltableextra,0,mmm,lll,kkk,jjj,iii,incol-FIRSTEXTRAin)=value; // assumes extra's are ordered in sequence
      }
      else{
        // incol is different then since extras start with 4 temperature quantities and then rest of normal extras
        if(incol>=FIRSTEXTRAin && incol<=LASTEXTRAin) EOSMAC(eosfulltableextra,0,mmm,lll,kkk,jjj,iii,incol-FIRSTEXTRAin)=value; // assumes extra's are ordered in sequence
        if(incol==TEMPUin) EOSMAC(eosfulltableextratemp,UTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP)=value;
        if(incol==TEMPPin) EOSMAC(eosfulltableextratemp,PTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP)=value;
        if(incol==TEMPCHIin) EOSMAC(eosfulltableextratemp,CHIDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP)=value;
        if(incol==TEMPSin) EOSMAC(eosfulltableextratemp,STOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP)=value;
      }
    }
    else{
      // degen tables the same but with different name
      if(whichtablevar==FULLTABLE){
        if(incol==UTOTOFFSETin) EOSMAC(eosfulltabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN)=value;
        if(incol==PTOTOFFSETin) EOSMAC(eosfulltabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN)=value;
        if(incol==CHIOFFSETin) EOSMAC(eosfulltabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN)=value;
        if(incol==STOTOFFSETin) EOSMAC(eosfulltabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN)=value;

        if(utotdegencut[whichtablevar]<=DEGENCUTLASTOLDVERSION){
        }
        else{ // utotdegencut[whichtablevar]>DEGENCUTLASTOLDVERSION
          if(incol==UTOTINin) EOSMAC(eosfulltabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN)=value;
          if(incol==PTOTINin) EOSMAC(eosfulltabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN)=value;
          if(incol==CHIINin) EOSMAC(eosfulltabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN)=value;
          if(incol==STOTINin) EOSMAC(eosfulltabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN)=value;

          if(incol==UTOTOUTin) EOSMAC(eosfulltabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN)=value;
          if(incol==PTOTOUTin) EOSMAC(eosfulltabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN)=value;
          if(incol==CHIOUTin) EOSMAC(eosfulltabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN)=value;
          if(incol==STOTOUTin) EOSMAC(eosfulltabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN)=value;
        }
      }
      else{
        if(incol==UTOTOFFSETin) EOSMAC(eosfulltableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN)=value;
        if(incol==PTOTOFFSETin) EOSMAC(eosfulltableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN)=value;
        if(incol==CHIOFFSETin) EOSMAC(eosfulltableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN)=value;
        if(incol==STOTOFFSETin) EOSMAC(eosfulltableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN)=value;

        if(utotdegencut[whichtablevar]<=DEGENCUTLASTOLDVERSION){
        }
        else{ // utotdegencut[whichtablevar]>DEGENCUTLASTOLDVERSION
          if(incol==UTOTINin) EOSMAC(eosfulltableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN)=value;
          if(incol==PTOTINin) EOSMAC(eosfulltableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN)=value;
          if(incol==CHIINin) EOSMAC(eosfulltableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN)=value;
          if(incol==STOTINin) EOSMAC(eosfulltableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN)=value;

          if(incol==UTOTOUTin) EOSMAC(eosfulltableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN)=value;
          if(incol==PTOTOUTin) EOSMAC(eosfulltableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN)=value;
          if(incol==CHIOUTin) EOSMAC(eosfulltableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN)=value;
          if(incol==STOTOUTin) EOSMAC(eosfulltableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN)=value;
        }
      }
    }
  }
#endif
#if(ALLOWSIMPLETABLE==1)
  else if(whichtablevar==SIMPLETABLE || whichtablevar==SIMPLETABLEEXTRA){
    if(whichdegen==ISNOTDEGENTABLE){

      if(whichtablevar==SIMPLETABLE){
        if(incol==PofRHOUin) EOSMAC(eossimpletablestandard,0,mmm,lll,kkk,jjj,iii,PofRHOU-FIRSTEOSSTANDARD)=value;
        if(incol==CS2ofRHOUin) EOSMAC(eossimpletablestandard,0,mmm,lll,kkk,jjj,iii,CS2ofRHOU-FIRSTEOSSTANDARD)=value;

        if(incol==UofRHOPin) EOSMAC(eossimpletableguess,0,mmm,lll,kkk,jjj,iii,UofRHOP-FIRSTEOSGUESS)=value;

        if(incol==UofRHOSin) EOSMAC(eossimpletablediss,0,mmm,lll,kkk,jjj,iii,UofRHOS-FIRSTEOSDISS)=value;

        if(incol==DPDRHOofRHOUin) EOSMAC(eossimpletabledp,0,mmm,lll,kkk,jjj,iii,DPDRHOofRHOU-FIRSTEOSDP)=value;
        if(incol==DPDUofRHOUin) EOSMAC(eossimpletabledp,0,mmm,lll,kkk,jjj,iii,DPDUofRHOU-FIRSTEOSDP)=value;

        if(incol==SofRHOUin) EOSMAC(eossimpletablesden,0,mmm,lll,kkk,jjj,iii,SofRHOU-FIRSTEOSSDEN)=value;
        if(incol==DSDRHOofRHOUin) EOSMAC(eossimpletablesden,0,mmm,lll,kkk,jjj,iii,DSDRHOofRHOU-FIRSTEOSSDEN)=value;
        if(incol==DSDUofRHOUin) EOSMAC(eossimpletablesden,0,mmm,lll,kkk,jjj,iii,DSDUofRHOU-FIRSTEOSSDEN)=value;

        if(incol==SSofRHOCHIin) EOSMAC(eossimpletablesspec,0,mmm,lll,kkk,jjj,iii,SSofRHOCHI-FIRSTEOSSSPEC)=value;
        if(incol==DSSDRHOofRHOCHIin) EOSMAC(eossimpletablesspec,0,mmm,lll,kkk,jjj,iii,DSSDRHOofRHOCHI-FIRSTEOSSSPEC)=value;
        if(incol==DSSDCHIofRHOCHIin) EOSMAC(eossimpletablesspec,0,mmm,lll,kkk,jjj,iii,DSSDCHIofRHOCHI-FIRSTEOSSSPEC)=value;

        if(incol==PofRHOCHIin) EOSMAC(eossimpletablepofchi,0,mmm,lll,kkk,jjj,iii,PofRHOCHI-FIRSTEOSPOFCHI)=value;
        if(incol==IDRHO0DPin) EOSMAC(eossimpletablepofchi,0,mmm,lll,kkk,jjj,iii,IDRHO0DP-FIRSTEOSPOFCHI)=value;
        if(incol==IDCHIDPin) EOSMAC(eossimpletablepofchi,0,mmm,lll,kkk,jjj,iii,IDCHIDP-FIRSTEOSPOFCHI)=value;

        if(incol==TEMPUin) EOSMAC(eossimpletabletemp,UTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP)=value;
        if(incol==TEMPPin) EOSMAC(eossimpletabletemp,PTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP)=value;
        if(incol==TEMPCHIin) EOSMAC(eossimpletabletemp,CHIDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP)=value;
        if(incol==TEMPSin) EOSMAC(eossimpletabletemp,STOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP)=value;

        // note that if WHICHDATATYPEGENERAL==4, below is just not used
        if(incol>=FIRSTEXTRAin && incol<=LASTEXTRAin) EOSMAC(eossimpletableextra,0,mmm,lll,kkk,jjj,iii,incol-FIRSTEXTRAin)=value; // assumes extra's are ordered in sequence
      }
      else{
        // incol is different then since extras start with 4 temperature quantities and then rest of normal extras
        if(incol>=FIRSTEXTRAin && incol<=LASTEXTRAin) EOSMAC(eossimpletableextra,0,mmm,lll,kkk,jjj,iii,incol-FIRSTEXTRAin)=value; // assumes extra's are ordered in sequence
        if(incol==TEMPUin) EOSMAC(eossimpletableextratemp,UTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP)=value;
        if(incol==TEMPPin) EOSMAC(eossimpletableextratemp,PTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP)=value;
        if(incol==TEMPCHIin) EOSMAC(eossimpletableextratemp,CHIDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP)=value;
        if(incol==TEMPSin) EOSMAC(eossimpletableextratemp,STOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP)=value;
      }
    }
    else{
      // degen tables the same but with different name
      if(whichtablevar==SIMPLETABLE){
        if(incol==UTOTOFFSETin) EOSMAC(eossimpletabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN)=value;
        if(incol==PTOTOFFSETin) EOSMAC(eossimpletabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN)=value;
        if(incol==CHIOFFSETin) EOSMAC(eossimpletabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN)=value;
        if(incol==STOTOFFSETin) EOSMAC(eossimpletabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN)=value;

        if(utotdegencut[whichtablevar]<=DEGENCUTLASTOLDVERSION){
        }
        else{ // utotdegencut[whichtablevar]>DEGENCUTLASTOLDVERSION
          if(incol==UTOTINin) EOSMAC(eossimpletabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN)=value;
          if(incol==PTOTINin) EOSMAC(eossimpletabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN)=value;
          if(incol==CHIINin) EOSMAC(eossimpletabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN)=value;
          if(incol==STOTINin) EOSMAC(eossimpletabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN)=value;

          if(incol==UTOTOUTin) EOSMAC(eossimpletabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN)=value;
          if(incol==PTOTOUTin) EOSMAC(eossimpletabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN)=value;
          if(incol==CHIOUTin) EOSMAC(eossimpletabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN)=value;
          if(incol==STOTOUTin) EOSMAC(eossimpletabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN)=value;
        }
      }
      else{
        if(incol==UTOTOFFSETin) EOSMAC(eossimpletableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN)=value;
        if(incol==PTOTOFFSETin) EOSMAC(eossimpletableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN)=value;
        if(incol==CHIOFFSETin) EOSMAC(eossimpletableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN)=value;
        if(incol==STOTOFFSETin) EOSMAC(eossimpletableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN)=value;

        if(utotdegencut[whichtablevar]<=DEGENCUTLASTOLDVERSION){
        }
        else{ // utotdegencut[whichtablevar]>DEGENCUTLASTOLDVERSION
          if(incol==UTOTINin) EOSMAC(eossimpletableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN)=value;
          if(incol==PTOTINin) EOSMAC(eossimpletableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN)=value;
          if(incol==CHIINin) EOSMAC(eossimpletableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN)=value;
          if(incol==STOTINin) EOSMAC(eossimpletableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN)=value;

          if(incol==UTOTOUTin) EOSMAC(eossimpletableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN)=value;
          if(incol==PTOTOUTin) EOSMAC(eossimpletableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN)=value;
          if(incol==CHIOUTin) EOSMAC(eossimpletableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN)=value;
          if(incol==STOTOUTin) EOSMAC(eossimpletableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN)=value;
        }
      }
    }
  }
#endif
#if(ALLOWSIMPLEZOOMTABLE==1)
  else if(whichtablevar==SIMPLEZOOMTABLE || whichtablevar==EXTRASIMPLEZOOMTABLE){
    if(whichdegen==ISNOTDEGENTABLE){

      if(whichtablevar==SIMPLEZOOMTABLE){
        if(incol==PofRHOUin) EOSMAC(eossimplezoomtablestandard,0,mmm,lll,kkk,jjj,iii,PofRHOU-FIRSTEOSSTANDARD)=value;
        if(incol==CS2ofRHOUin) EOSMAC(eossimplezoomtablestandard,0,mmm,lll,kkk,jjj,iii,CS2ofRHOU-FIRSTEOSSTANDARD)=value;

        if(incol==UofRHOPin) EOSMAC(eossimplezoomtableguess,0,mmm,lll,kkk,jjj,iii,UofRHOP-FIRSTEOSGUESS)=value;

        if(incol==UofRHOSin) EOSMAC(eossimplezoomtablediss,0,mmm,lll,kkk,jjj,iii,UofRHOS-FIRSTEOSDISS)=value;

        if(incol==DPDRHOofRHOUin) EOSMAC(eossimplezoomtabledp,0,mmm,lll,kkk,jjj,iii,DPDRHOofRHOU-FIRSTEOSDP)=value;
        if(incol==DPDUofRHOUin) EOSMAC(eossimplezoomtabledp,0,mmm,lll,kkk,jjj,iii,DPDUofRHOU-FIRSTEOSDP)=value;

        if(incol==SofRHOUin) EOSMAC(eossimplezoomtablesden,0,mmm,lll,kkk,jjj,iii,SofRHOU-FIRSTEOSSDEN)=value;
        if(incol==DSDRHOofRHOUin) EOSMAC(eossimplezoomtablesden,0,mmm,lll,kkk,jjj,iii,DSDRHOofRHOU-FIRSTEOSSDEN)=value;
        if(incol==DSDUofRHOUin) EOSMAC(eossimplezoomtablesden,0,mmm,lll,kkk,jjj,iii,DSDUofRHOU-FIRSTEOSSDEN)=value;

        if(incol==SSofRHOCHIin) EOSMAC(eossimplezoomtablesspec,0,mmm,lll,kkk,jjj,iii,SSofRHOCHI-FIRSTEOSSSPEC)=value;
        if(incol==DSSDRHOofRHOCHIin) EOSMAC(eossimplezoomtablesspec,0,mmm,lll,kkk,jjj,iii,DSSDRHOofRHOCHI-FIRSTEOSSSPEC)=value;
        if(incol==DSSDCHIofRHOCHIin) EOSMAC(eossimplezoomtablesspec,0,mmm,lll,kkk,jjj,iii,DSSDCHIofRHOCHI-FIRSTEOSSSPEC)=value;

        if(incol==PofRHOCHIin) EOSMAC(eossimplezoomtablepofchi,0,mmm,lll,kkk,jjj,iii,PofRHOCHI-FIRSTEOSPOFCHI)=value;
        if(incol==IDRHO0DPin) EOSMAC(eossimplezoomtablepofchi,0,mmm,lll,kkk,jjj,iii,IDRHO0DP-FIRSTEOSPOFCHI)=value;
        if(incol==IDCHIDPin) EOSMAC(eossimplezoomtablepofchi,0,mmm,lll,kkk,jjj,iii,IDCHIDP-FIRSTEOSPOFCHI)=value;

        if(incol==TEMPUin) EOSMAC(eossimplezoomtabletemp,UTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP)=value;
        if(incol==TEMPPin) EOSMAC(eossimplezoomtabletemp,PTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP)=value;
        if(incol==TEMPCHIin) EOSMAC(eossimplezoomtabletemp,CHIDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP)=value;
        if(incol==TEMPSin) EOSMAC(eossimplezoomtabletemp,STOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP)=value;

        // note that if WHICHDATATYPEGENERAL==4, below is just not used
        if(incol>=FIRSTEXTRAin && incol<=LASTEXTRAin) EOSMAC(eossimplezoomtableextra,0,mmm,lll,kkk,jjj,iii,incol-FIRSTEXTRAin)=value; // assumes extra's are ordered in sequence
      }
      else{
        // incol is different then since extras start with 4 temperature quantities and then rest of normal extras
        if(incol>=FIRSTEXTRAin && incol<=LASTEXTRAin) EOSMAC(eossimplezoomtableextra,0,mmm,lll,kkk,jjj,iii,incol-FIRSTEXTRAin)=value; // assumes extra's are ordered in sequence
        if(incol==TEMPUin) EOSMAC(eossimplezoomtableextratemp,UTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP)=value;
        if(incol==TEMPPin) EOSMAC(eossimplezoomtableextratemp,PTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP)=value;
        if(incol==TEMPCHIin) EOSMAC(eossimplezoomtableextratemp,CHIDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP)=value;
        if(incol==TEMPSin) EOSMAC(eossimplezoomtableextratemp,STOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP)=value;
      }
    }
    else{
      // degen tables the same but with different name
      if(whichtablevar==SIMPLEZOOMTABLE){
        if(incol==UTOTOFFSETin) EOSMAC(eossimplezoomtabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN)=value;
        if(incol==PTOTOFFSETin) EOSMAC(eossimplezoomtabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN)=value;
        if(incol==CHIOFFSETin) EOSMAC(eossimplezoomtabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN)=value;
        if(incol==STOTOFFSETin) EOSMAC(eossimplezoomtabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN)=value;

        if(utotdegencut[whichtablevar]<=DEGENCUTLASTOLDVERSION){
        }
        else{ // utotdegencut[whichtablevar]>DEGENCUTLASTOLDVERSION
          if(incol==UTOTINin) EOSMAC(eossimplezoomtabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN)=value;
          if(incol==PTOTINin) EOSMAC(eossimplezoomtabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN)=value;
          if(incol==CHIINin) EOSMAC(eossimplezoomtabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN)=value;
          if(incol==STOTINin) EOSMAC(eossimplezoomtabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN)=value;

          if(incol==UTOTOUTin) EOSMAC(eossimplezoomtabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN)=value;
          if(incol==PTOTOUTin) EOSMAC(eossimplezoomtabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN)=value;
          if(incol==CHIOUTin) EOSMAC(eossimplezoomtabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN)=value;
          if(incol==STOTOUTin) EOSMAC(eossimplezoomtabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN)=value;
        }
      }
      else{
        if(incol==UTOTOFFSETin) EOSMAC(eossimplezoomtableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN)=value;
        if(incol==PTOTOFFSETin) EOSMAC(eossimplezoomtableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN)=value;
        if(incol==CHIOFFSETin) EOSMAC(eossimplezoomtableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN)=value;
        if(incol==STOTOFFSETin) EOSMAC(eossimplezoomtableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN)=value;

        if(utotdegencut[whichtablevar]<=DEGENCUTLASTOLDVERSION){
        }
        else{ // utotdegencut[whichtablevar]>DEGENCUTLASTOLDVERSION
          if(incol==UTOTINin) EOSMAC(eossimplezoomtableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN)=value;
          if(incol==PTOTINin) EOSMAC(eossimplezoomtableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN)=value;
          if(incol==CHIINin) EOSMAC(eossimplezoomtableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN)=value;
          if(incol==STOTINin) EOSMAC(eossimplezoomtableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN)=value;

          if(incol==UTOTOUTin) EOSMAC(eossimplezoomtableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN)=value;
          if(incol==PTOTOUTin) EOSMAC(eossimplezoomtableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN)=value;
          if(incol==CHIOUTin) EOSMAC(eossimplezoomtableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN)=value;
          if(incol==STOTOUTin) EOSMAC(eossimplezoomtableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN)=value;
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
// and last function parameter goes from "FTYPEEOS value" -> "FTYPEEOS *value"
//
// input jjj should be 0 if inputting degen table data
// only used by read_setup_eostable()
// no longer refer to macro label "inextra" since put into canonical order and position before calling set_ or get_ functions
static void get_arrays_eostable(int whichdegen, int whichtablevar, int mmm, int lll, int kkk, int jjj, int iii, int incol, FTYPEEOS *value)
{







  // overrides:
  // then modify indices since different sub tables have different dimensions.  This will cause repeated-reads to put into same locations in memory.
  if(WHICHDATATYPEGENERAL==4){
    if(
       ((whichtablevar==FULLTABLE || whichtablevar==SIMPLETABLE || whichtablevar==SIMPLEZOOMTABLE) && (incol<FIRSTEXTRAin || incol>LASTEXTRAin)) || // for all non-extras in normal table
       ((whichtablevar==FULLTABLEEXTRA || whichtablevar==SIMPLETABLEEXTRA || whichtablevar==EXTRASIMPLEZOOMTABLE) && (incol<FIRSTEXTRAin || incol>LASTEXTRAin)) || // for temperature in extra table
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
  else if(whichtablevar==FULLTABLE || whichtablevar==FULLTABLEEXTRA){
    if(whichdegen==ISNOTDEGENTABLE){

      if(whichtablevar==FULLTABLE){
        if(incol==PofRHOUin) *value=EOSMAC(eosfulltablestandard,0,mmm,lll,kkk,jjj,iii,PofRHOU-FIRSTEOSSTANDARD);
        if(incol==CS2ofRHOUin) *value=EOSMAC(eosfulltablestandard,0,mmm,lll,kkk,jjj,iii,CS2ofRHOU-FIRSTEOSSTANDARD);

        if(incol==UofRHOPin) *value=EOSMAC(eosfulltableguess,0,mmm,lll,kkk,jjj,iii,UofRHOP-FIRSTEOSGUESS);

        if(incol==UofRHOSin) *value=EOSMAC(eosfulltablediss,0,mmm,lll,kkk,jjj,iii,UofRHOS-FIRSTEOSDISS);

        if(incol==DPDRHOofRHOUin) *value=EOSMAC(eosfulltabledp,0,mmm,lll,kkk,jjj,iii,DPDRHOofRHOU-FIRSTEOSDP);
        if(incol==DPDUofRHOUin) *value=EOSMAC(eosfulltabledp,0,mmm,lll,kkk,jjj,iii,DPDUofRHOU-FIRSTEOSDP);

        if(incol==SofRHOUin) *value=EOSMAC(eosfulltablesden,0,mmm,lll,kkk,jjj,iii,SofRHOU-FIRSTEOSSDEN);
        if(incol==DSDRHOofRHOUin) *value=EOSMAC(eosfulltablesden,0,mmm,lll,kkk,jjj,iii,DSDRHOofRHOU-FIRSTEOSSDEN);
        if(incol==DSDUofRHOUin) *value=EOSMAC(eosfulltablesden,0,mmm,lll,kkk,jjj,iii,DSDUofRHOU-FIRSTEOSSDEN);

        if(incol==SSofRHOCHIin) *value=EOSMAC(eosfulltablesspec,0,mmm,lll,kkk,jjj,iii,SSofRHOCHI-FIRSTEOSSSPEC);
        if(incol==DSSDRHOofRHOCHIin) *value=EOSMAC(eosfulltablesspec,0,mmm,lll,kkk,jjj,iii,DSSDRHOofRHOCHI-FIRSTEOSSSPEC);
        if(incol==DSSDCHIofRHOCHIin) *value=EOSMAC(eosfulltablesspec,0,mmm,lll,kkk,jjj,iii,DSSDCHIofRHOCHI-FIRSTEOSSSPEC);

        if(incol==PofRHOCHIin) *value=EOSMAC(eosfulltablepofchi,0,mmm,lll,kkk,jjj,iii,PofRHOCHI-FIRSTEOSPOFCHI);
        if(incol==IDRHO0DPin) *value=EOSMAC(eosfulltablepofchi,0,mmm,lll,kkk,jjj,iii,IDRHO0DP-FIRSTEOSPOFCHI);
        if(incol==IDCHIDPin) *value=EOSMAC(eosfulltablepofchi,0,mmm,lll,kkk,jjj,iii,IDCHIDP-FIRSTEOSPOFCHI);

        if(incol==TEMPUin) *value=EOSMAC(eosfulltabletemp,UTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP);
        if(incol==TEMPPin) *value=EOSMAC(eosfulltabletemp,PTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP);
        if(incol==TEMPCHIin) *value=EOSMAC(eosfulltabletemp,CHIDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP);
        if(incol==TEMPSin) *value=EOSMAC(eosfulltabletemp,STOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP);

        // note that if WHICHDATATYPEGENERAL==4, below is just not used
        if(incol>=FIRSTEXTRAin && incol<=LASTEXTRAin) *value=EOSMAC(eosfulltableextra,0,mmm,lll,kkk,jjj,iii,incol-FIRSTEXTRAin); // assumes extra's are ordered in sequence
      }
      else{
        // incol is different then since extras start with 4 temperature quantities and then rest of normal extras
        if(incol>=FIRSTEXTRAin && incol<=LASTEXTRAin) *value=EOSMAC(eosfulltableextra,0,mmm,lll,kkk,jjj,iii,incol-FIRSTEXTRAin); // assumes extra's are ordered in sequence
        if(incol==TEMPUin) *value=EOSMAC(eosfulltableextratemp,UTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP);
        if(incol==TEMPPin) *value=EOSMAC(eosfulltableextratemp,PTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP);
        if(incol==TEMPCHIin) *value=EOSMAC(eosfulltableextratemp,CHIDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP);
        if(incol==TEMPSin) *value=EOSMAC(eosfulltableextratemp,STOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP);
      }
    }
    else{
      // degen tables the same but with different name
      if(whichtablevar==FULLTABLE){
        if(incol==UTOTOFFSETin) *value=EOSMAC(eosfulltabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN);
        if(incol==PTOTOFFSETin) *value=EOSMAC(eosfulltabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN);
        if(incol==CHIOFFSETin) *value=EOSMAC(eosfulltabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN);
        if(incol==STOTOFFSETin) *value=EOSMAC(eosfulltabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN);

        if(utotdegencut[whichtablevar]<=DEGENCUTLASTOLDVERSION){
        }
        else{ // utotdegencut[whichtablevar]>DEGENCUTLASTOLDVERSION
          if(incol==UTOTINin) *value=EOSMAC(eosfulltabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN);
          if(incol==PTOTINin) *value=EOSMAC(eosfulltabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN);
          if(incol==CHIINin) *value=EOSMAC(eosfulltabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN);
          if(incol==STOTINin) *value=EOSMAC(eosfulltabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN);

          if(incol==UTOTOUTin) *value=EOSMAC(eosfulltabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN);
          if(incol==PTOTOUTin) *value=EOSMAC(eosfulltabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN);
          if(incol==CHIOUTin) *value=EOSMAC(eosfulltabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN);
          if(incol==STOTOUTin) *value=EOSMAC(eosfulltabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN);
        }
      }
      else{
        if(incol==UTOTOFFSETin) *value=EOSMAC(eosfulltableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN);
        if(incol==PTOTOFFSETin) *value=EOSMAC(eosfulltableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN);
        if(incol==CHIOFFSETin) *value=EOSMAC(eosfulltableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN);
        if(incol==STOTOFFSETin) *value=EOSMAC(eosfulltableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN);

        if(utotdegencut[whichtablevar]<=DEGENCUTLASTOLDVERSION){
        }
        else{ // utotdegencut[whichtablevar]>DEGENCUTLASTOLDVERSION
          if(incol==UTOTINin) *value=EOSMAC(eosfulltableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN);
          if(incol==PTOTINin) *value=EOSMAC(eosfulltableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN);
          if(incol==CHIINin) *value=EOSMAC(eosfulltableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN);
          if(incol==STOTINin) *value=EOSMAC(eosfulltableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN);

          if(incol==UTOTOUTin) *value=EOSMAC(eosfulltableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN);
          if(incol==PTOTOUTin) *value=EOSMAC(eosfulltableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN);
          if(incol==CHIOUTin) *value=EOSMAC(eosfulltableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN);
          if(incol==STOTOUTin) *value=EOSMAC(eosfulltableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN);
        }
      }
    }
  }
#endif
#if(ALLOWSIMPLETABLE==1)
  else if(whichtablevar==SIMPLETABLE || whichtablevar==SIMPLETABLEEXTRA){
    if(whichdegen==ISNOTDEGENTABLE){

      if(whichtablevar==SIMPLETABLE){
        if(incol==PofRHOUin) *value=EOSMAC(eossimpletablestandard,0,mmm,lll,kkk,jjj,iii,PofRHOU-FIRSTEOSSTANDARD);
        if(incol==CS2ofRHOUin) *value=EOSMAC(eossimpletablestandard,0,mmm,lll,kkk,jjj,iii,CS2ofRHOU-FIRSTEOSSTANDARD);

        if(incol==UofRHOPin) *value=EOSMAC(eossimpletableguess,0,mmm,lll,kkk,jjj,iii,UofRHOP-FIRSTEOSGUESS);

        if(incol==UofRHOSin) *value=EOSMAC(eossimpletablediss,0,mmm,lll,kkk,jjj,iii,UofRHOS-FIRSTEOSDISS);

        if(incol==DPDRHOofRHOUin) *value=EOSMAC(eossimpletabledp,0,mmm,lll,kkk,jjj,iii,DPDRHOofRHOU-FIRSTEOSDP);
        if(incol==DPDUofRHOUin) *value=EOSMAC(eossimpletabledp,0,mmm,lll,kkk,jjj,iii,DPDUofRHOU-FIRSTEOSDP);

        if(incol==SofRHOUin) *value=EOSMAC(eossimpletablesden,0,mmm,lll,kkk,jjj,iii,SofRHOU-FIRSTEOSSDEN);
        if(incol==DSDRHOofRHOUin) *value=EOSMAC(eossimpletablesden,0,mmm,lll,kkk,jjj,iii,DSDRHOofRHOU-FIRSTEOSSDEN);
        if(incol==DSDUofRHOUin) *value=EOSMAC(eossimpletablesden,0,mmm,lll,kkk,jjj,iii,DSDUofRHOU-FIRSTEOSSDEN);

        if(incol==SSofRHOCHIin) *value=EOSMAC(eossimpletablesspec,0,mmm,lll,kkk,jjj,iii,SSofRHOCHI-FIRSTEOSSSPEC);
        if(incol==DSSDRHOofRHOCHIin) *value=EOSMAC(eossimpletablesspec,0,mmm,lll,kkk,jjj,iii,DSSDRHOofRHOCHI-FIRSTEOSSSPEC);
        if(incol==DSSDCHIofRHOCHIin) *value=EOSMAC(eossimpletablesspec,0,mmm,lll,kkk,jjj,iii,DSSDCHIofRHOCHI-FIRSTEOSSSPEC);

        if(incol==PofRHOCHIin) *value=EOSMAC(eossimpletablepofchi,0,mmm,lll,kkk,jjj,iii,PofRHOCHI-FIRSTEOSPOFCHI);
        if(incol==IDRHO0DPin) *value=EOSMAC(eossimpletablepofchi,0,mmm,lll,kkk,jjj,iii,IDRHO0DP-FIRSTEOSPOFCHI);
        if(incol==IDCHIDPin) *value=EOSMAC(eossimpletablepofchi,0,mmm,lll,kkk,jjj,iii,IDCHIDP-FIRSTEOSPOFCHI);

        if(incol==TEMPUin) *value=EOSMAC(eossimpletabletemp,UTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP);
        if(incol==TEMPPin) *value=EOSMAC(eossimpletabletemp,PTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP);
        if(incol==TEMPCHIin) *value=EOSMAC(eossimpletabletemp,CHIDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP);
        if(incol==TEMPSin) *value=EOSMAC(eossimpletabletemp,STOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP);

        // note that if WHICHDATATYPEGENERAL==4, below is just not used
        if(incol>=FIRSTEXTRAin && incol<=LASTEXTRAin) *value=EOSMAC(eossimpletableextra,0,mmm,lll,kkk,jjj,iii,incol-FIRSTEXTRAin); // assumes extra's are ordered in sequence
      }
      else{
        // incol is different then since extras start with 4 temperature quantities and then rest of normal extras
        if(incol>=FIRSTEXTRAin && incol<=LASTEXTRAin) *value=EOSMAC(eossimpletableextra,0,mmm,lll,kkk,jjj,iii,incol-FIRSTEXTRAin); // assumes extra's are ordered in sequence
        if(incol==TEMPUin) *value=EOSMAC(eossimpletableextratemp,UTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP);
        if(incol==TEMPPin) *value=EOSMAC(eossimpletableextratemp,PTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP);
        if(incol==TEMPCHIin) *value=EOSMAC(eossimpletableextratemp,CHIDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP);
        if(incol==TEMPSin) *value=EOSMAC(eossimpletableextratemp,STOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP);
      }
    }
    else{
      // degen tables the same but with different name
      if(whichtablevar==SIMPLETABLE){
        if(incol==UTOTOFFSETin) *value=EOSMAC(eossimpletabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN);
        if(incol==PTOTOFFSETin) *value=EOSMAC(eossimpletabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN);
        if(incol==CHIOFFSETin) *value=EOSMAC(eossimpletabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN);
        if(incol==STOTOFFSETin) *value=EOSMAC(eossimpletabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN);

        if(utotdegencut[whichtablevar]<=DEGENCUTLASTOLDVERSION){
        }
        else{ // utotdegencut[whichtablevar]>DEGENCUTLASTOLDVERSION
          if(incol==UTOTINin) *value=EOSMAC(eossimpletabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN);
          if(incol==PTOTINin) *value=EOSMAC(eossimpletabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN);
          if(incol==CHIINin) *value=EOSMAC(eossimpletabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN);
          if(incol==STOTINin) *value=EOSMAC(eossimpletabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN);

          if(incol==UTOTOUTin) *value=EOSMAC(eossimpletabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN);
          if(incol==PTOTOUTin) *value=EOSMAC(eossimpletabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN);
          if(incol==CHIOUTin) *value=EOSMAC(eossimpletabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN);
          if(incol==STOTOUTin) *value=EOSMAC(eossimpletabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN);
        }
      }
      else{
        if(incol==UTOTOFFSETin) *value=EOSMAC(eossimpletableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN);
        if(incol==PTOTOFFSETin) *value=EOSMAC(eossimpletableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN);
        if(incol==CHIOFFSETin) *value=EOSMAC(eossimpletableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN);
        if(incol==STOTOFFSETin) *value=EOSMAC(eossimpletableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN);

        if(utotdegencut[whichtablevar]<=DEGENCUTLASTOLDVERSION){
        }
        else{ // utotdegencut[whichtablevar]>DEGENCUTLASTOLDVERSION
          if(incol==UTOTINin) *value=EOSMAC(eossimpletableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN);
          if(incol==PTOTINin) *value=EOSMAC(eossimpletableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN);
          if(incol==CHIINin) *value=EOSMAC(eossimpletableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN);
          if(incol==STOTINin) *value=EOSMAC(eossimpletableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN);

          if(incol==UTOTOUTin) *value=EOSMAC(eossimpletableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN);
          if(incol==PTOTOUTin) *value=EOSMAC(eossimpletableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN);
          if(incol==CHIOUTin) *value=EOSMAC(eossimpletableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN);
          if(incol==STOTOUTin) *value=EOSMAC(eossimpletableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN);
        }
      }
    }
  }
#endif
#if(ALLOWSIMPLEZOOMTABLE==1)
  else if(whichtablevar==SIMPLEZOOMTABLE || whichtablevar==EXTRASIMPLEZOOMTABLE){
    if(whichdegen==ISNOTDEGENTABLE){

      if(whichtablevar==SIMPLEZOOMTABLE){
        if(incol==PofRHOUin) *value=EOSMAC(eossimplezoomtablestandard,0,mmm,lll,kkk,jjj,iii,PofRHOU-FIRSTEOSSTANDARD);
        if(incol==CS2ofRHOUin) *value=EOSMAC(eossimplezoomtablestandard,0,mmm,lll,kkk,jjj,iii,CS2ofRHOU-FIRSTEOSSTANDARD);

        if(incol==UofRHOPin) *value=EOSMAC(eossimplezoomtableguess,0,mmm,lll,kkk,jjj,iii,UofRHOP-FIRSTEOSGUESS);

        if(incol==UofRHOSin) *value=EOSMAC(eossimplezoomtablediss,0,mmm,lll,kkk,jjj,iii,UofRHOS-FIRSTEOSDISS);

        if(incol==DPDRHOofRHOUin) *value=EOSMAC(eossimplezoomtabledp,0,mmm,lll,kkk,jjj,iii,DPDRHOofRHOU-FIRSTEOSDP);
        if(incol==DPDUofRHOUin) *value=EOSMAC(eossimplezoomtabledp,0,mmm,lll,kkk,jjj,iii,DPDUofRHOU-FIRSTEOSDP);

        if(incol==SofRHOUin) *value=EOSMAC(eossimplezoomtablesden,0,mmm,lll,kkk,jjj,iii,SofRHOU-FIRSTEOSSDEN);
        if(incol==DSDRHOofRHOUin) *value=EOSMAC(eossimplezoomtablesden,0,mmm,lll,kkk,jjj,iii,DSDRHOofRHOU-FIRSTEOSSDEN);
        if(incol==DSDUofRHOUin) *value=EOSMAC(eossimplezoomtablesden,0,mmm,lll,kkk,jjj,iii,DSDUofRHOU-FIRSTEOSSDEN);

        if(incol==SSofRHOCHIin) *value=EOSMAC(eossimplezoomtablesspec,0,mmm,lll,kkk,jjj,iii,SSofRHOCHI-FIRSTEOSSSPEC);
        if(incol==DSSDRHOofRHOCHIin) *value=EOSMAC(eossimplezoomtablesspec,0,mmm,lll,kkk,jjj,iii,DSSDRHOofRHOCHI-FIRSTEOSSSPEC);
        if(incol==DSSDCHIofRHOCHIin) *value=EOSMAC(eossimplezoomtablesspec,0,mmm,lll,kkk,jjj,iii,DSSDCHIofRHOCHI-FIRSTEOSSSPEC);

        if(incol==PofRHOCHIin) *value=EOSMAC(eossimplezoomtablepofchi,0,mmm,lll,kkk,jjj,iii,PofRHOCHI-FIRSTEOSPOFCHI);
        if(incol==IDRHO0DPin) *value=EOSMAC(eossimplezoomtablepofchi,0,mmm,lll,kkk,jjj,iii,IDRHO0DP-FIRSTEOSPOFCHI);
        if(incol==IDCHIDPin) *value=EOSMAC(eossimplezoomtablepofchi,0,mmm,lll,kkk,jjj,iii,IDCHIDP-FIRSTEOSPOFCHI);

        if(incol==TEMPUin) *value=EOSMAC(eossimplezoomtabletemp,UTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP);
        if(incol==TEMPPin) *value=EOSMAC(eossimplezoomtabletemp,PTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP);
        if(incol==TEMPCHIin) *value=EOSMAC(eossimplezoomtabletemp,CHIDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP);
        if(incol==TEMPSin) *value=EOSMAC(eossimplezoomtabletemp,STOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP);

        // note that if WHICHDATATYPEGENERAL==4, below is just not used
        if(incol>=FIRSTEXTRAin && incol<=LASTEXTRAin) *value=EOSMAC(eossimplezoomtableextra,0,mmm,lll,kkk,jjj,iii,incol-FIRSTEXTRAin); // assumes extra's are ordered in sequence
      }
      else{
        // incol is different then since extras start with 4 temperature quantities and then rest of normal extras
        if(incol>=FIRSTEXTRAin && incol<=LASTEXTRAin) *value=EOSMAC(eossimplezoomtableextra,0,mmm,lll,kkk,jjj,iii,incol-FIRSTEXTRAin); // assumes extra's are ordered in sequence
        if(incol==TEMPUin) *value=EOSMAC(eossimplezoomtableextratemp,UTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP);
        if(incol==TEMPPin) *value=EOSMAC(eossimplezoomtableextratemp,PTOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP);
        if(incol==TEMPCHIin) *value=EOSMAC(eossimplezoomtableextratemp,CHIDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP);
        if(incol==TEMPSin) *value=EOSMAC(eossimplezoomtableextratemp,STOTDIFF,mmm,lll,kkk,jjj,iii,TEMPGEN-FIRSTEOSTEMP);
      }
    }
    else{
      // degen tables the same but with different name
      if(whichtablevar==SIMPLEZOOMTABLE){
        if(incol==UTOTOFFSETin) *value=EOSMAC(eossimplezoomtabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN);
        if(incol==PTOTOFFSETin) *value=EOSMAC(eossimplezoomtabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN);
        if(incol==CHIOFFSETin) *value=EOSMAC(eossimplezoomtabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN);
        if(incol==STOTOFFSETin) *value=EOSMAC(eossimplezoomtabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN);

        if(utotdegencut[whichtablevar]<=DEGENCUTLASTOLDVERSION){
        }
        else{ // utotdegencut[whichtablevar]>DEGENCUTLASTOLDVERSION
          if(incol==UTOTINin) *value=EOSMAC(eossimplezoomtabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN);
          if(incol==PTOTINin) *value=EOSMAC(eossimplezoomtabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN);
          if(incol==CHIINin) *value=EOSMAC(eossimplezoomtabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN);
          if(incol==STOTINin) *value=EOSMAC(eossimplezoomtabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN);

          if(incol==UTOTOUTin) *value=EOSMAC(eossimplezoomtabledegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN);
          if(incol==PTOTOUTin) *value=EOSMAC(eossimplezoomtabledegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN);
          if(incol==CHIOUTin) *value=EOSMAC(eossimplezoomtabledegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN);
          if(incol==STOTOUTin) *value=EOSMAC(eossimplezoomtabledegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN);
        }
      }
      else{
        if(incol==UTOTOFFSETin) *value=EOSMAC(eossimplezoomtableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN);
        if(incol==PTOTOFFSETin) *value=EOSMAC(eossimplezoomtableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN);
        if(incol==CHIOFFSETin) *value=EOSMAC(eossimplezoomtableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN);
        if(incol==STOTOFFSETin) *value=EOSMAC(eossimplezoomtableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOFFSET-FIRSTEOSDEGEN);

        if(utotdegencut[whichtablevar]<=DEGENCUTLASTOLDVERSION){
        }
        else{ // utotdegencut[whichtablevar]>DEGENCUTLASTOLDVERSION
          if(incol==UTOTINin) *value=EOSMAC(eossimplezoomtableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN);
          if(incol==PTOTINin) *value=EOSMAC(eossimplezoomtableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN);
          if(incol==CHIINin) *value=EOSMAC(eossimplezoomtableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN);
          if(incol==STOTINin) *value=EOSMAC(eossimplezoomtableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSIN-FIRSTEOSDEGEN);

          if(incol==UTOTOUTin) *value=EOSMAC(eossimplezoomtableextradegen,UTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN);
          if(incol==PTOTOUTin) *value=EOSMAC(eossimplezoomtableextradegen,PTOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN);
          if(incol==CHIOUTin) *value=EOSMAC(eossimplezoomtableextradegen,CHIDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN);
          if(incol==STOTOUTin) *value=EOSMAC(eossimplezoomtableextradegen,STOTDIFF,mmm,lll,kkk,jjj,iii,EOSOUT-FIRSTEOSDEGEN);
        }
      }
    }
  }
#endif

}





// determine whether should do log interpolation at the input time
static int get_dologinterp_in(int degentable, int whichfun)
{
  int whichinterp1,whichinterp2,loginterp;
  int qi;



#if(DOLOGINTERP)

  if(degentable){
    loginterp=1;
  }
  else{
    // GODMARK: Can make array that stores this info, looked up by whichfun as index
    // functions (F) F(rho0,u)
    whichinterp1=(whichfun==PofRHOCHIin||whichfun==UofRHOPin||whichfun==TEMPPin||whichfun==PofRHOUin||whichfun==CS2ofRHOUin||whichfun==SofRHOUin||whichfun==SSofRHOCHIin||(whichfun>=FIRSTEXTRAin && whichfun<=LASTEXTRAin)||whichfun==TEMPUin||whichfun==TEMPCHIin||whichfun==UofRHOSin||whichfun==TEMPSin||whichfun==UofRHOSin);
    // functions (F) F(rho0,p)
    whichinterp2=(whichfun==DPDRHOofRHOUin||whichfun==DPDUofRHOUin||whichfun==DSDRHOofRHOUin||whichfun==DSDUofRHOUin||whichfun==DSSDRHOofRHOCHIin||whichfun==DSSDCHIofRHOCHIin||whichfun==IDRHO0DPin||whichfun==IDCHIDPin);
    
    //dualfprintf(fail_file,"whichfun=%d whichinterp1=%d whichinterp2=%d\n",whichfun,whichinterp1,whichinterp2);
    
    if(whichinterp1) loginterp=1;
    else if(whichinterp2) loginterp=0;
    else{
      dualfprintf(fail_file,"Undefined whichfun=%d in get_eos_fromlookup_linear(): degentable=%d\n",whichfun, degentable);
      //      for(qi=1;qi<=NUMINDEPDIMENS+1;qi++) dualfprintf(fail_file,"%d : vartypearray=%d indexarray=%d\n",qi,vartypearray[qi],indexarray[qi]);
      myexit(3287623);
    }
  }
#else
  loginterp=0;
#endif



  return(loginterp);

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
