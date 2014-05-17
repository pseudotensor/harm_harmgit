////////////////////////////
//
// EOS from KAZ FULL
//
////////////////////////////


// PIPETODO:
// 1) Streamline existing lookup code to use pointer referneces to functions to avoid conditionals.
// 2) Only go to 1E-10 or something instead of 1E-14 for inversion




// non-dependent and dependent macros
//#include "kazfulleos.global.h" // included at top-level now

// small-sized definitions
#include "kazfulleos.defs.h"



// memory and pointers for EOS tables
#include "kazfulleos.eostablesdefs.h"


///////////////
//
// some globally (within this file) accessed functions
//
///////////////
static int getsingle_eos_fromtable(int whichfun, int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2, FTYPEEOS *answer);
static int get_eos_fromtable(int whichtablesubtype, int *iffun, int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2, FTYPEEOS *answers, int *badlookups);
static int get_eos_fromtable_noextrapolation(int whichtablesubtype, int *iffun, int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2, FTYPEEOS *answers, int *badlookups);
static int which_eostable(int whichtablesubtype, int ifdegencheck, int whichindep, int *vartypearraylocal, FTYPE *qarray, int *whichtablevar);
static void eos_lookup_degen(int begin, int end, int skip, int whichtablevar, int whichindep, int *vartypearraylocal, FTYPE *qarray, FTYPEEOS *indexarrayvar);
static void eos_lookup_degen_utotdegencut23(int begin, int end, int skip, int whichtablevar, int whichindep, int *vartypearraylocal, FTYPE *qarray, FTYPEEOS *indexarrayvar);
static void eos_lookup_prepost_degen(int whichdegen, int whichtablevar, int whichindep, int *vartypearraylocal, FTYPE *qarray, FTYPEEOS *indexarrayvar);
static int get_whichindep_fromwhichd(int whichd, int *whichindep);
static void eos_lookup_modify_qarray(int whichdegen, FTYPEEOS *myanswers, int whichtablevar, int whichindep, int *vartypearraylocal, FTYPE *qarray, FTYPEEOS *indexarrayvar);

static void lookup_yespecial(int whichtablevar, int vartypearraylocal, FTYPE qarray, FTYPEEOS *indexarrayvar);
static void lookup_log(int whichtablevar, int vartypearraylocal, FTYPE qarray, FTYPEEOS *indexarrayvar);

static void reset_repeatedfun(int isextratype, int whichd, int numcols, int firsteos, int *repeatedfun);

static void pre_truncate_qarray(int ifdegencheck, int whichtabletry, int *vartypearraylocal, FTYPEEOS *qarray);



// include C files for simplicity
#include "kazfulleos_set_arrays.c"
#include "kazfulleos_read_setup_eostable.c"
#include "kazfulleos_eosfunctions.c"
#include "kazfulleos_lookuptable_interpolation.c"







// called during read_setup_eostable()
void initeos_kazfulleos(void)
{
#if(ALLOWKAZEOS)
  
  // initialize repeated qarray's
  int qi,whichd,extrai;
  for(whichd=0;whichd<NUMEOSDEGENQUANTITIESMEM1;whichd++){
    for(extrai=0;extrai<NUMEXTRATABLETYPES;extrai++) for(qi=FIRSTINDEPDIMEN;qi<=LASTINDEPDIMENUSED;qi++) qoldarray[extrai][whichd][qi]=-BIG;
    //    kaziiowhichd[whichd]=kazjjowhichd[whichd]=kazkkowhichd[whichd]=kazllowhichd[whichd]=kazmmowhichd[whichd]=INITKAZINDEX;
  }
  for(qi=FIRSTINDEPDIMEN;qi<=LASTINDEPDIMENUSED;qi++) qoldarrayextras[qi]=-BIG; // all same UEOS independent variable, so no need for [whichd] dependence
  doallextrasold=-1;
#endif
  
}





// see if need to extrapolate values before doing lookup
static int get_eos_fromtable(int whichtablesubtype, int *iffun, int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2, FTYPEEOS *answers, int *badlookups)
{
  int doextrap;
  FTYPE quant1temp1,quant1temp2;
  FTYPEEOS answers1[MAXEOSPIPELINE];
  int badlookups1[MAXEOSPIPELINE];
  FTYPEEOS answers2[MAXEOSPIPELINE];
  int badlookups2[MAXEOSPIPELINE];
  int numcols,coli;
  int failreturn;
  int tableiter;



  // no failures by default
  failreturn=0;
  // no extrapolation by default
  doextrap=0;


  if(EXTRAPOLATEHIGHRHOU){


    // check if need to extrapolate
    if(quant1>rhoupperlimit){
      // then need to extrapolate instead of assuming constancy as will happen if just pass to lookup
      quant1temp1=0.9*rhoupperlimit;
      quant1temp2=0.99*rhoupperlimit;
      doextrap=1;
    }
    else{
      quant1temp1=quant1;
      quant1temp2=quant1;
    }

    
    // do extrapolation if required
    if(doextrap){
      failreturn=get_eos_fromtable_noextrapolation(whichtablesubtype, iffun, whichd, EOSextra, quant1temp1, quant2, answers1, badlookups1);
      failreturn+=get_eos_fromtable_noextrapolation(whichtablesubtype, iffun, whichd, EOSextra, quant1temp2, quant2, answers2, badlookups2);
    }
  }




  if(doextrap && EXTRAPOLATEHIGHRHOU==1){
    // then do extrapolation from good lookups
    // get final extrapolated value using bi-linear interpolation with only 2 values
    numcols = numcolintablesubtype[whichtablesubtype]; // ok use of numcolintablesubtype

    FTYPE divisor=1.0/(quant1temp2-quant1temp1);

    // now get extrapolation
    for(coli=0;coli<numcols;coli++){
      if(iffun[coli]){
        if(badlookups1[coli]==0 && badlookups2[coli]==0){
          badlookups[coli]=0;
          answers[coli] = answers1[coli]*(quant1temp2-quant1) + answers2[coli]*(quant1-quant1temp1);
          // apply divisor
          answers[coli] =  answers[coli]*divisor;
        }
        else badlookups[coli]=1; // indicate one of extrapolated versions was bad lookup
      }// end iffun
    }// end over coli

  }
  else{
    // normal non-extrapolated lookup
    failreturn=get_eos_fromtable_noextrapolation(whichtablesubtype, iffun, whichd, EOSextra, quant1, quant2, answers, badlookups);
  }


  return(failreturn);
}






// non-pipelined way to get single value from pipelined lookup
// calls pipelined version with 1 value requested
// so automatically does extrapolation
static int getsingle_eos_fromtable(int whichfun, int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2, FTYPEEOS *answer)
{
  FTYPEEOS answers[MAXEOSPIPELINE];
  int badlookups[MAXEOSPIPELINE];
  int iffun[MAXEOSPIPELINE];
  int numcols,coli;
  int whichtablesubtype=whichtablesubtypeinquantity[whichfun];
  int whichcol=whichcolinquantity[whichfun];

  ///////////
  //
  // can pipeline up to numcols quantities
  //
  ///////////
  numcols = numcolintablesubtype[whichtablesubtype]; // ok use of numcolintablesubtype
  // default iffun=0
  for(coli=0;coli<numcols;coli++) iffun[coli]=0;
  // correct iffun
  iffun[whichcol]=1;

  get_eos_fromtable(whichtablesubtype,iffun,whichd,EOSextra,quant1,quant2,answers,badlookups);

  *answer = answers[whichcol];

  return(badlookups[whichcol]); // can do this for just 1 quantity

}









////////////////////////////
//
// Choose behavior of loop depending upon whether including degen offset
//
////////////////////////////
// GODMARK: Could set this 
#if(ALLOWDEGENOFFSET)
// then need to lookup degen part first
#define whichdegenstart 1
#define whichdegenend 0
#else
#define whichdegenstart 0
#define whichdegenend 0
#endif


// for now assume TDYNORYE is locked to timestep that is order 1E-6 seconds, so can use lowest leos
// which = which function looking up
// q1 = rhob
// q2 = u or P or \chi
// notice that q1-q2 are pass by value, so internally changed but externally remains unchanged

// iffun[0...] = 0 or 1 depending upon if want that quantity in that subtable
// answers[0...] is filled with answer if iffun==1
// for now assume can only do one whichd type per call, so quantities within table subtypes have to agree on their independent variables
static int get_eos_fromtable_noextrapolation(int whichtablesubtype, int *iffun, int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2, FTYPEEOS *answers, int *badlookups)
{
  int whichindep;
  int whichdegen;
  FTYPEEOS myanswers[2][MAXEOSPIPELINE]; // for [2]: 0 = non-degen answer  1=degen answer
  FTYPE qarray[NUMINDEPDIMENSMEM];
  FTYPE qfloor[NUMINDEPDIMENSMEM];
  int i,qi,coli;
  int repeatedeos;
  int numcols;
  int failreturn=0;
  int vartypearraylocal[NUMINDEPDIMENSMEM];
  int localiffun[MAXEOSPIPELINE];
  int localbadlookups[MAXEOSPIPELINE];
  int localwhichtablesubtype;
  int isextratype;
  int firsteos;
  int localnumcols;




  ///////////////////////////////
  //
  // Ensure that setup table
  //
  ///////////////////////////////
#if(PRODUCTION==0)
  if(didsetupkazeos==0){
    dualfprintf(fail_file,"Requested WHICHEOS=KAZEOS but didn't initialize kaz table by calling read_setup_eostable()\n");
    myexit(151068); 
  }

  // For repeatedeos, don't want coincidence where (e.g.) utotdiff==stotdiff and think repeated.  Need independents to be same.  Now controlled with using [whichd] on qoldarray
  if(whichd!=whichdintablesubtype[whichtablesubtype] && whichdintablesubtype[whichtablesubtype]!=NOSUCHDIFF ){
    dualfprintf(fail_file,"lookup not setup for different whichd=%d in this whichtablesubtype=%d with whichd=%d\n",whichd,whichtablesubtype,whichdintablesubtype[whichtablesubtype]);
  }
#endif



 
  //////////////////////
  //
  // Determine which independent variable to use for "temperature": utot, ptot, or chi (or utotdiff, ptotdiff, chidiff for degen case)
  //
  //////////////////////
  get_whichindep_fromwhichd(whichd,&whichindep);
  // used to translate q1-q5 to INDEP0-INDEP6, and used with limits and sizes of tables
  //  for(qi=FIRSTINDEPDIMEN;qi<=LASTINDEPDIMEN;qi++){
  for(qi=FIRSTINDEPDIMEN;qi<=LASTINDEPDIMEN;qi++){
    vartypearraylocal[qi]=vartypearray[qi]; // use of global, which is ok since unchanging
  }
  // OVERRIDE global value for temperature-like quantity
  vartypearraylocal[TEMPLIKEINDEP]=whichindep;



  ///////////
  //
  // can pipeline up to numcols quantities
  //
  ///////////
  numcols = numcolintablesubtype[whichtablesubtype]; // ok use of numcolintablesubtype
  isextratype=isextraintablesubtype[whichtablesubtype];
  firsteos=firsteosintablesubtype[whichtablesubtype];

  for(coli=0;coli<numcols;coli++){
    if(iffun[coli]) badlookups[coli]=1; // all bad by default.  Using this to check if need to get this coli during processing
    else badlookups[coli]=0; // unwanted are not bad
    // SO DO NOT CHANGE iffun[] here or in subfunctions.  Also do NOT default badlookups to 1 in subfunctions since already done here and use this to see what other functions yet to do.
    // So can only change badlookups to 0, can't set to 1
  }


  ///////////////////////////////
  //
  // nuclear offset (as consistent with HELM code to obtain correct energy per baryon)
  //
  ///////////////////////////////
  offsetquant2(+1.0,whichd, EOSextra, quant1, quant2, &quant2);



  ////////////////////////////
  //
  // Set array of independent EOS quantities from qarray and EOSextra[]
  // Note that input is not qarray[], so can change qarray[] as required during EOS routine (e.g. for truncation, etc.)
  //
  ////////////////////////////
  qarray[RHOINDEP]=quant1;
  qarray[TEMPLIKEINDEP]=quant2;
  // qarray[TEMPLIKEINDEP+?] are always stored in EOSextra
  for(qi=TEMPLIKEINDEP+1;qi<=LASTINDEPDIMENUSED;qi++){
    qarray[qi] = EOSextra[vartypeeosextraarray[qi]];
  }


  //////////////////////////////////////////////////////
  //
  // check if already got here with these q values (done after assignment of q3-q5, so distinct -- done before floors on q since they make less distinct anyways)
  //
  // Note that repeatedeos uses only qarray and whichd, but not based upon whichtablesubtype
  // This is so that in principle input independents could be same but switch to another subtable so don't have to re-get indexarray[]
  //
  // repeated entry can occur as long as:
  // 1) Same whichd
  // 2) Same qarray
  // 3) Resolves to same type of table (i.e. extra version or not extra version) since this is case where "whichd" and "qarray" and "whichtable" are not enough information for a given single function or whichtablesubtype.
  //
  // This is done instead of having (e.g.) qoldarray[MAXEOSPIPELINE] since many subtypes have same whichd and are on same range and want to be able to repeat across table subtypes with same "whichd".  extra generally on different range, which is why it has its own array index
  //
  // Then have old result stored that came from some table as stored in whichtable[]
  //
  // Result is in repeatedfun[] and resultold[] that *are* [MAXEOSPIPELINE] since don't care where answer came from as long as have one in storage
  //
  //////////////////////////////////////////////////////


  repeatedeos=1;
  for(qi=FIRSTINDEPDIMEN;qi<=LASTINDEPDIMENUSED;qi++){ // only need to repeat used independent variables, not all
    repeatedeos*=(fabs(qarray[qi]-qoldarray[isextratype][whichd][qi])<OLDTOLERANCE);
  }




  //////////////////////////////////////////////////////
  //
  // Compute result or retrieve already computed result
  //
  //////////////////////////////////////////////////////
  if(repeatedeos){

    //////////
    //
    // grab old result if exists and reset badlookups->0 since no longer need to get that quantity
    //
    //////////
    for(coli=0;coli<numcols;coli++){
      // check if already got the function wanted
      if(repeatedfun[firsteos+coli] && iffun[coli]){
        // choice of table irrelevant if already found solution
        answers[coli]=resultold[firsteos+coli];
        badlookups[coli]=0; // no longer need to get this one
      }
    }


    //////////
    //
    // check if still quantities to get
    //
    //////////
    int stilltoget=0; for(coli=0;coli<numcols;coli++) if(iffun[coli] && badlookups[coli]==1) stilltoget++;


    //////////
    //
    // get quantities not stored in repeated array but grid index is same as previous since repeatedeos
    //
    //////////
    if(stilltoget){
      // then get result

      // whatever is done after degen table lookup and normal table lookup
      whichdegen=whichdegenend;
      
      // check if previously did get table result and previously used degeneracy and non-degeneracy tables fell within same table type (full, simple, simplezoom)
      if(whichtable[isextratype][whichd][ISNOTDEGENTABLE]==NOTABLE || whichtable[isextratype][whichd][ISDEGENTABLE]!=whichtable[isextratype][whichd][ISNOTDEGENTABLE]){
        diag_eosfaillookup((int)EOSextra[IGLOBAL],(int)EOSextra[JGLOBAL],(int)EOSextra[KGLOBAL]);
        if(debugfail>=2){
          dualfprintf(fail_file,"Stilltoget, but tables not consistent or bad: i=%d j=%d k=%d notdegen=%d degen=%d\n",(int)EOSextra[IGLOBAL],(int)EOSextra[JGLOBAL],(int)EOSextra[KGLOBAL],whichtable[isextratype][whichd][ISDEGENTABLE],whichtable[isextratype][whichd][ISNOTDEGENTABLE]);
          for(coli=0;coli<numcols;coli++){
            if(iffun[coli] && badlookups[coli]==1){
              dualfprintf(fail_file,"Stilltoget: whichtablesubtype=%d coli=%d\n",whichtablesubtype,coli);
            }
          }
        }
        // badlookups[all coli] will remain for all bad ones
        // return now since nothing else to do
        return(0); // but not hard failure
      }
      else{
        // DEBUG:
        //if(debugfail>=2) dualfprintf(fail_file,"REPEATFROMLOOKUP\n"); // DEBUG
        //if(debugfail>=2) dualfprintf(fail_file,"CHECKLOOKUP:(repeated) %d %d %d : %d %d\n",(int)EOSextra[IGLOBAL],(int)EOSextra[JGLOBAL],(int)EOSextra[KGLOBAL],whichdegen,whichfun);


        // tell get_eos_fromlookup() which functions to try and get.  Better than having to look at badlookups inside get_eos_fromlookup()
        for(coli=0;coli<numcols;coli++) if(iffun[coli]&&badlookups[coli]) localiffun[coli]=1; else localiffun[coli]=0;

        // do lookup directly using old indexarray that is valid since independents unchanged
        failreturn += get_eos_fromlookup(repeatedeos,WHICHEOSDIMEN,whichdegen, whichtable[isextratype][whichd][whichdegen], whichtablesubtype, localiffun, whichindep, quant1, vartypearraylocal, indexarray[isextratype][whichd], myanswers[whichdegen],badlookups);

        // see if good lookup
        for(coli=0;coli<numcols;coli++){
          if(localiffun[coli]){ // only look at new coli's
            if(badlookups[coli]){
              diag_eosfaillookup((int)EOSextra[IGLOBAL],(int)EOSextra[JGLOBAL],(int)EOSextra[KGLOBAL]);
              if(debugfail>=2) dualfprintf(fail_file,"Repeated lookup got no valid lookup in table: i=%d j=%d k=%d\n",(int)EOSextra[IGLOBAL],(int)EOSextra[JGLOBAL],(int)EOSextra[KGLOBAL]); // actually bad failure
            } // end if bad lookup
            else{
              // setup return of answer
              answers[coli]=myanswers[whichdegen][coli];

              // now save result if got result so can use next time
              resultold[firsteos+coli]=myanswers[whichdegen][coli];
              repeatedfun[firsteos+coli]=1;
            }// end else if good lookup
          }// end localiffun==1
        }// end over coli

      }// end else if correct table and did previously get a result from table
    }// end if still quantity to get
    else{

      // then done since nothing left to get -- it was all inside the repeated array

    } // end else if done





  }// end if repreatedeos
  else{ // else if not repeatedeos, so need to do full lookup





    /////////////
    //
    // setup old values for next time since doing new lookup
    //
    /////////////
    for(qi=FIRSTINDEPDIMEN;qi<=LASTINDEPDIMENUSED;qi++) qoldarray[isextratype][whichd][qi]=qarray[qi]; // should be qarray *before* any degeneracy subtraction on qarray[TEMPLIKEINDEP]
    // once values of independent variables changes, reset so that none of quantities in subtable are set as repeated, no matter iffun==0 or 1
    reset_repeatedfun(isextratype,whichd,numcols,firsteos,repeatedfun);
    // if not repeated, then reset whichtable="prior table used" to access the subtabletype
    whichtable[isextratype][whichd][ISNOTDEGENTABLE]=whichtable[isextratype][whichd][ISDEGENTABLE]=NOTABLE;




    ////////////////////////////
    //
    // Obtain interpolated function
    //
    // 1) Get whichtable for degen table -- which determines ieos,keos,leos
    // 2) Determine jeos after have degen offset from interpolation of ieos,keos,leos
    //
    ////////////////////////////
    for(whichdegen=whichdegenstart;whichdegen>=whichdegenend;whichdegen--){


      //////////////
      //
      // tell get_eos_fromlookup() which functions to try and get.  Better than having to look at badlookups inside get_eos_fromlookup()
      // if degen type lookup, get *all* degens
      //
      //////////////
      if(whichdegen==ISNOTDEGENTABLE){
        localnumcols=numcols;
        for(coli=0;coli<localnumcols;coli++) if(iffun[coli]&&badlookups[coli]){ localwhichtablesubtype=whichtablesubtype; localiffun[coli]=1; localbadlookups[coli]=1;} else { localwhichtablesubtype=whichtablesubtype; localiffun[coli]=0; localbadlookups[coli]=0;}
      }
      else{
        localnumcols=get_numcols(whichtable[isextratype][whichd][whichdegen],SUBTYPEDEGEN);
        for(coli=0;coli<localnumcols;coli++) {localwhichtablesubtype=SUBTYPEDEGEN; localiffun[coli]=1; localbadlookups[coli]=1;}
      }


      

      // Once looked-up ieos,keos,leos, don't need to do full lookup again (but do need to do full interpolation "again" that is on different set of quantities)
      // that is, only jeos changes
      // see if within table, and if so lookup interpolated result

      // can't predict that when taking subtraction on utot -> utotdiff that will end up within same table, so always get both tables
      failreturn += which_eostable(whichtablesubtype,whichdegen,whichindep, vartypearraylocal, qarray, &whichtable[isextratype][whichd][whichdegen]);

      // check if got a table, and if so then use it, otherwise return(1) and assume using off-table values
      if(whichtable[isextratype][whichd][whichdegen]==NOTABLE){
        if(debugfail>=3){
          dualfprintf(fail_file,"subtype=%d whichtable=%d : whichd=%d whichdegen=%d :: ig=%d jg=%d kg=%d\n",whichtablesubtype,whichtable[isextratype][whichd][whichdegen],whichd,whichdegen,(int)EOSextra[IGLOBAL],(int)EOSextra[JGLOBAL],(int)EOSextra[KGLOBAL]);
          for(qi=FIRSTINDEPDIMEN;qi<=LASTINDEPDIMENUSED;qi++) dualfprintf(fail_file,"qarray[%d]=%21.15g\n",qi,qarray[qi]);
          for(coli=0;coli<localnumcols;coli++) if(iffun[coli]) dualfprintf(fail_file,"coli=%d bad=%d\n",coli,badlookups[coli]);
        }
        diag_eosfaillookup((int)EOSextra[IGLOBAL],(int)EOSextra[JGLOBAL],(int)EOSextra[KGLOBAL]);
        if(debugfail>=2) dualfprintf(fail_file,"qarray not within table: i=%d j=%d k=%d\n",(int)EOSextra[IGLOBAL],(int)EOSextra[JGLOBAL],(int)EOSextra[KGLOBAL]);
        // badlookups will be used since didn't change to 0
        // return now since nothing else to do
        return(0); // but not hard failure

      }
      else{

        // check that degen and normal table resolved to same table, since can't use degen value from one table to get values from another table
        if(whichdegen==whichdegenend){

          if(whichtable[isextratype][whichd][ISNOTDEGENTABLE]==NOTABLE || whichtable[isextratype][whichd][ISNOTDEGENTABLE]!=whichtable[isextratype][whichd][ISDEGENTABLE]){
            diag_eosfaillookup((int)EOSextra[IGLOBAL],(int)EOSextra[JGLOBAL],(int)EOSextra[KGLOBAL]);
            if(debugfail>=2) dualfprintf(fail_file,"Degen and normal table selections different: normal=%d degen=%d i=%d j=%d k=%d\n",whichtable[isextratype][whichd][ISNOTDEGENTABLE],whichtable[isextratype][whichd][ISDEGENTABLE],(int)EOSextra[IGLOBAL],(int)EOSextra[JGLOBAL],(int)EOSextra[KGLOBAL]);
            // badlookups will be used since didn't change to 0
            // return now since nothing else to do
            return(0); // but not hard failure
          }
        }



        ////////////
        //
        // Lookup the positions: ieos,jeos,keos,leos for a given table, independent variable definition, and indep values "qarray"
        //
        // for whichdegen==1, below does set indexarray[][2]=0 for accessing degen table
        //
        ////////////
        eos_lookup_prepost_degen(whichdegen, whichtable[isextratype][whichd][whichdegen], whichindep, vartypearraylocal, qarray, indexarray[isextratype][whichd]);



        ////////////////////////////
        //
        // now compute result
        //
        // For whichdegen==1, myanswers[whichdegen==1] then contains either {utot/ptot/chitot/stot}offset for utotdegencut==0,1 OR contains 3 quantities corressponding to {utot,ptot,chitot,stot}{0,in,out} used by eos_lookup_modify_qarray() to obtain indexed version of qarray[TEMPLIKEINDEP]
        //
        // For whichdegen==0, presume qarray[TEMPLIKEINDEP] correct for lookup and so lookup quantities within subtable=whichtablesubtype and lookup those quantities listed as iffun[whichfun]=1
        //
        ////////////////////////////
 

        //////////////
        //
        // do either non-degen or degen lookup
        //
        //////////////
        failreturn=get_eos_fromlookup(repeatedeos,WHICHEOSDIMEN,whichdegen, whichtable[isextratype][whichd][whichdegen], localwhichtablesubtype, localiffun, whichindep, quant1, vartypearraylocal, indexarray[isextratype][whichd],myanswers[whichdegen],localbadlookups);


        ////////////
        //
        // Check if bad lookup
        //
        ////////////
        for(coli=0;coli<localnumcols;coli++){
          if(localiffun[coli] && localbadlookups[coli]==1){
            if(debugfail>=3){
              dualfprintf(fail_file,"Looked up, but no answer within table: coli=%d out of %d\n",coli,localnumcols);
              dualfprintf(fail_file,"subtype=%d whichtable=%d : whichd=%d whichdegen=%d :: ig=%d jg=%d kg=%d\n",whichtablesubtype,whichtable[isextratype][whichd][whichdegen],whichd,whichdegen,(int)EOSextra[IGLOBAL],(int)EOSextra[JGLOBAL],(int)EOSextra[KGLOBAL]);
              for(qi=FIRSTINDEPDIMEN;qi<=LASTINDEPDIMENUSED;qi++) dualfprintf(fail_file,"qarray[%d]=%21.15g\n",qi,qarray[qi]);
            }
            diag_eosfaillookup((int)EOSextra[IGLOBAL],(int)EOSextra[JGLOBAL],(int)EOSextra[KGLOBAL]);
            if(debugfail>=2) dualfprintf(fail_file,"Looked up, but no answer within table: coli=%d out of %d i=%d j=%d k=%d\n",coli,localnumcols,(int)EOSextra[IGLOBAL],(int)EOSextra[JGLOBAL],(int)EOSextra[KGLOBAL]);

          }
        }


        // transfer over badlookups if final non-degen quantity
        if(whichdegen==0) for(coli=0;coli<numcols;coli++) if(iffun[coli]&&badlookups[coli] && localbadlookups[coli]==0) badlookups[coli]=0;


        if(whichdegen==1){
          ////////////////////////////
          //
          // translate qarray[TEMPLIKEINDEP] from original value to value required to do final temperature-like quantity lookup
          //
          // myanswers[whichdegen] should have multiple values with utotdegencut>=2
          //
          ////////////////////////////
          eos_lookup_modify_qarray(whichdegen, myanswers[whichdegen],whichtable[isextratype][whichd][whichdegen],whichindep,vartypearraylocal,qarray,indexarray[isextratype][whichd]);
        }


      } // end if whichtable!=NOTABLE and so qarray[] within some table's independent variables
    }// end loop over degenerate and then normal table (or if ALLOWDEGENOFFSET==0, then only for normal table)
    

    // finally, normal table lookup gives answer over *all* originally desired coli's
    for(coli=0;coli<numcols;coli++){
      if(iffun[coli] && badlookups[coli]==0){
        answers[coli]=myanswers[ISNOTDEGENTABLE][coli];
        // setup repeated for next time
        resultold[firsteos+coli]=answers[coli];
        repeatedfun[firsteos+coli]=1;

        // DEBUG:
        // if(whichtablesubtype==SUBTYPEEXTRA){
        //   dualfprintf(fail_file,"coli=%d answers=%21.15g\n",coli,answers[coli]);
        // }


      }
    }


    // DEBUG:
    //    if(whichtablesubtype==SUBTYPEEXTRA){
    //      for(qi=FIRSTINDEPDIMEN;qi<=LASTINDEPDIMENUSED;qi++) dualfprintf(fail_file,"qarray[%d]=%21.15g index[ise=%d][whichd=%d][qi=%d]=%21.15g\n",qi,qarray[qi],isextratype,whichd,qi,indexarray[isextratype][whichd][qi]);
    //    }



  }// end else if not repeated






  ///////////////////////////////
  //
  // apply nuclear offset to functions that have been stored in table with offset so that HARM uses *true* function values
  //
  // Note that this offset is applied after any repeated function values are stored.  So must be done whether repeated lookup or not repeated lookup.
  //
  ///////////////////////////////
  int whichdfake;
  for(coli=0;coli<numcols;coli++){
    if(iffun[coli] && badlookups[coli]==0){
      if(needspostoffset(whichd,whichtablesubtype,coli,&whichdfake)){
        offsetquant2(-1.0,whichdfake, EOSextra, quant1, answers[coli], &answers[coli]);
      }
    }
  }




  /////////////////
  //
  // Report failure status
  //
  /////////////////
  return(failreturn);


}





// reset repeatedfun to zero for all tables associated with both same isextratype and whichd
// use of repeated function depends upon same whichd and isextratype
static void reset_repeatedfun(int isextratype, int whichd, int numcols, int firsteos, int *repeatedfunlocal)
{
  int coli;


  if(isextratype){

    // note that one has no access to extra degen table

    if(whichd==UTOTDIFF) for(coli=0;coli<numcols;coli++) repeatedfunlocal[firsteos+coli]=0;
    else{
      dualfprintf(fail_file,"No such whichd=%d for isextratype=%d\n",whichd,isextratype);
      myexit(346262121);
    }
  }
  else{
    int numcolslocal,firsteoslocal;
    int subtype;
    for(subtype=0;subtype<NUMTABLESUBTYPES;subtype++){
      if(whichd==whichdintablesubtype[subtype]){
        // then this subtype should be reset
      
        numcolslocal = numcolintablesubtype[subtype]; // ok use of numcolintablesubtype
        firsteoslocal=firsteosintablesubtype[subtype];
        for(coli=0;coli<numcolslocal;coli++) repeatedfunlocal[firsteoslocal+coli]=0;
      }
    }
  }


}





// which_eostable() is written for each type of dataset(s) to used -- too complicated to make general

// check if density and pressure-like quantity are within table
// all quantities are real code units (i.e. linear mapping instead of log10 mapping)
// here we only presume to care about density and internal energy, while H and T are presumed to be truncated rather than extended with some alternative
// ifdegencheck: 0 = normal table check  1 = ignores q2 (u,p,chi) since generally q2 is largest range possible and later will restrict/check if within the full table
// whichindep = which independent variable (based upon which function looking up)
static int which_eostable(int whichtablesubtype, int ifdegencheck, int whichindep, int *vartypearraylocal, FTYPE *qarray, int *whichtablevar)
{
  int whichtabletry;
  // don't assume tables are setup so q1 and q2 always have same range
  // assume all quantities depend on limits in same way, so "which" doesn't matter
  // assume HEOS and YEEOS are always withing range for now


  /////////////////////////////
  //
  // Note that q3,q4,q5 ((TDYN or YE),(TDYN or YNU), Hcm) are computed such that is forced to be within FULLTABLE limits so consistently used
  //
  // Note that for q2, for utotdegencut>=2 that qarray[TEMPLIKEINDEP] is actually an fractional index "i/N"=lutotdiff, and lineartablelimits[] is correctly that range of "i/N"=lutotdiff from 0..1.0
  /////////////////////////////



#if(ALLOWFULLTABLE)
  if(whichtablesubtype==SUBTYPEEXTRA && WHICHDATATYPEGENERAL==4){
    whichtabletry=FULLTABLEEXTRA;
  }
  else whichtabletry=FULLTABLE;


  // pre-truncate some quantities (should change qarray upon returning!)
  pre_truncate_qarray(ifdegencheck, whichtabletry, vartypearraylocal, qarray);
      
  if(qarray[1]>=lineartablelimits[whichtabletry][vartypearraylocal[1]][0] && qarray[1]<=lineartablelimits[whichtabletry][vartypearraylocal[1]][1]
     &&
     (ifdegencheck || qarray[TEMPLIKEINDEP]>=lineartablelimits[whichtabletry][vartypearraylocal[TEMPLIKEINDEP]][0] && qarray[TEMPLIKEINDEP]<=lineartablelimits[whichtabletry][vartypearraylocal[TEMPLIKEINDEP]][1])
     ){
    *whichtablevar=whichtabletry;
    return(0);
  }
#endif
#if(ALLOWSIMPLETABLE)
  if(whichtablesubtype==SUBTYPEEXTRA && WHICHDATATYPEGENERAL==4){
    whichtabletry=SIMPLETABLEEXTRA;
  }
  else whichtabletry=SIMPLETABLE;


  // pre-truncate some quantities (should change qarray upon returning!)
  pre_truncate_qarray(ifdegencheck, whichtabletry, vartypearraylocal, qarray);

    
  if(qarray[1]>=lineartablelimits[whichtabletry][vartypearraylocal[1]][0] && qarray[1]<=lineartablelimits[whichtabletry][vartypearraylocal[1]][1]
     &&
     (ifdegencheck || qarray[TEMPLIKEINDEP]>=lineartablelimits[whichtabletry][vartypearraylocal[TEMPLIKEINDEP]][0] && qarray[TEMPLIKEINDEP]<=lineartablelimits[whichtabletry][vartypearraylocal[TEMPLIKEINDEP]][1])
     ){
    *whichtablevar=whichtabletry;
    return(0);
  }
#endif


  if(debugfail>=2){ // DEBUG: was turned on when debugging EOS
    dualfprintf(fail_file,"NOT IN LOOKUP: ifdegencheck=%d whichindep=%d qarray[1]=%21.15g qarray[TEMPLIKEINDEP]=%21.15g\n",ifdegencheck,whichindep,qarray[1],qarray[TEMPLIKEINDEP]);
  }
  *whichtablevar=NOTABLE;
  return(1);

  // GODMARK:
  // alternative to this function is that we simply place a floor on the linear values of rho, u, H, and T so always within table values
  // Use: lineartablelimits[ii][0,1]
  //
  // problem is that for inversion, truncation still will be an if check so as slow as above anyways, so can just truncate here if wanted and then more general.  Only slows down simpler single-type calls to EOS.

}




// truncate some quantities before doing lookup of which table will use to do lookup
static void pre_truncate_qarray(int ifdegencheck, int whichtabletry, int *vartypearraylocal, FTYPEEOS *qarray)
{
  int qi;



  for(qi=1;qi<=NUMINDEPDIMENS;qi++){

    //////////////////////////////////
    // always fully truncate Y_e and Y^0_\nu
    //////////////////////////////////
    if(qi==YEINDEP || qi==YNUINDEP){

      if(qarray[qi]<1.000001*lineartablelimits[whichtabletry][vartypearraylocal[qi]][0]) qarray[qi]=1.000001*lineartablelimits[whichtabletry][vartypearraylocal[qi]][0];
      if(qarray[qi]>0.999999*lineartablelimits[whichtabletry][vartypearraylocal[qi]][1]) qarray[qi]=0.999999*lineartablelimits[whichtabletry][vartypearraylocal[qi]][1];
    }
  }


  if(TRUNCATEHIGHRHOU){
    // SUPERNOTE: Note we truncate rho to be within table if at high densities!  This is to best approximate behavior at very high densities should they form instead of reverting to (bah) ideal gas EOS or such things.
    // SUPERNOTE: We also truncate temperature-like quantities also since also bad to revert to ideal gas for most quantities.
    // SUPERNOTE: Otherwise, would have to extend simple table beyond 1E15 where no nuclear data.
    if(qarray[RHOINDEP]>lineartablelimits[whichtabletry][vartypearraylocal[RHOINDEP]][1]) qarray[RHOINDEP]=lineartablelimits[whichtabletry][vartypearraylocal[RHOINDEP]][1];
    if(ifdegencheck==0 && qarray[TEMPLIKEINDEP]>lineartablelimits[whichtabletry][vartypearraylocal[TEMPLIKEINDEP]][1]) qarray[TEMPLIKEINDEP]=lineartablelimits[whichtabletry][vartypearraylocal[TEMPLIKEINDEP]][1];
    
  }


}




// rho0,u,H,T have independent lookups for mapping to i,j,k,l and already done in eos_lookup_degen()
// u here stands for utotdiff,ptotdiff,and chidiff depending upon whichindep
// here only need to look up jeos associated with u
// notice that qarray and vartypearraylocal start at 1 not 0
static void eos_lookup_degen(int begin, int end, int skip, int whichtablevar, int whichindep, int *vartypearraylocal, FTYPE *qarray, FTYPEEOS *indexarrayvar)
{
  int qi;

  // is this expensive?
  // log_base(R-R0)
  // avoid nan's by choosing minimal density

  // old simple log10 way:
  // i = (x - x0)*[(N-1)/(x1-x0)] such that if x=x0, then i=0, if x=x1, then i=N-1
  
  // new generalized log way:
  // i = (log_base(R-R0) - x_in)/dx

  for(qi=begin;qi<=end;qi++){
    if(qi==skip) continue;

    if(qi==YEINDEP && whichyelooptype[whichtablevar]==YELOOPTYPESPECIAL){
      lookup_yespecial(whichtablevar,vartypearraylocal[qi],qarray[qi],&indexarrayvar[qi]);
    }
    else{
      lookup_log(whichtablevar,vartypearraylocal[qi],qarray[qi],&indexarrayvar[qi]);
    }


    // DEBUG:
    // dualfprintf(fail_file,"whichtablevar=%d qi=%d var=%d q=%21.15g index=%21.15g\n",whichtablevar,qi,vartypearraylocal[qi],qarray[qi],indexarrayvar[qi]);


  }


}



// get i(value)
static void lookup_log(int whichtablevar, int vartypearraylocal, FTYPE qarray, FTYPEEOS *indexarrayvar)
{
  FTYPEEOS logq;
  FTYPEEOS prelogq;


  prelogq = MAX(qarray    -lineartablelimits[whichtablevar][vartypearraylocal]  [3],SMALL);
  logq    = log10(prelogq)*lineartablelimits[whichtablevar][vartypearraylocal][2];
  
  *indexarrayvar = (logq   -tablelimits[whichtablevar][vartypearraylocal]  [0])*tablelimits[whichtablevar][vartypearraylocal]  [3];

  // DEBUG:
  //  dualfprintf(fail_file,"prelogq=%21.15g logq=%21.15g\n",prelogq,logq);
  //  dualfprintf(fail_file,"others: %21.15g %21.15g %21.15g %21.15g\n",lineartablelimits[whichtablevar][vartypearraylocal]  [3],lineartablelimits[whichtablevar][vartypearraylocal][2],tablelimits[whichtablevar][vartypearraylocal]  [0],tablelimits[whichtablevar][vartypearraylocal]  [3]);

}


// Get i(ye)
// see yetable.nb
// SUPERNOTE: Note we truncate Y_e to be within table!
static void lookup_yespecial(int whichtablevar, int vartypearraylocal, FTYPE qarray, FTYPEEOS *indexarrayvar)
{
  FTYPEEOS num = tablesize[whichtablevar][YEEOS];
  FTYPEEOS yei = lineartablelimits[whichtablevar][YEEOS][0];
  FTYPEEOS yef = lineartablelimits[whichtablevar][YEEOS][1];
  FTYPEEOS yegrid1 = eosyegrid1[whichtablevar];
  FTYPEEOS yegrid2 = eosyegrid2[whichtablevar];
  FTYPEEOS xgrid1 = eosxgrid1[whichtablevar];
  FTYPEEOS xgrid2 = eosxgrid2[whichtablevar];
  FTYPEEOS ye=qarray;

  // first range:
  if(ye<=yei){
    *indexarrayvar=0.0; // truncate
  }
  else if(ye<=yegrid1){
    *indexarrayvar = (num-1.0)*xgrid1*log10(ye/yei)/log10(yegrid1/yei);
  }
  else if(ye<=yegrid2){
    *indexarrayvar = (num-1.0)*(ye*(xgrid1-xgrid2) + xgrid2*yegrid1 - xgrid1*yegrid2)/(yegrid1-yegrid2);
  }
  else if(ye<=yef){
    *indexarrayvar = (num-1.0)*(ye*(xgrid2-1.0) + yegrid2 - xgrid2*yef)/(yegrid2-yef);
  }
  else{
    *indexarrayvar=num-1.0; // truncate
  }

}



// for this function, presume qarray[TEMPLIKEINDEP] is in lutotdiff="i/N" form, so simpler conversion to "i" units
static void eos_lookup_degen_utotdegencut23(int begin, int end, int skip, int whichtablevar, int whichindep, int *vartypearraylocal, FTYPE *qarray, FTYPEEOS *indexarrayvar)
{
  int qi;

  for(qi=begin;qi<=end;qi++){
    if(qi==skip) continue;
    
    indexarrayvar[qi] = (qarray[qi])*tablelimits[whichtablevar][vartypearraylocal[qi]]  [3]; // tablelimits[][][3] is N[qi], so gives index

  }


}


// whichdegen: 1 = do degen lookup  0 = do normal lookup without energy indepenent variable lookup
static void eos_lookup_prepost_degen(int whichdegen, int whichtablevar, int whichindep, int *vartypearraylocal, FTYPE *qarray, FTYPEEOS *indexarrayvar)
{

  void eos_lookup_degen(int begin, int end, int skip, int whichtablevar, int whichindep, int *vartypearraylocal, FTYPE *qarray, FTYPEEOS *indexarrayvar);
  int skip;
  int qi;




  if(whichdegen==ISNOTDEGENTABLE){
    // if normal table after degen table, then only need to get q2 location
    skip = NONEINDEP;// nonsense value to avoid skipping

    if(utotdegencut[whichtablevar]>DEGENCUTLASTOLDVERSION){
      eos_lookup_degen_utotdegencut23(TEMPLIKEINDEP,TEMPLIKEINDEP,skip, whichtablevar, whichindep, vartypearraylocal, qarray, indexarrayvar);
    }
    else{
      eos_lookup_degen(TEMPLIKEINDEP,TEMPLIKEINDEP,skip, whichtablevar, whichindep, vartypearraylocal, qarray, indexarrayvar);
    }
  }
  else{
    // if doing degentable, lookall up except q=UEOS/PEOS/CHIEOS/SEOS
    skip=TEMPLIKEINDEP; 
    eos_lookup_degen(FIRSTINDEPDIMEN,LASTINDEPDIMENUSED,skip,whichtablevar, whichindep, vartypearraylocal, qarray, indexarrayvar);
    indexarrayvar[TEMPLIKEINDEP]=0.0; // set TEMPLIKEINDEP as 0 since degentable has this as 0 (may be floating -> integer issue in how this floating value is used)
  }



  
  /////////////////
  //
  // enforce floor on temperature-like quantity
  // assume table has low density and internal energy and don't want to go much below this since then implied temperature becomes undefined
  // if density is outside table, can go to lower temperatures than in table
  //
  ///////////////////
#if(ALLOWDEGENOFFSET)
  if(indexarrayvar[TEMPLIKEINDEP]<0) indexarrayvar[TEMPLIKEINDEP]=0;

  /////////
  //
  // truncate upper temperatures.
  // have to truncate degen type u here since table checker doesn't know limits of degen table
  //
  /////////
  if(TRUNCATEHIGHRHOU){
    qi=TEMPLIKEINDEP;
    if(indexarrayvar[qi]>0.9999999*(tablesize[whichtablevar][vartypearraylocal[qi]]-1.0)){
      indexarrayvar[qi]=0.9999999*(tablesize[whichtablevar][vartypearraylocal[qi]]-1.0);
    }
  }
#endif



  //////////
  //
  // assume bottom of table is default for other untouched dimensions unless changed here
  //
  /////////
  for(qi=LASTINDEPDIMENUSED+1;qi<=NUMINDEPDIMENS;qi++){
    indexarrayvar[qi]=0.0;
  }



}




// If utotdegencut<=1: Translate qarray[TEMPLIKEINDEP] from {utot,ptot,chi,stot} -> {utotdiff,ptotdiff,chidiff,stotdiff}
// If utotdegencut>=2: Translate qarray[TEMPLIKEINDEP] from {utot,ptot,chi,stot} -> {lutotdiff,lptotdiff,lchidiff,lstotdiff} that actually correspond to "i/N" : the fractional index of the grid
// vartypearraylocal[TEMPLIKEINDEP] resolves to one of UEOS,PEOS,CHIEOS,SEOS as set globally when starting lookup
static void eos_lookup_modify_qarray(int whichdegen, FTYPEEOS *myanswers, int whichtablevar, int whichindep, int *vartypearraylocal, FTYPE *qarray, FTYPEEOS *indexarrayvar)
{

  /////////////////////////
  //      
  // subtract offset from actual code value to see where within table we are in terms of the offset
  // when doing normal table this change won't matter (i.e. q2 not used again, and meaning of q2 is undefined when mixing function with independent variable, but avoid if statement)
  // if whichdegen==1, after below line then q2 contains utotoffset, ptotoffset, or chioffset consistent with independent variables used for EOSQUANTITIES



  if(utotdegencut[whichtablevar]<=1){
    // always true that ALLOWDEGENOFFSET if here with whichdegen==1
    // don't allow u<utotdiff that would imply T<~0
    // here myanswers[EOSOFFSET] = utotoffset
    FTYPEEOS UTOTOFFSET=myanswers[EOSOFFSET];
    qarray[TEMPLIKEINDEP] = qarray[TEMPLIKEINDEP] - UTOTOFFSET;
    qarray[TEMPLIKEINDEP] = MAX(qarray[TEMPLIKEINDEP],lineartablelimits[whichtablevar][vartypearraylocal[TEMPLIKEINDEP]][0]); // lowest value of utotdiff = utot-utotoffset = utotoffset-utotoffset = 0.0 or minimum of table


  }
  else{

    // if here, then myanswers[EOSOFFSET]=U0, myanswers[EOSIN]=UIN, myanswers[EOSOUT]=UOUT, so apply inversion to "i/N"
    FTYPEEOS UTOT0=myanswers[EOSOFFSET];
    FTYPEEOS UTOTIN=myanswers[EOSIN];
    FTYPEEOS UTOTOUT=myanswers[EOSOUT];
    FTYPEEOS UTOT=qarray[TEMPLIKEINDEP];
    // e.g. from eos_extract.m: lutotdiff = log10( (utot - utotoffset)./(utotin - utotoffset))./log10( (utotout - utotoffset)./(utotin - utotoffset));

    qarray[TEMPLIKEINDEP] = log10( (UTOT-UTOT0)/(UTOTIN-UTOT0) )/log10( (UTOTOUT-UTOT0)/(UTOTIN-UTOT0) );

    // min of utotdegencut>=2 table format is lutotdiff="i/N"=0 and at this point qarray[TEMPLIKEINDEP]=lutotdiff
    qarray[TEMPLIKEINDEP] = MAX(qarray[TEMPLIKEINDEP],lineartablelimits[whichtablevar][vartypearraylocal[TEMPLIKEINDEP]][0]);

    // presume ok to let qarray[TEMPLIKEINDEP] be beyond "i/N" after which ideal gas EOS used

    //    dualfprintf(fail_file,"UTOT0=%21.15g UTOTIN=%21.15g UTOTOUT=%21.15g UTOT=%21.15g qarray[2]=%21.15g\n",UTOT0,UTOTIN,UTOTOUT,UTOT,qarray[TEMPLIKEINDEP]);

  }
  



}









// obtain whichd from whichindep
static int get_whichindep_fromwhichd(int whichd, int *whichindep)
{

  if(whichd==UTOTDIFF) *whichindep=UEOS;
  else if(whichd==PTOTDIFF) *whichindep=PEOS;
  else if(whichd==CHIDIFF) *whichindep=CHIEOS;
  else if(whichd==STOTDIFF) *whichindep=SEOS;
  else{
    dualfprintf(fail_file,"Undefined whichd=%d in get_whichindep_fromwhichd()\n",whichd);
    myexit(21937525);
  }

  return(0);

}







