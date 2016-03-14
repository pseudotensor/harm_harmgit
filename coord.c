
#include "decs.h"

/*! \file coord.c
    \brief User coordinates and other coordinates stuff
    
// this file contains all the coordinate dependent
// parts of the code, except the initial and boundary
// conditions 

// static variables with global scope to this file
// could make any or all of these true global if want to change them in, say, init.c

// Notes:
// 1) BCtype[] used below, but should presume BCtype[] could change even if storing X and V.  So far use of R0SING and DISKSURFACE is ok since not changing to/from those types.


// OPENMPMARK: Note that previously had many static (so global to this file) variables.  However, many were not thread safe since depended upon position or time.
// So need to ensure all statics here are set ONLY in set_coord_() type functions and NOT set in bl_coord() or dxdxp() or other functions that are called by multiple threads.
// Hence all these statics must be constant values in time and space.
// Note that JCM did have the variables set in the right location previously, just didn't have variables defined in right place since wasn't trying to be thread safe.
// Note that Intel Thread Checker didn't catch multiple thread use of myhslope.

*/


static void blcoord_singfixes(FTYPE *X, FTYPE *V);


///////////////////////
//
// These static variables can only be set in set_coord_parms() type functions, not bl_coord() or dxdxp().  if need to set inside bl_coord() or dxdxp(), then move variable there!
//
///////////////////////
// for defcoord==COMPLEX1TH
static FTYPE der0=9;
//static FTYPE Ri=8.0; // too unresolved for 256x128
static FTYPE Ri=20.0;

// for defcoord==COMPLEX2TH
static FTYPE x2trans;
static FTYPE m2,d2,c2,m3,b3,thetatores;

// for defcoord==JET1COORDS
static FTYPE h0,hf,rh0,myrout,dmyhslope1dr,dmyhslope2dx1,x1in,x1out;
static FTYPE npow;

// for defcoord=JET2COORDS
static FTYPE r1jet,njet,rpjet;

// for defcoord=JET3COORDS
static FTYPE r0jet,rsjet,Qjet;

// for SJETCOORDS
static FTYPE fracphi;
static FTYPE npow2;
static FTYPE cpow2;
static FTYPE rbr;
static FTYPE x1br;
static FTYPE fracdisk;
static FTYPE fracjet;
static FTYPE rsjet;
static FTYPE r0grid;
static FTYPE r0jet;
static FTYPE rjetend;
static FTYPE r0disk;
static FTYPE rdiskend;
static FTYPE jetnu;
static FTYPE x10;
static FTYPE x20;
#define USESJETLOGHOVERR 1
#if(USESJETLOGHOVERR)
static FTYPE torusrmax; // was extern and original located in init.sashatorus.c, but should be inverted as now is so extern is in init.sashatorus.c.
#endif
static FTYPE torusrmax_loc;



// for defcoord=JET6COORDS
static FTYPE ntheta,htheta,rsjet2,r0jet2,rsjet3,r0jet3; // and rs,r0
static FTYPE cpow3;

// for defcoord=BPTHIN1
static FTYPE bp_npow,bp_r1jet,bp_njet1,bp_njet,bp_r0jet,bp_rsjet,bp_Qjet, bp_ntheta,bp_htheta,bp_rsjet2,bp_r0jet2,bp_rsjet3,bp_r0jet3, bp_rs, bp_r0,bp_rsinner,bp_r0inner,bp_npow2,bp_cpow2,bp_rbr,bp_x1br, bp_h0; 

// for defcoord=PULSARCOORDS
static FTYPE hinner,houter;

// for defcoord=JET4COORDS
static FTYPE rs,r0;

// UNI2LOG similar to LOGSINTH with new features and simpler startx/dx and grid growth factor
static int Nstar;
static FTYPE Rstar,Afactor;

// for defcoord=JET5COORDS
static FTYPE AAAA,AAA,BBB,DDD,CCCC,Rj;
static FTYPE ii0;

// for defcoord=SJETCOORDS
static void vofx_sjetcoords( FTYPE *X, FTYPE *V );  //original coordinates
static void vofx_cylindrified( FTYPE *Xin, void (*vofx)(FTYPE*, FTYPE*), FTYPE *Vout ); //coordinate "cylindrifier"

// for defcoord=JET6COORDSTHIN
static FTYPE th_npow,th_r1jet,th_njet1,th_njet,th_r0jet,th_rsjet,th_Qjet, th_ntheta,th_htheta,th_rsjet2,th_r0jet2,th_rsjet3,th_r0jet3, th_rs, th_r0,th_npow2,th_cpow2,th_rbr,th_x1br, th_h0; 



/// can call when no dependencies
void set_coord_parms(int defcoordlocal)
{

  set_coord_parms_nodeps(defcoordlocal);
  set_coord_parms_deps(defcoordlocal);


}



/// Things to set that only depend upon defcoord and nothing else
/// NOTEMARK: By nothing else, that means things like hslope, R0, Rin, or anything else that user might set in init_grid()
/// Otherwise, even if set hslope here for example and then use hslope to set something else, hslope could change and that other parameter would be wrong and coordinates would be mismatched.
void set_coord_parms_nodeps(int defcoordlocal)
{


#if(USEOPENMP)
  if(omp_in_parallel()){
    dualfprintf(fail_file,"set_coord_parms_nodeps() called in parallel region\n");
    myexit(784653446);
  }
#endif



  // assumes R0, Rin, Rout, and hslope are so general that are set in init.c
  if (defcoordlocal == USERCOORD) {
    extern void set_coord_parms_nodeps_user(int defcoordlocal);
    set_coord_parms_nodeps_user(defcoordlocal);
  }
  else if (defcoordlocal == LOGRSINTH) {
  }
  else if (defcoordlocal == REBECCAGRID) {
  }
  else if (defcoordlocal == COMPLEX0TH) {
  }
  else if(defcoordlocal == UNIRSINTH || defcoordlocal == UNIRSINTH2){
  }
  else if (defcoordlocal == EQMIRROR) {
  }
  else if(defcoordlocal == COMPLEX1TH) {
  }
  else if(defcoordlocal == COMPLEX2TH) {
    x2trans=0.1; // user settable, must be same as below in dxdxp
  }
  else if (defcoordlocal == LOGRUNITH) { // uniform theta and log in radius
  }
  else if (defcoordlocal == JET1COORDS) {
    // optimal is npow=10 R0=-3
    npow=1.0;
    //R0=0.0;

    // must be same as in dxdxp()
    hf=2.0-0.22;
    rh0=40.0;

  }
  else if (defcoordlocal == JET2COORDS) {
    npow=1.0;

    // must be same as in dxdxp()
    if(1){
      r1jet=16.0;
      njet=0.3;
      rpjet=0.9;
    }
    else{
      r1jet=9.0;
      njet=0.3;
      rpjet=.9;
    }
  }
  else if (defcoordlocal == JET3COORDS) {
    npow=1.0;

    // must be same as in dxdxp()
    if(0){ // first attempt
      r1jet=2.8;
      njet=0.3;
      r0jet=7.0;
      rsjet=21.0;
      Qjet=1.7;
    }
    else if(0){ // chosen to resolve disk then resolve jet
      r1jet=2.8;
      njet=0.3;
      r0jet=20.0;
      rsjet=80.0;
      Qjet=1.8;
    }
    else if(1){
      r1jet=2.8;
      njet=0.3;
      r0jet=20.0;
      rsjet=80.0;
      Qjet=1.3; // chosen to help keep jet resolved even within disk region
    }
  }
  else if (defcoordlocal == SJETCOORDS) {
  }
  else if (defcoordlocal == JET6COORDS) {

    // see jet3coords_checknew.nb
    npow=1.0;

    /////////////////////
    // RADIAL GRID SETUP
    /////////////////////
    npow=1.0;  //don't change it, essentially equivalent to changing cpow2

    //radial hyperexponential grid

    //power exponent
    npow2=6.0; // WALD: 6.0->4.0

    cpow2=1.0; //exponent prefactor (the larger it is, the more hyperexponentiation is)
    //    cpow3=0.01;
    cpow3=1.0;
    //radius at which hyperexponentiation kicks in
    //    rbr = 1E3;
    //    rbr = 5E2; // WALD 5E2->5E7 // SUPER-EDD: 2E3
    rbr = 2E3; // SUPERMADNEW


    // must be same as in dxdxp()
    // GODMARK: Note njet here is overwritten by njet later, but could have been different values if setup variable names differently.
    if(0){ // first attempt
      r1jet=2.8;
      njet=0.3;
      r0jet=7.0;
      rsjet=21.0;
      Qjet=1.7;
    }
    else if(0){ // chosen to resolve disk then resolve jet
      r1jet=2.8;
      njet=0.3;
      r0jet=20.0;
      rsjet=80.0;
      Qjet=1.8;
    }
    else if(0){
      r1jet=2.8;
      njet=0.3;
      r0jet=15.0;
      rsjet=40.0;
      Qjet=1.3; // chosen to help keep jet resolved even within disk region
      //      Qjet=1.7; // chosen to help keep jet resolved even within disk region
    }
    else if(1){ // SUPERMADNEW
      r1jet=30.0;
      njet=0.7;
      r0jet=30.0;
      rsjet=40.0;
      Qjet=1.6;
    }

    // for switches from normal theta to ramesh theta
    if(0){
      rs=40.0; // shift
      r0=20.0; // divisor
      r0jet3=20.0; // divisor
      rsjet3=0.0; // subtractor
    }
    else{// SUPERMADNEW
      rs=40.0; // shift
      r0=40.0; // divisor
      r0jet3=40.0; // divisor
      rsjet3=0.0; // subtractor
    }

    // for theta1
    //    hslope=0.3 ; // resolve inner-radial region near equator

    // for theta2
    //h0=0.3; // inner-radial "hslope" for theta2
    h0=0.2; // inner-radial "hslope" for theta2
    //h0=0.1; // inner-radial "hslope" for theta2 // for thinner disks, change this.
    // GODMARK: Note that this overwrites above njet!
    // power \theta_j \propto r^{-njet}
    if(0){
      njet=1.0; // WALD: 1.0->0.0
    }
    else{//SUPERMADNEW
      // no change
    }


    // see fix_3dpoledtissue.nb
    if(0){
#if(0)
      ntheta=21.0;
      htheta=0.15;
      rsjet2=5.0;
      r0jet2=2.0;
#else
      ntheta=5.0;
      htheta=0.15;
      rsjet2=5.0;
      r0jet2=2.0;
#endif
    }
    else{ // SUPERMADNEW
      ntheta=5.0;
      if(a<0.4){
        htheta=0.15;
        rsjet2=5.0;
        r0jet2=2.0;
      }
      else{
        htheta=0.15;
        rsjet2=8.0;
        r0jet2=3.0;
      }
    }
  }
  else if (defcoordlocal == JET6COORDSTHIN) {

    // see jet3coords_checknew.nb
    th_npow=1.0;

    /////////////////////
    // RADIAL GRID SETUP
    /////////////////////
    th_npow=1.0;  //don't change it, essentially equivalent to changing cpow2

    //radial hyperexponential grid
    //    npow2=4.0; //power exponent
    th_npow2=4.0; //power exponent
    th_cpow2=1.0; //exponent prefactor (the larger it is, the more hyperexponentiation is)
    //    rbr = 1E3;  //radius at which hyperexponentiation kicks in
    th_rbr = 1E2;  //radius at which hyperexponentiation kicks in



    // must be same as in dxdxp()
    // GODMARK: Note njet here is overwritten by njet later, but could have been different values if setup variable names differently.
    if(0){ // first attempt
      th_r1jet=2.8;
      th_njet=0.3;
      th_r0jet=7.0;
      th_rsjet=21.0;
      th_Qjet=1.7;
    }
    else if(0){ // chosen to resolve disk then resolve jet
      th_r1jet=2.8;
      th_njet=0.3;
      th_r0jet=20.0;
      th_rsjet=80.0;
      th_Qjet=1.8;
    }
    else if(1){
      th_r1jet=2.8;
      th_njet=0.3;
      th_r0jet=15.0;
      th_rsjet=40.0;
      th_Qjet=2.0-0.05; // chosen to help keep jet resolved even within disk region
    }

    // for switches from normal theta to ramesh theta
    th_rs=60.0; // shift
    th_r0=20.0; // divisor
 
    // for theta1
    //    hslope=0.3 ; // resolve inner-radial region near equator
    th_r0jet3=20.0; // divisor
    th_rsjet3=0.0; // subtractor

    // for theta2
    th_h0=0.05; // inner-radial "hslope" for theta2
    //h0=0.1; // inner-radial "hslope" for theta2 // for thinner disks, change this.
    // GODMARK: Note that this overwrites above njet!
    th_njet=0.0; // power \theta_j \propto r^{-njet}


    // see fix_3dpoledtissue.nb
#if(0)
    th_ntheta=21.0;
    th_htheta=0.15;
    th_rsjet2=5.0;
    th_r0jet2=2.0;
#else
    th_ntheta=5.0;
    th_htheta=0.02;
    th_rsjet2=5.0;
    th_r0jet2=2.0;
#endif

  }
  else if (defcoordlocal == BPTHIN1) {

    // see jet3coords_checknew.nb
    bp_npow=1.0;

    /////////////////////
    // RADIAL GRID SETUP
    /////////////////////
    bp_npow=1.0;  //don't change it, essentially equivalent to changing cpow2

    //radial hyperexponential grid
    //    npow2=4.0; //power exponent
    bp_npow2=5.0; //MAVARANOTE must be odd now unless I add a sign explicitly to power component this contributes to sum in exponent //10.0; //5.0; // 10.0;    // MARKNOTE set to 10.0 before using BP values //power exponent
    bp_cpow2=1.0; //exponent prefactor (the larger it is, the more hyperexponentiation is)
    //    rbr = 1E3;  //radius at which hyperexponentiation kicks in
    bp_rbr = 200.0;  //radius at which hyperexponentiation kicks in



    // must be same as in dxdxp()
    // GODMARK: Note njet here is overwritten by njet later, but could have been different values if setup variable names differently.
    if(0){ // first attempt
      bp_r1jet=2.8;
      bp_njet1=0.3;
      bp_r0jet=7.0;
      bp_rsjet=21.0;
      bp_Qjet=1.7;
    }
    else if(0){ // chosen to resolve disk then resolve jet
      bp_r1jet=2.8;
      bp_njet1=0.3;
      bp_r0jet=20.0;
      bp_rsjet=80.0;
      bp_Qjet=1.8;
    }
    else if(1){
      bp_r1jet=2.8;
      bp_njet1=0.1; // MARKNOTE set to 0.3 before using BP values
      bp_r0jet=35.0;
      bp_rsjet=30.0;
      bp_Qjet=1.9;//-hslope; // chosen to help keep jet resolved even within disk region
    }

    // for switches from normal theta to ramesh theta
    bp_rs=200.0; // shift
    bp_r0=60.0; // divisor
 
    // for switches from innermost region of disk (inside horizon) to regular disk to increase timestep set by smallest vertical cell size
    bp_rsinner=2.0;//5.6;//4.0*Rin; //MAVARACHANGE changed from 4. and added the Rin so that the ratio bp_rsinner/r doesn't grow too large or too fast. maybe make that ratio **.5 to be even safer?
    bp_r0inner=1.33; //maybe 1.0 is too quick? not really same problem as outer radii I suppose since it just flattens off;

    // for theta1
    //    hslope=0.3 ; // resolve inner-radial region near equator
    bp_r0jet3=200.0; // divisor
    bp_rsjet3=0.0; //MAVARANOTE0.0; // subtractor

    // for theta2
    bp_h0=0.1; // inner-radial "hslope" for theta2
    // GODMARK: Note that this overwrites above njet!
    bp_njet=0.5;  // MARKNOTE set to 1.0 before using BP values // power \theta_j \propto r^{-njet}


    // see fix_3dpoledtissue.nb
#if(0)//HIGHRES // MAVARACHANGE I choose this because bp_ntheta 5 is less than the 0 used for the thin regime for the bp study. so, 5 note extreme enough.
    bp_ntheta=21.0; //13.0; // MAVARANOTE only use 21 for high res, use 15 for mid-res, non for low-res
    bp_htheta=0.45; // changed from .15 to be in line with my own additions for theta-flaring
    bp_rsjet2=5.0;
    bp_r0jet2=2.0;
#endif
#if(1) //MIDRES    note that lowres doesn't use polefix code
    bp_ntheta=15.0;
    bp_htheta=0.15;
    bp_rsjet2=5.0;
    bp_r0jet2=2.0;
#endif


  }
  else if (defcoordlocal == JET5COORDS) {
    // exp grid merged with exp-exp grid
    // parameters solved using hyperexp_gridnew.nb
    // Depends upon Rin, Rout, Rj, totalsize[1], and for below we used Rin=1.2, Rout=10^(10), Rj=200, TS1=256, and ii0,CC,Rj for radial arctan
    // Rj probably doesn't have to be the same thing, but for now it is since this is where grid changes alot and after which much lower resolution

#define JET5TOTALSIZE (256)
#define JET5RIN (1.2)
#define JET5ROUT (1.0E10)

    // checks
    if(totalsize[1]!=JET5TOTALSIZE){
      dualfprintf(fail_file,"Current version of JET5COORDS requires totalsize[1]=256\n");
      myexit(348766346);
    }

    // probably don't need to set Rout, but should in case user expects it
    Rin=JET5RIN;
    Rout=JET5ROUT;

    AAAA=1999999.0/10000000.0; // 0.1999999
    R0=AAAA; // effectively R0 is AAAA

    AAA=0.0413459589685779052930351140071389811117796472908765122327766247871075306910922595355681060060416677474341974954736231119642058094;
    //         \
    //      4691814961939384683077670140242359180355488020296128748293771170061841869426340268505040612342717691948841149166838622798123171255523798596 \
    //      5818680547438536476798449141070248313113472199351567812172169767872353912078416440520778774394376979127646837398673038048093220394452697865 \
    //      6270959185899435937659309684785579314134506823471357528404980034204759236451791458247221099942310718563615360919275492961171913096250029921 \
    //      4602829374292931981963109018352727784709700476977587800651816760158953266217495724111120779750291873137970716049214552447910902507302350084 \
    //      8366635098055440220680071046750709751515603946356328254368704921793804765397498748680467740257796496333857261456594221499641719912265168350 \
    //      4446757128527854665947791196023009608509876847107243924551220166086198337782848685465542023953471925358173682394567187686793903610577481426 \
    //      443229809938992518536215813767712488;

    BBB=-11.730265173318042629015657514515818843547237290015385234914265620733433049284290050282184485;
    //         \
    //      2224409622873658975481691490225985175297668826368377571878745806585146061321107420949292473465725043868512923396688827607703073903863397993 \
    //      2646469619777851403629789176127509226495862609935541778298333029485043643870849780090683673135501928153917249217151527805740858293313649525 \
    //      1577278440925240711864457756461606512526439401436022717818816632223921127138564528631741386676062985954517161211597241100034689634177608896 \
    //      2609005653930210512763437351860547338180817425783109360814416873897105930531521949202423966791692201348394578245024979983828552773565452313 \
    //      9427391374740569005171185204456464084827678985364372511780199067938494425840793867427758167352923180341529476553568176436448946131066359190 \
    //      5940692681151957860911707070494249851557079085990630079038126862565401104450140890585318688339870170091825436768403083368348398427545899255 \
    //      91796942327692774760083202686018064712337746785361430187509415618384926325;

    DDD=0.055717934049496306640561701541245682756775321774122925900528950779091605796182691888079948813285317249415181430716632182425357598150;
    //         \
    //      9878315509631578902141939586183299923885176395546296181175926730413185041822815279504848229461737912619728960820136753407636343528158363579 \
    //      3015197668283754748946354820375719585870499599501931419253237151291548191762172574322092394412311930543039905230389647665970462514296542459 \
    //      6284119291308204987329361603064178657563059514526318343174490238136571828705510375083642956288576409110703492926390263521348433706396587719 \
    //      4076567127699214794076208495287633220730415894184618889390695112123895752908866508866230483057346798669374527001159485781544411226064994838 \
    //      9334249843823760573294461519240571810219944708047592001849464422681300782954876398308748055099942704435601529920984600254710048152168797114 \
    //      1500470907761667436213106445853175905130588447192393215135706743963572902579082652204804896110286228839660989062139576983695311500579346405 \
    //      67924814417848136324354859648264319;

    // control radial arctan
    ii0=(FTYPE)(totalsize[1])*0.5;
    CCCC=5.0;
    Rj=200.0;

   

    // control \theta's arctan
    r1jet=2.8;
    njet=0.3;
    r0jet=20.0;
    rsjet=80.0;
    Qjet=1.3; // chosen to help keep jet resolved even within disk region

  }
  else if (defcoordlocal == PULSARCOORDS) {

    if(0){// pulsar in force free
      // pulsar_grid.nb for theta part and for the radial part:
      // see pulsar_gridnew.nb
      // for Rout=10^6 and R0=0.786*Rin Rin=4.84, npow=10 gives same dr/r as npow=1 R0=0.9*Rin at r=Rin
      npow=1.0;
      
      // must be same as in dxdxp()
      r0jet=5.0; // spread in radius over which hslope changes
      rsjet=18.0; // location of current sheet beginning for NS pulsar
    }
    else if(1){ // NS-pulsar in GRMHD
      npow=10.0;
      r0jet=5.0; // spread in radius over which hslope changes
      rsjet=15.0; // location of current sheet beginning for NS pulsar
    }

  }
  else if (defcoordlocal == UNIFORMCOORDS) {
    //uniform grid for Cartesian coordinates
  }
  else if (defcoordlocal == BILOGCYLCOORDS) {
    npow=10.0; // exponential rate
  }
  else if (defcoordlocal == RAMESHCOORDS || defcoordlocal == RAMESHCOORDS_HALFDISK) {
    //    myhslope=pow( (*r-rsjet)/r0jet , njet);
    npow=10.0;
    //npow=3.0;

    r0jet=2.0; // divisor
    njet=0.34; // power \theta_j \propto r^{-njet}
    //njet=1.0;
    rsjet=0.5; // subtractor
  }
  else if (defcoordlocal == JET4COORDS ) {
    // see net_jet_grid.nb
    // for small Rout, should use R0~0 (i.e. instead  use R0~-3) or hslope>1


    // this coordinate system uses:  R0 and npow for radius , hslope for theta1 , rsjet and r0 for switch and switchi , h0, rs, r0, njet for theta2 (as in JET3COORDS)

    // npow, R0, rs, r0, hslope, h0, r0jet, rsjet, njet

    // for radial grid
    npow=1.0;
    // npow=10.0;
    //npow=3.0;
    R0 = -3.0;
   
    // for switches
    rs=15.0;
    r0=25.0;
 
    // for theta1
    hslope=0.3 ; // resolve inner-radial region near equator
    // below 2 not used right now
    r0jet=15.0; // divisor
    rsjet=0.0; // subtractor

    // for theta2
    njet=0.34; // power \theta_j \propto r^{-njet}
  }
  else if (defcoordlocal == UNI2LOG) {

    if(1){
      Nstar = 20; // # of cells between Rin and Rstar (probably Rin=0)
      Afactor = 1000.0; // roughly Rout/Rstar
    }
    else{
      // GODMARK
      Nstar = 0; 
      Afactor = 1.01;
    }
    
  }
  else{
    dualfprintf(fail_file,"Shouldn't reach end of set_coord_parms: You set defcoordlocal=%d\n",defcoordlocal);
    myexit(1);
  }

}


/// stuff that depends upon ANYTHING external that might be set by user in init_grid() or by system in init_defgrid(), like R0, Rin, hslope, h_over_r, etc.
void set_coord_parms_deps(int defcoordlocal)
{

#if(USEOPENMP)
  if(omp_in_parallel()){
    dualfprintf(fail_file,"set_coord_parms_deps() called in parallel region\n");
    myexit(784653446);
  }
#endif

  // assumes R0, Rin, Rout, and hslope are so general that are set in init.c
  if (defcoordlocal == USERCOORD) {
    extern void set_coord_parms_deps_user(int defcoordlocal);
    set_coord_parms_deps_user(defcoordlocal);
  }
  else if (defcoordlocal == LOGRSINTH) {
  }
  else if (defcoordlocal == REBECCAGRID) {
  }
  else if (defcoordlocal == COMPLEX0TH) {
  }
  else if(defcoordlocal == UNIRSINTH || defcoordlocal == UNIRSINTH2){
  }
  else if (defcoordlocal == EQMIRROR) {
  }
  else if(defcoordlocal == COMPLEX1TH) {
  }
  else if(defcoordlocal == COMPLEX2TH) {
    thetatores=2.5*h_over_r;

    // fixed coefficients
    m2=(3.*(-2.*thetatores + M_PI))/(2.*x2trans) + (4.*thetatores)/(-1. + 2.*x2trans);
    d2=(2.*thetatores - M_PI + 2.*M_PI*x2trans)/(-2.*pow(x2trans,3.) + 4.*pow(x2trans,4.));
    c2=(6.*thetatores - 3.*M_PI + 6.*M_PI*x2trans)/(2.*pow(x2trans,2.) - 4.*pow(x2trans, 3.));
    m3=(2.*thetatores)/(1. - 2.*x2trans);
    b3=M_PI/2. + thetatores/(-1. + 2.*x2trans);
  }
  else if (defcoordlocal == LOGRUNITH) { // uniform theta and log in radius
  }
  else if (defcoordlocal == JET1COORDS) {
    h0=hslope;
    myrout=Rout;
    dmyhslope1dr = (hf-h0)/(myrout-rh0);
    dmyhslope2dx1=(hf-h0)/(x1out-x1in);
    x1in=log(Rin-R0);
    x1out=log(Rout-R0);
  }
  else if (defcoordlocal == JET2COORDS) {
  }
  else if (defcoordlocal == JET3COORDS) {
  }
  else if (defcoordlocal == SJETCOORDS) {   // AKMARK
    /////////////////////
    // RADIAL GRID SETUP
    /////////////////////
    npow=global_npow;  //don't change it, essentially equivalent to changing cpow2

    //radial hyperexponential grid
    npow2=global_npow2; //power exponent
    cpow2=global_cpow2; //exponent prefactor (the larger it is, the more hyperexponentiation is)
    rbr = global_rbr;  //radius at which hyperexponentiation kicks in
    x1br = log( rbr - R0 ) / npow;  //the corresponding X[1] value

    /////////////////////
    //ANGULAR GRID SETUP
    /////////////////////

    x10 = global_x10;
    x20 = global_x20;
     
    //transverse resolution fraction devoted to different components
    //(sum should be <1)
    fracdisk = global_fracdisk;
    fracjet = global_fracjet;

    jetnu = global_jetnu;  //the nu-parameter that determines jet shape

    //subtractor, controls the size of the last few cells close to axis:
    //if rsjet = 0, then no modification <- *** default for use with grid cylindrification
    //if rsjet ~ 0.5, the grid is nearly vertical rather than monopolar,
    //                which makes the timestep larger
    rsjet = global_rsjet; 

    //distance at which theta-resolution is *exactly* uniform in the jet grid -- want to have this at BH horizon;
    //otherwise, near-uniform near jet axis but less resolution (much) further from it
    //the larger r0grid, the larger the thickness of the jet 
    //to resolve
    r0grid = global_r0grid;    

    //distance at which jet part of the grid becomes monopolar
    //should be the same as r0disk to avoid cell crowding at the interface of jet and disk grids
    r0jet = global_r0jet;
    
    //distance after which the jet grid collimates according to the usual jet formula
    //the larger this distance, the wider is the jet region of the grid
    rjetend = global_rjetend;
    
    //distance at which disk part of the grid becomes monopolar
    //the larger r0disk, the larger the thickness of the disk 
    //to resolve
    r0disk = global_r0disk;

    //distance after which the disk grid collimates to merge with the jet grid
    //should be roughly outer edge of the disk
    rdiskend = global_rdiskend;
#if(USESJETLOGHOVERR)
    torusrmax_loc = torusrmax;
#else
    torusrmax_loc = 0.; //if not used, fill with dummy value
#endif


    /////////////////////
    //PHI GRID SETUP
    /////////////////////
    if( dofull2pi ) {
      fracphi = 1;
    }
    else {
      fracphi = global_fracphi;  //phi-extent measured in units of 2*PI, i.e. 0.25 means PI/2
    }   
  }
  else if (defcoordlocal == JET6COORDS) {
    x1br = log( rbr - R0 ) / npow;  //the corresponding X[1] value
  }
  else if (defcoordlocal == BPTHIN1) {

    /////////////////////
    //PHI GRID SETUP
    /////////////////////
    if( dofull2pi ) {
      fracphi = 1;
    }
    else {
      fracphi = global_fracphi;  //phi-extent measured in units of 2*PI, i.e. 0.25 means PI/2
    }   
    bp_x1br = log( bp_rbr - R0 ) / bp_npow;  //the corresponding X[1] value
  }
  else if (defcoordlocal == JET6COORDSTHIN) {
    th_x1br = log( th_rbr - R0 ) / th_npow;  //the corresponding X[1] value
  }
  else if (defcoordlocal == JET5COORDS) {
  }
  else if (defcoordlocal == PULSARCOORDS) {

    if(0){// pulsar in force free
      hinner=hslope; // hslope specifies inner hslope
      houter=hslope*0.05; // reduce by some arbitrary factor (currently 1/20)
    }
    else if(1){ // NS-pulsar in GRMHD
      hinner=1.9*hslope; // hslope specifies inner hslope
      //houter=hslope*0.001; // reduce by some arbitrary factor (currently 1/20)
      houter=hslope*1.5; // increase houter up to 2.0
    }

  }
  else if (defcoordlocal == UNIFORMCOORDS) {
  }
  else if (defcoordlocal == BILOGCYLCOORDS) {
  }
  else if (defcoordlocal == RAMESHCOORDS || defcoordlocal == RAMESHCOORDS_HALFDISK) {
  }
  else if (defcoordlocal == JET4COORDS ) {
    // for theta2
    h0=hslope; // inner-radial "hslope" for theta2
  }
  else if (defcoordlocal == UNI2LOG) {

    if(1){
      Rstar = 10.0*1E5/Lunit; // 10km
    }
    else{
      // GODMARK
      Rstar = Rin;
    }
    
    if(Nstar==0){
      if(fabs(Rstar-Rin)>SMALL){
        dualfprintf(fail_file,"If Nstar=0 then Rstar=Rin must be set\n");
        myexit(9279);
      }
    }
    
    trifprintf("Rstar = %26.20g Nstar=%d Afactor=%26.20g\n",Rstar,Nstar,Afactor);
  }
  else{
    dualfprintf(fail_file,"Shouldn't reach end of set_coord_parms: You set defcoordlocal=%d\n",defcoordlocal);
    myexit(1);
  }

}



/// write coordinate parameters
void write_coord_parms(int defcoordlocal)
{
  FILE *out;
  int dimen;

#if(USEOPENMP)
  if(omp_in_parallel()){
    dualfprintf(fail_file,"write_coord_parms_parms() called in parallel region\n");
    myexit(784653446);
  }
#endif


  if(myid==0){
    if((out=fopen("coordparms.dat","wt"))==NULL){
      dualfprintf(fail_file,"Couldn't write coordparms.dat file\n");
      myexit(1);
    }
    else{

      // same for all coords (notice no carraige return)
      fprintf(out,"%26.20g %26.20g %26.20g %26.20g %d ",R0,Rin,Rout,hslope,dofull2pi);

      if (defcoordlocal == USERCOORD) {
        extern void write_coord_parms_user(int defcoordlocal, FILE *out);
        write_coord_parms_user(defcoordlocal,out);
      }
      else if (defcoordlocal == LOGRSINTH) {
      }
      else if (defcoordlocal == REBECCAGRID) {
      }
      else if (defcoordlocal == COMPLEX0TH) {
      }
      else if(defcoordlocal == UNIRSINTH || defcoordlocal == UNIRSINTH2){
      }
      else if (defcoordlocal == EQMIRROR) {
      }
      else if(defcoordlocal == COMPLEX1TH) {
      }
      else if(defcoordlocal == COMPLEX2TH) {
        fprintf(out,"%26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g\n",x2trans,thetatores,m2,d2,c2,m3,b3,h_over_r);
      }
      else if (defcoordlocal == LOGRUNITH) { // uniform theta and log in radius
        DIMENLOOP(dimen) fprintf(out,"%26.20g ",Rin_array[dimen]);
        DIMENLOOP(dimen) fprintf(out,"%26.20g ",Rout_array[dimen]);
        fprintf(out,"\n");
      }
      else if (defcoordlocal == JET1COORDS) {
        fprintf(out,"%26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g\n",npow,h0,hf,rh0,myrout,dmyhslope1dr,dmyhslope2dx1,x1in,x1out);
      }
      else if (defcoordlocal == JET2COORDS) {
        fprintf(out,"%26.20g %26.20g %26.20g %26.20g\n",npow,r1jet,njet,rpjet);
      }
      else if (defcoordlocal == JET3COORDS) {
        fprintf(out,"%26.20g %26.20g %26.20g %26.20g %26.20g %26.20g\n",npow,r1jet,njet,r0jet,rsjet,Qjet);
      }
      else if (defcoordlocal == SJETCOORDS) {
        fprintf(out,"%26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g\n",npow,r1jet,njet,r0grid,r0jet,rjetend,rsjet,Qjet,fracphi,npow2,cpow2,rbr,x1br,fracdisk,fracjet,r0disk,rdiskend,torusrmax_loc,jetnu,x10,x20);
      }
      else if (defcoordlocal == JET6COORDS) {
        fprintf(out,"%26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g\n",npow,r1jet,njet,r0jet,rsjet,Qjet,ntheta,htheta,rsjet2,r0jet2,rsjet3,r0jet3,rs,r0,npow2,cpow2,rbr,x1br,cpow3);
      }
      else if (defcoordlocal == JET6COORDSTHIN) {
        fprintf(out,"%26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g\n",th_npow,th_r1jet,th_njet,th_r0jet,th_rsjet,th_Qjet,th_ntheta,th_htheta,th_rsjet2,th_r0jet2,th_rsjet3,th_r0jet3,th_rs,th_r0,th_npow2,th_cpow2,th_rbr,th_x1br,th_rbr,th_h0,th_njet1);
      }
      else if (defcoordlocal == BPTHIN1) {
        fprintf(out,"%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",bp_npow,bp_r1jet,bp_njet,bp_r0jet,bp_rsjet,bp_Qjet,bp_ntheta,bp_htheta,bp_rsjet2,bp_r0jet2,bp_rsjet3,bp_r0jet3,bp_rs,bp_r0,bp_rsinner,bp_r0inner,bp_npow2,bp_cpow2,bp_rbr,bp_x1br,fracphi);   // MARKTODO   add bp_h0? and add bp_njet1
      }
      else if (defcoordlocal == JET5COORDS) {
        fprintf(out,"%26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g\n",AAAA,AAA,BBB,DDD,ii0,CCCC,Rj);
        fprintf(out,"%26.20g %26.20g %26.20g %26.20g %26.20g\n",r1jet,njet,r0jet,rsjet,Qjet);
      }
      else if (defcoordlocal == PULSARCOORDS) {
        fprintf(out,"%26.20g %26.20g %26.20g %26.20g %26.20g\n",npow,hinner,houter,r0jet,rsjet);
      }
      else if (defcoordlocal == UNIFORMCOORDS) {
        //uniform grid for Cartesian coordinates
        DIMENLOOP(dimen) fprintf(out,"%26.20g ",Rin_array[dimen]);
        DIMENLOOP(dimen) fprintf(out,"%26.20g ",Rout_array[dimen]);
        fprintf(out,"\n");
      }
      else if (defcoordlocal == BILOGCYLCOORDS) {
        fprintf(out,"%26.20g\n",npow);
      }
      else if (defcoordlocal == RAMESHCOORDS || defcoordlocal == RAMESHCOORDS_HALFDISK) {
        fprintf(out,"%26.20g %26.20g %26.20g %26.20g\n",npow,r0jet,njet,rsjet);
      }
      else if (defcoordlocal == JET4COORDS) {
        // npow, rs, r0, h0, r0jet, njet, rsjet
        fprintf(out,"%26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g\n",npow,rs,r0,h0,r0jet,njet,rsjet);
      }
      else if (defcoordlocal == UNI2LOG) {
        fprintf(out,"%d %26.20g %26.20g\n",Nstar,Rstar,Afactor);
      }
      else{
        dualfprintf(fail_file,"Shouldn't reach end of write_coord_parms: You set defcoordlocal=%d\n",defcoordlocal);
        myexit(1);
      }

      fclose(out);
    }
  }
}


/// read coordinate parameters
void read_coord_parms(int defcoordlocal)
{
  FILE *in;
  FTYPE ftemp;
  int dimen;


#if(USEOPENMP)
  if(omp_in_parallel()){
    dualfprintf(fail_file,"read_coord_parms_parms() called in parallel region\n");
    myexit(784653446);
  }
#endif
  
  if(myid==0){
    in=fopen("coordparms.dat","rt");
    if(in==NULL){
      dualfprintf(fail_file,"Couldn't read coordparms.dat file.  I'll assume coded coordinates and let restart header overwrite any global restart parameters\n");
      set_coord_parms(defcoord);
    }
    else{
      // don't want to overwrite since restart file sets this
      //      fscanf(in,HEADER5IN,&ftemp,&ftemp,&ftemp,&ftemp,&ftemp);
      // NO: jon_interp.c requires read these in, so assume restart file has equal values to coordparms.dat file
      fscanf(in,HEADER4IN,&R0,&Rin,&Rout,&hslope);
      fscanf(in,"%d",&dofull2pi);



      if (defcoordlocal == USERCOORD) {
        extern void read_coord_parms_user(int defcoordlocal, FILE *in);
        read_coord_parms_user(defcoordlocal, in);
      }
      else if (defcoordlocal == LOGRSINTH) {
      }
      else if (defcoordlocal == REBECCAGRID) {
      }
      else if (defcoordlocal == COMPLEX0TH) {
      }
      else if(defcoordlocal == UNIRSINTH || defcoordlocal == UNIRSINTH2){
      }
      else if (defcoordlocal == EQMIRROR) {
      }
      else if(defcoordlocal == COMPLEX1TH) {
      }
      else if(defcoordlocal == COMPLEX2TH) {
        fscanf(in,HEADER8IN,&x2trans,&thetatores,&m2,&d2,&c2,&m3,&b3,&h_over_r);
      }
      else if (defcoordlocal == LOGRUNITH) { // uniform theta and log in radius
        DIMENLOOP(dimen) fscanf(in,HEADERONEIN,&Rin_array[dimen]);
        DIMENLOOP(dimen) fscanf(in,HEADERONEIN,&Rout_array[dimen]);
      }
      else if (defcoordlocal == JET1COORDS) {
        fscanf(in,HEADER9IN,&npow,&h0,&hf,&rh0,&myrout,&dmyhslope1dr,&dmyhslope2dx1,&x1in,&x1out);
      }
      else if (defcoordlocal == JET2COORDS) {
        fscanf(in,HEADER4IN,&npow,&r1jet,&njet,&rpjet);
      }
      else if (defcoordlocal == JET3COORDS) {
        fscanf(in,HEADER6IN,&npow,&r1jet,&njet,&r0jet,&rsjet,&Qjet);
      }
      else if (defcoordlocal == SJETCOORDS) {
        fscanf(in,HEADER9IN,&npow,&r1jet,&njet,&r0grid,&r0jet,&rjetend,&rsjet,&Qjet,&fracphi);
        fscanf(in,HEADER9IN,&npow2,&cpow2,&rbr,&x1br,&fracdisk,&fracjet,&r0disk,&rdiskend,&torusrmax_loc);
        fscanf(in,HEADER3IN,&jetnu,&x10,&x20);
      }
      else if (defcoordlocal == JET6COORDS) {
        fscanf(in,HEADER19IN,&npow,&r1jet,&njet,&r0jet,&rsjet,&Qjet,&ntheta,&htheta,&rsjet2,&r0jet2,&rsjet3,&r0jet3,&rs,&r0,&npow2,&cpow2,&rbr,&x1br,&cpow3);
      }
      else if (defcoordlocal == JET6COORDSTHIN) {
        fscanf(in,HEADER21IN,&th_npow,&th_r1jet,&th_njet,&th_r0jet,&th_rsjet,&th_Qjet,&th_ntheta,&th_htheta,&th_rsjet2,&th_r0jet2,&th_rsjet3,&th_r0jet3,&th_rs,&th_r0,&th_npow2,&th_cpow2,&th_rbr,&th_x1br,&th_rbr,&th_h0,&th_njet1);
      }
      else if (defcoordlocal == BPTHIN1) {
        fscanf(in,HEADER21IN,&bp_npow,&bp_r1jet,&bp_njet,&bp_r0jet,&bp_rsjet,&bp_Qjet,&bp_ntheta,&bp_htheta,&bp_rsjet2,&bp_r0jet2,&bp_rsjet3,&bp_r0jet3,&bp_rs,&bp_r0,&bp_rsinner,&bp_r0inner,&bp_npow2,&bp_cpow2,&bp_rbr,&bp_x1br,&fracphi);
      }
      else if (defcoordlocal == JET5COORDS) {
        fscanf(in,HEADER7IN,&AAAA,&AAA,&BBB,&DDD,&ii0,&CCCC,&Rj);
        fscanf(in,HEADER5IN,&r1jet,&njet,&r0jet,&rsjet,&Qjet);
      }
      else if (defcoordlocal == PULSARCOORDS) {
        fscanf(in,HEADER5IN,&npow,&hinner,&houter,&r0jet,&rsjet);
      }
      else if (defcoordlocal == UNIFORMCOORDS) {
        //uniform grid for Cartesian coordinates
        DIMENLOOP(dimen) fscanf(in,HEADERONEIN,&Rin_array[dimen]);
        DIMENLOOP(dimen) fscanf(in,HEADERONEIN,&Rout_array[dimen]);
      }
      else if (defcoordlocal == BILOGCYLCOORDS) {
        fscanf(in,HEADERONEIN,&npow);
      }
      else if (defcoordlocal == RAMESHCOORDS|| defcoordlocal == RAMESHCOORDS_HALFDISK) {
        fscanf(in,HEADER4IN,&npow,&r0jet,&njet,&rsjet);
      }
      else if (defcoordlocal == JET4COORDS) {
        fscanf(in,HEADER7IN,&npow,&rs,&r0,&h0,&r0jet,&njet,&rsjet);
        // npow, rs, r0, h0, r0jet, njet, rsjet
      }
      else if (defcoordlocal == UNI2LOG) {
        fscanf(in,"%d",&Nstar);
        fscanf(in,HEADERONEIN,&Rstar);
        fscanf(in,HEADERONEIN,&Afactor);
      }
      else{
        dualfprintf(fail_file,"Shouldn't reach end of read_coord_parms: You set defcoordlocal=%d\n",defcoordlocal);
        myexit(1);
      }

      fclose(in);
    }
  }

#if(USEMPI)
  // broadcast
  MPI_Bcast(&R0, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&Rin, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&Rout, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&hslope, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&dofull2pi, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);

  
  if (defcoordlocal == USERCOORD) {
    extern void read_coord_parms_mpi_user(int defcoordlocal);
    read_coord_parms_mpi_user(defcoordlocal);
  }
  else if (defcoordlocal == LOGRSINTH) {
  }
  else if (defcoordlocal == REBECCAGRID) {
  }
  else if (defcoordlocal == COMPLEX0TH) {
  }
  else if(defcoordlocal == UNIRSINTH || defcoordlocal == UNIRSINTH2){
  }
  else if (defcoordlocal == EQMIRROR) {
  }
  else if(defcoordlocal == COMPLEX1TH) {
  }
  else if(defcoordlocal == COMPLEX2TH) {
    MPI_Bcast(&x2trans, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&thetatores, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&m2, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&d2, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&c2, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&m3, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&b3, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    //  MPI_Bcast(&h_over_r, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD); // set by pre_init_specific_init() in init.c
  }
  else if (defcoordlocal == LOGRUNITH) { // uniform theta and log in radius
    DIMENLOOP(dimen) MPI_Bcast(&Rin_array[dimen], 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    DIMENLOOP(dimen) MPI_Bcast(&Rout_array[dimen], 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  }
  else if (defcoordlocal == JET1COORDS) {
    MPI_Bcast(&npow, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&h0, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&hf, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&rh0, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&myrout, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&dmyhslope1dr, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&dmyhslope2dx1, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&x1in, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&x1out, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  }
  else if (defcoordlocal == JET2COORDS) {
    MPI_Bcast(&npow, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&r1jet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&njet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&rpjet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  }
  else if (defcoordlocal == JET3COORDS) {
    MPI_Bcast(&npow, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&r1jet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&njet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&r0jet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&rsjet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&Qjet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  }
  else if (defcoordlocal == SJETCOORDS) {
    MPI_Bcast(&npow, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&r1jet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&njet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&r0grid, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&r0jet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&rjetend, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&rsjet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&Qjet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    //new params
    MPI_Bcast(&fracphi, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&npow2, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&cpow2, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&rbr, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&x1br, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&fracdisk, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&fracjet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&r0disk, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&rdiskend, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&torusrmax_loc, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&jetnu, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&x10, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&x20, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  }
  else if (defcoordlocal == JET6COORDS) {
    MPI_Bcast(&npow, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&r1jet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&njet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&r0jet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&rsjet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&Qjet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&ntheta, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&htheta, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&rsjet2, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&r0jet2, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&rsjet3, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&r0jet3, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&rs, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&r0, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&npow2, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&cpow2, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&rbr, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&x1br, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&cpow3, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  }
  else if (defcoordlocal == JET6COORDSTHIN) {
    MPI_Bcast(&th_npow, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&th_r1jet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&th_njet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&th_r0jet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&th_rsjet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&th_Qjet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&th_ntheta, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&th_htheta, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&th_rsjet2, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&th_r0jet2, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&th_rsjet3, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&th_r0jet3, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&th_rs, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&th_r0, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&th_npow2, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&th_cpow2, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&th_rbr, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&th_x1br, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&th_rbr, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&th_h0, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&th_njet1, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  }
  else if (defcoordlocal == BPTHIN1) {
    MPI_Bcast(&bp_npow, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&bp_r1jet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&bp_njet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&bp_r0jet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&bp_rsjet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&bp_Qjet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&bp_ntheta, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&bp_htheta, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&bp_rsjet2, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&bp_r0jet2, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&bp_rsjet3, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&bp_r0jet3, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&bp_rs, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&bp_r0, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&bp_rsinner, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&bp_r0inner, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&bp_npow2, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&bp_cpow2, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&bp_rbr, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&bp_x1br, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&fracphi, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  }
  else if (defcoordlocal == JET5COORDS) {
    MPI_Bcast(&AAAA, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&AAA, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&BBB, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&DDD, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&ii0, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&CCCC, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&Rj, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);

    MPI_Bcast(&r1jet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&njet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&r0jet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&rsjet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&Qjet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  }
  else if (defcoordlocal == PULSARCOORDS) {
    MPI_Bcast(&npow, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&hinner, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&houter, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&r0jet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&rsjet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  }
  else if (defcoordlocal == UNIFORMCOORDS) {
    //uniform grid for Cartesian coordinates
    DIMENLOOP(dimen) MPI_Bcast(&Rin_array[dimen], 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    DIMENLOOP(dimen) MPI_Bcast(&Rout_array[dimen], 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  }
  else if (defcoordlocal == BILOGCYLCOORDS) {
    MPI_Bcast(&npow, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  }
  else if (defcoordlocal == RAMESHCOORDS|| defcoordlocal == RAMESHCOORDS_HALFDISK) {
    MPI_Bcast(&npow, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&r0jet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&njet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&rsjet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  }
  else if (defcoordlocal == JET4COORDS) {
    MPI_Bcast(&npow, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&rs, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&r0, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&h0, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&r0jet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&njet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&rsjet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  }
  else if (defcoordlocal == UNI2LOG) {
    MPI_Bcast(&Nstar, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&Rstar, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&Afactor, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  }
  else{
    dualfprintf(fail_file,"Shouldn't reach end of read_coord_parms MPI stuff: You set defcoordlocal=%d\n",defcoordlocal);
    myexit(1);
  }
  
#endif

}



/// Returns boyer-lindquist coordinte of point
void bl_coord(FTYPE *X, FTYPE *V)
{
  extern FTYPE mysin(FTYPE th);
  FTYPE myx2;
  FTYPE mysign,ts1,fnstar,myNrat;
  FTYPE BB,CC;
  FTYPE myhslope1,myhslope2,myhslope;
  FTYPE flip1,flip2;
  FTYPE th0,th0toprint,th2,switch0,switch2,switchinner0,switchinner2,switchrad0,switchrad2,thetasign,x2temp;
  FTYPE r,dtheta2dx1,dtheta2dx2,dtheta0dx2,dtheta0dx1,dswitch0dr,dswitch2dr;
  //  FTYPE th0,th2,switch0,switch2;
  //  FTYPE r,dtheta2dx1,dtheta2dx2,dtheta0dx2,dtheta0dx1,dswitch0dr,dswitch2dr;
  FTYPE X0;
  // for defcoord=JET5COORDS
  FTYPE ii,logform,radialarctan,thetaarctan; // temp vars
  //for SJETCOORDS
  FTYPE theexp,theexp1,theexp2;


  // AKMARK: coordinates defined, in particular, phi wedge (e.g., V[3]=2.0*M_PI*X[3])
  // below will give correct dxdxp[1][1], etc.
  V[0]=X[0]; // assume time = X[0] means time now and negative means time in past and positive means future

  // in spherical polar coords: t=V[0] r=V[1] th=V[2] phi=V[3]


  if(defcoord == USERCOORD) {
    extern void blcoord_user(FTYPE *X, FTYPE *V);
    blcoord_user(X,V);
  }
  else if (defcoord == LOGRSINTH) {
#if(1)
    if(BCtype[X1DN]==R0SING){
      if(R0>=0.0){
        dualfprintf(fail_file,"With log grid and R0SING must have R0<0 instead of %26.20g\n",R0);
        myexit(8274);
      }
      X0 = log(-R0);
      if(X[1]>X0) V[1] = R0+exp(X[1]) ;
      else V[1] = -(R0+R0*R0*exp(-X[1])) ;
      //      dualfprintf(fail_file,"X0=%26.20g V[1]=%26.20g\n",X0,V[1]);
    }
    else{
      V[1] = R0+exp(X[1]) ;
      // if V[1]=r<0 here, the presume only where interpolation boundaries are, not evolved quantities, and not extending so far negative radius that reach beyond, e.g. light cylinder so that velocities will be undefined with simple extrapolation
    }
#else

    V[1] = R0+exp(X[1]) ;
#endif


    //    V[1] = Rin*exp(X[1]) ;
    //V[1] = Rin * exp(X[1]);
    if(X[2]<0.5){
      V[2] = M_PI * X[2] + ((1. - hslope) / 2.) * mysin(2. * M_PI * X[2]);
    }
    else{
      //      V[2] = 0.5*M_PI + M_PI * fabs(X[2]-0.5) + ((1. - hslope) / 2.) * (-mysin(2. * M_PI * (1.0-X[2])));
      V[2] = M_PI - (M_PI * (1.0-X[2])) + ((1. - hslope) / 2.) * (-mysin(2. * M_PI * (1.0-X[2])));
    }

    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];

  }
  else if (defcoord == REBECCAGRID) {
    //    V[1] = Rin*exp(X[1]) ;
    V[1] = R0+exp(X[1]) ;
    //V[1] = Rin * exp(X[1]);
    //  if(X[2]<0.5){
    //  V[2] = M_PI * X[2] + ((1. - hslope) / 2.) * mysin(2. * M_PI * X[2]);
    // }
    // else{
    //      V[2] = 0.5*M_PI + M_PI * fabs(X[2]-0.5) + ((1. - hslope) / 2.) * (-mysin(2. * M_PI * (1.0-X[2])));
    // V[2] = M_PI - (M_PI * (1.0-X[2])) + ((1. - hslope) / 2.) * (-mysin(2. * M_PI * (1.0-X[2])));
    // }
    V[2]=(hslope*((X[2]-0.5)/0.5) + (1-hslope)*pow((X[2]-0.5)/0.5, 7.0)+1.)*M_PI/2.;
    // default is uniform \phi grid
    V[3]=0.5*M_PI*X[3];

  }
  else if (defcoord == COMPLEX0TH) {
    V[1] = R0+Rin * exp(X[1] * log(Rout / Rin));
    V[2] =
      ((-49. * hslope + 60. * M_PI) * X[2]) / 12. +
      ((247. * hslope - 240. * M_PI) * pow(X[2],2)) / 12. +
      ((-83. * hslope + 80. * M_PI) * pow(X[2],3)) / 2. -
      (5. * (-25. * hslope + 24. * M_PI) * pow(X[2], 4)) / 3. +
      (2. * (-25. * hslope + 24. * M_PI) * pow(X[2], 5)) / 3.;
    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];
  }
  else if(defcoord==UNIRSINTH){
    V[1]=X[1];
    V[2] = M_PI* X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);
    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];
  }
  else if(defcoord==UNIRSINTH2){
    V[1]=Rin + (Rout-Rin)*X[1];
    V[2] = M_PI* X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);
    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];
  }
  else if (defcoord == EQMIRROR) {
    // MIRROR at equator, equator is outer theta edge
    V[1] = R0+exp(X[1]) ;
    V[2] = M_PI * X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);
    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];
  }
  else if(defcoord == COMPLEX1TH) {
    V[1] = R0+exp(X[1]) ;

    V[2] = (der0*X[2]*(-32.*pow(-1. + X[2],3.)*pow(X[2],2.)*(-1. + 2.*X[2]) - 
                       Ri*(-1. + X[2])*pow(-1. + 2.*X[2],3.)*
                       (-1. + 7.*(-1. + X[2])*X[2])) + 
            M_PI*Ri*pow(X[2],3.)*(70. + 
                                  3.*X[2]*(-105. + 2.*X[2]*(91. + 10.*X[2]*(-7. + 2.*X[2])))))/Ri;

    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];
  }
  else if(defcoord == COMPLEX2TH) {

    V[1] = R0+exp(X[1]) ;

    // now assign values
    if(X[2]<0.5){ myx2=X[2]; flip1=0.0; flip2=1.0;}
    else{ myx2=1.0-X[2]; flip1=M_PI; flip2=-1.0;}

    if(myx2<=x2trans){
      V[2] = flip1+flip2*(d2*pow(myx2,3.0)+c2*pow(myx2,2.0)+m2*myx2);
    }
    else{
      V[2] = flip1+flip2*(m3*myx2+b3);
    }
    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];

  }
  else if (defcoord == LOGRUNITH) { // uniform theta and log in radius
    V[1] = R0+exp(X[1]) ;
    V[2] = Rin_array[2] + X[2] * ( Rout_array[2] - Rin_array[2] );
    V[3] = Rin_array[3] + X[3] * ( Rout_array[3] - Rin_array[3] );
  }
  else if (defcoord == JET1COORDS) {
    V[1] = R0+exp(pow(X[1],npow)) ;

    myhslope1 = h0+dmyhslope1dr*((V[1])-rh0);
    myhslope2 = h0+(hf-h0)*(X[1]-x1in)/(x1out-x1in);

    myhslope=(myhslope1+myhslope2)*0.5;
    V[2] = M_PI * X[2] + ((1. - myhslope) / 2.) * sin(2. * M_PI * X[2]);
    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];
  }
  else if (defcoord == JET2COORDS) {
    V[1] = R0+exp(pow(X[1],npow)) ;

    myhslope=2.0-pow(V[1]/r1jet,njet*(-1.0+exp(1.0)/exp(V[1]+rpjet)));

    V[2] = M_PI * X[2] + ((1. - myhslope) / 2.) * sin(2. * M_PI * X[2]);
    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];
  }
  else if (defcoord == JET3COORDS) {
    V[1] = R0+exp(pow(X[1],npow)) ;

    myhslope=2.0-Qjet*pow(V[1]/r1jet,-njet*(0.5+1.0/M_PI*atan(V[1]/r0jet-rsjet/r0jet)));

    if(X[2]<0.5){
      V[2] = M_PI * X[2] + ((1. - myhslope) / 2.) * mysin(2. * M_PI * X[2]);
    }
    else{
      V[2] = M_PI - (M_PI * (1.0-X[2])) + ((1. - myhslope) / 2.) * (-mysin(2. * M_PI * (1.0-X[2])));
    }
    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];
  }
  else if (defcoord == SJETCOORDS) {
#if(0)
    //use original coordinates
    vofx_sjetcoords( X, V );
#else
    //apply cylindrification to original coordinates
    //[this internally calls vofx_sjetcoords()] 
    vofx_cylindrified( X, vofx_sjetcoords, V );
#endif
    
  }
  else if (defcoord == JET6COORDS) {

#if(0) // no change in exponentiation
    // JET3COORDS-like radial grid
    V[1] = R0+exp(pow(X[1],npow)) ;
#elif(1)

    theexp = npow*X[1];
    if( X[1] > x1br ) {
      theexp += cpow2 * pow(X[1]-x1br,npow2);
    }
    V[1] = R0+exp(theexp);


    //    FTYPE npowtrue,npowlarger=10.0;
    //    FTYPE npowrs=1E3;
    //    FTYPE npowr0=2E2;
    //    npowtrue = npow + (npowlarger-npow)*(0.5+1.0/M_PI*atan((V[1]-npowrs)/npowr0));
    //    V[1] = R0+exp(pow(X[1],npowtrue)) ;
#elif(0)
    // avoid jump in grid at rbr
    // determine switches
    FTYPE r0rbr=rbr/2.0;
    switch0 = 0.5+1.0/M_PI*atan((V[1]-rbr)/r0rbr); // 1 for outer r

    FTYPE V1 = R0+exp(npow*X[1]);
    FTYPE V2 = R0+exp(npow*X[1] + cpow2 * pow(cpow3*(X[1]-x1br*1.0),npow2));

    V[1] = V1*(1.0-switch0) + V2*switch0;

#endif



    FTYPE theta1,theta2,arctan2;


#if(0)
    // JET3COORDS-based:
    myhslope=2.0-Qjet*pow(V[1]/r1jet,-njet*(0.5+1.0/M_PI*atan(V[1]/r0jet-rsjet/r0jet)));
    theta1 = M_PI * X[2] + ((1. - myhslope) * 0.5) * mysin(2. * M_PI * X[2]);
#else
    // RAMESH BASED
    // myhslope here is h2 in MCAF paper
    //    // h0 here is h3 in MCAF paper
    //FTYPE njetvsr;
    //if(V[1]<rbr) njetvsr=njet;
    //    else njetvsr=njet/(V[1])*rbr;
    //else njetvsr=
    //njetvsr=njet;

    FTYPE localrbr=rbr; //500.0; // rbr;
//    FTYPE localrbrr0=MAX(100.0,localrbr/2.0);
    FTYPE localrbrr0=100.0;

    switch0 = 0.5+1.0/M_PI*atan((V[1]-localrbr)/localrbrr0); // large r
    switch2 = 1.0-switch0; // small r


    FTYPE myhslope1=h0 + pow( (V[1]-rsjet3)/r0jet3 , njet);
    FTYPE myhslope2=h0 + pow( (localrbr-rsjet3)/r0jet3 , njet);
    myhslope = myhslope1*switch2 + myhslope2*switch0;

    // determine theta2
    if(X[2]>1.0) myx2=2.0-X[2];
    else if(X[2]<0.0) myx2=-X[2];
    else myx2=X[2];

    th2 = 0.5*M_PI*(1.0 + atan(myhslope*(myx2-0.5))/atan(myhslope*0.5));

    if(X[2]>1.0) th2=2.0*M_PI-th2;
    else if(X[2]<0.0) th2=-th2;

    // determine theta0
    // JET3COORDS-based:
    myhslope1=2.0-Qjet*pow(V[1]/r1jet,-njet*(0.5+1.0/M_PI*atan(V[1]/r0jet-rsjet/r0jet)));
    myhslope2=2.0-Qjet*pow(localrbr/r1jet,-njet*(0.5+1.0/M_PI*atan(localrbr/r0jet-rsjet/r0jet)));
    myhslope = myhslope1*switch2 + myhslope2*switch0;

    if(0){
      // myhslope here is h0 in MCAF paper
      th0 = M_PI * X[2] + ((1. - myhslope) * 0.5) * mysin(2. * M_PI * X[2]);
    }
    else{  // SUPERMADNEW
      // poly grid
      FTYPE xi=((1. - myhslope) * 0.5);
      th0 = M_PI * .5 * (myhslope*(2.0*X[2]-1.0) + (1.0-myhslope)*pow(2.0*X[2]-1.0,3.0)+1.);
      //      dualfprintf(fail_file,"myhslope=%g th0=%g\n",myhslope,th0);
    }

    // determine switches (only function of radius and not x2 or theta)
    switch0 = 0.5+1.0/M_PI*atan((V[1]-rs)/r0); // switch in .nb file // switch0->1 as r->infinity
    switch2 = 1.0-switch0; // for inner radial region

    // this works because all functions are monotonic, so final result is monotonic.  Also, th(x2=1)=Pi and th(x2=0)=0 as required
    theta1 = th0*switch2 + th2*switch0; // th0 is activated for small V[1] and th2 is activated at large radii.  Notice that sum of switch2+switch0=1 so normalization correct.

#endif
    
    // fix_3dpoledtissue.nb based:
    theta2 = M_PI*0.5*(htheta*(2.0*X[2]-1.0)+(1.0-htheta)*pow(2.0*X[2]-1.0,ntheta)+1.0);

    // generate interpolation factor
    arctan2 = 0.5 + 1.0/M_PI*(atan( (V[1]-rsjet2)/r0jet2) );

    // now interpolate between them
    V[2] = theta2 + arctan2*(theta1-theta2);
    


    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];
  }
  else if (defcoord == BPTHIN1) {

#if(0) // no change in exponentiation
    // JET3COORDS-like radial grid
    V[1] = R0+exp(pow(X[1],bp_npow)) ;
#else

#if(0)
    theexp = bp_npow*X[1];
    if( X[1] > bp_x1br ) {
      theexp += bp_cpow2 * pow(X[1]-bp_x1br,bp_npow2);
    }
    V[1] = R0+exp(theexp);
#else
    /*if(X[1]<bp_x1br){
      switchrad0=0.0;
      switchrad2=1.0;
    }
    else
    {
	switchrad0=1.0;
	switchrad2=0.0;
	}*/
    switchrad0 = 0.5+1.0/M_PI*atan((X[1]-bp_x1br)*5.*10./dx[1]/totalsize[1]); //dx[1]/N1);//totalsize[1]); // switch in .nb file
    switchrad2 = 0.5-1.0/M_PI*atan((X[1]-bp_x1br)*5.*10./dx[1]/totalsize[1]); //dx[1]~.03151/N1);//totalsize[1]); // switchi in .nb file

    theexp1 = bp_npow*X[1];
    theexp2 = theexp1+bp_cpow2 * pow(X[1]-bp_x1br,bp_npow2);
    V[1] = (R0+exp(theexp1))*switchrad2 + (R0+exp(theexp2))*switchrad0; 
#endif

    //    FTYPE npowtrue,npowlarger=10.0;
    //    FTYPE npowrs=1E3;
    //    FTYPE npowr0=2E2;
    //    npowtrue = npow + (npowlarger-npow)*(0.5+1.0/M_PI*atan((V[1]-npowrs)/npowr0));
    //    V[1] = R0+exp(pow(X[1],npowtrue)) ;
#endif



    FTYPE theta1,theta2,arctan2;


#if(0)
    // JET3COORDS-based:
    myhslope=2.0-bp_Qjet*pow(V[1]/bp_r1jet,-bp_njet*(0.5+1.0/M_PI*atan(V[1]/bp_r0jet-bp_rsjet/bp_r0jet)));
    theta1 = M_PI * X[2] + ((1. - myhslope) * 0.5) * mysin(2. * M_PI * X[2]);
#else
    // RAMESH BASED
    // myhslope here is h2 in MCAF paper
    // h0 here is h3 in MCAF paper
    //if(V[1] > bp_rsjet3){
    myhslope=bp_h0 + pow( (V[1]-bp_rsjet3)/bp_r0jet3 , bp_njet);
    //}
    /*else {
    myhslope=bp_h0;
    }*/
    // determine theta2
    if(X[2]>1.0) myx2=2.0-X[2];
    else if(X[2]<0.0) myx2=-X[2];
    else myx2=X[2];

    th2 = M_PI*.5*(.2*(2.0*myx2-1.0) + (1.0-.2)*pow(2.0*myx2-1,3.0)+1.0);

    if(X[2]>1.0) th2=2.0*M_PI-th2;
    else if(X[2]<0.0) th2=-th2;

    // determine theta0
    // JET3COORDS-based:
    //    myhslope=2.0-bp_Qjet*pow(V[1]/bp_r1jet,-bp_njet1*(0.5+1.0/M_PI*atan(V[1]/bp_r0jet-bp_rsjet/bp_r0jet)));
    myhslope=hslope;
    // myhslope here is h0 in MCAF paper


    //    th0 = M_PI * .5 * (1. + (1.-((1. - myhslope) * 0.5))*(2.*X[2]-1.) + ((1. - myhslope) * 0.5)*pow(2.*X[2]-1.,9) ) ; // MARKTODODONE  switched to poly type from Noble+ 2010 on June 10, 2013
    FTYPE xi=((1. - myhslope) * 0.5);
    //  th0 = M_PI * .5 * (1. + (1.-xi)*(2.*X[2]-1.) + xi*pow(2.*X[2]-1.,9) ) ; // MARKTODODONE  switched to poly type from Noble+ 2010 on June 10, 2013
    //switchinner0 = 0.5+1.0/M_PI*atan((V[1]-bp_rsinner)/bp_r0inner); // switch in .nb file
    //switchinner2 = 0.5-1.0/M_PI*atan((V[1]-bp_rsinner)/bp_r0inner); // switchi in .nb file
    
    //th0 = M_PI * .5 * (.2*(2.0*X[2]-1.0) + (1.0-.2)*pow(2.0*X[2]-1.0,9.0)+1.) ;
    /*
    if(X[2]>=.5){
      thetasign=+1.0;
      x2temp=X[2];
	}
    else {
      thetasign=-1.0;
      x2temp=1.0-X[2];
    }
    */
    th0 = M_PI * .5 * (.11875*(1.+(bp_rsinner/V[1]))*(2.0*X[2]-1.0) +(1.0-.11875*(1.+(bp_rsinner/V[1])))*pow(2.0*X[2]-1.0,9.0)+1.) ; // .1096=.1425/(1+6/20) -- .11875 is .1425/(1+4/20) --- .17 is .2/(1.+4/15.)
    //    if(X[2]>=0.5 && (mycpupos[2]==ncpux2/2 && ncpux2>1 || ncpux2==1)) printf("at radius %21.15g and X[2] = %21.15g the diff is %21.15e\n",V[1],X[2],th0toprint);

    // th0 = M_PI * .5 * (.2*(2.0*X[2]-1.0) + (1.0-.2)*pow(2.0*X[2]-1.0,9.0)+1.) ;
    
    // determine switches (only function of radius and not x2 or theta)
    switch0 = 0.5+1.0/M_PI*atan((V[1]-bp_rs)/bp_r0); // switch in .nb file
    switch2 = 0.5-1.0/M_PI*atan((V[1]-bp_rs)/bp_r0); // switchi in .nb file

    // this works because all functions are monotonic, so final result is monotonic.  Also, th(x2=1)=Pi and th(x2=0)=0 as required
    theta1 = th0*switch2 + th2*switch0; // th0 is activated for small V[1] and th2 is activated at large radii.  Notice that sum of switch2+switch0=1 so normalization correct.

#endif

#if(1)    
    // fix_3dpoledtissue.nb based:
    theta2 = M_PI*0.5*(.11875*(1.+(bp_rsinner/V[1]))*(2.0*X[2]-1.0)+(1.0-.11875*(1.+(bp_rsinner/V[1])))*pow(2.0*X[2]-1.0,bp_ntheta)+1.0);

    // generate interpolation factor
    arctan2 = 0.5 + 1.0/M_PI*(atan( (V[1]-bp_rsjet2)/bp_r0jet2) ); // MAVARA: outside a certain radius this switches the v2 dependence from theta2 to theta1....known as interpolation. this interpolation fixes pole issue. previous involving switch0/2 involves more fundamental difference in vertical distribution.

    // now interpolate between them
    V[2] = theta2 + arctan2*(theta1-theta2);
#else
    V[2] = theta1;
#endif


    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];
  }
  else if (defcoord == JET5COORDS) {

    // radial grid
    ii=X[1]*((FTYPE)totalsize[1]); // assume X[1]=0..1 for active grid going from Rin to Rout

    // radial arctan
    radialarctan=(1.0/2.0) + (1.0/M_PI)*atan((ii-ii0)/CCCC);

    // merge of exp grid and exp-exp grid

    // first term is normal exp, second term is exp-exp
    logform = AAA*ii + exp(BBB+DDD*ii)*radialarctan;

    // final radius
    V[1] = AAAA+exp(logform);

    


    /////////////////
    // \theta grid

    // theta arctan
    thetaarctan=(1.0/2.0) + (1.0/M_PI)*atan((V[1]-rsjet)/r0jet);

    // h(r)
    myhslope=2.0-Qjet*pow(V[1]/r1jet,-njet*thetaarctan);

    // final theta
    if(X[2]<0.5){
      V[2] = M_PI * X[2] + ((1. - myhslope) / 2.) * mysin(2. * M_PI * X[2]);
    }
    else{
      V[2] = M_PI - (M_PI * (1.0-X[2])) + ((1. - myhslope) / 2.) * (-mysin(2. * M_PI * (1.0-X[2])));
    }

    //////////////////
    // \phi grid
    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];

  }
  else if (defcoord == JET6COORDSTHIN) {

#if(0) // no change in exponentiation
    // JET3COORDS-like radial grid
    V[1] = R0+exp(pow(X[1],th_npow)) ;
#else

    theexp = th_npow*X[1];
    if( X[1] > th_x1br ) {
      theexp += th_cpow2 * pow(X[1]-th_x1br,th_npow2);
    }
    V[1] = R0+exp(theexp);


    //    FTYPE npowtrue,npowlarger=10.0;
    //    FTYPE th_npowrs=1E3;
    //    FTYPE th_npowr0=2E2;
    //    npowtrue = th_npow + (npowlarger-th_npow)*(0.5+1.0/M_PI*atan((V[1]-th_npowrs)/th_npowr0));
    //    V[1] = R0+exp(pow(X[1],npowtrue)) ;
#endif



    FTYPE theta1,theta2,arctan2;


#if(0)
    // JET3COORDS-based:
    myhslope=2.0-th_Qjet*pow(V[1]/th_r1jet,-th_njet*(0.5+1.0/M_PI*atan(V[1]/th_r0jet-th_rsjet/th_r0jet)));
    theta1 = M_PI * X[2] + ((1. - myhslope) * 0.5) * mysin(2. * M_PI * X[2]);
#else
    // RAMESH BASED
    // myhslope here is h2 in MCAF paper
    //    // h0 here is h3 in MCAF paper
    //FTYPE njetvsr;
    //if(V[1]<th_rbr) njetvsr=th_njet;
    //    else njetvsr=th_njet/(V[1])*th_rbr;
    //else njetvsr=
    //njetvsr=th_njet;
    
    if(V[1]<th_rbr){
      myhslope=th_h0 + pow( (V[1]-th_rsjet3)/th_r0jet3 , th_njet);
    }
    else myhslope=th_h0 + pow( (th_rbr-th_rsjet3)/th_r0jet3 , th_njet);

    // determine theta2
    if(X[2]>1.0) myx2=2.0-X[2];
    else if(X[2]<0.0) myx2=-X[2];
    else myx2=X[2];

    th2 = 0.5*M_PI*(1.0 + atan(myhslope*(myx2-0.5))/atan(myhslope*0.5));

    if(X[2]>1.0) th2=2.0*M_PI-th2;
    else if(X[2]<0.0) th2=-th2;

    // determine theta0
    // JET3COORDS-based:
    if(V[1]<th_rbr){
      myhslope=2.0-th_Qjet*pow(V[1]/th_r1jet,-th_njet*(0.5+1.0/M_PI*atan(V[1]/th_r0jet-th_rsjet/th_r0jet)));
    }
    else myhslope=2.0-th_Qjet*pow(th_rbr/th_r1jet,-th_njet*(0.5+1.0/M_PI*atan(th_rbr/th_r0jet-th_rsjet/th_r0jet)));


    if(0){
      // myhslope here is h0 in MCAF paper
      th0 = M_PI * X[2] + ((1. - myhslope) * 0.5) * mysin(2. * M_PI * X[2]);
    }
    if(1){
      // poly grid
      FTYPE xi=((1. - myhslope) * 0.5);
      th0 = M_PI * .5 * (myhslope*(2.0*X[2]-1.0) + (1.0-myhslope)*pow(2.0*X[2]-1.0,9.0)+1.);
    }

    // determine switches (only function of radius and not x2 or theta)
    switch0 = 0.5+1.0/M_PI*atan((V[1]-th_rs)/th_r0); // switch in .nb file
    switch2 = 0.5-1.0/M_PI*atan((V[1]-th_rs)/th_r0); // switchi in .nb file

    // this works because all functions are monotonic, so final result is monotonic.  Also, th(x2=1)=Pi and th(x2=0)=0 as required
    theta1 = th0*switch2 + th2*switch0; // th0 is activated for small V[1] and th2 is activated at large radii.  Notice that sum of switch2+switch0=1 so normalization correct.

#endif
    
    // fix_3dpoledtissue.nb based:
    theta2 = M_PI*0.5*(th_htheta*(2.0*X[2]-1.0)+(1.0-th_htheta)*pow(2.0*X[2]-1.0,th_ntheta)+1.0);

    // generate interpolation factor
    arctan2 = 0.5 + 1.0/M_PI*(atan( (V[1]-th_rsjet2)/th_r0jet2) );

    // now interpolate between them
    V[2] = theta2 + arctan2*(theta1-theta2);
    


    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];
  }
  else if (defcoord == JET5COORDS) {

    // radial grid
    ii=X[1]*((FTYPE)totalsize[1]); // assume X[1]=0..1 for active grid going from Rin to Rout

    // radial arctan
    radialarctan=(1.0/2.0) + (1.0/M_PI)*atan((ii-ii0)/CCCC);

    // merge of exp grid and exp-exp grid

    // first term is normal exp, second term is exp-exp
    logform = AAA*ii + exp(BBB+DDD*ii)*radialarctan;

    // final radius
    V[1] = AAAA+exp(logform);

    


    /////////////////
    // \theta grid

    // theta arctan
    thetaarctan=(1.0/2.0) + (1.0/M_PI)*atan((V[1]-rsjet)/r0jet);

    // h(r)
    myhslope=2.0-Qjet*pow(V[1]/r1jet,-njet*thetaarctan);

    // final theta
    if(X[2]<0.5){
      V[2] = M_PI * X[2] + ((1. - myhslope) / 2.) * mysin(2. * M_PI * X[2]);
    }
    else{
      V[2] = M_PI - (M_PI * (1.0-X[2])) + ((1. - myhslope) / 2.) * (-mysin(2. * M_PI * (1.0-X[2])));
    }

    //////////////////
    // \phi grid
    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];

  }
  else if (defcoord == PULSARCOORDS) {
    V[1] = R0+exp(pow(X[1],npow)) ;

    myhslope=(0.5+1.0/M_PI*atan((V[1]-rsjet)/r0jet))*(houter-hinner)+hinner;

    if(X[2]<0.5){
      V[2] = M_PI * X[2] + ((1. - myhslope) / 2.) * mysin(2. * M_PI * X[2]);
    }
    else{
      V[2] = M_PI - (M_PI * (1.0-X[2])) + ((1. - myhslope) / 2.) * (-mysin(2. * M_PI * (1.0-X[2])));
    }
    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];
  }
  else if (defcoord == UNIFORMCOORDS) {
    //uniform grid for Cartesian coordinates
    V[1] = Rin_array[1] + X[1] * ( Rout_array[1] - Rin_array[1] );
    V[2] = Rin_array[2] + X[2] * ( Rout_array[2] - Rin_array[2] );
    V[3] = Rin_array[3] + X[3] * ( Rout_array[3] - Rin_array[3] );
  }
  else if (defcoord == BILOGCYLCOORDS) {
    // R : cylindrical radius, assumes X[1]=0..1
    // exponential grid
    V[1] = ((Rout-Rin)*exp(npow*X[1])+Rin*exp(npow)-Rout)/(exp(npow)-1.0);
    
    // z : cylindrical height, assumes X[2]=-1..1
    // bi-exponential grid
    // here the grid goes from Zin to Zout in a bi-log way, and X[2]=0 is Z=0
    if(X[2]>0.0) V[2] = ((Zout-0)*exp(npow*fabs(X[2])) + 0*exp(npow)-Zout)/(exp(npow)-1.0);
    else V[2] = ((Zin-0)*exp(npow*fabs(X[2])) + 0*exp(npow)-Zin)/(exp(npow)-1.0);
    
  } 
  else if (defcoord == RAMESHCOORDS|| defcoord == RAMESHCOORDS_HALFDISK) {
    V[1] = R0+exp(pow(X[1],npow)) ;

    myhslope=pow( (V[1]-rsjet)/r0jet , njet);

    if(X[2]>1.0) myx2=2.0-X[2];
    else if(X[2]<0.0) myx2=-X[2];
    else myx2=X[2];

    V[2] = 0.5*M_PI*(1.0 + atan(myhslope*(myx2-0.5))/atan(myhslope*0.5));

    if(X[2]>1.0) V[2]=2.0*M_PI-V[2];
    else if(X[2]<0.0) V[2]=-V[2];

    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];


  }
  else if (defcoord == JET4COORDS) {


#if(0)
    // combines RAMESHCOORDS with original simple SINTH grid

    if(BCtype[X1DN]==R0SING){
      if(R0>=0.0){
        dualfprintf(fail_file,"With log grid and R0SING must have R0<0 instead of %26.20g\n",R0);
        myexit(8274);
      }
      X0 = log(-R0);
      if(X[1]>X0) V[1] = R0+exp(X[1]) ;
      else V[1] = -(R0+R0*R0*exp(-X[1])) ;
      //      dualfprintf(fail_file,"X0=%26.20g V[1]=%26.20g\n",X0,V[1]);
    }
    else{
      V[1] = R0+exp(X[1]) ;
      // if V[1]=r<0 here, the presume only where interpolation boundaries are, not evolved quantities, and not extending so far negative radius that reach beyond, e.g. light cylinder so that velocities will be undefined with simple extrapolation
    }
#else
    // JET3COORDS-like radial grid
    V[1] = R0+exp(pow(X[1],npow)) ;
#endif



    myhslope=h0 + pow( (V[1]-rsjet)/r0jet , njet);

    // determine theta2
    if(X[2]>1.0) myx2=2.0-X[2];
    else if(X[2]<0.0) myx2=-X[2];
    else myx2=X[2];

    th2 = 0.5*M_PI*(1.0 + atan(myhslope*(myx2-0.5))/atan(myhslope*0.5));

    if(X[2]>1.0) th2=2.0*M_PI-th2;
    else if(X[2]<0.0) th2=-th2;

    // determine theta0
    myhslope=hslope;

    if(X[2]<0.5){
      th0 = M_PI * X[2] + ((1. - myhslope) / 2.) * mysin(2. * M_PI * X[2]);
    }
    else{
      th0 = M_PI - (M_PI * (1.0-X[2])) + ((1. - myhslope) / 2.) * (-mysin(2. * M_PI * (1.0-X[2])));
    }

    // determine switches (only function of radius and not x2 or theta)
    switch0 = 0.5+1.0/M_PI*atan((V[1]-rs)/r0); // switch in .nb file
    switch2 = 0.5-1.0/M_PI*atan((V[1]-rs)/r0); // switchi in .nb file

    FTYPE theta1,theta2,arctan2;

    // this works because all functions are monotonic, so final result is monotonic.  Also, th(x2=1)=Pi and th(x2=0)=0 as required
    theta1 = th0*switch2 + th2*switch0; // th0 is activated for small V[1] and th2 is activated at large radii.  Notice that sum of switch2+switch0=1 so normalization correct.

    // fix_3dpoledtissue.nb based:
    theta2 = M_PI*0.5*(htheta*(2.0*X[2]-1.0)+(1.0-htheta)*pow(2.0*X[2]-1.0,ntheta)+1.0);

    // generate interpolation factor
    arctan2 = 0.5 + 1.0/M_PI*(atan( (V[1]-rsjet2)/r0jet2) );

    // now interpolate between them
    V[2] = theta2 + arctan2*(theta1-theta2);




    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];


  }
  else if (defcoord == UNI2LOG) {


    if(BCtype[X1DN]==R0SING && X[1]<0.0){
      mysign=-1.0;
    }
    else{
      mysign=1.0;
    }
    ts1 = (FTYPE)totalsize[1];
    fnstar = ((FTYPE)Nstar);
    myNrat = ts1/(fnstar+SMALL);


    if( Nstar==0 || fabs(X[1])>=1.0/myNrat ){ // same as if(fabs(i)>=Nstar)
      // log grid (similar to UNIRSINTH, but now such that startx=0 and dx=1/totalsize[1] and starts at Rstar and arbitrary growth factor Afactor)
      V[1] = mysign*(Rstar + (Rout-Rstar)*(pow(Afactor, (fabs(X[1])*ts1-fnstar)/(ts1-fnstar))-1.0)/(Afactor-1.0));
    }
    else{
      // see UNI2LOGgrid.nb in mathematica
      //
      // uniform grid (sign is correct already with uniform grid)
      //      V[1] = Rin + (Rstar-Rin)*myNrat*X[1];  // purely uniform grid leads to big connection coefficient at i=Nstar-1
      // so use smoother transition
      // i/Nstar = x1*totalsize[1]/Nstar = x1*myNrat
      // sign is not correct for this grid, so fix it
      BB = ( (Rout-Rstar)*log(Afactor)/( (Afactor-1.0)*(myNrat-1.0) ) );
      CC = BB + (Rin-Rstar);
      V[1] = mysign*( Rstar + BB*(fabs(X[1])*myNrat-1.0) + CC*(fabs(X[1])*myNrat-1.0)*(fabs(X[1])*myNrat-1.0) );
    }

    // fully express r=0 coordinate singularity
    // V[1]\sim 1E-14 if don't do this, and then gdet!=0 at r=0 and then flux through r=0 causes MBH and $a$ to be computed wrong.
    if(BCtype[X1DN]==R0SING && fabs(V[1]-10.0*NUMEPSILON)<10.0*NUMEPSILON ){
      V[1] = 0.0; // so coordinate singularity is fully represented in metric, etc.
    }

    


    // x2-direction



    if(X[2]<0.5){
      V[2] = M_PI * X[2] + ((1. - hslope) / 2.) * mysin(2. * M_PI * X[2]);
    }
    else{
      //      V[2] = 0.5*M_PI + M_PI * fabs(X[2]-0.5) + ((1. - hslope) / 2.) * (-mysin(2. * M_PI * (1.0-X[2])));
      V[2] = M_PI - (M_PI * (1.0-X[2])) + ((1. - hslope) / 2.) * (-mysin(2. * M_PI * (1.0-X[2])));
    }

    // x3-direction


    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];

  }
  else{
    dualfprintf(fail_file,"Shouldn't reach end of bl_coord: defcoord=%d\n",defcoord);
    myexit(1);
  }








  blcoord_singfixes(X,V);


}

static void blcoord_singfixes(FTYPE *X, FTYPE *V)
{

  //////
  //
  // don't allow to be smaller to avoid singularity
  // noted this caused problems with jon_interp in calculating jacobian
  //
  /////
  if(POSDEFMETRIC){
    if(V[2]<0) V[2] = -V[2];
    if(V[2]>M_PI) V[2]=2.0*M_PI-V[2];

    if(V[3]<0) V[3] = -V[3];
    if(V[3]>2.0*M_PI) V[3]=4.0*M_PI-V[3];

  }


#if( COORDSINGFIXCYL )   //SUPERSASMARK fix the singularity for the cylinrical coordinates
  // NOTEMARK: just shifting (e.g.) i=0 cell up a bit, nothing else to do. Assume only 1 grid cell (in "i") is there within such tolerance of SINGSMALL
  if(fabs(V[1]-0.0)<SINGSMALL) V[1]=SINGSMALL;
#endif


  if(ISSPCMCOORDNATIVE(MCOORD) && COORDSINGFIX){
    // avoid polar axis if SPC.  Also used for CTSTAG approach so can evolve B2
#if(1)
    // So use X[2] only -- closer to using j itself that we don't have available.
    if(BCtype[X2DN]==POLARAXIS && fabs(startx[TH]-X[TH])<SINGSMALL) V[TH]=SINGSMALL;
    FTYPE endx2=startx[TH]+totalsize[TH]*dx[TH];
    if(BCtype[X2UP]==POLARAXIS && fabs(endx2-X[TH])<SINGSMALL) V[TH]=M_PI-SINGSMALL;
#endif
#if(0)
    // OK, but worry about large radii where \theta is small towards axis
    if (BCtype[X2DN]==POLARAXIS && fabs(V[TH]) < SINGSMALL) V[TH]+=SINGSMALL;
    else if (BCtype[X2UP]==POLARAXIS && fabs(M_PI-V[TH]) < SINGSMALL)  V[TH]-=SINGSMALL;
#endif
#if(0)
    // WRONG!
    if (BCtype[X2DN]==POLARAXIS && fabs(V[TH]) < SINGSMALL){
      if(V[TH]>=0) V[TH]=SINGSMALL;
      if(V[TH]<0) V[TH]=-SINGSMALL;
    }
    if (BCtype[X2UP]==POLARAXIS && fabs(M_PI-V[TH]) < SINGSMALL){
      if(V[TH]>=M_PI) V[TH]=M_PI+SINGSMALL;
      if(V[TH]<M_PI) V[TH]=M_PI-SINGSMALL;
    }
#endif

  }

}


/// special v(x) for sjet coordinates
void vofx_sjetcoords( FTYPE *X, FTYPE *V )
{
  //for SJETCOORDS
  FTYPE theexp;
  FTYPE Ftrgen( FTYPE x, FTYPE xa, FTYPE xb, FTYPE ya, FTYPE yb );
  FTYPE limlin( FTYPE x, FTYPE x0, FTYPE dx, FTYPE y0 );
  FTYPE minlin( FTYPE x, FTYPE x0, FTYPE dx, FTYPE y0 );
  FTYPE mins( FTYPE f1, FTYPE f2, FTYPE df );
  FTYPE maxs( FTYPE f1, FTYPE f2, FTYPE df );
  FTYPE thetaofx2(FTYPE x2, FTYPE ror0nu);
  FTYPE  fac, faker, ror0nu;
  FTYPE fakerdisk, fakerjet;
  FTYPE rbeforedisk, rinsidedisk, rinsidediskmax, rafterdisk;
    
#define DOIMPROVEJETCOORDS 1
#if(DOIMPROVEJETCOORDS)
  FTYPE ror0nudisk, ror0nujet, thetadisk, thetajet;
#endif

  V[0] = X[0];

  theexp = npow*X[1];

  if( X[1] > x1br ) {
    theexp += cpow2 * pow(X[1]-x1br,npow2);
  }
  V[1] = R0+exp(theexp);

#if(0) //JON's method
  myhslope=2.0-Qjet*pow(V[1]/r1jet,-njet*(0.5+1.0/M_PI*atan(V[1]/r0grid-rsjet/r0grid)));

  if(X[2]<0.5){
    V[2] = M_PI * X[2] + ((1. - myhslope) / 2.) * mysin(2. * M_PI * X[2]);
  }
  else{
    V[2] = M_PI - (M_PI * (1.0-X[2])) + ((1. - myhslope) / 2.) * (-mysin(2. * M_PI * (1.0-X[2])));
  }
#elif(1) //SASHA's
  fac = Ftrgen( fabs(X[2]), fracdisk, 1-fracjet, 0, 1 );

  //faker = fac*V[1] + (1 - fac)*limlin(V[1],r0disk,0.5*r0disk,r0disk)*minlin(V[1],rdiskend,0.5*rdiskend,r0disk)/r0disk - rsjet*Rin;
    
  rbeforedisk = mins( V[1], r0disk, 0.5*r0disk );
#if(USESJETLOGHOVERR)
  //rinsidedisk = 1 for r < torusrmax_loc and increases logarithmically while r <= rdiskend, after which it  
  //levels off to the value = rinsidediskmax
  rinsidedisk = pow( 1. + 0.5*log10(mins(maxs(1,V[1]/torusrmax_loc,0.5),rdiskend/torusrmax_loc,0.5*rdiskend/torusrmax_loc)), 2./jetnu );
  rinsidediskmax = pow( 1. + 0.5*log10(rdiskend/torusrmax_loc), 2./jetnu);
#else
  rinsidedisk = 1.;
  rinsidediskmax = 1.;
#endif
  rafterdisk = maxs( 1, 1 + (V[1]-rdiskend)*r0jet/(rjetend*r0disk*rinsidediskmax), 0.5*rdiskend*r0jet/(rjetend*r0disk*rinsidediskmax) );
    
  fakerdisk = rbeforedisk * rinsidedisk * rafterdisk;
    
  fakerjet = mins( V[1], r0jet, 0.5*r0jet ) * maxs( 1, V[1]/rjetend, 0.5 );
    
#if( DOIMPROVEJETCOORDS )
  ror0nudisk = pow( (fakerdisk - rsjet*Rin)/r0grid, jetnu/2 );
  ror0nujet = pow( (fakerjet - rsjet*Rin)/r0grid, jetnu/2 );
  thetadisk = thetaofx2( X[2], ror0nudisk );
  thetajet = thetaofx2( X[2], ror0nujet );
  V[2] = fac*thetajet + (1 - fac)*thetadisk; 
#else
  faker = fac*fakerjet + (1 - fac)*fakerdisk - rsjet*Rin;
  ror0nu = pow( faker/r0grid, jetnu/2 );
  V[2] = thetaofx2( X[2], ror0nu );
#endif
  
  
  

#else
  //if((1+X[2])/2.<0.5){
  //  V[2] = M_PI * (1+X[2])/2. + ((1. - hslope) / 2.) * mysin(2. * M_PI * (1+X[2])/2.);
  //}
  //else{
  //  //      V[2] = 0.5*M_PI + M_PI * fabs(X[2]-0.5) + ((1. - hslope) / 2.) * (-mysin(2. * M_PI * (1.0-X[2])));
  //  V[2] = M_PI - (M_PI * (1.0-(1+X[2])/2.)) + ((1. - hslope) / 2.) * (-mysin(2. * M_PI * (1.0-(1+X[2])/2.)));
  //}
  V[2] = M_PI_2l * (1.0+ X[2]); 
#endif

  // default is uniform \phi grid
  V[3]=2.0*M_PI*X[3];
}


/// theta(x2) special
FTYPE thetaofx2(FTYPE x2, FTYPE ror0nu)
{
  FTYPE theta;
  if( x2 < -0.5 ) {
    theta = 0       + atan( tan((x2+1)*M_PI_2l)/ror0nu );
  }
  else if( x2 >  0.5 ) {
    theta = M_PI    + atan( tan((x2-1)*M_PI_2l)/ror0nu );
  }
  else {
    theta = M_PI_2l + atan( tan(x2*M_PI_2l)*ror0nu );
  }
  return(theta);
}  


/// Jacobian for dx uniform per dx nonuniform (dx/dr / dx/dr')
/// i.e. Just take d(bl-coord)/d(ksp uniform coord)
/// e.g. dr/dx1 d\theta/dx2
///
/// take note of the ordering of indicies
/// dxdxp[j][k]=dxdxp[mu][nu]=(dx^\mu_{BL}/dx^\nu_{KSP uni})
///
/// should make this numerical like connection, then to conserve CPU, would need all over grid
void dxdxprim(FTYPE *X, FTYPE *V, FTYPE (*dxdxp)[NDIM])
{
  void dxdxp_numerical(FTYPE *X, FTYPE (*dxdxp)[NDIM]);
  void dxdxp_analytic(FTYPE *X, FTYPE *V, FTYPE (*dxdxp)[NDIM]);

  if(defcoord<ANALYTICSWITCH){ // then have analytic dxdxp
    dxdxp_analytic(X,V,dxdxp);
  }
  else{
    dxdxp_numerical(X,dxdxp);
  }

  if(ISSPCMCOORD(MCOORD)){
    // below is because \int_0^\theta sin\theta d\theta d\phi is approximated as locally instead of finite volume
    if(totalsize[2]==1 && FIXGDETSPC_WHEN_1DRADIAL){
      dxdxp[2][2] = 2.0/dx[2]; // so that d\theta = 2
    }

#if(1)
    // dxdxp[2][1] changes sign across axis if consistently decollimate/collimate grid.  Yet, this will cause gv311 [primecoords] to be asymmetric for any given theta,phi when theta goes negative or beyond \pi.
    // To have correct symmetries in true SPC (not just extended \theta meanginless domain) must change sign and use FLIPU3AXIS 0 , FLIPB3AXIS 0 , FLIPU2AXIS 0 , and FLIPB2AXIS 0 .
    // This is the only way (e.g.) uu1 will be properly symmetric at respective theta,\phi positions when including ghost zones.  EMF needs v1=uu1/u0, so symmetry requires this to be correct.
    // This forces symmetry on primitives across pole, but may be bad for interpolation unless flip sign when doing interpolation itself.
    FTYPE endx2=startx[TH]+totalsize[TH]*dx[TH];
    if(X[2]<startx[2] || X[2]>endx2){
      dxdxp[2][0]*=-1.0;
      dxdxp[2][1]*=-1.0;
      dxdxp[2][3]*=-1.0;
    }
#endif

  }

}





/// should make this numerical like connection, then to conserve CPU, would need all over grid
void dxdxp_analytic(FTYPE *X, FTYPE *V, FTYPE (*dxdxp)[NDIM])
{
  int j,k;
  extern FTYPE mycos(FTYPE th);
  extern FTYPE mysin(FTYPE th);
  FTYPE myhslope1,myhslope2,myhslope;
  FTYPE dmyhslopedr,dmyhslopedx1;
  FTYPE myx2;
  FTYPE flip1,flip2;
  FTYPE th0,th2,switch0,switch2;
  FTYPE r,dtheta2dx1,dtheta2dx2,dtheta0dx2,dtheta0dx1,dswitch0dr,dswitch2dr;
  FTYPE X0;
  // for defcoord=JET5COORDS
  FTYPE ii,logform,radialarctan,thetaarctan; // temp vars



  // default identity transformation
  DLOOP(j,k) dxdxp[j][k]=0.0;
  DLOOPA(j) dxdxp[j][j]=1.0;

  if (defcoord == USERCOORD) {
    extern void dxdxp_analytic_user(FTYPE *X, FTYPE *V, FTYPE (*dxdxp)[NDIM]);
    dxdxp_analytic_user(X,V,dxdxp);
  }
  else if (defcoord == LOGRSINTH) {
    dxdxp[1][1] = V[1]-R0;
    if(X[2]<0.5){
      dxdxp[2][2] = M_PI + (1. - hslope) * M_PI * mycos(2. * M_PI * X[2]);
    }
    else{
      dxdxp[2][2] = M_PI + (1. - hslope) * M_PI * mycos(2. * M_PI * (1.0-X[2]) );
    }
    dxdxp[3][3] = 2.0*M_PI;
  }
  else if (defcoord == REBECCAGRID) {
    dxdxp[1][1] = V[1]-R0;
    dxdxp[2][2]=M_PI/2.*( (hslope * X[2])/0.5 + (7. * (1-hslope))/(pow(0.5, 7.)) *pow(X[2]-0.5, 6.)) ;

    // if(X[2]<0.5){
    //  dxdxp[2][2] = M_PI + (1. - hslope) * M_PI * mycos(2. * M_PI * X[2]);
    // }
    // else{
    //  dxdxp[2][2] = M_PI + (1. - hslope) * M_PI * mycos(2. * M_PI * (1.0-X[2]) );
    //  }
    dxdxp[3][3] = 0.5*M_PI;
  }
  else if (defcoord == COMPLEX0TH) {
    dxdxp[1][1] = (V[1]-R0) * log(Rout / Rin);

    dxdxp[2][2] = (-49. * hslope + 60. * M_PI) / 12. +
      ((247. * hslope - 240. * M_PI) * X[2]) / 6. +
      (3. * (-83. * hslope + 80. * M_PI) * pow(X[2], 2)) / 2. -
      (20. * (-25. * hslope + 24. * M_PI) * pow(X[2], 3)) / 3. +
      (10. * (-25. * hslope + 24. * M_PI) * pow(X[2], 4)) / 3.;

    dxdxp[3][3] = 2.0*M_PI;

  } else if(defcoord == UNIRSINTH || defcoord == UNIRSINTH2){
    dxdxp[2][2] = (totalsize[2]==1) ? (2.0) : (M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2]));
    dxdxp[3][3] = 2.0*M_PI;
  }
  else if (defcoord == EQMIRROR) {
    dxdxp[1][1] = V[1]-R0;
    dxdxp[2][2] = M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2]);
    dxdxp[3][3] = 2.0*M_PI;
  } 
  else if (defcoord == COMPLEX1TH) {
    dxdxp[1][1] = V[1]-R0;
    dxdxp[2][2] = (210.*M_PI*Ri*pow(1. - 2.*X[2],2.)*pow(-1. + X[2],2.)*
                   pow(X[2],2.) + der0*
                   (-32.*pow(-1. + X[2],2.)*pow(X[2],2.)*
                    (3. + 14.*(-1. + X[2])*X[2]) - 
                    Ri*pow(1. - 2.*X[2],2.)*
                    (-1. + 2.*(-1. + X[2])*X[2]*(2. + 49.*(-1. + X[2])*X[2]))))/Ri;
    dxdxp[3][3] = 2.0*M_PI;
  } 
  else if(defcoord == COMPLEX2TH) {
    dxdxp[1][1] = V[1]-R0;

    // now assign values
    if(X[2]<0.5){ myx2=X[2]; flip1=0.0; flip2=1.0;}
    else{ myx2=1.0-X[2]; flip1=M_PI; flip2=-1.0;}

    if(myx2<=x2trans){
      dxdxp[2][2] = (3.0*d2*pow(myx2,2.0)+2.0*c2*pow(myx2,1.0)+m2);
    }
    else{
      dxdxp[2][2] = (m3);
    }


    dxdxp[3][3] = 2.0*M_PI;


  }
  else if (defcoord == LOGRUNITH) {
    dxdxp[1][1] = V[1]-R0;
    dxdxp[2][2] = Rout_array[2] - Rin_array[2];
    dxdxp[3][3] = Rout_array[3] - Rin_array[3];
  }
  else if (defcoord == JET1COORDS) {

    //dxdxp[1][1] = npow*(V[1]-R0)*pow(log(V[1]-R0),(npow-1.0)/npow);
    dxdxp[1][1] = npow*(V[1]-R0)*pow(X[1],npow-1.0);

    myhslope1 = h0+dmyhslope1dr*(V[1]-rh0);
    myhslope2 = h0+dmyhslope2dx1*(X[1]-x1in);

    dmyhslopedx1=0.5*(dmyhslope1dr*dxdxp[1][1]+dmyhslope2dx1);
    myhslope=0.5*(myhslope1+myhslope2);

    dxdxp[2][2] = M_PI + (1. - myhslope) * M_PI * cos(2. * M_PI * X[2]);
    // d\theta/dx1 not 0
    // d\theta/dx1 = (d\theta/dr)*(dr/dx1)  (is this generally true or just when V[1](x1)?
    dxdxp[2][1] = -0.5*dmyhslopedx1* sin(2. * M_PI * X[2]);

    dxdxp[3][3] = 2.0*M_PI;
  }
  else if (defcoord == JET2COORDS) {
    // drdx1
    dxdxp[1][1] = npow*(V[1]-R0)*pow(X[1],npow-1.0);


    myhslope=2.0-pow(V[1]/r1jet,njet*(-1.0+exp(1.0)/exp(V[1]+rpjet)));

    dmyhslopedr=-((pow(exp(1.0),-V[1] - rpjet)*myhslope*njet*(-exp(1.0) + pow(exp(1.0),V[1] + rpjet) + exp(1.0)*V[1]*log(V[1]/r1jet)))/V[1]);
    dmyhslopedx1=dmyhslopedr*dxdxp[1][1];

    dxdxp[2][2] = M_PI + (1. - myhslope) * M_PI * cos(2. * M_PI * X[2]);
    dxdxp[2][1] = -0.5*dmyhslopedx1* sin(2. * M_PI * X[2]);
    dxdxp[3][3] = 2.0*M_PI;
  }
  else if(defcoord == JET3COORDS){
    // drdx1
    dxdxp[1][1] = npow*(V[1]-R0)*pow(X[1],npow-1.0);


    myhslope=2.0-Qjet*pow(V[1]/r1jet,-njet*(0.5+1.0/M_PI*atan(V[1]/r0jet-rsjet/r0jet)));

    dmyhslopedr=-((Qjet*(-((njet*(0.5 + atan(V[1]/r0jet - rsjet/r0jet)/M_PI))/V[1]) - (njet*r0jet*log(V[1]/r1jet))/(M_PI*(pow(r0jet,2) + pow(V[1] - rsjet,2)))))/pow(V[1]/r1jet,njet*(0.5 + atan(V[1]/r0jet - rsjet/r0jet)/M_PI)));

    dmyhslopedx1=dmyhslopedr*dxdxp[1][1];

    if(X[2]<0.5){
      dxdxp[2][2] = M_PI + (1. - myhslope) * M_PI * mycos(2. * M_PI * X[2]);
      dxdxp[2][1] = -0.5*dmyhslopedx1* mysin(2. * M_PI * X[2]);
    }
    else{
      dxdxp[2][2] = M_PI + (1. - myhslope) * M_PI * mycos(2. * M_PI * (1.0-X[2]));
      dxdxp[2][1] = -0.5*dmyhslopedx1* (-mysin(2. * M_PI * (1.0-X[2])));
    }


    dxdxp[3][3] = 2.0*M_PI;

    

  }
  else if(defcoord == JET6COORDS){
    dualfprintf(fail_file,"Should not be computing JET6COORDS analytically\n");
    myexit(34698346);
    dxdxp[3][3] = 2.0*M_PI;    
  }
  else if(defcoord == JET6COORDSTHIN){
    dualfprintf(fail_file,"Should not be computing JET6COORDSTHIN analytically\n");
    myexit(34698346);
    dxdxp[3][3] = 2.0*M_PI;    
  }
  else if(defcoord == BPTHIN1){
    dualfprintf(fail_file,"Should not be computing BPTHIN1 analytically\n");
    myexit(34698346);
    dxdxp[3][3] = 2.0*M_PI;    
  }
  else if(defcoord == JET5COORDS){
    dualfprintf(fail_file,"Should not be computing JET5COORDS analytically\n");
    myexit(34698346);
    dxdxp[3][3] = 2.0*M_PI;
  }
  else if(defcoord == PULSARCOORDS){
    // drdx1
    dxdxp[1][1] = npow*(V[1]-R0)*pow(X[1],npow-1.0);

    myhslope=(0.5+1.0/M_PI*atan((V[1]-rsjet)/r0jet))*(houter-hinner)+hinner;
    dmyhslopedr=(houter-hinner)*r0jet/(M_PI*(r0jet*r0jet+(V[1]-rsjet)*(V[1]-rsjet)));
    dmyhslopedx1=dmyhslopedr*dxdxp[1][1];

    dxdxp[2][2] = M_PI + (1. - myhslope) * M_PI * cos(2. * M_PI * X[2]);
    dxdxp[2][1] = -0.5*dmyhslopedx1* sin(2. * M_PI * X[2]);



    dxdxp[3][3] = 2.0*M_PI;



    

  }
  else if (defcoord == UNIFORMCOORDS) {
    //uniform grid for Cartesian coordinates
    dxdxp[1][1] = Rout_array[1] - Rin_array[1];
    dxdxp[2][2] = Rout_array[2] - Rin_array[2];
    dxdxp[3][3] = Rout_array[3] - Rin_array[3];
    //    dualfprintf(fail_file,"COORD.c: dxdxp[1][1]=%26.20g\n",dxdxp[1][1]);
  }
  else if (defcoord == BILOGCYLCOORDS) {
    myexit(6666);
  }
  else if(defcoord==RAMESHCOORDS|| defcoord == RAMESHCOORDS_HALFDISK){


    //    V[1] = R0+exp(pow(X[1],npow)) ;
    //    myhslope=pow( (*r-rsjet)/r0jet , njet);
    //    V[2] = 0.5*M_PI*(1.0 + atan(myhslope*(X[2]-0.5))/atan(myhslope*0.5));


    // drdx1
    dxdxp[1][1] = npow*(V[1]-R0)*pow(X[1],npow-1.0);
    // drdx2 = 0


    myhslope=pow( (V[1]-rsjet)/r0jet , njet);
    dmyhslopedr=(njet/r0jet)*pow( (V[1]-rsjet)/r0jet , njet-1.0);

    if(!finite(dmyhslopedr)){
      dualfprintf(fail_file,"Problem with dmyhslopedr=%g\n",dmyhslopedr);
      dualfprintf(fail_file,"njet=%g r=%g rsjet=%g r0jet=%g\n",njet,V[1],rsjet,r0jet);
      myexit(1);
    }

    // dhslope/dx1
    dmyhslopedx1=dmyhslopedr*dxdxp[1][1];


    if(X[2]>1.0) myx2=2.0-X[2];
    else if(X[2]<0.0) myx2=-X[2];
    else myx2=X[2];


    // d\theta/dx2
    // (2*Pi*h(r))/(ArcTan(h(r)/2.)*(4 + Power(1 - 2*x2,2)*Power(h(r),2)))
    
    dxdxp[2][2] = (2.0*M_PI*myhslope)/(atan(myhslope*0.5)*(4.0 + pow(1.0 - 2.0*myx2,2.0)*pow(myhslope,2.0)));


    // d\theta/dr
    //(Pi*((-4*ArcTan((-0.5 + x2)*h(r)))/(4 + Power(h(r),2)) + 
    //       (4*(-1 + 2*x2)*ArcTan(h(r)/2.))/
    //        (4 + Power(1 - 2*x2,2)*Power(h(r),2)))*
    //     Derivative(1)(h)(r))/(4.*Power(ArcTan(h(r)/2.),2))


    // d\theta/dx1  = d\theta/dr dr/dx1
    dxdxp[2][1] =     (M_PI*dmyhslopedx1*(
                                          (-4.0*atan((-0.5 + myx2)*myhslope))/(4. + pow(myhslope,2.)) + 
                                          (4.*(-1. + 2.*myx2)*atan(myhslope/2.))/(4. + pow(1. - 2.*myx2,2.)*pow(myhslope,2.))
                                          )
                       )/(4.*pow(atan(myhslope/2.),2.));


    if(X[2]>1.0) dxdxp[2][1]*=-1.0;
    if(X[2]<0.0) dxdxp[2][1]*=-1.0;
    

    dxdxp[3][3] = 2.0*M_PI;



  }
  else if(defcoord == JET4COORDS){


    r=V[1];

    /////////////////////
    //
    // determine theta2
    myhslope=h0 + pow( (V[1]-rsjet)/r0jet , njet);

    if(X[2]>1.0) myx2=2.0-X[2];
    else if(X[2]<0.0) myx2=-X[2];
    else myx2=X[2];

    th2 = 0.5*M_PI*(1.0 + atan(myhslope*(myx2-0.5))/atan(myhslope*0.5));

    if(X[2]>1.0) th2=2.0*M_PI-th2;
    else if(X[2]<0.0) th2=-th2;

    //////////////////////////
    //
    // determine theta0
    myhslope=hslope;

    if(X[2]<0.5){
      th0 = M_PI * X[2] + ((1. - myhslope) / 2.) * mysin(2. * M_PI * X[2]);
    }
    else{
      th0 = M_PI - (M_PI * (1.0-X[2])) + ((1. - myhslope) / 2.) * (-mysin(2. * M_PI * (1.0-X[2])));
    }

    // determine switches (only function of radius and not x2 or theta)
    switch0 = 0.5+1.0/M_PI*atan((V[1]-rs)/r0); // switch in .nb file
    switch2 = 0.5-1.0/M_PI*atan((V[1]-rs)/r0); // switchi in .nb file


    // now get derivatives

    // drdx1
    dxdxp[1][1] = npow*(V[1]-R0)*pow(X[1],npow-1.0);
    // drdx2 = 0


    ////////////////////////////////////
    //
    // for theta0
    //
    /////////////////////////////////////
    if(X[2]<0.5){
      dtheta0dx2 = M_PI + (1. - hslope) * M_PI * mycos(2. * M_PI * X[2]);
    }
    else{
      dtheta0dx2 = M_PI + (1. - hslope) * M_PI * mycos(2. * M_PI * (1.0-X[2]) );
    }

    dtheta0dx1 = 0.0; // for this simple grid with fixed hslope


    ////////////////////////////////////
    //
    // for theta2
    //
    /////////////////////////////////////
    myhslope=pow( (V[1]-rsjet)/r0jet , njet);
    dmyhslopedr=(njet/r0jet)*pow( (V[1]-rsjet)/r0jet , njet-1.0);

    if(!finite(dmyhslopedr)){
      dualfprintf(fail_file,"Problem with dmyhslopedr=%g\n",dmyhslopedr);
      dualfprintf(fail_file,"njet=%g r=%g rsjet=%g r0jet=%g\n",njet,V[1],rsjet,r0jet);
      myexit(1);
    }

    // dhslope/dx1
    dmyhslopedx1=dmyhslopedr*dxdxp[1][1];


    if(X[2]>1.0) myx2=2.0-X[2];
    else if(X[2]<0.0) myx2=-X[2];
    else myx2=X[2];


    
    dtheta2dx2 = (2.0*M_PI*myhslope)/(atan(myhslope*0.5)*(4.0 + pow(1.0 - 2.0*myx2,2.0)*pow(myhslope,2.0)));

    // d\theta/dx1  = d\theta/dr dr/dx1
    dtheta2dx1 =     (M_PI*dmyhslopedx1*(
                                         (-4.0*atan((-0.5 + myx2)*myhslope))/(4. + pow(myhslope,2.)) + 
                                         (4.*(-1. + 2.*myx2)*atan(myhslope/2.))/(4. + pow(1. - 2.*myx2,2.)*pow(myhslope,2.))
                                         )
                      )/(4.*pow(atan(myhslope/2.),2.));

    if(X[2]>1.0) dtheta2dx1*=-1.0;
    if(X[2]<0.0) dtheta2dx1*=-1.0;


    // switches
    dswitch0dr=1.0/(M_PI*r0*(1.0 + (r - rs)*(r - rs)/(r0*r0)));
    dswitch2dr=-dswitch0dr;

    // actual final derivatives
    dxdxp[2][2] = dtheta0dx2*switch2 + dtheta2dx2*switch0;

    // assumes r doesn't depend on x2
    dxdxp[2][1] =  dtheta0dx1*switch2 + th0*dswitch2dr*dxdxp[1][1] + dtheta2dx1*switch0+th2*dswitch0dr*dxdxp[1][1] ;






    dxdxp[3][3] = 2.0*M_PI;



  }

  else{
    dualfprintf(fail_file,"Shouldn't reach end of dxdxp: defcoord=%d\n",defcoord);
    myexit(1);
  }




}



#define GAMMIEDERIVATIVE 0
#define NUMREC 1  // at present this is too slow since dxdxp is called many times.  Could setup permenant memory space, but kinda stupid

// Note that can't use NUMREC for connection if using numerical dxdxp.  Any other form of connection and then any dxdxp can be used.

// For example, if both connection and dxdxp are computed using GAMMIEDERIVATIVE, seems to work fine as a decent numerical approximation

// Also, one can use NUMREC for dxdxp if using analytic connection or GAMMIEDERIVATIVE connection.

//#define DXDERTYPE DIFFNUMREC 
#define DXDERTYPE DIFFGAMMIE

// see conn_func() for notes
#if((REALTYPE==DOUBLETYPE)||(REALTYPE==FLOATTYPE))
#define DXDELTA (MY1EM5)
//#define DXDELTA (dx[1])
#elif(REALTYPE==LONGDOUBLETYPE)
#define DXDELTA (MY1EM6)
#endif

// finite volume form (nice if one could express so exact cancellation occurs in equations, but some errors vanish so that cannot identify how close to continuous solution one is since here one solves particular discrete equations)
//#define GENDXDELTA(k) (k==TT ? DXDELTA*Diffx[k] : 0.5*dx[k])

// finite difference form (nice to have "analytical" solution for dxdxp (and connection) so errors result that can be used to indicate how close to continuous solution one is)
// Diffx[k] enters because DXDELTA is relative to total scale involved, which is Diffx[k], not always 1.0
#define GENDXDELTA(k) (DXDELTA*Diffx[k])

/// numerical derivative of V w.r.t. X
void dxdxp_numerical(FTYPE *X, FTYPE (*dxdxp)[NDIM])
{
  int j,k,l;
  FTYPE Xh[NDIM], Xl[NDIM];
  FTYPE Vh[NDIM],Vl[NDIM];
  FTYPE blcoordsimple(struct of_geom *ptrgeom, FTYPE*X, int i, int j);
  extern int dfridr(FTYPE (*func)(struct of_geom *,FTYPE*,int,int), struct of_geom *ptrgeom, FTYPE *X,int ii, int jj, int kk, FTYPE *ans);
  void donothing(FTYPE *temp);
  FTYPE temp;
  FTYPE dxmachine[NDIM];
  struct of_geom geom;
  struct of_geom *ptrgeom;
  int failreturn;


  // setup dummy grid location since dxdxp doesn't need to know if on grid since don't store dxdxp (needed for dfridr())
  ptrgeom=&geom;
  ptrgeom->i=0;
  ptrgeom->j=0;
  ptrgeom->k=0;
  ptrgeom->p=NOWHERE;


  //  for(k=0;k<NDIM;k++) for(j=0;j<NDIM;j++){
  DLOOP(j,k){
    failreturn=0;
    if(DXDERTYPE==DIFFNUMREC){
      failreturn=dfridr(blcoordsimple,ptrgeom,X,0,j,k,&(dxdxp[j][k]));
    }
    if(DXDERTYPE==DIFFGAMMIE || failreturn==1){


      // I setup X and V relationship for time and phi to be correct now
      // Was usind dxdxp[3][3]=1 when V[3]=2.0*M_PI*X[3], so that was incorrect -- a bug
      /*
        if((j==TT)||(k==TT)){
        // assume no transformation of time coordinate and no mixing of t-coordinate with other coordinates (except what already in metric)
        if(j!=k) dxdxp[j][k]=0.0;
        else dxdxp[j][k]=1.0;
        }
        else if((j==PH)||(k==PH)){
        // assume no transformation of phi coordinate and no mixing of phi coordinate with other coordinates (at least no additional to existing metric)
        if(j!=k) dxdxp[j][k]=0.0;
        else dxdxp[j][k]=1.0;
        }
        else{
      */
      for(l=0;l<NDIM;l++){
        Xl[l]=Xh[l]=X[l]; // location of derivative
        temp = X[l]-GENDXDELTA(l);
        donothing(&temp);
        X[l] = temp+GENDXDELTA(l);
        temp = X[l]+GENDXDELTA(l);
        donothing(&temp);
        dxmachine[l] = temp-X[l];
      }
   
      //   Xh[k]+=dxmachine[k]; // shift up
      //   Xl[k]-=dxmachine[k]; // shift down

      Xh[k]+=GENDXDELTA(k); // shift up
      Xl[k]-=GENDXDELTA(k); // shift down

      //   dualfprintf(fail_file,"k=%d del=%g\n",k,GENDXDELTA(k));

      // below 2 lines redundant because gets both coordinates, but ok
      bl_coord(Xh, Vh);
      bl_coord(Xl, Vl);
      dxdxp[j][k] = (Vh[j] - Vl[j]) / (Xh[k] - Xl[k]);

      // GODMARK: unless N is a power of 2, X causes V to not be machine representable
      // GODMARK: Also, not only Xh-Xl, but each Xl and Xh must be machine representable

      // So even for a uniform grid dxdxp can vary near machine level
      //      dualfprintf(fail_file,"Vh=%26.20g Vl=%26.20g Xh=%2.15g Xl=%26.20g DX=%26.20g\n",Vh[j],Vl[j],Xh[k],Xl[k],GENDXDELTA(k));
      //      dualfprintf(fail_file,"(Vh[%d] - Vl[%d])=%26.20g (Xh[%d] - Xl[%d])=%26.20g\n",j,j,(Vh[j] - Vl[j]),k,k,(Xh[k] - Xl[k]));
      // }

      if(j==k && fabs(dxdxp[j][k])<NUMEPSILON){
        dualfprintf(fail_file,"dxdxp[%d][%d]=%g is too small.  Ensure SINGSMALL=%g < %g\n",j,k,dxdxp[j][k],SINGSMALL,(Xh[k] - Xl[k]));
      }

    }
  }
}


/// was using volatile, but not thread safe
void donothing(FTYPE *temp)
{
  *temp=*temp;

}


/// get V(X)
FTYPE blcoordsimple(struct of_geom *ptrgeom, FTYPE*X, int i, int j) // i not used
{
  FTYPE V[NDIM];
  
  // "dummy linear relationship" is right way and now setup when calling bl_coord()
  //  if((j==TT)||(j==PH)) return(X[j]); // dummy linear relationship
  //  else{
  bl_coord(X, V);
  return(V[j]);
  //  }

}



///////////////////////////////////////////////////////////////////
// 
// Below set X uniform grid -- usually doesn't change.
// Can usually force startx[1]=startx[2]=0. and dx[1]=1./N1 dx[2]=1./N2
// 
///////////////////////////////////////////////////////////////////


/// some grid location, dxs
/// could find this by root finding.  Needed if no obvious bounds
/// alternatively, could always define grid so x1=0..1 and x2=0..1 (likely more reasonable!)
void set_points()
{
  int jj;

  //for SJETCOORDS
  FTYPE x1max0, x1max,dxmax;
  int iter;
  const FTYPE RELACC = NUMEPSILON*100.0;
  const int ITERMAX = 100;

  // below 2 things not used since we set X[0] directly using t
  startx[0]=0;
  dx[0]=dt;

  if (defcoord == USERCOORD) {
    extern void set_points_user(void);
    set_points_user();
  }
  else if (defcoord == LOGRSINTH) {
    startx[1] = log(Rin-R0);
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = log((Rout-R0)/(Rin-R0)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 1.0/totalsize[3];
  }
  else if (defcoord == REBECCAGRID) {
    startx[1] = log(Rin-R0);
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = log((Rout-R0)/(Rin-R0)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 1.0/totalsize[3];
  }
  else if (defcoord == COMPLEX0TH) {
    startx[1] = 0.;
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = 1. / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 1.0/totalsize[3];
  } else if (defcoord == UNIRSINTH) {
    startx[1] = Rin;
    //    startx[2] = 0.324;
    startx[2] = 0.0;
    startx[3] = 0.;
    dx[1] = (Rout-Rin) / totalsize[1];
    //    dx[2] = (1.0-2*.324) / totalsize[2];
    dx[2] = 1.0 / totalsize[2];
    dx[3] = 1.0/totalsize[3];
  } else if (defcoord == UNIRSINTH2) {
    startx[1] = 0.0;
    startx[2] = 0.0;
    startx[3] = 0.0;
    dx[1] = 1.0 / totalsize[1];
    dx[2] = 1.0 / totalsize[2];
    dx[3] = 1.0 / totalsize[3];
  }
  else if (defcoord == EQMIRROR) {
    startx[1] = log(Rin-R0);
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = log((Rout-R0)/(Rin-R0)) / totalsize[1];
    dx[2] = 0.5 / totalsize[2];
    dx[3] = 1.0/totalsize[3];
  }
  else if (defcoord == COMPLEX1TH) {
    startx[1] = log(Rin-R0);
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = log((Rout-R0)/(Rin-R0)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 1.0/totalsize[3];
  }
  else if (defcoord == COMPLEX2TH) {
    startx[1] = log(Rin-R0);
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = log((Rout-R0)/(Rin-R0)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 1.0/totalsize[3];
  }
  else if (defcoord == LOGRUNITH) {
    startx[1] = log(Rin-R0);
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = log((Rout-R0)/(Rin-R0)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 1.0/totalsize[3];
  } 
  else if (defcoord == JET1COORDS) {
    startx[1] = pow(log(Rin-R0),1.0/npow);
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = (pow(log(Rout-R0),1.0/npow)-pow(log(Rin-R0),1.0/npow)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 1.0/totalsize[3];
  } 
  else if (defcoord == JET2COORDS) {
    startx[1] = pow(log(Rin-R0),1.0/npow);
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = (pow(log(Rout-R0),1.0/npow)-pow(log(Rin-R0),1.0/npow)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 1.0/totalsize[3];
  } 
  else if (defcoord == JET3COORDS) {
    startx[1] = pow(log(Rin-R0),1.0/npow);
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = (pow(log(Rout-R0),1.0/npow)-pow(log(Rin-R0),1.0/npow)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 1.0/totalsize[3];
  } 
  else if (defcoord == SJETCOORDS) {
    startx[1] = log(Rin-R0)/npow;

    if( Rout < rbr ) {
      x1max = log(Rout-R0)/npow;
    }
    else {
      x1max0 = 1;
      x1max = 2;

      //find the root via iterations
      for( iter = 0; iter < ITERMAX; iter++ ) {
        if( fabs((x1max - x1max0)/x1max) < RELACC ) {
          break;
        }
        x1max0 = x1max;
        dxmax= (pow( (log(Rout-R0) - npow*x1max0)/cpow2, 1./npow2 ) + x1br) - x1max0;

        // need a slight damping factor
        FTYPE dampingfactor=0.5;
        x1max = x1max0 + dampingfactor*dxmax;
      }

      if( iter == ITERMAX ) {
        trifprintf( "Error: iteration procedure for finding x1max has not converged: x1max = %g, dx1max/x1max = %g, iter = %d\n",
                    x1max, (x1max-x1max0)/x1max, iter );
        exit(1);
      }
      else {
        trifprintf( "x1max = %g (dx1max/x1max = %g, itno = %d)\n", x1max, (x1max-x1max0)/x1max, iter );
      }
    }

    startx[2] = -1.;
    startx[3] = 0.;
    dx[1] = ( x1max - startx[1] ) /totalsize[1];
    dx[2] = 2. / totalsize[2];
    dx[3] = fracphi/totalsize[3];
  } 
  else if (defcoord == JET6COORDS) { // same as JET3COORDS since radial grid same
    startx[1] = pow(log(Rin-R0),1.0/npow);
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = (pow(log(Rout-R0),1.0/npow)-pow(log(Rin-R0),1.0/npow)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 1.0/totalsize[3];

#if(1)
    startx[1] = log(Rin-R0)/npow;

    trifprintf( "ITERATIVE dx1: Rout=%26.20g R0=%26.20g npow=%26.20g cpow2=%26.20g cpow3=%26.20g npow2=%26.20g x1br=%26.20g rbr=%26.20g\n",Rout,R0,npow,cpow2,cpow3,npow2,x1br,rbr);

    if( Rout < rbr ) {
      x1max = log(Rout-R0)/npow;
    }
    else {
      x1max0 = 1.1*x1br;
      x1max = 1.2*x1br;

      //find the root via iterations
      for( iter = 0; iter < ITERMAX; iter++ ) {

        // trifprintf( "iter=%d x1max=%26.20g x2max0=%26.20g\n",iter,x1max0,x1max);

        if( fabs((x1max - x1max0)/x1max) < RELACC ) {
          break;
        }
        x1max0 = x1max;

        if(1){
          dxmax= (pow( (log(Rout-R0) - npow*x1max0)/cpow2, 1./npow2 ) + x1br*1.0) - x1max0;
        }
        else{
          // f-f0 = (x-x0)*dfdx -> if f=Rout -> x = (Rout-f0)/dfdx+x0
          
          FTYPE dVdx1=(npow + cpow2*npow2*cpow3*pow(cpow3*(x1max0-x1br*1.0),npow2-1.0)) * exp(npow*x1max0 + cpow2 * pow(cpow3*(x1max0-x1br*1.0),npow2));
          FTYPE V0 = R0 + exp(npow*x1max0 + cpow2 * pow(cpow3*(x1max0-x1br*1.0),npow2));
          
          dxmax=(Rout-V0)/dVdx1; // x-x0

          dualfprintf(fail_file,"dVdx1=%g V0=%g dxmax=%g x1max=%g x1max0=%g\n",dVdx1,V0,dxmax,x1max,x1max0);
        }

        // need a slight damping factor
        FTYPE dampingfactor=0.5;
        x1max = x1max0 + dampingfactor*dxmax;


      }

      if( iter == ITERMAX ) {
        trifprintf( "Error: iteration procedure for finding x1max has not converged: x1max = %g, dx1max/x1max = %g, iter = %d\n",
                    x1max, (x1max-x1max0)/x1max, iter );
        exit(1);
      }
      else {
        trifprintf( "x1max = %g (dx1max/x1max = %g, itno = %d)\n", x1max, (x1max-x1max0)/x1max, iter );
      }
    }

    dx[1] = ( x1max - startx[1] ) /totalsize[1];
#endif


  } 
  else if (defcoord == JET6COORDSTHIN) {
    startx[1] = pow(log(Rin-R0),1.0/th_npow);
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = (pow(log(Rout-R0),1.0/th_npow)-pow(log(Rin-R0),1.0/th_npow)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 1.0/totalsize[3];

#if(1)
    startx[1] = log(Rin-R0)/th_npow;

    trifprintf( "ITERATIVE dx1: Rout=%26.20g R0=%26.20g npow=%26.20g cpow2=%26.20g npow2=%26.20g x1br=%26.20g rbr=%26.20g\n",Rout,R0,th_npow,th_cpow2,th_npow2,th_x1br,th_rbr);

    if( Rout < th_rbr ) {
      x1max = log(Rout-R0)/th_npow;
    }
    else {
      x1max0 = 1;
      x1max = 2;

      //find the root via iterations
      for( iter = 0; iter < ITERMAX; iter++ ) {

        // trifprintf( "iter=%d x1max=%26.20g x2max0=%26.20g\n",iter,x1max0,x1max);

        if( fabs((x1max - x1max0)/x1max) < RELACC ) {
          break;
        }
        x1max0 = x1max;
        dxmax= (pow( (log(Rout-R0) - th_npow*x1max0)/th_cpow2, 1./th_npow2 ) + th_x1br) - x1max0;

        // need a slight damping factor
        FTYPE dampingfactor=0.5;
        x1max = x1max0 + dampingfactor*dxmax;

      }

      if( iter == ITERMAX ) {
        trifprintf( "Error: iteration procedure for finding x1max has not converged: x1max = %g, dx1max/x1max = %g, iter = %d\n",
                    x1max, (x1max-x1max0)/x1max, iter );
        exit(1);
      }
      else {
        trifprintf( "x1max = %g (dx1max/x1max = %g, itno = %d)\n", x1max, (x1max-x1max0)/x1max, iter );
      }
    }

    dx[1] = ( x1max - startx[1] ) /totalsize[1];
#endif


  } 
  else if (defcoord == BPTHIN1) { // same as JET3COORDS since radial grid same
    startx[1] = pow(log(Rin-R0),1.0/bp_npow);
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = (pow(log(Rout-R0),1.0/bp_npow)-pow(log(Rin-R0),1.0/bp_npow)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = fracphi/totalsize[3];

#if(1)
    startx[1] = log(Rin-R0)/bp_npow;

    trifprintf( "ITERATIVE dx1: Rout=%21.15g R0=%21.15g bp_npow=%21.15g bp_cpow2=%21.15g bp_npow2=%21.15g bp_x1br=%21.15g bp_rbr=%21.15g\n",Rout,R0,bp_npow,bp_cpow2,bp_npow2,bp_x1br,bp_rbr);

    if( Rout < bp_rbr ) {
      x1max = log(Rout-R0)/bp_npow;
    }
    else {
      x1max0 = 1;
      x1max = 2;

      //find the root via iterations
      for( iter = 0; iter < ITERMAX; iter++ ) {

        //      trifprintf( "iter=%d x1max=%21.15g x2max0=%21.15g\n",iter,x1max0,x1max);

        if( fabs((x1max - x1max0)/x1max) < RELACC ) {
          break;
        }
        x1max0 = x1max;
        dxmax= (pow( (log(Rout-R0) - bp_npow*x1max0)/bp_cpow2, 1./bp_npow2 ) + bp_x1br) - x1max0;

        // need a slight damping factor
        FTYPE dampingfactor=0.5;
        x1max = x1max0 + dampingfactor*dxmax;

      }

      if( iter == ITERMAX ) {
        trifprintf( "Error: iteration procedure for finding x1max has not converged: x1max = %g, dx1max/x1max = %g, iter = %d\n",
                    x1max, (x1max-x1max0)/x1max, iter );
        exit(1);
      }
      else {
        trifprintf( "x1max = %g (dx1max/x1max = %g, itno = %d)\n", x1max, (x1max-x1max0)/x1max, iter );
      }
    }

    dx[1] = ( x1max - startx[1] ) /totalsize[1];
#endif


  } 
  else if (defcoord == JET5COORDS) {
    startx[1] = 0.0;
    startx[2] = 0.0;
    startx[3] = 0.0;
    dx[1] = 1.0 / totalsize[1];
    dx[2] = 1.0 / totalsize[2];
    dx[3] = 1.0 / totalsize[3];
  } 
  else if (defcoord == PULSARCOORDS) {
    startx[1] = pow(log(Rin-R0),1.0/npow);
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = (pow(log(Rout-R0),1.0/npow)-pow(log(Rin-R0),1.0/npow)) / totalsize[1];
    dx[2] = 1. / totalsize[2];

    dx[3] = 1.0 / totalsize[3];
    //    dx[3] = 2.0*M_PI;
  } 
  else if (defcoord == BILOGCYLCOORDS) {
    startx[1] = 0.;
    startx[2] = -1.0;
    startx[3] = 0.0;
    dx[1] = 1.0 / totalsize[1];
    dx[2] = 2.0 / totalsize[2];
    dx[3] = 1.0 / totalsize[3];
  }
  else if (defcoord == UNIFORMCOORDS) {
    //uniform grid for Cartesian coordinates
    startx[1] = 0.;
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = 1. / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 1. / totalsize[3];
  }
  else if (defcoord == RAMESHCOORDS) {
    startx[1] = pow(log(Rin-R0),1.0/npow);
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = (pow(log(Rout-R0),1.0/npow)-pow(log(Rin-R0),1.0/npow)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 1. / totalsize[3];
  } 
  else if (defcoord == RAMESHCOORDS_HALFDISK) {
    startx[1] = pow(log(Rin-R0),1.0/npow);

    if(BCtype[X2DN]==DISKSURFACE){
      startx[2] = 0.5;
    }
    else if(BCtype[X2UP]==DISKSURFACE){
      startx[2] = 0.0;
    }

    startx[3] = 0.;
    dx[1] = (pow(log(Rout-R0),1.0/npow)-pow(log(Rin-R0),1.0/npow)) / totalsize[1];
    dx[2] = 0.5 / totalsize[2];

    dx[3] = 1. / totalsize[3];
  } 
  else if (defcoord == JET4COORDS) {
    startx[1] = pow(log(Rin-R0),1.0/npow);
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = (pow(log(Rout-R0),1.0/npow)-pow(log(Rin-R0),1.0/npow)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 1. / totalsize[3];
  } 
  else if (defcoord == UNI2LOG) {
    startx[1] = 0.0;
    startx[2] = 0.0;
    startx[3] = 0.0;
    dx[1] = 1.0 / totalsize[1];
    dx[2] = 1.0 / totalsize[2];
    dx[3] = 1.0 / totalsize[3];
  }
  else{
    dualfprintf(fail_file,"Shouldn't reach end of set_points: defcoord=%d\n",defcoord);
    myexit(1);
  }

  
  dV = dx[1] * dx[2] * dx[3]; // computational volume 
  dVF = dV ; // full 3d volume

  // determine endx[]
  DLOOPA(jj) endx[jj] = startx[jj] + totalsize[jj]*dx[jj];
  // overwrite TT term so Diffx[TT] = DXDELTA, which is used
  jj=TT; endx[TT] = startx[TT] + DXDELTA;

  // total endx-startx
  // fabs() used since only ever used as absolute value
  DLOOPA(jj) Diffx[jj] = fabs(endx[jj] - startx[jj]);

  //  DLOOPA(jj) fprintf(stderr,"jj=%d Diffx=%26.20g endx=%26.20g startx=%26.20g\n",jj,Diffx[jj],endx[jj],startx[jj]);


}



#undef DXDERTYPE
#undef DXDELTA



#define MAXIHOR MAXBND
#define FRACN1 (0.1)
#define ADJUSTFRACT (0.25)


/// set i(horizon)
int setihor(void)
{
  // set to smaller of either totalsize[1]*0.1 or MAXIHOR
  if(totalsize[1]*FRACN1>MAXIHOR) return((int)MAXIHOR);
  else return((int)((FTYPE)totalsize[1]*(FTYPE)FRACN1));
}


/// there's probably a way to do this in general
/// probably can root find to get this
/// set Rin so horizon exactly on FACE1 at i=ihor
FTYPE setRin(int ihor)
{
 
  FTYPE ftemp;
  FTYPE ihoradjust;


  ihoradjust=((FTYPE)ihor)+ADJUSTFRACT; // can't have grid edge exactly on horizon due to ucon_calc()

  //  fprintf(stderr,"ihoradjust = %26.20g\n",ihoradjust);

  if(defcoord == USERCOORD){
    extern FTYPE setRin_user(int ihor, FTYPE ihoradjust);
    return(setRin_user(ihor,ihoradjust));
  }
  else if(defcoord == LOGRSINTH){
    ftemp=ihoradjust/(FTYPE)totalsize[1];
    return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
  }
  else if(defcoord == REBECCAGRID){
    ftemp=ihoradjust/(FTYPE)totalsize[1];
    return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
  }
  else if(defcoord == COMPLEX0TH){ // even though form appears different in X1, same Rin results
    ftemp=ihoradjust/(FTYPE)totalsize[1];
    return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
  }
  else if(defcoord == UNIRSINTH || defcoord == UNIRSINTH2){ // uniform
    ftemp=ihoradjust/(FTYPE)totalsize[1];
    return((Rhor-ftemp*Rout)/(1.0-ftemp));
  }
  else if(defcoord == EQMIRROR){ // as defcoord=0
    ftemp=ihoradjust/(FTYPE)totalsize[1];
    return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
  }
  else if(defcoord == COMPLEX1TH){ // as above
    ftemp=ihoradjust/(FTYPE)totalsize[1];
    return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
  }
  else if(defcoord == COMPLEX2TH){ // as above
    ftemp=ihoradjust/(FTYPE)totalsize[1];
    return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
  }
  else if(defcoord == LOGRUNITH){ // as above
    ftemp=ihoradjust/(FTYPE)totalsize[1];
    return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
  }
  else if(defcoord == JET1COORDS){
    if(npow==1.0){
      ftemp=ihoradjust/(FTYPE)totalsize[1];
      return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
    }
    else{
      // SUPERGODMARK : need to change for npow>1.0
      return(1.1);
    }
  }
  else if(defcoord == JET2COORDS){
    if(npow==1.0){
      ftemp=ihoradjust/(FTYPE)totalsize[1];
      return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
    }
    else{
      // SUPERGODMARK : need to change for npow>1.0
      return(1.1);
    }
  }
  else if(defcoord == JET3COORDS){
    // see jet3coords_checknew.nb to have chosen Rin and ihor and compute required R0
    if(npow==1.0){
      ftemp=ihoradjust/(FTYPE)totalsize[1];
      return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
    }
    else{
      dualfprintf(fail_file,"ihoradjust=%26.20g totalsize[1]=%d Rhor=%26.20g R0=%26.20g npow=%26.20g Rout=%26.20g\n",ihoradjust,totalsize[1],Rhor,R0,npow,Rout);
      return(R0+exp( pow((totalsize[1]*pow(log(Rhor-R0),1.0/npow) - ihoradjust*pow(log(Rout-R0),1.0/npow))/(totalsize[1]-ihoradjust),npow)));
    }
  }
  else if(defcoord == SJETCOORDS){
    dualfprintf( fail_file, "setRin(): not implemented for SJETCOORDS\n" );
    myexit(1);
  }
  else if(defcoord == JET6COORDS){
    // see jet3coords_checknew.nb (and fix_3dpolestissue.nb) to have chosen Rin and ihor and compute required R0
    if(npow==1.0){
      ftemp=ihoradjust/(FTYPE)totalsize[1];
      return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
    }
    else if(npow2>0){
      return(1.2);
    }
    else{
      dualfprintf(fail_file,"ihoradjust=%26.20g totalsize[1]=%d Rhor=%26.20g R0=%26.20g npow=%26.20g Rout=%26.20g\n",ihoradjust,totalsize[1],Rhor,R0,npow,Rout);
      return(R0+exp( pow((totalsize[1]*pow(log(Rhor-R0),1.0/npow) - ihoradjust*pow(log(Rout-R0),1.0/npow))/(totalsize[1]-ihoradjust),npow)));
    }
  }
  else if(defcoord == JET6COORDSTHIN){
    // see jet3coords_checknew.nb (and fix_3dpolestissue.nb) to have chosen Rin and ihor and compute required R0
    if(th_npow==1.0){
      ftemp=ihoradjust/(FTYPE)totalsize[1];
      return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
    }
    else if(th_npow2>0){
      return(1.2);
    }
    else{
      dualfprintf(fail_file,"ihoradjust=%26.20g totalsize[1]=%d Rhor=%26.20g R0=%26.20g npow=%26.20g Rout=%26.20g\n",ihoradjust,totalsize[1],Rhor,R0,th_npow,Rout);
      return(R0+exp( pow((totalsize[1]*pow(log(Rhor-R0),1.0/th_npow) - ihoradjust*pow(log(Rout-R0),1.0/th_npow))/(totalsize[1]-ihoradjust),th_npow)));
    }
  }
  else if(defcoord == BPTHIN1){
    // see jet3coords_checknew.nb (and fix_3dpolestissue.nb) to have chosen Rin and ihor and compute required R0
    if(bp_npow==1.0){
      ftemp=ihoradjust/(FTYPE)totalsize[1];
      return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
    }
    else if(bp_npow2>0){
      return(1.2);
    }
    else{
      dualfprintf(fail_file,"ihoradjust=%21.15g totalsize[1]=%d Rhor=%21.15g R0=%21.15g bp_npow=%21.15g Rout=%21.15g\n",ihoradjust,totalsize[1],Rhor,R0,bp_npow,Rout);
      return(R0+exp( pow((totalsize[1]*pow(log(Rhor-R0),1.0/bp_npow) - ihoradjust*pow(log(Rout-R0),1.0/bp_npow))/(totalsize[1]-ihoradjust),bp_npow)));
    }
  }
  else if(defcoord == JET5COORDS){
    return(JET5RIN); // not free!  Must be consistent with how determine coordinates
  }
  else if(defcoord == PULSARCOORDS){
    if(npow==1.0){
      ftemp=ihoradjust/(FTYPE)totalsize[1];
      return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
    }
    else{
      // SUPERGODMARK : need to change for npow>1.0
      return(1.1);
    }
  }
  else if (defcoord == UNIFORMCOORDS) {
    //uniform grid for Cartesian coordinates
    //no horizon in Cartesian coords -> do not set anything
    return( 1.0 );
  }
  else if(defcoord==RAMESHCOORDS || defcoord==RAMESHCOORDS_HALFDISK){
    if(npow==1.0){
      ftemp=ihoradjust/(FTYPE)totalsize[1];
      return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
    }
    else{
      // SUPERGODMARK : need to change for npow>1.0
      return(1.1);
    }
  }
  else if(defcoord == JET4COORDS){
    if(npow==1.0){
      ftemp=ihoradjust/(FTYPE)totalsize[1];
      return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
    }
    else{
      // SUPERGODMARK : need to change for npow>1.0
      return(1.1);
    }
  }
  else if(defcoord == UNI2LOG){
    dualfprintf(fail_file,"This grid's setRin is not setup defcoord=%d\n",defcoord);
    myexit(1);
    return(-1);
  }

  dualfprintf(fail_file,"Shouldn't reach end of setRin: defcoord=%d\n",defcoord);
  myexit(1);
  return(-1);

}



/// get X(i,j,k)
/// CENT is default position in degenerate cases
void coord(int i, int j, int k, int loc, FTYPE *X)
{
  // as in get_geometry(), these ?global's are used by other routines as a global indicator of where in position space we are so don't have to pass to all subfunctions
  //iglobal=i;
  //jglobal=j;
  //kglobal=k;
  //pglobal=loc;


  if (loc == CENT) {
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];
    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];
    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
  }
  else if (loc == FACE1) {
#if(N1>1)
    X[1] = startx[1] + (i + startpos[1]      ) * dx[1];
#else
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];
#endif

    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];
    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
  }
  else if (loc == FACE2) {
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];

#if(N2>1)
    X[2] = startx[2] + (j + startpos[2]      ) * dx[2];
#else 
    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];
#endif

    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
  }
  else if (loc == FACE3) {
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];
    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];

#if(N3>1)
    X[3] = startx[3] + (k + startpos[3]      ) * dx[3];
#else
    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
#endif
  }
  else if (loc == CORN1) {
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];

#if(N2>1)
    X[2] = startx[2] + (j + startpos[2]      ) * dx[2];
#else
    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];
#endif
#if(N3>1)
    X[3] = startx[3] + (k + startpos[3]      ) * dx[3];
#else
    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
#endif
  }
  else if (loc == CORN2) {
#if(N1>1)
    X[1] = startx[1] + (i + startpos[1]      ) * dx[1];
#else
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];
#endif

    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];

#if(N3>1)
    X[3] = startx[3] + (k + startpos[3]      ) * dx[3];
#else
    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
#endif
  }
  else if (loc == CORN3) {
#if(N1>1)
    X[1] = startx[1] + (i + startpos[1]      ) * dx[1];
#else
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];
#endif
#if(N2>1)
    X[2] = startx[2] + (j + startpos[2]      ) * dx[2];
#else
    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];
#endif

    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
  }
  else if (loc == CORNT) {
#if(N1>1)
    X[1] = startx[1] + (i + startpos[1]      ) * dx[1];
#else
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];
#endif
#if(N2>1)
    X[2] = startx[2] + (j + startpos[2]      ) * dx[2];
#else
    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];
#endif
#if(N3>1)
    X[3] = startx[3] + (k + startpos[3]      ) * dx[3];
#else
    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
#endif
  }

  // present time is always denoted by X[0]=0 so that past times are negative X[0]
  //  X[0]=0.0;
  // below only applies for longsteps method
  //  X[0]=t; // use X[0] to indicate time in standard temporal coordinates

  // present time of metric quantities
  // assume Xmetricnew setup in set_grid before this function (coord) is called
  X[TT] = Xmetricnew[TT];

  return;
}


/// Get X(i,j,k)
/// like coord() but does not constrain X when in reduced dimensionality
/// used for infinitesimal differences such as when numerically computing connection or dxdxp's
/// at the moment this function not used since X is defined directly rather than based upon grid
void coord_free(int i, int j, int k, int loc, FTYPE *X)
{

  //iglobal=i;
  //jglobal=j;
  //kglobal=k;
  //pglobal=loc;

  if (loc == CENT) {
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];
    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];
    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
  }
  else if (loc == FACE1) {
    X[1] = startx[1] + (i + startpos[1]      ) * dx[1];
    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];
    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
  }
  else if (loc == FACE2) {
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];
    X[2] = startx[2] + (j + startpos[2]      ) * dx[2];
    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
  }
  else if (loc == FACE3) {
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];
    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];
    X[3] = startx[3] + (k + startpos[3]      ) * dx[3];
  }
  else if (loc == CORN1) {
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];
    X[2] = startx[2] + (j + startpos[2]      ) * dx[2];
    X[3] = startx[3] + (k + startpos[3]      ) * dx[3];
  }
  else if (loc == CORN2) {
    X[1] = startx[1] + (i + startpos[1]      ) * dx[1];
    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];
    X[3] = startx[3] + (k + startpos[3]      ) * dx[3];
  }
  else if (loc == CORN3) {
    X[1] = startx[1] + (i + startpos[1]      ) * dx[1];
    X[2] = startx[2] + (j + startpos[2]      ) * dx[2];
    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
  }
  else if (loc == CORNT) {
    X[1] = startx[1] + (i + startpos[1]      ) * dx[1];
    X[2] = startx[2] + (j + startpos[2]      ) * dx[2];
    X[3] = startx[3] + (k + startpos[3]      ) * dx[3];
  }

  // present time is always denoted by X[0]=0 so that past times are negative X[0]
  //  X[0]=0.0;
  // below only applies for longsteps method
  //  X[0]=t; // use X[0] to indicate time in standard temporal coordinates

  // present time of metric quantities
  // assume Xmetricnew setup in set_grid before this function (coord) is called
  X[TT] = Xmetricnew[TT];

  return;
}



/// identical to coord except declaration using FTYPEs
void coordf(FTYPE i, FTYPE j, FTYPE k, int loc, FTYPE *X)
{
  //iglobal=ROUND2INT(i);
  //jglobal=ROUND2INT(j);
  //kglobal=ROUND2INT(k);
  //pglobal=loc;

  if (loc == CENT) {
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];
    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];
    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
  }
  else if (loc == FACE1) {
#if(N1>1)
    X[1] = startx[1] + (i + startpos[1]      ) * dx[1];
#else
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];
#endif

    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];
    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
  }
  else if (loc == FACE2) {
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];

#if(N2>1)
    X[2] = startx[2] + (j + startpos[2]      ) * dx[2];
#else 
    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];
#endif

    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
  }
  else if (loc == FACE3) {
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];
    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];

#if(N3>1)
    X[3] = startx[3] + (k + startpos[3]      ) * dx[3];
#else
    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
#endif
  }
  else if (loc == CORN1) {
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];

#if(N2>1)
    X[2] = startx[2] + (j + startpos[2]      ) * dx[2];
#else
    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];
#endif
#if(N3>1)
    X[3] = startx[3] + (k + startpos[3]      ) * dx[3];
#else
    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
#endif
  }
  else if (loc == CORN2) {
#if(N1>1)
    X[1] = startx[1] + (i + startpos[1]      ) * dx[1];
#else
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];
#endif

    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];

#if(N3>1)
    X[3] = startx[3] + (k + startpos[3]      ) * dx[3];
#else
    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
#endif
  }
  else if (loc == CORN3) {
#if(N1>1)
    X[1] = startx[1] + (i + startpos[1]      ) * dx[1];
#else
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];
#endif
#if(N2>1)
    X[2] = startx[2] + (j + startpos[2]      ) * dx[2];
#else
    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];
#endif

    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
  }
  else if (loc == CORNT) {
#if(N1>1)
    X[1] = startx[1] + (i + startpos[1]      ) * dx[1];
#else
    X[1] = startx[1] + (i + startpos[1] + 0.5) * dx[1];
#endif
#if(N2>1)
    X[2] = startx[2] + (j + startpos[2]      ) * dx[2];
#else
    X[2] = startx[2] + (j + startpos[2] + 0.5) * dx[2];
#endif
#if(N3>1)
    X[3] = startx[3] + (k + startpos[3]      ) * dx[3];
#else
    X[3] = startx[3] + (k + startpos[3] + 0.5) * dx[3];
#endif
  }

  // present time is always denoted by X[0]=0 so that past times are negative X[0]
  //  X[0]=0.0;
  // below only applies for longsteps method
  //  X[0]=t; // use X[0] to indicate time in standard temporal coordinates

  // present time of metric quantities
  // assume Xmetricnew setup in set_grid before this function (coord) is called
  X[TT] = Xmetricnew[TT];

  return;
}


/// get i,j,k(X)
void icoord(FTYPE *X,int loc, int *i, int *j, int *k)
{
  if(loc == CENT){
    *i = (int)((X[1]-startx[1])/dx[1] - 0.5) ;
    *j = (int)((X[2]-startx[2])/dx[2] - 0.5) ;
    *k = (int)((X[3]-startx[3])/dx[3] - 0.5) ;
  }

  if(startpos[1]+ *i>=totalsize[1]+N1BND) *i=totalsize[1]-1+N1BND;
  if(startpos[2]+ *j>=totalsize[2]+N2BND) *j=totalsize[2]-1+N2BND;
  if(startpos[3]+ *k>=totalsize[3]+N3BND) *k=totalsize[3]-1+N3BND;
  
  if(startpos[1]+ *i<-N1BND) *i=-N1BND;
  if(startpos[2]+ *j<-N2BND) *j=-N2BND;
  if(startpos[3]+ *k<-N3BND) *k=-N3BND;

}

#define INCLUDEROUND 0 // default

#ifdef WIN32
#undef INCLUDEROUND
#define INCLUDEROUND 1
#endif

// roundl is not part of long double library
#if(REALTYPE==LONGDOUBLETYPE)
#undef INCLUDEROUND
#define INCLUDEROUND 1
#endif

#if(INCLUDEROUND)
/// get rounded value of x
FTYPE round(FTYPE x)
{
  FTYPE xfloor,xceil;

  xfloor=floor(x);
  xceil=ceil(x);
  if(fabs(x-xfloor)>fabs(x-xceil)) return(xceil);
  else return(xfloor);
}

/// get rounded value of x
long int lrint(FTYPE x)
{
  return((long int)round(x));
}
#endif



/// get rounded value of x
FTYPE myround(FTYPE x)
{

#if(REALTYPE==FLOATTYPE)
  return(roundf(x));
#elif(REALTYPE==DOUBLETYPE || REALTYPE==LONGDOUBLETYPE)
  // ceill will be used automatically by conversion to long doubles
  return(round(x));
#endif

}


#if(0)
/// get rounded value of x
long int myround2int(FTYPE x)
{

#if(REALTYPE==FLOATTYPE)
  return(lrintf(x));
#elif(REALTYPE==DOUBLETYPE)
  return(lrint(x));
#elif(REALTYPE==LONGDOUBLETYPE)
  return(lrintl(x));
#endif
}
#else
#define myround2int(x) (ROUND2INT(x))
#endif

/// get i,j,k(X) rounded to integer
void icoord_round(FTYPE *X,int loc, int *i, int *j, int *k)
{
  FTYPE myround(FTYPE x);
  //  long int myround2int(FTYPE x);

  if(loc == CENT){
    *i = myround2int((X[1]-startx[1])/dx[1] - 0.5) ;
    *j = myround2int((X[2]-startx[2])/dx[2] - 0.5) ;
    *k = myround2int((X[3]-startx[3])/dx[3] - 0.5) ;
  }

  if(startpos[1]+ *i>=totalsize[1]+N1BND) *i=totalsize[1]-1+N1BND;
  if(startpos[2]+ *j>=totalsize[2]+N2BND) *j=totalsize[2]-1+N2BND;
  if(startpos[3]+ *k>=totalsize[3]+N3BND) *k=totalsize[3]-1+N3BND;
  
  if(startpos[1]+ *i<-N1BND) *i=-N1BND;
  if(startpos[2]+ *j<-N2BND) *j=-N2BND;
  if(startpos[3]+ *k<-N3BND) *k=-N3BND;

}


/// see if point inside surface
/// dir=X1UP,X1DN,etc.
int is_inside_surface(int dir, int ii, int jj, int kk, int pp)
{
  int is_on_surface(int dir, int ii, int jj, int kk, int pp);
  int ijk[NDIM];
  int dimen;
  int dirsign;
  int isonsurf;

  ijk[1]=ii;
  ijk[2]=jj;
  ijk[3]=kk;


  dimen=DIMEN(dir);
  dirsign=DIRSIGN(dir);
  
  isonsurf=is_on_surface(dir,ii,jj,kk,pp);

  if(isonsurf){
    return(1);
  }
  else if(
          (startpos[dimen]+ijk[dimen]<0 && dirsign==-1)
          || 
          (startpos[dimen]+ijk[dimen]>totalsize[dimen]-1 && dirsign==1)
          ){
    return(1);
  }
  else return(0);

}

/// see if point is on surface
int is_on_surface(int dir, int ii, int jj, int kk, int pp)
{
  int ijk[NDIM];
  int dimen;
  int dirsign;

  ijk[1]=ii;
  ijk[2]=jj;
  ijk[3]=kk;


  dimen=DIMEN(dir);
  dirsign=DIRSIGN(dir);

  if(
     ((dimen==1)&&(dirsign==-1)&&(startpos[dimen]+ii==0)                 && (pp==FACE1 || pp==CORN2 || pp==CORN3 || pp==CORNT))||
     ((dimen==1)&&(dirsign==1) &&(startpos[dimen]+ii==totalsize[dimen])  && (pp==FACE1 || pp==CORN2 || pp==CORN3 || pp==CORNT))||
     ((dimen==2)&&(dirsign==-1)&&(startpos[dimen]+jj==0)                 && (pp==FACE2 || pp==CORN1 || pp==CORN3 || pp==CORNT))||
     ((dimen==2)&&(dirsign==1) &&(startpos[dimen]+jj==totalsize[dimen])  && (pp==FACE2 || pp==CORN1 || pp==CORN3 || pp==CORNT))||
     ((dimen==3)&&(dirsign==-1)&&(startpos[dimen]+kk==0)                 && (pp==FACE3 || pp==CORN2 || pp==CORN1 || pp==CORNT))||
     ((dimen==3)&&(dirsign==1) &&(startpos[dimen]+kk==totalsize[dimen])  && (pp==FACE3 || pp==CORN2 || pp==CORN1 || pp==CORNT))
     ){
    return(1);
  }
  else return(0);

}



////////// BELOW ARE ACCESSING STORE POSITION INFORMATION
// GODMARK: Should improve this.  Right now if I call coord_ijk and then immediately bl_coord_ijk and first time called, then expensive since coord() called twice
// NEVER CALL *_ijk_*() type functions if input i,j,k is meant to be related to total grid that would access beyond stored arrays


/// normally-used dxdxp[i,j,k] = dV/dX
void dxdxprim_ijk(int i, int j, int k, int loc, FTYPE (*dxdxp)[NDIM])
{
  int jj,kk;
  FTYPE X[NDIM],V[NDIM];

  //iglobal=i;
  //jglobal=j;
  //kglobal=k;
  //pglobal=loc;


  if(didstorepositiondata && loc!=NOWHERE){
    DLOOP(jj,kk) dxdxp[jj][kk] = GLOBALMETMACP1A2(dxdxpstore,loc,i,j,k,jj,kk);
#if(PRODUCTION==0)
    if(i<-N1BND || i>N1-1+SHIFT1+N1BND || j<-N2BND || j>N2-1+SHIFT2+N2BND || k<-N3BND || k>N3-1+SHIFT3+N3BND){
      dualfprintf(fail_file,"Beyond stored location: %d %d %d\n",i,j,k);
      myexit(10365262);
    }
#endif
  }
  else if(loc!=NOWHERE){
    bl_coord_ijk_2(i, j, k, loc, X, V);
    dxdxprim(X,V,dxdxp);
  }
  else{
    dualfprintf(fail_file,"dxdxprim_ijk(): No X to compute V or dxdxp ijk=%d %d %d loc=%d\n",i,j,k,loc);
    myexit(8813);
  }

}


/// normally-used dxdxp[i,j,k]
void dxdxprim_ijk_2(struct of_geom *ptrgeom, FTYPE *X, FTYPE *V, FTYPE (*dxdxp)[NDIM])
{
  int i,j,k,loc;
  int jj,kk;



  i=ptrgeom->i;
  j=ptrgeom->j;
  k=ptrgeom->k;
  loc=ptrgeom->p;

  //iglobal=i;
  //jglobal=j;
  //kglobal=k;
  //pglobal=loc;


  if(didstorepositiondata && loc!=NOWHERE){
    DLOOP(jj,kk) dxdxp[jj][kk] = GLOBALMETMACP1A2(dxdxpstore,loc,i,j,k,jj,kk);
#if(PRODUCTION==0)
    if(i<-N1BND || i>N1-1+SHIFT1+N1BND || j<-N2BND || j>N2-1+SHIFT2+N2BND || k<-N3BND || k>N3-1+SHIFT3+N3BND){
      dualfprintf(fail_file,"Beyond stored location: %d %d %d\n",i,j,k);
      myexit(10365263);
    }
#endif
  }
  else if(loc!=NOWHERE){
    bl_coord_ijk_2(i, j, k, loc, X, V);
    dxdxprim(X,V,dxdxp);
  }
  else{
    dxdxprim(X,V,dxdxp);
  }

}

/// Get (dV/dX)^{-1}
void idxdxprim(FTYPE (*dxdxp)[NDIM], FTYPE (*idxdxp)[NDIM])
{

  matrix_inverse(PRIMECOORDS, dxdxp,idxdxp);

}




/// normally-used dxdxp[i,j,k]
void idxdxprim_ijk(int i, int j, int k, int loc, FTYPE (*idxdxp)[NDIM])
{
  int jj,kk;
  FTYPE X[NDIM],V[NDIM];
  FTYPE dxdxp[NDIM][NDIM];

  //iglobal=i;
  //jglobal=j;
  //kglobal=k;
  //pglobal=loc;


  if(didstorepositiondata && loc!=NOWHERE){
    DLOOP(jj,kk) idxdxp[jj][kk] = GLOBALMETMACP1A2(idxdxpstore,loc,i,j,k,jj,kk);
#if(PRODUCTION==0)
    if(i<-N1BND || i>N1-1+SHIFT1+N1BND || j<-N2BND || j>N2-1+SHIFT2+N2BND || k<-N3BND || k>N3-1+SHIFT3+N3BND){
      dualfprintf(fail_file,"Beyond stored location: %d %d %d\n",i,j,k);
      myexit(10365264);
    }
#endif
  }
  else if(loc!=NOWHERE){
    bl_coord_ijk_2(i, j, k, loc, X, V);
    dxdxprim(X,V,dxdxp);
    idxdxprim(dxdxp,idxdxp);
  }
  else{
    dualfprintf(fail_file,"idxdxprim_ijk(): No X, V, or dxdxp to compute idxdxp\n");
    myexit(8814);
  }

}


/// normally-used dxdxp[i,j,k]
/// {X,V} -> idxdxp only if ptrgoem->loc=NOWHERE
/// else {ptrgeom->{i,j,k}} -> {X,V,idxdxp}
void idxdxprim_ijk_2(struct of_geom *ptrgeom, FTYPE *X, FTYPE *V, FTYPE (*idxdxp)[NDIM])
{
  int i,j,k,loc;
  int jj,kk;
  FTYPE dxdxp[NDIM][NDIM];


  i=ptrgeom->i;
  j=ptrgeom->j;
  k=ptrgeom->k;
  loc=ptrgeom->p;

  //iglobal=i;
  //jglobal=j;
  //kglobal=k;
  //pglobal=loc;


  if(didstorepositiondata && loc!=NOWHERE){
    // if here, don't assume X and V already set, so set them both:
    bl_coord_ijk_2(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,X,V);
    DLOOP(jj,kk) idxdxp[jj][kk] = GLOBALMETMACP1A2(idxdxpstore,loc,i,j,k,jj,kk);
#if(PRODUCTION==0)
    if(i<-N1BND || i>N1-1+SHIFT1+N1BND || j<-N2BND || j>N2-1+SHIFT2+N2BND || k<-N3BND || k>N3-1+SHIFT3+N3BND){
      dualfprintf(fail_file,"Beyond stored location: %d %d %d\n",i,j,k);
      myexit(10365265);
    }
#endif
  }
  else if(loc!=NOWHERE){
    // then assume ptrgeom->{i,j,k} is good for ptrgeom->p
    bl_coord_ijk_2(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,X,V);
    dxdxprim(X,V,dxdxp);
    idxdxprim(dxdxp,idxdxp);
  }
  else{
    // if idxdxprim_ijk_2 is called, assume fed in X and V already
    dxdxprim(X,V,dxdxp);
    idxdxprim(dxdxp,idxdxp);
  }

}



/// normally-used bl_coord[i,j,k]
/// {i,j,k} -> V only
void bl_coord_ijk(int i, int j, int k, int loc, FTYPE *V)
{
  int jj;
  FTYPE X[NDIM];

  //iglobal=i;
  //jglobal=j;
  //kglobal=k;
  //pglobal=loc;

  

  if(didstorepositiondata && loc!=NOWHERE){
    DLOOPA(jj) V[jj] = GLOBALMACP1A1(Vstore,loc,i,j,k,jj);
#if(PRODUCTION==0)
    if(i<-N1BND-SHIFT1 || i>N1-1+SHIFT1*2+N1BND || j<-N2BND-SHIFT2 || j>N2-1+SHIFT2*2+N2BND || k<-N3BND-SHIFT3 || k>N3-1+SHIFT3*2+N3BND){
      dualfprintf(fail_file,"Beyond stored location: %d %d %d\n",i,j,k);
      myexit(10365266);
    }
#endif
    V[TT]=Xmetricnew[TT];
  }
  else if(loc!=NOWHERE){
    coord(i, j, k, loc, X);
    bl_coord(X,V);
  }
  else{
    dualfprintf(fail_file,"bl_coord_ijk(): No X to compute V : ijk=%d %d %d loc=%d\n",i,j,k,loc);
    myexit(8815);
  }


}

/// normally-used bl_coord[i,j,k]
/// {i,j,k} -> {X,V} both if loc!=NOWHERE
/// X -> V if loc==NOWHERE (so reduces to bl_coord(X,V)
void bl_coord_ijk_2(int i, int j, int k, int loc, FTYPE *X, FTYPE *V)
{
  int jj;

  //iglobal=i;
  //jglobal=j;
  //kglobal=k;
  //pglobal=loc;

  
  //  dualfprintf(fail_file,"i=%d didstorepositiondata=%d loc=%d\n",i,didstorepositiondata,loc);

  if(didstorepositiondata && loc!=NOWHERE){

    coord_ijk(i,j,k,loc,X); // get X

    DLOOPA(jj) V[jj] = GLOBALMACP1A1(Vstore,loc,i,j,k,jj);
#if(PRODUCTION==0)
    if(i<-N1BND-SHIFT1 || i>N1-1+SHIFT1*2+N1BND || j<-N2BND-SHIFT2 || j>N2-1+SHIFT2*2+N2BND || k<-N3BND-SHIFT3 || k>N3-1+SHIFT3*2+N3BND){
      dualfprintf(fail_file,"Beyond stored location: %d %d %d\n",i,j,k);
      myexit(10365267);
    }
#endif
    V[TT]=Xmetricnew[TT];  // not really true, but dont use V[TT] yet
  }
  else if(loc!=NOWHERE){
    // Assume X is output, not input, if using bl_coord_ijk_?.  If didn't know i,j,k, then would have been forced to use bl_coord(X,V) anyways.
    coord(i,j,k,loc,X);
    bl_coord(X,V);
  }
  else{
    // then must use X directly
    bl_coord(X,V);
    // if really NOWHERE, then not associated with i,j,k at all, so can't use coord() to get X, so shouldn't be using bl_coord_ijk_2() in first place.  But assume X is already set
  }


}

/// normally-used bl_coord[i,j,k]
/// returns both X and V given i,j,k,loc and does NOT try to use stored memory -- presumes could input arbitrary i,j,k not on a single CPU
void bl_coord_coord(int i, int j, int k, int loc, FTYPE *X, FTYPE *V)
{
  
  coord(i,j,k,loc,X);
  bl_coord(X,V);

}


/// normally-used bl_coord[i,j,k]
/// {i,j,k} -> X only
void coord_ijk(int i, int j, int k, int loc, FTYPE *X)
{
  int jj;

  //iglobal=i;
  //jglobal=j;
  //kglobal=k;
  //pglobal=loc;

  
  if(didstorepositiondata && loc!=NOWHERE){
    DLOOPA(jj) X[jj] = GLOBALMACP1A1(Xstore,loc,i,j,k,jj);
#if(PRODUCTION==0)
    if(i<-N1BND-SHIFT1 || i>N1-1+SHIFT1*2+N1BND || j<-N2BND-SHIFT2 || j>N2-1+SHIFT2*2+N2BND || k<-N3BND-SHIFT3 || k>N3-1+SHIFT3*2+N3BND){
      dualfprintf(fail_file,"Beyond stored location: %d %d %d\n",i,j,k);
      myexit(10365268);
    }
#endif
    // time is not constant
    X[TT] = Xmetricnew[TT];
  }
  else if(loc!=NOWHERE){
    // coord must use i,j,k, so always safe for any loc?
    coord(i, j, k, loc, X);
  }
  else{
    dualfprintf(fail_file,"coord_ijk() NOWHERE: Can't use loc=%d to compute X\n",loc);
    myexit(8816);
  }


}



/// smooth step function: 
/// Ftr = 0 if x < 0, Ftr = 1 if x > 1 and smoothly interps. in btw.
FTYPE Ftr( FTYPE x )
{
  FTYPE res;

  if( x <= 0 ) {
    res = 0;
  }
  else if( x >= 1 ) {
    res = 1;
  }
  else {
    res = (64 + cos(5*M_PI*x) + 70*sin((M_PI*(-1 + 2*x))/2.) + 5*sin((3*M_PI*(-1 + 2*x))/2.))/128.;
  }

  return( res );
}

FTYPE Ftrgenlin( FTYPE x, FTYPE xa, FTYPE xb, FTYPE ya, FTYPE yb )
{
  FTYPE Ftr( FTYPE x );
  FTYPE res;

  res = (x*ya)/xa + (-((x*ya)/xa) + ((x - xb)*(1 - yb))/(1 - xb) + yb)*Ftr((x - xa)/(-xa + xb));
  
  return( res );
}

FTYPE Ftrgen( FTYPE x, FTYPE xa, FTYPE xb, FTYPE ya, FTYPE yb )
{
  FTYPE Ftr( FTYPE x );
  FTYPE res;

  res = ya + (yb-ya)*Ftr( (x-xa)/(xb-xa) );
  
  return( res );
}

FTYPE Fangle( FTYPE x )
{
  FTYPE res;

  if( x <= -1 ) {
    res = 0;
  }
  else if( x >= 1 ) {
    res = x;
  }
  else {
    res = (1 + x + (-140*sin((M_PI*(1 + x))/2.) + (10*sin((3*M_PI*(1 + x))/2.))/3. + (2*sin((5*M_PI*(1 + x))/2.))/5.)/(64.*M_PI))/2.;
  }

  return( res );
  
}

FTYPE limlin( FTYPE x, FTYPE x0, FTYPE dxx, FTYPE y0 )
{
  FTYPE Fangle( FTYPE x );
  return( y0 - dxx * Fangle(-(x-x0)/dxx) );
}

FTYPE minlin( FTYPE x, FTYPE x0, FTYPE dxx, FTYPE y0 )
{
  FTYPE Fangle( FTYPE x );
  return( y0 + dxx * Fangle((x-x0)/dxx) );
}

FTYPE mins( FTYPE f1, FTYPE f2, FTYPE df )
{
  FTYPE limlin( FTYPE x, FTYPE x0, FTYPE dxx, FTYPE y0 );
  return( limlin(f1, f2, df, f2) );
}

FTYPE maxs( FTYPE f1, FTYPE f2, FTYPE df )
{
  FTYPE mins( FTYPE f1, FTYPE f2, FTYPE df );
  return( -mins(-f1, -f2, df) );
}


/// =mins if dir < 0
/// =maxs if dir >= 0
FTYPE minmaxs( FTYPE f1, FTYPE f2, FTYPE df, FTYPE dir )
{
  FTYPE mins( FTYPE f1, FTYPE f2, FTYPE df );
  FTYPE maxs( FTYPE f1, FTYPE f2, FTYPE df );
  if( dir>=0 ) {
    return( maxs(f1, f2, df) );
  }
  
  return( mins(f1, f2, df) );
}

static FTYPE sinth0( FTYPE *X0, FTYPE *X, void (*vofx)(FTYPE*, FTYPE*) );
static FTYPE sinth1in( FTYPE *X0, FTYPE *X, void (*vofx)(FTYPE*, FTYPE*) );
static FTYPE th2in( FTYPE *X0, FTYPE *X, void (*vofx)(FTYPE*, FTYPE*) );
static void to1stquadrant( FTYPE *Xin, FTYPE *Xout, int *ismirrored );
static FTYPE func1( FTYPE *X0, FTYPE *X,  void (*vofx)(FTYPE*, FTYPE*) );
static FTYPE func2( FTYPE *X0, FTYPE *X,  void (*vofx)(FTYPE*, FTYPE*) );

/// for sjet coords
///Converts copies Xin to Xout and converts
///but sets Xout[2] to lie in the 1st quadrant, i.e. Xout[2] \in [-1,0])
///if the point had to be mirrored
void to1stquadrant( FTYPE *Xin, FTYPE *Xout, int *ismirrored )
{
  FTYPE ntimes;
  int dim;
  
  DLOOPA(dim) Xout[dim] = Xin[dim];
  
  //bring the angle variables to -2..2 (for X) and -2pi..2pi (for V)
  ntimes = floor( (Xin[2]+2.0)/4.0 );
  //this forces -2 < Xout[2] < 2
  Xout[2] -= 4 * ntimes;
  
  *ismirrored = 0;
  
  if( Xout[2] > 0. ) {
    Xout[2] = -Xout[2];
    *ismirrored = 1-*ismirrored;
  }    

  //now force -1 < Xout[2] < 0
  if( Xout[2] < -1. ) {
    Xout[2] = -2. - Xout[2];
    *ismirrored = 1-*ismirrored;
  }
}

/// for sjet coords
FTYPE sinth0( FTYPE *X0, FTYPE *X, void (*vofx)(FTYPE*, FTYPE*) )
{
  FTYPE V0[NDIM];
  FTYPE Vc0[NDIM];
  FTYPE Xc0[NDIM];
  int dim;
  
  //X1 = {0, X[1], X0[1], 0}
  DLOOPA(dim) Xc0[dim] = X[dim];
  Xc0[2] = X0[2];
  
  vofx( Xc0, Vc0 );
  vofx( X0, V0 );
  
  
  return( V0[1] * sin(V0[2]) / Vc0[1] );
}

/// for sjet coords
FTYPE sinth1in( FTYPE *X0, FTYPE *X, void (*vofx)(FTYPE*, FTYPE*) )
{
  FTYPE V[NDIM];
  FTYPE V0[NDIM];
  FTYPE V0c[NDIM];
  FTYPE X0c[NDIM];
  int dim;
  
  //X1 = {0, X[1], X0[1], 0}
  DLOOPA(dim) X0c[dim] = X0[dim];
  X0c[2] = X[2];
  
  vofx( X, V );
  vofx( X0c, V0c );
  vofx( X0, V0 );
  
  return( V0[1] * sin(V0c[2]) / V[1] );
}


/// for sjet coords
FTYPE th2in( FTYPE *X0, FTYPE *X, void (*vofx)(FTYPE*, FTYPE*) )
{
  FTYPE V[NDIM];
  FTYPE V0[NDIM];
  FTYPE Vc0[NDIM];
  FTYPE Xc0[NDIM];
  FTYPE Xcmid[NDIM];
  FTYPE Vcmid[NDIM];
  int dim;
  FTYPE res;
  FTYPE th0;
  
  DLOOPA(dim) Xc0[dim] = X[dim];
  Xc0[2] = X0[2];
  vofx( Xc0, Vc0 ); 
  
  DLOOPA(dim) Xcmid[dim] = X[dim];
  Xcmid[2] = 0;
  vofx( Xcmid, Vcmid ); 

  vofx( X0, V0 ); 
  vofx( X, V ); 
  
  th0 = asin( sinth0(X0, X, vofx) );
  
  res = (V[2] - Vc0[2])/(Vcmid[2] - Vc0[2]) * (Vcmid[2]-th0) + th0;
  
  return( res );
}

/// for sjet coords
///Adjusts V[2]=theta so that a few innermost cells around the pole
///become cylindrical
///ASSUMES: poles are at 
///            X[2] = -1 and +1, which correspond to
///            V[2] = 0 and pi
void vofx_cylindrified( FTYPE *Xin, void (*vofx)(FTYPE*, FTYPE*), FTYPE *Vout )
{
  FTYPE npiovertwos;
  FTYPE X[NDIM], V[NDIM];
  FTYPE Vin[NDIM];
  FTYPE X0[NDIM], V0[NDIM];
  FTYPE Xtr[NDIM], Vtr[NDIM];
  FTYPE f1, f2, dftr;
  FTYPE sinth, th;
  int dim, ismirrored;
  
  vofx( Xin, Vin );

  // BRING INPUT TO 1ST QUADRANT:  X[2] \in [-1 and 0]
  to1stquadrant( Xin, X, &ismirrored );  
  vofx( X, V );

  //initialize X0: cylindrify region
  //X[1] < X0[1] && X[2] < X0[2] (value of X0[3] not used)
  X0[0] = Xin[0];
  X0[1] = x10;
  X0[2] = x20;
  X0[3] = 0;
  vofx( X0, V0 );
  
  //{0, roughly midpoint between grid origin and x10, -1, 0}
  DLOOPA(dim) Xtr[dim] = X[dim];
  Xtr[1] = log( 0.5*( exp(X0[1])+exp(startx[1]) ) );   //always bound to be between startx[1] and X0[1]
  vofx( Xtr, Vtr );

  f1 = func1( X0, X, vofx );
  f2 = func2( X0, X, vofx );
  dftr = func2( X0, Xtr, vofx ) - func1( X0, Xtr, vofx );
      
  // Compute new theta
  sinth = maxs( V[1]*f1, V[1]*f2, Vtr[1]*fabs(dftr)+SMALL ) / V[1];
  
  th = asin( sinth ); 
  
  //initialize Vout with the original values
  DLOOPA(dim) Vout[dim] = Vin[dim];

  //apply change in theta in the original quadrant
  if( 0 == ismirrored ) {
    Vout[2] = Vin[2] + (th - V[2]);
  }
  else {
    //if mirrrored, flip the sign 
    Vout[2] = Vin[2] - (th - V[2]);
  }
}

/// for sjet coords
FTYPE func1( FTYPE *X0, FTYPE *X,  void (*vofx)(FTYPE*, FTYPE*) )
{
  FTYPE V[NDIM];
  
  vofx( X, V );
  
  return( sin(V[2]) ); 
}

/// for sjet coords
FTYPE func2( FTYPE *X0, FTYPE *X,  void (*vofx)(FTYPE*, FTYPE*) )
{
  FTYPE V[NDIM];
  FTYPE Xca[NDIM];
  FTYPE func2;
  int dim;
  FTYPE sth1in, sth2in, sth1inaxis, sth2inaxis;
  
  //{0, X[1], -1, 0}
  DLOOPA(dim) Xca[dim] = X[dim];
  Xca[2] = -1;
  
  vofx( X, V );
  
  sth1in = sinth1in( X0, X, vofx );
  sth2in = sin( th2in(X0, X, vofx) );
  
  sth1inaxis = sinth1in( X0, Xca, vofx );
  sth2inaxis = sin( th2in(X0, Xca, vofx) );
  
  func2 = minmaxs( sth1in, sth2in, fabs(sth2inaxis-sth1inaxis)+SMALL, X[1] - X0[1] );
  
  return( func2 ); 
}




















































