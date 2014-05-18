
/*! \file mpi_fileio.c
     \brief Various file input/output using MPI functions
// Notes:
//
// mpicombinetype is global initial user choice for combinetype for MPI
// truempicombinetype is global system-forced choice for mpicombinetype set during initialization (and so used rest of the time except new initialization where we start fresh).  So "mpiio_final()" also uses truempicombinetype since related to last choice
// OPENMPMARK: Assume all calls to MPI file I/O routines are not done by multiple threads.
*/

#include "decs.h"





/// initialize mpiio
void mpiio_init(int bintxt, int sorted, FILE ** fpptr, long headerbytesize, int which, char *filename, int numcolumns,
                MPI_Datatype datatype, void **jonioptr, void **writebufptr)
{


  logsfprintf("\nmpiio_init begin\n");


  /////////////
  //
  // this check covers combine and seperate
  //
  /////////////
  if(!sorted){
    if(bintxt==BINARYOUTPUT){
      dualfprintf(fail_file,"No such thing as binary unsorted output\n");
      myexit(1);
    }
  }


  //////////////////////////
  //
  // clean up prior call (internally handle mpicombinetype)
  //
  //////////////////////////
  mpiio_final(bintxt, sorted, fpptr, headerbytesize, which, filename, numcolumns, datatype, jonioptr, writebufptr);


  //////////////////////////
  //
  // Initialize new call
  //
  //////////////////////////

  if(mpicombinetype==MPICOMBINEROMIO){
    mpiioromio_init_combine(INITROMIO, which, headerbytesize, filename,numcolumns, datatype,writebufptr,(void*)0x0);
  }
  else{
    mpiios_init(bintxt, sorted, fpptr, which, headerbytesize, filename, numcolumns, datatype,
                jonioptr, writebufptr);
  }

  logsfprintf("mpiio_init end\n");

}




/// used to cleanup file writing in non-blocking way
/// should be called before writing if previously wrote and also called at end of calculation to close last file (use fake write to do so)
void mpiio_final(int bintxt, int sorted, FILE ** fpptr, long headerbytesize, int which, char *filename, int numcolumns,
                 MPI_Datatype datatype, void **jonioptr, void **writebufptr)
{
  static int priorinit=0;



  if(priorinit){

    logsfprintf("\nmpiio_final begin\n");

    // "true" value same
    // doesn't use jonioptr or bintxt or sorted or fpptr
    // writebuf not set yet for combine or seperate, have to allocate at address writebufptr first.
    // then assume previously wrote file that now needs to be fully ended before one can write another
    mpiioromio_init_combine(WRITEENDROMIO, which, headerbytesize, filename,numcolumns, datatype,writebufptr,*writebufptr);
  
    logsfprintf("mpiio_final end\n");

  }


  //////////////////////////
  //
  // priorinit indicates to *next call* of this function that mpiio_final(priorinit) should end/close/cleanup previously written file if was using ROMIO and writing in this call
  // otherwise assume file writing call was fully completed already before
  // below assumes mpicombinetype==truempicombinetype when doing ROMIO
  //
  //////////////////////////

  if(mpicombinetype==MPICOMBINEROMIO && which==WRITEFILE){
    priorinit=1;
  }
  else priorinit=0;

  // if numcolumns==-1, then separate finish of ROMIO call
  if(numcolumns==-1) priorinit=0;

}




void mpiio_combine(int bintxt, int sorted,
                   int numcolumns, MPI_Datatype datatype,
                   FILE ** fpptr, void *jonio, void *writebuf)
{
#if(USEMPI)

  logsfprintf("mpiio start combine\n");

  if(USEROMIO){
    // doesn't use jonioptr or bintxt or sorted or fpptr
    // address not used, just value of writebuf
    // headerbytesize no longer needed
    mpiioromio_init_combine(WRITECLOSEROMIO, WRITEFILE,  0, "", numcolumns, datatype,&writebuf,writebuf);
  }
  else{
    if(truempicombinetype==MPICOMBINESIMPLE){
      
      if (sorted) {
        mpiios_combine(bintxt, datatype, numcolumns,fpptr, jonio, writebuf);
      }
      else{
        mpiiotu_combine(datatype, numcolumns, fpptr, writebuf);
      }
    }
    else if(truempicombinetype==MPICOMBINEMINMEM){
      mpiiomin_final(numcolumns,fpptr, jonio, writebuf);
    }
  }
#endif
  logsfprintf("mpiio end combine\n");

}

#define DEBUGMINMEM 0


/// this can be made to work without mpi for super-efficient single cpu mode, but not important.  Can still do 1 cpu with MPI and use this!
void mpiio_minmem(int readwrite, int whichdump, int i, int j, int k, int bintxt, int sorted,
                  int numcolumns, MPI_Datatype datatype,
                  FILE ** fp, void *jonio, void *writebuf)
{
  static struct blink * blinkptr;
  static struct blink * cpulinkptr;
  static long long int datasofar0,datasofarc0;
  static int done0=0;
  static long long int datagot;
  long long int bfi;
  long long int sii,uii;// sorted and unsorted index sometimes
#if(USEMPI)
  MPI_Request request;
  MPI_Request request0;
#endif
  long long int joniooffset;
  void *joniosubmit;

  unsigned char *jonio1;
  float *jonio4;
  double *jonio8;
  long double *jonio16;
  int *jonio4i;
  long long int *jonio8i;

  long long int mapjoniosorted,mapjoniounsorted;
  long long int datatoget0,datatogetc0,datasent0;
  int doset;
  int dofull;
  int dolocal;
  int doing0;
  long long int gi,gj,gk;//global i,j,k starting from first cpu reference and adding the sorted index
  long long int lastnum;
  static int masterdone;
  static int thisslavedone;

  int numfiles;
  int sizeofdatatype;
  long long int mypos;
  unsigned short short4char;


#if(USEMPI)

  sizeofdatatype=getsizeofdatatype(datatype);

#if(DEBUGMINMEM)
  dualfprintf(fail_file,"got here 0: i=%d j=%d k=%d\n",i,j,k); fflush(fail_file);
#endif

  if(sorted==UNSORTED){
    dualfprintf(fail_file,"Not setup to do unsorted with this method\n");
    myexit(10001);
  }


  /////////////
  //
  // very quickly determine whether at required point to begin read or write of data, or still busy with buffer.
  //
  ///////////
  if((i==0)&&(j==0)&&(k==0)){
    // then first time here for this cpu for this dump
    blinkptr=blinkptr0[whichdump];
    datagot=0;
    if(myid==0) masterdone=0;
    thisslavedone=0;
  }
  if(masterdone&&(myid==0)) return;
  if(thisslavedone&&(myid!=0)) return;


  if(readwrite==WRITEFILE){
    // i+1 since called inside loop before i++, but really has done i++ things
    bfi=((long long int)k*(long long int)N1*(long long int)N2+(long long int)j*(long long int)N1+(long long int)(i+1))*(long long int)numcolumns; // recall that we are always at col=0 since we always complete a column
  }
  else if(readwrite==READFILE){
    bfi=((long long int)k*(long long int)N1*(long long int)N2+(long long int)j*(long long int)N1+(long long int)i)*(long long int)numcolumns;
  }


  if(blinkptr==NULL) dolocal=0; // end of local list
  else{
    // see if completed a single node in the list yet
    if(readwrite==WRITEFILE){
      if(bfi==((long long int)datagot+(long long int)blinkptr->num)) dolocal=1; // at a completion point, need to send to cpu=0
      else   dolocal=0; // not yet completed
    }
    else if(readwrite==READFILE){
      if(bfi==(long long int)datagot) dolocal=1; // at a completion point, need to recv from cpu=0
      else   dolocal=0; // not yet completed
    }
  }

  ///////////////////////
  //
  // ANY CPU Send/Recv a buffer to cpu=0
  //
  ////////////////////////



  if(dolocal){
#if(DEBUGMINMEM)
    dualfprintf(fail_file,"got here 0.5: %lld %lld %lld %d\n",bfi,datagot+blinkptr->num,blinkptr->num,dolocal); fflush(fail_file);
#endif
    // We are ready to read or write from/to cpu=0
#if(DEBUGMINMEM)
    dualfprintf(fail_file,"got here 0.6\n"); fflush(fail_file);
#endif


    ///////////////////////
    //
    // WRITEFILE, so other CPU sends to cpu=0
    //
    ////////////////////////
    
    if(readwrite==WRITEFILE){
      // means we have put what we want into writebuf, now send to cpu=0
#if(DEBUGMINMEM)
      dualfprintf(fail_file,"got here 0.65 : %lld %lld %d %d\n",writebuf,blinkptr->num,MPIid[0],myid); fflush(fail_file);
#endif
      // must keep myid>0 cpus stalled until cpu=0 needs their data
      // that is, Wait below continues if data copied out of buffer, but we want pseudo-blocking call here
      if(myid>0) MPI_Issend(writebuf,blinkptr->num,datatype,MPIid[0],myid, MPI_COMM_GRMHD,&request);
      // can't stall cpu=0 since only below has recv, but ok since cpu=0 stuck below until cpu=0 data needed again
      else MPI_Isend(writebuf,blinkptr->num,datatype,MPIid[0],myid, MPI_COMM_GRMHD,&request);
#if(DEBUGMINMEM)
      dualfprintf(fail_file,"got here 0.66 : %lld\n",writebuf); fflush(fail_file);
#endif
      // adjust offset to keep standard writebuf BUFFERMAP format with just an offset (offset keeps true memory area always as written part)
      bufferoffset -= (long long int)blinkptr->num;// or =-(long long int)datagot+(long long int)blinkptr->num


    }
    else if(readwrite==READFILE){

      ///////////////////////
      //
      // READFILE, so other CPU recvs from cpu=0
      //
      ////////////////////////


      // means we need to fill writebuf with what cpu=0 will send us
      // recv's will Wait below until data is coming.
#if(DEBUGMINMEM)
      dualfprintf(fail_file,"got here 0.65 : %lld\n",writebuf); fflush(fail_file);
#endif
      MPI_Irecv(writebuf,blinkptr->num,datatype,MPIid[0],myid, MPI_COMM_GRMHD,&request);      
#if(DEBUGMINMEM)
      dualfprintf(fail_file,"got here 0.66 : %lld %lld\n",writebuf,blinkptr); fflush(fail_file);
#endif
      bufferoffset=(long long int)(-datagot);
    }// end if READFILE



    datagot += (long long int)(blinkptr->num);
    lastnum = (long long int)(blinkptr->num);
    // now iterate to next node in list
    blinkptr=(blinkptr->np);
#if(DEBUGMINMEM)
    dualfprintf(fail_file,"got here 0.67 : %lld\n",blinkptr); fflush(fail_file);
#endif
    if(blinkptr==NULL){
      thisslavedone=1;
      // then done with list, check to see if everything adds up
      if(
         ((readwrite==WRITEFILE)&&(-bufferoffset!=(long long int)N1*(long long int)N2*(long long int)N3*(long long int)numcolumns))
         ||((readwrite==READFILE)&&(-bufferoffset!=(long long int)N1*(long long int)N2*(long long int)N3*(long long int)numcolumns-(long long int)lastnum))
         )
        {
          dualfprintf(fail_file,"local total doesn't equal expected value\n got=%d demand=%d\n",-bufferoffset,(long long int)N1*(long long int)N2*(long long int)N3*(long long int)numcolumns);
          myexit(10002);
        }
    }
#if(DEBUGMINMEM)
    dualfprintf(fail_file,"got here 0.7\n"); fflush(fail_file);
#endif

  }

#if(DEBUGMINMEM)
  dualfprintf(fail_file,"got here .8\n"); fflush(fail_file);
#endif





  
  // for cpu=0, don't wait on local send/recv since cpu=0 needs to setup it's receives first (since cpu=0 sends to itself)
  // only do this if cpu=0 is just done with local data send/recv(i.e. dolocal==1) or cpu=0 is done (done0==1)
  // no way to decouple cpu=0 local from master, so master waits for local
  if(myid==0){
    // first time hit, get initial node in list and reset counters
    if((i==0)&&(j==0)&&(k==0)){
      cpulinkptr=cpulinkptr0[whichdump];
      datasofar0=0;
      datasofarc0=0;
    }
  }




  // check to see if done with cpu=0 data (we use cpu=0 to continue grabbing data as needed)
  if(readwrite==WRITEFILE){
    // wait until all data is grabbed from cpu=0 to send to disk before holding up here.
    if((myid==0)&&(bfi==(long long int)N1*(long long int)N2*(long long int)N3*(long long int)numcolumns)){
      done0=1;
    }
    else done0=0;
  }
  else if(readwrite==READFILE){
    // wait until last read done.  Still can't use writebuf since will do last local cpu=0 assignment using writebuf after done with distribution of rest of data using cpu=0 to other cpus
    if((myid==0)&&(thisslavedone)){
      done0=1; // then will be done after final read and grab for cpu=0 (then cpu=0's writebuf will be written to after all cpus done)
    }
    else done0=0;
  }







  ///////////////////////
  //
  // CPU==0 stuff
  //
  ////////////////////////


  if((myid==0)&&((dolocal==1)||(done0==1))){
    // We are ready to read or write from/to other cpus

    if (datatype == MPI_UNSIGNED_CHAR) jonio1 = (unsigned char *) jonio;
    else if (datatype == MPI_FLOAT) jonio4 = (float *) jonio;
    else if (datatype == MPI_DOUBLE) jonio8 = (double *) jonio;
    else if (datatype == MPI_LONG_DOUBLE) jonio16 = (long double *) jonio;
    else if (datatype == MPI_INT) jonio4i = (int *) jonio;
    else if (datatype == MPI_LONG_LONG_INT) jonio8i = (long long int *) jonio;



    // we need to loop over list until hit a marked end-cpu.  This will give us a chunk of data to sort for direct output
    // this is the offset to the sorted part of jonio
    joniooffset=(long long int)joniosize/((long long int)2);

#if(DEBUGMINMEM)
    dualfprintf(fail_file,"got here 2: joniooffset=%lld\n",joniooffset); fflush(fail_file);
#endif


    /////////////////////////////
    //
    // WRITEFILE
    //
    /////////////////////////////


    if(readwrite==WRITEFILE){
      dofull=1;
      while(dofull){


#if(DEBUGMINMEM)
        dualfprintf(fail_file,"got here 3\n"); fflush(fail_file);
#endif
        /////////
        // determine amount of data to get from cpu group
 
        // limit the amount of data once doing last grab
        if((long long int)datasofar0+(long long int)joniosize/((long long int)2)>(long long int)totalsize[1]*(long long int)totalsize[2]*(long long int)totalsize[3]*(long long int)numcolumns) datatogetc0=-(long long int)datasofar0+(long long int)totalsize[1]*(long long int)totalsize[2]*(long long int)totalsize[3]*(long long int)numcolumns;
        else datatogetc0=joniosize/((long long int)2);
        // new amount of data (that will be read total)
        datasofarc0 += (long long int)datatogetc0;


        doset=1;
        doing0=0;
        // every receieved dataset is kept in first half of jonio from start of jonio, sorted, then next data set received.
        datatoget0=0;
        while(doset){
#if(DEBUGMINMEM)
          dualfprintf(fail_file,"got here 3.4 : cpu=%d: datasofar0=%lld  datasofarc0=%lld\n",cpulinkptr->cpu,datasofar0,datasofarc0); fflush(fail_file);
#endif
          datatoget0 += (long long int)(cpulinkptr->num);
          if(cpulinkptr->cpu==0) doing0=1;
#if(DEBUGMINMEM)
          dualfprintf(fail_file,"got here 3.5: %lld %d %d\n",jonio,cpulinkptr->num,cpulinkptr->cpu); fflush(fail_file);
          dualfprintf(fail_file,"got here 3.51: %lld %lld %lld %lld\n",datatogetc0,totalsize[1],totalsize[2],totalsize[3]); fflush(fail_file);
#endif
          MPI_Irecv(jonio,cpulinkptr->num,datatype,MPIid[cpulinkptr->cpu],cpulinkptr->cpu,MPI_COMM_GRMHD,&request0);
#if(DEBUGMINMEM)
          dualfprintf(fail_file,"got here 3.515\n"); fflush(fail_file);
#endif
          // have myid wait before continuing to make sure receive is complete
          MPI_Wait(&request0,&mpichstatus);
          // easiest algorithm/mapping is done using loop over full sorted size, reduced per cpu by if statement and checked
          //   for(sii=0;sii<cpulinkptr->num;sii++){

#if(DEBUGMINMEM)
          dualfprintf(fail_file,"got here 3.52: datatogetc0=%lld\n",datatogetc0); fflush(fail_file);
#endif

#if(DEBUGMINMEM)
          dualfprintf(fail_file,"got here 3.54: %d %d %d\n",cpulinkptr->ri,cpulinkptr->rj,cpulinkptr->rk); fflush(fail_file);
#endif


          // repeat this loop until really get all of datatogetc0 from (multiple) CPUs
          // The if(gi,gj,gk) below constrains the filling of jonio so that fills in global array with required data into correct array positions
          uii=0;
          for(sii=0;sii<datatogetc0;sii++){
            // uii: iteration count for accessing memory from received data

            // sii: iteration count for accessing (with joniooffset) where to store some of received data

            // ri : starting global (over all CPUs) i-location
            // rj : "" for j-location
            // rk : "" for k-location

            // sii : sorted index that iterates to include column data

            // mypos: number of grid positions (not including columns) iterated beyond starting (ri,rj,rk) position

            // gi : global (over all CPUs) i-location
            // gj : "" for j-location
            // gk : "" for k-location

            mypos=(long long int)(sii/((long long int)numcolumns)) + (long long int)(cpulinkptr->ri) + (long long int)(cpulinkptr->rj)*(long long int)totalsize[1] + (long long int)(cpulinkptr->rk)*(long long int)totalsize[1]*(long long int)totalsize[2];

            gi=(long long int)mypos%((long long int)totalsize[1]);
            gj=(long long int)((mypos%((long long int)totalsize[1]*(long long int)totalsize[2]))/((long long int)totalsize[1]));
            gk=(long long int)(mypos/((long long int)totalsize[1]*(long long int)totalsize[2]));

#if(DEBUGMINMEM)
            dualfprintf(fail_file,"got here 3.55: sii=%lld mypos=%lld  gi=%lld gj=%lld gk=%lld\n",sii, mypos, gi,gj,gk); fflush(fail_file);
#endif

            if(
               (gi>=startpos0[1][cpulinkptr->cpu])&&
               (gi<=  endpos0[1][cpulinkptr->cpu])&&
               (gj>=startpos0[2][cpulinkptr->cpu])&&
               (gj<=  endpos0[2][cpulinkptr->cpu])&&
               (gk>=startpos0[3][cpulinkptr->cpu])&&
               (gk<=  endpos0[3][cpulinkptr->cpu])
               ){

#if(DEBUGMINMEM)
              dualfprintf(fail_file,"got here 3.56: did assign: sii=%lld joniooffset=%lld uii=%lld\n",sii,joniooffset,uii); fflush(fail_file);
#endif

              if (datatype == MPI_UNSIGNED_CHAR) jonio1[sii+joniooffset]=jonio1[uii++];
              else if (datatype == MPI_FLOAT) jonio4[sii+joniooffset]=jonio4[uii++];
              else if (datatype == MPI_DOUBLE) jonio8[sii+joniooffset]=jonio8[uii++];
              else if (datatype == MPI_LONG_DOUBLE) jonio16[sii+joniooffset]=jonio16[uii++];
              else if (datatype == MPI_INT) jonio4i[sii+joniooffset]=jonio4i[uii++];
              else if (datatype == MPI_LONG_LONG_INT) jonio8i[sii+joniooffset]=jonio8i[uii++];
            }
          }// end over sii


#if(DEBUGMINMEM)
          dualfprintf(fail_file,"got here 3.6\n"); fflush(fail_file);
#endif

          /////////////////
          //
          // check!
          // see if local data total is as expected
          //
          ////////////////
          if(uii!=(long long int)(cpulinkptr->num)){
            dualfprintf(fail_file,"uii and num for this cpu not same, algorithm failure: uii=%d num=%d datatogetc0=%d\n",uii,cpulinkptr->num,datatogetc0);
            myexit(10003);
          }
          if(cpulinkptr->end) doset=0;
          cpulinkptr=cpulinkptr->np;


        }// end while(doset)



        // see if total data is as expected
        datasofar0 += (long long int)datatoget0; // diagnostic
        if((long long int)datasofar0 != (long long int)datasofarc0){
          dualfprintf(fail_file,"cumulative data received via MPI and expected data is different: got=%d expected=%d\n",datasofar0,datasofarc0);
          myexit(10004);
        }



        // now take sorted part and write it to disk
        // now write out collected data using CPU=0


        /////////////
        //
        // write data to file
        //
        ///////////////
        if (datatype == MPI_UNSIGNED_CHAR) joniosubmit=(void*) (jonio1+joniooffset);
        else if (datatype == MPI_FLOAT) joniosubmit=(void*) (jonio4+joniooffset);
        else if (datatype == MPI_DOUBLE) joniosubmit=(void*) (jonio8+joniooffset);
        else if (datatype == MPI_LONG_DOUBLE) joniosubmit=(void*) (jonio16+joniooffset);
        else if (datatype == MPI_INT) joniosubmit=(void*) (jonio4i+joniooffset);
        else if (datatype == MPI_LONG_LONG_INT) joniosubmit=(void*) (jonio8i+joniooffset);

        if(docolsplit){
          numfiles=numcolumns;
        }
        else numfiles=1;

        if(bintxt==BINARYOUTPUT){
          if(numfiles==1) fwrite(joniosubmit, sizeofdatatype,datatoget0, *fp);
          else{
            for(sii=0;sii<datatoget0;sii++){
              if (datatype == MPI_UNSIGNED_CHAR) fwrite(&jonio1[sii+joniooffset], sizeofdatatype,1, fp[sii%numfiles]);
              else if (datatype == MPI_FLOAT) fwrite(&jonio4[sii+joniooffset], sizeofdatatype,1, fp[sii%numfiles]);
              else if (datatype == MPI_DOUBLE) fwrite(&jonio8[sii+joniooffset], sizeofdatatype,1, fp[sii%numfiles]);
              else if (datatype == MPI_LONG_DOUBLE) fwrite(&jonio16[sii+joniooffset], sizeofdatatype,1, fp[sii%numfiles]);
              else if (datatype == MPI_INT) fwrite(&jonio4i[sii+joniooffset], sizeofdatatype,1, fp[sii%numfiles]);
              else if (datatype == MPI_LONG_LONG_INT) fwrite(&jonio8i[sii+joniooffset], sizeofdatatype,1, fp[sii%numfiles]);
            }
          }
        }
        else if(bintxt==TEXTOUTPUT){ // properly ordered, so just dump it
          for(sii=0;sii<datatoget0;sii++){
            if (datatype == MPI_UNSIGNED_CHAR) fprintf(fp[sii%numfiles],"%c",jonio1[sii+joniooffset]);
            else if (datatype == MPI_FLOAT) fprintf(fp[sii%numfiles],"%15.7g",jonio4[sii+joniooffset]);
            else if (datatype == MPI_DOUBLE) fprintf(fp[sii%numfiles],"%21.15g",jonio8[sii+joniooffset]);
            else if (datatype == MPI_LONG_DOUBLE) fprintf(fp[sii%numfiles],"%31.25Lg",jonio16[sii+joniooffset]);
            else if (datatype == MPI_INT) fprintf(fp[sii%numfiles],"%d",jonio4i[sii+joniooffset]);
            else if (datatype == MPI_LONG_LONG_INT) fprintf(fp[sii%numfiles],"%lld",jonio8i[sii+joniooffset]);
            if(numfiles==1){
              if((sii+1)%numcolumns) fprintf(fp[sii%numfiles]," ");
              else fprintf(fp[sii%numfiles],"\n");
            }
            else fprintf(fp[sii%numfiles],"\n");
          }
        }
        if((done0==0)&&(doing0==1)) dofull=0; // cpu==0 still needs to deal with it's own data
        // if wasn't locally doing cpu=0 by master, then can't exit yet, continue till cpu=0 local data needed by master
        if(cpulinkptr==NULL){
          dofull=0; // end of list
          masterdone=1;
          if(datasofar0 != (long long int)totalsize[1]*(long long int)totalsize[2]*(long long int)totalsize[3]*(long long int)numcolumns){
            dualfprintf(fail_file,"write: total data written not equal to expected amount: wrote=%lld expected=%lld\n",datasofar0,(long long int)totalsize[1]*(long long int)totalsize[2]*(long long int)totalsize[3]*(long long int)numcolumns);
            myexit(10005);
          }
        }
        // otherwise continue
      }



    }
    else if(readwrite==READFILE){


      ///////////////////////////
      //
      // READFILE
      //
      ////////////////////////////


      dofull=1;
      while(dofull){
#if(DEBUGMINMEM)
        dualfprintf(fail_file,"got here 4\n"); fflush(fail_file);
#endif
        // we read data into 2nd half of jonio (sorted part), then desort for one cpu, then continue for next cpu 
 
        /////////
        // determine amount of data to get from file


        // limit the amount of data once doing last grab
        if(datasofar0+(long long int)joniosize/((long long int)2) > (long long int)totalsize[1]*(long long int)totalsize[2]*(long long int)totalsize[3]*(long long int)numcolumns) datatogetc0 = -(long long int)datasofar0 + (long long int)totalsize[1]*(long long int)totalsize[2]*(long long int)totalsize[3]*(long long int)numcolumns;
        else datatogetc0=joniosize/2;
        // new amount of data (that will be read total)
        datasofarc0+=datatogetc0;      

#if(DEBUGMINMEM)
        dualfprintf(fail_file,"got here1 : %lld %lld\n",datatogetc0,datasofarc0); fflush(fail_file);
#endif

        /////////////
        // get data from file
        if (datatype == MPI_UNSIGNED_CHAR) joniosubmit=(void*) (jonio1+joniooffset);
        else if (datatype == MPI_FLOAT) joniosubmit=(void*) (jonio4+joniooffset);
        else if (datatype == MPI_DOUBLE) joniosubmit=(void*) (jonio8+joniooffset);
        else if (datatype == MPI_LONG_DOUBLE) joniosubmit=(void*) (jonio16+joniooffset);
        else if (datatype == MPI_INT) joniosubmit=(void*) (jonio4i+joniooffset);
        else if (datatype == MPI_LONG_LONG_INT) joniosubmit=(void*) (jonio8i+joniooffset);
 

        if(docolsplit){
          numfiles=numcolumns;
        }
        else numfiles=1;

        if(bintxt==BINARYOUTPUT){
          // first let cpu=0 read data
          if(numfiles==1) fread(joniosubmit, sizeofdatatype,datatogetc0,*fp);
          else{
            for(sii=0;sii<datatoget0;sii++){
              if (datatype == MPI_UNSIGNED_CHAR) fread(&jonio1[sii+joniooffset], sizeofdatatype,1, fp[sii%numfiles]);
              else if (datatype == MPI_FLOAT) fread(&jonio4[sii+joniooffset], sizeofdatatype,1, fp[sii%numfiles]);
              else if (datatype == MPI_DOUBLE) fread(&jonio8[sii+joniooffset], sizeofdatatype,1, fp[sii%numfiles]);
              else if (datatype == MPI_LONG_DOUBLE) fread(&jonio16[sii+joniooffset], sizeofdatatype,1, fp[sii%numfiles]);
              else if (datatype == MPI_INT) fread(&jonio4i[sii+joniooffset], sizeofdatatype,1, fp[sii%numfiles]);
              else if (datatype == MPI_LONG_LONG_INT) fread(&jonio8i[sii+joniooffset], sizeofdatatype,1, fp[sii%numfiles]);
            }
          }
        }
        else if(bintxt==TEXTOUTPUT){ // properly ordered, so just grab it
          for(sii=0;sii<datatogetc0;sii++){
            if (datatype == MPI_UNSIGNED_CHAR){
              fscanf(fp[sii%numfiles],"%hu",&short4char);
              jonio1[sii+joniooffset]=short4char; // convert to unsigned char
            }
            else if (datatype == MPI_FLOAT) fscanf(fp[sii%numfiles],"%f",&jonio4[sii+joniooffset]);
            else if (datatype == MPI_DOUBLE) fscanf(fp[sii%numfiles],"%lf",&jonio8[sii+joniooffset]);
            else if (datatype == MPI_LONG_DOUBLE) fscanf(fp[sii%numfiles],"%Lf",&jonio16[sii+joniooffset]);
            else if (datatype == MPI_INT) fscanf(fp[sii%numfiles],"%d",&jonio4i[sii+joniooffset]);
            else if (datatype == MPI_LONG_LONG_INT) fscanf(fp[sii%numfiles],"%lld",&jonio8i[sii+joniooffset]);
          }
        }// end if TEXTOUTPUT
#if(DEBUGMINMEM)
        dualfprintf(fail_file,"got here2\n"); fflush(fail_file);
#endif



 
        /////////////////
        //
        // send data to other cpus (while(doset))
        //
        //////////////////
 
        // now that data is read into 2nd half of jonio, need to desort data into 1st half per cpu in a loop
        datasent0=0;
        datatoget0=0;
        doset=1;
        doing0=0;

        while(doset){
#if(DEBUGMINMEM)
          dualfprintf(fail_file,"got here 2.4\n"); fflush(fail_file);
#endif

          if(cpulinkptr->cpu==0) doing0=1;
          datatoget0 += (long long int)(cpulinkptr->num);

          // copy from 2nd half of jonio to first half the node cpu's data
          uii=0;
          for(sii=0;sii<datatogetc0;sii++){
#if(DEBUGMINMEM)
            dualfprintf(fail_file,"got here2.5: %lld\n",sii); fflush(fail_file);
#endif

            mypos=(long long int)((long long int)sii/((long long int)numcolumns)) + (long long int)(cpulinkptr->ri) + (long long int)(cpulinkptr->rj)*(long long int)totalsize[1] + ((long long int)cpulinkptr->rk)*(long long int)totalsize[1]*(long long int)totalsize[2];

            gi=((long long int)mypos)%((long long int)totalsize[1]);
            gj=(long long int)((((long long int)mypos)%((long long int)totalsize[1]*(long long int)totalsize[2]))/((long long int)totalsize[1]));
            gk=(long long int)(((long long int)mypos)/((long long int)totalsize[1]*(long long int)totalsize[2]));

#if(DEBUGMINMEM)
            dualfprintf(fail_file,"got here2.6: %lld %lld\n",gi,gj); fflush(fail_file);
#endif
            if(
               (gi>=startpos0[1][cpulinkptr->cpu])&&
               (gi<=endpos0[1][cpulinkptr->cpu])&&
               (gj>=startpos0[2][cpulinkptr->cpu])&&
               (gj<=endpos0[2][cpulinkptr->cpu])&&
               (gk>=startpos0[3][cpulinkptr->cpu])&&
               (gk<=endpos0[3][cpulinkptr->cpu])
               ){
#if(DEBUGMINMEM)
              dualfprintf(fail_file,"got here2.7 %lld\n",cpulinkptr->cpu); fflush(fail_file);
#endif

              if (datatype == MPI_UNSIGNED_CHAR) jonio1[uii++]=jonio1[sii+joniooffset];
              else if (datatype == MPI_FLOAT) jonio4[uii++]=jonio4[sii+joniooffset];
              else if (datatype == MPI_DOUBLE) jonio8[uii++]=jonio8[sii+joniooffset];
              else if (datatype == MPI_LONG_DOUBLE) jonio16[uii++]=jonio16[sii+joniooffset];
              else if (datatype == MPI_INT) jonio4i[uii++]=jonio16[sii+joniooffset];
              else if (datatype == MPI_LONG_LONG_INT) jonio8i[uii++]=jonio16[sii+joniooffset];
            }// end if within this CPUs data
          }// end over sii
#if(DEBUGMINMEM)
          dualfprintf(fail_file,"got here3\n"); fflush(fail_file);
#endif


          /////////////////////
          //
          // check!
          if(uii != (long long int)(cpulinkptr->num)){
            dualfprintf(fail_file,"uii and num for this cpu not same, algorithm failure: uii=%d num=%d datatogetc0=%d\n",uii,cpulinkptr->num,datatogetc0);
            myexit(10006);
          }
#if(DEBUGMINMEM)
          dualfprintf(fail_file,"got here4\n"); fflush(fail_file);   
#endif
          // jonio is the unsorted bit here starting at index=0 (for all cpus)
          MPI_Isend(jonio,cpulinkptr->num,datatype,MPIid[cpulinkptr->cpu],cpulinkptr->cpu,MPI_COMM_GRMHD,&request0);
          // have myid wait before continuing to make sure receive is complete
          // alternative is to have many sends, but then need to desort all at once which isn't easy since cycling through link list one at a time only once (could store starter pointer, etc.)
          MPI_Wait(&request0,&mpichstatus);
   
          datasent0 += (long long int)(cpulinkptr->num); // diagnostic

          if(cpulinkptr->end) doset=0;
          cpulinkptr=cpulinkptr->np;
#if(DEBUGMINMEM)
          dualfprintf(fail_file,"got here5\n"); fflush(fail_file);   
#endif

        }// end while doset





        ///////////////////
        //
        // check to see if total data is correct
        //
        //////////////////

#if(DEBUGMINMEM)
        dualfprintf(fail_file,"got here6\n"); fflush(fail_file);   
#endif

        datasofar0 += (long long int)datatoget0; // diagnostic
        if((long long int)datasofar0 != (long long int)datasofarc0){
          dualfprintf(fail_file,"cumulative data received via MPI and expected data is different: got=%d expected=%d\n",datasofar0,datasofarc0);
          myexit(10007);
        }
        if((long long int)datasent0 != (long long int)datatoget0){
          dualfprintf(fail_file,"data sent doesn't match data read\n");
          myexit(10008);
        }
        if((done0==0)&&(doing0==1)) dofull=0; // cpu==0 still needs to deal with more reads for it's own data
        if(cpulinkptr==NULL){
          masterdone=1; // i.e. don't come back
          dofull=0; // end of list
          if(datasofar0 != (long long int)totalsize[1]*(long long int)totalsize[2]*(long long int)totalsize[3]*(long long int)numcolumns){
            dualfprintf(fail_file,"read: total data written not equal to expected amount: wrote=%lld expected=%lld\n",datasofar0,(long long int)totalsize[1]*(long long int)totalsize[2]*(long long int)totalsize[3]*(long long int)numcolumns);
            myexit(10009);
          }
        }// end if NULL




        // otherwise continue
#if(DEBUGMINMEM)
        dualfprintf(fail_file,"got here7\n"); fflush(fail_file);   
#endif



      }// end while dofull




#if(DEBUGMINMEM)
      dualfprintf(fail_file,"got here8\n"); fflush(fail_file);   
#endif



    }// end if READFILE





  }// end if myid==0 && ((dolocal==0)||(doone==1))



#if(DEBUGMINMEM)
  dualfprintf(fail_file,"got here 6\n"); fflush(fail_file);
#endif

  // have myid wait before continuing so buffer can be released for more writing to
  if(dolocal) MPI_Wait(&request,&mpichstatus); // myid==0 handled specially since also master, and can't wait, even if master done first in this function, since master may not require cpu=0 at all for its current node set 

  //  logsfprintf("end mpiminmem_read\n");
#if(DEBUGMINMEM)
  dualfprintf(fail_file,"got here9\n"); fflush(fail_file);   
#endif


#endif
}






void mpiio_seperate(int bintxt, int sorted, int stage,
                    int numcolumns, MPI_Datatype datatype,
                    FILE ** fpptr, void *jonio, void *writebuf)
{
#if(USEMPI)

  logsfprintf("mpiio begin seperate\n");

  if(truempicombinetype==MPICOMBINEROMIO){
    // doesn't use jonioptr or bintxt or sorted
    // headerbytesize no longer needed
    if(stage==STAGE1) mpiioromio_init_combine(READROMIO, READFILE,  0, "", numcolumns, datatype,&writebuf,writebuf);
    if(stage==STAGE2) mpiioromio_init_combine(READFREEROMIO, READFILE,  0, "", numcolumns, datatype,&writebuf,writebuf);
  }
  else{

    if(truempicombinetype==MPICOMBINESIMPLE){
      
      if(sorted){
        mpiios_seperate(bintxt, stage, datatype, numcolumns, fpptr, jonio, writebuf);
      }
      else{
        mpiiotu_seperate(stage, datatype, numcolumns, fpptr, writebuf);
      }
      
      
    }
    else if(truempicombinetype==MPICOMBINEMINMEM){
      // do nothing on STAGE1
      if(stage==STAGE2) mpiiomin_final(numcolumns,fpptr, jonio, writebuf);
    }
  }

#endif
  logsfprintf("mpiio end seperate\n");

}




/// NICS Nautilus hints:
/// http://www.nics.tennessee.edu/user-support/mpi-tips-for-cray-xt5
/// http://www.nics.tennessee.edu/io-tips
/// filename only needed for INITROMIO
/// operationtype=INITROMIO or WRITECLOSEROMIO
/// mpi io using ROMIO
/// which=WRITEFILE or READFILE
void mpiioromio_init_combine(int operationtype, int which,  long headerbytesize, char *filename, int numcolumns,MPI_Datatype datatype, void **writebufptr,void *writebuf)
{
  int i;

  static long long int sizeofmemory;
  static int sizeofdatatype;

  static long double **writebuf16;
  static double **writebuf8;
  static float **writebuf4;
  static unsigned char **writebuf1;
  static int **writebuf4i;
  static long long int **writebuf8i;

  static int numfiles;
  static int romiocolumns;
  static char realfilename[MAXFILENAME];
#if(USEMPI&&USEROMIO)
  static MPI_Datatype newtype;
  static MPI_File fh;
  static MPI_Status status;
  static MPI_Request request;
#endif
  static int ndims, array_of_gsizes[4], array_of_distribs[4];
  static int order, len;
  static int array_of_dargs[4], array_of_psizes[4];
  static int bufcount, array_size;






  if(operationtype==INITROMIO){

    sizeofdatatype=getsizeofdatatype(datatype);

    logsfprintf("mpiioromio_init begin\n");


    //////////////////////////
    //
    // set dimensionality

    // see if splitting into individual columns
    if(docolsplit){
      numfiles=numcolumns;
      romiocolumns=1;
    }
    else{
      numfiles=1;
      romiocolumns=numcolumns;
    }


#if(USEMPI&&USEROMIO)
    if((COMPDIM==3)&&(romiocolumns>1)){
      //create the distributed array filetype
      // this still works if any reduced dimensionality
      ndims = 4;
      order = MPI_ORDER_C;
    
      array_of_gsizes[3] = romiocolumns;
      array_of_gsizes[2] = totalsize[1];
      array_of_gsizes[1] = totalsize[2];
      array_of_gsizes[0] = totalsize[3];

      sizeofmemory = (long long int)N1*(long long int)N2*(long long int)N3*(long long int)romiocolumns*(long long int)sizeofdatatype;
    
      array_of_distribs[3] = MPI_DISTRIBUTE_BLOCK;
      array_of_distribs[2] = MPI_DISTRIBUTE_BLOCK;
      array_of_distribs[1] = MPI_DISTRIBUTE_BLOCK;
      array_of_distribs[0] = MPI_DISTRIBUTE_BLOCK;
    
      array_of_dargs[3] = MPI_DISTRIBUTE_DFLT_DARG;
      array_of_dargs[2] = MPI_DISTRIBUTE_DFLT_DARG;
      array_of_dargs[1] = MPI_DISTRIBUTE_DFLT_DARG;
      array_of_dargs[0] = MPI_DISTRIBUTE_DFLT_DARG;
    
      array_of_psizes[3]=1;
      array_of_psizes[2]=ncpux1;
      array_of_psizes[1]=ncpux2;
      array_of_psizes[0]=ncpux3;
    }
    else if((COMPDIM==2)&&(romiocolumns>1)){
      //create the distributed array filetype
      ndims = 3;
      order = MPI_ORDER_C;
    
      array_of_gsizes[2] = romiocolumns;
      array_of_gsizes[1] = totalsize[1];
      array_of_gsizes[0] = totalsize[2];

      sizeofmemory = (long long int)N1*(long long int)N2*(long long int)romiocolumns*(long long int)sizeofdatatype;
    
      array_of_distribs[2] = MPI_DISTRIBUTE_BLOCK;
      array_of_distribs[1] = MPI_DISTRIBUTE_BLOCK;
      array_of_distribs[0] = MPI_DISTRIBUTE_BLOCK;
    
      array_of_dargs[2] = MPI_DISTRIBUTE_DFLT_DARG;
      array_of_dargs[1] = MPI_DISTRIBUTE_DFLT_DARG;
      array_of_dargs[0] = MPI_DISTRIBUTE_DFLT_DARG;
    
      array_of_psizes[2]=1;
      array_of_psizes[1]=ncpux1;
      array_of_psizes[0]=ncpux2;
    }
    else if((COMPDIM==1)&&(romiocolumns>1)){
      //create the distributed array filetype
      ndims = 1;
      order = MPI_ORDER_C;
    
      array_of_gsizes[1] = romiocolumns;
      array_of_gsizes[0] = totalsize[1];

      sizeofmemory = (long long int)N1*(long long int)romiocolumns*(long long int)sizeofdatatype;
    
      array_of_distribs[1] = MPI_DISTRIBUTE_BLOCK;
      array_of_distribs[0] = MPI_DISTRIBUTE_BLOCK;
    
      array_of_dargs[1] = MPI_DISTRIBUTE_DFLT_DARG;
      array_of_dargs[0] = MPI_DISTRIBUTE_DFLT_DARG;
    
      array_of_psizes[1]=1;
      array_of_psizes[0]=ncpux1;
    }
    else if((COMPDIM==3)&&(romiocolumns==1)){
      //create the distributed array filetype
      // this still works if any reduced dimensionality
      ndims=3;
      order = MPI_ORDER_C;
      
      array_of_gsizes[2] = totalsize[1];
      array_of_gsizes[1] = totalsize[2];
      array_of_gsizes[0] = totalsize[3];

      sizeofmemory = (long long int)N1*(long long int)N2*(long long int)N3*(long long int)sizeofdatatype;
      
      array_of_distribs[2] = MPI_DISTRIBUTE_BLOCK;
      array_of_distribs[1] = MPI_DISTRIBUTE_BLOCK;
      array_of_distribs[0] = MPI_DISTRIBUTE_BLOCK;
      
      array_of_dargs[2] = MPI_DISTRIBUTE_DFLT_DARG;
      array_of_dargs[1] = MPI_DISTRIBUTE_DFLT_DARG;
      array_of_dargs[0] = MPI_DISTRIBUTE_DFLT_DARG;
      
      array_of_psizes[2]=ncpux1;
      array_of_psizes[1]=ncpux2;
      array_of_psizes[0]=ncpux3;
    }
    else  if((COMPDIM==2)&&(romiocolumns==1)){
      //create the distributed array filetype
      ndims=2;
      order = MPI_ORDER_C;
      
      array_of_gsizes[1] = totalsize[1];
      array_of_gsizes[0] = totalsize[2];

      sizeofmemory = (long long int)N1*(long long int)N2*(long long int)sizeofdatatype;
      
      array_of_distribs[1] = MPI_DISTRIBUTE_BLOCK;
      array_of_distribs[0] = MPI_DISTRIBUTE_BLOCK;
      
      array_of_dargs[1] = MPI_DISTRIBUTE_DFLT_DARG;
      array_of_dargs[0] = MPI_DISTRIBUTE_DFLT_DARG;
      
      array_of_psizes[1]=ncpux1;
      array_of_psizes[0]=ncpux2;
    }
    else  if((COMPDIM==1)&&(romiocolumns==1)){
      //create the distributed array filetype
      ndims=1;
      order = MPI_ORDER_C;
      
      array_of_gsizes[0] = totalsize[1];

      sizeofmemory = (long long int)N1*(long long int)sizeofdatatype;
      
      array_of_distribs[0] = MPI_DISTRIBUTE_BLOCK;
      
      array_of_dargs[0] = MPI_DISTRIBUTE_DFLT_DARG;
      
      array_of_psizes[0]=ncpux1;
    }
    else if(romiocolumns==0){
      // write nothing actually
      ndims=1;
      order = MPI_ORDER_C;
      
      array_of_gsizes[0] = 1;

      sizeofmemory = 0; // 1*(long long int)sizeofdatatype;
      
      array_of_distribs[0] = MPI_DISTRIBUTE_BLOCK;
      
      array_of_dargs[0] = MPI_DISTRIBUTE_DFLT_DARG;
      
      array_of_psizes[0]=1;
    }
    else{
      dualfprintf(fail_file,"Shouldn't reach to end of ROMIO selection\n");
      myexit(17872);
    }
#endif
      
    // initialize filename
    // must be same filename as in dump_gen() when opening files here, so header is included at top of the dump
    if(docolsplit&&(numcolumns>1)){
      sprintf(realfilename,"%s-col%04d",filename,romiocoliter);
    }
    else strcpy(realfilename,filename);
     


    //////////////////////////
    //
    // initialize grid and file handler for ROMIO

#if(USEMPI&&USEROMIO)
    // NOT using mapping MPIid[myid] below.  Using GRMHD layout rank "myid" that ensures that ROMIO writes correctly spatially no matter what MPI rank is.
    MPI_Type_create_darray(numprocs, myid, ndims, array_of_gsizes, 
                           array_of_distribs, array_of_dargs,
                           array_of_psizes, order, datatype, &newtype);
    MPI_Type_commit(&newtype);

#if(0)
    MPI_Type_size(newtype, &bufcount); // ULTRASUPERGODMARK: bufcount is int limited to 2GB, while should be allowed to be larger.  Couldn't find info online about how to fix this and suggestions that 2GB is max buffer size of any MPI communication!  That's stupid, since ruins ROMIO ability.
    sizeofmemory=bufcount; //includes type
    bufcount = bufcount/sizeofdatatype; // number of elements of type
    // could avoid use of MPI_Type_size() and set sizeofmemory and bufcount directly, but still bufcount must be integer when used in MPI functions below, so only helps by sizeofdatatype in size
#else
    // then sizeofmemory already set
    bufcount=(long long int)sizeofmemory/(long long int)sizeofdatatype;
#endif


    // fail if MPI functions below can't handle buffer size
    if((long long int)N1*(long long int)N2*(long long int)N3*(long long int)romiocolumns*(long long int)sizeofdatatype>=(long long int)(2L*1024L*1024L*1024L) && sizeof(int)<=4){
      dualfprintf(fail_file,"JCM couldn't figure out how to modify ROMIO so would work when sizeof(int)==4 and buffer size is >2GB\n");
      myexit(867546243);
    }

    // setup file handler
    // setup to append, in case wanted binary header at front
    if(which==WRITEFILE){
      MPI_File_open(MPI_COMM_GRMHD, realfilename, MPI_MODE_APPEND | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
    }
    else if(which==READFILE){
      MPI_File_open(MPI_COMM_GRMHD, realfilename, MPI_MODE_RDONLY , MPI_INFO_NULL, &fh);
    }

    // this sets the distributed nature of the file writing
    // http://www.mpi-forum.org/docs/mpi-20-html/node198.htm#Node198 : "naive" "internal" "external32"
    // GODMARK: rather than "native" can use MPI_Register_datarep() to create a conversion between "native" and other types for reading/writing
    // but apparently can't produce formatted text output
    // MPI_Datarep_conversion_function() used too
    // http://www.mpi-forum.org/docs/mpi-20-html/node201.htm#Node201

    // Note that for 64-bit  machines "native" and "external32" are not compatible

    // GODMARK: For now use "native" since "external32" not supported even as of 07/21/08 by MPICH2 or OpenMPI
    MPI_File_set_view(fh, headerbytesize, datatype, newtype, "native", MPI_INFO_NULL);
    // now use "external32" so files are exactly well-defined for any system as in:
    // http://parallel.ru/docs/Parallel/mpi2/node200.html#Node200
    // perhaps my bin2txt program should use these MPI types too
    //
    // GODMARK: Avoid the below good feature since not supported in general
    //    MPI_File_set_view(fh, headerbytesize, datatype, newtype, "external32", MPI_INFO_NULL);
    // all that needs to be done now is initialize and later fill writebuf with the data
#endif

    //////////////////////////
    //
    // allocate memory for writebuf
   

    
    if(datatype==MPI_UNSIGNED_CHAR){ writebuf1=(unsigned char **)writebufptr; }
    else if(datatype==MPI_FLOAT){  writebuf4=(float **)writebufptr; }
    else if(datatype==MPI_DOUBLE){  writebuf8=(double**)writebufptr; }
    else if(datatype==MPI_LONG_DOUBLE){  writebuf16=(long double**)writebufptr; }
    else if(datatype==MPI_INT){  writebuf4i=(int **)writebufptr; }
    else if(datatype==MPI_LONG_LONG_INT){  writebuf8i=(long long int **)writebufptr; }


    if(datatype==MPI_UNSIGNED_CHAR) *writebuf1=(unsigned char*)malloc(sizeofmemory);
    else if(datatype==MPI_FLOAT) *writebuf4=(float*)malloc(sizeofmemory);
    else if(datatype==MPI_DOUBLE) *writebuf8 =(double*)malloc(sizeofmemory);
    else if(datatype==MPI_LONG_DOUBLE) *writebuf16 =(long double*)malloc(sizeofmemory);
    else if(datatype==MPI_INT) *writebuf4i=(int*)malloc(sizeofmemory);
    else if(datatype==MPI_LONG_LONG_INT) *writebuf8i=(long long int*)malloc(sizeofmemory);
    if(
       ((datatype==MPI_UNSIGNED_CHAR)&&(*writebuf1 == NULL)) ||
       ((datatype==MPI_FLOAT)&&(*writebuf4 == NULL)) ||
       ((datatype==MPI_DOUBLE)&&(*writebuf8 == NULL)) ||
       ((datatype==MPI_LONG_DOUBLE)&&(*writebuf16 == NULL)) ||
       ((datatype==MPI_INT)&&(*writebuf4i == NULL)) ||
       ((datatype==MPI_LONG_LONG_INT)&&(*writebuf8i == NULL)) 
       ){
      dualfprintf(fail_file,"Can't initialize writebuf memory for mpiioromio\n");
      myexit(10010);
    }

    logsfprintf("mpiioromio_init end\n");

  }// end if operationtype==INITROMIO


  // ROMIO examples
  // http://www.sesp.cse.clrc.ac.uk/Publications/paraio/paraio/paraio.html
  // http://www.sesp.cse.clrc.ac.uk/Publications/paraio/paraio/node53.html

  // GODMARK: Could use Asynchronous IO if using MPI-2
  // GODMARK: Could use non-blocking MPI_File_iread() and MPI_File_iwrite() with MPIO_Wait() if wanted to write while continuing processes

  //
  // http://www.nersc.gov/nusers/resources/software/libs/io/mpiio.php
  // Collective: MPI_File_iread_all() ?
  // OR:
  // http://www-unix.mcs.anl.gov/mpi/www/www3/MPI_File_read_all_begin.html
  // int MPI_File_read_all_begin(MPI_File fh, void *buf, int count, MPI_Datatype datatype)
  ///
  // http://www-unix.mcs.anl.gov/mpi/www/www3/MPI_File_read_all_end.html
  // int MPI_File_read_all_end(MPI_File fh, void *buf, MPI_Status *status)
  //
  // better description:
  // http://mpi.deino.net/mpi_functions/MPI_File_read_all_begin.html


#if(USEMPI&&USEROMIO)
  if(operationtype==READROMIO){
    logsfprintf("mpiioromio_seperate begin\n");

    if(BARRIERROMIOPRE==1) MPI_Barrier(MPI_COMM_GRMHD); // force barrier before begin writing so avoids large number of unexpected buffer space required.

#if(MPIVERSION>=2)
    // non-blocking but need data from read immediately, so not taking advantage
    MPI_File_read_all_begin(fh, writebuf, bufcount, datatype);
    MPI_File_read_all_end(fh, writebuf, &status);
    MPI_File_close(&fh);
#else
    MPI_File_read_all(fh, writebuf, bufcount, datatype, &status);
    MPI_File_close(&fh);
#endif

    logsfprintf("mpiioromio_seperate end\n");
  }
  else if(operationtype==READFREEROMIO){
    logsfprintf("mpiioromio_seperate_free begin\n");
    free(writebuf);
    MPI_Type_free(&newtype);
    logsfprintf("mpiioromio_seperate_free end\n");
  }
  else  if(operationtype==WRITECLOSEROMIO){
    logsfprintf("mpiioromio_combine begin\n");

    if(BARRIERROMIOPRE==1) MPI_Barrier(MPI_COMM_GRMHD); // force barrier before begin writing so avoids large number of unexpected buffer space required.

#if(MPIVERSION>=2)
    MPI_File_write_all_begin(fh, writebuf, bufcount, datatype);
    // wait till later to end, close, free writebuf, and free newtype
#else
    // now write the buffer:
    MPI_File_write_all(fh, writebuf, bufcount, datatype, &status);
    MPI_File_close(&fh);
    // free buffer and type
    free(writebuf);
    MPI_Type_free(&newtype);
#endif

    logsfprintf("mpiioromio_combine end\n");
  
  }
  else if(operationtype==WRITEENDROMIO){
#if(MPIVERSION>=2)
    logsfprintf("mpiioromio_combine writeend begin\n");
    // now write the buffer:
    MPI_File_write_all_end(fh, writebuf, &status);
    //    logsfprintf("mpiioromio_combine 1\n");
    MPI_File_close(&fh);
    //    logsfprintf("mpiioromio_combine 2: writebuf=%d\n",writebuf);
    // free buffer and type
    free(writebuf);
    //    logsfprintf("mpiioromio_combine 3\n");
    MPI_Type_free(&newtype);
    logsfprintf("mpiioromio_combine writeend end\n");
#else
    // nothing to do

#endif
  }

#endif // end if USEMPI&&USEROMIO


}










/// size of buffer for MPI writing/reading
int set_sizeofmemory(int numbuff, int sizeofdatatype, int numcolumns, long long int *sizeofmemory)
{
  // can't trust that (int)ceil() along will be upper integer
  *sizeofmemory = (long long int)sizeofdatatype * (long long int)ROUND2LONGLONGINT(ceil((FTYPE)N1 * (FTYPE)N2 * (FTYPE)N3 * (FTYPE)NUMBUFFERS/(FTYPE)numcolumns))*(long long int)numcolumns * (long long int)numbuff ; // default

  return(0);
}

int set_maxnumsize(int numcolumns, long long int *maxnumsize)
{

  *maxnumsize=(long long int)(ROUND2LONGLONGINT(ceil(ROUND2LONGLONGINT(ceil((FTYPE)(N1*N2*N3*NUMBUFFERS)/(FTYPE)numcolumns))*(FTYPE)(numcolumns))));

  return(0);
}




int set_numbuffers(int numcolumns, int *numbuffers)
{
  int myval;

  myval=(int)ROUND2INT(ceil((FTYPE)numcolumns/((FTYPE)N1*(FTYPE)N2*(FTYPE)N3)));

  if(N1*N2*N3<numcolumns) *numbuffers=myval;
  else *numbuffers=1;

  return(0);
}

long long int gcountmod(int numcolumns)
{
  long long int myval;
  
  myval=(long long int)ROUND2LONGLONGINT(ceil((FTYPE)((FTYPE)N1*(FTYPE)N2*(FTYPE)N3*(FTYPE)NUMBUFFERS/(FTYPE)numcolumns)))*(long long int)numcolumns;

  return(myval);
}





#define DEBUGSINIT 0


/// initializes memory buffers for MPI combining for all MPI combine types(simple, minmem, romio)
/// mpi io sorted or not
void mpiios_init(int bintxt, int sorted, FILE ** fp, int which, int headerbytesize, char *filename, int numcolumns,
                 MPI_Datatype datatype, void **jonioptr, void **writebufptr)
{
  FTYPE fakevar;
  int fakei;

  int i;
  
  long long int sizeofmemory;
  int sizeofdatatype;

  long double **jonio16;
  double **jonio8;
  float **jonio4;
  unsigned char **jonio1;
  int **jonio4i;
  long long int **jonio8i;

  long double **writebuf16;
  double **writebuf8;
  float **writebuf4;
  unsigned char **writebuf1;
  int **writebuf4i;
  long long int **writebuf8i;

  int numfiles;
  char newfilename[MAXFILENAME];

  int trygettingfile;


#if(USEMPI)

  sizeofdatatype=getsizeofdatatype(datatype);

  if(datatype==MPI_UNSIGNED_CHAR){ jonio1=(unsigned char **)jonioptr; writebuf1=(unsigned char **)writebufptr; }
  else if(datatype==MPI_FLOAT){ jonio4=(float **)jonioptr; writebuf4=(float **)writebufptr; }
  else if(datatype==MPI_DOUBLE){ jonio8=(double**)jonioptr; writebuf8=(double**)writebufptr; }
  else if(datatype==MPI_LONG_DOUBLE){ jonio16=(long double**)jonioptr; writebuf16=(long double**)writebufptr; }
  else if(datatype==MPI_INT){ jonio4i=(int **)jonioptr; writebuf4i=(int **)writebufptr; }
  else if(datatype==MPI_LONG_LONG_INT){ jonio8i=(long long int **)jonioptr; writebuf8i=(long long int **)writebufptr; }
  
  logsfprintf("mpiios_init begin\n");


  /////////////////////
  //
  //  open files for non-ROMIO writing
  //
  /////////////////////
  
  if(myid==0){  // total on CPU=0, always, since mpicombine=1
    if(docolsplit){
      numfiles=numcolumns;
    }
    else numfiles=1;

    for(i=0;i<numfiles;i++){
      if(docolsplit&&(numfiles>1)){
        sprintf(newfilename,"%s-col%04d",filename,i); // must be same name as in dump_gen()
      }
      else{
        sprintf(newfilename,"%s",filename);
      }


#define NUMTRYGETTINGFILE 3
      for(trygettingfile=0;trygettingfile<NUMTRYGETTINGFILE;trygettingfile++){

        if (which == WRITEFILE){
          if(bintxt==BINARYOUTPUT)         fp[i] = fopen(newfilename, "a");
          else if(bintxt==TEXTOUTPUT)      fp[i] = fopen(newfilename, "at");
          else{
            dualfprintf(fail_file,"No such bintxt=%d\n",bintxt);
            myexit(62061);
          }
          // file pointer set correctly upon append
        }
        else if (which == READFILE){
          if(bintxt==BINARYOUTPUT)     fp[i] = fopen(newfilename, "rb");
          else if(bintxt==TEXTOUTPUT)  fp[i] = fopen(newfilename, "rt");
          else{
            dualfprintf(fail_file,"No such bintxt=%d\n",bintxt);
            myexit(62062);
          }
        }
        if (fp[i] == NULL) {
          dualfprintf(fail_file, "error opening file: %s\n", newfilename);
          if(trygettingfile==NUMTRYGETTINGFILE) myexit(62676);
          else{
            dualfprintf(fail_file,"Pausing for disk to get into synch: %d\n",trygettingfile);
            if(MPIAVOIDFORK){
              // compute a bit to fake a sleep
              fakevar=0.8;
#define FAKENUM 100 // number of square roots to take to pause
              for(fakei=i;fakei<=FAKENUM;fakei++){
                fakevar=sqrt(fakevar);
              }       
            }
            else{
              system("sleep 1"); // or use something to pause for 1 second
            }
          }
        }
        else{
          // then done with loop (ensure no inner loop and so set correct loop to end)
          trygettingfile=NUMTRYGETTINGFILE;
          break;
        }
      }// end if trygettingfile
 

      // must set file pointer to after header
      if (which == READFILE){
        fseek(fp[i],headerbytesize,SEEK_SET);
      }

    }// end over numfiles
  }// end if myid==0
  



  if ( (sorted==SORTED)&&(myid == 0) ){  // total on CPU=0 is stored in memory somehow for sorting before output

    /////////////
    //
    // determine the memory needed for the mpicombinetype for cpu=0
    //
    ///////////////

    truempicombinetype=mpicombinetype; // initial action
    if(truempicombinetype==MPICOMBINESIMPLE){
      sizeofmemory = (long long int)sizeofdatatype * (long long int)totalsize[1] * (long long int)totalsize[2] * (long long int)totalsize[3] * (long long int)numcolumns;
    }
    else if(truempicombinetype==MPICOMBINEMINMEM){
      // check this calculation against setuplinklist()'s cpulist0 array size!
      // 2 needed since need to read sequentially, then order it into the other buffer for writing
      set_sizeofmemory(2,sizeofdatatype, numcolumns, &sizeofmemory);
      if(sizeofmemory>(long long int)sizeofdatatype * (long long int)totalsize[1] * (long long int)totalsize[2] * (long long int)totalsize[3] * (long long int)numcolumns){
        sizeofmemory = (long long int)sizeofdatatype * (long long int)totalsize[1] * (long long int)totalsize[2] * (long long int)totalsize[3] * (long long int)numcolumns;
        truempicombinetype=MPICOMBINESIMPLE; // then switch to simple method
      }
      // need memory to be at least larger than number of columns (X2 for 2 buffers)
      // don't want to work with chunks smaller than # of columns, and all chunks should come in # of column chunks times some integer multiple
      if(sizeofmemory<(long long int)sizeofdatatype*(long long int)numcolumns*(long long int)2) sizeofmemory=(long long int)sizeofdatatype*(long long int)numcolumns*(long long int)2;
      if(sizeofmemory<(long long int)2*(long long int)numcolumns){
        dualfprintf(fail_file,"problem, sizeofmemory=%lld < %lld=2*numcolumns\n",sizeofmemory,(long long int)2*(long long int)numcolumns);
        myexit(102000);
      }

    }
    joniosize=sizeofmemory/((long long int)sizeofdatatype);

#if(DEBUGMINMEM)
    dualfprintf(fail_file,"jonio sizeofmemory=%lld sizeofdatatype=%d\n",sizeofmemory,sizeofdatatype);
#endif

    if(datatype==MPI_UNSIGNED_CHAR) *jonio1=(unsigned char*)malloc(sizeofmemory);
    else if(datatype==MPI_FLOAT) *jonio4=(float*)malloc(sizeofmemory);
    else if(datatype==MPI_DOUBLE) *jonio8 =(double*)malloc(sizeofmemory);
    else if(datatype==MPI_LONG_DOUBLE) *jonio16 =(long double*)malloc(sizeofmemory);
    else if(datatype==MPI_INT) *jonio4i=(int*)malloc(sizeofmemory);
    else if(datatype==MPI_LONG_LONG_INT) *jonio8i=(long long int*)malloc(sizeofmemory);
    if(
       (datatype==MPI_UNSIGNED_CHAR)&&(*jonio1 == NULL) ||
       (datatype==MPI_FLOAT)&&(*jonio4 == NULL) ||
       (datatype==MPI_DOUBLE)&&(*jonio8 == NULL) ||
       (datatype==MPI_LONG_DOUBLE)&&(*jonio16 == NULL) ||
       (datatype==MPI_INT)&&(*jonio4i == NULL) ||
       (datatype==MPI_LONG_LONG_INT)&&(*jonio8i == NULL) 
       ){
      dualfprintf(fail_file, "Can't initialize jonio memory for mpiios_init\n");
      myexit(16010);
    }
#if(DEBUGSINIT)
    /*
      for(i=0;i<sizeofmemory/datatype;i++){
      if(datatype==sizeof(long double)) (*jonio16)[i]=-10000.0000;
      if(datatype==sizeof(double)) (*jonio8)[i]=-10000.0000;
      if(datatype==sizeof(float)) (*jonio4)[i]=-10000.0000;
      if(datatype==sizeof(unsigned char)) (*jonio1)[i]=100;
      }
      dualfprintf(fail_file,"got here1\n");
    */
#endif

    logsfprintf("mpiios_init jonio: sizeofmemory=%lld sizeofdatatype=%d\n",sizeofmemory,sizeofdatatype);

  }


  //////////////////////
  //
  // need to tell all CPUS the status of changed global vars
  //
  ///////////////////////

  MPI_Bcast(&truempicombinetype,1,MPI_INT,MPIid[0], MPI_COMM_GRMHD);


  ///////////////////////////////
  //
  // open per CPU writebuf
  //
  /////////////////////////////////

  if(truempicombinetype==MPICOMBINESIMPLE){
    sizeofmemory = (long long int)sizeofdatatype * (long long int)N1 * (long long int)N2 * (long long int)N3 * (long long int)numcolumns ;
  }
  else if(truempicombinetype==MPICOMBINEMINMEM){
    // maximum cpu=0 could require under any case
    set_sizeofmemory(1,sizeofdatatype, numcolumns, &sizeofmemory);
    if(sizeofmemory>(long long int)sizeofdatatype*(long long int)N1*(long long int)N2*(long long int)N3*(long long int)numcolumns) sizeofmemory=(long long int)sizeofdatatype*(long long int)N1*(long long int)N2*(long long int)N3*(long long int)numcolumns; // can't ask for more!
    if(sizeofmemory<(long long int)sizeofdatatype*(long long int)numcolumns) sizeofmemory=(long long int)sizeofdatatype*(long long int)numcolumns; // minimum, chunk at minimum of number of columns
    if(sizeofmemory<(long long int)numcolumns){
      dualfprintf(fail_file,"problem, sizeofmemory=%lld < %lld=numcolumns\n",sizeofmemory,(long long int)numcolumns);
      myexit(6900);
    }
  }

  writebufsize=sizeofmemory/((long long int)sizeofdatatype); // never used currently

#if(DEBUGMINMEM)
  dualfprintf(fail_file,"writebuf sizeofmemory=%lld sizeofdatatype=%d\n",sizeofmemory,sizeofdatatype);
#endif

  if(datatype==MPI_UNSIGNED_CHAR) *writebuf1=(unsigned char*)malloc(sizeofmemory);
  else if(datatype==MPI_FLOAT) *writebuf4=(float*)malloc(sizeofmemory);
  else if(datatype==MPI_DOUBLE) *writebuf8 =(double*)malloc(sizeofmemory);
  else if(datatype==MPI_LONG_DOUBLE) *writebuf16 =(long double*)malloc(sizeofmemory);
  else if(datatype==MPI_INT) *writebuf4i=(int*)malloc(sizeofmemory);
  else if(datatype==MPI_LONG_LONG_INT) *writebuf8i=(long long int*)malloc(sizeofmemory);
  if(
     (datatype==MPI_UNSIGNED_CHAR)&&(*writebuf1 == NULL) ||
     (datatype==MPI_FLOAT)&&(*writebuf4 == NULL) ||
     (datatype==MPI_DOUBLE)&&(*writebuf8 == NULL) ||
     (datatype==MPI_LONG_DOUBLE)&&(*writebuf16 == NULL) ||
     (datatype==MPI_INT)&&(*writebuf4i == NULL) ||
     (datatype==MPI_LONG_LONG_INT)&&(*writebuf8i == NULL) 
     ){
    dualfprintf(fail_file, "Can't initialize writebuf memory for mpiios_init: datatype=%d sizeofmemory=%d\n",datatype,sizeofmemory);
    myexit(86726);
  }
#if(DEBUGSINIT)
  /*
    for(i=0;i<sizeofmemory/sizeofdatatype;i++){
    if(datatype==sizeof(long double)) (*writebuf16)[i]=-10000.0000;
    if(datatype==sizeof(double)) (*writebuf8)[i]=-10000.0000;
    if(datatype==sizeof(float)) (*writebuf4)[i]=-10000.0000;
    if(datatype==sizeof(unsigned char)) (*writebuf1)[i]=100;
    }
    dualfprintf(fail_file,"got here2\n");
  */
#endif
  logsfprintf("mpiios_init end: sizeofmemory=%lld sizeofdatatype=%d\n",sizeofmemory,sizeofdatatype);
#endif

}



void mpiiomin_final(int numcolumns,FILE **fp, void *jonio, void *writebuf)
{
  int i;
  int numfiles;

  free(writebuf);  // writebuf used by each CPU
  if(myid==0){
    free(jonio);  // used by CPU=0
    if(docolsplit){
      numfiles=numcolumns;
    }
    else numfiles=1;
    for(i=0;i<numfiles;i++){
      fclose(fp[i]);
      fp[i] = NULL;
    }
  }
}

#if(USEMPI)


void mpiios_combine(int bintxt, MPI_Datatype datatype, int numcolumns,
                    FILE ** fp, void *jonio, void *writebuf)
{
  long long int i, j, k, l, col, mapvaluejonio, mapvaluetempbuf;
#if(USEMPI)
  MPI_Request rrequest;
  MPI_Request srequest;
#endif
  int othercpupos[COMPDIM + 1];

  unsigned char *jonio1;
  float *jonio4;
  double *jonio8;
  long double *jonio16;
  int *jonio4i;
  long long int *jonio8i;

  unsigned char *writebuf1;
  float *writebuf4;
  double *writebuf8;
  long double *writebuf16;
  int *writebuf4i;
  long long int *writebuf8i;
  
  int numfiles;
  int sizeofdatatype;
  
  logsfprintf("mpiios begin combine\n");

  sizeofdatatype=getsizeofdatatype(datatype);

  if (datatype == MPI_UNSIGNED_CHAR) jonio1 = (unsigned char *) jonio;
  else if (datatype == MPI_FLOAT) jonio4 = (float *) jonio;
  else if (datatype == MPI_DOUBLE) jonio8 = (double *) jonio;
  else if (datatype == MPI_LONG_DOUBLE) jonio16 = (long double *) jonio;
  else if (datatype == MPI_INT) jonio4i = (int *) jonio;
  else if (datatype == MPI_LONG_LONG_INT) jonio8i = (long long int *) jonio;
  
  if (datatype == MPI_UNSIGNED_CHAR) writebuf1 = (unsigned char *) writebuf;
  else if (datatype == MPI_FLOAT) writebuf4 = (float *) writebuf;
  else if (datatype == MPI_DOUBLE) writebuf8 = (double *) writebuf;
  else if (datatype == MPI_LONG_DOUBLE) writebuf16 = (long double *) writebuf;
  else if (datatype == MPI_INT) writebuf4i = (int *) writebuf;
  else if (datatype == MPI_LONG_LONG_INT) writebuf8i = (long long int *) writebuf;


#if(USEMPI)
  // no need for tempbuf, works since first write to jonio is CPU=0's writebuf
  if(myid!=0) MPI_Isend(writebuf, N1 * N2 * N3 * numcolumns, datatype, MPIid[0], myid, MPI_COMM_GRMHD, &srequest);
  if (myid == 0) {
    for (l = 0; l < numprocs; l++) {
      logsfprintf("on myid==0: mpiios combine: %d of numprocs=%d data=%lld\n",l,numprocs,(long long int) (N1 * N2 * N3 * numcolumns) );
      if(l!=0){
        MPI_Irecv(writebuf, N1 * N2 * N3 * numcolumns, datatype, MPIid[l], l, MPI_COMM_GRMHD, &rrequest);
        MPI_Wait(&rrequest, &mpichstatus);
      }

      othercpupos[1]=l%ncpux1;
      othercpupos[2]=(int)((l%(ncpux1*ncpux2))/ncpux1);
      othercpupos[3]=(int)(l/(ncpux1*ncpux2));

      // now fill jonio with proper sequence (i.e. tiled mapping)
      for (k = 0; k < N3; k++) for (j = 0; j < N2; j++)
                                 for (i = 0; i < N1; i++)
                                   for (col = 0; col < numcolumns; col++) {
                                     // mapvaluejonio is global single-dimensional index for position in total CPU space - in C order for "array[k][j][i]"
                                     // Note that this storage mapping function forces disk to have i as fastest regardless of order stored in original array.
                                     // This is fully compatible with any ORDERSTORAGE choice since global arrays are still properly accessed and this mapping just forces writebuf to be filled or read with fixed mapping function as desirable, so that user can change ORDERSTORAGE without the files written being changed.
                                     mapvaluejonio =
                                       + (long long int)ncpux1 * (long long int)N1 * (long long int)ncpux2 * (long long int)N2 * (long long int)numcolumns * ((long long int)k + (long long int)othercpupos[3] * (long long int)N3)
                                       + (long long int)ncpux1 * (long long int)N1 * (long long int)numcolumns * ((long long int)j + (long long int)othercpupos[2] * (long long int)N2)
                                       + (long long int)numcolumns * ((long long int)i + (long long int)othercpupos[1] * (long long int)N1)
                                       + (long long int)col;

                                     // mapvaluetempbuf is a single-buffer single-dimensional index for the position in the buffer in C-order
                                     mapvaluetempbuf =
                                       + (long long int)k * (long long int)N1 * (long long int)N2 * (long long int)numcolumns
                                       + (long long int)j * (long long int)N1 * (long long int)numcolumns
                                       + (long long int)i * (long long int)numcolumns + (long long int)col;
     
                                     if (datatype == MPI_UNSIGNED_CHAR) jonio1[mapvaluejonio] = writebuf1[mapvaluetempbuf];
                                     else if (datatype == MPI_FLOAT) jonio4[mapvaluejonio] = writebuf4[mapvaluetempbuf];
                                     else if (datatype == MPI_DOUBLE) jonio8[mapvaluejonio] = writebuf8[mapvaluetempbuf];
                                     else if (datatype == MPI_LONG_DOUBLE) jonio16[mapvaluejonio] = writebuf16[mapvaluetempbuf];
                                     else if (datatype == MPI_INT) jonio4i[mapvaluejonio] = writebuf4i[mapvaluetempbuf];
                                     else if (datatype == MPI_LONG_LONG_INT) jonio8i[mapvaluejonio] = writebuf8i[mapvaluetempbuf];
                                   }
    }
  }
  if(myid!=0) MPI_Wait(&srequest, &mpichstatus);
  free(writebuf);  // writebuf used by each CPU

  if (myid == 0) {
    logsfprintf("on myid==0: mpiios combine: writing: %lld\n",(long long int)totalsize[1] * (long long int)totalsize[2] * (long long int)totalsize[3] * (long long int)numcolumns);
    // now write out collected data using CPU=0
    if(docolsplit){
      numfiles=numcolumns;
    }
    else numfiles=1;

    if(bintxt==BINARYOUTPUT){
      // GODMARK: Below may be too large of a size for fwrite() to output
      if(numfiles==1) fwrite(jonio, sizeofdatatype,totalsize[1] * totalsize[2] * totalsize[3] * numcolumns, *fp);
      else{
        for(i=0;i<totalsize[1]*totalsize[2]*totalsize[3]*numcolumns;i++){
          if (datatype == MPI_UNSIGNED_CHAR) fwrite(&jonio1[i], sizeofdatatype,1, fp[i%numfiles]);
          else if (datatype == MPI_FLOAT) fwrite(&jonio4[i], sizeofdatatype,1, fp[i%numfiles]);
          else if (datatype == MPI_DOUBLE) fwrite(&jonio8[i], sizeofdatatype,1, fp[i%numfiles]);
          else if (datatype == MPI_LONG_DOUBLE) fwrite(&jonio16[i], sizeofdatatype,1, fp[i%numfiles]);
          else if (datatype == MPI_INT) fwrite(&jonio4i[i], sizeofdatatype,1, fp[i%numfiles]);
          else if (datatype == MPI_LONG_LONG_INT) fwrite(&jonio8i[i], sizeofdatatype,1, fp[i%numfiles]);
        }
      }
    }
    else if(bintxt==TEXTOUTPUT){ // properly ordered, so just dump it
      for(i=0;i<totalsize[1]*totalsize[2]*totalsize[3]*numcolumns;i++){
        if (datatype == MPI_UNSIGNED_CHAR) fprintf(fp[i%numfiles],"%c",jonio1[i]);
        else if (datatype == MPI_FLOAT) fprintf(fp[i%numfiles],"%15.7g",jonio4[i]);
        else if (datatype == MPI_DOUBLE) fprintf(fp[i%numfiles],"%21.15g",jonio8[i]);
        else if (datatype == MPI_LONG_DOUBLE) fprintf(fp[i%numfiles],"%31.25Lg",jonio16[i]);
        else if (datatype == MPI_INT) fprintf(fp[i%numfiles],"%d",jonio4i[i]);
        else if (datatype == MPI_LONG_LONG_INT) fprintf(fp[i%numfiles],"%lld",jonio8i[i]);
        if(numfiles==1){
          if((i+1)%numcolumns) fprintf(*fp," ");
          else fprintf(*fp,"\n");
        }
        else  fprintf(fp[i%numfiles],"\n");
      }
    }
    logsfprintf("on myid==0: mpiios combine: freeing\n");
    free(jonio);  // used by CPU=0
    for(i=0;i<numfiles;i++){
      fclose(fp[i]);
      fp[i] = NULL;
    }
  }// end if myid==0
#endif

  logsfprintf("mpiios end combine\n");

}



void mpiios_seperate(int bintxt, int stage, MPI_Datatype datatype, int numcolumns,
                     FILE ** fp, void *jonio,
                     void *writebuf)
{
  long long int i, j, k, l, col, mapvaluejonio, mapvaluetempbuf;
#if(USEMPI)
  MPI_Request rrequest;
  MPI_Request srequest;
#endif
  int othercpupos[COMPDIM + 1];

  int sizeofdatatype;

  unsigned char *jonio1;
  float *jonio4;
  double *jonio8;
  long double *jonio16;
  int *jonio4i;
  long long int *jonio8i;

  unsigned char *writebuf1;
  float *writebuf4;
  double *writebuf8;
  long double *writebuf16;
  int *writebuf4i;
  long long int *writebuf8i;
  unsigned short short4char;

  int numfiles;


  logsfprintf("mpiios begin seperate\n");

  sizeofdatatype=getsizeofdatatype(datatype);

  if (datatype == MPI_UNSIGNED_CHAR) jonio1 = (unsigned char *) jonio;
  else if (datatype == MPI_FLOAT) jonio4 = (float *) jonio;
  else if (datatype == MPI_DOUBLE) jonio8 = (double *) jonio;
  else if (datatype == MPI_LONG_DOUBLE) jonio16 = (long double *) jonio;
  else if (datatype == MPI_INT) jonio4i = (int *) jonio;
  else if (datatype == MPI_LONG_LONG_INT) jonio8i = (long long int *) jonio;


  if (datatype == MPI_UNSIGNED_CHAR) writebuf1 = (unsigned char *) writebuf;
  else if (datatype == MPI_FLOAT) writebuf4 = (float *) writebuf;
  else if (datatype == MPI_DOUBLE) writebuf8 = (double *) writebuf;
  else if (datatype == MPI_LONG_DOUBLE) writebuf16 = (long double *) writebuf;
  else if (datatype == MPI_INT) writebuf4i = (int *) writebuf;
  else if (datatype == MPI_LONG_LONG_INT) writebuf8i = (long long int *) writebuf;


#if(USEMPI)

  if (stage == STAGE1) {
    if (myid == 0) {
      if(docolsplit){
        numfiles=numcolumns;
      }
      else numfiles=1;
      
      if(bintxt==BINARYOUTPUT){
        // first let cpu=0 read data
        // GODMARK: Below may be too big for fread() to read
        if(numfiles==1) fread(jonio, sizeofdatatype,totalsize[1] * totalsize[2] * totalsize[3] * numcolumns, *fp);
        else{
          for(i=0;i<totalsize[1]*totalsize[2]*totalsize[3]*numcolumns;i++){
            if (datatype == MPI_UNSIGNED_CHAR) fread(&jonio1[i], sizeofdatatype,1, fp[i%numfiles]);
            else if (datatype == MPI_FLOAT) fread(&jonio4[i], sizeofdatatype,1, fp[i%numfiles]);
            else if (datatype == MPI_DOUBLE)  fread(&jonio8[i], sizeofdatatype,1, fp[i%numfiles]);
            else if (datatype == MPI_LONG_DOUBLE) fread(&jonio16[i], sizeofdatatype,1, fp[i%numfiles]);
            else if (datatype == MPI_INT) fread(&jonio4i[i], sizeofdatatype,1, fp[i%numfiles]);
            else if (datatype == MPI_LONG_LONG_INT) fread(&jonio8i[i], sizeofdatatype,1, fp[i%numfiles]);
          }
        }
      }
      else if(bintxt==TEXTOUTPUT){ // properly ordered, so just dump it
        for(i=0;i<totalsize[1]*totalsize[2]*totalsize[3]*numcolumns;i++){
          if (datatype == MPI_UNSIGNED_CHAR){
            fscanf(fp[i%numfiles],"%hu",&short4char);
            jonio1[i]=short4char; // convert
          }
          else if (datatype == MPI_FLOAT) fscanf(fp[i%numfiles],"%f",&jonio4[i]);
          else if (datatype == MPI_DOUBLE) fscanf(fp[i%numfiles],"%lf",&jonio8[i]);
          else if (datatype == MPI_LONG_DOUBLE) fscanf(fp[i%numfiles],"%Lf",&jonio16[i]);
          else if (datatype == MPI_INT) fscanf(fp[i%numfiles],"%d",&jonio4i[i]);
          else if (datatype == MPI_LONG_LONG_INT) fscanf(fp[i%numfiles],"%lld",&jonio8i[i]);
        }
      }
    }
    // writebuf is CPU=0's tempbuf for each CPU, including CPU=0, which is done last
    if (myid == 0) {
      for (l = numprocs-1 ; l >=0; l--) {

        othercpupos[1]=l%ncpux1;
        othercpupos[2]=(int)((l%(ncpux1*ncpux2))/ncpux1);
        othercpupos[3]=(int)(l/(ncpux1*ncpux2));

        // now unfill jonio with proper sequence (i.e. tiled mapping)
        for (k = 0; k < N3; k++) for (j = 0; j < N2; j++) {
            for (i = 0; i < N1; i++) {
              for (col = 0; col < numcolumns; col++) {

                // mapvaluejonio is global single-dimensional index for position in total CPU space - in C order for "array[k][j][i]"
                mapvaluejonio =
                  + (long long int)ncpux1 * (long long int)N1 * (long long int)ncpux2 * (long long int)N2 * (long long int)numcolumns * ((long long int)k + (long long int)othercpupos[3] * (long long int)N3)
                  + (long long int)ncpux1 * (long long int)N1 * (long long int)numcolumns * ((long long int)j + (long long int)othercpupos[2] * (long long int)N2)
                  + (long long int)numcolumns * ((long long int)i + (long long int)othercpupos[1] * (long long int)N1)
                  + (long long int)col;

                // mapvaluetempbuf is a single-buffer single-dimensional index for the position in the buffer in C-order
                mapvaluetempbuf =
                  + (long long int)k * (long long int)N1 * (long long int)N2 * (long long int)numcolumns
                  + (long long int)j * (long long int)N1 * (long long int)numcolumns
                  + (long long int)i * (long long int)numcolumns + (long long int)col;
     
                // debug check
                if((mapvaluetempbuf<0)||(mapvaluetempbuf>=(long long int)N1*N2*N3*numcolumns)){
                  dualfprintf(fail_file,"mapvaluetempbuf out of range: %d\n",mapvaluetempbuf);
                  myexit(96726);
                }
                if((mapvaluejonio<0)||(mapvaluejonio>=(long long int)totalsize[1]*(long long int)totalsize[2]*(long long int)totalsize[3]*(long long int)numcolumns)){
                  dualfprintf(fail_file,"mapvaluejonio out of range: %d\n",mapvaluejonio);
                  myexit(6726);
                }
                if (datatype == MPI_UNSIGNED_CHAR) writebuf1[mapvaluetempbuf] = jonio1[mapvaluejonio];
                if (datatype == MPI_FLOAT) writebuf4[mapvaluetempbuf] = jonio4[mapvaluejonio];
                if (datatype == MPI_DOUBLE) writebuf8[mapvaluetempbuf] = jonio8[mapvaluejonio];
                if (datatype == MPI_LONG_DOUBLE) writebuf16[mapvaluetempbuf] = jonio16[mapvaluejonio];
                if (datatype == MPI_INT) writebuf4i[mapvaluetempbuf] = jonio4i[mapvaluejonio];
                if (datatype == MPI_LONG_LONG_INT) writebuf8i[mapvaluetempbuf] = jonio8i[mapvaluejonio];
              }
            }
          }
        if(l!=0){
   
          MPI_Isend(writebuf, N1 * N2 * N3 * numcolumns, datatype, MPIid[l], l, MPI_COMM_GRMHD, &srequest);
          MPI_Wait(&srequest, &mpichstatus);
        }
      }
      free(jonio); // done with jonio after loop
    }
    else{
      // chosen CPU to receive data from CPU=0
      
      MPI_Irecv(writebuf, N1 * N2 * N3 * numcolumns, datatype, MPIid[0], myid, MPI_COMM_GRMHD, &rrequest);
      MPI_Wait(&rrequest, &mpichstatus); // writebuf used until      
    }
  } else if (stage == STAGE2) {
    free(writebuf);
    if (myid == 0) {
      fclose(*fp);
      *fp = NULL;
    }
  }
#endif

  logsfprintf("mpiios end seperate\n");



}


void mpiiotu_combine(MPI_Datatype datatype, int numcolumns, FILE ** fp, void *writebuf)
{
  long long int i, j, k, l, col;
#if(USEMPI)
  MPI_Request rrequest;
  MPI_Request srequest;
#endif
  int othercpupos[COMPDIM + 1];

  int sizeofdatatype;

  unsigned char *writebuf1;
  float *writebuf4;
  double *writebuf8;
  long double *writebuf16;
  int *writebuf4i;
  long long int *writebuf8i;
  
  long long int bufferwritemap;
  int numfiles;

  sizeofdatatype=getsizeofdatatype(datatype);

  if (datatype == MPI_UNSIGNED_CHAR) writebuf1 = (unsigned char *) writebuf;
  else if (datatype == MPI_FLOAT) writebuf4 = (float *) writebuf;
  else if (datatype == MPI_DOUBLE) writebuf8 = (double *) writebuf;
  else if (datatype == MPI_LONG_DOUBLE) writebuf16 = (long double *) writebuf;
  else if (datatype == MPI_INT) writebuf4i = (int *) writebuf;
  else if (datatype == MPI_LONG_LONG_INT) writebuf8i = (long long int *) writebuf;
  

#if(USEMPI)
  if(myid!=0) MPI_Isend(writebuf, N1 * N2 * N3 * numcolumns, datatype, MPIid[0], myid, MPI_COMM_GRMHD, &srequest);
  if (myid == 0) {   
    if(docolsplit){
      numfiles=numcolumns;
    }
    else numfiles=1;

    // done in forward order, no need to use tempbuf since CPU=0's writebuf is first out
    for (l = 0; l <numprocs; l++) {
      if(l!=0){
        MPI_Irecv(writebuf, N1 * N2 * N3 * numcolumns, datatype, MPIid[l], l, MPI_COMM_GRMHD, &rrequest);
        MPI_Wait(&rrequest, &mpichstatus);
      }
      // now write writebuf
      DUMPGENLOOP{
        for(col=0;col<numcolumns;col++){

          bufferwritemap=(long long int)col + (long long int)numcolumns*((long long int)i+(long long int)N1*(long long int)j+(long long int)N1*(long long int)N2*(long long int)k);

          if(datatype== MPI_UNSIGNED_CHAR){
            fprintf(fp[(bufferwritemap)%numfiles],"%c ",writebuf1[bufferwritemap]);
          }
          else if(datatype==MPI_FLOAT){
            fprintf(fp[(bufferwritemap)%numfiles],"%15.7g ",writebuf4[bufferwritemap]);
          }
          else if(datatype==MPI_DOUBLE){
            fprintf(fp[(bufferwritemap)%numfiles],"%21.15g ",writebuf8[bufferwritemap]);
          }
          else if(datatype==MPI_LONG_DOUBLE){
            fprintf(fp[(bufferwritemap)%numfiles],"%31.25Lg ",writebuf16[bufferwritemap]);
          }
          else if(datatype==MPI_INT){
            fprintf(fp[(bufferwritemap)%numfiles],"%d ",writebuf4i[bufferwritemap]);
          }
          else if(datatype==MPI_LONG_LONG_INT){
            fprintf(fp[(bufferwritemap)%numfiles],"%lld ",writebuf8i[bufferwritemap]);
          }
          if(numfiles>1) fprintf(fp[(bufferwritemap)%numfiles],"\n");
        }
        if(numfiles==1) fprintf(*fp,"\n");
      }
    }    
  }
  if(myid!=0) MPI_Wait(&srequest, &mpichstatus);
  free(writebuf);  // writebuf used by each CPU

  if (myid == 0) {
    for(i=0;i<numfiles;i++){
      fclose(fp[i]);
      fp[i] = NULL;
    }
  }
#endif

}

/// fill writebuf with each cpu's data set,using CPU=0 to process the file
void mpiiotu_seperate(int stage, MPI_Datatype datatype, int numcolumns,
                      FILE ** fp,void *writebuf)
{
  long long int i, j, k, l, col;
#if(USEMPI)
  MPI_Request rrequest;
  MPI_Request srequest;
#endif
  int othercpupos[COMPDIM + 1];

  void *tempbuf;
  void *sendbuf;

  int numfiles;

  int sizeofdatatype;

  unsigned char *writebuf1;
  float *writebuf4;
  double *writebuf8;
  long double *writebuf16;
  int *writebuf4i;
  long long int *writebuf8i;

  unsigned char *tempbuf1;
  float *tempbuf4;
  double *tempbuf8;
  long double *tempbuf16;
  int *tempbuf4i;
  long long int *tempbuf8i;

  unsigned char *sendbuf1;
  float *sendbuf4;
  double *sendbuf8;
  long double *sendbuf16;
  int *sendbuf4i;
  long long int *sendbuf8i;

  unsigned short short4char;

  long long int bufferwritemap;



  sizeofdatatype=getsizeofdatatype(datatype);

  if (datatype == MPI_UNSIGNED_CHAR) writebuf1 = (unsigned char *) writebuf;
  else if (datatype == MPI_FLOAT) writebuf4 = (float *) writebuf;
  else if (datatype == MPI_DOUBLE) writebuf8 = (double *) writebuf;
  else if (datatype == MPI_LONG_DOUBLE) writebuf16 = (long double *) writebuf;
  else if (datatype == MPI_INT) writebuf4i = (int *) writebuf;
  else if (datatype == MPI_LONG_LONG_INT) writebuf8i = (long long int *) writebuf;

  if(myid==0){
    if((tempbuf=malloc(sizeofdatatype*N1*N2*N3*numcolumns))==NULL){
      dualfprintf(fail_file,"Can't open tempbuf in gammieio_sep\n");
      myexit(7566);
    }

    if (datatype == MPI_UNSIGNED_CHAR) tempbuf1 = (unsigned char *) tempbuf;
    else if (datatype == MPI_FLOAT) tempbuf4 = (float *) tempbuf;
    else if (datatype == MPI_DOUBLE) tempbuf8 = (double *) tempbuf;
    else if (datatype == MPI_LONG_DOUBLE) tempbuf16 = (long double *) tempbuf;
    else if (datatype == MPI_INT) tempbuf4i = (int *) tempbuf;
    else if (datatype == MPI_LONG_LONG_INT) tempbuf8i = (long long int *) tempbuf;

  }

#if(USEMPI)

  if (stage == 1) {
    if (myid == 0) {
      if(docolsplit){
        numfiles=numcolumns;
      }
      else numfiles=1;

      for (l = 0; l < numprocs; l++) {
        if(l==0){
          sendbuf=writebuf;
          sendbuf1=writebuf1;
          sendbuf4=writebuf4;
          sendbuf8=writebuf8;
          sendbuf4i=writebuf4i;
          sendbuf8i=writebuf8i;
        }
        else{
          sendbuf=tempbuf;
          sendbuf1=tempbuf1;
          sendbuf4=tempbuf4;
          sendbuf8=tempbuf8;
          sendbuf4i=tempbuf4i;
          sendbuf8i=tempbuf8i;
        }

        DUMPGENLOOP for (col = 0; col < numcolumns; col++) {

          bufferwritemap=(long long int)col+(long long int)numcolumns*((long long int)i+(long long int)N1*(long long int)j+(long long int)N1*(long long int)N2*(long long int)k);

          if(datatype==MPI_UNSIGNED_CHAR){
            fscanf(fp[(bufferwritemap)%numfiles],"%hu",&short4char);
            sendbuf1[bufferwritemap]=short4char;
          }
          else if(datatype==MPI_FLOAT) fscanf(fp[(bufferwritemap)%numfiles],"%f",&sendbuf4[bufferwritemap]);
          else if(datatype==MPI_DOUBLE) fscanf(fp[(bufferwritemap)%numfiles],"%lf",&sendbuf8[bufferwritemap]);
          else if(datatype==MPI_LONG_DOUBLE) fscanf(fp[(bufferwritemap)%numfiles],"%Lf",&sendbuf16[bufferwritemap]);
          else if(datatype==MPI_INT) fscanf(fp[(bufferwritemap)%numfiles],"%d",&sendbuf4i[bufferwritemap]);
          else if(datatype==MPI_LONG_LONG_INT) fscanf(fp[(bufferwritemap)%numfiles],"%lld",&sendbuf8i[bufferwritemap]);
        }
        if(l!=0){
          MPI_Isend(sendbuf, N1 * N2 * N3 * numcolumns, datatype, MPIid[l], l,  MPI_COMM_GRMHD, &srequest);
          // have to wait before filling sendbuf buffer again for next CPU
          MPI_Wait(&srequest, &mpichstatus);
        }
      }
      free(tempbuf);
    }
    else{
      MPI_Irecv(writebuf, N1 * N2 * N3 * numcolumns, datatype, MPIid[0], myid, MPI_COMM_GRMHD, &rrequest);
      MPI_Wait(&rrequest, &mpichstatus); // writebuf used until
    }
  } else if (stage == 2) {
    free(writebuf);
    if (myid == 0) {
      for(i=0;i<numfiles;i++){
        fclose(fp[i]);
        fp[i] = NULL;
      }
    }
  }
#endif



}




#endif


int getsizeofdatatype(MPI_Datatype datatype)
{
  int sizeofdatatype;

  if (datatype == MPI_UNSIGNED_CHAR) sizeofdatatype=sizeof(unsigned char);
  else if (datatype == MPI_FLOAT) sizeofdatatype=sizeof(float);
  else if (datatype == MPI_DOUBLE) sizeofdatatype=sizeof(double);
  else if (datatype == MPI_LONG_DOUBLE) sizeofdatatype=sizeof(long double);
  else if (datatype == MPI_INT) sizeofdatatype=sizeof(int);
  else if (datatype == MPI_LONG) sizeofdatatype=sizeof(long int); // sometimes different than (int)
  else if (datatype == MPI_LONG_LONG_INT) sizeofdatatype=sizeof(long long int);
  else{
    dualfprintf(fail_file,"No such datatype=%d\n",datatype);
    myexit(8676);
  }

  return(sizeofdatatype);
}


/// a simple max, assumes local cpu max already found
/// sends results back to all cpus
void mpimax(SFTYPE*maxptr)
{
  SFTYPE send;
  
#if(USEMPI)
  send = *maxptr;
  MPI_Allreduce(&send, maxptr, 1, MPI_SFTYPE, MPI_MAX,MPI_COMM_GRMHD);
#endif
}


/// a simple max, assumes local cpu max already found
/// sends results back to all cpus
void mpiimax(int*maxptr)
{
  int send;
  
#if(USEMPI)
  send = *maxptr;
  MPI_Allreduce(&send, maxptr, 1, MPI_INT, MPI_MAX,MPI_COMM_GRMHD);
#endif
}

/// a simple max, assumes local cpu max already found
/// sends results back to all cpus
void mpiisum(int*sumptr)
{
  int send;
  
#if(USEMPI)
  send = *sumptr;
  MPI_Allreduce(&send, sumptr, 1, MPI_INT, MPI_SUM,MPI_COMM_GRMHD);
#endif
}

/// a simple max, assumes local cpu max already found
/// sends results back to all cpus
void mpiFTYPEsum(FTYPE*sumptr)
{
  FTYPE send;
  
#if(USEMPI)
  send = *sumptr;
  MPI_Allreduce(&send, sumptr, 1, MPI_FTYPE, MPI_SUM,MPI_COMM_GRMHD);
#endif
}


/// as above, but only myid=recvid gets result
void mpiisum0(int*sumptr, int recvid)
{
  int send;
  
#if(USEMPI)
  send = *sumptr;
  MPI_Reduce(&send, sumptr, 1, MPI_INT, MPI_SUM,MPIid[recvid], MPI_COMM_GRMHD);
#endif
}

/// as above, but only myid=recvid gets result
void mpildsum0(long int*sumptr, int recvid)
{
  long int send;
  
#if(USEMPI)
  send = *sumptr;
  MPI_Reduce(&send, sumptr, 1, MPI_LONG, MPI_SUM,MPIid[recvid], MPI_COMM_GRMHD);
#endif
}

/// a simple max, assumes local cpu max already found
/// sends results back to all cpus
void mpimin(SFTYPE*minptr)
{
  SFTYPE send;
  
#if(USEMPI)
  send = *minptr;
  MPI_Allreduce(&send, minptr, 1, MPI_SFTYPE, MPI_MIN,MPI_COMM_GRMHD);
#endif
}

/// a simple max, assumes local cpu max already found
/// sends results back to all cpus
void mpifmin(FTYPE*minptr)
{
  FTYPE send;
  
#if(USEMPI)
  send = *minptr;
  MPI_Allreduce(&send, minptr, 1, MPI_FTYPE, MPI_MIN,MPI_COMM_GRMHD);
#endif
}


/// only used for images to get min max over whole dumped domain.
void prminmaxsum(FTYPE (*p)[NSTORE2][NSTORE3][NPR], int start,int nmemb, FTYPE *maxptr, FTYPE*minptr,FTYPE*sumptr)
{
  long long int i,j,k,pl;
  FTYPE maxsend,minsend,sumsend;
  int domin,domax,dosum;

  if(maxptr==NULL) domax=0; else domax=1;
  if(minptr==NULL) domin=0; else domin=1;
  if(sumptr==NULL) dosum=0; else dosum=1;
  
  for(pl=start;pl<start+nmemb;pl++){
    
    if(domin) minptr[pl]=1E30;
    if(domax) maxptr[pl]=-1E30;
    if(dosum) sumptr[pl]=0;
  }
  DUMPGENLOOP{
    for(pl=start;pl<start+nmemb;pl++){
      if(domax) if (MACP0A1(p,i,j,k,pl) > maxptr[pl]) maxptr[pl] = MACP0A1(p,i,j,k,pl);
      if(domin) if (MACP0A1(p,i,j,k,pl) < minptr[pl]) minptr[pl] = MACP0A1(p,i,j,k,pl);
      if(dosum) sumptr[pl]+=MACP0A1(p,i,j,k,pl);
    }
  }
#if(USEMPI)
  for(pl=start;pl<start+nmemb;pl++){    
    if(domax){
      maxsend = maxptr[pl];
      MPI_Allreduce(&maxsend, &maxptr[pl], 1, MPI_FTYPE, MPI_MAX, MPI_COMM_GRMHD);
    }
    if(domin){
      minsend = minptr[pl];
      MPI_Allreduce(&minsend, &minptr[pl], 1, MPI_FTYPE, MPI_MIN, MPI_COMM_GRMHD);
    }
    if(dosum){
      sumsend = sumptr[pl];
      MPI_Allreduce(&sumsend, &sumptr[pl], 1, MPI_FTYPE, MPI_SUM, MPI_COMM_GRMHD);
    }
  }
#endif
}





/// write to file or MPI buffer
void myfwrite(int bintxt, MPI_Datatype datatype, void *ptr, int start, int nmemb, int i, int j, int k, FILE**stream,void*writebuf)
{
  long long int pl;
  void *voidbuf;
  int sizeofdatatype;

  static long double *writebuf16;
  static double *writebuf8;
  static float *writebuf4;
  static unsigned char *writebuf1;
  static int *writebuf4i;
  static long long int *writebuf8i;

  static long double *ptr16;
  static double *ptr8;
  static float *ptr4;
  static unsigned char *ptr1;
  static int *ptr4i;
  static long long int *ptr8i;

  long long int streamnum;



  sizeofdatatype=getsizeofdatatype(datatype);

  if (datatype == MPI_UNSIGNED_CHAR) writebuf1 = (unsigned char *) writebuf;
  else if (datatype == MPI_FLOAT) writebuf4 = (float *) writebuf;
  else if (datatype == MPI_DOUBLE) writebuf8 = (double *) writebuf;
  else if (datatype == MPI_LONG_DOUBLE) writebuf16 = (long double *) writebuf;
  else if (datatype == MPI_INT) writebuf4i = (int *) writebuf;
  else if (datatype == MPI_LONG_LONG_INT) writebuf8i = (long long int *) writebuf;
  else{
    dualfprintf(fail_file,"No such datatype=%d\n",datatype);
    myexit(1);
  }

  if (datatype == MPI_UNSIGNED_CHAR) ptr1 = (unsigned char *) ptr;
  else if (datatype == MPI_FLOAT) ptr4 = (float *) ptr;
  else if (datatype == MPI_DOUBLE) ptr8 = (double *) ptr;
  else if (datatype == MPI_LONG_DOUBLE) ptr16 = (long double *) ptr;
  else if (datatype == MPI_INT) ptr4i = (int *) ptr;
  else if (datatype == MPI_LONG_LONG_INT) ptr8i = (long long int *) ptr;
  else{
    dualfprintf(fail_file,"No such datatype=%d\n",datatype);
    myexit(1);
  }

  if(mpicombine==0){
    if(bintxt==BINARYOUTPUT){ // binaryoutput
      if(!docolsplit){
        fwrite((char*)ptr+(sizeofdatatype*start), sizeofdatatype, nmemb, stream[0]);
        nextbuf+=nmemb;
      }
      else{
        for(pl=start;pl<start+nmemb;pl++){
          streamnum=nextbuf;
          fwrite((char*)ptr+(sizeofdatatype*pl), sizeofdatatype, 1, stream[streamnum]);
          nextbuf++;
        }
      }
    }
    else{ // text output
      for(pl=start;pl<start+nmemb;pl++){
        if(docolsplit) streamnum=nextbuf;
        else streamnum=0;
        if (datatype == MPI_UNSIGNED_CHAR) fprintf(stream[streamnum],"%c ",ptr1[pl]);
        else if (datatype == MPI_FLOAT) fprintf(stream[streamnum],"%15.7g ",ptr4[pl]);
        else if (datatype == MPI_DOUBLE) fprintf(stream[streamnum],"%21.15g ",ptr8[pl]);
        else if (datatype == MPI_LONG_DOUBLE) fprintf(stream[streamnum],"%31.25Lg ",ptr16[pl]);
        else if (datatype == MPI_INT) fprintf(stream[streamnum],"%d ",ptr4i[pl]);
        else if (datatype == MPI_LONG_LONG_INT) fprintf(stream[streamnum],"%lld ",ptr8i[pl]);
        nextbuf++;
      }
    }
  }
  else{ // mpicombine==1
    if(docolsplit&&USEROMIO){ // column splitting with ROMIO
      for(pl=start;pl<start+nmemb;pl++){
        if(nextbuf==romiocoliter){ // only write if doing that column
          // BUFFERMAP2 only depends on i,j, not column number
          if (datatype == MPI_UNSIGNED_CHAR) writebuf1[BUFFERMAP2] = ptr1[pl];
          if (datatype == MPI_FLOAT) writebuf4[BUFFERMAP2] = ptr4[pl];
          if (datatype == MPI_DOUBLE) writebuf8[BUFFERMAP2] = ptr8[pl];
          if (datatype == MPI_LONG_DOUBLE) writebuf16[BUFFERMAP2] = ptr16[pl];
          if (datatype == MPI_INT) writebuf4i[BUFFERMAP2] = ptr4i[pl];
          if (datatype == MPI_LONG_LONG_INT) writebuf8i[BUFFERMAP2] = ptr8i[pl];
        }
        nextbuf++;
      }
    }
    else{ // no ROMIO column splitting, just normal MPI buffer writing
      for(pl=start;pl<start+nmemb;pl++){
        if (datatype == MPI_UNSIGNED_CHAR) writebuf1[BUFFERMAP] = ptr1[pl];
        if (datatype == MPI_FLOAT) writebuf4[BUFFERMAP] = ptr4[pl];
        if (datatype == MPI_DOUBLE) writebuf8[BUFFERMAP] = ptr8[pl];
        if (datatype == MPI_LONG_DOUBLE) writebuf16[BUFFERMAP] = ptr16[pl];
        if (datatype == MPI_INT) writebuf4i[BUFFERMAP] = ptr4i[pl];
        if (datatype == MPI_LONG_LONG_INT) writebuf8i[BUFFERMAP] = ptr8i[pl];
      }
    }
  }
}

/// same kind of process as myfwrite, see comments there
void myfread(int bintxt, MPI_Datatype datatype, void *ptr, int start, int nmemb, int i, int j, int k, FILE**stream,void*writebuf)
{
  long long int pl;
  void *voidbuf;
  int sizeofdatatype;

  static long double *writebuf16;
  static double *writebuf8;
  static float *writebuf4;
  static unsigned char *writebuf1;
  static int *writebuf4i;
  static long long int *writebuf8i;

  static long double *ptr16;
  static double *ptr8;
  static float *ptr4;
  static unsigned char *ptr1;
  static int *ptr4i;
  static long long int *ptr8i;

  unsigned short short4char;

  int streamnum;



  sizeofdatatype=getsizeofdatatype(datatype);

  if (datatype == MPI_UNSIGNED_CHAR) writebuf1 = (unsigned char *) writebuf;
  else if (datatype == MPI_FLOAT) writebuf4 = (float *) writebuf;
  else if (datatype == MPI_DOUBLE) writebuf8 = (double *) writebuf;
  else if (datatype == MPI_LONG_DOUBLE) writebuf16 = (long double *) writebuf;
  else if (datatype == MPI_INT) writebuf4i = (int *) writebuf;
  else if (datatype == MPI_LONG_LONG_INT) writebuf8i = (long long int *) writebuf;
  else{
    dualfprintf(fail_file,"No such datatype=%d\n",datatype);
    myexit(1);
  }

  if (datatype == MPI_UNSIGNED_CHAR) ptr1 = (unsigned char *) ptr;
  else if (datatype == MPI_FLOAT) ptr4 = (float *) ptr;
  else if (datatype == MPI_DOUBLE) ptr8 = (double *) ptr;
  else if (datatype == MPI_LONG_DOUBLE) ptr16 = (long double *) ptr;
  else if (datatype == MPI_INT) ptr4i = (int *) ptr;
  else if (datatype == MPI_LONG_LONG_INT) ptr8i = (long long int *) ptr;
  else{
    dualfprintf(fail_file,"No such datatype=%d\n",datatype);
    myexit(1);
  }

  if(mpicombine==0){
    if(bintxt==BINARYOUTPUT){
      if(!docolsplit){
        fread((char*)ptr+(sizeofdatatype*start), sizeofdatatype, nmemb, stream[0]);
        nextbuf+=nmemb;
      }
      else{
        for(pl=start;pl<start+nmemb;pl++){
          streamnum=nextbuf;
          fread((char*)ptr+(sizeofdatatype*pl), sizeofdatatype, 1, stream[streamnum]);
          nextbuf++;
        }
      }
    }
    else{
      for(pl=start;pl<start+nmemb;pl++){
        if(docolsplit) streamnum=nextbuf;
        else streamnum=0;
        if (datatype == MPI_UNSIGNED_CHAR){
          fscanf(stream[streamnum],"%hu",&short4char);
          ptr1[pl]=short4char;
        }
        else if (datatype == MPI_FLOAT) fscanf(stream[streamnum],"%f",&ptr4[pl]);
        else if (datatype == MPI_DOUBLE) fscanf(stream[streamnum],"%lf",&ptr8[pl]);
        else if (datatype == MPI_LONG_DOUBLE) fscanf(stream[streamnum],"%Lf",&ptr16[pl]);
        else if (datatype == MPI_INT) fscanf(stream[streamnum],"%d",&ptr4i[pl]);
        else if (datatype == MPI_LONG_LONG_INT) fscanf(stream[streamnum],"%lld",&ptr8i[pl]);
        nextbuf++;
      }
    }
  }
  else{
    if(docolsplit&&USEROMIO){
      for(pl=start;pl<start+nmemb;pl++){
        if(nextbuf==romiocoliter){
          if (datatype == MPI_UNSIGNED_CHAR) ptr1[pl]=writebuf1[BUFFERMAP2]; 
          if (datatype == MPI_FLOAT) ptr4[pl]=writebuf4[BUFFERMAP2]; 
          if (datatype == MPI_DOUBLE) ptr8[pl]=writebuf8[BUFFERMAP2]; 
          if (datatype == MPI_LONG_DOUBLE) ptr16[pl]=writebuf16[BUFFERMAP2]; 
          if (datatype == MPI_INT) ptr4i[pl]=writebuf4i[BUFFERMAP2]; 
          if (datatype == MPI_LONG_LONG_INT) ptr8i[pl]=writebuf8i[BUFFERMAP2]; 
        }
        nextbuf++;
      }
    }
    else for(pl=start;pl<start+nmemb;pl++){
        if (datatype == MPI_UNSIGNED_CHAR) ptr1[pl]=writebuf1[BUFFERMAP]; 
        if (datatype == MPI_FLOAT) ptr4[pl]=writebuf4[BUFFERMAP]; 
        if (datatype == MPI_DOUBLE) ptr8[pl]=writebuf8[BUFFERMAP]; 
        if (datatype == MPI_LONG_DOUBLE) ptr16[pl]=writebuf16[BUFFERMAP]; 
        if (datatype == MPI_INT) ptr4i[pl]=writebuf4i[BUFFERMAP]; 
        if (datatype == MPI_LONG_LONG_INT) ptr8i[pl]=writebuf8i[BUFFERMAP]; 
      }
  }
}


/// sets values between 2 pointers, typically to cumulate values into writebuf array for later use.
void myset(MPI_Datatype datatype, void *ptr, int start, int nmemb, void*writebuf)
{
  long long int pl;

  static long double *writebuf16;
  static double *writebuf8;
  static float *writebuf4;
  static unsigned char *writebuf1;
  static int *writebuf4i;
  static long long int *writebuf8i;

  static long double *ptr16;
  static double *ptr8;
  static float *ptr4;
  static unsigned char *ptr1;
  static int *ptr4i;
  static long long int *ptr8i;

  if (datatype == MPI_UNSIGNED_CHAR) writebuf1 = (unsigned char *) writebuf;
  else if (datatype == MPI_FLOAT) writebuf4 = (float *) writebuf;
  else if (datatype == MPI_DOUBLE) writebuf8 = (double *) writebuf;
  else if (datatype == MPI_LONG_DOUBLE) writebuf16 = (long double *) writebuf;
  else if (datatype == MPI_INT) writebuf4i = (int *) writebuf;
  else if (datatype == MPI_LONG_LONG_INT) writebuf8i = (long long int *) writebuf;
  else{
    dualfprintf(fail_file,"No such datatype=%d\n",datatype);
    myexit(1);
  }

  if (datatype == MPI_UNSIGNED_CHAR) ptr1 = (unsigned char *) ptr;
  else if (datatype == MPI_FLOAT) ptr4 = (float *) ptr;
  else if (datatype == MPI_DOUBLE) ptr8 = (double *) ptr;
  else if (datatype == MPI_LONG_DOUBLE) ptr16 = (long double *) ptr;
  else if (datatype == MPI_INT) ptr4i = (int *) ptr;
  else if (datatype == MPI_LONG_LONG_INT) ptr8i = (long long int *) ptr;
  else{
    dualfprintf(fail_file,"No such datatype=%d\n",datatype);
    myexit(1);
  }

  for(pl=start;pl<start+nmemb;pl++){
    // nextcol is the iterator here, a global variable currently
    if (datatype == MPI_UNSIGNED_CHAR) writebuf1[nextcol++] = ptr1[pl];
    else if (datatype == MPI_FLOAT) writebuf4[nextcol++] = ptr4[pl];
    else if (datatype == MPI_DOUBLE) writebuf8[nextcol++] = ptr8[pl];
    else if (datatype == MPI_LONG_DOUBLE) writebuf16[nextcol++] = ptr16[pl];
    else if (datatype == MPI_INT) writebuf4i[nextcol++] = ptr4i[pl];
    else if (datatype == MPI_LONG_LONG_INT) writebuf8i[nextcol++] = ptr8i[pl];
  }
}



/// very similar to myset, just switched assignments
void myget(MPI_Datatype datatype, void *ptr, int start, int nmemb, void*writebuf)
{
  long long int pl;

  static long double *writebuf16;
  static double *writebuf8;
  static float *writebuf4;
  static unsigned char *writebuf1;
  static int *writebuf4i;
  static long long int *writebuf8i;

  static long double *ptr16;
  static double *ptr8;
  static float *ptr4;
  static unsigned char *ptr1;
  static int *ptr4i;
  static long long int *ptr8i;

  if (datatype == MPI_UNSIGNED_CHAR) writebuf1 = (unsigned char *) writebuf;
  else if (datatype == MPI_FLOAT) writebuf4 = (float *) writebuf;
  else if (datatype == MPI_DOUBLE) writebuf8 = (double *) writebuf;
  else if (datatype == MPI_LONG_DOUBLE) writebuf16 = (long double *) writebuf;
  else if (datatype == MPI_INT) writebuf4i = (int *) writebuf;
  else if (datatype == MPI_LONG_LONG_INT) writebuf8i = (long long int *) writebuf;
  else{
    dualfprintf(fail_file,"No such datatype=%d\n",datatype);
    myexit(1);
  }

  if (datatype == MPI_UNSIGNED_CHAR) ptr1 = (unsigned char *) ptr;
  else if (datatype == MPI_FLOAT) ptr4 = (float *) ptr;
  else if (datatype == MPI_DOUBLE) ptr8 = (double *) ptr;
  else if (datatype == MPI_LONG_DOUBLE) ptr16 = (long double *) ptr;
  else if (datatype == MPI_INT) ptr4i = (int *) ptr;
  else if (datatype == MPI_LONG_LONG_INT) ptr8i = (long long int *) ptr;
  else{
    dualfprintf(fail_file,"No such datatype=%d\n",datatype);
    myexit(1);
  }

  for(pl=start;pl<start+nmemb;pl++){
    // nextcol is the iterator here, a global variable currently
    if (datatype == MPI_UNSIGNED_CHAR) ptr1[pl] = writebuf1[nextcol++];
    else if (datatype == MPI_FLOAT) ptr4[pl] = writebuf4[nextcol++];
    else if (datatype == MPI_DOUBLE) ptr8[pl] = writebuf8[nextcol++];
    else if (datatype == MPI_LONG_DOUBLE) ptr16[pl] = writebuf16[nextcol++];
    else if (datatype == MPI_INT) ptr4i[pl] = writebuf4i[nextcol++];
    else if (datatype == MPI_LONG_LONG_INT) ptr8i[pl] = writebuf8i[nextcol++];
  }
}




int init_linklists(void)
{
  struct blink * blinkptr;
  struct blink * cpulinkptr;
  int i;
  long long int numlists;
  long long int numcells;
  int maxnumcolumns;




  /////////////////////
  //
  // Below shouldn't be modified for user purposes
  //
  // setup number of buffers
  //
  ///////////////////////

  maxnumcolumns=0;
  for(i=0;i<NUMDUMPTYPES;i++){
    if(maxnumcolumns<dnumcolumns[i]) maxnumcolumns=dnumcolumns[i];
  }
  // buffer must at least hold maxcolumns of data, and since buffer is only N1*N2*N3 big, make sure that at least NUMBUFFERS*N1*N2*N3>maxnumcolumns
  set_numbuffers(maxnumcolumns,&NUMBUFFERS);


  for(i=0;i<NUMDUMPTYPES;i++) if(dnumcolumns[i]>0) setuplinklist(dnumcolumns[i],i);


  trifprintf("end setuplinklists: %d\n",NUMDUMPTYPES);



  ///////////////////////////
  //
  // setup link lists for setuplinklist()
  //
  ///////////////////////////

  trifprintf("start per cpu lists\n");
  // check link lists
  for(i=0;i<NUMDUMPTYPES;i++){
    if(dnumcolumns[i]>0){
      logfprintf("i=%d\n",i);
      blinkptr=blinkptr0[i];
      numlists=0;
      numcells=0;
      while(blinkptr!=NULL){
        numcells+=blinkptr->num;
        //      logfprintf("i=%d num=%d, numtotal=%d\n",i,blinkptr->num,numcells);
        numlists++;
        blinkptr=blinkptr->np; // next one
      }
      logfprintf("i=%d numlists=%lld numcells=%lld\n",i,numlists,numcells);
      numlists=0;
    }
  }

  // check cpu=0 link list
  if(myid==0){
    trifprintf("start cpu==0 lists\n");
    for(i=0;i<NUMDUMPTYPES;i++){
      if(dnumcolumns[i]>0){
        logfprintf("i=%d\n",i);
        cpulinkptr=cpulinkptr0[i];
        numlists=0;
        numcells=0;
        while(cpulinkptr!=NULL){
          numcells+=cpulinkptr->num;
          // logfprintf("i=%d num=%d, cpu=%d, li=%d, lj=%d, lk=%d, col=%d, numtotal=%lld\n",i,cpulinkptr->num,cpulinkptr->cpu,cpulinkptr->i,cpulinkptr->j,cpulinkptr->k,cpulinkptr->col,numcells);
          numlists++;
          cpulinkptr=cpulinkptr->np; // next one
        }
        logfprintf("i=%d numlists=%lld numcells=%lld\n",i,numlists,numcells);
        numlists=0;
      }
    }
  }



  return(0);






}


/// setuplinklist() not a user function
/// use for mpiminio() functions (see init_mpi.c)
/// Had to convert to long long int for these variables and the blink structure in global.structs.h
/// This was required to deal with a large number of cpus, N1,N2,N3, and number of columns when totalsize[1]*totalsize[2]*totalsize[3]*numcolumns>2GB that is limit for size of signed integer that is default size on many systems
int setuplinklist(int numcolumns,int which)
{
  long long int gcount,lcount,numlinks;
  long long int i,j,k,col,li,lj,lk,pi,pj,pk,pid,firstlink;
  struct blink * clinkptr0, *clinkptr;
  struct blink * linkptr0for0, *linkptrfor0;
  long long int *lcountfor0;
  long long int firstlinkfor0;
  long long int *firstlijk,*li0,*lj0,*lk0,*lcol0;
  long long int ri,rj,rk,rcol;
  long long int *cpulist0;
  long long int numcpusinlist0,lcpu,itercpu,buffersize;
  long long int maxnumsize;


  set_maxnumsize(numcolumns,&maxnumsize);




  if(myid==0){
    // cpulist0's size is maximum possible number of cpus in a list due to buffer size
    //    buffersize=(int)(ceil(ceil((FTYPE)(N1*N2*N3*NUMBUFFERS)/(FTYPE)numcolumns)*(FTYPE)(numcolumns)/(FTYPE)N1));
    buffersize=numprocs;
    stderrfprintf("max cpus in a list=%lld\n",buffersize); fflush(stderr);
    if((cpulist0=(long long int*)malloc(sizeof(long long int)*buffersize))==NULL){
      dualfprintf(fail_file,"can't allocate cpulist0\n");
      myexit(10012);
    }
    if((lcountfor0=(long long int*)malloc(sizeof(long long int)*numprocs))==NULL){
      dualfprintf(fail_file,"can't allocate lcountfor0\n");
      myexit(10013);
    }
    if((firstlijk=(long long int*)malloc(sizeof(long long int)*numprocs))==NULL){
      dualfprintf(fail_file,"can't allocate firstlijk\n");
      myexit(10014);
    }
    if((li0=(long long int*)malloc(sizeof(long long int)*numprocs))==NULL){
      dualfprintf(fail_file,"can't allocate li0\n");
      myexit(10015);
    }
    if((lj0=(long long int*)malloc(sizeof(long long int)*numprocs))==NULL){
      dualfprintf(fail_file,"can't allocate lj0\n");
      myexit(10016);
    }
    if((lk0=(long long int*)malloc(sizeof(long long int)*numprocs))==NULL){
      dualfprintf(fail_file,"can't allocate lk0\n");
      myexit(10017);
    }
    if((lcol0=(long long int*)malloc(sizeof(long long int)*numprocs))==NULL){
      dualfprintf(fail_file,"can't allocate lcol0\n");
      myexit(10018);
    }
    for(i=0;i<buffersize;i++){
      cpulist0[i]=0;
    }
    for(i=0;i<numprocs;i++){
      lcountfor0[i]=firstlijk[i]=li0[i]=lj0[i]=lk0[i]=lcol0[i]=0;
    }
  }



  numcpusinlist0=0;

  clinkptr0=NULL;
  gcount=0;
  lcount=0;
  numlinks=0;
  firstlink=1;
  if(myid==0){
    for(itercpu=0;itercpu<numprocs;itercpu++){  firstlijk[itercpu]=1; }
    linkptr0for0=NULL;
    firstlinkfor0=1;
  }




  /////////////////////////
  // general loop
  for(k=0;k<ncpux3*N3;k++)  for(j=0;j<ncpux2*N2;j++)  for(i=0;i<ncpux1*N1;i++) for(col=0;col<numcolumns;col++){
          // relative local index
          li=i%N1;
          lj=j%N2;
          lk=k%N3;
          // cpu position number
          pi=(i/(long long int)N1);
          pj=(j/(long long int)N2);
          pk=(k/(long long int)N3);
          // cpu id for this data
          pid=pk*ncpux2*ncpux1+pj*ncpux1+pi;
          if(myid==pid) lcount++;
          if(myid==0){
            lcountfor0[pid]++;
            // below is if we have need this cpu's data (pid) and need to mark starting point on full grid
            if(firstlijk[pid]){
              cpulist0[numcpusinlist0++]=pid;
              li0[pid]=i;
              lj0[pid]=j;
              lk0[pid]=k;
              lcol0[pid]=col;
              if(col!=0){
                dualfprintf(fail_file,"col!=0 col=%lld, so chunking bad: numcolumns=%d which=%d\n",col,numcolumns,which);
                myexit(10019);
              }
              firstlijk[pid]=0;
            }
          }
          gcount++;
          //    if(myid==0){
          //  dualfprintf(fail_file,"%lld %lld %lld %lld\n",numcpusinlist0,gcount,pid,cpulist0[numcpusinlist0]); fflush(fail_file);
          // }
          //    logfprintf("%lld %lld %lld %lld %lld %lld %lld %lld\n",li,lj,lk,pi,pj,pk,pid,lcount,gcount);
          // 1st below if is to catch every buffer amount, while 2nd if part is needed to account for when the number of buffers is such that the last buffer isn't completely needed
          // this should work for any numcolumns or NUMBUFFERS, even at very last zone no matter what
          // chunk in minimum size of numcolumns
          if((gcount%(gcountmod(numcolumns))==0)||(gcount==(long long int)totalzones*(long long int)numcolumns)){
            // ok, so numcolumns can't exceed the buffer size, highly unlikely to happen, and checked for!
            if(myid==0){
              // must do in order determined to have data, not numerical order
              for(itercpu=0;itercpu<numcpusinlist0;itercpu++){
                lcpu=cpulist0[itercpu];
                if(lcountfor0[lcpu]>0){
                  if(itercpu==0){ // first cpu in list
                    ri=li0[lcpu];
                    rj=lj0[lcpu];
                    rk=lk0[lcpu];
                    rcol=lcol0[lcpu];
                  }
                  if(firstlinkfor0){
                    linkptrfor0=linkptr0for0=addlink(NULL);
                    firstlinkfor0=0;
                  }
                  else{
                    linkptrfor0=addlink(linkptrfor0);
                  }
                  linkptrfor0->cpu=lcpu;
                  linkptrfor0->num=lcountfor0[lcpu];
                  linkptrfor0->i=li0[lcpu];
                  linkptrfor0->j=lj0[lcpu];
                  linkptrfor0->k=lk0[lcpu];
                  linkptrfor0->col=lcol0[lcpu];
                  linkptrfor0->ri=ri;
                  linkptrfor0->rj=rj;
                  linkptrfor0->rk=rk;
                  linkptrfor0->rcol=rcol;
                  linkptrfor0->end=0;
     
                  lcountfor0[lcpu]=0; // reset counter for this id
                  firstlijk[lcpu]=1; // reset starting value
                }
                else{
                  dualfprintf(fail_file,"wtf: shoudn't be here.  Maybe passed more CPUs to batch system (mpirun) than passed to code?\n");
                  myexit(10020);
                }
              }
              // the last link is here identified as the last in the series of cpus to communicate with.  There's at least one new link here!
              linkptrfor0->end=1;
              numcpusinlist0=0; // reset list of cpus for this list
            }
            if(lcount>0){
              logfprintf("numcolumns=%d lcount=%lld\n",numcolumns,lcount); 
              // initialize another structure
              // set previous structure value to this structure, set this next one to NULL
              if(firstlink){
                clinkptr=clinkptr0=addlink(NULL);
                clinkptr->num=lcount;
                firstlink=0;
              }
              else{
                clinkptr=addlink(clinkptr);
                clinkptr->num=lcount;
              }
              lcount=0;
            }
          }// otherwise continue
        }      // now we have a link list for each cpu that determines how long each next buffer is that needs to be sent to cpu=0
  blinkptr0[which]=clinkptr0;
  cpulinkptr0[which]=linkptr0for0;

  return(0);
}

/// add link for forward-only link list
struct blink * addlink(struct blink * clinkptr)
{
  struct blink *pb;

  pb=(struct blink *)malloc(sizeof(struct blink));
  pb->np=NULL; // terminate list
  // set last link's pointer to this new structure
  if(clinkptr!=NULL) clinkptr->np=pb;

  return(pb);
}

