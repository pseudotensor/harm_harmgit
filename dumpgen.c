#include "decs.h"

/*! \file dumpgen.c
  \brief General Data dump functions

  Data dump functions controlling how header + block of 3D data is outputted
*/


/// General dumping function
/// Used for all dump types for block data output
int dump_gen(int readwrite, long dump_cnt, int bintxt, int whichdump, MPI_Datatype datatype, char *fileprefix, char *fileformat, char *filesuffix, int (*headerfun) (int whichdump, int whichdumpversion, int numcolumns, int bintxt, FILE*headerptr),int (*setgetcontent) (int i, int j, int k, MPI_Datatype datatype, void*setbuf))
{
  int i = 0, j = 0, k = 0, l = 0, col = 0;
  FILE **fpp;
  char dfnam[MAXFILENAME]={'\0'};
  char dfnamreal[MAXFILENAME]={'\0'};
  char localfileformat[MAXFILENAME]={'\0'};
  void *jonio;
  static void *writebuf; // must be static so ROMIO non-blocking has permenant pointer address for this
  char truemyidtxt[MAXFILENAME]={'\0'};
  char filerw[MAXFILENAME]={'\0'};
  FILE *headerptr;
  int numfiles,coliter;
  void *setbuf;
  int sizeofdatatype;
  int romiocloopend;
  int headerbintxt;
  int dumpbintxt;
  fpos_t headerendpos;
  long int uptodatabytesize;
  int binextension;
  int checkstatus;
  int gopastlinebreak(FILE *stream);
  int whichdumpversion;



  // see if want to only do ROMIO finish-up when doing MPIVERSION>=2
  if(whichdump==FAKEDUMPTYPE){
    // only writebufptr and which  realy used, so just ignore other values
    // -1 indicates separate finish romio call, so that next call to mpiio_final() will not think some new file is ready to be finished.
    mpiio_final(0, 0, NULL, 0, MPICOMBINEROMIO, NULL, -1, 0, NULL, &writebuf);
    // done with finishing-up prior ROMIO combine if using ROMIO
    return(0);
  }


  ////////////
  //
  // setup file format for header and dump
  //
  ////////////

  if(readwrite==READFILE){
    strcpy(filerw,"rb");  //atch
    //if(bintxt==BINARYOUTPUT)
    //else strcpy(filerw,"rt");
  }
  else if(readwrite==WRITEFILE){
    strcpy(filerw,"wb");  //atch
    //if(bintxt==BINARYOUTPUT) strcpy(filerw,"w");
    //else strcpy(filerw,"wt");
  }

  if(bintxt==BINARYOUTPUT) headerbintxt=dumpbintxt=BINARYOUTPUT;
  else if(bintxt==TEXTOUTPUT) headerbintxt=dumpbintxt=TEXTOUTPUT;
  else if(bintxt==MIXEDOUTPUT){
    headerbintxt=TEXTOUTPUT;
    dumpbintxt=BINARYOUTPUT;
  }



  numcolumns=dnumcolumns[whichdump];
  whichdumpversion=dnumversion[whichdump];
  docolsplit=DOCOLSPLIT[whichdump]; // docolsplit global var for now


  ////////////////////
  //
  // See if enough HD space
  //
  ////////////////////

  if(mpicombine){
    if(dumpbintxt==BINARYOUTPUT){
      if(myid==0) isenoughfreespace((unsigned long long)(sizeof(FTYPE))*(unsigned long long)numcolumns*(unsigned long long)(totalsize[1])*(unsigned long long)(totalsize[2])*(unsigned long long)(totalsize[3]) );
      else isenoughfreespace(0);
    }
    else{// text
      if(myid==0) isenoughfreespace((unsigned long long)(22)*(unsigned long long)numcolumns*(unsigned long long)(totalsize[1])*(unsigned long long)(totalsize[2])*(unsigned long long)(totalsize[3]) );
      else isenoughfreespace(0);
    }
  }
  else{
    if(dumpbintxt==BINARYOUTPUT){
      isenoughfreespace((unsigned long long)(sizeof(FTYPE))*(unsigned long long)numcolumns*(unsigned long long)(N1)*(unsigned long long)(N2)*(unsigned long long)(N3) );
    }
    else{// text
      isenoughfreespace((unsigned long long)(21)*(unsigned long long)numcolumns*(unsigned long long)(N1)*(unsigned long long)(N2)*(unsigned long long)(N3) );
    }
  }


  /////////////////////
  //
  // Allocate memory for setting up setbuf
  //
  //////////////////////

  sizeofdatatype=getsizeofdatatype(datatype);
  setbuf=malloc(numcolumns*sizeofdatatype);
  if(setbuf==NULL){
    dualfprintf(fail_file,"cannot allocate memory to setbuf in %s %s with numcolumns=%d and sizeofdatatype=%d\n",fileprefix,filesuffix,numcolumns,sizeofdatatype);
    myexit(927656247);
  }


  //  trifprintf("numcolumns=%d sizeofdatatype=%d setbuf=%d total=%d\n",numcolumns,sizeofdatatype,setbuf,numcolumns*sizeofdatatype);

  //////////////////////////////////
  //
  //  Set up DOCOLSPLIT for normal and ROMIO loop
  //
  ///////////////////////////////////



  if(docolsplit){
    numfiles=numcolumns;
    if(mpicombine&&USEMPI&&USEROMIO) romiocloopend=numfiles;
    else romiocloopend=1;
  }
  else{
    numfiles=1;
    romiocloopend=1;
  }


  //////////////////////////////////
  //
  //  Define file output and open it
  //
  ///////////////////////////////////

  // say whether .bin is allowed or not if binary
  if(fileprefix[0]=='i') binextension=0; // images don't need binary extension
  else binextension=1;


  // sometimes all CPUs need to know filename (e.g. ROMIO)
  // setup file suffix
  if((dumpbintxt==BINARYOUTPUT)&&(binextension)){
    if(USEMPI&&(mpicombine==0)&&(numprocs>1)) sprintf(truemyidtxt,".bin.%04d",myid);
    else strcpy(truemyidtxt,".bin");
  }
  else{
    if(USEMPI&&(mpicombine==0)&&(numprocs>1)) sprintf(truemyidtxt,".%04d",myid);
    else strcpy(truemyidtxt,"");
  }

  // setup filename
  if(dump_cnt>=0){
    strcpy(localfileformat,"%s");
    strcat(localfileformat,fileformat);
    strcat(localfileformat,"%s");
    strcat(localfileformat,"%s");
    sprintf(dfnam, localfileformat, fileprefix, dump_cnt, filesuffix, truemyidtxt);
  }
  else{ // then no file number wanted (i.e. for gdump())
    sprintf(dfnam, "%s%s%s", fileprefix, filesuffix, truemyidtxt);
  }



  ////////////////
  //
  // open files, or open files for header if mpicombine==1, which for mpicombine==1 gets reopened later by MPI routines
  //
  ///////////////

  checkstatus=0; // no error so far
  int problemloadingfile=0;
  if((USEMPI&&(myid==0)&&(mpicombine==1))||(mpicombine==0)){// for mpicombine==1 even with ROMIO, real filename and header not needed

    // only one CPU does header if mpicombine==1, header+dump done in all CPUs if mpicombine==0
    // create files for each column, or each column's header if mpicombine==1
    if((fpp=(FILE**)malloc(sizeof(FILE*)*numfiles))==NULL){
      dualfprintf(fail_file,"couldn't open fpp in dump()\n");
      myexit(836565474);
    }// now fpp[i] indexes a list of file pointers


    // setup each file corresponding to each column
    COLLOOP(coliter){
      if(docolsplit&&(numfiles>1)){
        sprintf(dfnamreal,"%s-col%04d",dfnam,coliter);
      }
      else strcpy(dfnamreal,dfnam);
      
      if ((fpp[coliter] = fopen(dfnamreal, filerw)) == NULL) {
        dualfprintf(fail_file, "error opening %s %s (fullname=%s) file\n",fileprefix,filesuffix,dfnamreal);
        dualfprintf(fail_file, "Check if disk full or have correct permissions\n");
        problemloadingfile=1;
      }
    }// end COLLOOP
  }// end if myid or mpicombine==0 (split the loop so that can check for file existence first)


  // need to broadcast whether got all files or had problem
#if(USEMPI)
  MPI_Bcast(&problemloadingfile,1,MPI_INT,MPIid[0], MPI_COMM_GRMHD);
#endif
  if(problemloadingfile){
    // indicate failure, but one may wish to set some defaults if no file, so don't hard fail.
    return(FILENOTFOUND);
  }



  checkstatus=0; // no error so far
  if((USEMPI&&(myid==0)&&(mpicombine==1))||(mpicombine==0)){// for mpicombine==1 even with ROMIO, real filename and header not needed
    // setup each file corresponding to each column
    COLLOOP(coliter){

      //////////////////////////////////
      //
      //  read or write header: header is read/written in whatever style chosen to the top of each dump file created
      //
      ///////////////////////////////////
      if(headerfun!=NULL){
        headerfun(whichdump,whichdumpversion,numcolumns,headerbintxt,fpp[coliter]); // outputs header to each column file (or just one file, or all CPU files, etc.)
      }

      ////////////////////////////
      //
      // check that file one is reading is in right format (no need to check writing format)
      // assumed puts stream back to location before entered function
      //
      ////////////////////////////
      if(readwrite==READFILE){
        checkstatus=check_fileformat(readwrite, bintxt, whichdump, numcolumns, docolsplit, mpicombine, sizeofdatatype, fpp[coliter]);
      }

      // deal with transition between header and data
      if(readwrite==READFILE){
        if(bintxt==TEXTOUTPUT || bintxt==MIXEDOUTPUT){
          if(headerfun!=NULL){
            // now move past \n
            if(gopastlinebreak(fpp[coliter])) checkstatus=1;
          }
        }
      }
      // get position that would start real data
      uptodatabytesize=ftell(fpp[coliter]);


    }
    
    // don't close if mpicombine==0, since used in a moment, else mpicombine==1 it's reopened by MPI routines
    if(USEMPI&&(myid==0)&&(mpicombine==1)) COLLOOP(coliter) fclose(fpp[coliter]); // will get reopened later by MPI routines


  }



  // need to broadcast the header size to other CPUs for ROMIO
#if(USEMPI&&USEROMIO)
  MPI_Bcast(&uptodatabytesize,1,MPI_LONG,MPIid[0], MPI_COMM_GRMHD);
#endif


  ///////////////////////////////////////////////////////////
  //
  // Check file format status
  //
  //
  ////////////////////////////////////////////////////////////
  if(failed==0) failed=checkstatus;
  error_check(ERRORCODEBELOWCLEANFINISH+100); // number should be >ERRORCODEBELOWCLEANFINISH for myexit to avoid dumping
  


    
  ///////////////////////////////////////////////////////////
  //
  // loop over columns for per-column buffer ROMIO dump
  //
  //
  ////////////////////////////////////////////////////////////


  ROMIOCOLLOOP(romiocoliter) { // only loop if mpicombine&&USEMPI&&USEROMIO&&docolsplit==1
    if(romiocloopend>1) trifprintf("romiocoliter=%d of romiocloopend=%d\n",romiocoliter,romiocloopend);


    
    // setup MPI buffer if mpicombine==1
    if( mpicombine == 0 ) { // then one file per CPU if USEMPI or just normal file writing on 1CPU
      writebuf=NULL;
    }
    else mpiio_init(dumpbintxt,sortedoutput, fpp, uptodatabytesize, readwrite, dfnam, numcolumns, datatype, &jonio, &writebuf);
    // if USEROMIO==1 then numcolumns interpreted properly for docolsplit


    if(readwrite==READFILE){
      //////////////////////////////////
      //
      // read DUMP 
      //
      //////////////////////////////////
      
      if (mpicombine == 1) {
#if(USEMPI)
        mpiio_seperate(binaryoutput,sortedoutput, STAGE1, numcolumns, datatype, fpp, jonio, writebuf);
#endif
      }
    }



    //////////////////
    //
    // DUMP LOOP
    //
    //////////////////



    if(readwrite==READFILE){
      BUFFERINIT0;
      DUMPGENLOOP { // diagnostic loop
        // buffer init starts the parallel index
        BUFFERINIT;
        // initialize to 0th column
        COLINIT;

        ///////////////////////////////////////
        //
        // READFILE
        //
        //////////////////////
        if((mpicombine)&&(truempicombinetype==MPICOMBINEMINMEM)) mpiio_minmem(READFILE,whichdump,i,j,k,dumpbintxt,sortedoutput,numcolumns,datatype, fpp,jonio,writebuf);

 
        // read all at once
        myfread(dumpbintxt,datatype,setbuf,0,numcolumns,i,j,k,fpp,writebuf);
 
        // check
        if(nextbuf!=numcolumns){
          dualfprintf(fail_file,"Number of columns (numcolumns=%d) isn't equal to number of columns/buffers attempted (nextbuf=%lld)\n",numcolumns,nextbuf);
          myexit(932736466);
        }
 
        // get the content of 1 row
        setgetcontent(i,j,k,datatype,setbuf);
 
        // check
        if(nextcol!=numcolumns){
          dualfprintf(fail_file,"Number of columns (numcolumns=%d) isn't equal to number of columns attempted (nextcol=%d)\n",numcolumns,nextcol);
          myexit(836745613);
        }
      }// end DUMPGENLOOP
    }// end readwrite==READFILE
    else if(readwrite==WRITEFILE){
      BUFFERINIT0;
      DUMPGENLOOP { // diagnostic loop



        // buffer init starts the parallel index
        BUFFERINIT;
        // initialize to 0th column
        COLINIT;
        ///////////////////////////////////////
        //
        // WRITEFILE
        //
        //////////////////////

        // set the content of 1 row
        setgetcontent(i,j,k,datatype,setbuf);

        // check
        if(nextcol!=numcolumns){
          dualfprintf(fail_file,"Number of columns (numcolumns=%d) isn't equal to number of columns attempted (nextcol=%d)\n",numcolumns,nextcol);
          myexit(19785566);
        }

        // write all at once
        myfwrite(dumpbintxt,datatype,setbuf,0,numcolumns,i,j,k,fpp,writebuf);

 
        // check
        if(nextbuf!=numcolumns){
          dualfprintf(fail_file,"Number of columns (numcolumns=%d) isn't equal to number of columns/buffers attempted (nextbuf=%d)\n",numcolumns,nextbuf);
          myexit(94675455);
        }


        // finish up this row
        if((mpicombine==0)&&(dumpbintxt==TEXTOUTPUT)) COLLOOP(coliter) fprintf(fpp[coliter],"\n");
        if((mpicombine)&&(truempicombinetype==MPICOMBINEMINMEM)) mpiio_minmem(WRITEFILE,whichdump,i,j,k,dumpbintxt,sortedoutput,numcolumns,datatype, fpp, jonio,writebuf);



      }// end DUMPGENLOOP
    }//end readwrite==WRITEFILE
  



    //////////////////
    //
    // Close dump file and write/close file if mpicombine==1
    //
    //////////////////


    if (mpicombine == 0){
      COLLOOP(coliter) if (fpp[coliter] != NULL) fclose(fpp[coliter]);
    }
    else{
#if(USEMPI)
      if(readwrite==WRITEFILE) mpiio_combine(dumpbintxt, sortedoutput, numcolumns, datatype, fpp, jonio, writebuf);
      else if(readwrite==READFILE) mpiio_seperate(binaryoutput,sortedoutput, STAGE2, numcolumns, datatype, fpp, jonio, writebuf);
#endif
    }

  }// end column loop for ROMIO&&docolsplit



  // free the set/get buffer
  if(setbuf!=NULL) free(setbuf);



  return (0);
}









/// Output header for any dumping type
/// In reading/writing any header one has binary/text format
/// Need single function that read/write in binary/text so have consistent input/output format always
/// also need to Bcast it sometimes
///
/// examples:
/// binary read
///    fread(&idum1, sizeof(int), 1, headerptr);
///
/// text read:
///    fscanf(headerptr,"%ld",&DTr);
///
/// Bcast example:
/// MPI_Bcast(&avgscheme[1],1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
///
/// binary write:
///    fwrite(&totalsize[1], sizeof(int), 1, headerptr);
///
/// text write:
///    fprintf(headerptr,"%ld",DTr);
///
///
/// format assumed to have no space at end and will add that if writing
/// assume root=MPIid[0] and MPI_COMM_GRMHD for Bcast
/// readwrite:
///#define WRITEHEAD 0
///#define READHEAD 1
/// bintxt: BINARYOUTPUT TEXTOUTPUT only choices -- if mixed then choose text
int header1_gen(int accessmemory, int readwrite, int bintxt, int bcasthead, void *ptr, size_t size, char *format, size_t nmemb, MPI_Datatype datatype, FILE *stream)
{
  void *ptrlocal;
  //
  unsigned char *ptr1;
  float *ptr4;
  double *ptr8;
  long double *ptr16;
  int *ptr4i;
  long *ptr4l;
  long long int *ptr8i;
  //
  unsigned char var1;
  float var4;
  double var8;
  long double var16;
  int var4i;
  long var4l;
  long long int var8i;
  //
  char formatwithspace[MAXFILENAME]; // really format length, not file name
  size_t ii;
  char largedumbspace[200]; // GODMARK: Should be larger than long doubles
  int jj;

  
  if(accessmemory){
    ptrlocal=ptr;
  }
  else{

    // emulates argument passing to this function through assignment of pointers
    var1=0;
    var4=0.0;
    var8=0.0;
    var16=0.0;
    var4i=0;
    var4l=0;
    var8i=0;

    if (datatype == MPI_UNSIGNED_CHAR) ptrlocal = (unsigned char *) &var1;
    else if (datatype == MPI_FLOAT) ptrlocal = (float *) &var4;
    else if (datatype == MPI_DOUBLE) ptrlocal = (double *) &var8;
    else if (datatype == MPI_LONG_DOUBLE) ptrlocal = (long double *) &var16;
    else if (datatype == MPI_INT) ptrlocal = (int *) &var4i;
    else if (datatype == MPI_LONG) ptrlocal = (long *) &var4l;
    else if (datatype == MPI_LONG_LONG_INT) ptrlocal = (long long int *) &var8i;
    else{
      dualfprintf(fail_file,"No such datatype=%d\n",datatype);
      myexit(76293623);
    }

  }


  
  ///////////////////////
  //
  // resolve data type as necessary for fprintf
  //
  ///////////////////////
  if (datatype == MPI_UNSIGNED_CHAR) ptr1 = (unsigned char *) ptrlocal;
  else if (datatype == MPI_FLOAT) ptr4 = (float *) ptrlocal;
  else if (datatype == MPI_DOUBLE) ptr8 = (double *) ptrlocal;
  else if (datatype == MPI_LONG_DOUBLE) ptr16 = (long double *) ptrlocal;
  else if (datatype == MPI_INT) ptr4i = (int *) ptrlocal;
  else if (datatype == MPI_LONG) ptr4l = (long *) ptrlocal;
  else if (datatype == MPI_LONG_LONG_INT) ptr8i = (long long int *) ptrlocal;
  else{
    dualfprintf(fail_file,"No such datatype=%d\n",datatype);
    myexit(76293623);
  }




  if(readwrite==READHEAD){

    if(bintxt==BINARYOUTPUT){
      if(accessmemory) fread(ptrlocal,size,nmemb,stream);
      else{
        for(ii=0;ii<nmemb;ii++){
          fread(ptrlocal,size,1,stream); // just repeatedly fill same dummy space
        }
      }
    }
    else if(bintxt==TEXTOUTPUT || bintxt==MIXEDOUTPUT){

      for(ii=0;ii<nmemb;ii++){

        if(accessmemory) jj=ii;
        else jj=0; // repeatedly read into same dummy space

        if (datatype == MPI_UNSIGNED_CHAR) fscanf(stream,format,&ptr1[jj]);
        else if (datatype == MPI_FLOAT) fscanf(stream,format,&ptr4[jj]);
        else if (datatype == MPI_DOUBLE) fscanf(stream,format,&ptr8[jj]);
        else if (datatype == MPI_LONG_DOUBLE) fscanf(stream,format,&ptr16[jj]);
        else if (datatype == MPI_INT) fscanf(stream,format,&ptr4i[jj]);
        else if (datatype == MPI_LONG) fscanf(stream,format,&ptr4l[jj]);
        else if (datatype == MPI_LONG_LONG_INT) fscanf(stream,format,&ptr8i[jj]);
        else{
          dualfprintf(fail_file,"No such datatype=%d\n",datatype);
          myexit(40968321);
        }
      }
    }
    else{
      dualfprintf(fail_file,"No such bintxt=%d in readwrite=%d\n",bintxt,readwrite);
      myexit(249684962);
    }

  }
  else if(readwrite==WRITEHEAD){

    if(bintxt==BINARYOUTPUT){
      if(accessmemory) fwrite(ptrlocal,size,nmemb,stream);
      else{
        for(ii=0;ii<nmemb;ii++){
          fwrite(ptrlocal,size,1,stream); // repeatedly write from same dummy space (value=0)
        }
      }
    }
    else if(bintxt==TEXTOUTPUT || bintxt==MIXEDOUTPUT){
      sprintf(formatwithspace,"%s ",format);
      
      for(ii=0;ii<nmemb;ii++){

        if(accessmemory) jj=ii;
        else jj=0; // repeatedly read into same dummy space

        if (datatype == MPI_UNSIGNED_CHAR) fprintf(stream,formatwithspace,ptr1[jj]);
        else if (datatype == MPI_FLOAT) fprintf(stream,formatwithspace,ptr4[jj]);
        else if (datatype == MPI_DOUBLE) fprintf(stream,formatwithspace,ptr8[jj]);
        else if (datatype == MPI_LONG_DOUBLE) fprintf(stream,formatwithspace,ptr16[jj]);
        else if (datatype == MPI_INT) fprintf(stream,formatwithspace,ptr4i[jj]);
        else if (datatype == MPI_LONG) fprintf(stream,formatwithspace,ptr4l[jj]);
        else if (datatype == MPI_LONG_LONG_INT) fprintf(stream,formatwithspace,ptr8i[jj]);
        else{
          dualfprintf(fail_file,"No such datatype=%d\n",datatype);
          myexit(98346834);
        }
      }
    }
    else{
      dualfprintf(fail_file,"No such bintxt=%d in readwrite=%d\n",bintxt,readwrite);
      myexit(24934963);
    }

  }
  else if(readwrite==NOTHINGHEAD){
    if(!bcasthead){
      dualfprintf(fail_file,"Entered header1_gen() with nothing do to\n");
      myexit(24672672);
    }
  }
  else{
    dualfprintf(fail_file,"Entered header1_gen() with nothing do to version2\n");
    myexit(24672673);
  }




  // assume only broadcast when reading header information
  // only do if accessing memory (not dummy memory)
  if(bcasthead && accessmemory){
    // bintxt doesn't matter
    // assume root=MPIid[0] and MPI_COMM_GRMHD
#if(USEMPI)
    MPI_Bcast(ptr, (int)nmemb, datatype, MPIid[0], MPI_COMM_GRMHD);
#endif
    
  }



  //  return(0); // indicates no failure
  return(nmemb);

}






/// check that file read is in right format to avoid error in setup of restart header or data sizes
int check_fileformat(int readwrite, int bintxt, int whichdump, int numcolumnsvar, int docolsplitvar, int mpicombinevar, int sizeofdatatype, FILE *stream)
{
  long onlyheaderbytesize;
  long uptodatabytesize;
  long withintransitionbytesize;
  long totalbytesize;
  //
  int truenumcolumns;
  long long int datawordnumber;
  long long int totaldatasize;
  long long int wordtotal;
  long long int badwordtotal;
  long long int databytesize;
  long long int fullheaderbytesize;
  int get_word_count(long long int databytesize, long long int *wordtotal, FILE *stream);
  unsigned char mychar;
  int gopastlinebreak(FILE *stream);


  // get position of stream, which indicates size of header in bytes
  onlyheaderbytesize=ftell(stream);


  // find transition between header and data
  if(readwrite==READFILE){
    if(bintxt==TEXTOUTPUT || bintxt==MIXEDOUTPUT){
      gopastlinebreak(stream);
    }
  }
  // up to and including '\n' minus just beyond header
  // will be same as onlyheaderbytesize if nothing between header and data or if binary
  uptodatabytesize=ftell(stream);

  // get number of bytes up to first \n
  if(readwrite==READFILE){
    if(bintxt==TEXTOUTPUT || bintxt==MIXEDOUTPUT){
      fseek(stream,0,SEEK_SET);
      gopastlinebreak(stream);
    }
  }
  fullheaderbytesize=ftell(stream);


  // see if header is really up to first \n and not multiple \n
  if(fullheaderbytesize!=uptodatabytesize){
    dualfprintf(fail_file,"restart read found \\n mismatch: fullheaderbytesize=%lld uptodatabytesize=%lld\n",fullheaderbytesize,uptodatabytesize);
    return(1);
  }



  // go back to where was before getting up to first \n
  fseek(stream,uptodatabytesize,SEEK_SET);
  
  // byte size of transition region
  withintransitionbytesize = uptodatabytesize - onlyheaderbytesize;

  // DEBUG:
  //  dualfprintf(fail_file,"onlyheaderbytesize=%ld uptodatabytesize=%ld withintransitionbytesize=%ld\n",onlyheaderbytesize,uptodatabytesize,withintransitionbytesize);


  // go to end of file
  fseek(stream,0,SEEK_END);
  // get byte size
  totalbytesize=ftell(stream);


  // get bytes in data region
  databytesize = totalbytesize-uptodatabytesize;


  // determine number of columns within a file
  if(docolsplitvar) truenumcolumns=1;
  else truenumcolumns=numcolumnsvar;

  // determine number of words in data section
  if(mpicombinevar) datawordnumber=totalsize[1]*totalsize[2]*totalsize[3]*truenumcolumns;
  else datawordnumber=N1*N2*N3*truenumcolumns;

  // determine total bytes in data section
  totaldatasize=datawordnumber*sizeofdatatype;



  

  
  // only have checks for this case so far
  if(readwrite==READFILE && whichdump==RESTARTDUMPTYPE){


    ///////////////////////
    //
    // Check header
    //
    ///////////////////////
    if(bintxt==TEXTOUTPUT || bintxt==MIXEDOUTPUT){

      // first start back where ended header that will be first byte of data section
      fseek(stream,onlyheaderbytesize,SEEK_SET);

      // get word count from onlyheaderbytesize up to '\n'
      get_word_count(withintransitionbytesize, &badwordtotal, stream);

      if(badwordtotal>0){
        dualfprintf(fail_file,"restart read found extra words (badwordtotal=%lld) in header or could be that reading of header passed into data section\n",badwordtotal);
        return(1);
      }

    }



    ///////////////////////
    //
    // Check data
    //
    ///////////////////////
    if(bintxt==BINARYOUTPUT || bintxt==MIXEDOUTPUT){
      // header is binary for BINARYOUTPUT and text for MIXEDOUTPUT
      // data is binary

      // just check that data section is right size
      // this is easier since data section has a single data type unlike header
      // header has to be right size for this to work out
      if(databytesize != totaldatasize){
        dualfprintf(fail_file,"restart read binary header/data found databytesize=%d and totaldatasize=%d\n",databytesize,totaldatasize);
        return(1);
      }
      

    }
    else if(bintxt==TEXTOUTPUT){
      // header is text
      // data is text

      // in this case we don't grab \n just treating it as space and counting words as normal.  This is a more strict test that restart file is in correct format

      // for data as text could count words somehow (using wc and system) but many clusters don't allow system()
      // so have to do it manually

      // first start back where ended header that will be first byte of data section
      // don't use uptodatabytesize since want to catch errors in word count (want data to be exactly expected word count -- this checks that header is not too long)
      //fseek(stream,onlyheaderbytesize,SEEK_SET);
      fseek(stream,uptodatabytesize,SEEK_SET);  //ATCH: otherwise, would miss the last word in dump file since databytesize definition uses uptodatabytesize as data start pos

      // get word count
      get_word_count(databytesize, &wordtotal, stream);


      if(wordtotal!=datawordnumber){
        dualfprintf(fail_file,"restart read text data found wordtotal=%lld while datawordnumber=%lld\n",wordtotal,datawordnumber);
        dualfprintf(fail_file,"onlyheaderbytesize=%lld ,totalbytesize=%lld,databytesize=%lld,truenumcolumns=%lld,datawordnumber=%lld,totaldatasize=%lld\n",onlyheaderbytesize,totalbytesize,databytesize,truenumcolumns,datawordnumber,totaldatasize);
        return(1);
      }


    }
  }

  //////////////////////////
  //
  // whatever we did to the stream, put it back to just after header function called
  //
  //////////////////////////
  fseek(stream,onlyheaderbytesize,SEEK_SET);


  return(0);

}



/// move file pointer past line break
int gopastlinebreak(FILE *stream)
{
  int mychar;

  // then need to grab up to '\n'
  while(1){
    mychar=fgetc(stream);
    if(mychar=='\n' || mychar=='\r') break;
    if(feof(stream)){
      dualfprintf(fail_file,"Reached end of file while seeking \\n or \\r in header\n");
      return(1);
    }
  }
 
  return(0);
}

/// Get word count for data file to ensure correct size instead of assuming correct.
int get_word_count(long long int databytesize, long long int *wordtotal, FILE *stream)
{
  unsigned char mychar;
  int wordchar,spacechar;
  long long int i;
  int thischarisword;



  // now run through text data section
  // a word is defined as some characters not including space, \n, \r, etc. and being bounded at least on one side by such characters
  wordchar=0;
  spacechar=0;
  *wordtotal=0;
  // go through each byte
  for(i=0;i<databytesize;i++){
    // read-in a byte (character) at a time
    mychar=(unsigned char)fgetc(stream);


    if(feof(stream)){
      dualfprintf(fail_file,"Something is wrong with databytesize or loop since databytesize=%lld but found EOF\n",databytesize);
      return(1);
    }

    thischarisword=-1;
    if(mychar=='\n' || mychar=='\r' || mychar==' ' || mychar=='\t' || mychar=='\v'){
      // then got word-break (delimiter) type character
      spacechar++;
      thischarisword=0;
    }
    else{
      wordchar++;
      thischarisword=1;
    }


    // check for word at start
    if(i==0 && wordchar==1 && thischarisword){
      // then detected new word appearing and this is what we are counting (instead of ends of words)
      (*wordtotal)++;
    }

    // check for word entering as "space word"
    if(wordchar==1 && spacechar>0  && thischarisword){
      // this defines presence of word we just entered
      (*wordtotal)++;

      // reset spacechar
      spacechar=0;
    }
 
    // reset wordchar if within space region
    if(spacechar>0){
      wordchar=0;
    }

    // DEBUG:
    //    dualfprintf(fail_file,"i=%lld mychar=%c spacechar=%d wordchar=%d wordtotal=%lld\n",i,mychar,spacechar,wordchar,*wordtotal);


  }// end over all bytes


  return(0);
}
