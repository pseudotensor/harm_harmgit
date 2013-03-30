#include "decs.h"


// directly from interp.c

///////////////////////////////////////
//
// MEMORY RELATED FUNCTIONS AND IMAGE WRITE FUNCTION
//
////////////////////////////////////////

unsigned char **cmatrix(nrl,nrh,ncl,nch)
     int nrl,nrh,ncl,nch;
{
  unsigned char **m;
  int i ;

  m=(unsigned char **)malloc((unsigned) (nrh-nrl+1)*sizeof(unsigned char *));
  if (!m) exit(20) ;
  m -= nrl;

  for(i=nrl;i<=nrh;i++) {
    m[i]=(unsigned char *)malloc((unsigned) (nch-ncl+1)*sizeof(unsigned char));
    if (!m[i]) exit(20) ;
    m[i] -= ncl;
  }
  return m;
}

FTYPE **fmatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
  FTYPE **m;
  int i ;

  m=(FTYPE **)malloc((unsigned) (nrh-nrl+1)*sizeof(FTYPE *));
  if (!m) exit(20) ;
  m -= nrl;

  for(i=nrl;i<=nrh;i++) {
    m[i]=(FTYPE *)malloc((unsigned) (nch-ncl+1)*sizeof(FTYPE));
    if (!m[i]) exit(20) ;
    m[i] -= ncl;
  }
  return m;
}
FTYPE ***f3matrix(nzl,nzh,nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch,nzl,nzh;
{
  FTYPE ***m;
  int i,j ;

  m=(FTYPE ***)malloc((unsigned) (nzh-nzl+1)*sizeof(FTYPE **));
  if (!m) exit(20) ;
  m -= nzl;

  for(i=nzl;i<=nzh;i++){
    m[i]=(FTYPE **)malloc((unsigned) (nrh-nrl+1)*sizeof(FTYPE *));
    if (!m[i]) exit(20) ;
    m[i] -= nrl;

    for(j=nrl;j<=nrh;j++) {
      m[i][j]=(FTYPE *)malloc((unsigned) (nch-ncl+1)*sizeof(FTYPE));
      if (!m[i][j]) exit(20) ;
      m[i][j] -= ncl;
    }
    // m[i][j][k]=m[z][r][c]

  }

  return m;
}
FTYPE ****f4matrix(nql,nqh,nzl,nzh,nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch,nzl,nzh,nql,nqh;
{
  FTYPE ****m;
  int k,i,j ;

  m=(FTYPE ****)malloc((unsigned) (nqh-nql+1)*sizeof(FTYPE ***));
  if (!m) exit(20) ;
  m -= nql;

  for(k=nql;k<=nqh;k++){

    m[k]=(FTYPE ***)malloc((unsigned) (nzh-nzl+1)*sizeof(FTYPE **));
    if (!m[k]) exit(20) ;
    m[k] -= nzl;
   
    for(i=nzl;i<=nzh;i++){
      m[k][i]=(FTYPE **)malloc((unsigned) (nrh-nrl+1)*sizeof(FTYPE *));
      if (!m[k][i]) exit(20) ;
      m[k][i] -= nrl;
     
      for(j=nrl;j<=nrh;j++) {
        m[k][i][j]=(FTYPE *)malloc((unsigned) (nch-ncl+1)*sizeof(FTYPE));
        if (!m[k][i][j]) exit(20) ;
        m[k][i][j] -= ncl;
      }
      // m[h][i][j][k]=m[q][z][r][c]
    }
   
  }
  return m;
}

FTYPE *****f5matrix(ncoll,ncolh,nql,nqh,nzl,nzh,nrl,nrh,ncl,nch)
int ncoll,ncolh,nrl,nrh,ncl,nch,nzl,nzh,nql,nqh;
{
  FTYPE *****m;
  int k,i,j,l ;

  m=(FTYPE *****)malloc((unsigned) (ncolh-ncoll+1)*sizeof(FTYPE ****));
  if (!m) exit(20) ;
  m -= ncoll;

  for(k=ncoll;k<=ncolh;k++){

    m[k]=(FTYPE ****)malloc((unsigned) (nqh-nql+1)*sizeof(FTYPE ***));
    if (!m[k]) exit(20) ;
    m[k] -= nql;
   
    for(i=nql;i<=nqh;i++){
      m[k][i]=(FTYPE ***)malloc((unsigned) (nzh-nzl+1)*sizeof(FTYPE **));
      if (!m[k][i]) exit(20) ;
      m[k][i] -= nzl;
     
      for(j=nzl;j<=nzh;j++) {
        m[k][i][j]=(FTYPE **)malloc((unsigned) (nrh-nrl+1)*sizeof(FTYPE*));
        if (!m[k][i][j]) exit(20) ;
        m[k][i][j] -= nrl;
       
        for(l=nrl;l<=nrh;l++) {
          m[k][i][j][l]=(FTYPE *)malloc((unsigned) (nch-ncl+1)*sizeof(FTYPE));
          if (!m[k][i][j][l]) exit(20) ;
          m[k][i][j][l] -= ncl;
        }
        // m[col][h][i][j][k]=m[col][q][z][r][c]=m[k][i][j][l][.]
      }
    }
  }
  return m;
}

unsigned char ***c3matrix(nzl,nzh,nrl,nrh,ncl,nch)
     int nrl,nrh,ncl,nch,nzl,nzh;
{
  unsigned char ***m;
  int i,j ;

  m=(unsigned char ***)malloc((unsigned) (nzh-nzl+1)*sizeof(unsigned char **));
  if (!m) exit(20) ;
  m -= nzl;

  for(i=nzl;i<=nzh;i++){
    m[i]=(unsigned char **)malloc((unsigned) (nrh-nrl+1)*sizeof(unsigned char *));
    if (!m[i]) exit(20) ;
    m[i] -= nrl;

    for(j=nrl;j<=nrh;j++) {
      m[i][j]=(unsigned char *)malloc((unsigned) (nch-ncl+1)*sizeof(unsigned char));
      if (!m[i][j]) exit(20) ;
      m[i][j] -= ncl;
    }
    // m[i][j][k]=m[z][r][c]

  }

  return m;
}
unsigned char ****c4matrix(nql,nqh,nzl,nzh,nrl,nrh,ncl,nch)
     int nrl,nrh,ncl,nch,nzl,nzh,nql,nqh;
{
  unsigned char ****m;
  int i,j,k ;


  m=(unsigned char ****)malloc((unsigned) (nqh-nql+1)*sizeof(unsigned char ***));
  if (!m) exit(20) ;
  m -= nql;
 
  for(k=nql;k<=nqh;k++){
   
    m[k]=(unsigned char ***)malloc((unsigned) (nzh-nzl+1)*sizeof(unsigned char **));
    if (!m[k]) exit(20) ;
    m[k] -= nzl;
   
    for(i=nzl;i<=nzh;i++){
      m[k][i]=(unsigned char **)malloc((unsigned) (nrh-nrl+1)*sizeof(unsigned char *));
      if (!m[k][i]) exit(20) ;
      m[k][i] -= nrl;
     
      for(j=nrl;j<=nrh;j++) {
        m[k][i][j]=(unsigned char *)malloc((unsigned) (nch-ncl+1)*sizeof(unsigned char));
        if (!m[k][i][j]) exit(20) ;
        m[k][i][j] -= ncl;
      }
      // m[h][i][j][k]=m[q][z][r][c]

    }
  }
 
  return m;
}
unsigned char *****c5matrix(ncoll,ncolh,nql,nqh,nzl,nzh,nrl,nrh,ncl,nch)
     int ncoll,ncolh,nrl,nrh,ncl,nch,nzl,nzh,nql,nqh;
{
  unsigned char *****m;
  int i,j,k,l ;


  m=(unsigned char *****)malloc((unsigned) (ncolh-ncoll+1)*sizeof(unsigned char ****));
  if (!m) exit(20) ;
  m -= ncoll;
 
  for(k=ncoll;k<=ncolh;k++){
   
    m[k]=(unsigned char ****)malloc((unsigned) (nqh-nql+1)*sizeof(unsigned char ***));
    if (!m[k]) exit(20) ;
    m[k] -= nql;
   
    for(i=nql;i<=nqh;i++){
      m[k][i]=(unsigned char ***)malloc((unsigned) (nzh-nzl+1)*sizeof(unsigned char **));
      if (!m[k][i]) exit(20) ;
      m[k][i] -= nzl;
     
      for(j=nzl;j<=nzh;j++) {
        m[k][i][j]=(unsigned char **)malloc((unsigned) (nrh-nrl+1)*sizeof(unsigned char*));
        if (!m[k][i][j]) exit(20) ;
        m[k][i][j] -= nrl;

        for(l=nrl;l<=nrh;l++) {
          m[k][i][j][l]=(unsigned char *)malloc((unsigned) (nch-ncl+1)*sizeof(unsigned char));
          if (!m[k][i][j][l]) exit(20) ;
          m[k][i][j][l] -= ncl;
        }
      }
      // m[col][h][i][j][k]=m[col][q][z][r][c]=m[k][i][j][l][.]

    }
  }
 
  return m;
}

void free_cmatrix(unsigned char **m, long nrl, long nrh, long ncl, long nch)
{
  long i;

  for (i = nrh; i >= nrl; i--)
    free((unsigned char *) (m[i] + ncl));
  free((unsigned char *) (m + nrl));
}

void free_fmatrix(FTYPE **m, long nrl, long nrh, long ncl, long nch)
{
  long i;

  for (i = nrh; i >= nrl; i--)
    free((FTYPE *) (m[i] + ncl));
  free((FTYPE *) (m + nrl));
}

void free_c3matrix(unsigned char ***m, long nzl, long nzh, long nrl, long nrh, long ncl, long nch)
{
  long i,j;

  for(i=nzh;i>=nzl;i--){
    for (j = nrh; j >= nrl; j--){
      free((unsigned char *) (m[i][j] + ncl));
    }
    free((unsigned char *) (m[i] + nrl));
  }
  free((unsigned char *) (m + nzl));
}

void free_f3matrix(FTYPE ***m, long nzl, long nzh, long nrl, long nrh, long ncl, long nch)
{
  long i,j;

  for(i=nzh;i>=nzl;i--){
    for (j = nrh; j >= nrl; j--){
      free((FTYPE *) (m[i][j] + ncl));
    }
    free((FTYPE *) (m[i] + nrl));
  }
  free((FTYPE *) (m + nzl));
}

void free_c4matrix(unsigned char ****m, long nql, long nqh, long nzl, long nzh, long nrl, long nrh, long ncl, long nch)
{
  long i,j,k;

  for(k=nqh;k>=nql;k--){
    for(i=nzh;i>=nzl;i--){
      for (j = nrh; j >= nrl; j--){
        free((unsigned char *) (m[k][i][j] + ncl));
      }
      free((unsigned char *) (m[k][i] + nrl));
    }
    free((unsigned char *) (m[k] + nzl));
  }
  free((unsigned char *) (m + nql));
}

void free_f4matrix(FTYPE ****m, long nql, long nqh, long nzl, long nzh, long nrl, long nrh, long ncl, long nch)
{
  long i,j,k;

  for(k=nqh;k>=nql;k--){
    for(i=nzh;i>=nzl;i--){
      for (j = nrh; j >= nrl; j--){
        free((FTYPE *) (m[k][i][j] + ncl));
      }
      free((FTYPE *) (m[k][i] + nrl));
    }
    free((FTYPE *) (m[k] + nzl));
  }
  free((FTYPE *) (m + nql));
}

void free_c5matrix(unsigned char *****m, long ncoll, long ncolh, long nql, long nqh, long nzl, long nzh, long nrl, long nrh, long ncl, long nch)
{
  long i,j,k,l;

  for(l=ncolh;l>=ncoll;l--){
    for(k=nqh;k>=nql;k--){
      for(i=nzh;i>=nzl;i--){
        for (j = nrh; j >= nrl; j--){
          free((unsigned char *) (m[l][k][i][j] + ncl));
        }
        free((unsigned char *) (m[l][k][i] + nrl));
      }
      free((unsigned char *) (m[l][k] + nzl));
    }
    free((unsigned char *) (m[l] + nql));
  }
  free((unsigned char *) (m + ncoll));
}

void free_f5matrix(FTYPE *****m, long ncoll, long ncolh, long nql, long nqh, long nzl, long nzh, long nrl, long nrh, long ncl, long nch)
{
  long i,j,k,l;

  for(l=ncolh;l>=ncoll;l--){
    for(k=nqh;k>=nql;k--){
      for(i=nzh;i>=nzl;i--){
        for (j = nrh; j >= nrl; j--){
          free((FTYPE *) (m[l][k][i][j] + ncl));
        }
        free((FTYPE *) (m[l][k][i] + nrl));
      }
      free((FTYPE *) (m[l][k] + nzl));
    }
    free((FTYPE *) (m[l] + nql));
  }
  free((FTYPE *) (m + ncoll));
}





