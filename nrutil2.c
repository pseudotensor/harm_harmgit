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

void free_f3matrix(FTYPE ***m, long nzl, long nzh, long nrl, long nrh, long ncl, long nch)
{
  long i,j;

  for(i=nzh;i>=nzl;i--){
    for (j = nrh; j >= nrl; j--){
      free((FTYPE *) (m[j][i] + ncl));
    }
    free((FTYPE *) (m[i] + nrl));
  }
  free((FTYPE *) (m + nzl));
}


