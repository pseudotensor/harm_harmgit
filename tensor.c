#include "decs.h"

// JCM: removed (unsigned) inside malloc since (unsigned)=(unsigned int) and argument was long int -- makes non sense and will truncate if asking for more than what integer number is

extern void nrerror(char error_text[]);

#define NR_END 1
#define FREE_ARG char*

FTYPE ***dtensor(long nrl, long nrh, long ncl, long nch, long ndl,
                 long ndh)
/* allocate a FTYPE 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] 
 */
{
  long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep =
    ndh - ndl + 1;
  FTYPE ***t;

  /* allocate pointers to pointers to rows */
  t = (FTYPE ***)
    malloc( ((nrow + NR_END) * sizeof(FTYPE **)));
  if (!t)
    nrerror("allocation failure 1 in tensor()");
  t += NR_END;
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl] = (FTYPE **)
    malloc( ((nrow * ncol + NR_END) *
             sizeof(FTYPE *)));
  if (!t[nrl])
    nrerror("allocation failure 2 in tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl] = (FTYPE *)
    malloc( ((nrow * ncol * ndep + NR_END) *
             sizeof(FTYPE)));
  if (!t[nrl][ncl])
    nrerror("allocation failure 3 in tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for (j = ncl + 1; j <= nch; j++)
    t[nrl][j] = t[nrl][j - 1] + ndep;
  for (i = nrl + 1; i <= nrh; i++) {
    t[i] = t[i - 1] + ncol;
    t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
    for (j = ncl + 1; j <= nch; j++)
      t[i][j] = t[i][j - 1] + ndep;
  }

  /* return pointer to array of pointers to rows */
  return t;
}

void free_dtensor(FTYPE ***t, long nrl, long nrh, long ncl, long nch,
                  long ndl, long ndh)
/* free a FTYPE tensor allocated by tensor() */
{
  free((FREE_ARG) (t[nrl][ncl] + ndl - NR_END));
  free((FREE_ARG) (t[nrl] + ncl - NR_END));
  free((FREE_ARG) (t + nrl - NR_END));
}
