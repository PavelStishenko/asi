#ifndef blacsutils_h
#define blacsutils_h

#include <complex>

#define    DTYPE_             0                   /* Descriptor Type */
#define    CTXT_              1                     /* BLACS context */
#define    M_                 2             /* Global Number of Rows */
#define    N_                 3          /* Global Number of Columns */
#define    MB_                4                 /* Row Blocking Size */
#define    NB_                5              /* Column Blocking Size */
#define    RSRC_              6            /* Starting Processor Row */
#define    CSRC_              7         /* Starting Processor Column */
#define    LLD_               8           /* Local Leading Dimension */
#define    DLEN_              9                 /* Descriptor Length */

#define EXTERN extern "C"
#define double_complex_t std::complex<double>


EXTERN void blacs_pinfo_(int *mypnum, int *nprocs);
EXTERN void blacs_get_(int*, int*, int*);
EXTERN void blacs_gridinit_(int*, char*, int*, int*); // BLACS_GRIDINIT( ICONTXT, ORDER, NPROW, NPCOL )
EXTERN void blacs_barrier_(int *icontxt, char *scope);

EXTERN void dgesd2d_(int *ConTxt, int *m, int *n, double *A, int *lda, int *rdest, int *cdest);
EXTERN void dgerv2d_ (int *ConTxt, int *m, int *n, double *A, int *lda, int *rsrc, int *csrc);

// https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top/scalapack-routines/scalapack-redistribution-copy-routines/p-gemr2d.html
// http://www.netlib.org/scalapack/slug/node164.html
// http://www.netlib.org/scalapack/slug/node168.html
EXTERN void pdgemr2d_(int *m , int *n , double *a , int *ia , int *ja , int *desca , double *b , int *ib , int *jb , int *descb , int *ictxt);
EXTERN void pzgemr2d_(int *m , int *n , double_complex_t *a , int *ia , int *ja , int *desca , double_complex_t *b , int *ib , int *jb , int *descb , int *ictxt);

template <typename T>
void pXgemr2d(int *m , int *n , T *a , int *ia , int *ja , int *desca , T *b , int *ib , int *jb , int *descb , int *ictxt);

template <>
void pXgemr2d<double>          (int *m , int *n ,           double *a , int *ia , int *ja , int *desca ,           double *b , int *ib , int *jb , int *descb , int *ictxt)
{
  assert(*m>0);
  assert(*n>0);
  pdgemr2d_(m , n , a , ia , ja ,desca ,b ,ib ,jb ,descb ,ictxt);
}
template <>
void pXgemr2d<double_complex_t>(int *m , int *n , double_complex_t *a , int *ia , int *ja , int *desca , double_complex_t *b , int *ib , int *jb , int *descb , int *ictxt)
{
  assert(*m>0);
  assert(*n>0);
  pzgemr2d_(m , n , a , ia , ja ,desca ,b ,ib ,jb ,descb ,ictxt);
}

EXTERN void descinit_(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*); // DESC, M, N, MB, NB, IRSRC, ICSRC, ICTXT, LLD, INFO

EXTERN void blacs_gridinfo_(int *, int *, int *, int *, int *); // blacs_gridinfo (icontxt, nprow, npcol, myrow, mycol);
EXTERN int numroc_(int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);

EXTERN void ztpttr_(char *UPLO, int *N, double_complex_t *AP, double_complex_t *A, int *LDA, int *INFO);
EXTERN void dtpttr_(char *UPLO, int *N, double           *AP, double           *A, int *LDA, int *INFO);

template <typename T>
void Xtpttr(char *UPLO, int *N, T *AP, T *A, int *LDA, int *INFO);

template <>
void Xtpttr(char *UPLO, int *N, double *AP, double *A, int *LDA, int *INFO)
{
  dtpttr_(UPLO, N, AP, A, LDA, INFO);
}

template <>
void Xtpttr(char *UPLO, int *N, double_complex_t *AP, double_complex_t *A, int *LDA, int *INFO)
{
  ztpttr_(UPLO, N, AP, A, LDA, INFO);
}

int get_system_context(int blacs_context)
{
  int what=10, val;
  blacs_get_(&blacs_context, &what, &val); // get system context of the blacs_context
  return val;
}

int get_default_system_context()
{
  int system_context, zero=0;
  blacs_get_(&zero, &zero, &system_context); // get default system context
  return system_context;
}

/*
  init_context: system context
  Returns: new BLACS context
*/
int make_blacs_context(int init_context, int MP, int NP)
{
  char order = 'R';
  blacs_gridinit_(&init_context, &order, &MP, &NP);
 
  return init_context;
}

void blacs_desc_init(int M, int N, int blacs_context, int *desc)
{
  assert(M>0);
  assert(N>0);
  if (blacs_context  == -1)
  {
    memset(desc, 0, sizeof(int)*DLEN_);
    desc[DTYPE_] = 1;
    desc[CTXT_] = -1;
    return;    // stub for non-participating process
  }
  
  int nprow, npcol, myrow, mycol;
  blacs_gridinfo_(&blacs_context, &nprow, &npcol, &myrow, &mycol);  
  
  int MB=1, NB=1, rsrc=0, csrc=0;
  int info;
  
  int lld = numroc_(&M, &MB, &myrow, &rsrc, &nprow);
  if (lld < 1) 
  {
    lld = 1;
  }

  descinit_(desc, &M, &N, &MB, &NB,  &rsrc,  &csrc, &blacs_context,  &lld, &info); // DESC, M, N, MB, NB, IRSRC, ICSRC, ICTXT, LLD, INFO
  assert(!info);
}


int locrow(int *desc)
{
  int nprow, npcol, myrow, mycol;
  blacs_gridinfo_(&desc[CTXT_], &nprow, &npcol, &myrow, &mycol);
  
  return numroc_(&desc[M_], &desc[MB_], &myrow, &desc[RSRC_], &nprow);
}

int loccol(int *desc)
{
  int nprow, npcol, myrow, mycol;
  blacs_gridinfo_(&desc[CTXT_], &nprow, &npcol, &myrow, &mycol);
  
  return numroc_(&desc[N_], &desc[NB_], &mycol, &desc[CSRC_], &npcol);
}

#endif // blacsutils_h
