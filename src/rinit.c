/* Registration of C routines */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void viterbi(double*, double*, double*, int*, int*, int*, int*, double*, double*, double*, double*, int*, double*);
void viterbi2(double*, double*, double*, int*, int*, int*, int*, double*, double*, double*, double*, double*, int*, double*);

#if _MSC_VER >= 1000
__declspec(dllexport)
#endif

static const R_CMethodDef R_CDef[] = {
  {"viterbi", (DL_FUNC)&viterbi,13},
  {"viterbi2", (DL_FUNC)&viterbi2,14},
  {NULL, NULL, 0},
};

void R_init_Test(DllInfo *info)
{
  R_registerRoutines(info,R_CDef,NULL,NULL,NULL);
}
