#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .C calls */
extern void forward_backward(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void simulate_markov(void *, void *, void *, void *, void *, void *);
extern void viterbi(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"forward_backward", (DL_FUNC) &forward_backward, 18},
  {"simulate_markov",  (DL_FUNC) &simulate_markov,   6},
  {"viterbi",          (DL_FUNC) &viterbi,          13},
  {NULL, NULL, 0}
};

void R_init_hhsmm(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
