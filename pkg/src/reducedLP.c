#include "R.h"
#include "R_ext/BLAS.h"
#include "lp_lib.h"
#include "lpSolveAPI.h"
#include "reducedLP.h"

void reducedLP(double* A, int* pldA, int* pp, double* btilde, double *ub,
               double* beta, int* pPivotingRule, int* pStatus)
{
  int ldA = *pldA, p = *pp, i = 0, j = 0, pivotingRule = -1, status = -1;
  const int IONE = 1;
  double tolerance = 1e-9;
  double lpsInfinity = 0.0;
  double DONE = 1.0;
  lprec* lp = NULL;
  double* swapRow = (double*) R_alloc(p, sizeof(double));
  int* iwork = (int*) R_alloc(ldA-p-1, sizeof(int));

  pivotingRule = *pPivotingRule;

  for(i = 0; i < 2*p; i++)
    beta[i] = 0.0;

  for(i = 0; i < ldA-p; i++)
    iwork[i] = i+1;

  for(j = 0; j < p; j++) {
    i = j;
    while(fabs(A[i+j*ldA]) < tolerance) i++;
    if(i > j) {
      if(i >= ldA)
        error("zero column encountered during free variable pivoting");

      F77_CALL(dcopy)(&p, A+i+j*ldA, &ldA, swapRow, &IONE);
      F77_CALL(dcopy)(&p, A+j+j*ldA, &ldA, A+i+j*ldA, &ldA);
      F77_CALL(dcopy)(&p, swapRow, &IONE, A+j+j*ldA, &ldA);

    }
    pivot(A, ldA, 2*p, j, j, 0, 0);
  }

  lp = _make_lp(ldA-p-1, p);

  lpsInfinity = _get_infinite(lp);
  _set_verbose(lp, NEUTRAL);
  _set_pivoting(lp, pivotingRule);

  for(j = 1; j <= p; j++) {
    _set_columnex(lp, j, ldA-p-1, A+p+(j+p-1)*ldA, iwork);
    _set_obj(lp, j, A[(j+p)*ldA-1]);

    if(ub[j-1] == R_PosInf)
      ub[j-1] = lpsInfinity;
    _set_upbo(lp, j, ub[j-1]);
  }

  _set_rh_vec(lp, btilde);

  for(i = 1; i <= ldA-p-1; i++)
    _set_constr_type(lp, i, LE);

/*print_lp(lp);*/

  _solve(lp);
  status = _get_status(lp);
  
  if(status == 0) {
    _get_variables(lp, beta+p);

    for(j = p; j > 0; j--) {
      i = 2*p - j;
      beta[j-1] = -1.0 * F77_CALL(ddot)(&i, beta+j, &IONE, A+(j-1)+j*ldA, &ldA);
    }
  }

  *pStatus = status;

  _delete_lp(lp);
  lp = NULL;
}



/* void subst(double* A, int* pldA, int* pp, double* beta)
{
  int ldA = *pldA, p = *pp, i = 0, j = 0;
  const int IONE = 1;
  double tolerance = 1e-9;
  double DONE = 1.0;
  lprec* lp = NULL;
  double* swapRow = (double*) R_alloc(p, sizeof(double));
  int* iwork = (int*) R_alloc(ldA-p-1, sizeof(int));

  for(i = 0; i < 2*p; i++)
    beta[i] = 0.0;

  for(i = 0; i < ldA-p; i++)
    iwork[i] = i+1;

  for(j = 0; j < p; j++) {
    i = j;
    while(fabs(A[i+j*ldA]) < tolerance) i++;
    if(i > j) {
      if(i >= ldA)
        error("zero column encountered during freee variable pivoting");

      F77_CALL(dcopy)(&p, A+i+j*ldA, &ldA, swapRow, &IONE);
      F77_CALL(dcopy)(&p, A+j+j*ldA, &ldA, A+i+j*ldA, &ldA);
      F77_CALL(dcopy)(&p, swapRow, &IONE, A+j+j*ldA, &ldA);

    }
    pivot(A, ldA, 2*p, j, j, 0, 0);
  }

  lp = make_lp(ldA-p-1, ldA-1);
  set_verbose(lp, NEUTRAL);
  set_pivoting(lp, PRICER_FIRSTINDEX);

  for(j = 1; j <= p; j++) {
    set_columnex(lp, j, ldA-p-1, A+p+(j+p-1)*ldA, iwork);
    set_obj(lp, j, A[(j+p)*ldA-1]);
  }

  for(j = 1; j <= ldA-p-1; j++) {
    set_columnex(lp, j+p, 1, &DONE, &j);
  }

  for(j = 1; j <= ldA - 1; j++)
    set_upbo(lp, j, 1.0);

  for(j = 0; j <= ldA-p-2; j++)
    iwork[j] = -(j+p+1);
  set_basis(lp, iwork, FALSE);

//print_lp(lp);

  solve(lp);

  get_variables(lp, beta+p);

  for(j = p; j > 0; j--) {
    i = 2*p - j;
    beta[j-1] = -1.0 * F77_CALL(ddot)(&i, beta+j, &IONE, A+(j-1)+j*ldA, &ldA);
  }

  delete_lp(lp);
  lp = NULL;
}*/


void pivot(double* A, int ldA, int twop, int row, int col, int i, int j)
{
  int k = 0, m = 0;
  double alpha = 0.0, Aij = 1.0;

  i = i + row;
  j = j + col;

  Aij = A[i+j*ldA];

  for(k = col; k < j; k++)
    A[i+k*ldA] = A[i+k*ldA] / Aij;

  for(k = j+1; k < twop; k++) {
    A[i+k*ldA] = A[i+k*ldA] / Aij;
  }

  A[i+j*ldA] = 1.0;

  m = twop - col;

  for(k = row; k < i; k++) {
    alpha = -1.0 * A[k+j*ldA];
    F77_CALL(daxpy)(&m, &alpha, A+i+col*ldA, &ldA, A+k+col*ldA, &ldA);
  }

  for(k = i+1; k < ldA; k++) {
    alpha = -1.0 * A[k+j*ldA];
    F77_CALL(daxpy)(&m, &alpha, A+i+col*ldA, &ldA, A+k+col*ldA, &ldA);
  }
}
