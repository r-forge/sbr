#include "R.h"
#include "lp_lib.h"
#include "lpSolveAPI.h"
#include "linprog.h"

void linprog(double* obj, double* A, int* pldA, int* pp, int* peq, double* b,
             double* lb, double* ub, int* pnInts, int* intvec,
             int* pPivotingRule, int* pSimplexType, int* basis,
             int* pepslevel, double* eps, int* ppresolve,
             double* pOptimalSolution, double* x, int* pStatus)
{
  lprec* lp = NULL;
  int i = 0, j = 0, ldA = -1, p = -1, eq = -1, nInts = -1, pivotingRule = -1;
  int simplexType = -1, epslevel = -1, presolve = -1, ret = -1023;
  unsigned char check = FALSE;
  double lpsInfinity = 0.0;
  int* rowno = NULL;

  ldA = *pldA;
  p = *pp;
  eq = *peq;
  nInts = *pnInts;
  pivotingRule = *pPivotingRule;
  simplexType = *pSimplexType;
  epslevel = *pepslevel;
  presolve = *ppresolve;

  rowno = (int*) R_alloc(ldA, sizeof(int));
  for(i = 0; i < ldA; i++)
    rowno[i] = i+1;
  
  lp = _make_lp(ldA, p);
  if(lp == NULL)
    error("unable to create linear program");

  _set_verbose(lp, 0);

  for(j = 0; j < p; j++) {
    _set_columnex(lp, j+1, ldA, A+j*ldA, rowno);
    _set_obj(lp, j+1, obj[j]);
  }    

  _set_rh_vec(lp, b);

  for(i = 1; i < eq; i++)
    _set_constr_type(lp, i, LE);

  for(i = eq; i <= ldA; i++)
    _set_constr_type(lp, i, EQ);

  lpsInfinity = _get_infinite(lp);
  for(j = 0; j < p; j++) {
    if(lb[j] == R_NegInf && ub[j] == R_PosInf)
      _set_unbounded(lp, j+1);
    else {
      if(ub[j] == R_PosInf)
        ub[j] = lpsInfinity;
      else if(lb[j] == R_NegInf)
        lb[j] = -1.0 * lpsInfinity;
      _set_bounds(lp, j+1, lb[j], ub[j]);
    }
  }

  for(j = 0; j < nInts; j++)
    _set_int(lp, intvec[j], TRUE);

  if(pivotingRule >= 0)
    _set_pivoting(lp, pivotingRule);

  if(simplexType >= 0)
    _set_simplextype(lp, simplexType);

  if(basis[0] == 0) {
    check = _set_basis(lp, basis, FALSE);
    if(check == FALSE)
      warning("user supplied basis was invalid - using phase 1 simplex");
  }

  if(epslevel >= 0) {
    _set_epslevel(lp, epslevel);
  }
  else {
    if(eps[0] > 0.0)
      _set_epsb(lp, eps[0]);
    if(eps[1] > 0.0)
      _set_epsd(lp, eps[1]);
    if(eps[2] > 0.0)
      _set_epsel(lp, eps[2]);
    if(eps[3] > 0.0)
      _set_epsint(lp, eps[3]);
    if(eps[4] > 0.0)
      _set_epsperturb(lp, eps[4]);
    if(eps[5] > 0.0)
      _set_epspivot(lp, eps[5]);
  }

  if(presolve > 0)
    _set_presolve(lp, presolve, _get_presolveloops(lp));

  ret = _solve(lp);

  eps[0] = _get_epsb(lp);
  eps[1] = _get_epsd(lp);
  eps[2] = _get_epsel(lp);
  eps[3] = _get_epsint(lp);
  eps[4] = _get_epsperturb(lp);
  eps[5] = _get_epspivot(lp);

  *pOptimalSolution = _get_objective(lp);
  _get_variables(lp, x);
  *pStatus = _get_status(lp);

  _delete_lp(lp);
}


