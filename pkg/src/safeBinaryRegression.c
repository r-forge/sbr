#include "R.h"
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include "lp_lib.h"
#include "safeBinaryRegression.h"
#include "linprog.h"
#include "reducedLP.h"

make_lp_func *_make_lp;
set_verbose_func *_set_verbose;
set_columnex_func *_set_columnex;
set_obj_func *_set_obj;
set_rh_vec_func *_set_rh_vec;
set_constr_type_func *_set_constr_type;
get_infinite_func *_get_infinite;
set_unbounded_func *_set_unbounded;
set_bounds_func *_set_bounds;
set_upbo_func *_set_upbo;
set_int_func *_set_int;
set_pivoting_func *_set_pivoting;
set_simplextype_func *_set_simplextype;
set_basis_func *_set_basis;
set_epslevel_func *_set_epslevel;
set_epsb_func *_set_epsb;
set_epsd_func *_set_epsd;
set_epsel_func *_set_epsel;
set_epsint_func *_set_epsint;
set_epsperturb_func *_set_epsperturb;
set_epspivot_func *_set_epspivot;
set_presolve_func *_set_presolve;
get_presolveloops_func *_get_presolveloops;
solve_func *_solve;
get_epsb_func *_get_epsb;
get_epsd_func *_get_epsd;
get_epsel_func *_get_epsel;
get_epsint_func *_get_epsint;
get_epsperturb_func *_get_epsperturb;
get_epspivot_func *_get_epspivot;
get_objective_func *_get_objective;
get_variables_func *_get_variables;
get_status_func *_get_status;
delete_lp_func *_delete_lp;


void R_init_safeBinaryRegression(DllInfo *info)
{
  const char lpSolveAPI[] = "lpSolveAPI";

  R_NativePrimitiveArgType linprogAgrs[19] = {REALSXP, REALSXP, INTSXP, INTSXP,
                                              INTSXP, REALSXP, REALSXP, REALSXP,
                                              INTSXP, INTSXP, INTSXP, INTSXP,
                                              INTSXP, INTSXP, REALSXP, INTSXP,
                                              REALSXP, REALSXP, INTSXP};

  R_NativePrimitiveArgType reducedLpAgrs[8] = {REALSXP, INTSXP, INTSXP, REALSXP,
                                               REALSXP, REALSXP, INTSXP,
                                               INTSXP};

  void reducedLP(double* A, int* pldA, int* pp, double* btilde, double *ub,
               double* beta, int* pPivotingRule, int* pStatus);


  R_CMethodDef cMethods[] = {
    {"linprog", (DL_FUNC) &linprog, 19, linprogAgrs},
    {"reducedLP", (DL_FUNC) &reducedLP, 8, reducedLpAgrs},
    {NULL, NULL, 0}
  };

  R_registerRoutines(info, cMethods, NULL, NULL, NULL);

  _make_lp = (make_lp_func*) R_GetCCallable(lpSolveAPI, "make_lp");
  _set_verbose = (set_verbose_func*) R_GetCCallable(lpSolveAPI, "set_verbose");
  _set_columnex = (set_columnex_func*) R_GetCCallable(lpSolveAPI, "set_columnex");
  _set_obj = (set_obj_func*) R_GetCCallable(lpSolveAPI, "set_obj");
  _set_rh_vec = (set_rh_vec_func*) R_GetCCallable(lpSolveAPI, "set_rh_vec");
  _set_constr_type = (set_constr_type_func*) R_GetCCallable(lpSolveAPI, "set_constr_type");
  _get_infinite = (get_infinite_func*) R_GetCCallable(lpSolveAPI, "get_infinite");
  _set_unbounded = (set_unbounded_func*) R_GetCCallable(lpSolveAPI, "set_unbounded");
  _set_bounds = (set_bounds_func*) R_GetCCallable(lpSolveAPI, "set_bounds");
  _set_upbo = (set_upbo_func*) R_GetCCallable(lpSolveAPI, "set_upbo");
  _set_int = (set_int_func*) R_GetCCallable(lpSolveAPI, "set_int");
  _set_pivoting = (set_pivoting_func*) R_GetCCallable(lpSolveAPI, "set_pivoting");
  _set_simplextype = (set_simplextype_func*) R_GetCCallable(lpSolveAPI, "set_simplextype");
  _set_basis = (set_basis_func*) R_GetCCallable(lpSolveAPI, "set_basis");
  _set_epslevel = (set_epslevel_func*) R_GetCCallable(lpSolveAPI, "set_epslevel");
  _set_epsb = (set_epsb_func*) R_GetCCallable(lpSolveAPI, "set_epsb");
  _set_epsd = (set_epsd_func*) R_GetCCallable(lpSolveAPI, "set_epsd");
  _set_epsel = (set_epsel_func*) R_GetCCallable(lpSolveAPI, "set_epsel");
  _set_epsint = (set_epsint_func*) R_GetCCallable(lpSolveAPI, "set_epsint");
  _set_epsperturb = (set_epsperturb_func*) R_GetCCallable(lpSolveAPI, "set_epsperturb");
  _set_epspivot = (set_epspivot_func*) R_GetCCallable(lpSolveAPI, "set_epspivot");
  _set_presolve = (set_presolve_func*) R_GetCCallable(lpSolveAPI, "set_presolve");
  _get_presolveloops = (get_presolveloops_func*) R_GetCCallable(lpSolveAPI, "get_presolveloops");
  _solve = (solve_func*) R_GetCCallable(lpSolveAPI, "solve");
  _get_epsb = (get_epsb_func*) R_GetCCallable(lpSolveAPI, "get_epsb");
  _get_epsd = (get_epsd_func*) R_GetCCallable(lpSolveAPI, "get_epsd");
  _get_epsel = (get_epsel_func*) R_GetCCallable(lpSolveAPI, "get_epsel");
  _get_epsint = (get_epsint_func*) R_GetCCallable(lpSolveAPI, "get_epsint");
  _get_epsperturb = (get_epsperturb_func*) R_GetCCallable(lpSolveAPI, "get_epsperturb");
  _get_epspivot = (get_epspivot_func*) R_GetCCallable(lpSolveAPI, "get_epspivot");
  _get_objective = (get_objective_func*) R_GetCCallable(lpSolveAPI, "get_objective");
  _get_variables = (get_variables_func*) R_GetCCallable(lpSolveAPI, "get_variables");
  _get_status = (get_status_func*) R_GetCCallable(lpSolveAPI, "get_status");
  _delete_lp = (delete_lp_func*) R_GetCCallable(lpSolveAPI, "delete_lp");
}


