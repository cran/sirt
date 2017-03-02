#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP sirt_calc_copula_itemcluster_C(SEXP, SEXP);
extern SEXP sirt_calccounts_pcm_groups_C(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_choppin_rowaveraging(SEXP, SEXP, SEXP);
extern SEXP sirt_eigenvaluesDsirt(SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_evm_aux(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_evm_comp_matrix_poly(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_firsteigenvalsirt2(SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_gooijer_csn_table(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_ia_optim_lambda(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_ia_optim_nu(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_interval_index_C(SEXP, SEXP);
extern SEXP sirt_isop_tests_C(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_md_pattern_csource(SEXP);
extern SEXP sirt_mle_pcm_C(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_mle_pcm_group_C(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_mlnormal_proc_variance_shortcut_XY_restructure(SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_mlnormal_proc_variance_shortcut_Z_restructure(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_mlnormal_update_beta_rcpp_helper(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_mlnormal_update_V_rcpp_helper(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_MML2_CALCPOST_V1(SEXP, SEXP, SEXP);
extern SEXP sirt_MML2_CALCPOST_V2(SEXP, SEXP, SEXP);
extern SEXP sirt_MML2_CALCPOST_V3(SEXP, SEXP, SEXP);
extern SEXP sirt_MML2_RASCHTYPE_COUNTS(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_monoreg_rowwise_Cpp(SEXP, SEXP);
extern SEXP sirt_noharm_estFcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_noharm_estPcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_noharm_estPsicpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_parameters_jackknife(SEXP);
extern SEXP sirt_polychoric2_aux_rcpp(SEXP, SEXP, SEXP);
extern SEXP sirt_probs_pcm_groups_C(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_probs_pcm_nogroups_C(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_rm_arraymult1(SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_RM_CALCPOST(SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_rm_facets_calcprobs_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_rm_probraterfct1(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_rowCumsums2_source(SEXP);
extern SEXP sirt_rowKSmallest_C(SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_rowMaxsCPP_source(SEXP);
extern SEXP sirt_rowmins2_bundle_C(SEXP, SEXP);
extern SEXP sirt_SMIRT_CALCPOST(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_SMIRT_CALCPROB_COMP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_SMIRT_CALCPROB_NONCOMP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_SMIRT_CALCPROB_PARTCOMP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sirt_tetrachoric2_rcpp_aux(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"sirt_calc_copula_itemcluster_C",                      (DL_FUNC) &sirt_calc_copula_itemcluster_C,                       2},
    {"sirt_calccounts_pcm_groups_C",                        (DL_FUNC) &sirt_calccounts_pcm_groups_C,                         7},
    {"sirt_choppin_rowaveraging",                           (DL_FUNC) &sirt_choppin_rowaveraging,                            3},
    {"sirt_eigenvaluesDsirt",                               (DL_FUNC) &sirt_eigenvaluesDsirt,                                4},
    {"sirt_evm_aux",                                        (DL_FUNC) &sirt_evm_aux,                                         6},
    {"sirt_evm_comp_matrix_poly",                           (DL_FUNC) &sirt_evm_comp_matrix_poly,                            9},
    {"sirt_firsteigenvalsirt2",                             (DL_FUNC) &sirt_firsteigenvalsirt2,                              4},
    {"sirt_gooijer_csn_table",                              (DL_FUNC) &sirt_gooijer_csn_table,                               7},
    {"sirt_ia_optim_lambda",                                (DL_FUNC) &sirt_ia_optim_lambda,                                 8},
    {"sirt_ia_optim_nu",                                    (DL_FUNC) &sirt_ia_optim_nu,                                    11},
    {"sirt_interval_index_C",                               (DL_FUNC) &sirt_interval_index_C,                                2},
    {"sirt_isop_tests_C",                                   (DL_FUNC) &sirt_isop_tests_C,                                    5},
    {"sirt_md_pattern_csource",                             (DL_FUNC) &sirt_md_pattern_csource,                              1},
    {"sirt_mle_pcm_C",                                      (DL_FUNC) &sirt_mle_pcm_C,                                       8},
    {"sirt_mle_pcm_group_C",                                (DL_FUNC) &sirt_mle_pcm_group_C,                                 9},
    {"sirt_mlnormal_proc_variance_shortcut_XY_restructure", (DL_FUNC) &sirt_mlnormal_proc_variance_shortcut_XY_restructure,  4},
    {"sirt_mlnormal_proc_variance_shortcut_Z_restructure",  (DL_FUNC) &sirt_mlnormal_proc_variance_shortcut_Z_restructure,   8},
    {"sirt_mlnormal_update_beta_rcpp_helper",               (DL_FUNC) &sirt_mlnormal_update_beta_rcpp_helper,                7},
    {"sirt_mlnormal_update_V_rcpp_helper",                  (DL_FUNC) &sirt_mlnormal_update_V_rcpp_helper,                  11},
    {"sirt_MML2_CALCPOST_V1",                               (DL_FUNC) &sirt_MML2_CALCPOST_V1,                                3},
    {"sirt_MML2_CALCPOST_V2",                               (DL_FUNC) &sirt_MML2_CALCPOST_V2,                                3},
    {"sirt_MML2_CALCPOST_V3",                               (DL_FUNC) &sirt_MML2_CALCPOST_V3,                                3},
    {"sirt_MML2_RASCHTYPE_COUNTS",                          (DL_FUNC) &sirt_MML2_RASCHTYPE_COUNTS,                           6},
    {"sirt_monoreg_rowwise_Cpp",                            (DL_FUNC) &sirt_monoreg_rowwise_Cpp,                             2},
    {"sirt_noharm_estFcpp",                                 (DL_FUNC) &sirt_noharm_estFcpp,                                 16},
    {"sirt_noharm_estPcpp",                                 (DL_FUNC) &sirt_noharm_estPcpp,                                 16},
    {"sirt_noharm_estPsicpp",                               (DL_FUNC) &sirt_noharm_estPsicpp,                               16},
    {"sirt_parameters_jackknife",                           (DL_FUNC) &sirt_parameters_jackknife,                            1},
    {"sirt_polychoric2_aux_rcpp",                           (DL_FUNC) &sirt_polychoric2_aux_rcpp,                            3},
    {"sirt_probs_pcm_groups_C",                             (DL_FUNC) &sirt_probs_pcm_groups_C,                              6},
    {"sirt_probs_pcm_nogroups_C",                           (DL_FUNC) &sirt_probs_pcm_nogroups_C,                            5},
    {"sirt_rm_arraymult1",                                  (DL_FUNC) &sirt_rm_arraymult1,                                   4},
    {"sirt_RM_CALCPOST",                                    (DL_FUNC) &sirt_RM_CALCPOST,                                     4},
    {"sirt_rm_facets_calcprobs_cpp",                        (DL_FUNC) &sirt_rm_facets_calcprobs_cpp,                        12},
    {"sirt_rm_probraterfct1",                               (DL_FUNC) &sirt_rm_probraterfct1,                                5},
    {"sirt_rowCumsums2_source",                             (DL_FUNC) &sirt_rowCumsums2_source,                              1},
    {"sirt_rowKSmallest_C",                                 (DL_FUNC) &sirt_rowKSmallest_C,                                  4},
    {"sirt_rowMaxsCPP_source",                              (DL_FUNC) &sirt_rowMaxsCPP_source,                               1},
    {"sirt_rowmins2_bundle_C",                              (DL_FUNC) &sirt_rowmins2_bundle_C,                               2},
    {"sirt_SMIRT_CALCPOST",                                 (DL_FUNC) &sirt_SMIRT_CALCPOST,                                  6},
    {"sirt_SMIRT_CALCPROB_COMP",                            (DL_FUNC) &sirt_SMIRT_CALCPROB_COMP,                             6},
    {"sirt_SMIRT_CALCPROB_NONCOMP",                         (DL_FUNC) &sirt_SMIRT_CALCPROB_NONCOMP,                          6},
    {"sirt_SMIRT_CALCPROB_PARTCOMP",                        (DL_FUNC) &sirt_SMIRT_CALCPROB_PARTCOMP,                         7},
    {"sirt_tetrachoric2_rcpp_aux",                          (DL_FUNC) &sirt_tetrachoric2_rcpp_aux,                           3},
    {NULL, NULL, 0}
};

void R_init_sirt(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
