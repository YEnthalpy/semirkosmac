// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// km
List km(NumericVector e, NumericVector d, NumericVector p);
RcppExport SEXP _semirkosmac_km(SEXP eSEXP, SEXP dSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type e(eSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d(dSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(km(e, d, p));
    return rcpp_result_gen;
END_RCPP
}
// gehan_ns
arma::mat gehan_ns(arma::mat x, arma::vec y, arma::uvec d, arma::vec beta, arma::vec p, int n);
RcppExport SEXP _semirkosmac_gehan_ns(SEXP xSEXP, SEXP ySEXP, SEXP dSEXP, SEXP betaSEXP, SEXP pSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(gehan_ns(x, y, d, beta, p, n));
    return rcpp_result_gen;
END_RCPP
}
// gehan_mtg
arma::mat gehan_mtg(arma::mat x, arma::vec d, arma::vec e, arma::uvec ind);
RcppExport SEXP _semirkosmac_gehan_mtg(SEXP xSEXP, SEXP dSEXP, SEXP eSEXP, SEXP indSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type e(eSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ind(indSEXP);
    rcpp_result_gen = Rcpp::wrap(gehan_mtg(x, d, e, ind));
    return rcpp_result_gen;
END_RCPP
}
// gehan_smth
arma::mat gehan_smth(const arma::mat& x, const arma::vec& y, const arma::vec& d, const arma::vec& p, const arma::vec& b, const int& n);
RcppExport SEXP _semirkosmac_gehan_smth(SEXP xSEXP, SEXP ySEXP, SEXP dSEXP, SEXP pSEXP, SEXP bSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(gehan_smth(x, y, d, p, b, n));
    return rcpp_result_gen;
END_RCPP
}
// gehan_s_mtg
arma::mat gehan_s_mtg(const arma::mat& x, const arma::vec& y, const arma::vec& d, const arma::vec& p, const arma::vec& bt, const arma::uvec ind, const int& n);
RcppExport SEXP _semirkosmac_gehan_s_mtg(SEXP xSEXP, SEXP ySEXP, SEXP dSEXP, SEXP pSEXP, SEXP btSEXP, SEXP indSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type bt(btSEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type ind(indSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(gehan_s_mtg(x, y, d, p, bt, ind, n));
    return rcpp_result_gen;
END_RCPP
}
// gehan_s_jaco
arma::mat gehan_s_jaco(const arma::mat& x, const arma::vec& y, const arma::vec& d, const arma::vec& p, const arma::vec& b, const int& n);
RcppExport SEXP _semirkosmac_gehan_s_jaco(SEXP xSEXP, SEXP ySEXP, SEXP dSEXP, SEXP pSEXP, SEXP bSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(gehan_s_jaco(x, y, d, p, b, n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_semirkosmac_km", (DL_FUNC) &_semirkosmac_km, 3},
    {"_semirkosmac_gehan_ns", (DL_FUNC) &_semirkosmac_gehan_ns, 6},
    {"_semirkosmac_gehan_mtg", (DL_FUNC) &_semirkosmac_gehan_mtg, 4},
    {"_semirkosmac_gehan_smth", (DL_FUNC) &_semirkosmac_gehan_smth, 6},
    {"_semirkosmac_gehan_s_mtg", (DL_FUNC) &_semirkosmac_gehan_s_mtg, 7},
    {"_semirkosmac_gehan_s_jaco", (DL_FUNC) &_semirkosmac_gehan_s_jaco, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_semirkosmac(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
