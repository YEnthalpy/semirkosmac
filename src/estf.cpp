#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
List km(NumericVector e, NumericVector d,
        NumericVector p)
{
  int t = e.length();
  // int l = 0;
  // while (l < t)
  // {
  //   if (d[l] == 1)
  //   {
  //     break;
  //   }
  //   l++;
  // }
  NumericVector s(t, 1.0);
  // NumericVector eu(t);
  double den = 0;
  int i = t - 1;
  // int k = t - 1;
  while (i >= 0)
  {
    double ei = e[i];
    int j = i;
    double nut = 0;
    while (j >= 0)
    {
      if (e[j] < ei)
      {
        break;
      }
      else
      {
        double pj = p[j];
        den += 1 / pj;
        if (d[j] == 1)
        {
          nut += 1 / pj;
        }
      }
      j--;
    }
    // i = j;
    // if (nut == 0)
    // {
    //   continue;
    // }
    double tmp = 1 - nut / den;
    s[j + 1] = tmp;
    i = j;
    // s[k] = tmp;
    // eu[k] = e[i + 1];
    // k--;
  }
  // s.erase(0, k + 1);
  // eu.erase(0, k + 1);
  // NumericVector out(t - k - 1);
  NumericVector out(t);
  double start = 1;
  for (int i = 0; i < t; i++)
  {
    start = start * s[i];
    out[i] = start;
  }
  // if (l != 0)
  // {
  //   out.push_front(1);
  //   eu.push_front(e[l - 1]);
  // }
  List out1 = List::create(e, out);
  return out1;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat gehan_ns(arma::mat x, arma::vec y, arma::uvec d,
                   arma::vec beta, arma::vec p, int n)
{
  arma::vec e = y - x * beta;
  int r = y.n_elem;
  int m = x.n_cols;
  arma::mat out(r, m);
  for (int i = 0; i < r; i++)
  {
    if (d[i] == 0)
    {
      continue;
    }
    arma::rowvec nut(m);
    for (int j = 0; j < r; j++)
    {
      if (e[i] <= e[j])
      {
        nut += (x.row(i) - x.row(j)) / p[j];
      }
    }
    out.row(i) = nut / p[i] / n / r / r;
  }
  return out;
}



// Binary search find the first element that isn't smaller than x
// arry is sorted

int f_loc(arma::vec arry, double x, int low, int high)
{
  if (arry[0] >= x)
  {
    return 0;
  }
  int mid = 0;
  while (low <= high)
  {
    mid = (low + high) / 2;
    if (arry[mid] >= x)
    {
      if (arry[mid - 1] < x)
      {
        return mid;
      }
      high = mid - 1;
    }
    else
    {
      low = mid + 1;
    }
  }
  return mid;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat gehan_mtg(arma::mat x, arma::vec d, arma::vec e,
                    arma::uvec ind)
{
  int n = x.n_rows;
  int m = x.n_cols;
  int r0 = ind.n_elem;
  arma::mat out(n, m);
  // Define the pilot subsample
  arma::mat x_pt = x.rows(ind);
  arma::vec d_pt = d.elem(ind);
  arma::vec e_pt = e.elem(ind);
  // Calculate Zi and Ni in the pilot subsample
  arma::mat xsum(r0, m);
  arma::rowvec csum(r0);
  arma::rowvec xtmp(m);
  double ctmp = 0;
  int k = r0 - 1;
  while (k >= 0)
  {
    ctmp += 1;
    xtmp += x_pt.row(k);
    xsum.row(k) = xtmp;
    csum[k] = ctmp;
    int l = k + 1;
    // Dealing with ties
    while (e_pt[l] == e_pt[k])
    {
      xsum.row(l) = xsum.row(k);
      csum[l] = csum[k];
      l++;
    }
    k--;
  }
  // Calculate the contribution of the ith observation
  for (int i = 0; i < n; i++)
  {
    // Find the location of e[i] in e_pt (O{log(r_0)})
    // By binary search
    int loc = f_loc(e_pt, e[i], 0, r0 - 1);
    arma::rowvec cont_i(m);
    arma::rowvec xi = x.row(i);
    if (d[i] == 1)
    {
      cont_i = csum[loc] * xi - xsum.row(loc);
    }
    // Dealing with ties
    int nloc = loc;
    while (e[i] == e_pt[nloc])
    {
      nloc++;
    }
    if (e[i] > e_pt[nloc])
    {
      nloc++;
    }
    nloc--;
    for (int j = 0; j <= nloc; j++)
    {
      if (d_pt[j] == 1)
      {
        cont_i += xsum.row(j) / csum[j] - xi;
      }
    }
    out.row(i) = cont_i / r0;
  }
  return out;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat gehan_smth(const arma::mat &x, const arma::vec &y,
                     const arma::vec &d, const arma::vec &p,
                     const arma::vec &b, const int &n)
{
  int r = y.n_elem;
  int m = b.n_elem;
  arma::vec e = y - x * b;
  arma::mat out(r, m);
  for (int i = 0; i < r - 1; i++)
  {
    arma::mat xj = x.rows(i + 1, r - 1);
    arma::mat xdif = x.row(i) - xj.each_row();
    arma::vec edif = e.subvec(i + 1, r - 1) - e[i];
    arma::vec dij = d.subvec(i + 1, r - 1);
    arma::vec rij = sqrt(sum(xdif % xdif, 1));
    arma::vec pij = p.subvec(i + 1, r - 1);
    arma::vec phij = arma::normcdf(sqrt(r) * edif / rij);
    phij.replace(arma::datum::nan, 1);
    out.row(i) = sum(xdif.each_col() % ((phij % (d[i] + dij) - dij) / pij) / p[i], 0) / n / n / r / r;
  }
  return out;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat gehan_s_mtg(const arma::mat &x, const arma::vec &y, 
                      const arma::vec &d, const arma::vec &p,
                      const arma::vec &bt, const arma::uvec ind,
                      const int &n)
{
  int r = x.n_rows;
  int m = x.n_cols;
  arma::vec e = y - x * bt;
  int r0 = ind.n_elem;
  // Define the pilot subsample
  arma::mat x_pt = x.rows(ind);
  arma::vec d_pt = d.elem(ind);
  arma::vec e_pt = e.elem(ind);
  arma::vec p_pt = p.elem(ind);
  // Calculate Zi and Ni in the pilot subsample
  arma::mat xpsum(r, m);
  arma::vec psum(r);
  for (int i = 0; i < r; i++)
  {
    arma::mat xdif = x_pt.each_row() - x.row(i);
    arma::vec edif = e_pt - e[i];
    arma::vec rij = sqrt(arma::sum(xdif % xdif, 1));
    arma::vec phij = arma::normcdf(sqrt(r) * edif / rij);
    phij.replace(arma::datum::nan, 0);
    psum[i] = arma::accu(phij / p_pt) / r0;
    xpsum.row(i) = arma::sum(x_pt.each_col() % (phij / p_pt), 0) / r0;
  }
  arma::mat xpsum_pt = xpsum.rows(ind);
  arma::vec psum_pt = psum.elem(ind);
  arma::mat tmp = xpsum_pt.each_col() / psum_pt;
  // Calculate the estimating function
  arma::mat out(r, m);
  arma::mat out1 = x.each_col() % psum - xpsum;
  arma::mat out2(r, m);
  for (int j = 0; j < r0; j++)
  {
    if (d_pt[j] == 0 || std::isnan(arma::accu(tmp.row(j))))
    {
      continue;
    }
    arma::mat xdif = x.each_row() - x_pt.row(j);
    arma::vec edif = e - e_pt[j];
    arma::vec rij = sqrt(arma::sum(xdif % xdif, 1));
    arma::vec phij = arma::normcdf(sqrt(r) * edif / rij);
    phij.replace(arma::datum::nan, 0);
    out2 += (x.each_row() - tmp.row(j)).each_col() % phij / p_pt[j] / r0;
  }
  out = (out1.each_col() % (d / p)  - out2.each_col() / p) / r0 / n / n;
  return out;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat gehan_s_jaco(const arma::mat &x, const arma::vec &y,
                       const arma::vec &d, const arma::vec &p,
                       const arma::vec &b, const int &n)
{
  int r = y.n_elem;
  int m = b.n_elem;
  arma::vec e = y - x * b;
  arma::mat out(m, m);
  for (int i = 0; i < r - 1; i++)
  {
    arma::mat xj = x.rows(i + 1, r - 1);
    arma::mat xdif = x.row(i) - xj.each_row();
    arma::vec edif = e.subvec(i + 1, r - 1) - e[i];
    arma::vec dij = d.subvec(i + 1, r - 1) + d[i];
    arma::vec rij = sqrt(sum(xdif % xdif, 1));
    arma::vec pij = p.subvec(i + 1, r - 1);
    arma::vec phij = arma::normpdf(sqrt(r) * edif / rij);
    phij.replace(arma::datum::nan, 1);
    arma::mat xdif2 = xdif.each_col() % (phij % dij / rij / pij);
    xdif2.replace(arma::datum::nan, 0);
    out += xdif.t() * xdif2 / p[i] / n / n / r / r * sqrt(r);
  }
  return out;
}


// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// arma::mat logrank_ns(arma::mat x, arma::vec y, arma::uvec d,
//                      arma::vec beta, arma::vec p)
// {
//   arma::vec e = y - x * beta;
//   int r = y.n_elem;
//   int m = x.n_cols;
//   arma::mat out(r, m);
//   for (int i = 0; i < r; i++)
//   {
//     if (d[i] == 0)
//     {
//       continue;
//     }
//     arma::rowvec nut(m);
//     double den = 0;
//     for (int j = 0; j < r; j++)
//     {
//       double pj = p[j];
//       if (e[i] <= e[j])
//       {
//         nut += x.row(j) / pj;
//         den += 1 / pj;
//       }
//     }
//     out.row(i) = (x.row(i) - nut / den) / p[i] / r;
//   }
//   return out;
// }
// 
// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// arma::mat logrank_smth(arma::mat x, arma::vec y, arma::uvec d,
//                        arma::vec beta, arma::vec p)
// {
//   arma::vec e = y - x * beta;
//   int r = y.n_elem;
//   int m = x.n_cols;
//   arma::mat out(r, m);
//   for (int i = 0; i < r; i++)
//   {
//     if (d[i] == 0)
//     {
//       continue;
//     }
//     arma::rowvec nut(m);
//     double den = 0;
//     for (int j = 0; j < r; j++)
//     {
//       double rij = arma::norm(x.row(i) - x.row(j));
//       double phij = arma::normcdf((e[j] - e[i]) / rij) / p[j];
//       nut += x.row(j) * phij;
//       den += phij;
//     }
//     out.row(i) = (x.row(i) - nut / den) / p[i] / r;
//   }
//   return out;
// }

// arma::vec groyden(const arma::mat &x, const arma::vec &y,
//                   const arma::vec &d, const arma::vec &p,
//                   const arma::vec &b0, const int &n,
//                   const int &maxit, const double &tol)
// {
//   int m = b0.n_elem;
//   arma::mat jki = arma::inv(gehan_s_jaco(x, y, d, p, b0, n));
//   arma::vec fk = gehan_smth(x, y, d, p, b0, n);
//   arma::vec dk1 = -jki * fk;
//   arma::vec bk1 = b0 + dk1;
//   arma::vec out(m);
//   int k = 0;
//   while (k < maxit)
//   {
//     arma::vec fk1 = gehan_smth(x, y, d, p, bk1, n);
//     arma::vec dfk1 = fk1 - fk;
//     arma::vec tmp = sum(jki.each_row() % dfk1.t(), 1);
//     arma::mat jk1i = jki + (dk1 - tmp) * (sum(jki.each_col() % dk1, 0)) / sum(tmp % dk1);
//     arma::vec dk_new = -jk1i * fk1;
//     arma::vec bk_new = bk1 + dk1;
//     double ero = arma::norm(dk_new);
//     if (ero < tol)
//     {
//       out = dk_new;
//       break;
//     }
//     arma::vec dk1 = dk_new;
//     arma::vec bk1 = bk_new;
//     arma::vec fk = fk1;
//     arma::mat jki = jk1i;
//     k++;
//   }
//   return (out);
// }
// 
// double up_alpf(double &alph, double &f, double &fnew,
//                double &tau_min, double &tau_max)
// {
//   double alpht = alph * alph * f / (fnew + f * (2 * alph - 1));
//   double out = 0;
//   if (alpht < tau_min * alph)
//     out = tau_min * alph;
//   else if (alpht > tau_max * alph)
//     out = tau_max * alph;
//   else
//     out = alpht;
//   return out;
// }
// 
// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// arma::vec bb(const arma::mat &x, const arma::vec &y,
//              const arma::vec &d, const arma::vec &p,
//              const arma::vec &b0, const int &n,
//              const int &maxit, const double &tol)
// {
//   // Set constant values
//   double sig_min = 1e-10, sig_max = 1e10;
//   double tau_min = 0.1, tau_max = 0.5, gam = 1e-4;
//   int M = 10;
//   int m = b0.n_elem;
//   // First iteration
//   arma::vec bk = b0;
//   arma::vec Fk = gehan_smth(x, y, d, p, bk, n);
//   double nF0 = arma::norm(Fk);
//   double fk = nF0 * nF0;
//   double sig_k = 1.0;
//   // Initiate some vectors
//   arma::vec Fbark(M, arma::fill::zeros);
//   arma::vec sk(m, arma::fill::zeros);
//   arma::vec yk(m, arma::fill::zeros);
//   arma::vec out(m);
//   // Start iteration
//   int k = 0;
//   while (k <= maxit)
//   {
//     double alphap = 1.0, alpham = 1.0;
//     if (k > 0)
//       sig_k = arma::accu(sk % sk) / arma::accu(sk % yk);
//     if (abs(sig_k) > sig_max || abs(sig_k) < sig_min)
//     {
//       if (fk > 1)
//         sig_k = 1;
//       else if (fk <= 1 & fk >= 1e-10)
//         sig_k = 1 / sqrt(fk);
//       else if (fk < 1e-10)
//         sig_k = 1e5;
//     }
//     if (k == 0 && (nF0 > 1))
//       sig_k = 1 / nF0;
//     // Update Fbk and get new fmk
//     Fbark.shed_row(M - 1);
//     Fbark.insert_rows(0, 1);
//     Fbark[0] = fk;
//     double fmk = Fbark.max();
//     // Update dk
//     arma::vec dk = -sig_k * Fk;
//     // line searching method
//     double etak = nF0 / (1 + k) / (1 + k);
//     while (1)
//     {
//       arma::vec bkp = bk + alphap * dk;
//       arma::vec Fkp = gehan_smth(x, y, d, p, bkp, n);
//       double fkp = arma::accu(Fkp % Fkp);
//       if (fkp <= fmk + etak - gam * alphap * alphap * fk)
//       {
//         yk = Fkp - Fk;
//         Fk = Fkp;
//         sk = bkp - bk;
//         bk = bkp;
//         break;
//       }
//       arma::vec bkm = bk - alpham * dk;
//       arma::vec Fkm = gehan_smth(x, y, d, p, bkm, n);
//       double fkm = arma::accu(Fkm % Fkm);
//       if (fkm <= fmk + etak - gam * alpham * alpham * fk)
//       {
//         yk = Fkm - Fk;
//         Fk = Fkm;
//         sk = bkm - bk;
//         bk = bkm;
//         break;
//       }
//       // Update alphap
//       alphap = up_alpf(alphap, fk, fkp, tau_min, tau_max);
//       // Update alpham
//       alpham = up_alpf(alpham, fk, fkm, tau_min, tau_max);
//     }
//     fk = arma::accu(Fk % Fk);
//     if (sqrt(fk) <= tol)
//     {
//       out = bk;
//       break;
//     }
//     else
//       k++;
//   }
//   if (k < maxit)
//   {
//     out.insert_rows(0, 2);
//     out[0] = 0;
//     out[1] = k;
//   }
//   if (k == maxit)
//   {
//     out.insert_rows(0, 2);
//     out[0] = 1;
//     out[1] = k;
//   }
//   return (out);
// }
