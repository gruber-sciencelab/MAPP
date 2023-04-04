/*
##############################################################################
#
#   Optimize MARA model parameters (find maximum likelihood estimates)
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 20-01-2020
#   LICENSE: GPL_3.0
#   LICENSE: https://www.gnu.org/licenses/gpl-3.0.html
#
##############################################################################
*/

#include <Rcpp.h>

#include "optim.c"
#include "zeroin.c"
using namespace Rcpp;

// Container for all shared data
class Box {
 public:
  NumericMatrix ft, t_prefactor, N, NAb, Nft, Nt;
  NumericVector A, b, c, sumS_ft;
  int e, m;
  Box(NumericMatrix ft_, NumericMatrix t_prefactor_, NumericMatrix N_,
      NumericMatrix NAb_, NumericVector A_, NumericVector b_, NumericVector c_,
      NumericVector sumS_ft_, NumericMatrix Nft_, NumericMatrix Nt_, int e_,
      int m_)
      : ft(ft_),
        t_prefactor(t_prefactor_),
        N(N_),
        A(A_),
        b(b_),
        c(c_),
        sumS_ft(sumS_ft_),
        e(e_),
        m(m_),
        NAb(NAb_),
        Nft(Nft_),
        Nt(Nt_) {}
};

// function to minimize in BFGS search
double negative_LL(int n, double *params, void *ex) {
  // cast the void pointer to the Box pointer for de-reference
  Box *ptr = static_cast<Box *>(ex);
  NumericMatrix ft = (*ptr).ft;
  NumericMatrix t_prefactor = (*ptr).t_prefactor;
  NumericMatrix N = (*ptr).N;
  NumericVector c = (*ptr).c;
  int m = (*ptr).m;

  int n_samples = ft.ncol();
  int n_exons = ft.nrow();
  double term1 = 0;
  double term2 = 0;
  double bcNA;
  for (int s = 0; s < n_samples; s++) {
    for (int e = 0; e < n_exons; e++) {
      bcNA = params[s + n_samples] + c[e] + N(e, m) * params[s];
      term1 += ft(e, s) * bcNA;
      term2 += t_prefactor(e, s) * log(1 + exp(bcNA));
    }
  }
  double negative_LL = -1 * (term1 - term2);
  return negative_LL;
}

// function fo minimize in BFGS search
void grad_negative_LL(int n, double *params, double *g, void *ex) {
  // cast the void pointer to the Box pointer for de-reference
  Box *ptr = static_cast<Box *>(ex);
  NumericMatrix ft = (*ptr).ft;
  NumericMatrix t_prefactor = (*ptr).t_prefactor;
  NumericMatrix N = (*ptr).N;
  NumericMatrix Nft = (*ptr).Nft;
  NumericMatrix Nt = (*ptr).Nt;
  NumericVector c = (*ptr).c;
  int m = (*ptr).m;

  int n_samples = ft.ncol();
  int n_exons = ft.nrow();
  double temp;
  double e_bcNA;

  for (int s = 0; s < n_samples; s++) {
    // d(-LL)/dA_i
    temp = 0.0;
    for (int e = 0; e < n_exons; e++) {
      e_bcNA = exp(params[s + n_samples] + c[e] + N(e, m) * params[s]);
      temp += Nft(e, s) - (Nt(e, s) * e_bcNA / (1 + e_bcNA));
    }
    g[s] = -1 * temp;
    // d(-LL)/b_i
    temp = 0.0;
    for (int e = 0; e < n_exons; e++) {
      e_bcNA = exp(params[s + n_samples] + c[e] + N(e, m) * params[s]);
      temp += ft(e, s) - (t_prefactor(e, s) * e_bcNA / (1 + e_bcNA));
    }
    g[s + n_samples] = -1 * temp;
  }
}

// function that checks the initial endpoints in bisection method
int sgn(double x) {
  if (x > 0.0) {
    return 1;
  } else if (x < 0.0) {
    return -1;
  } else {
    return 0;
  }
}

// bisection optimization: simplified version: N*A+b is one constant
double dLL_dc(double x, void *info) {
  // cast the void pointer to the Box pointer for de-reference
  Box *ptr = static_cast<Box *>(info);
  const int n_samples = (*ptr).ft.ncol();
  int e = (*ptr).e;
  int m = (*ptr).m;
  double e_bcNA;
  double value = (*ptr).sumS_ft[e];
  for (int s = 0; s < n_samples; s++) {
    e_bcNA = exp(x + (*ptr).NAb(e, s));
    value -= (*ptr).t_prefactor(e, s) * e_bcNA / (1 + e_bcNA);
  }
  return value;
}

// [[Rcpp::export]]
List fit_model_parameters(NumericMatrix i, NumericMatrix t, NumericMatrix N,
                          NumericVector t_crit_arg, bool average_expressions) {
  // prepare data structures for the parameters and results, also other
  // variables NumericVector and NumericMatrix is always initialised with zero
  const int n_samples = i.ncol();
  const int n_motifs = N.ncol();
  const int n_exons = i.nrow();
  NumericMatrix A_stdA_b_stdb(n_motifs, 4 * n_samples + 1);
  NumericMatrix NAb(n_exons, n_samples);
  NumericVector A(n_samples);
  NumericVector b(n_samples);
  NumericVector c(n_exons);
  NumericMatrix f(n_exons, n_samples);
  NumericMatrix t_prefactor(n_exons, n_samples);
  NumericMatrix ft(n_exons, n_samples);
  double mean_A;
  double mean_b;
  double x_upper;
  double x_lower;
  double fx_upper;
  double fx_lower;
  int uniroot_max_iter;
  double uniroot_eps;

  // calculate the inclusion fractions (PSI scores)
  for (int e = 0; e < n_exons; e++) {
    for (int s = 0; s < n_samples; s++) {
      f(e, s) = i(e, s) / t(e, s);
    }
  }

  // average TPM expressions over samples
  if (average_expressions){

    double t_crit = t_crit_arg(0);

    // calculate the shared t prefactor (including t_crit)
    double temp_for_row_avg;
    for (int e = 0; e < n_exons; e++) {
      temp_for_row_avg = mean(t(e, _));
      for (int s = 0; s < n_samples; s++) {
        t_prefactor(e, s) = temp_for_row_avg / (temp_for_row_avg + t_crit);
      }
    }

    // calculate the C,X,R and R*t_prefactor
    // the latter variables with t_prefactor have in fact t_prefactor*R (!)
    double C = sum(t);
    double X = sum(t_prefactor);
    double R = C / X;
    for (int e = 0; e < n_exons; e++) {
      for (int s = 0; s < n_samples; s++) {
        t_prefactor(e, s) = t_prefactor(e, s) * R;
      }
    }

  // operate on original per-sample TPM values
  } else {

    // calculate the shared t prefactor (including t_crit)
    for (int e = 0; e < n_exons; e++) {
      for (int s = 0; s < n_samples; s++) {
        t_prefactor(e, s) = t(e, s) / (t(e, s) + t_crit_arg(s));
      }
    }

    // calculate the C,X,R and R*t_prefactor
    // the latter variables with t_prefactor have in fact t_prefactor*R (!)
    NumericVector C = colSums(t);
    NumericVector X = colSums(t_prefactor);
    NumericVector R = C / X;
    for (int e = 0; e < n_exons; e++) {
      for (int s = 0; s < n_samples; s++) {
        t_prefactor(e, s) = t_prefactor(e, s) * R(s);
      }
    }
  }

  // calculate the shared term: f*t_prefactor (including R correction)
  for (int e = 0; e < n_exons; e++) {
    for (int s = 0; s < n_samples; s++) {
      ft(e, s) = f(e, s) * t_prefactor(e, s);
    }
  }

  // LL-related
  double LL_current;
  double LL_old;
  const double LL_eps = 0.01;
  const double LL_conv_test = -log(1 + LL_eps);

  // arguments for the vmmin C function
  int A_b_size = n_samples * 2;
  double *A_b = vect(A_b_size);
  for (int xx = 0; xx < A_b_size; xx++) A_b[xx] = 0;
  double zero_zero = 0.0;
  double *Fmin = &zero_zero;
  int maxit = 10000;
  int trace = 0;
  int *mask;
  mask = (int *)R_alloc(A_b_size, sizeof(int));
  for (int xx = 0; xx < A_b_size; xx++) mask[xx] = 1;
  double abstol = 1e-15;  // this one is not relevant in our case
  double reltol = 1e-15;
  int nREPORT = 10;
  int zero = 0;
  int *fncount = &zero;
  int *grcount = &zero;
  int *fail = &zero;

  // precomputed term for all motifs for all exons:
  NumericVector sumS_ft = rowSums(ft);

  // precomputed terms for a given motif
  NumericMatrix Nft(n_exons, n_samples);
  NumericMatrix Nt(n_exons, n_samples);
  NumericMatrix Hessian_shared_terms(n_exons, 4 * n_samples);
  NumericMatrix d2LL_db_s_2_matrix(n_exons, n_samples);
  NumericMatrix d2LL_dA_s_2_matrix(n_exons, n_samples);
  NumericMatrix d2LL_dA_s_db_s_matrix(n_exons, n_samples);
  NumericVector H_AA(n_samples);
  NumericVector H_Ab(n_samples);
  NumericVector H_bb(n_samples);

  // given A,b=0 MLE for c are given by:
  NumericVector c_init = rowSums(ft) / (rowSums(t_prefactor) - rowSums(ft));
  for (int e = 0; e < n_exons; e++) {
    c_init[e] = log(c_init[e]);
  }
  // treat it as const

  // create the container with data for BFGS optimization
  Box box(ft, t_prefactor, N, NAb, A, b, c, sumS_ft, Nft, Nt, -1, -1);
  void *void_box = &box;
  // from now on manipulate them only via the box object!

  // parameters for the std_A inference
  List Richardson_params = List::create(
      Named("eps") = 1e-4, _["zero.tol"] = sqrt(2.220446e-16 / 7e-7),
      _["r"] = 2, _["v"] = 2, _["d"] = 1);

  // for every motif:
  for (int m = 0; m < n_motifs; m++) {
    // initialize the precomputed terms for a given motif
    for (int e = 0; e < n_exons; e++) {
      for (int s = 0; s < n_samples; s++) {
        Nft(e, s) = N(e, m) * ft(e, s);
        Nt(e, s) = N(e, m) * t_prefactor(e, s);
      }
    }

    // A,b,c get overwritten after every motif so re-assign the initial values:
    // for all samples initial A,b = 0,
    // for all exons initial c = MLE given A,b = 0 (as above)
    for (int s = 0; s < n_samples; s++) {
      box.A[s] = 0;
      box.b[s] = 0;
    }
    for (int e = 0; e < n_exons; e++) {
      box.c[e] = c_init[e];
    }
    box.m = m;

    // Given A,b=0 and MLE for c calculate the initial negative LL
    double term1 = 0;
    double term2 = 0;
    for (int s = 0; s < n_samples; s++) {
      for (int e = 0; e < n_exons; e++) {
        term1 += ft(e, s) * c[e];
        term2 += t_prefactor(e, s) * log(1 + exp(c[e]));
      }
    }
    LL_current = -1 * (term1 - term2);

    // Expectation-Maximization:
    do {
      // remember the old likelihood and now move: i->i+1
      LL_old = LL_current;

      // create a vector of A+b parameters for BFGS optimization
      for (int xx = 0; xx < n_samples; ++xx) {
        A_b[xx] = box.A[xx];
        A_b[xx + n_samples] = box.b[xx];
      }

      // Calculate MLE for A and b for all samples, c is treated as input here
      // call the BFGS implemented in C
      Fmin = &zero_zero;
      vmmin(A_b_size, A_b, Fmin, *negative_LL, *grad_negative_LL, maxit, trace,
            mask, abstol, reltol, nREPORT, void_box, fncount, grcount, fail);

      // get the MLE back to the vectors
      for (int xx = 0; xx < n_samples; ++xx) {
        box.A[xx] = A_b[xx];
        box.b[xx] = A_b[xx + n_samples];
      }

      // Neutral move that realizes the constraints on A and b; c will adjust to
      // them: The following four lines center A, b around zero - constraints
      mean_A = mean(box.A);
      mean_b = mean(box.b);
      box.A = box.A - mean_A;
      box.b = box.b - mean_b;
      // The following line compensates the shift and brings LL back to maximum
      box.c = box.c + mean_b + N(_, m) * mean_A;

      // given MLE for A,b: optimize for every c separately: find root of
      // dLL/dc_e.

      // pre-calculate a term shared for all c just once: NAb = N*A+b
      for (int e = 0; e < n_exons; e++) {
        for (int s = 0; s < n_samples; s++) {
          box.NAb(e, s) = box.N(e, m) * box.A[s] + box.b[s];
        }
      }

      for (int e = 0; e < n_exons; e++) {
        // define the endpoints for bisective root finding based on the previous
        // MLE for every c
        x_upper = box.c[e] + 10;
        x_lower = box.c[e] - 10;

        // calculate the dLL/dc_e at the endpoints
        fx_upper = box.sumS_ft[e];
        for (int s = 0; s < n_samples; s++) {
          fx_upper -= (box.t_prefactor(e, s) * exp(x_upper + NAb(e, s)) /
                       (1 + exp(x_upper + NAb(e, s))));
        }
        fx_lower = box.sumS_ft[e];
        for (int s = 0; s < n_samples; s++) {
          fx_lower -= (box.t_prefactor(e, s) * exp(x_lower + NAb(e, s)) /
                       (1 + exp(x_lower + NAb(e, s))));
        }

        // stop if the signs at the endpoints are the same
        if (sgn(fx_upper) == sgn(fx_lower))
          Rcpp::stop("incorrect range for bisection");

        // UNIROOT:
        box.e = e;
        uniroot_eps = pow(10, -30);
        uniroot_max_iter = 10000;
        box.c[e] = R_zeroin2(x_lower, x_upper, fx_lower, fx_upper, *dLL_dc,
                             void_box, &uniroot_eps, &uniroot_max_iter);
      }

      // calculate the new -LL given the updated parameters
      LL_current = negative_LL(A_b_size, A_b, void_box);

      // iterave until the LL converge i.e. while L[n] / L[n-1] >= 1 + eps
    } while (LL_current - LL_old < LL_conv_test);

    // ====================== ANALYTICAL STANDARD DEVIATIONS AROUND A and B

    // pre-calculate shared terms for the 2nd order partial derrivatives of LL:
    // 4 blocks (e,s):
    //   x = c_e + b_s + N_e * A_s
    //   Prob. inclusion: Theta = e^x / [1+e^x]
    //   1 - Theta
    //   Theta * (1-Theta)
    for (int e = 0; e < n_exons; e++) {
      for (int s = 0; s < n_samples; s++) {
        Hessian_shared_terms(e, s) = box.c[e] + box.b[s] + box.A[s] * N(e, m);
        Hessian_shared_terms(e, n_samples + s) =
            exp(Hessian_shared_terms(e, s)) /
            (1 + exp(Hessian_shared_terms(e, s)));
        Hessian_shared_terms(e, 2 * n_samples + s) =
            1 - Hessian_shared_terms(e, n_samples + s);
        Hessian_shared_terms(e, 3 * n_samples + s) =
            Hessian_shared_terms(e, 2 * n_samples + s) *
            Hessian_shared_terms(e, n_samples + s);
      }
    }

    // pre-calculate the 2nd order partial derrivatives of LL:
    for (int e = 0; e < n_exons; e++) {
      for (int s = 0; s < n_samples; s++) {
        d2LL_db_s_2_matrix(e, s) =
            Hessian_shared_terms(e, 3 * n_samples + s) * t_prefactor(e, s);
        d2LL_dA_s_db_s_matrix(e, s) = d2LL_db_s_2_matrix(e, s) * N(e, m);
        d2LL_dA_s_2_matrix(e, s) = d2LL_dA_s_db_s_matrix(e, s) * N(e, m);
      }
    }
    H_AA = colSums(d2LL_dA_s_2_matrix);
    H_Ab = colSums(d2LL_dA_s_db_s_matrix);
    H_bb = colSums(d2LL_db_s_2_matrix);

    // move the A,b and calculate their standard deviations
    for (int s = 0; s < n_samples; s++) {
      A_stdA_b_stdb(m, s) = box.A[s];
      A_stdA_b_stdb(m, s + n_samples) =
          sqrt(H_bb[s] / (H_bb[s] * H_AA[s] - pow(H_Ab[s], 2)));
      A_stdA_b_stdb(m, s + 2 * n_samples) = box.b[s];
      A_stdA_b_stdb(m, s + 3 * n_samples) =
          sqrt(H_AA[s] / (H_bb[s] * H_AA[s] - pow(H_Ab[s], 2)));
    }

    // save the LogLikelihood score as well
    A_stdA_b_stdb(m, 4 * n_samples) = -1 * LL_current;
  }

  // return a List:
  // NumericMatrix with per-sample A_s, std_A_s, b_s, std_b_s, LL
  // NumericVector with c_e
  List return_list(2);
  return_list[0] = A_stdA_b_stdb;
  return_list[1] = box.c;
  return return_list;
}
