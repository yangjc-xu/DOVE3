#include <iostream>
#include <armadillo>
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
struct CoxPHInformation{
  double LogLikelihood;
  arma::vec ScoreFun;
  arma::mat InformationMat;
};

arma::rowvec BS(double t, arma::vec knots, bool constantVE){
  int npc;
  double temp;
  arma::rowvec B;

  if (constantVE) {
    npc = knots.n_elem;
    temp = t>knots(npc-1) ? t-knots(npc-1) : 0;
  } else {
    npc = knots.n_elem + 1;
    temp = 0;
  }

  B.set_size(npc);
  B(0) = t-temp;
  for (int i = 1; i < npc; ++i) {
    B(i) = (t>knots(i-1) ? t-knots(i-1) : 0)-temp;
  }
  //B = B*0.0329;
  return B;
}

arma::rowvec BS_infection(double t, arma::vec knots, bool constantVE){
  int npc;
  double temp;
  arma::rowvec B;

  if (constantVE) {
    npc = knots.n_elem;
    temp = t>knots(npc-1) ? t-knots(npc-1) : 0;
  } else {
    npc = knots.n_elem + 1;
    temp = 0;
  }

  B.set_size(npc);
  B(0) = 1;
  for (int i = 1; i < npc; ++i) {
    B(i) = (t>knots(i-1) ? t-knots(i-1) : 0)-temp;
  }
  //B = B*0.0329;
  return B;
}

arma::rowvec BS2(double Time, arma::vec V, arma::vec S, double FirstInf,
                 arma::field<arma::vec> knots, arma::vec dimension, bool constantVE,
                 int infection_ind, int interact_ind, bool interact){
  int npc = arma::sum(dimension);
  arma::vec loc = arma::cumsum(dimension);
  arma::rowvec B(npc); B.zeros();
  int n_v = V.n_elem;

  if(interact){
    for(int j = 0; j < n_v; ++j){
      int start = (V(j)-1 == 0) ? 0 : loc(V(j)-2);
      int end = loc(V(j)-1)-1;
      arma::vec knots_j = knots(V(j)-1);
      if(V(j) < interact_ind){
        if(Time > S(j) & S(j) <= FirstInf){ B.subvec(start, end) = BS(Time - S(j), knots_j, constantVE);}
      } else if(V(j) >= interact_ind & V(j) < infection_ind){
        if(Time > S(j) & S(j) > FirstInf){ B.subvec(start, end) = BS(Time - S(j), knots_j, constantVE);}
      } else{
        if(Time > S(j)){ B.subvec(start, end) += BS_infection(Time - S(j), knots_j, constantVE);}
      }
    }
  } else{
    for(int j = 0; j < n_v; ++j){
      int start = (V(j)-1 == 0) ? 0 : loc(V(j)-2);
      int end = loc(V(j)-1)-1;
      arma::vec knots_j = knots(V(j)-1);
      if(V(j) < infection_ind){
        if(Time > S(j)){ B.subvec(start, end) = BS(Time - S(j), knots_j, constantVE);}
      } else{
        if(Time > S(j)){ B.subvec(start, end) += BS_infection(Time - S(j), knots_j, constantVE);}
      }
    }
  }
  return B;
}

struct CoxPHInformation get_Info(arma::vec beta, arma::vec gamma, arma::vec d,
                                 arma::vec Time, arma::vec t, arma::vec Delta, arma::vec C, arma::field<arma::vec> V,
                                 arma::field<arma::vec> S, arma::mat X, arma::vec m,
                                 arma::vec VacTime, arma::vec FirstInfTime, arma::vec dimension, int infection_ind, int interact_ind,
                                 arma::mat Z, arma::vec Score_constant, arma::field<arma::vec> knots,
                                 bool constantVE, bool interact, double cutoff){
  struct CoxPHInformation res;
  int p1 = beta.n_elem;
  int p2 = gamma.n_elem;
  int K = t.n_elem;
  double LogLikelihood = 0;
  arma::vec Score(p1+p2); Score.zeros();
  arma::mat InfoMat(p1+p2,p1+p2); InfoMat.zeros();

  // compute Xbeta and Zgamma
  arma::vec BetaX = X * beta;
  arma::vec GammaZ = Z * gamma;

  // This k loop is for time points. K = 470.
  arma::vec S0(K);S0.zeros();
  arma::mat S1(p1+p2,K); S1.zeros();
  arma::cube S2(p1+p2,p1+p2,K); S2.zeros();
  for(int k = 0; k < K; ++k){
    double tk = t(k);
    double temp0 = 0;
    arma::vec temp1(p1+p2); temp1.zeros();
    arma::mat temp2(p1+p2,p1+p2); temp2.zeros();

    // compute log-likelihood, score and information matrix
    arma::uvec index = arma::find(C >= tk);
    // This j loop is for subjects. The index.n_elem is several millions.
    for(int j = 0; j < index.n_elem; ++j){
      arma::vec V_index_j = V(index(j));
      arma::vec S_index_j = S(index(j));
      double FirstInf_index_j = FirstInfTime(index(j));
      arma::rowvec Zk = BS2(tk, V_index_j, S_index_j, FirstInf_index_j, knots, dimension, constantVE, infection_ind, interact_ind, interact);
      arma::rowvec XZk = join_horiz(X.row(index(j)), Zk);
      double temp = exp(BetaX(index(j)) + arma::as_scalar(Zk*gamma));
      // remove subjects who get infections with 30 days out of risk set
      arma::uvec index_inf = arma::find(V_index_j >= infection_ind);
      if(index_inf.n_elem > 0){
        for(int q = 0; q < index_inf.n_elem; ++q){
          if((tk - S_index_j(index_inf(q))) <= cutoff & (tk - S_index_j(index_inf(q))) > 0){
            temp = 0;
          }
        }
      }
      if(tk <= VacTime(index(j))){
        temp = 0;
      }
      temp0 += temp;
      temp1 += temp * XZk.t();
      temp2 += temp * XZk.t() * XZk;
    }

    S0(k) = temp0;
    S1.col(k) = temp1;
    S2.slice(k) = temp2;
  }

  for(int k = 0; k < K; ++k){
    LogLikelihood -= d(k) * log(S0(k));
    Score -= d(k) * (S1.col(k) / S0(k));
    arma::vec S3 = S1.col(k) / S0(k);
    InfoMat += d(k) * (S2.slice(k) / S0(k) - S3 * S3.t());
  }

  Score += Score_constant;

  res.LogLikelihood = LogLikelihood;
  res.ScoreFun = Score;
  res.InformationMat = InfoMat;

  return res;
}


// [[Rcpp::export]]
Rcpp::List Cox_general(arma::vec Time, arma::vec t, arma::vec Delta, arma::vec C, arma::field<arma::vec> V,
                       arma::field<arma::vec> S, arma::mat X, arma::vec m,
                       arma::vec VacTime, arma::vec FirstInfTime, arma::vec dimension, int infection_ind, int interact_ind,
                       arma::field<arma::vec> knots, arma::mat ResMat, double eps, int MaxIter,
                       bool constantVE, bool interact, double cutoff){
  Rcpp::List ret;
  int tau = max(Time);
  int n_type = knots.size();
  arma::mat vh(tau, n_type); vh.zeros();
  arma::mat var_vh(tau, n_type); var_vh.zeros();
  arma::mat ci_upper(tau, n_type); ci_upper.zeros();
  arma::mat ci_lower(tau, n_type); ci_lower.zeros();
  int p1 = X.n_cols;
  int p2 = arma::sum(dimension);
  int n = X.n_rows;
  int n_all = arma::sum(m);
  arma::mat Z(n_all,p2); Z.zeros();
  int K = t.n_elem;
  arma::vec d(K); d.zeros();
  arma::vec m_sum(n); m_sum.zeros();
  arma::vec beta(p1); beta.zeros();
  arma::vec gamma(p2); gamma.zeros();
  arma::vec theta(p1+p2); theta.zeros();
  arma::vec theta_old(p1+p2); theta_old.zeros();
  arma::vec Score_constant(p1+p2); Score_constant.zeros();
  CoxPHInformation temp_old;
  CoxPHInformation temp_new;
  int Iter = 0;
  double error = 1;

  // compute d
  arma::uvec index_event = arma::find(Delta > 0);
  arma::vec Time_event = Time(index_event);
  for(int k = 0; k < K; ++k){
    arma::uvec index_tk = arma::find(Time_event == t(k));
    d(k) = index_tk.n_elem;
  }

  // compute m_sum
  for(int i = 0; i < n; ++i){
    double start = (i == 0) ? 0 : m_sum(i-1);
    m_sum(i) = start + m(i);
  }

  // compute Z_i(T_ij)
  for(int i = 0; i < n; ++i){
    int start_i = (i == 0) ? 0 : m_sum(i-1);
    int end_i = start_i + m(i) - 1;
    arma::vec Time_i = Time.subvec(start_i, end_i);
    arma::vec V_i = V(i);
    arma::vec S_i = S(i);
    double FirstInf_i = FirstInfTime(i);
    for(int j = 0; j < m(i); ++j){
      Z.row(start_i+j) = BS2(Time_i(j), V_i, S_i, FirstInf_i, knots, dimension, constantVE, infection_ind, interact_ind, interact);
    }
  }

  // compute Score_constant
  for(int i = 0; i < n; ++i){
    int start_i = (i == 0) ? 0 : m_sum(i-1);
    for(int j = 0; j < m(i); ++j){
      if(Delta(start_i+j) == 1){
        arma::rowvec XZ = join_horiz(X.row(i), Z.row(start_i+j));
        Score_constant += XZ.t();
      }
    }
  }

  temp_old = get_Info(theta_old.head(p1), theta_old.tail(p2), d,
                      Time, t, Delta, C, V, S, X, m,
                      VacTime, FirstInfTime, dimension, infection_ind, interact_ind,
                      Z, Score_constant, knots, constantVE, interact, cutoff);
  while(error > eps && Iter < MaxIter){
    theta = theta_old + ResMat.t() * arma::solve(ResMat * temp_old.InformationMat * ResMat.t(), ResMat * temp_old.ScoreFun);
    std::cout << "beta = " << theta.head(p1).t() << std::endl;
    std::cout << "gamma = " << theta.tail(p2).t() << std::endl;
    temp_new = get_Info(theta.head(p1), theta.tail(p2), d,
                        Time, t, Delta, C, V, S, X, m,
                        VacTime, FirstInfTime, dimension, infection_ind, interact_ind,
                        Z, Score_constant, knots, constantVE, interact, cutoff);
    error = arma::abs(theta.head(p1)-theta_old.head(p1)).max() + arma::abs(theta.tail(p2)-theta_old.tail(p2)).max();
    theta_old = theta;
    temp_old = temp_new;
    ++Iter;
    std::cout << "Iteration = " << Iter << ", Difference = " << error << "\n" << std::endl;
  }

  beta = theta.head(p1);
  gamma = theta.tail(p2);
  arma::mat cov = ResMat.t() * arma::inv(ResMat * temp_new.InformationMat * ResMat.t()) * ResMat;
  arma::mat cov_gamma = cov.submat(p1,p1,p1+p2-1,p1+p2-1);

  for(int k = 0; k < n_type; ++k){
    arma::vec knots_k = knots(k);
    int start_k = (k == 0) ? 0 : arma::sum(dimension.subvec(0,k-1));
    int end_k = start_k + dimension(k) - 1;
    arma::vec gamma_k = gamma.subvec(start_k, end_k);
    arma::mat cov_gamma_k = cov_gamma.submat(start_k, start_k, end_k, end_k);
    for(int l = 1; l <= tau; ++l){
      arma::rowvec B = (k >= infection_ind-1) ? BS_infection(l, knots_k, constantVE) : BS(l, knots_k, constantVE);
      vh(l-1,k) = 1 - exp(arma::as_scalar(B * gamma_k));
      var_vh(l-1,k) = exp(2*arma::as_scalar(B * gamma_k)) * arma::as_scalar(B * cov_gamma_k * B.t());
      double var_logh = arma::as_scalar(B * cov_gamma_k * B.t());
      ci_upper(l-1,k) = 1 - exp(arma::as_scalar(B * gamma_k) - sqrt(var_logh) * 1.96);
      ci_lower(l-1,k) = 1 - exp(arma::as_scalar(B * gamma_k) + sqrt(var_logh) * 1.96);
    }
  }

  ret["beta"] = beta;
  ret["gamma"] = gamma;
  ret["Covariance"] = cov;
  ret["vh"] = vh;
  ret["var_vh"] = var_vh;
  ret["ci_upper"] = ci_upper;
  ret["ci_lower"] = ci_lower;

  return ret;
}
