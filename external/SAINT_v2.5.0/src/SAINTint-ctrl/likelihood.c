#include "saint.h"

/***************************************************/
/*        computing likelihoods in log scale       */
/***************************************************/


/*************************/
/*    all interactions   */
/*************************/
float LRprop(PARAM *param, PRIOR *prior, DATA *data) {
  int i,j,id;
  float pos, neg, maxl;
  float lik_new, lik_old;
  lik_new = 0.0;
  lik_old = 0.0;
  for(i=0;i<data->nuinter;i++) {
    pos = 0.0; neg = 0.0;
    for(j=0;j<data->n_u2a[i];j++) { 
      id = data->u2a[i][j];
      if(data->ctrl[data->i2IP[id]] == 0) {
        pos += log_gaussian(data->d[id], (param->lambda_true[id]), param->eta[data->i2p[id]]);
        neg += log_gaussian(data->d[id], (param->lambda_false[id]), param->eta0[data->i2p[id]]);
      }
    }
    maxl = pos > neg ? pos : neg;
    pos = exp(pos - maxl);
    neg = exp(neg - maxl);
    lik_new += log(param->ptrue_tmp * pos + (1.0-param->ptrue_tmp) * neg);
    lik_old += log(param->ptrue * pos + (1.0-param->ptrue) * neg);
  }
  return lik_new - lik_old;
}

float loglik_all(PARAM *param, PRIOR *prior, DATA *data) {
  int i;
  float lambda;
  float lik = 0.0;
  for(i=0;i<data->ninter;i++) {
    if(param->Z[data->a2u[i]] && data->miss[i] == 0) {
      lambda = (param->lambda_true[i]);
      lik += log_gaussian(data->d[i], lambda, param->eta[data->i2p[i]]);
    }
    else if (param->Z[data->a2u[i]] == 0) {
      lambda = (param->lambda_false[i]);
      lik += log_gaussian(data->d[i], lambda, param->eta0[data->i2p[i]]);
    }
    else {}
  }
  return lik;
}

float loglik_all_class(PARAM *param, PRIOR *prior, DATA *data, int cl) {
  /* loglik by class */
  int i;
  float lambda;
  float lik = 0.0;
  for(i=0;i<data->ninter;i++) {
    if(cl) {
      if(param->Z[data->a2u[i]] && data->miss[i] == 0) {
        lambda = (param->lambda_true[i]);
        lik += log_gaussian(data->d[i], lambda, param->eta[data->i2p[i]]);
      }
    }
    else {
      if(param->Z[data->a2u[i]] == 0) {
        lambda = (param->lambda_false[i]);
        lik += log_gaussian(data->d[i], lambda, param->eta0[data->i2p[i]]);
      }
    }
  }
  return lik;
}

float loglik_all_class_tmp(PARAM *param, PRIOR *prior, DATA *data, int cl) {
  /* loglik by class */
  int i;
  float lambda;
  float lik = 0.0;
  for(i=0;i<data->ninter;i++) {
    if(cl) {
      if(param->Z[data->a2u[i]] && data->miss[i] == 0) {
        lambda = (param->lambda_true_tmp[i]);
        lik += log_gaussian(data->d[i], lambda, param->eta[data->i2p[i]]);
      }
    }
    else {
      if(param->Z[data->a2u[i]] == 0) {
        lambda = (param->lambda_false_tmp[i]);
        lik += log_gaussian(data->d[i], lambda, param->eta0[data->i2p[i]]);
      }
    }
  }
  return lik;
}

















