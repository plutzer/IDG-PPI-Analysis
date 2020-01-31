#include "saint.h"

/***************************************************/
/*        computing likelihoods in log scale       */
/***************************************************/

float log_poisson_prop(float N, float lambda) {
  float res = -lambda + N * log(lambda);
  return res;
}

float log_poisson_g_prop(float N, float lambda, float theta) {
  float lambda1, lambda2, out;
  lambda2 = 1.0 - 1.0 / sqrt(theta);
  lambda1 = lambda / sqrt(theta);
  out = log(lambda1) + (N - 1.0) * log(lambda1 + N * lambda2) - (lambda1 + N * lambda2);
  return out;
}


/*************************/
/*    all interactions   */
/*************************/
float LRprop(PARAM *param, PRIOR *prior, DATA *data) {
  int i,j,id;
  float pos, neg, maxl, tmp;
  float false_avg, true_avg;
  float lik_new, lik_old;
  lik_new = 0.0;
  lik_old = 0.0;
  for(i=0;i<data->nuinter;i++) {
    pos = 0.0; neg = 0.0;
    for(j=0;j<data->n_u2a[i];j++) { 
      id = data->u2a[i][j];
      if(1) {
        pos += log_poisson_g_prop(data->d[id], exp(param->lambda_true[id]), param->eta[data->i2p[id]]);
        neg += log_poisson_g_prop(data->d[id], exp(param->lambda_false[id]), param->eta0[data->i2p[id]]);
      }
    }
    maxl = pos > neg ? pos : neg;
    pos = exp(pos - maxl);
    neg = exp(neg - maxl);
    lik_new += log(param->ptrue_tmp * pos + (1.0-param->ptrue_tmp) * neg);
    lik_old += log(param->ptrue * pos + (1.0-param->ptrue) * neg);
  }


  for(i=0;i<data->nprey;i++) {
    if(data->preyFlag[i] == 0) {
      false_avg = 0.0;
      for(j=0;j<data->preyNinter[i];j++) {
        id = data->p2i[i][j];
        false_avg += param->lambda_false[id];
      }
      false_avg /= ((float) data->preyNinter[i]);
      true_avg = 0.0;
      for(j=0;j<data->preyNinter[i];j++) {
        id = data->p2i[i][j];
        true_avg += param->lambda_true[id];
      }
      true_avg /= ((float) data->preyNinter[i]);
      tmp = ((float) (data->nIP - data->preyNinter[i]));
      pos = tmp * log_poisson_g_prop(0.0, exp(true_avg), param->eta[i]);
      neg = tmp * log_poisson_g_prop(0.0, exp(false_avg), param->eta0[i]);
      maxl = pos > neg ? pos : neg;
      pos = exp(pos - maxl);
      neg = exp(neg - maxl);
      lik_new += log(param->ptrue_tmp * pos + (1.0-param->ptrue_tmp) * neg);
      lik_old += log(param->ptrue * pos + (1.0-param->ptrue) * neg);
    }
  }
  return lik_new - lik_old;
}

float loglik_all(PARAM *param, PRIOR *prior, DATA *data) {
  int i;
  float lambda;
  float lik = 0.0;
  for(i=0;i<data->ninter;i++) {
    if(param->Z[data->a2u[i]]) {
      lambda = exp(param->lambda_true[i]);
      lik += log_poisson_g_prop(data->d[i], lambda, param->eta[data->i2p[i]]);
    }
    else {
      lambda = exp(param->lambda_false[i]);
      lik += log_poisson_g_prop(data->d[i], lambda, param->eta0[data->i2p[i]]);
    }
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
      if(param->Z[data->a2u[i]]) {
        lambda = exp(param->lambda_true[i]);
        lik += log_poisson_g_prop(data->d[i], lambda, param->eta[data->i2p[i]]);
      }
    }
    else {
      if(param->Z[data->a2u[i]] == 0) {
        lambda = exp(param->lambda_false[i]);
        lik += log_poisson_g_prop(data->d[i], lambda, param->eta0[data->i2p[i]]);
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
      if(param->Z[data->a2u[i]]) {
        lambda = exp(param->lambda_true_tmp[i]);
        lik += log_poisson_g_prop(data->d[i], lambda, param->eta[data->i2p[i]]);
      }
    }
    else {
      if(param->Z[data->a2u[i]] == 0) {
        lambda = exp(param->lambda_false_tmp[i]);
        lik += log_poisson_g_prop(data->d[i], lambda, param->eta0[data->i2p[i]]);
      }
    }
  }
  return lik;
}

















