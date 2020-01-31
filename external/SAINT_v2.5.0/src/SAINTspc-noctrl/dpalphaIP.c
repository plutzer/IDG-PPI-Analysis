#include "saint.h"

/********* ALPHA_IP *********/

void DP_alpha_IP_gamma(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r) {
  int i,j;
  int wsum[_MAX_COMP_];
  int wrevsum[_MAX_COMP_];
  float gammap[_MAX_COMP_];
  for(i=0;i<_MAX_COMP_;i++) {
    wsum[i] = 0;
    wrevsum[i] = 0;
  }
  for(i=0;i<data->nIP;i++) {
    if(data->ctrl[i] == 0) (wsum[prior->w_alpha_IP[i]])++;
  }
  for(i=_MAX_COMP_-1;i>=0;i--) {
    for(j=i;j<_MAX_COMP_;j++) wrevsum[i] += wsum[j];
  }
  for(i=0;i<_MAX_COMP_-1;i++) gammap[i] = gsl_ran_beta(r, (1.0 + (double) wsum[i]), ((double) prior->rho_alpha_IP) + ((double) wrevsum[i+1]));
  gammap[_MAX_COMP_-1] = 1.0;
  prior->gamma_alpha_IP[0] = gammap[0];
  for(i=1;i<_MAX_COMP_;i++) {
    prior->gamma_alpha_IP[i] = gammap[i];
    for(j=0;j<i;j++) prior->gamma_alpha_IP[i] *= (1.0 - gammap[j]); 
  }
}

void DP_alpha_IP_w(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r, int pid) {
  int i,j,id;
  float cur_alpha_IP, tmp_lambda, maxl;
  float prob[_MAX_COMP_];
  for(i=0;i<_MAX_COMP_;i++) prob[i] = log(prior->gamma_alpha_IP[i]);
  cur_alpha_IP = param->alpha_IP[pid];
  for(i=0;i<_MAX_COMP_;i++) {
    for(j=0;j<data->IPNinter[pid];j++) {
      id = data->IP2i[pid][j];
      if(param->Z[data->a2u[id]] && data->d[id] > 0.0) {
        tmp_lambda = param->lambda_true[id] + prior->theta_alpha_IP[i] - cur_alpha_IP;
        prob[i] += log_poisson_g_prop(data->d[id], exp(tmp_lambda), param->eta[data->i2p[id]]);
      }
    }
  }
  maxl = vec_max(prob, _MAX_COMP_);
  for(i=0;i<_MAX_COMP_;i++) prob[i] -= maxl;
  for(i=0;i<_MAX_COMP_;i++) prob[i] = exp(prob[i]);
  prior->w_alpha_IP[pid] = ranMultinom(r, prob, _MAX_COMP_);
  param->alpha_IP[pid] = prior->theta_alpha_IP[prior->w_alpha_IP[pid]];

  for(j=0;j<data->IPNinter[pid];j++) {
    id = data->IP2i[pid][j];
    param->lambda_true[id] += param->alpha_IP[pid] - cur_alpha_IP;
  }
}

void DP_alpha_IP_theta(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r, int pid, int *inuse) {
  int i, j, id, accept;
  float Delta, mhratio, newval, scale;
  scale = prior->gamma_alpha_IP[pid] / (1.0 - prior->gamma_alpha_IP[pid]);
  if(inuse[pid] == 0) {
    newval = gsl_ran_gaussian(r, sqrt(prior->v_alpha_IP)) + prior->m_alpha_IP;
    Delta = newval - prior->theta_alpha_IP[pid];
    prior->theta_alpha_IP[pid] = newval;
  }
  else {
    /* metropolis-hastings */
    mhratio = 0.0;
    Delta = gsl_ran_gaussian(r, 0.25);
    for(i=0;i<data->nIP;i++) {
      if(prior->w_alpha_IP[i] == pid) {
        for(j=0;j<data->IPNinter[i];j++) {
          id = data->IP2i[i][j];
          if(param->Z[data->a2u[id]] && data->d[id] > 0.0) {
            param->lambda_true_tmp[id] = param->lambda_true[id] + Delta;
            mhratio += log_poisson_g_prop(data->d[id], exp(param->lambda_true_tmp[id]), param->eta[data->i2p[id]]) 
                     - log_poisson_g_prop(data->d[id], exp(param->lambda_true[id]), param->eta[data->i2p[id]]);
          }
        }
      }
    }
    mhratio += log_gaussian(prior->theta_alpha_IP[pid] + Delta, prior->m_alpha_IP, prior->v_alpha_IP) 
             - log_gaussian(prior->theta_alpha_IP[pid], prior->m_alpha_IP, prior->v_alpha_IP); 
    accept = gsl_ran_flat(r, 0.0, 1.0) <= GSL_MIN(1.0, exp(mhratio)) ? 1 : 0 ;

    /* if accepted, update param and lambda */
    if(accept) {
      prior->theta_alpha_IP[pid] += Delta;
      for(i=0;i<data->nIP;i++) {
        if(prior->w_alpha_IP[i] == pid && data->ctrl[i] == 0) {
          param->alpha_IP[i] += Delta;
          for(j=0;j<data->IPNinter[i];j++) {
            id = data->IP2i[i][j];
            param->lambda_true[id] += Delta; 
          }
        }
      }
    }
  }
}

void DP_alpha_IP(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r) {
  int i;
  float mean;
  int inuse[_MAX_COMP_];
  for(i=0;i<_MAX_COMP_;i++) inuse[i] = 0;
  for(i=0;i<data->nIP;i++) inuse[prior->w_alpha_IP[i]] = 1;

  DP_alpha_IP_gamma(param, prior, data, r);
  for(i=0;i<data->nIP;i++) {
    if(data->ctrl[i] == 0) DP_alpha_IP_w(param, prior, data, r, i);
  }
  for(i=0;i<_MAX_COMP_;i++) DP_alpha_IP_theta(param, prior, data, r, i, inuse);
  /* loglik update */

  mean = 0.0;
  for(i=0;i<_MAX_COMP_;i++) mean += prior->gamma_alpha_IP[i] * prior->theta_alpha_IP[i];
  for(i=0;i<_MAX_COMP_;i++) prior->theta_alpha_IP[i] -= mean;
  for(i=0;i<data->nIP;i++) param->alpha_IP[i] -= mean;
  param->beta0 += mean;

}


