#include "saint.h"

/********* ALPHA_IP *********/

void DP_mu_gamma(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r) {
  int i,j;
  int wsum[_MAX_COMP_];
  int wrevsum[_MAX_COMP_];
  float gammap[_MAX_COMP_];
  for(i=0;i<_MAX_COMP_;i++) {
    wsum[i] = 0;
    wrevsum[i] = 0;
  }
  for(i=0;i<data->nprey;i++) (wsum[prior->w_mu[i]])++;
  for(i=_MAX_COMP_-1;i>=0;i--) {
    for(j=i;j<_MAX_COMP_;j++) wrevsum[i] += wsum[j];
  }
  for(i=0;i<_MAX_COMP_-1;i++) gammap[i] = gsl_ran_beta(r, (1.0 + (double) wsum[i]), ((double) prior->rho_mu) + ((double) wrevsum[i+1]));
  gammap[_MAX_COMP_-1] = 1.0;
  prior->gamma_mu[0] = gammap[0];
  for(i=1;i<_MAX_COMP_;i++) {
    prior->gamma_mu[i] = gammap[i];
    for(j=0;j<i;j++) prior->gamma_mu[i] *= (1.0 - gammap[j]); 
  }
}

void DP_mu_w(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r, int pid) {
  int i,j,id;
  float cur_mu, tmp_lambda, maxl, tmp;
  float prob[_MAX_COMP_];
  for(i=0;i<_MAX_COMP_;i++) prob[i] = log(prior->gamma_mu[i]);
  cur_mu = param->mu[pid];
  for(i=0;i<_MAX_COMP_;i++) {
    for(j=0;j<data->preyNinter[pid];j++) {
      id = data->p2i[pid][j];
      if(data->ctrl[data->i2IP[id]]) {
        tmp_lambda = param->lambda_false[id] + prior->theta_mu[i] - cur_mu;
        tmp = data->d[id];
        prob[i] += log_gaussian(tmp, (tmp_lambda), param->eta0[data->i2p[id]]);
      }
    }
  }
  maxl = vec_max(prob, _MAX_COMP_);
  for(i=0;i<_MAX_COMP_;i++) prob[i] -= maxl;
  for(i=0;i<_MAX_COMP_;i++) prob[i] = exp(prob[i]);
  prior->w_mu[pid] = ranMultinom(r, prob, _MAX_COMP_);
  param->mu[pid] = prior->theta_mu[prior->w_mu[pid]];

  for(j=0;j<data->preyNinter[pid];j++) {
    id = data->p2i[pid][j];
    param->lambda_false[id] += param->mu[pid] - cur_mu;
  }
}

void DP_mu_theta(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r, int pid, int *inuse) {
  int i, j, id, accept;
  float Delta, mhratio, newval, scale, tmp;
  scale = prior->gamma_mu[pid] / (1.0 - prior->gamma_mu[pid]);
  if(inuse[pid] == 0) {
    newval = gsl_ran_gaussian(r, sqrt(prior->v_mu)) + prior->m_mu;
    Delta = newval - prior->theta_mu[pid];
    prior->theta_mu[pid] = newval;
  }
  else {
    /* metropolis-hastings */
    mhratio = 0.0;
    Delta = gsl_ran_gaussian(r, 0.25);
    for(i=0;i<data->nprey;i++) {
      if(prior->w_mu[i] == pid) {
        for(j=0;j<data->preyNinter[i];j++) {
          id = data->p2i[i][j];
          if(data->ctrl[data->i2IP[id]]) {
            tmp = data->d[id];
            param->lambda_false_tmp[id] = param->lambda_false[id] + Delta;
            mhratio += log_gaussian(tmp, (param->lambda_false_tmp[id]), param->eta0[i]) 
                     - log_gaussian(tmp, (param->lambda_false[id]), param->eta0[i]);
          }
        }
      }
    }
    mhratio += log_gaussian(prior->theta_mu[pid] + Delta, prior->m_mu, prior->v_mu) 
             - log_gaussian(prior->theta_mu[pid], prior->m_mu, prior->v_mu); 
    accept = gsl_ran_flat(r, 0.0, 1.0) <= GSL_MIN(1.0, exp(mhratio)) ? 1 : 0 ;

    /* if accepted, update param and lambda */
    if(accept) {
      prior->theta_mu[pid] += Delta;
      for(i=0;i<data->nprey;i++) {
        if(prior->w_mu[i] == pid) {
          param->mu[i] += Delta;
          for(j=0;j<data->preyNinter[i];j++) {
            id = data->p2i[i][j];
            param->lambda_false[id] += Delta; 
          }
        }
      }
    }
  }
}

void DP_mu(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r) {
  int i;
  float mean;
  int inuse[_MAX_COMP_];
  for(i=0;i<_MAX_COMP_;i++) inuse[i] = 0;
  for(i=0;i<data->nprey;i++) inuse[prior->w_mu[i]] = 1;

  DP_mu_gamma(param, prior, data, r);
  for(i=0;i<data->nprey;i++) DP_mu_w(param, prior, data, r, i);
  for(i=0;i<_MAX_COMP_;i++) DP_mu_theta(param, prior, data, r, i, inuse);
  /* loglik update */

  mean = 0.0;
  for(i=0;i<_MAX_COMP_;i++) mean += prior->gamma_mu[i] * prior->theta_mu[i];
  for(i=0;i<_MAX_COMP_;i++) prior->theta_mu[i] -= mean;
  for(i=0;i<data->nprey;i++) param->mu[i] -= mean;
  param->betac += mean;

}


