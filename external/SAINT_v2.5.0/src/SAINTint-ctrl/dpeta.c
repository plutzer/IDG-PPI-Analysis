#include "saint.h"

/********* ALPHA_prey *********/

void DP_eta_gamma(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r) {
  int i,j;
  int wsum[_MAX_COMP_];
  int wrevsum[_MAX_COMP_];
  float gammap[_MAX_COMP_];
  for(i=0;i<_MAX_COMP_;i++) {
    wsum[i] = 0;
    wrevsum[i] = 0;
  }
  for(i=0;i<data->nprey;i++) (wsum[prior->w_eta[i]])++;
  for(i=_MAX_COMP_-1;i>=0;i--) {
    for(j=i;j<_MAX_COMP_;j++) wrevsum[i] += wsum[j];
  }
  for(i=0;i<_MAX_COMP_-1;i++) gammap[i] = gsl_ran_beta(r, (1.0 + (double) wsum[i]), ((double) prior->rho_eta) + ((double) wrevsum[i+1]));
  gammap[_MAX_COMP_-1] = 1.0;
  prior->gamma_eta[0] = gammap[0];
  for(i=1;i<_MAX_COMP_;i++) {
    prior->gamma_eta[i] = gammap[i];
    for(j=0;j<i;j++) prior->gamma_eta[i] *= (1.0 - gammap[j]); 
  }
}

void DP_eta_w(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r, int pid) {
  int i,j,id;
  float cur_eta, tmp_lambda, maxl;
  float prob[_MAX_COMP_];
  for(i=0;i<_MAX_COMP_;i++) prob[i] = log(prior->gamma_eta[i]);
  cur_eta = param->eta[pid];
  for(i=0;i<_MAX_COMP_;i++) {
    for(j=0;j<data->preyNinter[pid];j++) {
      id = data->p2i[pid][j];
      if(param->Z[data->a2u[id]]) tmp_lambda = param->lambda_true[id];
      else tmp_lambda = param->lambda_false[id];
      prob[i] += log_gaussian(data->d[id], (tmp_lambda), prior->theta_eta[i]);
    }
  }
  maxl = vec_max(prob, _MAX_COMP_);
  for(i=0;i<_MAX_COMP_;i++) prob[i] -= maxl;
  for(i=0;i<_MAX_COMP_;i++) prob[i] = exp(prob[i]);
  prior->w_eta[pid] = ranMultinom(r, prob, _MAX_COMP_);
  param->eta[pid] = prior->theta_eta[prior->w_eta[pid]];
}

float log_exponential(float x, float mean) {
  float res = -log(mean) - x / mean;
  return res;
}

void DP_eta_theta(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r, int pid, int *inuse) {
  int i, j, id, accept, pass;
  float Delta, mhratio, newval, scale, tmp_lambda, tmp;
  scale = prior->gamma_eta[pid] / (1.0 - prior->gamma_eta[pid]);
  if(inuse[pid] == 0) {
    pass = 0;
    while(!pass) {
      newval = 1.0 / gsl_ran_gamma(r, 100.0, 1.0);
      if(newval < 2.0) pass = 1;
    }
    Delta = newval - prior->theta_eta[pid];
    prior->theta_eta[pid] = newval;
  }
  else {
    /* metropolis-hastings */
    mhratio = 0.0;
    Delta = gsl_ran_gaussian(r, 0.1);
    if(prior->theta_eta[pid] + Delta <= 0.0 || prior->theta_eta[pid] + Delta > 2.0) {
      accept = 0;
    }
    else {
      for(i=0;i<data->nprey;i++) {
        if(prior->w_eta[i] == pid) {
          for(j=0;j<data->preyNinter[i];j++) {
            id = data->p2i[i][j];
            if(param->Z[data->a2u[id]] && data->miss[id] == 0) {
              tmp_lambda = param->lambda_true[id];
              tmp = data->d[id];
              mhratio += log_gaussian(tmp, (tmp_lambda), prior->theta_eta[pid]+Delta)
                       - log_gaussian(tmp, (tmp_lambda), prior->theta_eta[pid]);
            }
          }
        }
      }
      mhratio += log_inv_gamma( (prior->theta_eta[pid]+Delta), prior->shape_eta, prior->scale_eta)
               - log_inv_gamma( prior->theta_eta[pid], prior->shape_eta, prior->scale_eta);
      accept = gsl_ran_flat(r, 0.0, 1.0) <= GSL_MIN(1.0, exp(mhratio)) ? 1 : 0 ;
    }

    /* if accepted, update param and lambda */
    if(accept) {
      prior->theta_eta[pid] += Delta;
      for(i=0;i<data->nprey;i++) {
        if(prior->w_eta[i] == pid) {
          param->eta[i] += Delta;
        }
      }
    }
  }
}

void DP_eta(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r) {
  int i;
  int inuse[_MAX_COMP_];
  for(i=0;i<_MAX_COMP_;i++) inuse[i] = 0;
  for(i=0;i<data->nprey;i++) inuse[prior->w_eta[i]] = 1;

  DP_eta_gamma(param, prior, data, r);
  for(i=0;i<data->nprey;i++) DP_eta_w(param, prior, data, r, i);
  for(i=0;i<_MAX_COMP_;i++) DP_eta_theta(param, prior, data, r, i, inuse);
  /* loglik update */

}


