#include "saint.h"

/********* ALPHA_prey *********/

void DP_eta0_gamma(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r) {
  int i,j;
  int wsum[_MAX_COMP_];
  int wrevsum[_MAX_COMP_];
  float gammap[_MAX_COMP_];
  for(i=0;i<_MAX_COMP_;i++) {
    wsum[i] = 0;
    wrevsum[i] = 0;
  }
  for(i=0;i<data->nprey;i++) (wsum[prior->w_eta0[i]])++;
  for(i=_MAX_COMP_-1;i>=0;i--) {
    for(j=i;j<_MAX_COMP_;j++) wrevsum[i] += wsum[j];
  }
  for(i=0;i<_MAX_COMP_-1;i++) gammap[i] = gsl_ran_beta(r, (1.0 + (double) wsum[i]), ((double) prior->rho_eta0) + ((double) wrevsum[i+1]));
  gammap[_MAX_COMP_-1] = 1.0;
  prior->gamma_eta0[0] = gammap[0];
  for(i=1;i<_MAX_COMP_;i++) {
    prior->gamma_eta0[i] = gammap[i];
    for(j=0;j<i;j++) prior->gamma_eta0[i] *= (1.0 - gammap[j]); 
  }
}

void DP_eta0_w(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r, int pid) {
  int i,j,id;
  float cur_eta, tmp_lambda, maxl;
  float prob[_MAX_COMP_];
  for(i=0;i<_MAX_COMP_;i++) prob[i] = log(prior->gamma_eta0[i]);
  cur_eta = param->eta0[pid];
  for(i=0;i<_MAX_COMP_;i++) {
    for(j=0;j<data->preyNinter[pid];j++) {
      id = data->p2i[pid][j];
      if(data->ctrl[data->i2IP[id]]) {
        tmp_lambda = param->lambda_false[id];
        prob[i] += log_gaussian(data->d[id], (tmp_lambda), prior->theta_eta0[i]);
      }
    }
  }
  maxl = vec_max(prob, _MAX_COMP_);
  for(i=0;i<_MAX_COMP_;i++) prob[i] -= maxl;
  for(i=0;i<_MAX_COMP_;i++) prob[i] = exp(prob[i]);
  prior->w_eta0[pid] = ranMultinom(r, prob, _MAX_COMP_);
  param->eta0[pid] = prior->theta_eta0[prior->w_eta0[pid]];
}


void DP_eta0_theta(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r, int pid, int *inuse) {
  int i, j, id, accept, pass;
  float Delta, mhratio, newval, scale, tmp_lambda;
  scale = prior->gamma_eta0[pid] / (1.0 - prior->gamma_eta0[pid]);
  if(inuse[pid] == 0) {
    pass = 0;
    while(!pass) {
      newval = 1.0 / gsl_ran_gamma(r, 100.0, 1.0);
      if(newval < 2.0) pass = 1;
    }
    Delta = newval - prior->theta_eta0[pid];
    prior->theta_eta0[pid] = newval;
  }
  else {
    /* metropolis-hastings */
    mhratio = 0.0;
    Delta = gsl_ran_gaussian(r, 0.1);
    if(prior->theta_eta0[pid] + Delta <= 0.0 || prior->theta_eta0[pid] + Delta > 2.0) {
      accept = 0;
    }
    else {
      for(i=0;i<data->nprey;i++) {
        if(prior->w_eta0[i] == pid) {
          for(j=0;j<data->preyNinter[i];j++) {
            id = data->p2i[i][j];
            if(data->ctrl[data->i2IP[id]]) {
              tmp_lambda = param->lambda_false[id];
              mhratio += log_gaussian(data->d[id], (tmp_lambda), prior->theta_eta0[pid]+Delta)
                       - log_gaussian(data->d[id], (tmp_lambda), prior->theta_eta0[pid]);
            }
          }
        }
      }
      mhratio += log_inv_gamma( (prior->theta_eta0[pid]+Delta), prior->shape_eta0, prior->scale_eta0)
               - log_inv_gamma( prior->theta_eta0[pid], prior->shape_eta0, prior->scale_eta0);
      accept = gsl_ran_flat(r, 0.0, 1.0) <= GSL_MIN(1.0, exp(mhratio)) ? 1 : 0 ;
    }

    /* if accepted, update param and lambda */
    if(accept) {
      prior->theta_eta0[pid] += Delta;
      for(i=0;i<data->nprey;i++) {
        if(prior->w_eta0[i] == pid) {
          param->eta0[i] += Delta;
        }
      }
    }
  }
}

void DP_eta0(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r) {
  int i;
  int inuse[_MAX_COMP_];
  for(i=0;i<_MAX_COMP_;i++) inuse[i] = 0;
  for(i=0;i<data->nprey;i++) inuse[prior->w_eta0[i]] = 1;

  DP_eta0_gamma(param, prior, data, r);
  for(i=0;i<data->nprey;i++) DP_eta0_w(param, prior, data, r, i);
  for(i=0;i<_MAX_COMP_;i++) DP_eta0_theta(param, prior, data, r, i, inuse);
  /* loglik update */

}


