#include "saint.h"

/********* ALPHA_prey *********/

void DP_alpha_prey_gamma(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r) {
  int i,j;
  int wsum[_MAX_COMP_];
  int wrevsum[_MAX_COMP_];
  float gammap[_MAX_COMP_];
  for(i=0;i<_MAX_COMP_;i++) {
    wsum[i] = 0;
    wrevsum[i] = 0;
  }
  for(i=0;i<data->nprey;i++) (wsum[prior->w_alpha_prey[i]])++;
  for(i=_MAX_COMP_-1;i>=0;i--) {
    for(j=i;j<_MAX_COMP_;j++) wrevsum[i] += wsum[j];
  }
  for(i=0;i<_MAX_COMP_-1;i++) gammap[i] = gsl_ran_beta(r, (1.0 + (double) wsum[i]), ((double) prior->rho_alpha_prey) + ((double) wrevsum[i+1]));
  gammap[_MAX_COMP_-1] = 1.0;
  prior->gamma_alpha_prey[0] = gammap[0];
  for(i=1;i<_MAX_COMP_;i++) {
    prior->gamma_alpha_prey[i] = gammap[i];
    for(j=0;j<i;j++) prior->gamma_alpha_prey[i] *= (1.0 - gammap[j]); 
  }
}

void DP_alpha_prey_w(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r, int pid) {
  int i,j,id;
  float cur_alpha_prey, tmp_lambda, maxl, tmp;
  float prob[_MAX_COMP_];
  for(i=0;i<_MAX_COMP_;i++) prob[i] = log(prior->gamma_alpha_prey[i]);
  cur_alpha_prey = param->alpha_prey[pid];
  for(i=0;i<_MAX_COMP_;i++) {
    for(j=0;j<data->preyNinter[pid];j++) {
      id = data->p2i[pid][j];
      if(param->Z[data->a2u[id]] && data->d[id] > 0.0) {
        tmp_lambda = param->lambda_true[id] + prior->theta_alpha_prey[i] - cur_alpha_prey;
        if(lowMode) tmp = GSL_MIN(_LM_, data->d[id]);
        else tmp = data->d[id];
        if(tmp > 0.0) prob[i] += log_poisson_g_prop(tmp, exp(tmp_lambda), param->eta[pid]);
      }
      else if (data->d[id] >= 2.0 && data->ctrlavg[pid] > 0.0 && (data->d[id] / data->ctrlavg[pid]) > 2.0) {
        tmp_lambda = param->lambda_true[id] + prior->theta_alpha_prey[i] - cur_alpha_prey;
        if(lowMode) tmp = GSL_MIN(_LM_, data->d[id]);
        else tmp = data->d[id];
        if(tmp > 0.0) prob[i] += log_poisson_g_prop(tmp, exp(tmp_lambda), param->eta[pid]);
      }
      else {}
    }
  }
  maxl = vec_max(prob, _MAX_COMP_);
  for(i=0;i<_MAX_COMP_;i++) prob[i] -= maxl;
  for(i=0;i<_MAX_COMP_;i++) prob[i] = exp(prob[i]);
  prior->w_alpha_prey[pid] = ranMultinom(r, prob, _MAX_COMP_);
  param->alpha_prey[pid] = prior->theta_alpha_prey[prior->w_alpha_prey[pid]];

  for(j=0;j<data->preyNinter[pid];j++) {
    id = data->p2i[pid][j];
    param->lambda_true[id] += param->alpha_prey[pid] - cur_alpha_prey;
  }
}

void DP_alpha_prey_theta(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r, int pid, int *inuse) {
  int i, j, id, accept;
  float Delta, mhratio, newval, scale;
  scale = prior->gamma_alpha_prey[pid] / (1.0 - prior->gamma_alpha_prey[pid]);
  if(inuse[pid] == 0) {
    newval = gsl_ran_gaussian(r, sqrt(prior->v_alpha_prey)) + prior->m_alpha_prey;
    Delta = newval - prior->theta_alpha_prey[pid];
    prior->theta_alpha_prey[pid] = newval;
  }
  else {
    /* metropolis-hastings */
    mhratio = 0.0;
    Delta = gsl_ran_gaussian(r, 1.0);
    for(i=0;i<data->nprey;i++) {
      if(prior->w_alpha_prey[i] == pid) {
        for(j=0;j<data->preyNinter[i];j++) {
          id = data->p2i[i][j];
          if(param->Z[data->a2u[id]] && data->d[id] > 0.0) {
            param->lambda_true_tmp[id] = param->lambda_true[id] + Delta;
            if(lowMode) {
              mhratio += log_poisson_g_prop(GSL_MIN(_LM_,data->d[id]), exp(param->lambda_true_tmp[id]), param->eta[i]) 
			   - log_poisson_g_prop(GSL_MIN(_LM_,data->d[id]), exp(param->lambda_true[id]), param->eta[i]);
            }
            else {
              mhratio += log_poisson_g_prop(data->d[id], exp(param->lambda_true_tmp[id]), param->eta[i]) 
			   - log_poisson_g_prop(data->d[id], exp(param->lambda_true[id]), param->eta[i]);
            }
          }
        }
      }
    }
    mhratio += log_gaussian(prior->theta_alpha_prey[pid] + Delta, prior->m_alpha_prey, prior->v_alpha_prey) 
             - log_gaussian(prior->theta_alpha_prey[pid], prior->m_alpha_prey, prior->v_alpha_prey); 
    accept = gsl_ran_flat(r, 0.0, 1.0) <= GSL_MIN(1.0, exp(mhratio)) ? 1 : 0 ;

    /* if accepted, update param and lambda */
    if(accept) {
      prior->theta_alpha_prey[pid] += Delta;
      for(i=0;i<data->nprey;i++) {
        if(prior->w_alpha_prey[i] == pid) {
          param->alpha_prey[i] += Delta;
          for(j=0;j<data->preyNinter[i];j++) {
            id = data->p2i[i][j];
            param->lambda_true[id] += Delta; 
          }
        }
      }
    }
  }
}

void DP_alpha_prey(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r) {
  int i;
  float mean;

  int inuse[_MAX_COMP_];
  for(i=0;i<_MAX_COMP_;i++) inuse[i] = 0;
  for(i=0;i<data->nprey;i++) {
    inuse[prior->w_alpha_prey[i]] = 1;
  }

  DP_alpha_prey_gamma(param, prior, data, r);
  for(i=0;i<data->nprey;i++) DP_alpha_prey_w(param, prior, data, r, i);
  for(i=0;i<_MAX_COMP_;i++) DP_alpha_prey_theta(param, prior, data, r, i, inuse);
  /* loglik update */

  mean = 0.0;
  for(i=0;i<_MAX_COMP_;i++) mean += prior->gamma_alpha_prey[i] * prior->theta_alpha_prey[i];
  for(i=0;i<_MAX_COMP_;i++) prior->theta_alpha_prey[i] -= mean;
  for(i=0;i<data->nprey;i++) param->alpha_prey[i] -= mean;
  param->beta0 += mean;

}


