#include "saint.h"

float log_gaussian(float x, float mu, float var) {
  float res = - .5 * pow(x-mu,2.0) / var - .5 * log(2.0 * M_PI * var);
  return res;
}

void sampleBeta0(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r) {
  int i,accept;
  float likratio, lr_new, lr_old;
  float diff = gsl_ran_gaussian(r,0.2);
  
  for(i=0;i<data->ninter;i++) param->lambda_true_tmp[i] = param->lambda_true[i] + diff;
  lr_new = loglik_all_class_tmp(param, prior, data, 1);
  lr_old = loglik_all_class(param, prior, data, 1);
  likratio = lr_new - lr_old; 
  likratio += log_gaussian(param->beta0 + diff, prior->m_beta, prior->v_beta) 
                     - log_gaussian(param->beta0, prior->m_beta, prior->v_beta);
  likratio = GSL_MIN(1.0, exp(likratio));
  accept = gsl_ran_flat(r,0.0,1.0) <= likratio ? 1 : 0;
  if(accept) {
    param->beta0 += diff;
    for(i=0;i<data->ninter;i++) param->lambda_true[i] = param->lambda_true_tmp[i];
    param->loglikTotal += (lr_new - lr_old);
  }  
}

void sampleBetac(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r) {
  int i,accept;
  float likratio, lr_new, lr_old;
  float diff = gsl_ran_gaussian(r,0.2);
  
  for(i=0;i<data->ninter;i++) param->lambda_false_tmp[i] = param->lambda_false[i] + diff;
  lr_new = loglik_all_class_tmp(param, prior, data, 0);
  lr_old = loglik_all_class(param, prior, data, 0);
  likratio = lr_new - lr_old; 
  likratio += log_gaussian(param->betac + diff, prior->m_beta, prior->v_beta) 
                     - log_gaussian(param->betac, prior->m_beta, prior->v_beta);
  likratio = GSL_MIN(1.0, exp(likratio));
  accept = gsl_ran_flat(r,0.0,1.0) <= likratio ? 1 : 0;
  if(accept) {
    param->betac += diff;
    for(i=0;i<data->ninter;i++) param->lambda_false[i] = param->lambda_false_tmp[i];
    param->loglikTotal += (lr_new - lr_old);
  }  
}


void sampleZ(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r) {
  /* Z and iZ */
  int i,j,id;
  int indiv, total;
  int isCtrl, isReverse, isMaxOne;
  float prob, maxl;
  float posi, negi;
  float pos, neg, tmp, tmp_lambda, tmp_neg;
  int cond1;
  for(i=0;i<data->nuinter;i++) {
    isCtrl = 0;
    for(j=0;j<data->n_u2a[i];j++) {
      id = data->u2a[i][j];
      if(data->ctrl[data->i2IP[id]] == 1) {
        isCtrl = 1;
        break;
      }
    }
    isReverse = 0;
    for(j=0;j<data->n_u2a[i];j++) {
      id = data->u2a[i][j];
      cond1 = data->d[id] < GSL_MAX(data->ctrlavg[data->i2p[id]], exp(param->lambda_false[id])) ? 1 : 0;
      if(cond1) {
        // isReverse = 1; 
        break;
      }
    }
    isMaxOne = 1;
    for(j=0;j<data->n_u2a[i];j++) {
      if(data->d2[data->u2a[i][j]] > 1.0) isMaxOne = 0;
    }

    if(isCtrl || isReverse || isMaxOne) {
      param->Z[i] = 0;
      for(j=0;j<data->n_u2a[i];j++) {
        id = data->u2a[i][j];
        param->iZ[id] = 0;
      }
    }
    else {
      pos = 0.0; neg = 0.0;
      for(j=0;j<data->n_u2a[i];j++) {
        id = data->u2a[i][j];
        tmp = data->d2[id];
        tmp_lambda = exp(param->lambda_true[id]);
        if(minFold) {
          if(tmp_lambda < GSL_MAX(exp(param->lambda_false[id]), data->ctrlavg[data->i2p[id]]) * _fold_) {
            tmp_lambda = GSL_MAX(exp(param->lambda_false[id]), data->ctrlavg[data->i2p[id]]) * _fold_;
          }
        }
    
        if(tmp_lambda < GSL_MAX(data->ctrlavg[data->i2p[id]], exp(param->lambda_false[id]))) {
          param->iZ[id] = 0;
        }
        else {
          tmp = data->d2[id];
          tmp_neg = GSL_MAX(0.1, exp(param->lambda_false[id]));
          if(tmp > tmp_lambda && tmp_lambda > tmp_neg) tmp = tmp_lambda;
          if(lowMode) tmp = GSL_MIN(_LM_, tmp);
          posi = log_poisson_g_prop(tmp, tmp_lambda, param->eta[data->i2p[id]]);
          negi = log_poisson_g_prop(tmp, tmp_neg, param->eta0[data->i2p[id]]);
          pos += posi;
          neg += negi;
          maxl = posi > negi ? posi : negi;
          posi -= maxl;
          negi -= maxl;
          prob = param->ptrue * exp(posi) / (param->ptrue * exp(posi) + (1.0-param->ptrue) * exp(negi));
          param->iZ[id] = gsl_ran_flat(r,0.0,1.0) <= prob ? 1 : 0;
          if(data->d2[id] == 0.0) param->iZ[id] = 0.0;
        }
      }
      /* Z */
      if(data->n_u2a[i] == 1) {
        id = data->u2a[i][0];
        param->Z[i] = param->iZ[id];	
      }
      else {
        maxl = pos > neg ? pos : neg;
        pos -= maxl;
        neg -= maxl;
        prob = param->ptrue * exp(pos) / (param->ptrue * exp(pos) + (1.0-param->ptrue) * exp(neg));
        param->Z[i] = gsl_ran_flat(r,0.0,1.0) <= prob ? 1 : 0;
        indiv = 0;
        total = data->n_u2a[i];
        for(j=0;j<data->n_u2a[i];j++) {
          id = data->u2a[i][j];
          if(param->iZ[id]) indiv++;
        }
        pos = ((double) indiv) / ((double) total);
        param->Z[i] = gsl_ran_bernoulli(r, pos);
      }
    }
  }
}

float logit(float x) {
  return log(x) - log(1-x);
}

float inverseLogit(float x) {
  return exp(x) / (1.0 + exp(x));
}

void sampleProportion(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r) {
  int accept;
  float mhratio;
  param->ptrue_tmp = inverseLogit(logit(param->ptrue) + gsl_ran_gaussian(r, 0.1));
  mhratio = LRprop(param, prior, data); 
  /* uniform prior, so no prior ratio, indep. symetric random walk, so no proposal ratio */
  accept = gsl_ran_flat(r,0.0,1.0) <= GSL_MIN(1.0, exp(mhratio)) ? 1 : 0;
  if(accept) param->ptrue = param->ptrue_tmp;  
}

/**************************************/
/*** Metropolis-Hastings with Gibbs ***/
/**************************************/
void mhgibbs(PARAM *param, PRIOR *prior, DATA *data, SUMMARY *summary, const gsl_rng *r, int updateSum) {
  if(gsl_ran_flat(r,0.0,1.0) <= 0.33) sampleBeta0(param, prior, data, r);
  if(gsl_ran_flat(r,0.0,1.0) <= 0.33) sampleBetac(param, prior, data, r);
  DP_alpha_prey(param, prior, data, r);
  /* DP_alpha_IP(param, prior, data, r); */
  DP_mu(param, prior, data, r); 
  //if(gsl_ran_flat(r,0.0,1.0) <= 0.33) DP_eta(param, prior, data, r); 
  //if(gsl_ran_flat(r,0.0,1.0) <= 0.33) DP_eta0(param, prior, data, r); 
  sampleZ(param, prior, data, r);
  compute_lambda_all(param, prior, data);
  if(gsl_ran_flat(r,0.0,1.0) <= 0.33) sampleProportion(param, prior, data, r);
}

void write_mcmc(PARAM *param, PRIOR *prior, DATA *data, FILE *fp1, FILE *fp2, FILE *fp3, int ct) {
  int i;
  fprintf(fp1, "%d\t", ct+1);
  for(i=0;i<data->nprey-1;i++) {
    fprintf(fp1, "%.3f\t", param->alpha_prey[i]);  
  }
  fprintf(fp1, "%.3f\n", param->alpha_prey[data->nprey-1]);

  fprintf(fp2, "%d\t", ct+1);
  for(i=0;i<data->nIP-1;i++) {
    fprintf(fp2, "%.3f\t", param->alpha_IP[i]);  
  }
  fprintf(fp2, "%.3f\n", param->alpha_IP[data->nIP-1]);

  fprintf(fp3, "%d\t", ct+1);
  for(i=0;i<data->nprey-1;i++) {
    fprintf(fp3, "%.3f\t", param->mu[i]);  
  }
  fprintf(fp3, "%.3f\n", param->mu[data->nprey-1]);
}



