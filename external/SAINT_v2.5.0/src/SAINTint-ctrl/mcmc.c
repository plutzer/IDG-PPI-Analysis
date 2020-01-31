#include "saint.h"

float gaussian(float x, float mu, float var) {
  float res = exp( - .5 * pow(x-mu,2.0) / var ) /  sqrt(2.0 * M_PI * var);
  return res;
}

float log_gaussian(float x, float mu, float var) {
  float res;
  if(x == GSL_NEGINF) {
    res = gsl_cdf_gaussian_P(dmin - log(2.0) - mu, sqrt(var));
    res = log(GSL_MAX(res, _tiny_));
  }  
  else res = - .5 * pow(x-mu,2.0) / var - .5 * log(2.0 * M_PI * var);
  return res;
}



float log_inv_gamma(float x, float sh, float sc) {
  float res = 0.0;
  res = sh * log(sc) - gsl_sf_lngamma(sh) - (sh + 1.0) * log(x) - sc / x;
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
  int isCtrl, isReverse;
  float prob, maxl;
  float posi, negi;
  float pos, neg, tmp, tmp_lambda, tmp_neg, scale;
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
      cond1 = data->d[id] < GSL_MAX(param->lambda_false[id], data->ctrlavg[data->i2p[id]]) ? 1 : 0;
      if(cond1) {
        isReverse = 1;
        break;
      }
    }

    if(isCtrl || isReverse) {
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
        tmp = data->d[id];
        tmp_lambda = (param->lambda_true[id]);
        scale = 1.0;
        if(tmp > tmp_lambda && tmp_lambda > param->lambda_false[id]) tmp = tmp_lambda;
        if(data->d[id] < data->ctrlavg[data->i2p[id]] + log(10.0) && data->ctrl_obs[data->i2p[id]] > 0 ) {
          tmp_lambda = data->ctrlavg[data->i2p[id]] + 2.0 * (data->ctrlavg[data->i2p[id]] + log(10.0) - data->d[id]);
          scale = 1.0;
        }
        posi = log_gaussian(tmp, tmp_lambda, param->eta[data->i2p[id]]);
        tmp_neg = (param->lambda_false[id]);
        tmp = GSL_MAX(data->d[id], tmp_neg);
        negi = log_gaussian(tmp, tmp_neg, param->eta0[data->i2p[id]] * scale);
        pos += posi;
        neg += negi;
        maxl = posi > negi ? posi : negi;
        posi -= maxl;
        negi -= maxl;
        prob = param->ptrue * exp(posi) / (param->ptrue * exp(posi) + (1.0-param->ptrue) * exp(negi));
        param->iZ[id] = gsl_ran_flat(r,0.0,1.0) <= prob ? 1 : 0;
        if(data->miss[id]) param->iZ[id] = 0;
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
        param->Z[i] = pos;
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


void updateMiss(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r) {
  int i,p,ctrl_status;
  float new_miss, mhratio, maxval, tmpeta, tmplambda;
  for(i=0;i<data->ninter;i++) {
    if(data->miss[i]) {
      ctrl_status = data->ctrl[data->i2IP[i]];
      p = data->i2p[i];
      maxval = data->dmin_ctrl[p];
      new_miss = data->d[i] + gsl_ran_gaussian(r, 1.0);
      if(ctrl_status) {
        tmplambda = param->lambda_false[i];
        tmpeta = param->eta0[p];
        mhratio = 0.0;
        if(data->ctrl_obs[p] > 0) {
          // mhratio += log_gaussian(new_miss, tmplambda, tmpeta);
          // mhratio -= log_gaussian(data->d[i], tmplambda, tmpeta);
        }
        mhratio += log_gaussian(new_miss, param->missMean, param->missVar);
        mhratio -= log_gaussian(data->d[i], param->missMean, param->missVar);
        mhratio = GSL_MIN(1.0, exp(mhratio));
        if(gsl_ran_flat(r,0.0,1.0) < mhratio) {
          data->d[i] = new_miss;
        }
      }
    }
  }
}

/**************************************/
/*** Metropolis-Hastings with Gibbs ***/
/**************************************/
void mhgibbs(PARAM *param, PRIOR *prior, DATA *data, SUMMARY *summary, const gsl_rng *r, int updateSum) {

  // updateMiss(param, prior, data, r);
  // set_ctrlavg(data);
  if(gsl_ran_flat(r,0.0,1.0) <= 0.2) sampleBeta0(param, prior, data, r);
  if(gsl_ran_flat(r,0.0,1.0) <= 0.2) sampleBetac(param, prior, data, r);
  DP_alpha_prey(param, prior, data, r);
  // DP_alpha_IP(param, prior, data, r);
  DP_mu(param, prior, data, r); 
  if(gsl_ran_flat(r,0.0,1.0) <= 0.2) DP_eta(param, prior, data, r); 
  if(gsl_ran_flat(r,0.0,1.0) <= 0.2) DP_eta0(param, prior, data, r); 
  sampleZ(param, prior, data, r);
  // fprintf(stderr, "%.3f\t(%.3f\t%.3f)\t(%.3f\t%.3f)\n", data->d[186], param->lambda_true[186], param->eta[data->i2p[186]], param->lambda_false[186], param->eta0[data->i2p[186]]);
  // fprintf(stderr, "%.3f\n", param->eta0[3]);
  compute_lambda_all(param, prior, data);
  if(gsl_ran_flat(r,0.0,1.0) <= 0.2) sampleProportion(param, prior, data, r);
}



void updateMissMean(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r) {
  int i;
  float new_miss, mhratio;
  new_miss = param->missMean + gsl_ran_gaussian(r, 0.1);

  mhratio = 0.0;
  for(i=0;i<data->ninter;i++) {
    if(data->miss[i]) {
      mhratio += log_gaussian(data->d[i], new_miss, param->missVar);
      mhratio -= log_gaussian(data->d[i], param->missMean, param->missVar);
    }
  }
  mhratio += log_gaussian(new_miss, 0.0, 10.0);
  mhratio -= log_gaussian(param->missMean, 0.0, 10.0);
  mhratio = GSL_MIN(1.0, exp(mhratio));
  if(gsl_ran_flat(r,0.0,1.0) < mhratio) {
    param->missMean = new_miss;
  }
}

void updateMissVar(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r) {
  int i;
  float new_miss, mhratio;
  new_miss = param->missVar + gsl_ran_gaussian(r, 0.1);
 
  if(new_miss > 0.0) {
    mhratio = 0.0;
    for(i=0;i<data->ninter;i++) {
      if(data->miss[i]) {
          mhratio += log_gaussian(data->d[i], param->missMean, new_miss);
          mhratio -= log_gaussian(data->d[i], param->missMean, param->missVar);
      }
    }
    mhratio -= log(gsl_ran_gamma_pdf(new_miss, 1.0, 1.0));
    mhratio += log(gsl_ran_gamma_pdf(param->missVar, 1.0, 1.0));
    mhratio = GSL_MIN(1.0, exp(mhratio));
    if(gsl_ran_flat(r,0.0,1.0) < mhratio) {
      param->missVar = new_miss;
    }
  }
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



