#include "saint.h"

/**************************************************************/
/*             initializing the model parameters              */
/**************************************************************/

void memory_param(PARAM *param, PRIOR *prior, DATA *data) {
  assert(param->loglik_prey = (float *) calloc(data->nprey, sizeof(float)));
  assert(param->loglik_IP = (float *) calloc(data->nIP, sizeof(float)));
  assert(param->alpha_prey = (float *) calloc(data->nprey, sizeof(float)));
  assert(param->alpha_IP = (float *) calloc(data->nIP, sizeof(float)));
  assert(param->mu = (float *) calloc(data->nprey, sizeof(float)));
  assert(param->eta = (float *) calloc(data->nprey, sizeof(float)));
  assert(param->eta0 = (float *) calloc(data->nprey, sizeof(float)));
  assert(param->iZ = (int *) calloc(data->ninter, sizeof(int)));
  assert(param->Z = (int *) calloc(data->nuinter, sizeof(int)));
  assert(param->lambda_true = (float *) calloc(data->ninter, sizeof(float)));
  assert(param->lambda_false = (float *) calloc(data->ninter, sizeof(float)));
  assert(param->lambda_true_tmp = (float *) calloc(data->ninter, sizeof(float)));
  assert(param->lambda_false_tmp = (float *) calloc(data->ninter, sizeof(float)));
}


void set_Z(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r) {
  /* Z and iZ */
  int i,j,id;
  int indiv, total;
  int isCtrl, isReverse;
  float prob, maxl;
  float posi, negi;
  float pos, neg, tmp;
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
      if(data->d[id] < param->lambda_false[id]) {
        isReverse = 0;
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
        posi = log_gaussian(tmp, param->lambda_true[id], param->eta[data->i2p[id]]);
        negi = log_gaussian(tmp, param->lambda_false[id], param->eta0[data->i2p[id]]);
        pos += posi;
        neg += negi;
        maxl = posi > negi ? posi : negi;
        posi -= maxl;
        negi -= maxl;
        prob = param->ptrue * exp(posi) / (param->ptrue * exp(posi) + (1.0-param->ptrue) * exp(negi));
        param->iZ[id] = gsl_ran_flat(r,0.0,1.0) <= prob ? 1 : 0;
      }
      /* Z */
      if(data->n_u2a[i] == 1) {
        id = data->u2a[i][0];
        param->Z[i] = param->iZ[id];	
      }
      else {
        /* maxl = pos > neg ? pos : neg;
        pos -= maxl;
        neg -= maxl;
        prob = param->ptrue * exp(pos) / (param->ptrue * exp(pos) + (1.0-param->ptrue) * exp(neg));
        param->Z[i] = gsl_ran_flat(r,0.0,1.0) <= prob ? 1 : 0; */
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

float vec_obs_mean(PARAM *param, DATA *data) {
  int i,ct;
  float m = 0.0;
  ct = 0;
  for(i=0;i<data->ninter;i++) {
    if(data->miss[i] == 0) {
      m += data->d[i];
      ct++;
    }
  }
  m /= ((float) ct);
  return m;
}

float vec_obs_var(PARAM *param, DATA *data, float m) {
  int i,ct;
  float v = 0.0;
  ct = 0;
  for(i=0;i<data->ninter;i++) {
    if(data->miss[i] == 0) {
      v += pow(data->d[i] - m, 2.0);
      ct++;
    }
  }
  v /= ((float) (ct-1));
  return v;
}

void initialize_param(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r) {
  int i,j,k;
  int isSource[data->nprey];  /* whether the prey is used for missing distribution or not */
  float tmp_mean, tmp_var;
  
  param->beta0 = 0.0;
  param->betac = 0.0;
  for(i=0;i<data->nprey;i++) param->alpha_prey[i] = prior->theta_alpha_prey[prior->w_alpha_prey[i]];
  for(i=0;i<data->nIP;i++) {
    param->alpha_IP[i] = 0.0;
    /* if(data->ctrl[i] == 0) param->alpha_IP[i] = prior->theta_alpha_IP[prior->w_alpha_IP[i]];
    else param->alpha_IP[i] = 0.0; */
  }
  for(i=0;i<data->nprey;i++) param->mu[i] = prior->theta_mu[prior->w_mu[i]];
  for(i=0;i<data->nprey;i++) param->eta[i] = prior->theta_eta[prior->w_eta[i]];
  for(i=0;i<data->nprey;i++) param->eta0[i] = prior->theta_eta0[prior->w_eta0[i]];

  compute_lambda_all(param, prior, data);
  set_Z(param, prior, data, r);
  /* param->loglikTotal = loglik_all(param, prior, data); */
  param->ptrue = 0.2;
  param->ptrue_tmp = 0.2;

  param->missMean = 0.0;
  j = 0;
  for(i=0;i<data->ninter;i++) {
    if(data->ctrl[data->i2IP[i]] && data->miss[i] == 0) {
      param->missMean += data->d[i];
      j++;
    }
  }
  param->missMean /= ((double) j);
  
  param->missVar = 0.0;
  for(i=0;i<data->ninter;i++) {
    if(data->ctrl[data->i2IP[i]] && data->miss[i] == 0) {
      param->missVar += pow(data->d[i] - param->missMean, 2.0);
    }
  }
  param->missVar /= ((double) (j-1));
  data->dvar = param->missVar;

  /*************************************************************/
  /* this is where the missing data distribution is determined */
  /* If the number of control is 1, then do the same as now    */
  /* Otherwise,                                                */
  /* 1. Identify prey proteins with any missing data           */
  /* 2. Use the observed intensities from these proteins       */
  /*    to estimate mean and variance                          */
  /*************************************************************/

  if(data->nctrl <= 1) {
    param->missMean = data->dmin__ctrl - 0.0 * sqrt(param->missVar);
    param->missVar *= 0.5;
  }
  else {
     for(i=0;i<data->nprey;i++) {
        isSource[i] = 0;
     }
     for(i=0;i<data->ninter;i++) {
        if(data->miss[i] == 1 && data->ctrl[data->i2IP[i]]) isSource[data->i2p[i]] = 1;
     }
     tmp_mean = 0.0;
     tmp_var = 0.0;
     k = 0;
     for(i=0;i<data->ninter;i++) {
        if(isSource[data->i2p[i]] == 1 && data->miss[i] == 0 && data->ctrl[data->i2IP[i]]) {
           tmp_mean += data->d[i];
           k++;
        }
     }
     tmp_mean /= ((float) k);
     for(i=0;i<data->ninter;i++) {
        if(isSource[data->i2p[i]] == 1 && data->miss[i] == 0 && data->ctrl[data->i2IP[i]]) {
           tmp_var += pow(data->d[i] - tmp_mean, 2.0);
           k++;
        }
     }
     tmp_var /= ((float) (k-1));
     param->missMean = tmp_mean;    
     param->missVar = tmp_var;

  }


  for(i=0;i<data->ninter;i++) {
    if(data->miss[i]) data->d[i] = gsl_ran_gaussian(r, sqrt(param->missVar)) + param->missMean;
  }
  //fprintf(stderr, "%.3f\n", param->missMean);
  //fprintf(stderr, "%.3f\n", param->missVar);
}

void set_param(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r) {
  memory_param(param, prior, data);
  initialize_param(param, prior, data, r);
}



