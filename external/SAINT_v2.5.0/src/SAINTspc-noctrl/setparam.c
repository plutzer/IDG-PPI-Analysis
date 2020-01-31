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
  int isCtrl, isReverse, isLarge;
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
      if(param->lambda_true[id] < param->lambda_false[id] || exp(param->lambda_false[id]) >= 5.0) {
        isReverse = 1;
        break;
      }
    }
    isLarge = 0;
    for(j=0;j<data->n_u2a[i];j++) {
      id = data->u2a[i][j];
      if(data->d[id] >= 5.0) {
        isLarge = 1;
        break;
      }
    }

    if(isLarge && !isCtrl) {
      param->Z[i] = 1;
      for(j=0;j<data->n_u2a[i];j++) {
        id = data->u2a[i][j];
        if(data->d[id] > 0.0) param->iZ[id] = 1;
      }
    }
    else if(isCtrl || isReverse) {
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
        posi = log_poisson_g_prop(tmp, exp(param->lambda_true[id]), param->eta[data->i2p[id]]);
        /* tmp = data->d[id] < exp(param->lambda_false[id]) ? exp(param->lambda_false[id]) : data->d[id]; */
        negi = log_poisson_g_prop(tmp, exp(param->lambda_false[id]), param->eta0[data->i2p[id]]);
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
        maxl = pos > neg ? pos : neg;
        pos -= maxl;
        neg -= maxl;
        prob = param->ptrue * exp(pos) / (param->ptrue * exp(pos) + (1.0-param->ptrue) * exp(neg));
        param->Z[i] = gsl_ran_flat(r,0.0,1.0) <= prob ? 1 : 0;
      }
    }
  }
}


void initialize_param(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r) {
  int i;

  param->beta0 = 0.0;
  param->betac = 0.0;
  for(i=0;i<data->nprey;i++) param->alpha_prey[i] = prior->theta_alpha_prey[prior->w_alpha_prey[i]];
  for(i=0;i<data->nIP;i++) {
    param->alpha_IP[i] = 0.0;
    /* if(data->ctrl[i] == 0) param->alpha_IP[i] = prior->theta_alpha_IP[prior->w_alpha_IP[i]];
    else param->alpha_IP[i] = 0.0; */
  }
  for(i=0;i<data->nprey;i++) param->mu[i] = prior->theta_mu[prior->w_mu[i]];
  if(param->modelvar) {
    for(i=0;i<data->nprey;i++) param->eta[i] = prior->theta_eta[prior->w_eta[i]];
    for(i=0;i<data->nprey;i++) param->eta0[i] = prior->theta_eta0[prior->w_eta0[i]];
  }
  else {
    for(i=0;i<data->nprey;i++) param->eta[i] = 1.0;
    for(i=0;i<data->nprey;i++) param->eta0[i] = 1.0;
  }

  compute_lambda_all(param, prior, data);
  set_Z(param, prior, data, r);
  /* param->loglikTotal = loglik_all(param, prior, data); */
  param->ptrue = 0.1;
  param->ptrue_tmp = 0.1;
}

void set_param(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r) {
  memory_param(param, prior, data);
  initialize_param(param, prior, data, r);
}



