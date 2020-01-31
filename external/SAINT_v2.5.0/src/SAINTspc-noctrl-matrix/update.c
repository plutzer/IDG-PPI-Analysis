#include "saint.h"

double logZIP(double x, double r, double lambda) {
  double res = 0.0;
  res = -exp(lambda) + x * lambda - gsl_sf_lngamma(x+1.0);
  res += (1.0-r) * exp(res);
  res += (x == 0 ? r : 0.0);
  return res;
}

double loglikAll(DATA *data, PRIOR *prior, PARAM *param) {
  int b, i, j;
  double tmp = 0.0;
  for(i=0;i<param->np;i++) {
    for(b=0;b<data->ninterRow[i];b++) {
      j = data->useRow[i][b];
      if(param->Y[i]==1 && param->Z[i][data->mtu[j]]==1) {
        tmp += ( data->d[i][j] * param->lambda_real[i][j] - exp(param->lambda_real[i][j]) );
      }
      else if(param->Y[i]==1 && param->Z[i][data->mtu[j]]==0) {
        tmp += logZIP(data->d[i][j], param->r0[j], param->lambda_real0[j]); 
	  }
      else {
        if(data->d[i][j] > 0.0) tmp += ( data->d[i][j] * param->lambda_cont[i][j] - exp(param->lambda_cont[i][j]) );
        else tmp += ( data->d[i][j] * param->lambda_real0[j] - exp(param->lambda_real0[j]) );
      }
    }
  }
  return tmp;
}

double loglikAll_realtmp(DATA *data, PRIOR *prior, PARAM *param) {
  int b, i, j;
  double tmp = 0.0;
  for(i=0;i<param->np;i++) {
    for(b=0;b<data->ninterRow[i];b++) {
      j = data->useRow[i][b];
      if(param->Y[i]==1 && param->Z[i][data->mtu[j]]==1) {
        tmp += ( data->d[i][j] * param->lambda_real_tmp[i][j] - exp(param->lambda_real_tmp[i][j]) );
      }
      else if(param->Y[i]==1 && param->Z[i][data->mtu[j]]==0) {
        tmp += logZIP(data->d[i][j], param->r0[j], param->lambda_real0[j]); 
      }
      else {
        if(data->d[i][j] > 0.0) tmp += ( data->d[i][j] * param->lambda_cont[i][j] - exp(param->lambda_cont[i][j]) );
        else tmp += ( data->d[i][j] * param->lambda_real0[j] - exp(param->lambda_real0[j]) );
      }
    }
  }
  return tmp;
}

double loglikAll_conttmp(DATA *data, PRIOR *prior, PARAM *param) {
  int b, i, j;
  double tmp = 0.0;
  for(i=0;i<param->np;i++) {
    for(b=0;b<data->ninterRow[i];b++) {
      j = data->useRow[i][b];
      if(param->Y[i]==1 && param->Z[i][data->mtu[j]]==1) {
        tmp += ( data->d[i][j] * param->lambda_real[i][j] - exp(param->lambda_real[i][j]) );
      }
      else if(param->Y[i]==1 && param->Z[i][data->mtu[j]]==0) {
        tmp += logZIP(data->d[i][j], param->r0[j], param->lambda_real0[j]); 
      }
      else {
        if(data->d[i][j] > 0.0) tmp += ( data->d[i][j] * param->lambda_cont_tmp[i][j] - exp(param->lambda_cont_tmp[i][j]) );
        else tmp += ( data->d[i][j] * param->lambda_real0[j] - exp(param->lambda_real0[j]) );
      }
    }
  }
  return tmp;
}


double loglikRow(DATA *data, PRIOR *prior, PARAM *param, int r) {
  int b, j;
  double tmp = 0.0;
  for(b=0;b<data->ninterRow[r];b++) {
    j = data->useRow[r][b];
    if(param->Y[r]==1 && param->Z[r][data->mtu[j]]==1) {
        tmp += ( data->d[r][j] * param->lambda_real[r][j] - exp(param->lambda_real[r][j]) );
    }
    else if(param->Y[r]==1 && param->Z[r][data->mtu[j]]==0) {
        tmp += logZIP(data->d[r][j], param->r0[j], param->lambda_real0[j]); 
    }
    else {
      if(data->d[r][j] > 0.0) tmp += ( data->d[r][j] * param->lambda_cont[r][j] - exp(param->lambda_cont[r][j]) );
      else tmp += ( data->d[r][j] * param->lambda_real0[j] - exp(param->lambda_real0[j]) );
    }
  }
  return tmp;
}

double loglikCol(DATA *data, PRIOR *prior, PARAM *param, int c) {
  int a,i;
  double tmp = 0.0;
  for(a=0;a<data->ninter[c];a++) {
	i = data->use[c][a];
    if(param->Y[i]==1 && param->Z[i][data->mtu[c]]==1) {
      tmp += ( data->d[i][c] * param->lambda_real[i][c] - exp(param->lambda_real[i][c]) );
    }
    else if(param->Y[i]==1 && param->Z[i][data->mtu[c]]==0) {
        tmp += logZIP(data->d[i][c], param->r0[c], param->lambda_real0[c]); 
    }
    else {
      if(data->d[i][c] > 0.0) tmp += ( data->d[i][c] * param->lambda_cont[i][c] - exp(param->lambda_cont[i][c]) );
      else tmp += ( data->d[i][c] * param->lambda_real0[c] - exp(param->lambda_real0[c]) );
    }
  }
  return tmp;
}

void calcLambdaRealRow(DATA *data, PRIOR *prior, PARAM *param, int r) {
  int b, j, pos;
  double tmp;
  for(b=0;b<data->ninterRow[r];b++) {
    j = data->useRow[r][b];
    tmp = param->beta[0];
    pos = 1;
    if(param->useAbun) {
      tmp += param->beta[pos] * data->preyAbundance[r];
      pos++;
    }
    if(param->useLen) {
      tmp += param->beta[pos] * data->preyLogLength[r];
      pos++;    
    }
    if(param->useCov) {
      tmp += param->beta[pos] * data->baitCoverage[j];
      pos++;    
    }
    tmp += ( param->alpha_prey[r] + param->alpha_bait[j] );
    param->lambda_real[r][j] = tmp;
  }  
}

void calcLambdaRealCol(DATA *data, PRIOR *prior, PARAM *param, int c) {
  int a, i, pos;
  double tmp;
  for(a=0;a<data->ninter[c];a++) {
    i = data->use[c][a];
    tmp = param->beta[0];
    pos = 1;
    if(param->useAbun) {
      tmp += param->beta[pos] * data->preyAbundance[i];
      pos++;
    }
    if(param->useLen) {
      tmp += param->beta[pos] * data->preyLogLength[i];
      pos++;    
    }
    if(param->useCov) {
      tmp += param->beta[pos] * data->baitCoverage[c];
      pos++;    
    }
    tmp += ( param->alpha_prey[i] + param->alpha_bait[c] ) ;
    param->lambda_real[i][c] = tmp;
  }
}

void calcLambdaReal(DATA *data, PRIOR *prior, PARAM *param) {
  int b, i, j, pos;
  double tmp;
  for(i=0;i<param->np;i++) {
    for(b=0;b<data->ninterRow[i];b++) {
      j = data->useRow[i][b];
      tmp = param->beta[0];
      pos = 1;
      if(param->useAbun) {
        tmp += param->beta[pos] * data->preyAbundance[i];
        pos++;
      }
      if(param->useLen) {
        tmp += param->beta[pos] * data->preyLogLength[i];
        pos++;    
      }
      if(param->useCov) {
        tmp += param->beta[pos] * data->baitCoverage[j];
        pos++;    
      }
      
      tmp += ( param->alpha_prey[i] + param->alpha_bait[j] ); 
      param->lambda_real[i][j] = tmp;
    }  
  }
}

void calcLambdaReal_tmp(DATA *data, PRIOR *prior, PARAM *param) {
  int b, i, j, pos;
  double tmp;
  for(i=0;i<param->np;i++) {
    for(b=0;b<data->ninterRow[i];b++) {
      j = data->useRow[i][b];
      tmp = param->beta_tmp[0];
      pos = 1;
      if(param->useAbun) {
        tmp += param->beta_tmp[pos] * data->preyAbundance[i];
        pos++;
      }
      if(param->useLen) {
        tmp += param->beta_tmp[pos] * data->preyLogLength[i];
        pos++;    
      }
      if(param->useCov) {
        tmp += param->beta_tmp[pos] * data->baitCoverage[j];
        pos++;    
      }
      
      tmp += ( param->alpha_prey[i] + param->alpha_bait[j] );
      param->lambda_real_tmp[i][j] = tmp;
    }  
  }
}

void calcLambdaContRow(DATA *data, PRIOR *prior, PARAM *param, int r) {
  int b, j, pos;
  double tmp;
  for(b=0;b<data->ninterRow[r];b++) {
    j = data->useRow[r][b];
    tmp = param->gamma[0];
    pos = 1;
    if(param->useAbun) {
      tmp += param->gamma[pos] * data->preyAbundance[r];
      pos++;
    }
    if(param->useLen) {
      tmp += param->gamma[pos] * data->preyLogLength[r];
      pos++;    
    }
    if(param->useCov) {
      tmp += param->gamma[pos] * data->baitCoverage[j];
      pos++;    
    }
    tmp += ( param->mu_prey[r]) ;      
    param->lambda_cont[r][j] = tmp;
  }  
}

void calcLambdaContCol(DATA *data, PRIOR *prior, PARAM *param, int c) {
  int a, i, pos;
  double tmp;
  for(a=0;a<data->ninter[c];a++) {
    i = data->use[c][a];
    tmp = param->gamma[0];
    pos = 1;
    if(param->useAbun) {
      tmp += param->gamma[pos] * data->preyAbundance[i];
      pos++;
    }
    if(param->useLen) {
      tmp += param->gamma[pos] * data->preyLogLength[i];
      pos++;    
    }
    if(param->useCov) {
      tmp += param->gamma[pos] * data->baitCoverage[c];
      pos++;    
    }
    tmp += (param->mu_prey[i]) ;
    param->lambda_cont[i][c] = tmp;
  }  
}

void calcLambdaCont(DATA *data, PRIOR *prior, PARAM *param) {
  int b, i,j, pos;
  double tmp;
  for(i=0;i<param->np;i++) {
    for(b=0;b<data->ninterRow[i];b++) {
      j = data->useRow[i][b];
      tmp = param->gamma[0];
      pos = 1;
      if(param->useAbun) {
        tmp += param->gamma[pos] * data->preyAbundance[i];
        pos++;
      }
      if(param->useLen) {
        tmp += param->gamma[pos] * data->preyLogLength[i];
        pos++;    
      }
      if(param->useCov) {
        tmp += param->gamma[pos] * data->baitCoverage[j];
        pos++;    
      }
      tmp += (param->mu_prey[i]) ;      
      param->lambda_cont[i][j] = tmp;
    }
  }  
}

void calcLambdaCont_tmp(DATA *data, PRIOR *prior, PARAM *param) {
  int b, i, j, pos;
  double tmp;
  for(i=0;i<param->np;i++) {
    for(b=0;b<data->ninterRow[i];b++) {
	  j = data->useRow[i][b];
      tmp = param->gamma_tmp[0];
      pos = 1;
      if(param->useAbun) {
        tmp += param->gamma_tmp[pos] * data->preyAbundance[i];
        pos++;
      }
      if(param->useLen) {
        tmp += param->gamma_tmp[pos] * data->preyLogLength[i];
        pos++;    
      }
      if(param->useCov) {
        tmp += param->gamma_tmp[pos] * data->baitCoverage[j];
        pos++;    
      }
      tmp += (param->mu_prey[i]) ;      
      param->lambda_cont_tmp[i][j] = tmp;
    }
  }  
}


double dMultGauss(double *x, double *mu, double *Sigma, int length) {
  int i,j,s;
  double temp[length];
  double dmult = 0.0;
  gsl_matrix *inv = gsl_matrix_alloc(length, length);
  gsl_matrix_view m = gsl_matrix_view_array(Sigma, length, length);
  gsl_permutation *p = gsl_permutation_alloc(length);
  gsl_linalg_LU_decomp(&m.matrix, p, &s);
  gsl_linalg_LU_invert(&m.matrix, p, inv);      
  for(j=0;j<length;j++) {
    temp[j] = 0.0;
    for(i=0;i<length;i++) temp[j] += (x[i]-mu[i]) * gsl_matrix_get(inv, i, j) ;
  }
  for(j=0;j<length;j++) dmult += temp[j] * (x[j]-mu[j]);
  gsl_permutation_free(p);
  gsl_matrix_free(inv);
  dmult *= -.5; 
  return dmult;
}

void updateBeta(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r, int k) {
  int i,j;
  double acceptProb; /* , rDraw; */
  
  /* Generate Proposal and 
     Copy Current Values of beta, lambda_real, lik to Temp */
  param->loglik = loglikAll(data, prior, param);
  for(i=0;i<param->nvar;i++) {
    param->beta_tmp[i] = param->beta[i];
  }
  param->beta_tmp[k] += gsl_ran_gaussian(r, _PSD_BETA_ * gsl_ran_gamma(r, 0.5, 0.5));
  calcLambdaReal_tmp(data, prior, param);
  param->loglik_tmp = loglikAll_realtmp(data, prior, param);
  /* fprintf(stderr, "likdiff=%f, oldbeta=%f, newbeta=%f\n", param->loglik_tmp - param->loglik, param->beta_tmp[k], param->beta[k]); */
      
  /* Calculate Accept-Reject Probability and Flip a coin */
  acceptProb = -param->loglik;
  if(param->nvar > 1) acceptProb -= dMultGauss(param->beta, prior->mean_beta, prior->var_beta, param->nvar);
  else acceptProb -= ( - pow(param->beta[0] - prior->mean_beta[0], 2.0) / (2.0 * prior->var_beta[0]) );
  acceptProb += param->loglik_tmp;
  if(param->nvar > 1) acceptProb += dMultGauss(param->beta_tmp, prior->mean_beta, prior->var_beta, param->nvar);
  else acceptProb += ( - pow(param->beta_tmp[0] - prior->mean_beta[0], 2.0) / (2.0 * prior->var_beta[0]) ); 
  acceptProb = exp(acceptProb);
  if(acceptProb > 1.0) acceptProb = 1.0;
  if(gsl_ran_flat(r, 0.0, 1.0) < acceptProb) {
    for(i=0;i<param->nvar;i++) param->beta[i] = param->beta_tmp[i];
    for(i=0;i<param->np;i++) {
      for(j=0;j<param->nb;j++) {
        param->lambda_real[i][j] = param->lambda_real_tmp[i][j];
      }
    }
  }
  else { }  
}


void updateAlphaPrey(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r) {
  int i,j,prey;
  double acceptProb, tmp;
  double mean;
  /* Generate Proposal and 
     Copy Current Values of beta, lambda_real, lik to Temp*/
  for(prey=0;prey<param->np;prey++) {
      param->loglikRow[prey] = loglikRow(data, prior, param, prey);
      for(i=0;i<param->np;i++) param->alpha_prey_tmp[i] = param->alpha_prey[i];
      tmp = gsl_ran_gaussian(r, _PSD_ALPHA_PREY_);
      param->alpha_prey[prey] += tmp;
      param->loglikRow_tmp[prey] = param->loglikRow[prey];

      for(j=0;j<param->nb;j++) {
        param->lambda_real_tmp[prey][j] = param->lambda_real[prey][j];
      }
      calcLambdaRealRow(data, prior, param, prey);
      param->loglikRow[prey] = loglikRow(data, prior, param, prey);
    
      /* Calculate Accept-Reject Probability and Flip a coin */
      acceptProb = param->loglikRow[prey];
      acceptProb += (-.5 * pow(param->alpha_prey[prey], 2.0) / prior->sigmasq_alpha_prey);
      acceptProb -= param->loglikRow_tmp[prey];
      acceptProb -= (-.5 * pow(param->alpha_prey_tmp[prey], 2.0) / prior->sigmasq_alpha_prey);
      acceptProb = exp(acceptProb);
      if(acceptProb > 1.0) acceptProb = 1.0;
      /* tmp = GSL_MAX(param->mu_prey[prey], param->mu_prey_flag[prey]);
      if(param->alpha_prey[prey] < tmp) acceptProb = 0.0; */
      if(gsl_ran_flat(r, 0.0, 1.0) < acceptProb) {
        /* mean = vec_mean(param->alpha_prey, param->np);
	for(i=0;i<param->np;i++) param->alpha_prey[i] -= mean;
	param->beta[0] += mean; */
      }
      else {
        param->alpha_prey[prey] = param->alpha_prey_tmp[prey]; 
        for(j=0;j<param->nb;j++) {
          param->lambda_real[prey][j] = param->lambda_real_tmp[prey][j];
        }
      }
  }
  mean = vec_mean(param->alpha_prey, param->np);
  for(prey=0;prey<param->np;prey++) param->alpha_prey[prey] -= mean;
  param->beta[0] += mean;
}

void updateAlphaDeltaBait(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r) {
  int i,bait;
  double acceptProb;
  double mean;  
  for(bait=0;bait<param->nb;bait++) {
      param->loglikCol[bait] = loglikCol(data, prior, param, bait);
        param->alpha_bait_tmp[bait] = param->alpha_bait[bait];
        param->alpha_bait[bait] += gsl_ran_gaussian(r, _PSD_ALPHA_BAIT_); 
      param->loglikCol_tmp[bait] = param->loglikCol[bait];
      for(i=0;i<param->np;i++) {
        param->lambda_real_tmp[i][bait] = param->lambda_real[i][bait];
      }
      calcLambdaRealCol(data, prior, param, bait);
      param->loglikCol[bait] = loglikCol(data, prior, param, bait);
    
      /* Calculate Accept-Reject Probability and Flip a coin */
      mean = vec_mean(param->alpha_bait, param->nb);
      acceptProb = param->loglikCol[bait];
      acceptProb += (-.5 * pow(param->alpha_bait[bait]-mean, 2.0) / prior->sigmasq_alpha_bait);
      acceptProb -= param->loglikCol_tmp[bait];
      acceptProb -= (-.5 * pow(param->alpha_bait_tmp[bait]-mean, 2.0) / prior->sigmasq_alpha_bait);
      acceptProb = exp(acceptProb);
      if(acceptProb > 1.0) acceptProb = 1.0;

      if(gsl_ran_flat(r, 0.0, 1.0) < acceptProb) { }
      else {
        param->alpha_bait[bait] = param->alpha_bait_tmp[bait];
        for(i=0;i<param->np;i++) {
          param->lambda_real[i][bait] = param->lambda_real_tmp[i][bait];
        }
      }
  }
  mean = vec_mean(param->alpha_bait, param->np);
  for(bait=0;bait<param->nb;bait++) param->alpha_bait[bait] -= mean;
  param->beta[0] += mean;
}

void updateGamma(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r, int k) {
  int i,j;
  double acceptProb; 
  
  param->loglik = loglikAll(data, prior, param);
  for(i=0;i<param->nvar;i++) {
    param->gamma_tmp[i] = param->gamma[i];
  }
  param->gamma_tmp[k] += gsl_ran_gaussian(r, _PSD_GAMMA_ * gsl_ran_gamma(r, 0.5,0.5));
  calcLambdaCont_tmp(data, prior, param);
  param->loglik_tmp = loglikAll_conttmp(data, prior, param);
    
  /* Calculate Accept-Reject Probability and Flip a coin */
  acceptProb = - param->loglik;
  if(param->nvar > 1) acceptProb -= dMultGauss(param->gamma, prior->mean_gamma, prior->var_gamma, param->nvar);
  else acceptProb -= ( - pow(param->gamma[0] - prior->mean_gamma[0], 2.0) / (2.0 * prior->var_gamma[0]) );
  acceptProb += param->loglik_tmp;
  if(param->nvar > 1) acceptProb += dMultGauss(param->gamma_tmp, prior->mean_gamma, prior->var_gamma, param->nvar);
  else acceptProb += ( - pow(param->gamma_tmp[0] - prior->mean_gamma[0], 2.0) / (2.0 * prior->var_gamma[0]) );
  acceptProb = exp(acceptProb); 
  if(acceptProb > 1.0) acceptProb = 1.0; 
  if(gsl_ran_flat(r, 0.0, 1.0) < acceptProb) {
    for(i=0;i<param->nvar;i++) param->gamma[i] = param->gamma_tmp[i];
    for(i=0;i<param->np;i++) {
      for(j=0;j<param->nb;j++) {
        param->lambda_cont[i][j] = param->lambda_cont_tmp[i][j];
      }
    }

  } 
  else {
  
  }       
}


void updateMuPrey(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r) {
  int i,j,prey;
  double acceptProb, tmp;
  double mean;
  /* Generate Proposal and 
     Copy Current Values of beta, lambda_real, lik to Temp*/
  for(prey=0;prey<param->np;prey++) {
      param->loglikRow[prey] = loglikRow(data, prior, param, prey);
      for(i=0;i<param->np;i++) param->mu_prey_tmp[i] = param->mu_prey[i];
      param->loglikRow_tmp[prey] = param->loglikRow[prey];
      for(j=0;j<param->nb;j++) {
        param->lambda_cont_tmp[prey][j] = param->lambda_cont[prey][j];
      }
      tmp = gsl_ran_gaussian(r, _PSD_MU_PREY_);    
      param->mu_prey[prey] += tmp;
      calcLambdaContRow(data, prior, param, prey);    
      param->loglikRow[prey] = loglikRow(data, prior, param, prey);
    
      /* Calculate Accept-Reject Probability and Flip a coin */
      acceptProb = param->loglikRow[prey];
      acceptProb += (-.5 * pow(param->mu_prey[prey], 2.0) / prior->sigmasq_mu_prey);
      acceptProb -= param->loglikRow_tmp[prey];
      acceptProb -= (-.5 * pow(param->mu_prey_tmp[prey], 2.0) / prior->sigmasq_mu_prey);
      acceptProb = exp(acceptProb);
      if(acceptProb > 1.0) acceptProb = 1.0;
      /* if(param->mu_prey[prey] > param->alpha_prey[prey]) acceptProb = 0.0; */
    
      if(gsl_ran_flat(r, 0.0, 1.0) < acceptProb) {  }
      else {
        param->mu_prey[prey] = param->mu_prey_tmp[prey];      
        for(j=0;j<param->nb;j++) {
          param->lambda_cont[prey][j] = param->lambda_cont_tmp[prey][j];
        }
      }
  }      
  mean = vec_mean(param->mu_prey, param->np);
  for(prey=0;prey<param->np;prey++) param->mu_prey[prey] -= mean;
  param->gamma[0] += mean;  
}
 
void updateLambdaReal0(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r) {
  int i,j;
  double shape, rate;
  double a, b;
  int n0;
  
  /* update r0's first */ 
  for(j=0;j<param->nb;j++) {
    n0 = param->np - data->ninter[j];
    a = ((double) n0) / (1.0 + exp(-exp(param->lambda_real0[j])));
    b = ((double) data->ninter[j]) + ((double) n0) *  exp(-exp(param->lambda_real0[j])) / (1 + exp(-exp(param->lambda_real0[j]))) ;
    param->r0[j] = a / (a + b);
  }
	  
  for(j=0;j<param->nb;j++) {
    n0 = param->np - data->ninter[j];
    shape = prior->epsilon_real[j]; 
    rate = prior->kappa_real[j]; 
    for(i=0;i<param->np;i++) {
      if(param->Y[i] == 1 && param->Z[i][data->mtu[j]] == 0 && data->use[j][i]) {
        shape += ((double) data->d[i][j]);
        rate += 1.0;
      }
    }
	rate += ((double) n0) *  exp(-param->lambda_real0[j]) / (1 + exp(-param->lambda_real0[j]));
    param->lambda_real0[j] = log(gsl_ran_gamma(r, shape, 1.0/rate));
  }
}
/*
void updateLambdaReal0(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r) {
  int i,j;
  double shape, rate;
  for(j=0;j<param->nb;j++) {
    shape = prior->epsilon_real[j] + 0.01 * ((double) param->np);
    rate = prior->kappa_real[j] + 0.002 * ((double) param->np);
    for(i=0;i<param->np;i++) {
      if(param->Y[i] == 1 && param->Z[i][data->mtu[j]] == 0) {  
        shape += ((double) data->d[i][j]);
        rate += 1.0;
      }
    }
    param->lambda_real0[j] = log(gsl_ran_gamma(r, shape, 1.0/rate));
  }
}
 */

void updatePriors(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r){
  updateSigmasqAlphaPrey(data, prior, param, r);
  updateSigmasqAlphaBait(data, prior, param, r);
  updateSigmasqMuPrey(data, prior, param, r);
}
 
void sampleY(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r) {
  int b, i,j,k,mfuk, tmpsum0, tmpsum, tmpmaxint, filterout1, filterout2, count;
  double comp1, comp2, ratio, tmp, dotprod;
  double tmp1, tmp2, tmp3, tmpmax;
  double tmp_prob[data->uniqueNum];
  double preal[2];
  int flagged[param->np];
  for(i=0;i<param->np;i++) {
      comp1 = param->pcont[1];
      comp2 = param->pcont[0];
      preal[0] = 0.0; preal[1] = 0.0;
      for(b=0;b<data->ninterRowUnique[i];b++) {
        j = data->useRowUnique[i][b];
        preal[0] += ((double) (1-param->Z[i][j])) ;
        preal[1] += ((double) param->Z[i][j]) ;
      }
      preal[0] += ((double) (data->uniqueNum - data->ninterRowUnique[i]));
      preal[0] /= ((double) data->uniqueNum);
      preal[1] /= ((double) data->uniqueNum);
      for(b=0;b<data->ninterRowUnique[i];b++) {
        j = data->useRowUnique[i][b];
        tmp1 = 0.0; tmp2 = 0.0; tmp3 = 0.0;
        for(k=0;k<data->uniqueSize[j];k++) {
	  mfuk = data->mfu[j][k];
          tmp = param->lambda_real[i][mfuk];
          tmp1 += -exp(param->lambda_real[i][mfuk]) + data->d[i][mfuk] * param->lambda_real[i][mfuk] - gsl_sf_lngamma(data->d[i][mfuk]+1.0);
          tmp2 += logZIP(data->d[i][mfuk], param->r0[mfuk], param->lambda_real0[mfuk]);
	  /* tmp2 += -exp(param->lambda_real0[mfuk]) + data->d[i][mfuk] * param->lambda_real0[mfuk]; */
	  tmp3 += -exp(param->lambda_cont[i][mfuk]) + data->d[i][mfuk] * param->lambda_cont[i][mfuk] - gsl_sf_lngamma(data->d[i][mfuk]+1.0);
	}
	tmpmax = GSL_MAX(tmp1, tmp2);
	tmpmax = GSL_MAX(tmpmax, tmp3);
	tmp1 -= tmpmax; tmp2 -= tmpmax; tmp3 -= tmpmax;
	tmp1 = exp(tmp1); tmp2 = exp(tmp2); tmp3 = exp(tmp3);
        if(tmp1 == 0.0 && tmp2 == 0.0) { tmp_prob[j] = 0.0; }
        else tmp_prob[j] = (preal[1] * tmp1) / ( preal[1] * tmp1 + preal[0] * tmp2);
	comp1 *= (param->Z[i][j] ? tmp1 : tmp2 ); 
        comp2 *= tmp3;   
      }
      if(comp1 == 0.0 && comp2 == 0.0) ratio = 0.0;
      else ratio = comp1 / (comp1 + comp2);
      if(isnan(ratio)) ratio = 0.0;
      if(!finite(ratio)) ratio = 1.0;
      if(gsl_ran_flat(r, 0.0, 1.0) < ratio) param->Y[i] = 1;
      else param->Y[i] = 0;
      /* if(vec_max(param->lambda_cont[i], param->nb) <= 2.0 && vec_max(param->lambda_real[i], param->nb) <= 2.0) param->Y[i] = 1; */
      tmpsum = 0; tmpsum0 = 0;
      for(b=0;b<data->ninterRowUnique[i];b++) {
        j = data->useRowUnique[i][b];
        tmpmaxint = 0;
        for(k=0;k<data->uniqueSize[j];k++) {
          if(((int) data->d[i][data->mfu[j][k]]) > tmpmaxint) tmpmaxint = ((int) data->d[i][data->mfu[j][k]]);
        }
        if(tmpmaxint >= 1) tmpsum0++; 
        if(param->Z[i][j]) tmpsum++;
      }

      filterout1 = ( (int) ( ((double) tmpsum) / ((double) data->uniqueNum)  >  param->ff_prop) );
      filterout2 = ( (int) ( ((double) tmpsum0) / ((double) data->uniqueNum)  <=  param->ff_prop) );

      if(filterout2) {
        param->Y[i] = 1;
        param->flagged[i] += 1.0;
      }

      flagged[i] = 0;
      if(filterout1) {
        flagged[i] = 1;
        param->Y[i] = 0;
      } 

      if(data->override[i] == 'C') param->Y[i] = 0;
      else if(data->override[i] == 'R') param->Y[i] = 1;
      else {  }
  }

  /* Calculate Frequency Distribution -- appearCont */
  count = 0;
  for(i=0;i<param->np;i++) {
    if(param->Y[i]==0) count++;
  }
  for(j=0;j<param->nb;j++) {
    param->appearCont[j] = 0.0;
  }
  for(j=0;j<param->nb;j++) {
    for(i=0;i<param->np;i++) {
      if(param->Y[i]==0 && data->d[i][j] > 0.0) param->appearCont[j] += (1.0 / ((double) count));
    }
  }

  /* Flag unflagged ones with high nonzero counts > ff_prop */  
  for(i=0;i<param->np;i++) {
    if(param->Y[i] && flagged[i] == 0) {}
    else {
      dotprod = 0.0;
      for(j=0;j<param->nb;j++) {
        if(data->d[i][j] > 0.0) dotprod += param->appearCont[j] / ((double) param->nb);
      }
      if(dotprod >= param->ff_prop) {
        param->Y[i] = 0;
      }
    }
  }
}

void sampleZ(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r) {
  int b, i,j,k,mfuk,anypos;
  double p1, p2;
  double tmp1, tmp2, tmp, tmpmax, prob, tp;
  for(i=0;i<param->np;i++) {
    for(b=0;b<data->ninterRowUnique[i];b++) {
      j = data->useRowUnique[i][b];
      if(data->uniqueSize[j] > 1) {
        anypos = 0;
        for(k=0;k<data->uniqueSize[j];k++) {
          if(data->d[i][data->mfu[j][k]] > 0.0) anypos = 1;
        }
        if(anypos == 0) {
          for(k=0;k<data->uniqueSize[j];k++) param->iZ[i][data->mfu[j][k]] = 0;
          param->Z[i][j] = 0;
        }
        else {
          p1 = p2 = 1.0;
          for(k=0;k<data->uniqueSize[j];k++) {
            mfuk = data->mfu[j][k];	
            tmp = data->d[i][mfuk] <= exp(param->lambda_real[i][mfuk]) ? param->lambda_real[i][mfuk] : log(data->d[i][mfuk]);  /* used to be masked -- nov-17-08 */
            tmp = param->lambda_real[i][mfuk]; 
            tmp1 = -exp(tmp) + data->d[i][mfuk] * tmp - gsl_sf_lngamma(data->d[i][mfuk]+1.0);
            tmp2 = logZIP(data->d[i][mfuk], param->r0[mfuk], param->lambda_real0[mfuk]);
	    tmpmax = GSL_MAX(tmp1, tmp2);
	    tmp1 -= tmpmax; tmp2 -= tmpmax;
	    p1 *= exp(tmp1); p2 *= exp(tmp2);
            /* make calls and calculate distance within replicate */
            tp = param->preal[1] * exp(tmp1) / (param->preal[1] * exp(tmp1) + param->preal[0] * exp(tmp2));
            param->iZ[i][mfuk] = gsl_ran_flat(r, 0.0, 1.0) < tp ? 1 : 0;
            if(data->d[i][mfuk] == 0.0) param->iZ[i][mfuk] = 0;
	  }
          p1 *= param->preal[1]; p2 *= param->preal[0];
          prob = p1 / (p1 + p2);
          if(gsl_ran_flat(r, 0.0, 1.0) < prob) param->Z[i][j] = 1;
          else param->Z[i][j] = 0;
        }
      }
      else {
        anypos = 0;
        if(data->d[i][data->mfu[j][0]] > 0.0) anypos = 1;
        if(anypos == 0) {
          param->iZ[i][data->mfu[j][0]] = 0;
          param->Z[i][j] = 0;
        }
        else {
          p1 = p2 = 1.0;
          mfuk = data->mfu[j][0];	
          /* tmp = data->d[i][mfuk] <= exp(param->lambda_real[i][mfuk]) ? param->lambda_real[i][mfuk] : log(data->d[i][mfuk]); */ /* used to be masked -- nov-17-08 */
          tmp = param->lambda_real[i][mfuk]; 
          tmp1 = -exp(tmp) + data->d[i][mfuk] * tmp - gsl_sf_lngamma(data->d[i][mfuk]+1.0);
          tmp2 = logZIP(data->d[i][mfuk], param->r0[mfuk], param->lambda_real0[mfuk]);

	  tmpmax = GSL_MAX(tmp1, tmp2);
	  tmp1 -= tmpmax; tmp2 -= tmpmax;
	  p1 *= exp(tmp1); p2 *= exp(tmp2);
            /* make calls and calculate distance within replicate */
          p1 *= param->preal[1]; p2 *= param->preal[0];
          prob = p1 / (p1 + p2);
          if(gsl_ran_flat(r, 0.0, 1.0) < prob) {
            param->Z[i][j] = 1;
            param->iZ[i][mfuk] = 1;
          }
          else {
            param->Z[i][j] = 0;
            param->iZ[i][mfuk] = 0;
          }
        }
      }
    }
  }
}

double calc_dist(int calls[], int num) {
  int i,j,count;
  double dist = 0.0;
  count = 0;
  for(i=0;i<(num-1);i++) {
    for(j=(i+1);j<num;j++) {
      dist += (calls[i] == calls[j] ? 0.0 : 1.0);
      count++; 
    }
  }
  dist /= ((double) count);
  return dist;
}

void updateSigmasqAlphaPrey(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r) {
  int i;
  double shape, rate;
  shape = 1.0;
  rate = 1.0;
  for(i=0;i<param->np;i++) {
    if(param->Y[i]) { 
      shape += 0.5;
      rate += 0.5 * pow(param->alpha_prey[i], 2.0);
    }
  }
  prior->sigmasq_alpha_prey = 1.0 / gsl_ran_gamma(r, shape, 1.0/rate); 
}

void updateSigmasqAlphaBait(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r) {
  int j;
  double shape, rate, mean;
  shape = 1.0;
  rate = 1.0;
  mean = vec_mean(param->alpha_bait, param->nb);
  for(j=0;j<param->nb;j++) {
    shape += 0.5;
    rate += 0.5 * pow(param->alpha_bait[j]-mean, 2.0);
  }
  prior->sigmasq_alpha_bait = 1.0 / gsl_ran_gamma(r, shape, 1.0/rate); 
}

void updateSigmasqMuPrey(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r) {
  int i;
  double shape, rate;
  shape = 1.0;
  rate = 1.0;
  for(i=0;i<param->np;i++) {
    if(param->Y[i]==0) { 
      shape += 0.5;
      rate += 0.5 * pow(param->mu_prey[i], 2.0);
    }
  }
  prior->sigmasq_mu_prey = 1.0 / gsl_ran_gamma(r, shape, 1.0/rate); 
}

void sampleP(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r) {
  int b, i, j;
  int n[2];
  double alpha[2];
  for(i=0;i<2;i++) n[i] = 0;
  for(i=0;i<param->np;i++) (n[param->Y[i]])++;
  for(i=0;i<2;i++) alpha[i] = prior->acont[i] + ((double) n[i]);
  gsl_ran_dirichlet(r, 2, alpha, param->pcont);

  for(i=0;i<2;i++) n[i] = 0;
  for(i=0;i<param->np;i++) {
    if(param->Y[i]) {
      for(b=0;b<data->ninterRowUnique[i];b++) {
        j = data->useRowUnique[i][b];
        (n[param->Z[i][j]])++;
      }
      n[0] += data->uniqueNum - data->ninterRowUnique[i];
    }
  }
  for(i=0;i<2;i++) alpha[i] = prior->areal[i] + ((double) n[i]);
  gsl_ran_dirichlet(r, 2, alpha, param->preal);
  /* fprintf(stderr, "C-%.3f R-%.3f NR-%.3f\n", param->pcont[0], param->preal[1], param->preal[0]); */
}

/**************************************/
/*** Metropolis-Hastings with Gibbs ***/
/**************************************/
void mhgibbs(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r) {
  int k;
  /* fprintf(stderr, "beta:");
  for(k=0;k<param->nvar;k++) fprintf(stderr, "%.5f ", param->beta[k]);
  fprintf(stderr, "gamma:");
  for(k=0;k<param->nvar;k++) fprintf(stderr, "%.5f ", param->gamma[k]);
  fprintf(stderr, "\n"); */
  for(k=0;k<param->nvar;k++) updateBeta(data, prior, param, r, k);
  updateAlphaPrey(data, prior, param, r);
  updateAlphaDeltaBait(data, prior, param, r);
  /* calcLambdaReal(data, prior, param); */
  updateLambdaReal0(data, prior, param, r);
  for(k=0;k<param->nvar;k++) updateGamma(data, prior, param, r, k);
  updateMuPrey(data, prior, param, r);
  updatePriors(data, prior, param, r);
  sampleY(data, prior, param, r);
  sampleZ(data, prior, param, r);    
  sampleP(data, prior, param, r);    
}

