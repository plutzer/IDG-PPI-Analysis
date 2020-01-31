#include "saint.h"

void init_data(DATA *data, int *p, int *q) {
  int i;
  data->nrow = *p;
  data->ncol = *q;
  assert(data->d = (double **) calloc(*p, sizeof(double *)));
  for(i=0; i<*p; i++) {
    assert(data->d[i] = (double *) calloc(*q, sizeof(double)));
  }
  assert(data->maxRow = (double *) calloc(*p, sizeof(double)));
  assert(data->preyORF = (char **) calloc(*p, sizeof(char *)));
  assert(data->override = (char *) calloc(*p, sizeof(char)));

  for(i=0; i<*p; i++) {
    assert(data->preyORF[i] = (char *) calloc(_MAX_NAME_, sizeof(char)));
  }
  assert(data->prey = (char **) calloc(*p, sizeof(char *))); 
  for(i=0; i<*p; i++) {
    assert(data->prey[i] = (char *) calloc(_MAX_NAME_, sizeof(char))); 
  }
  assert(data->experiment = (char **) calloc(*q, sizeof(char *)));
  for(i=0; i<*q; i++) {
    assert(data->experiment[i] = (char *) calloc(_MAX_NAME_, sizeof(char)));
  }
  assert(data->bait = (char **) calloc(*q, sizeof(char *))); 
  for(i=0; i<*q; i++) {
    assert(data->bait[i] = (char *) calloc(_MAX_NAME_, sizeof(char))); 
  }  
  assert(data->unique = (char **) calloc(*q, sizeof(char *))); 
  for(i=0; i<*q; i++) {
    assert(data->unique[i] = (char *) calloc(_MAX_NAME_, sizeof(char))); 
  }
  assert(data->uniqueSize = (int *) calloc(*q, sizeof(int)));
  assert(data->mtu = (int *) calloc(*q, sizeof(int)));
  assert(data->mfu = (int **) calloc(*q, sizeof(int *)));
  for(i=0; i<*q; i++) {
    assert(data->mfu[i] = (int *) calloc(_MAX_REPLICA_, sizeof(int)));
  }
  assert(data->baitCoverage = (double *) calloc(*q, sizeof(double)));
  assert(data->preyAbundance = (double *) calloc(*p, sizeof(double)));
  assert(data->preyLogLength = (double *) calloc(*p, sizeof(double)));
}

void free_data(DATA *data) {
  int i;
  for(i=0; i<data->nrow; i++) {
    free(data->d[i]);
    free(data->prey[i]); 
    free(data->preyORF[i]);
  }  
  free(data->d);
  free(data->maxRow);
  free(data->prey); 
  free(data->preyORF);
  for(i=0; i<data->ncol; i++) {
    free(data->bait[i]);
    free(data->experiment[i]);
  }
  for(i=0; i<data->ncol; i++) {
    free(data->unique[i]); 
    free(data->mfu[i]);
  }
  free(data->bait);
  free(data->experiment);
  free(data->unique);
  free(data->uniqueSize);
  free(data->mtu);
  free(data->mfu);
  free(data->baitCoverage);
  free(data->preyAbundance);
  free(data->preyLogLength);
}

void normalizeCoverage(double *cover, int *len) {
  int j;
  double tmp[*len];
  double max, med;
  for(j=0;j<*len;j++) cover[j] = cover[j] + 1.0;
  max = vec_max(cover, *len);
  for(j=0;j<*len;j++) cover[j] /= max;
  for(j=0;j<*len;j++) cover[j] *= 0.99;
  for(j=0;j<*len;j++) tmp[j] = gsl_cdf_gaussian_Pinv(cover[j], 1.0);
  med = vec_med(tmp, *len);
  for(j=0;j<*len;j++) cover[j] = tmp[j] - med;
}

void read_data(FILE *fp, DATA *data, int *p, int *q) {
  int i,j,k,cur,curU,exist;
  char buf[_MAX_BUF_];
  init_data(data, p, q); 
  double med;
  
  /* Read first two lines of info */
  for(j=0;j<4;j++) fscanf(fp,"%s",buf);
  for(j=0;j<data->ncol;j++) {
    fscanf(fp,"%s",buf);
    strcpy(data->experiment[j], buf);
  }

  for(j=0;j<4;j++) fscanf(fp,"%s",buf);
  for(j=0;j<data->ncol;j++) {
    fscanf(fp,"%s",buf);
    strcpy(data->bait[j], buf);
  }

  for(j=0;j<4;j++) fscanf(fp,"%s",buf);
  for(j=0;j<data->ncol;j++) {
    fscanf(fp,"%s",buf);
    data->baitCoverage[j] = atof(buf);
  }
  /* normalizeCoverage(data->baitCoverage, q);  */
  
  /* Read each row 
     First three elements are: 
        preyName, preyAbundance, preyLogLength */

  for(i=0;i<*p;i++) {
    fscanf(fp,"%s",buf);
    strcpy(data->prey[i], buf); 
    fscanf(fp,"%s",buf);
    data->preyAbundance[i] = atof(buf); 
    fscanf(fp,"%s",buf);
    data->preyLogLength[i] = atof(buf);
    fscanf(fp,"%s",buf);
    data->override[i] = buf[0];
    for(j=0;j<data->ncol;j++) {
      fscanf(fp,"%s",buf);
      data->d[i][j] = atof(buf);
    }
  }
  
  /* Unique Bait-Experiments Identification */
  cur = 0;
  strcpy(data->unique[cur], data->bait[cur]);
  cur++;
  for(j=1;j<*q;j++) {
    exist = 0;
    for(i=0;i<cur;i++) {
      if(strcmp(data->unique[i], data->bait[j]) == 0) exist = 1; 
    }
    if(exist == 0) {
      strcpy(data->unique[cur], data->bait[j]);
      cur++;
    }
  }
  data->uniqueNum = cur;

  for(i=0;i<data->uniqueNum;i++) {
    cur = 0;
    for(j=0;j<*q;j++) {
      if(strcmp(data->unique[i], data->bait[j]) == 0) {
        data->mtu[j] = i;
        data->mfu[i][cur] = j;
        cur++;
      }
    }
    data->uniqueSize[i] = cur;
  }

  /* use and useUnique indices */
  /* by column */
  assert(data->ninterUnique = (int *) calloc(data->uniqueNum, sizeof(int)));
  assert(data->useUnique = (int **) calloc(data->uniqueNum, sizeof(int *)));
  assert(data->ninter = (int *) calloc(*q, sizeof(int)));
  assert(data->use = (int **) calloc(*q, sizeof(int *)));
  /* by row */
  assert(data->ninterRowUnique = (int *) calloc(*p, sizeof(int)));
  assert(data->useRowUnique = (int **) calloc(*p, sizeof(int *)));
  assert(data->ninterRow = (int *) calloc(*p, sizeof(int)));
  assert(data->useRow = (int **) calloc(*p, sizeof(int *)));

  /* column */
  for(j=0;j<data->uniqueNum;j++) { /* Count first */ 
    data->ninterUnique[j] = 0;
    for(i=0;i<*p;i++) {
      exist = 0;
      for(k=0;k<data->uniqueSize[j];k++) {
        if(data->d[i][data->mfu[j][k]] > 0) exist = 1;
      }
      if(exist) { 
        (data->ninterUnique[j])++;
      }
    }
  }
  for(j=0;j<data->uniqueNum;j++) {
    for(k=0;k<data->uniqueSize[j];k++) data->ninter[data->mfu[j][k]] = data->ninterUnique[j];
  }
  
  for(k=0;k<*q;k++) {
     assert(data->use[k] = (int *) calloc(data->ninterUnique[data->mtu[k]], sizeof(int)));  
  }
  for(j=0;j<data->uniqueNum;j++) {
    assert(data->useUnique[j] = (int *) calloc(data->ninterUnique[j], sizeof(int)));  
  }
    
  for(j=0;j<data->uniqueNum;j++) {
    cur = 0;
    for(i=0;i<*p;i++) {
      exist = 0;
      for(k=0;k<data->uniqueSize[j];k++) {
        if(data->d[i][data->mfu[j][k]] > 0) exist = 1;
      }
      if(exist) { 
        for(k=0;k<data->uniqueSize[j];k++) data->use[data->mfu[j][k]][cur] = i;
        data->useUnique[j][cur] = i;
        cur++;
      }
    }
  }
  
  /* row - learn from column */
  for(i=0;i<*p;i++) {
    data->ninterRow[i] = 0;
    data->ninterRowUnique[i] = 0;
    for(j=0;j<data->uniqueNum;j++) {
      exist = 0;
      for(k=0;k<data->uniqueSize[j];k++) {
        if(data->d[i][data->mfu[j][k]] > 0) exist = 1;
      }
      if(exist) {
        data->ninterRow[i] += data->uniqueSize[j];
        (data->ninterRowUnique[i])++;
      }
    }
  }
  for(i=0;i<*p;i++) assert(data->useRow[i] = (int *) calloc(data->ninterRow[i], sizeof(int)));  
  for(i=0;i<*p;i++) assert(data->useRowUnique[i] = (int *) calloc(data->ninterRowUnique[i], sizeof(int)));  
  for(i=0;i<*p;i++) {
    curU = 0;
    cur = 0;
    for(j=0;j<data->uniqueNum;j++) {
      exist = 0;
      for(k=0;k<data->uniqueSize[j];k++) {
        if(data->d[i][data->mfu[j][k]] > 0) exist = 1;
      }
      if(exist) {
        data->useRowUnique[i][curU] = j;
        curU++;
        for(k=0;k<data->uniqueSize[j];k++) {
          data->useRow[i][cur] = data->mfu[j][k];
          cur++;
        }
      }
    }	
  }
  
  /* Normalize Count Information */
  for(j=0;j<data->ncol;j++) data->baitCoverage[j] = log(data->baitCoverage[j] + 1.0);
  med = vec_min(data->baitCoverage, data->ncol);
  for(j=0;j<data->ncol;j++) data->baitCoverage[j] -= med; 
  for(j=0;j<data->ncol;j++) data->baitCoverage[j] *= 1.0;
    
  for(i=0;i<data->nrow;i++) data->preyAbundance[i] = log( (data->preyAbundance[i] + 1.0) / (data->preyLogLength[i]) ); 
  med = vec_min(data->preyAbundance, data->nrow);
  for(i=0;i<data->nrow;i++) data->preyAbundance[i] -= med;  
  for(i=0;i<data->nrow;i++) data->preyAbundance[i] *= 1.0;

  for(i=0;i<data->nrow;i++) data->preyLogLength[i] = log(data->preyLogLength[i] + 1.0);
  med = vec_min(data->preyLogLength, data->nrow);
  for(i=0;i<data->nrow;i++) data->preyLogLength[i] -= med; 
  for(i=0;i<data->nrow;i++) data->preyLogLength[i] *= 1.0;

  for(i=0;i<data->nrow;i++) {
    data->maxRow[i] = 0.0;
    for(j=0;j<data->ncol;j++) {
      if(data->d[i][j] > data->maxRow[i]) {
        data->maxRow[i] = data->d[i][j];
      }
    }
  } 
}

void init_param(DATA *data, PARAM *param, int *p, int *q) {
  int i;
  param->np = *p;
  param->nb = *q;
  param->nvar = 1 + param->useAbun + param->useLen + param->useCov;
  assert(param->loglikRow = (double *) calloc(*p, sizeof(double)));
  assert(param->loglikCol = (double *) calloc(*q, sizeof(double)));
  assert(param->loglikRow_tmp = (double *) calloc(*p, sizeof(double)));
  assert(param->loglikCol_tmp = (double *) calloc(*q, sizeof(double)));
  
  /* Regression Parameters */
  assert(param->beta = (double *) calloc(param->nvar, sizeof(double )));
  assert(param->gamma = (double *) calloc(param->nvar, sizeof(double )));
  /* for(j=0;j<param->nb;j++) {
    assert(param->beta[j] = (double *) calloc(param->nvar, sizeof(double)));
    assert(param->gamma[j] = (double *) calloc(param->nvar, sizeof(double)));
  } */
  assert(param->alpha_prey = (double *) calloc(*p, sizeof(double)));
  assert(param->delta_prey = (double *) calloc(*p, sizeof(double)));  
  assert(param->alpha_bait = (double *) calloc(*q, sizeof(double)));
  assert(param->delta_bait = (double *) calloc(*q, sizeof(double)));  
  assert(param->mu_prey = (double *) calloc(*p, sizeof(double)));
  assert(param->mu_prey_flag = (double *) calloc(*p, sizeof(double)));
    
  assert(param->beta_tmp = (double *) calloc(param->nvar, sizeof(double )));
  assert(param->gamma_tmp = (double *) calloc(param->nvar, sizeof(double )));
  /* for(j=0;j<param->nb;j++) {
    assert(param->beta_tmp[j] = (double *) calloc(param->nvar, sizeof(double)));
    assert(param->gamma_tmp[j] = (double *) calloc(param->nvar, sizeof(double)));
  } */
  assert(param->alpha_prey_tmp = (double *) calloc(*p, sizeof(double)));
  assert(param->delta_prey_tmp = (double *) calloc(*p, sizeof(double)));  
  assert(param->alpha_bait_tmp = (double *) calloc(*q, sizeof(double)));
  assert(param->delta_bait_tmp = (double *) calloc(*q, sizeof(double)));  
  assert(param->mu_prey_tmp = (double *) calloc(*p, sizeof(double)));
  assert(param->mu_prey_flag_tmp = (double *) calloc(*p, sizeof(double)));
  
  /* Mixture Indicators */
  assert(param->iZ = (int **) calloc(*p, sizeof(int *)));
  for(i=0;i<*p;i++) assert(param->iZ[i] = (int *) calloc(*q, sizeof(int)));
  assert(param->Z = (int **) calloc(*p, sizeof(int *)));
  for(i=0;i<*p;i++) assert(param->Z[i] = (int *) calloc(data->uniqueNum, sizeof(int)));
  assert(param->Y = (int *) calloc(*p, sizeof(int)));
  
  /* Actual Mean Parameter for 
     Real Interactors and Contaminants */
  assert(param->lambda_real = (double **) calloc(*p, sizeof(double *)));
  for(i=0;i<*p;i++) assert(param->lambda_real[i] = (double *) calloc(*q, sizeof(double)));
  assert(param->lambda_cont = (double **) calloc(*p, sizeof(double *)));
  for(i=0;i<*p;i++) assert(param->lambda_cont[i] = (double *) calloc(*q, sizeof(double)));
  assert(param->r0 = (double *) calloc(*q, sizeof(double)));
  assert(param->lambda_real0 = (double *) calloc(*q, sizeof(double)));

  assert(param->lambda_real_tmp = (double **) calloc(*p, sizeof(double *)));
  for(i=0;i<*p;i++) assert(param->lambda_real_tmp[i] = (double *) calloc(*q, sizeof(double)));
  assert(param->lambda_cont_tmp = (double **) calloc(*p, sizeof(double *)));
  for(i=0;i<*p;i++) assert(param->lambda_cont_tmp[i] = (double *) calloc(*q, sizeof(double)));
  assert(param->lambda_real0_tmp = (double *) calloc(*q, sizeof(double)));
  
  assert(param->flagged = (double *) calloc(*p, sizeof(double)));
  assert(param->appearCont = (double *) calloc(*q, sizeof(double)));
}

void set_param(PARAM *param, PRIOR *prior, DATA *data, int iter, const gsl_rng *r) {
  int i,j,ct,anypos,k;
  double tmp;
  param->loglik = - ((double) (param->np * param->nb) );
  param->loglik_tmp = param->loglik;
  for(i=0;i<param->np;i++) param->loglikRow[i] = - ((double) param->np);
  for(j=0;j<param->nb;j++) param->loglikCol[j] = - ((double) param->nb);  
  /* Setting Initial Values for Poisson Regression Parameters */
  
  for(i=0;i<param->nvar;i++) param->beta[i] = 0.0; /* gsl_ran_gaussian(r, _PSD_BETA_); */
   
  for(i=0;i<param->np;i++) {
    tmp = vec_max(data->d[i], param->nb);
    param->alpha_prey[i] = 0.0;  
    param->delta_prey[i] = 0.0;      
    param->mu_prey[i] = 0.0;
    param->mu_prey_flag[i] = 0.0;
  }
  for(j=0;j<param->nb;j++) {
    param->alpha_bait[j] = 0.0;  
    param->delta_bait[j] = 0.0;  
  }
  for(i=0;i<param->nvar;i++) param->gamma[i] = 1.0; /* gsl_ran_gaussian(r, _PSD_GAMMA_); */
  param->loglik_tmp = - ((double) (param->np * param->nb) );
  for(i=0;i<param->np;i++) param->loglikRow_tmp[i] = - ((double) param->np);
  for(j=0;j<param->nb;j++) param->loglikCol_tmp[j] = - ((double) param->nb);  
  /* Setting Initial Values for Poisson Regression Parameters */
  for(i=0;i<param->nvar;i++) param->beta_tmp[i] = 0.0;
  for(i=0;i<param->np;i++) {
    param->alpha_prey_tmp[i] = 0.0;
    param->delta_prey_tmp[i] = 0.0; /* log(vec_med(data->d[i], data->ncol)+1.0); */
    param->mu_prey_tmp[i] = param->alpha_prey_tmp[i];
    param->mu_prey_flag_tmp[i] = param->mu_prey_tmp[i];
  }
  for(j=0;j<param->nb;j++) {
    param->alpha_bait_tmp[j] = 0.0;  
    param->delta_bait_tmp[j] = 0.0;  
  }
  for(i=0;i<param->nvar;i++) param->gamma_tmp[i] = 0.0;  
    
  /* Setting Initial Values for Mixture Indicators */
  for(i=0;i<param->np;i++) {
    ct = 0;
    if(gsl_ran_flat(r,0.0,1.0) < 0.5) param->Y[i] = 1;
    else param->Y[i] = 0;
    if(param->Y[i] == 1) {
      for(j=0;j<data->uniqueNum;j++) {
        anypos = 0;
        for(k=0;k<data->uniqueSize[j];k++) if(data->d[i][data->mfu[j][k]] > 0.0) anypos = 1;
        if((gsl_ran_flat(r,0.0,1.0) < 0.5) && anypos) param->Z[i][j] = 1;
        else param->Z[i][j] = 0;
        for(k=0;k<data->uniqueSize[j];k++) if(data->d[i][data->mfu[j][k]] > 20.0) param->Z[i][j] = 1;
      }
    }
    else {
      for(j=0;j<data->uniqueNum;j++) param->Z[i][j] = 0;
    }
    if(param->Y[i] == 1) {
      for(j=0;j<param->nb;j++) {
        if(data->d[i][j] > 3.0) param->iZ[i][j] = 1;
        else param->iZ[i][j] = 0;
      }
    }
    else {
      for(j=0;j<param->nb;j++) param->iZ[i][j] = 0;
    }
  }
  
  /* Mean Parameters for All Four Possible Mixture Components */
  for(j=0;j<param->nb;j++) {
    param->r0[j] = 0.95;
    param->lambda_real0[j] = log(1.0);
  }
  for(i=0;i<param->np;i++) {
    param->flagged[i] = 0.0;
  }
  for(j=0;j<param->nb;j++) {
    param->appearCont[j] = 0.0;
  }
  calcLambdaReal(data, prior, param);
  calcLambdaCont(data, prior, param);
  
  /* Mixture Proportions */
  param->pcont[0] = 0.2;
  param->pcont[1] = 0.8;
  param->preal[0] = 0.9;
  param->preal[1] = 0.1;
}

void free_param(PARAM *param) {
  int i;
  
  free(param->loglikRow);
  free(param->loglikCol);
  free(param->loglikRow_tmp);
  free(param->loglikCol_tmp);
  
  free(param->beta);
  free(param->gamma);
  free(param->alpha_prey);
  free(param->delta_prey);
  free(param->alpha_bait);
  free(param->delta_bait);
  free(param->mu_prey);  
  free(param->mu_prey_flag);

  free(param->beta_tmp);
  free(param->gamma_tmp);
  free(param->alpha_prey_tmp);
  free(param->delta_prey_tmp);
  free(param->alpha_bait_tmp);
  free(param->delta_bait_tmp);
  free(param->mu_prey_tmp);  
 
  for(i=0;i<param->np;i++) free(param->iZ[i]);
  free(param->iZ);
  for(i=0;i<param->np;i++) free(param->Z[i]);
  free(param->Z);
  free(param->Y);
  
  for(i=0;i<param->np;i++) free(param->lambda_real[i]);
  free(param->lambda_real);
  for(i=0;i<param->np;i++) free(param->lambda_cont[i]);
  free(param->lambda_cont);
  free(param->r0);
  free(param->lambda_real0);
  for(i=0;i<param->np;i++) free(param->lambda_real_tmp[i]);
  free(param->lambda_real_tmp);
  for(i=0;i<param->np;i++) free(param->lambda_cont_tmp[i]);
  free(param->lambda_cont_tmp);
  free(param->lambda_real0_tmp);  
}

void init_prior(PRIOR *prior, int *p, int *q) {
  prior->np = *p;
  prior->nb = *q;
  prior->nvar = 1 + prior->useAbun + prior->useLen + prior->useCov;
    
  assert(prior->mean_beta = (double *) calloc(prior->nvar, sizeof(double)));
  assert(prior->var_beta = (double *) calloc(prior->nvar * prior->nvar, sizeof(double)));
  assert(prior->mean_delta_bait = (double *) calloc(*p, sizeof(double)));
  assert(prior->var_delta_bait = (double *) calloc(*p, sizeof(double)));
  assert(prior->mean_alpha_prey = (double *) calloc(*p, sizeof(double)));
  assert(prior->var_alpha_prey = (double *) calloc(*p, sizeof(double)));
  assert(prior->mean_alpha_bait = (double *) calloc(*q, sizeof(double)));
  assert(prior->var_alpha_bait = (double *) calloc(*q, sizeof(double)));
  assert(prior->mean_gamma = (double *) calloc(prior->nvar, sizeof(double)));
  assert(prior->var_gamma = (double *) calloc(prior->nvar * prior->nvar, sizeof(double)));
  assert(prior->mean_mu_prey = (double *) calloc(*p, sizeof(double)));
  assert(prior->var_mu_prey = (double *) calloc(*p, sizeof(double)));
  assert(prior->epsilon_real = (double *) calloc(*q, sizeof(double)));
  assert(prior->epsilon_cont = (double *) calloc(*q, sizeof(double)));
  assert(prior->kappa_real = (double *) calloc(*q, sizeof(double)));
  assert(prior->kappa_cont = (double *) calloc(*q, sizeof(double)));
}

void set_prior(PRIOR *prior) {
  int i;
  /* Beta */
  for(i=0;i<prior->nvar;i++) prior->mean_beta[i] = 1.0;
  for(i=0;i<prior->nvar*prior->nvar;i++) prior->var_beta[i] = 0.0;
  for(i=0;i<prior->nvar;i++) prior->var_beta[i*prior->nvar + i] = 2.0;
  for(i=0;i<prior->nvar;i++) prior->mean_gamma[i] = 2.0;
  for(i=0;i<prior->nvar*prior->nvar;i++) prior->var_gamma[i] = 0.0;
  for(i=0;i<prior->nvar;i++) prior->var_gamma[i*prior->nvar + i] = 2.0;
  prior->mean_gamma[0] = 1.0;
  prior->var_gamma[0] = 0.01;
  
  for(i=0;i<prior->np;i++) {
    prior->mean_alpha_prey[i] = 0.0;
    prior->var_alpha_prey[i] = 2.00;
  }
  for(i=0;i<prior->np;i++) {
    prior->mean_delta_bait[i] = 0.0;
    prior->var_delta_bait[i] = 2.00;
  }
  for(i=0;i<prior->nb;i++) {
    prior->mean_alpha_bait[i] = 0.0;
    prior->var_alpha_bait[i] = 2.00;
  }
  for(i=0;i<prior->np;i++) {
    prior->mean_mu_prey[i] = 0.0;
    prior->var_mu_prey[i] = 2.00;
  }
  for(i=0;i<prior->nb;i++) {
    prior->epsilon_real[i] = ((double) prior->np) * 0.015;
    prior->kappa_real[i] = ((double) prior->np) * 0.01;
    prior->epsilon_cont[i] = 1.0;
    prior->kappa_cont[i] = 1.0; 
  }
  prior->acont[0] = 1.0;
  prior->acont[1] = 1.0;
  prior->areal[0] = 1.0;
  prior->areal[1] = 1.0;

  prior->sigmasq_alpha_prey = 1.0;
  prior->sigmasq_alpha_bait = 1.0;
  prior->sigmasq_mu_prey = 1.0;
  prior->shape_alpha_prey = 10.0 * (0.01 * ((double) prior->np));
  prior->rate_alpha_bait = 10.0 * (0.01 * ((double) prior->np));
  prior->shape_alpha_bait = 10.0 * (0.01 * ((double) prior->np));
  prior->rate_alpha_bait = 10.0 * (0.01 * ((double) prior->np));
  prior->shape_mu_prey = 10.0 * (0.01 * ((double) prior->np));
  prior->rate_mu_prey = 10.0 * (0.01 * ((double) prior->np));
 

}

void free_prior(PRIOR *prior) {
  free(prior->mean_beta);
  free(prior->var_beta);
  free(prior->mean_alpha_prey);
  free(prior->var_alpha_prey);
  free(prior->mean_delta_bait);
  free(prior->var_delta_bait);
  free(prior->mean_alpha_bait);
  free(prior->var_alpha_bait);
  free(prior->mean_gamma);
  free(prior->var_gamma);
  free(prior->mean_mu_prey);
  free(prior->var_mu_prey);
  free(prior->epsilon_real);
  free(prior->epsilon_cont);
  free(prior->kappa_real);
  free(prior->kappa_cont);
}

void init_summary(DATA *data, SUMMARY *summary, int *p, int *q) {
  int i;
  summary->np = *p;
  summary->nb = *q;
  assert(summary->iZ = (double **) calloc(*p, sizeof(double *)));
  for(i=0;i<*p;i++) assert(summary->iZ[i] = (double *) calloc(*q, sizeof(double)));
  assert(summary->Z = (double **) calloc(*p, sizeof(double *)));
  for(i=0;i<*p;i++) assert(summary->Z[i] = (double *) calloc(data->uniqueNum, sizeof(double)));
  /* assert(summary->lambda_real = (double **) calloc(*p, sizeof(double *)));
  for(i=0;i<*p;i++) assert(summary->lambda_real[i] = (double *) calloc(*q, sizeof(double)));
  assert(summary->lambda_cont = (double **) calloc(*p, sizeof(double *)));
  for(i=0;i<*p;i++) assert(summary->lambda_cont[i] = (double *) calloc(*q, sizeof(double))); */
  assert(summary->expect = (double **) calloc(*p, sizeof(double *)));
  for(i=0;i<*p;i++) assert(summary->expect[i] = (double *) calloc(*q, sizeof(double)));  
  /* assert(summary->lambda_real0 = (double *) calloc(*q, sizeof(double)));  */
  assert(summary->Y = (double *) calloc(*p, sizeof(double)));
  assert(summary->max_prob = (double *) calloc(*p, sizeof(double)));
}

void set_summary(DATA *data, SUMMARY *summary) {
  int i,j;
  /* Reproducibility score */
  for(i=0;i<summary->np;i++) {
    for(j=0;j<summary->nb;j++) {
      summary->iZ[i][j] = 0.0;
    }
  }

  /* for(j=0;j<summary->nb;j++) summary->lambda_real0[j] = 0.0; */

  for(i=0;i<summary->np;i++) {
    summary->Y[i] = 0.0;
    for(j=0;j<data->uniqueNum;j++) summary->Z[i][j] = 0.0;
    for(j=0;j<summary->nb;j++) {
      /* summary->lambda_real[i][j] = 0.0;
      summary->lambda_cont[i][j] = 0.0; */
      summary->expect[i][j] = 0.0;
    }
  }
}

void free_summary(SUMMARY *summary) {
  int i;
  for(i=0;i<summary->np;i++) {
    free(summary->iZ[i]);
    free(summary->Z[i]);
    /* free(summary->lambda_real[i]);
    free(summary->lambda_cont[i]); */
    free(summary->expect[i]);
  }
  free(summary->iZ);
  free(summary->Z);
  /* free(summary->lambda_real);
  free(summary->lambda_cont); */
  free(summary->expect);
  free(summary->Y);
  /* free(summary->lambda_real0); */
  free(summary->max_prob);
}

