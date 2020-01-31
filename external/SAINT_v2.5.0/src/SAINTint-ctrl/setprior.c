#include "saint.h"

void memory_prior(PARAM *param, PRIOR *prior, DATA *data) {
  assert(prior->w_alpha_prey = (int *) calloc(data->nprey, sizeof(int)));
  assert(prior->w_alpha_IP = (int *) calloc(data->nIP, sizeof(int)));
  assert(prior->w_mu = (int *) calloc(data->nprey, sizeof(int)));
  assert(prior->w_eta = (int *) calloc(data->nprey, sizeof(int)));
  assert(prior->w_eta0 = (int *) calloc(data->nprey, sizeof(int)));
}

void initialize_prior(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r) {
  int i,pass;
  float mean;
  float MAXC = ((float) _MAX_COMP_);
  prior->m_beta = 0.0;
  prior->v_beta = 100.0;
  prior->atrue = 0.1;
  prior->afalse = 1.0 - prior->atrue;

  prior->rho_alpha_prey = GSL_MAX(((float) data->nprey) * 0.01, 1.0);
  prior->m_alpha_prey = 0.0;
  prior->v_alpha_prey = 10.0;
  for(i=0;i<_MAX_COMP_;i++) prior->gamma_alpha_prey[i] = 1.0 / MAXC;  
  for(i=0;i<_MAX_COMP_;i++) prior->theta_alpha_prey[i] = gsl_ran_gaussian(r, 2.0) + 2.0;  
  for(i=0;i<data->nprey;i++) prior->w_alpha_prey[i] = ((int) gsl_ran_flat(r,0.0,MAXC)); 

  mean = 0.0;
  for(i=0;i<_MAX_COMP_;i++) mean += prior->gamma_alpha_prey[i] * prior->theta_alpha_prey[i];
  for(i=0;i<_MAX_COMP_;i++) prior->theta_alpha_prey[i] -= mean;

  prior->rho_alpha_IP = 1.0;
  prior->m_alpha_IP = 0.0;
  prior->v_alpha_IP = 10.0;
  for(i=0;i<_MAX_COMP_;i++) prior->gamma_alpha_IP[i] = 1.0 / MAXC;  
  /* for(i=0;i<_MAX_COMP_;i++) prior->theta_alpha_IP[i] = gsl_ran_gaussian(r, 2.0);  */
  for(i=0;i<_MAX_COMP_;i++) prior->theta_alpha_IP[i] = 0.0;  
  /* for(i=0;i<data->nIP;i++) prior->w_alpha_IP[i] = ((int) gsl_ran_flat(r,0.0,MAXC)); */
  for(i=0;i<data->nIP;i++) prior->w_alpha_IP[i] = 0; 
  mean = 0.0;
  for(i=0;i<_MAX_COMP_;i++) mean += prior->gamma_alpha_IP[i] * prior->theta_alpha_IP[i];
  for(i=0;i<_MAX_COMP_;i++) prior->theta_alpha_IP[i] -= mean;

  prior->rho_mu = GSL_MAX(((float) data->nprey) * 0.01, 1.0);
  prior->m_mu = 0.0;
  prior->v_mu = 10.0;
  for(i=0;i<_MAX_COMP_;i++) prior->gamma_mu[i] = 1.0 / MAXC;  
  for(i=0;i<_MAX_COMP_;i++) prior->theta_mu[i] = gsl_ran_gaussian(r, 2.0) - 2.0;  
  for(i=0;i<data->nprey;i++) prior->w_mu[i] = ((int) gsl_ran_flat(r,0.0,MAXC)); 
  mean = 0.0;
  for(i=0;i<_MAX_COMP_;i++) mean += prior->gamma_mu[i] * prior->theta_mu[i];
  for(i=0;i<_MAX_COMP_;i++) prior->theta_mu[i] -= mean;

  prior->rho_eta = GSL_MAX(((float) data->nprey) * 0.01, 1.0);
  prior->shape_eta = 1.0;
  prior->scale_eta = 1.0;
  for(i=0;i<_MAX_COMP_;i++) prior->gamma_eta[i] = 1.0 / MAXC;  
  for(i=0;i<_MAX_COMP_;i++) {
    pass = 0;
    while(!pass){
      prior->theta_eta[i] = gsl_ran_flat(r, 0.0, 1.0) ;  
      if(prior->theta_eta[i] < 1.0) pass = 1;
    }
  }
  for(i=0;i<data->nprey;i++) prior->w_eta[i] = ((int) gsl_ran_flat(r,0.0,MAXC)); 

  prior->rho_eta0 = GSL_MAX(((float) data->nprey) * 0.01, 1.0);
  prior->shape_eta0 = 1.0;
  prior->scale_eta0 = 1.0;
  for(i=0;i<_MAX_COMP_;i++) prior->gamma_eta0[i] = 1.0 / MAXC;  
  for(i=0;i<_MAX_COMP_;i++) {
    pass = 0;
    while(!pass){
      prior->theta_eta0[i] = gsl_ran_flat(r, 0.0, 1.0) ;  
      if(prior->theta_eta0[i] < 1.0) pass = 1;
    }
  }
  for(i=0;i<data->nprey;i++) prior->w_eta0[i] = ((int) gsl_ran_flat(r,0.0,MAXC)); 
}

void set_prior(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r) {
  memory_prior(param, prior, data);
  initialize_prior(param, prior, data, r);
}



