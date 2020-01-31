/*
    Copyright (C) <2011>  <Hyungwon Choi>
    For troubleshooting, contact hyung_won_choi@nuhs.edu.sg.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You can obtain a copy of the GNU General Public License from
    <http://www.gnu.org/licenses/>.
*/




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <ctype.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>

#define _MAX_PROT_       10000
#define _MAX_BUF_        2000
#define _MAX_NAME_       2000
#define _MAX_SAMPLE_      10000
#define _MAX_BAIT_        10000
#define _MAX_REPLICA_      100
#define _FF_PROP_        0.10
#define _SKIP_             1
#define _PRINT_FREQ_       100
#define _PSD_BETA_        1.0
#define _PSD_GAMMA_       1.0
#define _PSD_ALPHA_PREY_  2.0
#define _PSD_ALPHA_BAIT_  2.0
#define _PSD_DELTA_BAIT_  1.0
#define _PSD_MU_PREY_     1.0

typedef struct tagDATA {
  char **preyORF;
  char **prey;
  char **experiment;
  char **bait;
  double *baitCoverage;
  double *preyAbundance;
  double *preyLogLength;
  double **d;
  double *maxRow;  

  char *override;

  int nrow; 
  int ncol; 

  char **unique;
  int uniqueNum;
  int *uniqueSize;
  int *mtu; /* mapToUnique */
  int **mfu; /* mapFromUnique */

  int *ninter;
  int **use;
  int *ninterUnique;
  int **useUnique;
  int *ninterRow;
  int **useRow;
  int *ninterRowUnique;
  int **useRowUnique;
} DATA;

typedef struct tagPARAM{
  int np;
  int nb;
  int useAbun;
  int useLen;
  int useCov;
  int nvar;

  double loglik;
  double *loglikRow;
  double *loglikCol;
  double loglik_tmp;
  double *loglikRow_tmp;
  double *loglikCol_tmp;
  
  double *beta;
  double *alpha_prey;
  double *delta_prey;
  double *alpha_bait;
  double *delta_bait;
  double *gamma;
  double *mu_prey;
  double *mu_prey_flag;
  int **iZ;
  int **Z;
  int *Y;
  double pcont[2];
  double preal[2];
  double **lambda_real;
  double **lambda_cont;
  double *r0;
  double *lambda_real0;
  
  double *beta_tmp;
  double *alpha_prey_tmp;
  double *delta_prey_tmp;
  double *alpha_bait_tmp;
  double *delta_bait_tmp;
  double *gamma_tmp;
  double *mu_prey_tmp;
  double *mu_prey_flag_tmp;
  double **lambda_real_tmp;
  double **lambda_cont_tmp;
  double *lambda_real0_tmp;

  double ff_prop;
  double *flagged;
  double *appearCont;

} PARAM;

typedef struct tagPRIOR{
  int np;
  int nb;
  int useAbun;
  int useLen;
  int useCov;
  int nvar;
  double *mean_beta;
  double *var_beta;
  double *mean_alpha_prey;
  double *var_alpha_prey;
  double *mean_delta_bait;
  double *var_delta_bait;
  double *mean_alpha_bait;
  double *var_alpha_bait;
  double *mean_gamma;
  double *var_gamma;
  double *mean_mu_prey;
  double *var_mu_prey;
  double *epsilon_real;
  double *kappa_real;
  double *epsilon_cont;
  double *kappa_cont;
  double acont[2];
  double areal[2];
  
  double sigmasq_alpha_prey;
  double sigmasq_alpha_bait;
  double sigmasq_mu_prey;
  double shape_alpha_prey;
  double rate_alpha_prey;
  double shape_alpha_bait;
  double rate_alpha_bait;
  double shape_mu_prey;
  double rate_mu_prey;

} PRIOR;

typedef struct tagSUMMARY{
  int np;
  int nb;
  double **iZ;
  double **Z;
  double **lambda_real;
  double **lambda_cont;
  double *lambda_real0;
  double **expect;
  double *Y;
  double *max_prob;
} SUMMARY;


/*************/
/* functions */
/*************/

int nrow(FILE *fp);
int newlinechar(char *buf, int k);
int ncol(FILE *fp);

/****************   init_data.c  ***************************/
void init_data(DATA *data, int *p, int *q);
void free_data(DATA *data);
void read_data(FILE *fp, DATA *data, int *p, int *q);
void init_param(DATA *data, PARAM *param, int *p, int *q);
void set_param(PARAM *param, PRIOR *prior, DATA *data, int iter, const gsl_rng *r);
void free_param(PARAM *param);
void init_prior(PRIOR *prior, int *p, int *q);
void set_prior(PRIOR *prior);
void free_prior(PRIOR *prior);
void init_summary(DATA *data, SUMMARY *summary, int *p, int *q);
void set_summary(DATA *data, SUMMARY *summary);
void free_summary(SUMMARY *summary);
void normalizeCoverage(double *cover, int *len);

/****************   update.c   *****************************/
double logZIP(double x, double r, double lambda);
double loglikAll(DATA *data, PRIOR *prior, PARAM *param);
double loglikAll_realtmp(DATA *data, PRIOR *prior, PARAM *param);
double loglikAll_conttmp(DATA *data, PRIOR *prior, PARAM *param);
double loglikRow(DATA *data, PRIOR *prior, PARAM *param, int r);
double loglikCol(DATA *data, PRIOR *prior, PARAM *param, int c);
void calcLambdaRealRow(DATA *data, PRIOR *prior, PARAM *param, int r);
void calcLambdaRealCol(DATA *data, PRIOR *prior, PARAM *param, int c);
void calcLambdaReal(DATA *data, PRIOR *prior, PARAM *param);
void calcLambdaReal_tmp(DATA *data, PRIOR *prior, PARAM *param);
void calcLambdaContRow(DATA *data, PRIOR *prior, PARAM *param, int r);
void calcLambdaContCol(DATA *data, PRIOR *prior, PARAM *param, int c);
void calcLambdaCont(DATA *data, PRIOR *prior, PARAM *param);
void calcLambdaCont_tmp(DATA *data, PRIOR *prior, PARAM *param);
double dMultGauss(double *x, double *mu, double *Sigma, int length);
void updateBeta(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r, int k);
void updateAlphaDeltaPrey(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r);
void updateAlphaDeltaBait(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r);
void updateGamma(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r, int k);
void updateMuPrey(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r);
void updateLambdaReal0(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r);
void updatePriors(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r);

void updateSigmasqAlphaPrey(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r);
void updateSigmasqAlphaBait(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r);
void updateSigmasqMuPrey(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r);

void sampleY(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r);
void sampleZ(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r);
void sampleP(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r);
void mhgibbs(DATA *data, PRIOR *prior, PARAM *param, const gsl_rng *r);

double calc_dist(int calls[], int num);

/****************   mmath.c   ******************************/
double vec_sum(const double *vec, int len);
double vec_max(const double *vec, int len);
double vec_min(const double *vec, int len);
double vec_mean(const double *vec, int len);
double vec_var(const double *vec, int len);
double vec_med(const double *vec, int len);
double vec_mad(const double *vec, int len);
double poisson_unscaled_pdf(const double c, const double lambda);
