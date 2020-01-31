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
#include <unistd.h>
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

#define _MAX_BUF_        2000
#define _MAX_NAME_       2000
#define _MAX_COUNT_	250
#define _MAX_COMP_         15
#define _SKIP_             10
#define _PRINT_FREQ_       100
#define _HISTO_START_ -5.0
#define _HISTO_END_    5.0
#define _HISTO_BIN_   100
#define _HISTO_START2_ 0.5
#define _HISTO_END2_   20.5
#define _HISTO_BIN2_   200
#define _TRUNC_ 1000.0
#define _FR_ 0.01

typedef struct tagDATA {

  /*************/
  /* logistics */
  /*************/
  int ninter;
  int nuinter;
  int nprey;
  int nIP;
  int nbait;

  /**************************/
  /* interaction level data */
  /**************************/
  char **prey;
  char **bait;
  char **ip;   /* raw data, each row corresponds to one interaction, case-sensitive */
  float *d;
  float *d2;
  float *iprob;
  float *l;
  float *c;
  int *isCtrl;

  /*********************************/
  /* unique interaction level data */
  /*********************************/
  char **uprey;
  char **ubait;
  float *prob;
  int *isAnyCtrl;

  int *n_u2a; /* number of individual interactions per unique interactions */
  int **u2a;  /* unique interactions to individual interactions */  
  int *a2u;   /* individual interactions to unique interactions */ 
  /* crucial indicator for probability calculation */

  /***********************************/
  /* unique bait and prey level data */
  /***********************************/
  float *IPtotalAbundance;
  float *IPbaitCoverage;
  float *preyLength;
  char *preyOverride;

  char **PREY;  /* unique preys */
  char **PREYGENE;  /* unique preys */
  char **BAIT;  /* unique baits */
  char **IP;    /* unique IP #s */

  int nctrl;
  int ntest;
  int *ctrl;  /* index: control IPs or not: 'C' = control, 'T' = test */

  float *ctrlavg;

  int *preyNinter;  /* # interaction for prey */
  int *preyFlag;    /* if preyNinter is larger than freq */
  int *baitNinter;  /* # interaction for bait */
  int *IPNinter;    /* # interaction in an IP */
  int *baitNIP;     /* # IPs per bait         */

  /****************/
  /* mapping data */
  /****************/
  int *i2p;   /* index: interaction to prey */
  int *i2b;   /* index: interaction to bait */
  int *i2IP;  /* index: interaction to IP   */ 

  int **p2i;  /* index: prey to interaction */
  int **b2i;  /* index: bait to interaction */
  int **IP2i; /* index: IP to interaction   */ 

  int *ui2p;   /* index: unique interaction to prey */
  int *ui2b;   /* index: unique interaction to bait */
  /* no need to build reverse mapping for unique interactions */
  /* perhaps this mapping is unnecessary */

  int **b2IP; /* index: bait to IP */
  int *IP2b;  /* index: IP to bait */

} DATA;

typedef struct tagPARAM{
  float freq;
  float freqgroup;
  float modelvar;

  float loglikTotal;
  float *loglik_prey;
  float *loglik_IP;

  float beta0;
  float betac;
  float *alpha_prey;
  float *alpha_IP;
  float *mu;
  float *eta;
  float *eta0;

  int *iZ; /* individual interactions */
  int *Z; /* unique interactions */

  float ptrue;
  float ptrue_tmp;
  float *lambda_true;
  float *lambda_false;
  float *lambda_true_tmp;
  float *lambda_false_tmp;

} PARAM;

typedef struct tagPRIOR{

  /* parametric portion */
  float m_beta;  /* set to zero */
  float v_beta;
  float atrue, afalse;

  /* nonparametric portion */
  float rho_alpha_prey;
  float m_alpha_prey;
  float v_alpha_prey;
  int *w_alpha_prey;
  float gamma_alpha_prey[_MAX_COMP_];
  float theta_alpha_prey[_MAX_COMP_];

  float rho_alpha_IP;
  float m_alpha_IP;
  float v_alpha_IP;
  int *w_alpha_IP;
  float gamma_alpha_IP[_MAX_COMP_];
  float theta_alpha_IP[_MAX_COMP_];

  float rho_mu;
  float m_mu;
  float v_mu;
  int *w_mu;
  float gamma_mu[_MAX_COMP_];
  float theta_mu[_MAX_COMP_];

  float rho_eta;
  float mean_eta;
  int *w_eta;
  float gamma_eta[_MAX_COMP_];
  float theta_eta[_MAX_COMP_];

  float rho_eta0;
  float mean_eta0;
  int *w_eta0;
  float gamma_eta0[_MAX_COMP_];
  float theta_eta0[_MAX_COMP_];
} PRIOR;

typedef struct tagHISTOGRAM{
  float start[_HISTO_BIN_];
  float end[_HISTO_BIN_];
  float count[_HISTO_BIN_ + 2];
} HISTOGRAM;

typedef struct tagHISTOGRAM2{
  float start[_HISTO_BIN2_];
  float end[_HISTO_BIN2_];
  float count[_HISTO_BIN2_ + 2];
} HISTOGRAM2;


typedef struct tagSUMMARY{
  float freq;
  float *iZ;
  float *Z;
  float *alpha_prey;
  float *alpha_IP;
  float *mu;
  float *eta;
  float *eta0;
  float *lambda_true;
  float *lambda_false;
    float *FDR;
    
  HISTOGRAM hist_alpha_prey;
  HISTOGRAM hist_alpha_IP;
  HISTOGRAM hist_mu;
  HISTOGRAM2 hist_eta;
  HISTOGRAM2 hist_eta0;
} SUMMARY;


/* GLOBAL VARIABLES */
int burn;
int iter;
float freq;
float freqgroup;
int modelvar;
int NORMALIZE;

/*************/
/* functions */
/*************/

int nrow(FILE *fp);
int newlinechar(char *buf, int k);
int ncol(FILE *fp);
int commandLine(int argc, char **argv);
void print_DP(PRIOR *prior, DATA *data);

/* initdata.c */
void read_interaction_data(FILE *fpinter, DATA *data);
void find_unique_interaction(DATA *data);
int unique_elements(char **x, int *unique, int nx);
int count_unique_elements(char **x, int nx);
void centerData(float *x, int n, int takelog);
int mapPreyToData(DATA *data);
void read_prey_data(FILE *fpprey, DATA *data);
void mapIPtoBait(DATA *data);
int mapIPBaitToData(DATA *data);
void getIPinfo(DATA *data);
void read_bait_data(FILE *fpbait, DATA *data);

void set_ctrlavg(DATA *data);
void read_data(FILE *fpinter, FILE *fpprey, FILE *fpbait, DATA *data, float *freq, float *freqgroup);
void prey_flag(DATA *data, float *freq, float *freqgroup);

/* meancounts.c */
void compute_lambda_true_all(PARAM *param, PRIOR *prior, DATA *data);
void compute_lambda_false_all(PARAM *param, PRIOR *prior, DATA *data);
void compute_lambda_all(PARAM *param, PRIOR *prior, DATA *data);
/* void compute_lambda_true_prey(PARAM *param, PRIOR *prior, DATA *data, int pid);
void compute_lambda_false_prey(PARAM *param, PRIOR *prior, DATA *data, int pid);
void compute_lambda_prey(PARAM *param, PRIOR *prior, DATA *data, int pid);
void compute_lambda_true_IP(PARAM *param, PRIOR *prior, DATA *data, int ipid);
void compute_lambda_IP(PARAM *param, PRIOR *prior, DATA *data, int ipid); */

/* likelihood.c */
float log_poisson_prop(float N, float lambda);
float log_poisson_g_prop(float N, float lambda, float theta);
float LRprop(PARAM *param, PRIOR *prior, DATA *data);
float loglik_all(PARAM *param, PRIOR *prior, DATA *data);
float loglik_all_class(PARAM *param, PRIOR *prior, DATA *data, int cl);
float loglik_all_class_tmp(PARAM *param, PRIOR *prior, DATA *data, int cl);

/* mcmc.c */
float log_gaussian(float x, float mu, float var);
void sampleBeta0(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r);
void sampleBetac(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r);
void sampleZ(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r);
void contaminant(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r);
float logit(float x);
float inverseLogit(float x);
void sampleProportion(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r);
void mhgibbs(PARAM *param, PRIOR *prior, DATA *data, SUMMARY *summary, const gsl_rng *r, int updateSum);
void write_mcmc(PARAM *param, PRIOR *prior, DATA *data, FILE *fp1, FILE *fp2, FILE *fp3, int ct);

/* printmap.c */
void printInter(DATA *data);
void printUInter(DATA *data);
void printIP(DATA *data);
void printBait(DATA *data);
void printPrey(DATA *data);
void printMap(DATA *data);

/* setprior.c */
void memory_prior(PARAM *param, PRIOR *prior, DATA *data);
void initialize_prior(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r);
void set_prior(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r);

/* setparam.c */
void memory_param(PARAM *param, PRIOR *prior, DATA *data);
void initialize_param(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r);
void set_param(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r);
void set_Z(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r);

/* setsumamry.c */
void memory_summary(SUMMARY *summary, DATA *data);
void initialize_summary(SUMMARY *summary, DATA *data);
void initialize_histogram(HISTOGRAM *hist);
void initialize_histogram2(HISTOGRAM2 *hist);
void set_summary(SUMMARY *summary, DATA *data);
void calculateFDR(DATA *data, SUMMARY *summary);


void updateSummary(PARAM *param, PRIOR *prior, DATA *data, SUMMARY *summary);
void scaleSummary(SUMMARY *summary, DATA *data, int iter);
void updateHist_alpha_prey(HISTOGRAM *hist, PRIOR *prior);
void updateHist_alpha_IP(HISTOGRAM *hist, PRIOR *prior);
void updateHist_mu(HISTOGRAM *hist, PRIOR *prior);
void updateHist_eta(HISTOGRAM2 *hist, PRIOR *prior);
void updateHist_eta0(HISTOGRAM2 *hist, PRIOR *prior);
void updateHistogram(PARAM *param, PRIOR *prior, DATA *data, SUMMARY *summary);

/* dpalphaprey.c */
void DP_alpha_prey_gamma(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r);
void DP_alpha_prey_w(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r, int pid);
void DP_alpha_prey_theta(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r, int cid, int *inuse);
void DP_alpha_prey(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r);

/* dpalphaIP.c */
void DP_alpha_IP_gamma(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r);
void DP_alpha_IP_w(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r, int pid);
void DP_alpha_IP_theta(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r, int cid, int *inuse);
void DP_alpha_IP(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r);

/* dpmu.c */
void DP_mu_gamma(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r);
void DP_mu_w(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r, int pid);
void DP_mu_theta(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r, int cid, int *inuse);
void DP_mu(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r);

/* dpmu.c */
void DP_eta_gamma(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r);
void DP_eta_w(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r, int pid);
void DP_eta_theta(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r, int cid, int *inuse);
void DP_eta(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r);

/* dpmu.c */
void DP_eta0_gamma(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r);
void DP_eta0_w(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r, int pid);
void DP_eta0_theta(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r, int cid, int *inuse);
void DP_eta0(PARAM *param, PRIOR *prior, DATA *data, const gsl_rng *r);

/* result.c */
void write_interactions(DATA *data, SUMMARY *summary);
void write_unique_interactions(DATA *data, SUMMARY *summary);
void write_prey(DATA *data, SUMMARY *summary);
void write_IP(DATA *data, SUMMARY *summary);
void write_bait(DATA *data, SUMMARY *summary);
void write_histogram(FILE *fp, HISTOGRAM *hist);
void write_histogram2(FILE *fp, HISTOGRAM2 *hist);
void write_hyperprior(DATA *data, SUMMARY *summary);
void write_result(DATA *data, SUMMARY *summary);
void write_matrix_data(DATA *data, SUMMARY *summary);
void write_matrix_data2(DATA *data, SUMMARY *summary);


/****************   mmath.c   ******************************/
float vec_sum(const float *vec, int len);
float vec_max(const float *vec, int len);
float vec_min(const float *vec, int len);
float vec_mean(const float *vec, int len);
float vec_var(const float *vec, int len);
float vec_med(const float *vec, int len);
float vec_mad(const float *vec, int len);
float geometric_mean(float *x, int n);
int ranMultinom(const gsl_rng *r, float *p, int K);

