#include "saint.h"

int nrow(FILE *fp) {
  char buf[10000];
  int n = 0;
  while(fgets(buf, sizeof(buf), fp) != NULL) n++;
  return n;
}

int newlinechar(char *buf, int k) {
  int i;
  int found = 0;
  for(i=0;i<k;i++) {
    if(buf[i] == '\n') {
      found = 1;
      break;
    }
  }
  return found;
}

int ncol(FILE *fp) {
  char buf[10000];
  int i,cont = 0;
  fgets(buf, sizeof(buf), fp);
  for(i=0;i<10000;i++) {
    if(buf[i] == '\t' || buf[i] == ' ') cont++;
    if(buf[i] == '\0') break;
  }
  return cont;
}


void print_DP(PRIOR *prior, DATA *data) { 
  int i; 
  /* fprintf(stderr, "\nDP_alpha_IP\n");
  for(i=0;i<_MAX_COMP_;i++) fprintf(stderr, "%.2f\t", prior->gamma_alpha_IP[i]);
  fprintf(stderr, "\n");
  for(i=0;i<_MAX_COMP_;i++) fprintf(stderr, "%.2f\t", prior->theta_alpha_IP[i]);
  fprintf(stderr, "\n"); */

  fprintf(stderr, "\nDP_alpha_prey\n");
  for(i=0;i<_MAX_COMP_;i++) fprintf(stderr, "%.2f\t", prior->gamma_alpha_prey[i]);
  fprintf(stderr, "\n");
  for(i=0;i<_MAX_COMP_;i++) fprintf(stderr, "%.2f\t", prior->theta_alpha_prey[i]);
  fprintf(stderr, "\n");

  fprintf(stderr, "DP_mu\n");
  for(i=0;i<_MAX_COMP_;i++) fprintf(stderr, "%.2f\t", prior->gamma_mu[i]);
  fprintf(stderr, "\n");
  for(i=0;i<_MAX_COMP_;i++) fprintf(stderr, "%.2f\t", prior->theta_mu[i]);
  fprintf(stderr, "\n");

  fprintf(stderr, "DP_eta\n");
  for(i=0;i<_MAX_COMP_;i++) fprintf(stderr, "%.2f\t", prior->gamma_eta[i]);
  fprintf(stderr, "\n");
  for(i=0;i<_MAX_COMP_;i++) fprintf(stderr, "%.2f\t", prior->theta_eta[i]);
  fprintf(stderr, "\n");

  fprintf(stderr, "DP_eta0\n");
  for(i=0;i<_MAX_COMP_;i++) fprintf(stderr, "%.2f\t", prior->gamma_eta0[i]);
  fprintf(stderr, "\n");
  for(i=0;i<_MAX_COMP_;i++) fprintf(stderr, "%.2f\t", prior->theta_eta0[i]);
  fprintf(stderr, "\n"); 

  fprintf(stderr, "\n");
}

int commandLine(int argc, char **argv) {
  if (argc < 7) {
    fprintf(stderr, "\nusage: saint-spc-noctrl [interactionFile] [preyFile] [baitFile] [nburnin] [niter] [fthres] [fgroup] [var] [normalize] \n\n-burnin and iter: burn-in period and main iteration of MCMC\n-fthres: frequency threshold above which probability is set to 0\n-fgroup: frequency boundary dividing high and low frequency groups\n-var: binary [0/1] indicator for modeling variance of the count data distributions\n-normalize: whether to normalize the counts by total spectral counts\n\n");
    return 1;
  }
  /* interaction file: IPnumber \t bait \t prey \t spectralCount \n */
  /* prey file:        prey \t sequenceLength \n */
  /* bait file:        bait \t IPnumber \t isControl \n */

  FILE *fpinter = fopen(argv[1], "r");
  FILE *fpprey = fopen(argv[2], "r");
  FILE *fpbait = fopen(argv[3], "r");

  if(fpinter == NULL) { 
    fprintf(stderr, "Cannot locate interaction data %s.\n", argv[1]);
    return 1; 
  }

  if(fpprey == NULL) { 
    fprintf(stderr, "Cannot locate prey data %s.\n", argv[2]);
    return 1; 
  }

  if(fpbait == NULL) { 
    fprintf(stderr, "Cannot locate bait data %s.\n", argv[3]);
    return 1; 
  }
  
  if(argc < 5) { 
    fprintf(stderr, "The number of burnin was not provided. Set to 2,000.\n");
    burn = 2000;
  }
  else {
    burn = atoi(argv[4]);
  }

  if(argc < 6) { 
    fprintf(stderr, "The number of main interations was not provided. Set to 10,000.\n");
    iter = 10000;
  }
  else {
    iter = atoi(argv[5]);
  }

  if(argc < 7) { 
    fprintf(stderr, "The frequency threshold was not provided. Set to 0.1.\n");
    freq = 0.1;
  }
  else {
    freq = atof(argv[6]);
  }

  if(argc < 8) { 
    fprintf(stderr, "The frequency group boundary was not provided. Set to 0.01.\n");
    freqgroup = 0.01;
  }
  else {
    freqgroup = atof(argv[7]);
  }

  if(argc < 9) { 
    fprintf(stderr, "The indicator for variance modelling was not provided. Set to 0 (Yes).\n");
    modelvar = 0;
  }
  else {
    modelvar = atoi(argv[8]);
  }

  if(argc < 9) { 
    fprintf(stderr, "The indicator for normalization was not provided. Set to 0 (Yes).\n");
    NORMALIZE = 0;
  }
  else {
    NORMALIZE = atoi(argv[9]);
  }

  fclose(fpinter);
  fclose(fpprey);
  fclose(fpbait);
  return 0;
}


/***************************** MAIN ***************************/

int main(int argc, char **argv) {

  int i,ct;
  DATA data; 
  PARAM param;
  PRIOR prior;
  SUMMARY summary;
        
  const gsl_rng_type *T;
  gsl_rng *r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  /* Command Line */
  if(commandLine(argc, argv)) return 1;
  FILE *fpinter = fopen(argv[1], "r");
  FILE *fpprey = fopen(argv[2], "r");
  FILE *fpbait = fopen(argv[3], "r");

  /* Read interaction data, identify baits, preys, and IPs, 
     make unique interaction data frame, 
     identify the mapping between different levels of data */	
  system("mkdir LOG");      /* error logs */
  system("mkdir MAPPING");      /* mapping logs */
  system("mkdir MCMC");     /* posterior samples */
  system("mkdir RESULT");   /* posterior probabilities, other summaries */

  fprintf(stderr, "Reading data and mapping interactions\n");
  read_data(fpinter, fpprey, fpbait, &data, &freq, &freqgroup);
  printMap(&data);

  /* Set up model parameters and prior elicitation */
  set_prior(&param, &prior, &data, r);
  set_param(&param, &prior, &data, r);
  set_summary(&summary, &data);
  param.freq = freq;
  param.freqgroup = freqgroup;
  param.modelvar = modelvar;
  summary.freq = freq;
	
  /* updates and summary */	
  chdir("MCMC");
  FILE *fp1 = fopen("alpha_prey","w");
  FILE *fp2 = fopen("alpha_IP","w");
  FILE *fp3 = fopen("mu","w");

  /* burnin */
  ct = 0;
  fprintf(stderr, "Burn-in Period\n");

  for(i=0;i<burn;i++) {

    mhgibbs(&param, &prior, &data, &summary, r, 0);
    if((i+1) % _SKIP_ == 0) {
      write_mcmc(&param, &prior, &data, fp1, fp2, fp3, ct);
      ct++;
    }
    if((i+1) % _PRINT_FREQ_ == 0) {
      fprintf(stderr, "%d\t", i+1);
      /* fprintf(stderr, "%d (%.2f,%.2f,%.2f)\t", i+1, param.ptrue, param.beta0, param.betac); */
      /* print_DP(&prior, &data); */
    }
  }
  fprintf(stderr, "\n");

  /* main iterations */
  fprintf(stderr, "Main Iterations\n");
  ct = 0;
  for(i=0;i<iter;i++) {
    /* code for mh-gibbs */
    mhgibbs(&param, &prior, &data, &summary, r, 1);
    if((i+1) % _SKIP_ == 0) {
      write_mcmc(&param, &prior, &data, fp1, fp2, fp3, ct);
      updateSummary(&param, &prior, &data, &summary); 
      ct++;
    }
    if((i+1) % _PRINT_FREQ_ == 0) {
      fprintf(stderr, "%d\t", i+1);
      /* fprintf(stderr, "%d (%.2f,%.2f,%.2f)\t", i+1, param.ptrue, param.beta0, param.betac);
      print_DP(&prior, &data); */
    }
  }
  fprintf(stderr, "\n");
  scaleSummary(&summary, &data, ct);
    calculateFDR(&data, &summary);

  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  chdir("..");

  /* output */
  write_result(&data, &summary);  
  return 0;
}		



