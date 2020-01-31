#include "saint.h"

int nrow(FILE *fp) {
  char buf[100000];
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
  char buf[100000];
  int i,cont = 0;
  fgets(buf, sizeof(buf), fp);
  for(i=0;i<100000;i++) {
    if(buf[i] == '\t') cont++;
    if(buf[i] == '\0') break;
  }
  return cont;
}

int main(int argc, char **argv) {

  int b, i, j, k, l, m, p, q, burn, iter, total, tmp, u, tmpmaxint, count;
  double ff;  
  /* Data */

  char prob[_MAX_BUF_];
  DATA data; 
  PARAM param;
  PRIOR prior;
  SUMMARY summary;
  char fbait[_MAX_BUF_];
  char fprey[_MAX_BUF_];
  char fmu[_MAX_BUF_];
  char list[_MAX_BUF_];
  char iprob[_MAX_BUF_];
        
  const gsl_rng_type *T;
  gsl_rng *r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  param.useAbun = 1;
  param.useLen = 0;
  param.useCov = 1;
  prior.useAbun = 1;
  prior.useLen = 0;
  prior.useCov = 1;
  
  if ((argc < 6) || (argc > 9)) {
    fprintf(stderr, "usage: saint-spc-noctrl-matrix [interactomeData] [output file] [nburnin] [niter] [ff]\n       saint-spc-noctrl-matrix [interactomeData] [output file] [nburnin] [niter] [ff] [useAbundance(0/1)] [useLength(0/1)] [useCoverage(0/1)]\n");
    return 1;
  }  
  FILE *fpi = fopen(argv[1], "r");
  p = nrow(fpi)-3;
  q = ncol(fpi)-3;
  fclose(fpi);

  fprintf(stderr, "%d proteins, %d IPs\n", p, q);

  FILE *fp = fopen(argv[1], "r");
  FILE *fp_output = fopen(argv[2], "w");
  strcpy(prob, argv[2]);
  strcat(prob, "_prob");
  FILE *fp_outprob = fopen(prob, "w");
  strcpy(list, argv[2]);
  strcat(list, "_list");
  FILE *fp_list = fopen(list, "w");
  strcpy(fbait, argv[2]);
  strcat(fbait, "_alpha_bait");
  FILE *fpbait = fopen(fbait, "w");
  strcpy(fprey, argv[2]);
  strcat(fprey, "_alpha_prey");
  FILE *fpprey = fopen(fprey, "w");
  strcpy(fmu, argv[2]);
  strcat(fmu, "_mu_prey");
  FILE *fpmu = fopen(fmu, "w");
  strcpy(iprob, argv[2]);
  strcat(iprob, "_iprob");
  FILE *fp_iprob = fopen(iprob, "w");

  burn = atoi(argv[3]);
  iter = atoi(argv[4]);
  ff = atof(argv[5]);
  if(argc > 7) {
    tmp = atoi(argv[6]);
    if(&tmp != NULL) {
      param.useAbun = tmp;
      prior.useAbun = tmp;
    }
    tmp = atoi(argv[7]);
    if(&tmp != NULL) {
      param.useLen = tmp;
      prior.useLen = tmp;
    }
    tmp = atoi(argv[8]);
    if(&tmp != NULL) {
      param.useCov = tmp;
      prior.useCov = tmp;
    }  
  }
  
  if(fp == NULL) { 
    fprintf(stderr, "Interactome data %s does not exist.\n", argv[1]);
    return 1; 
  }
  if(&p == NULL) { 
    fprintf(stderr, "The number of prey proteins was not provided.\n");
    return 1; 
  }
  if(&q == NULL) { 
    fprintf(stderr, "The number of bait proteins was not provided.\n");
    return 1; 
  }
  if(&burn == NULL) { 
    fprintf(stderr, "The number of burnin iterations was not provided.\n");
    return 1; 
  }
  if(&iter == NULL) { 
    fprintf(stderr, "The number of main interations was not provided.\n");
    return 1; 
  }

/*  fprintf(stderr, "\n*************************************************\n");
  fprintf(stderr, "************     Welcome to SAInt     ***********\n");
  fprintf(stderr, "*** Significance Analysis of Interactome Data ***\n");
  fprintf(stderr, "*************************************************\n");
  
  if(param.useAbun) fprintf(stderr, "Background abundance is used in this model\n");
  if(param.useLen) fprintf(stderr, "Sequence length is used in this model\n");
  if(param.useCov) fprintf(stderr, "Bait Coverage is used in this model\n");
  fprintf(stderr, "*************************************************\n");
*/    
  /* Read Data */
  read_data(fp, &data, &p, &q);
  init_prior(&prior, &p, &q);
  init_param(&data, &param, &p, &q);
  param.ff_prop = ff;
  init_summary(&data, &summary, &p, &q);
  set_summary(&data, &summary);
  set_prior(&prior);
  set_param(&param, &prior, &data, iter, r);
  u = data.uniqueNum;  

  for(j=0;j<(p-1);j++) fprintf(fpprey, "%s\t", data.prey[j]);
  fprintf(fpprey, "%s\n", data.prey[p-1]);
  for(j=0;j<(q-1);j++) fprintf(fpbait, "%s\t", data.bait[j]);
  fprintf(fpbait, "%s\n", data.bait[q-1]);
  for(j=0;j<(p-1);j++) fprintf(fpmu, "%s\t", data.prey[j]);
  fprintf(fpmu, "%s\n", data.prey[p-1]);

  /* Estimation */
  fprintf(stderr, "Burnin:\n");
  for(i=0;i<burn;i++) {
    if((i+1) % _PRINT_FREQ_ == 0) {
      fprintf(stderr, "%d\t", i+1);
      /* fprintf(stderr, "C-%.3f R-%.3f NR-%.3f\n", param.pcont[0], param.preal[1], param.preal[0]); */
    }
    mhgibbs(&data, &prior, &param, r);
  }
  fprintf(stderr, "\n");
  fprintf(stderr, "Iterations: \n");
  total = 0;
    
  for(i=0;i<iter;i++) {
    if((i+1) % _PRINT_FREQ_ == 0) {
     fprintf(stderr, "%d\t", i+1);
     /* fprintf(stderr, "C-%.3f R-%.3f NR-%.3f\n", param.pcont[0], param.preal[1], param.preal[0]); */
    }
    mhgibbs(&data, &prior, &param, r);
    if((i+1) % _SKIP_ == 0) {
      for(j=0;j<p;j++) {
        summary.Y[j] += ((double) param.Y[j]);
        for(b=0;b<data.ninterRowUnique[j];b++) {
          k = data.useRowUnique[j][b];
          summary.Z[j][k] += ((double) param.Y[j]) * ((double) param.Z[j][k]);
        }
        for(b=0;b<data.ninterRow[j];b++) {
          k = data.useRow[j][b];
          summary.iZ[j][k] += ((double) param.Y[j]) * ((double) param.iZ[j][k]);
        } 
      }
      total++;
      for(j=0;j<(p-1);j++) fprintf(fpprey, "%.2f\t", param.alpha_prey[j]);
      fprintf(fpprey, "%.2f\n", param.alpha_prey[p-1]);
      for(j=0;j<(q-1);j++) fprintf(fpbait, "%.2f\t", param.alpha_bait[j]);
      fprintf(fpbait, "%.2f\n", param.alpha_bait[q-1]);
      for(j=0;j<(p-1);j++) fprintf(fpmu, "%.2f\t", param.mu_prey[j]);
      fprintf(fpmu, "%.2f\n", param.mu_prey[p-1]);
    }
  }
  for(j=0;j<p;j++) {
    summary.Y[j] /= ((double) total);
    for(b=0;b<data.ninterRowUnique[j];b++) {
      k = data.useRowUnique[j][b];
      summary.Z[j][k] /= ((double) total);
    }
    for(b=0;b<data.ninterRow[j];b++) {
      k = data.useRow[j][b];
      summary.iZ[j][k] /= ((double) total);
    } 
    summary.max_prob[j] = vec_max(summary.Z[j], u);
    param.flagged[j] /= ((double) (iter+burn));
  }

  fprintf(fp_list, "Experiment\tBait\tPrey\tSpec\tProb\tiProb\tContaminant\n");
  for(j=0;j<data.uniqueNum;j++) {
    for(b=0;b<data.ninterUnique[j];b++) {
      l = data.useUnique[j][b];
      for(k=0;k<data.uniqueSize[j];k++) {
        m = data.mfu[j][k];
        fprintf(fp_list, "%s\t%s\t%s\t%d\t%.3f\t%.3f\t%s\n", 
                        data.experiment[m], data.bait[m], data.prey[l],
                        (int) data.d[l][m], summary.Z[l][j], summary.iZ[l][m],
                        ((1.0 - summary.Y[l]) >= ff ? "yes" : "no")); 
      }
    }
  }

  /* Output */
  fprintf(stderr,"\nOutput summary..");
      
  fprintf(stderr, "..");  

  /* First line */
  fprintf(fp_outprob, "\t");
  for(j=0;j<(u-1);j++) fprintf(fp_outprob, "%s\t", data.unique[j]);
  fprintf(fp_outprob, "%s\n", data.unique[u-1]);

  /* First line */
  fprintf(fp_iprob, "\t");
  for(j=0;j<(q-1);j++) fprintf(fp_iprob, "%s\t", data.experiment[j]);
  fprintf(fp_iprob, "%s\n", data.experiment[q-1]);

  /* First line */
  fprintf(fp_output, "\t");
  for(j=0;j<(q-1);j++) fprintf(fp_output, "%s\t", data.experiment[j]);
  fprintf(fp_output, "%s\n", data.experiment[q-1]);

  /* Second line */
  fprintf(fp_output, "\t");
  for(j=0;j<(q-1);j++) fprintf(fp_output, "%s\t", data.bait[j]);
  fprintf(fp_output, "%s\n", data.bait[q-1]);

  for(j=0;j<p;j++) {
    count = 0;
    for(l=0;l<data.uniqueNum;l++) {
      tmpmaxint = 0;
      for(k=0;k<data.uniqueSize[l];k++) {
        if(((int) data.d[j][data.mfu[l][k]]) > tmpmaxint) tmpmaxint = ((int) data.d[j][data.mfu[l][k]]);
      }
      if(tmpmaxint >= 1) count++;
    }
    fprintf(fp_output, "%s\t", data.prey[j]);

    fprintf(fp_outprob, "%s\t", data.prey[j]);
    for(k=0;k<(u-1);k++) {
      fprintf(fp_outprob, "%.3f\t", summary.Z[j][k]);    
    }
    for(k=0;k<(q-1);k++) {
      fprintf(fp_output, "%d|%.3f|%.3f\t", ((int) data.d[j][k]), summary.Z[j][data.mtu[k]],summary.iZ[j][k]);
    }
    fprintf(fp_output, "%d|%.3f|%.3f\n", ((int) data.d[j][k]), summary.Z[j][data.mtu[q-1]], summary.iZ[j][q-1]);
    fprintf(fp_outprob, "%.3f\n", summary.Z[j][u-1]);
  }

  /* Matrix Output */
  for(j=0;j<p;j++) {
    fprintf(fp_iprob, "%s\t", data.prey[j]);
    for(k=0;k<(q-1);k++) {
      fprintf(fp_iprob, "%.3f\t", summary.iZ[j][k]);
    }
    fprintf(fp_iprob, "%.3f\n", summary.iZ[j][q-1]);
  }

  fprintf(stderr, "..\nThank you for using SAInt!\n");

  free_data(&data);
  free_prior(&prior);
  free_param(&param);
  free_summary(&summary);
  fclose(fp);
  fclose(fp_output);
  fclose(fp_outprob);
  fclose(fpbait);
  fclose(fpprey);
  fclose(fpmu);
  fclose(fp_iprob);
  fclose(fp_list);
  return 0;
}




