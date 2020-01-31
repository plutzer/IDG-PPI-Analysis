#include "saint.h"

/**************************************************/
/*             Writing Out the Results            */
/**************************************************/

void write_interactions(DATA *data, SUMMARY *summary) {
  int i;
  FILE *fp = fopen("interactions", "w");
  fprintf(fp, "IP\tBait\tPrey\tPreyGene\tIntensity\tProb\tiProb\n");
  for(i=0;i<data->ninter;i++) {
    if(data->ctrl[data->i2IP[i]] == 0 && data->miss[i] == 0) fprintf(fp, "%s\t%s\t%s\t%s\t%.3f\t%.3f\t%.3f\n", 
                   data->ip[i], data->bait[i], data->prey[i], data->PREYGENE[data->i2p[i]], data->d[i], 
				   summary->Z[data->a2u[i]], summary->iZ[i]);
  }
  fclose(fp);
}

void write_unique_interactions(DATA *data, SUMMARY *summary) {
  int i,j;
  int isCtrl;
  int id, ct, pid;
  float intsum;
  float maxp, avgp, geop, tmp;
  FILE *fp = fopen("unique_interactions", "w");
  fprintf(fp, "Bait\tPrey\tPreyGene\tIP\tIntensity\tIntensitySum\tNumRep\tProb\tiProb\tctrlIntensity\tAvgP\tMaxP\tFDR\n");
  
  for(i=0;i<data->nuinter;i++) {
    isCtrl = 0;
    for(j=0;j<data->n_u2a[i];j++) {
      id = data->u2a[i][j];
      if(data->ctrl[data->i2IP[id]]) isCtrl = 1;
    }
    if(isCtrl == 0) {
      fprintf(fp, "%s\t%s\t%s\t", data->ubait[i], data->uprey[i], data->PREYGENE[data->ui2p[i]]);

      for(j=0;j<(data->n_u2a[i]-1);j++) fprintf(fp, "%s|", data->ip[data->u2a[i][j]]);
      fprintf(fp, "%s\t", data->ip[data->u2a[i][data->n_u2a[i]-1]]);  

      intsum = 0.0;
      for(j=0;j<(data->n_u2a[i]-1);j++) {
        id = data->u2a[i][j];
        if(data->miss[id]==0) intsum += exp(data->d[id]);

        if(data->miss[id]==0) fprintf(fp, "%.3f", exp(data->d[id]));
        else fprintf(fp, ".");
        fprintf(fp, "|");
      }
      id = data->u2a[i][data->n_u2a[i]-1];
      if(data->miss[id]==0) intsum += exp(data->d[id]);

      if(data->miss[id]==0) fprintf(fp, "%.3f", exp(data->d[id]));  
      else fprintf(fp, ".");
      fprintf(fp, "\t");

      fprintf(fp, "%.3f\t%d\t", intsum, data->n_u2a[i]);

      fprintf(fp, "%.2f\t", summary->Z[i]);

      for(j=0;j<(data->n_u2a[i]-1);j++) {
        id = data->u2a[i][j];
        fprintf(fp, "%.2f", summary->iZ[id]);
        fprintf(fp, "|");
      }
      id = data->u2a[i][data->n_u2a[i]-1];
      fprintf(fp, "%.2f", summary->iZ[id]);  
      fprintf(fp, "\t");

      pid = data->ui2p[i];
      ct = 0;
      for(j=0;j<data->preyNinter[pid];j++) {
        id = data->p2i[pid][j];
        if(data->ctrl[data->i2IP[id]]) {
          if(data->miss[id]==0) fprintf(fp, "%.2f", exp(data->d[id]));
          else fprintf(fp, ".");

          if(ct < data->nctrl-1) fprintf(fp, "|");
          else fprintf(fp, "\t");
          ct++;
        }
      }

      /* maxp */
      maxp = 0.0;
      avgp = 0.0;
      geop = 0.0;
      for(j=0;j<data->n_u2a[i];j++) {
        id = data->u2a[i][j];
        if(summary->iZ[id] > maxp) maxp = summary->iZ[id];
        avgp += summary->iZ[id] / ((float) data->n_u2a[i]);
        tmp = data->d[id] == 0.0 ? 0.001 : summary->iZ[id];
        geop += log(tmp) / ((float) data->n_u2a[i]);
      } 
      geop = exp(geop);
      fprintf(fp, "%.4f\t%.4f\t%.4f\n", avgp, maxp, summary->FDR[i]);
      /* fprintf(fp, "%.2f\n", ((float) data->preyNinter[data->ui2p[i]]) / ((float) data->nIP)); */
    }
  }

  fclose(fp);
}

void write_prey(DATA *data, SUMMARY *summary) {
  int i;
  FILE *fp = fopen("preys", "w");
  fprintf(fp, "Prey\tAlpha_prey\tMu\n");
  for(i=0;i<data->nprey;i++) {
    fprintf(fp, "%s\t%.2f\t%.2f\n", data->PREY[i], summary->alpha_prey[i], summary->mu[i]);
  }
  fclose(fp);
}

void write_IP(DATA *data, SUMMARY *summary) {
  int i;
  FILE *fp = fopen("IPs", "w");
  fprintf(fp, "IP\tBait\tAlpha_IP\n");
  for(i=0;i<data->nIP;i++) {
    fprintf(fp, "%s\t%s\t%.2f\n", data->IP[i], data->BAIT[data->IP2b[i]], summary->alpha_IP[i]);
  }
  fclose(fp);
}

void write_bait(DATA *data, SUMMARY *summary) {
  int i,j;
  FILE *fp = fopen("baits", "w");
  fprintf(fp, "Bait\tIP\tAlpha_IP\n");
  for(i=0;i<data->nbait;i++) {
    fprintf(fp, "%s\t", data->BAIT[i]);    

    for(j=0;j<data->baitNIP[i]-1;j++) {
      fprintf(fp, "%s|", data->IP[data->b2IP[i][j]]);
    }
    fprintf(fp, "%s\t", data->IP[data->b2IP[i][data->baitNIP[i]-1]]);
    
    for(j=0;j<data->baitNIP[i]-1;j++) {
      fprintf(fp, "%.2f|", summary->alpha_IP[data->b2IP[i][j]]);
    }
    fprintf(fp, "%.2f\n", summary->alpha_IP[data->b2IP[i][data->baitNIP[i]-1]]);
  }
  fclose(fp);
}

void write_histogram(FILE *fp, HISTOGRAM *hist) {
  int i;
  fprintf(fp, "-inf\t%.2f\t%.2f\n", hist->start[0], hist->count[0]);
  for(i=0;i<_HISTO_BIN_;i++) {
    fprintf(fp, "%.2f\t%.2f\t%.2f\n", hist->start[i], hist->end[i], hist->count[i+1]);
  }
  fprintf(fp, "%.2f\tinf\t%.2f\n", hist->end[_HISTO_BIN_-1], hist->count[_HISTO_BIN_+1]);
}

void write_histogram2(FILE *fp, HISTOGRAM2 *hist) {
  int i;
  fprintf(fp, "-inf\t%.2f\t%.2f\n", hist->start[0], hist->count[0]);
  for(i=0;i<_HISTO_BIN2_;i++) {
    fprintf(fp, "%.2f\t%.2f\t%.2f\n", hist->start[i], hist->end[i], hist->count[i+1]);
  }
  fprintf(fp, "%.2f\tinf\t%.2f\n", hist->end[_HISTO_BIN2_-1], hist->count[_HISTO_BIN2_+1]);
}


void write_hyperprior(DATA *data, SUMMARY *summary) {
  FILE *fp1 = fopen("hist_alpha_prey", "w");
  FILE *fp2 = fopen("hist_alpha_IP", "w");
  FILE *fp3 = fopen("hist_mu", "w");
  FILE *fp4 = fopen("hist_eta", "w");
  FILE *fp5 = fopen("hist_eta0", "w");

  write_histogram(fp1, &(summary->hist_alpha_prey));
  write_histogram(fp2, &(summary->hist_alpha_IP));
  write_histogram(fp3, &(summary->hist_mu));
  write_histogram2(fp4, &(summary->hist_eta));
  write_histogram2(fp5, &(summary->hist_eta0));

  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  fclose(fp4);
  fclose(fp5);
}  

void write_result(DATA *data, SUMMARY *summary) {
  chdir("RESULT");
  write_interactions(data, summary);
  write_unique_interactions(data, summary);
  write_prey(data, summary);
  write_IP(data, summary);
  write_bait(data, summary);
  write_hyperprior(data, summary);  
  write_matrix_data(data, summary);
  write_matrix_data2(data, summary);
  chdir("..");
}

/******************************/

void write_matrix_data(DATA *data, SUMMARY *summary) {
  int i,j,k,id;
  int endLine, isMatch;
  FILE *fp = fopen("matrix_form","w"); 
  endLine = -1;
  for(j=0;j<data->nIP;j++) {
    if(data->ctrl[j] != 1) endLine = j;
  }

  /* header line 1 */
  fprintf(fp, "Bait\t");
  for(j=0;j<data->nIP;j++) {
    id = data->IP2b[j];
    if(data->ctrl[j]) {
	  fprintf(fp, "%s\t", data->BAIT[id]);
    }
  }
  for(j=0;j<data->nIP;j++) {
    id = data->IP2b[j];
    if(data->ctrl[j] == 0) {
      fprintf(fp, "%s", data->BAIT[id]);
      if(j==endLine) fprintf(fp, "\n");
      else fprintf(fp, "\t");
    }
  }

  /* header line 2 */
  fprintf(fp, "IP\t");
  for(j=0;j<data->nIP;j++) {
    if(data->ctrl[j]) {
      fprintf(fp, "%s\t", data->IP[j]);
    }
  }
  for(j=0;j<data->nIP;j++) {
    if(data->ctrl[j] == 0) {
	  fprintf(fp, "%s", data->IP[j]);
      if(j==endLine) fprintf(fp, "\n");
      else fprintf(fp, "\t");
    }
  }

  /* header line 3 */
  fprintf(fp, "Prey\t");
  for(j=0;j<data->nIP;j++) {
    if(data->ctrl[j]) {
	  fprintf(fp, "\t");
    }
  }
  for(j=0;j<data->nIP;j++) {
    if(data->ctrl[j] == 0) {
      fprintf(fp, "%.2f", summary->alpha_IP[j]);
      if(j==endLine) fprintf(fp, "\n");
      else fprintf(fp, "\t");
    }
  }

  for(i=0;i<data->nprey;i++) {
    fprintf(fp, "%s\t", data->PREY[i]);
    /* Control runs first */ 
    for(j=0;j<data->nIP;j++) {
      if(data->ctrl[j]) {
        /* find if prey wise data has this IP */
        /* if not, leave the space blank */
        /* else, biz as usual: (count | prob | lambda_s, lambda_ns) */
        isMatch = -1;
        for(k=0;k<data->preyNinter[i];k++) {
          id = data->p2i[i][k];
          if(strcmp(data->ip[id], data->IP[j]) == 0) isMatch = id;
        }
        if(isMatch == -1) {
          fprintf(fp, "\t");
        }
        else {
          if(data->miss[isMatch]==0) fprintf(fp, "%.2f|%.2f\t", data->d[isMatch], (summary->lambda_false[isMatch]));
          else fprintf(fp, "\t");
        }
      }
    }
    /* Rest of the IPs */
    for(j=0;j<data->nIP;j++) {
      if(data->ctrl[j] == 0) {
        isMatch = -1;
        for(k=0;k<data->preyNinter[i];k++) {
	  id = data->p2i[i][k];
	  if(strcmp(data->ip[id], data->IP[j]) == 0) isMatch = id;
        }
        if(isMatch == -1) {
          /* fprintf(fp, "\t"); */
        }
        else {
	  fprintf(fp, "%.2f|%.2f|%.2f|%.2f", data->d[isMatch], summary->Z[data->a2u[isMatch]], (summary->lambda_true[isMatch]), (summary->lambda_false[isMatch]));
        }

        if(j==endLine) fprintf(fp, "\n");
        else fprintf(fp, "\t");
      }
    }
  }

  fclose(fp);
}



void write_matrix_data2(DATA *data, SUMMARY *summary) {
  int i,j,k,id;
  int endLine, isMatch;
  FILE *fp = fopen("matrix_form_short","w"); 
  endLine = -1;
  for(j=0;j<data->nIP;j++) {
    if(data->ctrl[j] != 1) endLine = j;
  }

  /* header line 1 */
  fprintf(fp, "Bait\t");
  for(j=0;j<data->nIP;j++) {
    id = data->IP2b[j];
    if(data->ctrl[j]) {
	  fprintf(fp, "%s\t", data->BAIT[id]);
    }
  }
  for(j=0;j<data->nIP;j++) {
    id = data->IP2b[j];
    if(data->ctrl[j] == 0) {
      fprintf(fp, "%s", data->BAIT[id]);
      if(j==endLine) fprintf(fp, "\n");
      else fprintf(fp, "\t");
    }
  }

  /* header line 2 */
  fprintf(fp, "IP\t");
  for(j=0;j<data->nIP;j++) {
    if(data->ctrl[j]) {
      fprintf(fp, "%s\t", data->IP[j]);
    }
  }
  for(j=0;j<data->nIP;j++) {
    if(data->ctrl[j] == 0) {
	  fprintf(fp, "%s", data->IP[j]);
      if(j==endLine) fprintf(fp, "\n");
      else fprintf(fp, "\t");
    }
  }

  /* header line 3 */
  fprintf(fp, "Prey\t");
  for(j=0;j<data->nIP;j++) {
    if(data->ctrl[j]) {
	  fprintf(fp, "\t");
    }
  }
  for(j=0;j<data->nIP;j++) {
    if(data->ctrl[j] == 0) {
      fprintf(fp, "%.2f", summary->alpha_IP[j]);
      if(j==endLine) fprintf(fp, "\n");
      else fprintf(fp, "\t");
    }
  }

  for(i=0;i<data->nprey;i++) {
    fprintf(fp, "%s\t", data->PREY[i]);
    /* Control runs first */ 
    for(j=0;j<data->nIP;j++) {
      if(data->ctrl[j]) {
        /* find if prey wise data has this IP */
        /* if not, leave the space blank */
        /* else, biz as usual: (count | prob | lambda_s, lambda_ns) */
        isMatch = -1;
        for(k=0;k<data->preyNinter[i];k++) {
          id = data->p2i[i][k];
          if(strcmp(data->ip[id], data->IP[j]) == 0) isMatch = id;
        }
        if(isMatch == -1) {
          fprintf(fp, "\t");
        }
        else {
          fprintf(fp, "%.2f\t", data->d[isMatch]);
        }
      }
    }
    /* Rest of the IPs */
    for(j=0;j<data->nIP;j++) {
      if(data->ctrl[j] == 0) {
        isMatch = -1;
        for(k=0;k<data->preyNinter[i];k++) {
	  id = data->p2i[i][k];
	  if(strcmp(data->ip[id], data->IP[j]) == 0) isMatch = id;
        }
        if(isMatch == -1) {
          /* fprintf(fp, "\t"); */
        }
        else {
	  fprintf(fp, "(%.2f|%.2f|%.2f)", data->d[isMatch], summary->Z[data->a2u[isMatch]], summary->iZ[isMatch]);
        }

        if(j==endLine) fprintf(fp, "\n");
        else fprintf(fp, "\t");
      }
    }
  }

  fclose(fp);
}






