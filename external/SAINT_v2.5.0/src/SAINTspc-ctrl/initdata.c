#include "saint.h"

/*************************/
/* read interaction data */
/*************************/
void read_interaction_data(FILE *fpinter, DATA *data) {
  int i;
  char buf[100];

  data->ninter = nrow(fpinter);
  rewind(fpinter);

  assert(data->prey = (char **) calloc(data->ninter, sizeof(char *)));
  for(i=0;i<data->ninter;i++) assert(data->prey[i] = (char *) calloc(250, sizeof(char)));
  assert(data->bait = (char **) calloc(data->ninter, sizeof(char *)));
  for(i=0;i<data->ninter;i++) assert(data->bait[i] = (char *) calloc(250, sizeof(char)));
  assert(data->ip = (char **) calloc(data->ninter, sizeof(char *)));
  for(i=0;i<data->ninter;i++) assert(data->ip[i] = (char *) calloc(250, sizeof(char)));
  assert(data->d = (float *) calloc(data->ninter, sizeof(float)));
  assert(data->d2 = (float *) calloc(data->ninter, sizeof(float)));
  assert(data->iprob = (float *) calloc(data->ninter, sizeof(float)));
  assert(data->l = (float *) calloc(data->ninter, sizeof(float)));
  assert(data->c = (float *) calloc(data->ninter, sizeof(float)));

  assert(data->a2u = (int *) calloc(data->ninter, sizeof(int)));

  for(i=0;i<data->ninter;i++) {
    fscanf(fpinter, "%s", buf);
    strcpy(data->ip[i], buf);
    fscanf(fpinter, "%s", buf);
    strcpy(data->bait[i], buf);
    fscanf(fpinter, "%s", buf);
    strcpy(data->prey[i], buf);
    fscanf(fpinter, "%s", buf);
    data->d2[i] = atof(buf);
    if(data->d2[i] >= _TRUNC_) data->d[i] = _TRUNC_;
    else data->d[i] = data->d2[i];
    /* fprintf(stderr, "%s\t%s\t%s\t%f\n", data->ip[i], data->bait[i], data->prey[i], data->d[i]); */
  }
}


/***********************************************************************************************************/


/*********************************************************************************/
/* make unique interaction data and identify mapping between unique and all data */
/*********************************************************************************/
void find_unique_interaction(DATA *data) {
  int i,j,cur;
  int baitCompare, preyCompare;
  int isUnique[data->ninter];
  int nInstance[data->ninter]; /* this counts at the level of unique interactions */
  int counter[data->ninter];   /* same as above, used for mapping unique->individual */

  for(i=0;i<data->ninter;i++) {
    isUnique[i] = 1;
    nInstance[i] = 0;
    counter[i] = 0;
  }

  /* scan 1~n to mark unique interactions and count instances of each */
  cur = 0;
  for(i=0;i<(data->ninter-1);i++) {
    if(isUnique[i]) {
      (nInstance[cur])++;
      for(j=(i+1);j<data->ninter;j++) {
        if(isUnique[j]) {
          baitCompare = strcmp(data->bait[i], data->bait[j]);
          preyCompare = strcmp(data->prey[i], data->prey[j]);
          if(baitCompare == 0 && preyCompare == 0) {
            isUnique[j] = 0;
            (nInstance[cur])++;
          }
        }
      }
      cur++;
    }
  }
  if(isUnique[data->ninter-1]) {
    (nInstance[cur])++;
    cur++;
  }

  /* count # unique interactions */
  data->nuinter = 0;
  for(i=0;i<data->ninter;i++) {
    if(isUnique[i]) (data->nuinter)++;
  }

  /* memory business for unique interactions */
  assert(data->uprey = (char **) calloc(data->nuinter, sizeof(char *)));
  for(i=0;i<data->nuinter;i++) assert(data->uprey[i] = (char *) calloc(250, sizeof(char)));
  assert(data->ubait = (char **) calloc(data->nuinter, sizeof(char *)));
  for(i=0;i<data->nuinter;i++) assert(data->ubait[i] = (char *) calloc(250, sizeof(char)));
  assert(data->prob = (float *) calloc(data->nuinter, sizeof(float)));

  /* copy unique interactions */
  cur = 0;
  for(i=0;i<data->ninter;i++) {
    if(isUnique[i]) {
      strcpy(data->uprey[cur], data->prey[i]);
      strcpy(data->ubait[cur], data->bait[i]);
      data->prob[cur] = 0.0;
      cur++;
    }
  }
  if(data->nuinter > cur) fprintf(stderr, "Warning: possibly missed some unique interactions\n");
  else if(data->nuinter < cur) fprintf(stderr, "Warning: too many unique interactions, check mapping\n");
  else {}

  /* mapping between individual and unique interactions */
  assert(data->n_u2a = (int *) calloc(data->nuinter, sizeof(int)));
  assert(data->u2a = (int **) calloc(data->nuinter, sizeof(int *)));
  for(i=0;i<data->nuinter;i++) data->n_u2a[i] = nInstance[i];
  for(i=0;i<data->nuinter;i++) {
    assert(data->u2a[i] = (int *) calloc(data->n_u2a[i], sizeof(int)));
  }

  cur = 0; /* current index of unique */
  for(i=0;i<data->ninter;i++) {
    if(isUnique[i]) {
      data->a2u[i] = cur;
      data->u2a[cur][counter[cur]] = i;
      (counter[cur])++;
      for(j=(i+1);j<data->ninter;j++) {
        if(isUnique[j] == 0) {
          baitCompare = strcmp(data->bait[i], data->bait[j]);
          preyCompare = strcmp(data->prey[i], data->prey[j]);
          if(baitCompare == 0 && preyCompare == 0) {
            data->a2u[j] = cur;
            data->u2a[cur][counter[cur]] = j;
            (counter[cur])++;
          }
        }
      }
      cur++;
    }
  }
}

/***********************************************************************************************************/

/*****************************************************/
/* make indicators of uniqueness in character arrays */
/* returns the number of unique elements             */
/*****************************************************/
int unique_elements(char **x, int *unique, int nx) {
  int i,j;
  int nunique = nx;
  for(i=0;i<nx;i++) unique[i] = 1;
  for(i=0;i<(nx-1);i++) {
    if(unique[i]) {
	  for(j=(i+1);j<nx;j++) {
        if(strcmp(x[i], x[j]) == 0) {
		  unique[j] = 0;
          nunique--;
        }
      }
    }
  }
  return nunique;
}

int count_unique_elements(char **x, int nx) {
  int i,j;
  int unique[nx];
  int nunique = nx;
  for(i=0;i<nx;i++) unique[i] = 1;
  for(i=0;i<(nx-1);i++) {
    if(unique[i]) {
	  for(j=(i+1);j<nx;j++) {
        if(strcmp(x[i], x[j]) == 0) {
		  unique[j] = 0;
          nunique--;
        }
      }
    }
  }
  return nunique;
}

void centerData(float *x, int n, int takelog) {
  int i;
  float mean;
  if(takelog) {
    for(i=0;i<n;i++) x[i] = log(x[i] + 1.0);
  }
  mean = vec_mean(x, n);
  for(i=0;i<n;i++) x[i] -= mean;
}

/*********************************************************/
/* mapping between preys and interactions (unique & not) */
/*********************************************************/
int mapPreyToData(DATA *data) {
  int i,j;
  int cur;

  for(i=0;i<data->nprey;i++) data->preyNinter[i] = 0;

  assert(data->i2p = (int *) calloc(data->ninter, sizeof(int))); 
  for(i=0;i<data->ninter;i++) data->i2p[i] = -1;

  for(i=0;i<data->nprey;i++) {
    for(j=0;j<data->ninter;j++) {
      if(strcmp(data->PREY[i], data->prey[j]) == 0) {
        (data->preyNinter[i])++;
        data->i2p[j] = i;        
      }
    }
  }

  assert(data->p2i = (int **) calloc(data->nprey, sizeof(int *))); 
  for(i=0;i<data->nprey;i++) assert(data->p2i[i] = (int *) calloc(data->preyNinter[i], sizeof(int)));
  
  for(i=0;i<data->nprey;i++) {
    cur = 0;
    for(j=0;j<data->ninter;j++) {
      if(strcmp(data->PREY[i], data->prey[j]) == 0) {
        data->p2i[i][cur] = j;        
        cur++;
      }
      if(cur >= data->preyNinter[i]) break;
    }
  }

  assert(data->ui2p = (int *) calloc(data->nuinter, sizeof(int))); 
  for(i=0;i<data->nprey;i++) {
    for(j=0;j<data->nuinter;j++) {
      if(strcmp(data->PREY[i], data->uprey[j]) == 0) {
        data->ui2p[j] = i;        
      }
    }
  }
  
  /* report which prey in the prey file did not show up in the interaction file */  
  cur = 0;
  for(i=0;i<data->nprey;i++) {
    if(data->preyNinter[i] == 0) {
	  cur = 1;
      break;
    }
  }
  if(cur) {
    chdir("LOG");
    FILE *fptemp1 = fopen("PreysNotInData", "w");
    for(i=0;i<data->nprey;i++) {
      if(data->preyNinter[i] == 0) fprintf(fptemp1, "%s\n", data->PREY[i]);
    }
    fclose(fptemp1);
    chdir("..");
    return 1;
  }
  
  /* report which prey in the interaction file did not show up in the prey file */  
  cur = 0;
  for(i=0;i<data->ninter;i++) {
    if(data->i2p[i] == -1) {
	  cur = 1;
      break;
    }
  }
  if(cur) {
    chdir("LOG");
    FILE *fptemp2 = fopen("PreysNotInList", "w");
    for(i=0;i<data->ninter;i++) {
      if(data->i2p[i] == -1) fprintf(fptemp2, "%d\t%s\t%s\t%s\n", i+1, data->ip[i], data->bait[i], data->prey[i]);
    }
	fclose(fptemp2);
    chdir("..");
    return 1;
  }

  return 0;
}

/**************************************************************/
/* read prey data and check discrepancy with interaction data */
/**************************************************************/
void read_prey_data(FILE *fpprey, DATA *data) {
  int i, nprey;
  char buf[256];
  nprey = nrow(fpprey);
  rewind(fpprey);
  data->nprey = nprey;

  assert(data->PREY = (char **) calloc(nprey, sizeof(char *)));
  for(i=0;i<nprey;i++) assert(data->PREY[i] = (char *) calloc(250, sizeof(char)));
  assert(data->PREYGENE = (char **) calloc(nprey, sizeof(char *)));
  for(i=0;i<nprey;i++) assert(data->PREYGENE[i] = (char *) calloc(250, sizeof(char)));
  assert(data->preyLength = (float *) calloc(nprey, sizeof(float)));
  assert(data->preyNinter = (int *) calloc(nprey, sizeof(int)));
  assert(data->ctrlavg = (float *) calloc(nprey, sizeof(float)));

  for(i=0;i<nprey;i++) {
    fscanf(fpprey, "%s", buf);
    strcpy(data->PREY[i], buf);
    fscanf(fpprey, "%s", buf);
    data->preyLength[i] = atof(buf);
    /* fprintf(stderr, "%s\t%f\n", data->PREY[i], data->preyLength[i]); */
    data->ctrlavg[i] = 0.0;
    fscanf(fpprey, "%s", buf);
    strcpy(data->PREYGENE[i], buf);
  }
  centerData(data->preyLength, nprey, 1);
  mapPreyToData(data);
  
  for(i=0;i<data->ninter;i++) data->l[i] = data->preyLength[data->i2p[i]];
}



/***********************************************************************************************************/
void mapIPtoBait(DATA *data) {
  int i,j;
  int nbait, nIP, cur;
  char temp[data->nIP][256];
  int uniqueBaits[data->nIP];

  nIP = data->nIP;
  nbait = unique_elements(data->BAIT, uniqueBaits, nIP);
  data->nbait = nbait;  
  assert(data->baitNIP = (int *) calloc(nbait, sizeof(int)));
  for(i=0;i<data->nbait;i++) data->baitNIP[i] = 0;

  cur = 0;
  for(i=0;i<data->nIP;i++) {
     if(uniqueBaits[i]) {
       strcpy(temp[cur], data->BAIT[i]);
       cur++;
     }
  }
  if(cur != data->nbait) fprintf(stderr, "check bait-IP file\n");

  for(i=0;i<data->nbait;i++) {
    cur = 0;
    for(j=0;j<data->nIP;j++) {        
      if(strcmp(temp[i], data->BAIT[j]) == 0) {
        cur++;
      }
    }
    data->baitNIP[i] = cur;
  }

  assert(data->IP2b = (int *) calloc(data->nIP, sizeof(int)));
  assert(data->b2IP = (int **) calloc(data->nbait, sizeof(int *)));
  for(i=0;i<data->nbait;i++) assert(data->b2IP[i] = (int *) calloc(data->baitNIP[i], sizeof(int)));

  for(i=0;i<data->nbait;i++) {
    cur = 0;
    for(j=0;j<data->nIP;j++) {        
      if(strcmp(temp[i], data->BAIT[j]) == 0) {
        data->IP2b[j] = i;
        data->b2IP[i][cur] = j;
        cur++;
      }
    }
    data->baitNIP[i] = cur;
  }

  for(i=0;i<data->nbait;i++) strcpy(data->BAIT[i], temp[i]);

}


int mapIPBaitToData(DATA *data) {
  /* Part I: bait to data */
  int i,j;
  int cur;

  assert(data->baitNinter = (int *) calloc(data->nbait, sizeof(int)));

  for(i=0;i<data->nbait;i++) data->baitNinter[i] = 0;
  for(i=0;i<data->nIP;i++) data->IPNinter[i] = 0;

  assert(data->i2b = (int *) calloc(data->ninter, sizeof(int))); 
  for(i=0;i<data->ninter;i++) data->i2b[i] = -1;
  assert(data->i2IP = (int *) calloc(data->ninter, sizeof(int))); 
  for(i=0;i<data->ninter;i++) data->i2IP[i] = -1;

  for(i=0;i<data->nIP;i++) {
    for(j=0;j<data->ninter;j++) {
      if(strcmp(data->IP[i], data->ip[j]) == 0) {
        (data->IPNinter[i])++;
        data->i2IP[j] = i;        
      }
    }
  }

  for(i=0;i<data->nbait;i++) {
    for(j=0;j<data->ninter;j++) {
      if(strcmp(data->BAIT[i], data->bait[j]) == 0) {
        (data->baitNinter[i])++;
        data->i2b[j] = i;        
      }
    }
  }

  assert(data->IP2i = (int **) calloc(data->nIP, sizeof(int *))); 
  for(i=0;i<data->nIP;i++) assert(data->IP2i[i] = (int *) calloc(data->IPNinter[i], sizeof(int)));
  
  for(i=0;i<data->nIP;i++) {
    cur = 0;
    for(j=0;j<data->ninter;j++) {
      if(strcmp(data->IP[i], data->ip[j]) == 0) {
        data->IP2i[i][cur] = j;        
        cur++;
      }
      if(cur >= data->IPNinter[i]) break;
    }
  }

  assert(data->b2i = (int **) calloc(data->nbait, sizeof(int *))); 
  for(i=0;i<data->nbait;i++) assert(data->b2i[i] = (int *) calloc(data->baitNinter[i], sizeof(int)));
  
  for(i=0;i<data->nbait;i++) {
    cur = 0;
    for(j=0;j<data->ninter;j++) {
      if(strcmp(data->BAIT[i], data->bait[j]) == 0) {
        data->b2i[i][cur] = j;        
        cur++;
      }
      if(cur >= data->baitNinter[i]) break;
    }
  }

  /* from unique interactions to bait/IP */
  assert(data->ui2b = (int *) calloc(data->nuinter, sizeof(int))); 
  for(i=0;i<data->nbait;i++) {
    for(j=0;j<data->nuinter;j++) {
      if(strcmp(data->BAIT[i], data->ubait[j]) == 0) data->ui2b[j] = i;        
    }
  }
  
  /* report which bait/IP in the bait file did not show up in the interaction file */  
  cur = 0;
  for(i=0;i<data->nbait;i++) {
    if(data->IPNinter[i] == 0) {
	  cur = 1;
      break;
    }
  }
  if(cur) {
    chdir("LOG");
    FILE *fptemp1 = fopen("IPNotInData", "w");
    for(i=0;i<data->nIP;i++) {
      if(data->IPNinter[i] == 0) fprintf(fptemp1, "%s\t%s\n", data->IP[i], data->BAIT[i]);
    }
	fclose(fptemp1);
    chdir("..");
    return 1;
  }
  
  /* report which baits/IPs in the interaction file did not show up in the bait/IP file */  
  cur = 0;
  for(i=0;i<data->ninter;i++) {
    if(data->i2IP[i] == -1) {
	  cur = 1;
      break;
    }
  }
  if(cur) {
    chdir("LOG");
    FILE *fptemp2 = fopen("IPNotInList", "w");
    for(i=0;i<data->ninter;i++) {
      if(data->i2IP[i] == -1) fprintf(fptemp2, "%d\t%s\t%s\t%s\n", i+1, data->ip[i], data->bait[i], data->prey[i]);
    }
	fclose(fptemp2);
    chdir("..");
    return 1;
  }
  return 0;
}



void getIPinfo(DATA *data) {
  int i,j;
  int IPmatch, BAITmatch, PREYmatch;
  char buf[256];
  assert(data->IPbaitCoverage = (float *) calloc(data->nIP, sizeof(float)));
  assert(data->IPtotalAbundance = (float *) calloc(data->nIP, sizeof(float)));
  for(i=0;i<data->nIP;i++) {
    data->IPbaitCoverage[i] = 0.0;
    data->IPtotalAbundance[i] = 0.0;
  }

  for(i=0;i<data->nIP;i++) {
    strcpy(buf, data->BAIT[data->IP2b[i]]);
    for(j=0;j<data->ninter;j++) {
      IPmatch = strcmp(data->ip[j], data->IP[i]);
      BAITmatch = strcmp(data->bait[j], buf);
      PREYmatch = strcmp(data->prey[j], buf);
      if(IPmatch == 0) {
         data->IPtotalAbundance[i] += data->d[j];
         if(BAITmatch == 0 && PREYmatch == 0) data->IPbaitCoverage[i] = data->d[j] / data->preyLength[data->i2p[j]];
      }
    }
    /* if(data->IPbaitCoverage[i] == 0.0) {
      fprintf(stderr, "IP %s (bait %s) has no bait-bait interaction\n", data->IP[i], data->BAIT[data->IP2b[i]]);
    } */
  }


  /* for(i=0;i<data->nIP;i++) {
    fprintf(stderr, "%d: %s %s %.2f %.2f\n", i+1, data->IP[i], data->BAIT[data->IP2b[i]], data->IPbaitCoverage[i], data->IPtotalAbundance[i]);
  } */
  
}



/**************************************************************/
/* read bait data and check discrepancy with interaction data */
/**************************************************************/
void read_bait_data(FILE *fpbait, DATA *data) {
  int i, nbait, nIP;
  char buf[256];
  nIP = nrow(fpbait);
  rewind(fpbait);
  data->nIP = nIP;

  data->nctrl = 0;
  data->ntest = 0;

  assert(data->BAIT = (char **) calloc(nIP, sizeof(char *)));
  for(i=0;i<nIP;i++) assert(data->BAIT[i] = (char *) calloc(250, sizeof(char)));
  assert(data->IP = (char **) calloc(nIP, sizeof(char *)));
  for(i=0;i<nIP;i++) assert(data->IP[i] = (char *) calloc(250, sizeof(char)));
  assert(data->ctrl = (int *) calloc(nIP, sizeof(int)));
  assert(data->IPNinter = (int *) calloc(nIP, sizeof(int)));

  for(i=0;i<nIP;i++) {
    fscanf(fpbait, "%s", buf);
    strcpy(data->IP[i], buf);
    fscanf(fpbait, "%s", buf);
    strcpy(data->BAIT[i], buf);   /* not unique at this point */
    fscanf(fpbait, "%s", buf);
    if(buf[0] == 'C' || buf[0] == 'c') {
      data->ctrl[i] = 1;	/* note that control is marked as 1, test is as 0 */
      (data->nctrl)++;
    }
    else {
      data->ctrl[i] = 0;
      (data->ntest)++;
    }
    /* fprintf(stderr, "%s\t%s\t%d\n", data->IP[i], data->BAIT[i], data->ctrl[i]); */
  } 

  /* check whether IPs are unique or not */

  mapIPtoBait(data);
  nbait = data->nbait;
  mapIPBaitToData(data);
  getIPinfo(data);  /* bait coverage and total abundance */
  centerData(data->IPbaitCoverage, nIP, 1);
  centerData(data->IPtotalAbundance, nIP, 1);  /* these quantities are on log scale, mean centered now. */

  for(i=0;i<data->ninter;i++) {
    data->c[i] = data->IPtotalAbundance[data->i2IP[i]];
  }
}


/***********************************************************************************************************/

void set_ctrlavg(DATA *data) {
  int i;
  for(i=0;i<data->ninter;i++) {
    if(data->ctrl[data->i2IP[i]]) {
      /* if(data->ctrlavg[data->i2p[i]] < data->d[i]) data->ctrlavg[data->i2p[i]] = data->d[i]; */
      data->ctrlavg[data->i2p[i]] += data->d[i];
    }
  }
  for(i=0;i<data->nprey;i++) data->ctrlavg[i] /= ((float) data->nctrl);
}

/**************************************************************/
/*            master function for reading the data            */
/**************************************************************/
void read_data(FILE *fpinter, FILE *fpprey, FILE *fpbait, DATA *data) {
  read_interaction_data(fpinter, data);
  find_unique_interaction(data);
  read_prey_data(fpprey, data);
  read_bait_data(fpbait, data);
  
  /* make a function to filter out interactions with no matching preys and baits */

  set_ctrlavg(data);
}

