#include "saint.h"


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

  /* count # unique interactions */
  data->nuinter = 0;
  for(i=0;i<data->ninter;i++) {
    if(isUnique[i]) (data->nuinter)++;
  }

  /* memory business for unique interactions */
  assert(data->uprey = (char **) calloc(data->nuinter, sizeof(char *)));
  for(i=0;i<data->nuinter;i++) assert(data->uprey[i] = (char *) calloc(500, sizeof(char)));
  assert(data->ubait = (char **) calloc(data->nuinter, sizeof(char *)));
  for(i=0;i<data->nuinter;i++) assert(data->ubait[i] = (char *) calloc(500, sizeof(char)));
  assert(data->prob = (double *) calloc(data->nuinter, sizeof(double)));

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

int unique_elements_copy(char **x, char **uniq, int nx) {
  int i,j;
  int nunique;
  int unique[nx];

  nunique = nx;
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
  j = 0;
  for(i=0;i<nx;i++) {
    if(unique[i]) {
      strcpy(uniq[j], x[i]);
      j++;
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


/*********************************************************/
/* mapping between preys and interactions (unique & not) */
/*********************************************************/
int mapPreyToData(DATA *data) {
  int i,j;
  int cur;

  assert(data->preyNinter = (int *) calloc(data->nprey, sizeof(int)));
  for(i=0;i<data->nprey;i++) data->preyNinter[i] = 0;

  assert(data->i2p = (int *) calloc(data->ninter, sizeof(int))); 
  for(i=0;i<data->ninter;i++) data->i2p[i] = -1;

  for(j=0;j<data->ninter;j++) {
    for(i=0;i<data->nprey;i++) {
      if(strcmp(data->PREY[i], data->prey[j]) == 0) {
        (data->preyNinter[i])++;
        data->i2p[j] = i;        
        break;
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
      if(cur >= data->preyNinter[i]) {
	break;
      }
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
    chdir("reformat_log");
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
    chdir("reformat_log");
    FILE *fptemp2 = fopen("PreysNotInList", "w");
    fprintf(stderr, "Some prey(s) are missing from the prey file.\nCheck PreysNotInList file in reformat_log folder and insert them in the prey file.\nThere may be duplicates in this list.\n");
    for(i=0;i<data->ninter;i++) {
      if(data->i2p[i] == -1) fprintf(fptemp2, "%d\t%s\n", i+1, data->prey[i]);
    }
    fclose(fptemp2);
    chdir("..");
    return 1;
  }

  return 0;
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
    for(j=0;j<data->nIP;j++) {
      if(strcmp(temp[i], data->BAIT[j]) == 0) (data->baitNIP[i])++;
    }
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
  for(i=0;i<data->nIP;i++) {
    if(data->IPNinter[i] == 0) {
	  cur = 1;
      break;
    }
  }
  if(cur) {
    chdir("reformat_log");
    FILE *fptemp1 = fopen("IPNotInData", "w");
    for(i=0;i<data->nIP;i++) {
      if(data->IPNinter[i] == 0) fprintf(fptemp1, "%s\t%s\n", data->IP[i], data->BAIT[data->IP2b[i]]);
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
    chdir("reformat_log");
    FILE *fptemp2 = fopen("IPNotInList", "w");
    for(i=0;i<data->ninter;i++) {
      if(data->i2IP[i] == -1) fprintf(fptemp2, "%d\t%s\n", i+1, data->ip[i]);
    }
    fclose(fptemp2);
    chdir("..");
    return 1;
  }
  return 0;
}


/**************************************************************/
/* read bait data and check discrepancy with interaction data */
/**************************************************************/

/***********************************************************************************************************/




