#include "saint.h"

int takeMIN(int a, int b) {
  int res = a;
  if(b < a) res = b;
  return res;
}

/*************************/
/* read interaction data */
/*************************/
int read_all_data(FILE *fpinter, FILE *fpprey, FILE *fpbait, DATA *data) {
  int i,j;
  char buf[10000];
  int nIP, nprey, nuIP, nuprey, ninter, nuinter, cur;
  int *uPreyCount;
  int *uPreyLen;
  char **uniquePrey;
  char **uniquePreyGene;
  char **uniqueIP;
  int *uniqueInter;
  char **allInter;
  double resid;
  int *prey_appear;

  char **t_prey;
  char **t_bait;
  char **t_ip;
  double *t_d;

  char **t_preyname;
  int *t_preylen;

  /******************************/
  /* read interaction text file */
  /******************************/
  ninter = nrow(fpinter);
  rewind(fpinter);
  data->ninter = ninter;

  assert(allInter = (char **) calloc(ninter, sizeof(char *)));
  for(i=0;i<data->ninter;i++) assert(allInter[i] = (char *) calloc(500, sizeof(char)));
  assert(uniqueInter = (int *) calloc(ninter, sizeof(int)));
  /* assert(uniqueInter = (char **) calloc(ninter, sizeof(char *)));
  for(i=0;i<data->ninter;i++) assert(uniqueInter[i] = (char *) calloc(500, sizeof(char))); */

  assert(data->prey = (char **) calloc(data->ninter, sizeof(char *)));
  for(i=0;i<data->ninter;i++) assert(data->prey[i] = (char *) calloc(500, sizeof(char)));
  assert(data->bait = (char **) calloc(data->ninter, sizeof(char *)));
  for(i=0;i<data->ninter;i++) assert(data->bait[i] = (char *) calloc(500, sizeof(char)));
  assert(data->ip = (char **) calloc(data->ninter, sizeof(char *)));
  for(i=0;i<data->ninter;i++) assert(data->ip[i] = (char *) calloc(500, sizeof(char)));
  assert(data->d = (double *) calloc(data->ninter, sizeof(double)));
  assert(data->iprob = (double *) calloc(data->ninter, sizeof(double)));

  assert(data->a2u = (int *) calloc(data->ninter, sizeof(int)));

  assert(t_prey = (char **) calloc(data->ninter, sizeof(char *)));
  for(i=0;i<data->ninter;i++) assert(t_prey[i] = (char *) calloc(500, sizeof(char)));
  assert(t_bait = (char **) calloc(data->ninter, sizeof(char *)));
  for(i=0;i<data->ninter;i++) assert(t_bait[i] = (char *) calloc(500, sizeof(char)));
  assert(t_ip = (char **) calloc(data->ninter, sizeof(char *)));
  for(i=0;i<data->ninter;i++) assert(t_ip[i] = (char *) calloc(500, sizeof(char)));
  assert(t_d = (double *) calloc(data->ninter, sizeof(double)));

  for(i=0;i<data->ninter;i++) {
    fscanf(fpinter, "%s", buf);
    strcpy(data->ip[i], buf);
    strcpy(allInter[i], buf);
    strcat(allInter[i], " ");
    fscanf(fpinter, "%s", buf);
    strcpy(data->bait[i], buf);
    fscanf(fpinter, "%s", buf);
    strcpy(data->prey[i], buf);
    strcat(allInter[i], buf);
    fscanf(fpinter, "%s", buf);
    data->d[i] = atof(buf);
  }

  data->type = 0;
  for(i=0;i<data->ninter;i++) {
    resid = data->d[i] - ((double) ((int) data->d[i]));
    if(resid > 0.0) {
      data->type = 1;
      break;
    }
  }

  nuinter = unique_elements(allInter, uniqueInter, ninter);
  if(nuinter != ninter) {
    fprintf(stderr, "Duplicate Warning: IP-prey pair must be unique in the interaction file...");
    cur = 0;
    for(i=0;i<data->ninter;i++) {
      if(uniqueInter[i]) {
        strcpy(t_prey[cur], data->prey[i]);
        strcpy(t_bait[cur], data->bait[i]);
        strcpy(t_ip[cur], data->ip[i]);
        t_d[cur] = data->d[i];
        cur++;
      }
    }
    data->ninter = nuinter;
    for(i=0;i<data->ninter;i++) {
        strcpy(data->prey[i], t_prey[i]);
        strcpy(data->bait[i], t_bait[i]);
        strcpy(data->ip[i], t_ip[i]);
        data->d[i] = t_d[i];
    }
    fprintf(stderr, "fixed.\n");
  }

  for(i=0;i<ninter;i++) free(allInter[i]);
  free(allInter);
  free(uniqueInter);


  /***********************/
  /* read prey text file */
  /***********************/
  rewind(fpprey);
  nprey = nrow(fpprey);
  rewind(fpprey);
  data->nprey = nprey;    
  assert(uPreyCount = (int *) calloc(nprey, sizeof(int)));
  assert(uPreyLen = (int *) calloc(nprey, sizeof(int)));
  assert(data->PREY = (char **) calloc(nprey, sizeof(char *)));
  for(i=0;i<nprey;i++) assert(data->PREY[i] = (char *) calloc(500, sizeof(char)));
  assert(data->PREYGENE = (char **) calloc(nprey, sizeof(char *)));
  for(i=0;i<nprey;i++) assert(data->PREYGENE[i] = (char *) calloc(500, sizeof(char)));
  assert(uniquePrey = (char **) calloc(nprey, sizeof(char *)));
  for(i=0;i<nprey;i++) assert(uniquePrey[i] = (char *) calloc(500, sizeof(char)));
  assert(uniquePreyGene = (char **) calloc(nprey, sizeof(char *)));
  for(i=0;i<nprey;i++) assert(uniquePreyGene[i] = (char *) calloc(500, sizeof(char)));
  assert(data->preyLen = (int *) calloc(nprey, sizeof(int)));

  for(i=0;i<nprey;i++) uPreyCount[i] = 0; 
  
  for(i=0;i<nprey;i++) {
    fscanf(fpprey, "%s", buf);
    if(strlen(buf) > 500) {
      fprintf(stderr, "Prey name %s is longer than 500 characters.\n", buf);
      return 1;
    }
    strcpy(data->PREY[i], buf);

    if(data->type == 0) {
      fscanf(fpprey, "%s", buf);
      data->preyLen[i] = atoi(buf);   /* not unique at this point */
    }

    fscanf(fpprey, "%s", buf);
    strcpy(data->PREYGENE[i], buf);
  } 

  /* if preys are not unique, then write out prey file and reread. */
  nuprey = unique_elements_copy(data->PREY, uniquePrey, nprey);
  if(nuprey != nprey) {
    fprintf(stderr, "Duplicate Warning: Preys must be unique in the prey file...fixed.\n");
  }

  for(i=0;i<nprey;i++) { 
    for(j=0;j<nuprey;j++) {
      if(strcmp(data->PREY[i], uniquePrey[j]) == 0) {
        (uPreyCount[j])++;
        break;
      }
    }
  }     

  for(j=0;j<nuprey;j++) { 
    for(i=0;i<nprey;i++) {
      if(strcmp(uniquePrey[j], data->PREY[i]) == 0) {
        if(data->type == 0) {
          uPreyLen[j] = data->preyLen[i];
        }
        strcpy(uniquePreyGene[j], data->PREYGENE[i]);
        break;
      }
    }
  }    

  data->nprey = nuprey;
  assert(prey_appear = (int *) calloc(nuprey, sizeof(int)));
  for(i=0;i<nuprey;i++) {
    prey_appear[i] = 0;
    for(j=0;j<data->ninter;j++) {
      if(strcmp(uniquePrey[i], data->prey[j]) == 0) {
        prey_appear[i] = 1;
        break;
      }
    }
  }

  assert(t_preyname = (char **) calloc(nuprey, sizeof(char *)));
  for(i=0;i<nuprey;i++) assert(t_preyname[i] = (char *) calloc(500, sizeof(char)));
  assert(t_preylen = (int *) calloc(nuprey, sizeof(int)));

  FILE *fpp = fopen("prey.new", "w");
  cur = 0;
  for(j=0;j<nuprey;j++) {
    if(prey_appear[j]) {
      if(data->type) fprintf(fpp, "%s\t%s\n", uniquePrey[j], uniquePreyGene[j]);
      else fprintf(fpp, "%s\t%d\t%s\n", uniquePrey[j], uPreyLen[j], uniquePreyGene[j]);
      strcpy(t_preyname[cur], uniquePrey[j]);
      if(data->type == 0) t_preylen[cur] = uPreyLen[j];
      cur++;
    }
  }
  fclose(fpp);

  data->nprey = cur;
  for(i=0;i<cur;i++) {
    strcpy(data->PREY[i], t_preyname[i]);
    if(data->type == 0) data->preyLen[i] = t_preylen[i];    
  }

  free(prey_appear);

  for(i=0;i<nprey;i++) {
    free(uniquePrey[i]);
  }
  free(uniquePrey);
  free(uPreyCount);
  free(uPreyLen);
  
  /***********************/
  /* read bait text file */
  /***********************/
  nIP = nrow(fpbait);
  rewind(fpbait);
  data->nIP = nIP;

  data->nctrl = 0;
  data->ntest = 0;

  assert(data->BAIT = (char **) calloc(nIP, sizeof(char *)));
  for(i=0;i<nIP;i++) assert(data->BAIT[i] = (char *) calloc(500, sizeof(char)));
  assert(data->IP = (char **) calloc(nIP, sizeof(char *)));
  for(i=0;i<nIP;i++) assert(data->IP[i] = (char *) calloc(500, sizeof(char)));
  assert(uniqueIP = (char **) calloc(nIP, sizeof(char *)));
  for(i=0;i<nIP;i++) assert(uniqueIP[i] = (char *) calloc(500, sizeof(char)));
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
  } 

  /* if IP names are duplicated, return 1 */
  nuIP = unique_elements_copy(data->IP, uniqueIP, nIP);
  if(nuIP != nIP) {
    fprintf(stderr, "Duplicate Warning: IP names must be unique in the bait file. User must resolve this.\n");
    return 1;
  }

  for(i=0;i<nIP;i++) free(uniqueIP[i]);
  free(uniqueIP);


  return 0;
}

/******************************************************/
/* Reformatting the data: prey map, control reduction */
/******************************************************/
  
int reformat_data(DATA *data) {
  int i, j, id, cur;
  double spec[data->nprey][data->_K_];
  double allspec[data->nctrl];
  int ctrlCounts[data->nprey];

  // fprintf(stderr, "%d\n", data->nctrl);

  /* reformat control data: do this only if(K > 5) */
  for(i=0;i<data->nprey;i++) {
    for(j=0;j<data->_K_;j++) spec[i][j] = 0.0;
  }    
  if(data->nctrl > data->_K_) {  
    for(i=0;i<data->nprey;i++) {
      cur = 0;
      for(j=0;j<data->nctrl;j++) allspec[j] = 0;
      for(j=0;j<data->preyNinter[i];j++) {
        id = data->p2i[i][j];
        // fprintf(stderr, "%s %s\n", data->bait[id], data->prey[id]);
        if(data->ctrl[data->i2IP[id]]) { 
          allspec[cur] = data->d[id];
          cur++;
        }
      }
      ctrlCounts[i] = cur;
      if(cur > 0) {
        gsl_sort(allspec, 1, data->nctrl);
        for(j=0;j<data->nctrl;j++) {
           spec[i][j] = allspec[data->nctrl-1-j];    
           if(j >= data->_K_) break;
        }
      }
    }
  }

  /*   write interaction_intermediate */
  FILE *fpi = fopen("interaction.intermediate", "w");

  if(data->nctrl > data->_K_) {
    for(i=0;i<data->ninter;i++) {
      if(data->ctrl[data->i2IP[i]] == 0) {
        if(data->type) fprintf(fpi, "%s\t%s\t%s\t%f\n", data->ip[i], data->bait[i], data->prey[i], data->d[i]);
        else fprintf(fpi, "%s\t%s\t%s\t%d\n", data->ip[i], data->bait[i], data->prey[i], ((int) data->d[i]));
      }
    }
    for(i=0;i<data->nprey;i++) {
      if(ctrlCounts[i] > 0) {
	for(j=0;j<data->_K_;j++) {
          if(spec[i][j] > 0.0) {
            if(data->type) fprintf(fpi, "CTRL%d\tCTRL%d\t%s\t%.3f\n", j+1, j+1, data->PREY[i], spec[i][j]);
            else fprintf(fpi, "CTRL%d\tCTRL%d\t%s\t%d\n", j+1, j+1, data->PREY[i], (int) spec[i][j]);
          }
        }
      }
    }
  }
  else {
    for(i=0;i<data->ninter;i++) {
      if(data->type) fprintf(fpi, "%s\t%s\t%s\t%f\n", data->ip[i], data->bait[i], data->prey[i], data->d[i]);
      else fprintf(fpi, "%s\t%s\t%s\t%d\n", data->ip[i], data->bait[i], data->prey[i], ((int) data->d[i]));
    }
  }
  fclose(fpi);


  /*   write bait.new */
  FILE *fpb = fopen("bait.new", "w");
  if(data->nctrl > data->_K_) {
    for(i=0;i<data->nIP;i++) {    
      if(data->ctrl[i] == 0) {
        fprintf(fpb, "%s\t%s\tT\n", data->IP[i], data->BAIT[data->IP2b[i]]);
      }
    }
    for(j=0;j<data->_K_;j++) fprintf(fpb, "CTRL%d\tCTRL%d\tC\n", j+1, j+1);
  }
  else {
    for(i=0;i<data->nIP;i++) {
      fprintf(fpb, "%s\t%s\t%s\n", data->IP[i], data->BAIT[data->IP2b[i]], data->ctrl[i] ? "C" : "T");
    }
  }
  fclose(fpb);
  
  return 0;

}






