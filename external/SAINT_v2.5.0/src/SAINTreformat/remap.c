#include "saint.h"

/***********************************************************************************************************/

void reread_data(DATA *data) {
  int i;
  char buf[10000];
  int nIP, nprey, ninter;


  /***********************/
  /* read prey text file */
  /***********************/
  FILE *fpprey = fopen("prey.new", "r");
  nprey = nrow(fpprey);
  rewind(fpprey);
  data->nprey = nprey;    
  assert(data->PREY = (char **) calloc(nprey, sizeof(char *)));
  for(i=0;i<nprey;i++) assert(data->PREY[i] = (char *) calloc(500, sizeof(char)));
  assert(data->preyLen = (int *) calloc(nprey, sizeof(int)));
  assert(data->PREYGENE = (char **) calloc(nprey, sizeof(char *)));
  for(i=0;i<nprey;i++) assert(data->PREYGENE[i] = (char *) calloc(500, sizeof(char)));
 
  for(i=0;i<nprey;i++) {

    fscanf(fpprey, "%s", buf);
    if(strlen(buf) > 500) {
      fprintf(stderr, "Prey name %s is longer than 500 characters.\n", buf);
    }
    strcpy(data->PREY[i], buf);
    if(data->type == 0) {
      fscanf(fpprey, "%s", buf);
      data->preyLen[i] = atoi(buf);   /* not unique at this point */
    }
    fscanf(fpprey, "%s", buf);
    strcpy(data->PREYGENE[i], buf);
  } 
  fclose(fpprey);


  
  /***********************/
  /* read bait text file */
  /***********************/
  FILE *fpbait = fopen("bait.new", "r");
  nIP = nrow(fpbait);
  rewind(fpbait);
  data->nIP = nIP;
  data->nbait = data->nIP;

  data->nctrl = 0;
  data->ntest = 0;

  assert(data->BAIT = (char **) calloc(nIP, sizeof(char *)));
  for(i=0;i<nIP;i++) assert(data->BAIT[i] = (char *) calloc(500, sizeof(char)));
  assert(data->IP = (char **) calloc(nIP, sizeof(char *)));
  for(i=0;i<nIP;i++) assert(data->IP[i] = (char *) calloc(500, sizeof(char)));
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

  /******************************/
  /* read interaction text file */
  /******************************/
  FILE *fpinter = fopen("interaction.intermediate", "r");
  ninter = nrow(fpinter);
  rewind(fpinter);
  data->ninter = ninter;
  data->nuinter = data->ninter;
  assert(data->prey = (char **) calloc(data->ninter, sizeof(char *)));
  for(i=0;i<data->ninter;i++) assert(data->prey[i] = (char *) calloc(250, sizeof(char)));
  assert(data->bait = (char **) calloc(data->ninter, sizeof(char *)));
  for(i=0;i<data->ninter;i++) assert(data->bait[i] = (char *) calloc(250, sizeof(char)));
  assert(data->ip = (char **) calloc(data->ninter, sizeof(char *)));
  for(i=0;i<data->ninter;i++) assert(data->ip[i] = (char *) calloc(250, sizeof(char)));
  assert(data->d = (double *) calloc(data->ninter, sizeof(double)));
  assert(data->iprob = (double *) calloc(data->ninter, sizeof(double)));

  assert(data->a2u = (int *) calloc(data->ninter, sizeof(int)));

  for(i=0;i<data->ninter;i++) {
    fscanf(fpinter, "%s", buf);
    strcpy(data->ip[i], buf);
    fscanf(fpinter, "%s", buf);
    strcpy(data->bait[i], buf);
    fscanf(fpinter, "%s", buf);
    strcpy(data->prey[i], buf);
    fscanf(fpinter, "%s", buf);
    data->d[i] = atof(buf);
  }

  fclose(fpinter);
}


/*
void remap_data(DATA *data) {
  int map_well;
  find_unique_interaction(data);
  map_well = mapPreyToData(data);
  mapIPtoBait(data);
  mapIPBaitToData(data);
}


*/
