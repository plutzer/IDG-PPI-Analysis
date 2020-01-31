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
    if(buf[i] == '\t' || buf[i] == ' ') cont++;
    if(buf[i] == '\0') break;
  }
  return cont;
}



int commandLine(DATA *data, int argc, char **argv) {
  if (argc == 2 && strcmp(argv[1], "-v") == 0) {
    fprintf(stderr, "SAINT version 2.5.0\n");
    return 1;
  }

  if (!(argc == 4 || argc == 5)) {
    fprintf(stderr, "usage: saint-spc-reformat [interactionFile] [preyFile] [baitFile]\n");
    fprintf(stderr, "usage: saint-spc-reformat [interactionFile] [preyFile] [baitFile] [# Control]\n");
    return 1;
  }

  /* interaction file: IPnumber \t bait \t prey \t spectralCount \n */
  /* prey file:        prey \t sequenceLength \n */
  /* bait file:        IPnumber \t bait \t isControl \n */

  FILE *fpinter = fopen(argv[1], "r");
  FILE *fpprey = fopen(argv[2], "r");
  FILE *fpbait = fopen(argv[3], "r");
  if(argc == 5) data->_K_ = atoi(argv[4]);
  else data->_K_ = 100;

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
  fclose(fpprey);
  fclose(fpbait);
  fclose(fpinter);
  return 0;
}




/**************************************************************/
/*            master function for reading the data            */
/**************************************************************/
int read_data(FILE *fpinter, FILE *fpprey, FILE *fpbait, DATA *data, DATA *newdata) {
  int read_well, map_well;

  read_well = read_all_data(fpinter, fpprey, fpbait, data);

  if(read_well != 0) return 1;

  find_unique_interaction(data);
  map_well = mapPreyToData(data);
  mapIPtoBait(data);
  mapIPBaitToData(data);

  read_well = reformat_data(data); /* here I take the maximum K counts */

  newdata->type = data->type;
  reread_data(newdata);

  find_unique_interaction(newdata);
  map_well = mapPreyToData(newdata);
  mapIPtoBait(newdata);
  mapIPBaitToData(newdata);

  system("rm -rf interaction.intermediate");

  return 0;
}





/***************************** MAIN ***************************/

int main(int argc, char **argv) {
  int progress;
  DATA data; 
  DATA newdata;
  
  /* Command Line */
  if(commandLine(&data, argc, argv)) return 1;
  FILE *fpinter = fopen(argv[1], "r");
  FILE *fpprey = fopen(argv[2], "r");
  FILE *fpbait = fopen(argv[3], "r");

  /* Read interaction data, identify baits, preys, and IPs, 
     make unique interaction data frame, 
     identify the mapping between different levels of data */	

  system("mkdir reformat_log");      /* mapping logs */
  progress = read_data(fpinter, fpprey, fpbait, &data, &newdata);
  if(progress != 0) return 1;
  /* printMap(&data); */
  /* we can also check if preys are only from controls here */

 
  append(&newdata);


  fclose(fpinter);
  fclose(fpprey);
  fclose(fpbait);
  return 0;
}		



