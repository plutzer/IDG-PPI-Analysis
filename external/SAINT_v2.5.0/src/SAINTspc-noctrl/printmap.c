#include "saint.h"

void printInter(DATA *data) {
  int i;
  FILE *fp = fopen("interaction","w");
  fprintf(fp, "ip\tbait\tprey\tIP\tBAIT\tPREY\tubait\tuprey\n");
  for(i=0;i<data->ninter;i++) {
    fprintf(fp, "%s\t%s\t%s\t", data->ip[i], data->bait[i], data->prey[i]);  
    fprintf(fp, "%s\t%s\t%s\t", data->IP[data->i2IP[i]], data->BAIT[data->i2b[i]], data->PREY[data->i2p[i]]);  
    fprintf(fp, "%s\t%s\n", data->ubait[data->a2u[i]], data->uprey[data->a2u[i]]);  
  }
  fclose(fp);
}

void printUInter(DATA *data) {
  int i,j,k;
  FILE *fp = fopen("unique_interaction","w");
  fprintf(fp, "ubait\tuprey\tubait\tuprey\tip\tbait\tprey\n");
  for(i=0;i<data->nuinter;i++) {
    for(j=0;j<data->n_u2a[i];j++) {
      k = data->u2a[i][j];
      fprintf(fp, "%s\t%s\t", data->ubait[i], data->uprey[i]);  
      fprintf(fp, "%s\t%s\t%s\n", data->ip[k], data->bait[k], data->prey[k]);  
    }
  }

  fprintf(fp, "\n\n************************\n\n");
  
  fprintf(fp, "ubait\tuprey\tBAIT\tPREY\n");
  for(i=0;i<data->nuinter;i++) {
    fprintf(fp, "%s\t%s\t%s\t%s\n", data->ubait[i], data->uprey[i], data->BAIT[data->ui2b[i]], data->PREY[data->ui2p[i]]);
  }

  fclose(fp);
}

void printIP(DATA *data) {
  int i,j,k;
  /* IP to bait */
  FILE *fp = fopen("IP","w");
  fprintf(fp, "IP\tBAIT\n");
  for(i=0;i<data->nIP;i++) {
    fprintf(fp, "%s\t%s\n", data->IP[i], data->BAIT[data->IP2b[i]]);
  }

  fprintf(fp, "\n\n************************\n\n");
  
  /* IP to interactions */
  fprintf(fp, "IP\tip\tbait\tprey\n");
  for(i=0;i<data->nIP;i++) {
    for(j=0;j<data->IPNinter[i];j++) {
      k = data->IP2i[i][j];
      fprintf(fp, "%s\t%s\t%s\t%s\n", data->IP[i], data->ip[k], data->bait[k], data->prey[k]);
    }
  }
  fclose(fp);
}

void printBait(DATA *data) {
  int i,j,k;
  FILE *fp = fopen("bait","w");

  /* bait to IP */
  fprintf(fp, "BAIT\tIP\n");
  for(i=0;i<data->nbait;i++) {
    for(j=0;j<data->baitNIP[i];j++) {
      k = data->b2IP[i][j];
      fprintf(fp, "%s\t%s\n", data->BAIT[i], data->IP[k]); 
    }
  }

  fprintf(fp, "\n\n************************\n\n");

  /* bait to interaction */
  fprintf(fp, "BAIT\tip\tbait\tprey\n");
  for(i=0;i<data->nbait;i++) {
    for(j=0;j<data->baitNinter[i];j++) {
      k = data->b2i[i][j];
      fprintf(fp, "%s\t%s\t%s\t%s\n", data->BAIT[i], data->ip[k], data->bait[k], data->prey[k]);
    }
  }
  fclose(fp);
}

void printPrey(DATA *data) {
  int i,j,k;
  FILE *fp = fopen("prey","w");
  /* prey to interaction */
  for(i=0;i<data->nprey;i++) {
    for(j=0;j<data->preyNinter[i];j++) {
      k = data->p2i[i][j];
      fprintf(fp, "%s\t%s\t%s\t%s\n", data->PREY[i], data->ip[k], data->bait[k], data->prey[k]);
    }
  }
  fclose(fp);
}

void printMap(DATA *data) {
  chdir("MAPPING");
  printInter(data);
  printUInter(data);
  printIP(data);
  printBait(data);
  printPrey(data);
  chdir("..");
}


