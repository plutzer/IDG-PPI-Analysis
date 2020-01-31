#include "saint.h"

int intersect(int *x, int *y, int *res, int nx, int ny) {
  int i,j;
  int n = 0;
  int found[nx];
  for(i=0;i<nx;i++) found[i] = 0;

  for(i=0;i<nx;i++) {
    for(j=0;j<ny;j++) {
      if(x[i] == y[j]) {            
        found[i] = 1;
        break;
      }
    }
  }

  for(i=0;i<nx;i++) {
    if(found[i]) {
      res[n] = x[i];
      n++; 
    }
  }
   
  return n;
}

void append(DATA *data) {
  int i,j,k,l;
  int ninter, nIP;
  int isMatch;
  int isCtrl;
  int common[1000000];
  char IP[500][500];  /* should be enough */
  FILE *fp = fopen("interaction.new","w");

  for(j=0;j<data->nbait;j++) {
    nIP = data->baitNIP[j];
    for(i=0;i<nIP;i++) strcpy(IP[i], data->IP[data->b2IP[j][i]]);
    isCtrl = 0;
    for(i=0;i<nIP;i++) {
      if(data->ctrl[data->b2IP[j][i]] == 1) isCtrl = 1;
    }
    /* check if this is test bait or control bait */

    if(isCtrl == 0) {

      for(i=0;i<data->nprey;i++) {


        ninter = intersect(data->b2i[j], data->p2i[i], common, data->baitNinter[j], data->preyNinter[i]);
        if(ninter > 0) {

          for(k=0;k<nIP;k++) {
            isMatch = -1;
            for(l=0;l<ninter;l++) {
              if(strcmp(IP[k], data->ip[common[l]]) == 0) {
                isMatch = common[l];
                break;
              }
            }
            if(isMatch >= 0) {
               if(strcmp(data->BAIT[j], data->bait[isMatch]) != 0) fprintf(stderr, "match error - bait\n");
               if(strcmp(data->PREY[i], data->prey[isMatch]) != 0) fprintf(stderr, "match error - prey\n");

               if(data->type) fprintf(fp, "%s\t%s\t%s\t%f\n", data->ip[isMatch], data->bait[isMatch], data->prey[isMatch], data->d[isMatch]);
               else fprintf(fp, "%s\t%s\t%s\t%d\n", data->ip[isMatch], data->bait[isMatch], data->prey[isMatch], (int) data->d[isMatch]);

            }
            else {
              if(strcmp(data->BAIT[j], data->PREY[i]) != 0) {
                if(data->type) fprintf(fp, "%s\t%s\t%s\t%.1f\n", IP[k], data->BAIT[j], data->PREY[i], 0.0);
                else fprintf(fp, "%s\t%s\t%s\t%d\n", IP[k], data->BAIT[j], data->PREY[i], 0);
              }
            }
          }
        }

      }

    }
    else {
      for(i=0;i<data->nprey;i++) {
        ninter = intersect(data->b2i[j], data->p2i[i], common, data->baitNinter[j], data->preyNinter[i]);
        for(k=0;k<nIP;k++) {
          isMatch = -1;
          for(l=0;l<ninter;l++) {
            if(strcmp(IP[k], data->ip[common[l]]) == 0) {
              isMatch = common[l];
              break;
            }
          }
          if(isMatch >= 0) {
            if(strcmp(data->BAIT[j], data->bait[isMatch]) != 0) fprintf(stderr, "match error - bait\n");
            if(strcmp(data->PREY[i], data->prey[isMatch]) != 0) fprintf(stderr, "match error - prey\n");
            if(strcmp(data->bait[isMatch], data->prey[isMatch]) != 0) {
              if(data->type) fprintf(fp, "%s\t%s\t%s\t%f\n", data->ip[isMatch], data->bait[isMatch], data->prey[isMatch], data->d[isMatch]);
              else fprintf(fp, "%s\t%s\t%s\t%d\n", data->ip[isMatch], data->bait[isMatch], data->prey[isMatch], (int) data->d[isMatch]);
            }
          }
          else {
            if(strcmp(data->BAIT[j], data->PREY[i]) != 0) {
              if(data->type) fprintf(fp, "%s\t%s\t%s\t%.1f\n", IP[k], data->BAIT[j], data->PREY[i], 0.0);
              else fprintf(fp, "%s\t%s\t%s\t%d\n", IP[k], data->BAIT[j], data->PREY[i], 0);
            }
          }
        }
      }
    }
  }  
  fclose(fp);
}


