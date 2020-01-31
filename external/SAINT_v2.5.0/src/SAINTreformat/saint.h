
/*
    Copyright (C) <2011>  <Hyungwon Choi>
    For troubleshooting, contact hyung_won_choi@nuhs.edu.sg.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You can obtain a copy of the GNU General Public License from
    <http://www.gnu.org/licenses/>.
*/




#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <assert.h>
#include <gsl/gsl_sort.h>

#define _MAX_BUF_ 500

typedef struct tagDATA {
  int _K_;
  int type; /* 0: count data; 1: intensity data */

  /*************/
  /* logistics */
  /*************/
  int ninter;
  int nuinter;
  int nprey;
  int nIP;
  int nbait;

  /**************************/
  /* interaction level data */
  /**************************/
  char **prey;
  char **bait;
  char **ip;   /* raw data, each row corresponds to one interaction, case-sensitive */
  double *d;
  double *iprob;

  /*********************************/
  /* unique interaction level data */
  /*********************************/
  char **uprey;
  char **ubait;
  double *prob;

  int *n_u2a; /* number of individual interactions per unique interactions */
  int **u2a;  /* unique interactions to individual interactions */  
  int *a2u;   /* individual interactions to unique interactions */ 
  /* crucial indicator for probability calculation */

  /***********************************/
  /* unique bait and prey level data */
  /***********************************/

  char **PREY;  /* unique preys */
  char **PREYGENE;
  char **BAIT;  /* unique baits */
  char **IP;    /* unique IP #s */

  int nctrl;
  int ntest;
  int *ctrl;  /* index: control IPs or not: 'C' = control, 'T' = test */

  int *preyNinter;  /* # interaction for prey */
  int *baitNinter;  /* # interaction for bait */
  int *IPNinter;    /* # interaction in an IP */
  int *baitNIP;     /* # IPs per bait         */
  int *preyLen;

  /****************/
  /* mapping data */
  /****************/
  int *i2p;   /* index: interaction to prey */
  int *i2b;   /* index: interaction to bait */
  int *i2IP;  /* index: interaction to IP   */ 

  int **p2i;  /* index: prey to interaction */
  int **b2i;  /* index: bait to interaction */
  int **IP2i; /* index: IP to interaction   */ 

  int *ui2p;   /* index: unique interaction to prey */
  int *ui2b;   /* index: unique interaction to bait */
  /* no need to build reverse mapping for unique interactions */
  /* perhaps this mapping is unnecessary */

  int **b2IP; /* index: bait to IP */
  int *IP2b;  /* index: IP to bait */

} DATA;


/*************/
/* functions */
/*************/

int takeMIN(int a, int b);

int nrow(FILE *fp);
int newlinechar(char *buf, int k);
int ncol(FILE *fp);
int commandLine(DATA *data, int argc, char **argv);

int read_all_data(FILE *fpinter, FILE *fpprey, FILE *fpbait, DATA *data);
void find_unique_interaction(DATA *data);
int unique_elements(char **x, int *unique, int nx);
int unique_elements_copy(char **x, char **uniq, int nx);
int count_unique_elements(char **x, int nx);
int mapPreyToData(DATA *data);
void prey_data(DATA *data);
void mapIPtoBait(DATA *data);
int mapIPBaitToData(DATA *data);
void map_bait_data(DATA *data);
int read_data(FILE *fpinter, FILE *fpprey, FILE *fpbait, DATA *data, DATA *newdata);
int reformat_data(DATA *data);
void reread_data(DATA *newdata);

/* new addition of files here for remapping */


void find_unique_interaction_anew(DATA *data);
void mapPreyToData_anew(DATA *data);
void mapIPtoBait_anew(DATA *data);
void mapIPBaitToData_anew(DATA *data);
void remap_data(DATA *data);


void printInter(DATA *data);
void printUInter(DATA *data);
void printIP(DATA *data);
void printBait(DATA *data);
void printPrey(DATA *data);
void printMap(DATA *data);

int intersect(int *x, int *y, int *res, int nx, int ny);
void append(DATA *data);



