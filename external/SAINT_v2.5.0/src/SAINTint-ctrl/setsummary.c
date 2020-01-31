#include "saint.h"

/**************************************************************/
/*             initializing the model summaryeters            */
/**************************************************************/

void memory_summary(SUMMARY *summary, DATA *data) {
  assert(summary->iZ = (float *) calloc(data->ninter, sizeof(float)));
  assert(summary->Z = (float *) calloc(data->nuinter, sizeof(float)));
  assert(summary->alpha_prey = (float *) calloc(data->nprey, sizeof(float)));
  assert(summary->alpha_IP = (float *) calloc(data->nIP, sizeof(float)));
  assert(summary->mu = (float *) calloc(data->nprey, sizeof(float)));
  assert(summary->eta = (float *) calloc(data->nprey, sizeof(float)));
  assert(summary->eta0 = (float *) calloc(data->nprey, sizeof(float)));
  assert(summary->lambda_true = (float *) calloc(data->ninter, sizeof(float)));
  assert(summary->lambda_false = (float *) calloc(data->ninter, sizeof(float)));
  assert(summary->FDR = (float *) calloc(data->nuinter, sizeof(float)));
}

void initialize_summary(SUMMARY *summary, DATA *data) {
  int i;
  for(i=0;i<data->ninter;i++) summary->iZ[i] = 0.0;
  for(i=0;i<data->nuinter;i++) summary->Z[i] = 0.0;
  for(i=0;i<data->nprey;i++) summary->alpha_prey[i] = 0.0;
  for(i=0;i<data->nIP;i++) summary->alpha_IP[i] = 0.0;
  for(i=0;i<data->nprey;i++) summary->mu[i] = 0.0;
  for(i=0;i<data->nprey;i++) summary->eta[i] = 0.0;
  for(i=0;i<data->nprey;i++) summary->eta0[i] = 0.0;
  for(i=0;i<data->ninter;i++) summary->lambda_true[i] = 0.0;
  for(i=0;i<data->ninter;i++) summary->lambda_false[i] = 0.0;
  for(i=0;i<data->nuinter;i++) summary->FDR[i] = 0.0;
}

void initialize_histogram(HISTOGRAM *hist) {
  int i;
  float binsize = ((float) (_HISTO_END_ - _HISTO_START_)) / ((float) _HISTO_BIN_);
  for(i=0;i<_HISTO_BIN_;i++) {
    hist->start[i] = _HISTO_START_ + ((float) i) * binsize;
    hist->end[i] = _HISTO_START_ + ((float) (i+1)) * binsize;
  }
  for(i=0;i<(_HISTO_BIN_+2);i++) hist->count[i] = 0.0;
}

void initialize_histogram2(HISTOGRAM2 *hist) {
  int i;
  float binsize = ((float) (_HISTO_END2_ - _HISTO_START2_)) / ((float) _HISTO_BIN2_);
  for(i=0;i<_HISTO_BIN2_;i++) {
    hist->start[i] = _HISTO_START2_ + ((float) i) * binsize;
    hist->end[i] = _HISTO_START2_ + ((float) (i+1)) * binsize;
  }
  for(i=0;i<(_HISTO_BIN2_+2);i++) hist->count[i] = 0.0;
}

void set_summary(SUMMARY *summary, DATA *data) {
  memory_summary(summary, data);
  initialize_summary(summary, data);
  initialize_histogram(&(summary->hist_alpha_prey));
  initialize_histogram(&(summary->hist_alpha_IP));
  initialize_histogram(&(summary->hist_mu));
  initialize_histogram2(&(summary->hist_eta)); 
  initialize_histogram2(&(summary->hist_eta0)); 
}

void updateSummary(PARAM *param, PRIOR *prior, DATA *data, SUMMARY *summary) {
  int i;
  for(i=0;i<data->ninter;i++) summary->iZ[i] += ((float) param->iZ[i]);
  for(i=0;i<data->nuinter;i++) summary->Z[i] += ((float) param->Z[i]);
  for(i=0;i<data->nprey;i++) summary->alpha_prey[i] += param->alpha_prey[i];
  for(i=0;i<data->nIP;i++) summary->alpha_IP[i] += param->alpha_IP[i];
  for(i=0;i<data->nprey;i++) summary->mu[i] += param->mu[i];
  for(i=0;i<data->nprey;i++) summary->eta[i] += param->eta[i];
  for(i=0;i<data->nprey;i++) summary->eta0[i] += param->eta0[i];
  for(i=0;i<data->ninter;i++) summary->lambda_true[i] += param->lambda_true[i];
  for(i=0;i<data->ninter;i++) summary->lambda_false[i] += param->lambda_false[i];
  updateHistogram(param, prior, data, summary);
}

void scaleSummary(SUMMARY *summary, DATA *data, int iter) {
  int i;
  float scale = 1.0 / ((float) iter);
  for(i=0;i<data->ninter;i++) summary->iZ[i] *= scale;
  for(i=0;i<data->nuinter;i++) summary->Z[i] *= scale;
  for(i=0;i<data->nprey;i++) summary->alpha_prey[i] *= scale;
  for(i=0;i<data->nIP;i++) summary->alpha_IP[i] *= scale;
  for(i=0;i<data->nprey;i++) summary->mu[i] *= scale;
  for(i=0;i<data->nprey;i++) summary->eta[i] *= scale;
  for(i=0;i<data->nprey;i++) summary->eta0[i] *= scale;
  for(i=0;i<data->ninter;i++) summary->lambda_true[i] *= scale;
  for(i=0;i<data->ninter;i++) summary->lambda_false[i] *= scale;
}


void calculateFDR(DATA *data, SUMMARY *summary) {
    int i,j,id;
    float avgp;
    float numer, denom, thres;
    float tmp_avgp[data->nuinter];
    
    for(i=0;i<data->nuinter;i++) {
        avgp = 0.0;
        for(j=0;j<data->n_u2a[i];j++) {
            id = data->u2a[i][j];
            avgp += summary->iZ[id] / ((float) data->n_u2a[i]);
        }
        tmp_avgp[i] = avgp;
    }
    
    for(i=0;i<data->nuinter;i++) {
        numer = 0.0;
        denom = 0.0;
        thres = tmp_avgp[i];
        for(j=0;j<data->nuinter;j++) {
            if(tmp_avgp[j] >= tmp_avgp[i]) {
                numer += (1.0 - tmp_avgp[j]);
                denom += 1.0;
            }
            
        }
        if(denom == 0.0) summary->FDR[i] = 0.0;
        else summary->FDR[i] = numer / denom;
    }
    
}



/*************************************/
/**        Histogram updates        **/
/*************************************/
void updateHist_alpha_prey(HISTOGRAM *hist, PRIOR *prior) {
  int i,j;
  for(i=0;i<_MAX_COMP_;i++) {
    if(prior->theta_alpha_prey[i] < hist->start[0]) {
	  hist->count[0] += prior->gamma_alpha_prey[i];
    }
    else if(prior->theta_alpha_prey[i] >= hist->end[_HISTO_BIN_-1]) {
      hist->count[_HISTO_BIN_ + 1] += prior->gamma_alpha_prey[i];
    }
    else {
      for(j=0;j<_HISTO_BIN_;j++) {
        if(prior->theta_alpha_prey[i] >= hist->start[j] && prior->theta_alpha_prey[i] < hist->end[j]) {
          hist->count[j+1] += prior->gamma_alpha_prey[i];
          break;
        }
      }
    }
  }
}

void updateHist_alpha_IP(HISTOGRAM *hist, PRIOR *prior) {
  int i,j;
  for(i=0;i<_MAX_COMP_;i++) {
    if(prior->theta_alpha_IP[i] < hist->start[0]) {
      hist->count[0] += prior->gamma_alpha_IP[i];
    }
    else if(prior->theta_alpha_IP[i] >= hist->end[_HISTO_BIN_-1]) {
      hist->count[_HISTO_BIN_ + 1] += prior->gamma_alpha_IP[i];
    }
    else {
      for(j=0;j<_HISTO_BIN_;j++) {
        if(prior->theta_alpha_IP[i] >= hist->start[j] && prior->theta_alpha_IP[i] < hist->end[j]) {
          hist->count[j+1] += prior->gamma_alpha_IP[i];
          break;
        }
      }
    }
  }
}

void updateHist_mu(HISTOGRAM *hist, PRIOR *prior) {
  int i,j;
  for(i=0;i<_MAX_COMP_;i++) {
    if(prior->theta_mu[i] < hist->start[0]) {
	  hist->count[0] += prior->gamma_mu[i];
    }
    else if(prior->theta_mu[i] >= hist->end[_HISTO_BIN_-1]) {
      hist->count[_HISTO_BIN_ + 1] += prior->gamma_mu[i];
    }
    else {
      for(j=0;j<_HISTO_BIN_;j++) {
        if(prior->theta_mu[i] >= hist->start[j] && prior->theta_mu[i] < hist->end[j]) {
          hist->count[j+1] += prior->gamma_mu[i];
          break;
        }
      }
    }
  }
}

void updateHist_eta(HISTOGRAM2 *hist, PRIOR *prior) {
  int i,j;
  for(i=0;i<_MAX_COMP_;i++) {
    if(prior->theta_eta[i] < hist->start[0]) {
	  hist->count[0] += prior->gamma_eta[i];
    }
    else if(prior->theta_eta[i] >= hist->end[_HISTO_BIN2_-1]) {
      hist->count[_HISTO_BIN2_ + 1] += prior->gamma_eta[i];
    }
    else {
      for(j=0;j<_HISTO_BIN2_;j++) {
        if(prior->theta_eta[i] >= hist->start[j] && prior->theta_eta[i] < hist->end[j]) {
          hist->count[j+1] += prior->gamma_eta[i];
          break;
        }
      }
    }
  }
}

void updateHist_eta0(HISTOGRAM2 *hist, PRIOR *prior) {
  int i,j;
  for(i=0;i<_MAX_COMP_;i++) {
    if(prior->theta_eta0[i] < hist->start[0]) {
	  hist->count[0] += prior->gamma_eta0[i];
    }
    else if(prior->theta_eta0[i] >= hist->end[_HISTO_BIN2_-1]) {
      hist->count[_HISTO_BIN2_ + 1] += prior->gamma_eta0[i];
    }
    else {
      for(j=0;j<_HISTO_BIN2_;j++) {
        if(prior->theta_eta0[i] >= hist->start[j] && prior->theta_eta0[i] < hist->end[j]) {
          hist->count[j+1] += prior->gamma_eta0[i];
          break;
        }
      }
    }
  }
}


void updateHistogram(PARAM *param, PRIOR *prior, DATA *data, SUMMARY *summary) {
  updateHist_alpha_prey(&(summary->hist_alpha_prey), prior);  
  updateHist_alpha_IP(&(summary->hist_alpha_IP), prior);  
  updateHist_mu(&(summary->hist_mu), prior);  
  updateHist_eta(&(summary->hist_eta), prior);  
  updateHist_eta0(&(summary->hist_eta0), prior);  
}



