#include "saint.h"

float vec_sum(const float *vec, int len) {
  int i;
  float res;
  res=vec[0];
  for(i=1;i<len;i++) res+=vec[i];
  return res;
}

float vec_max(const float *vec, int len) {
  int i;
  float res;
  res=vec[0];
  for(i=1;i<len;i++) {
    if(res<vec[i]) res=vec[i];
  }
  return res;
}

float vec_min(const float *vec, int len) {
  int i;
  float res;
  res=vec[0];
  for(i=1;i<len;i++) {
    if(res>vec[i]) res=vec[i];
  }
  return res;
}

float vec_mean(const float *vec, int len) {
  float tmp=0.0;
  int i;
  for(i=0;i<len;i++) tmp+=vec[i];
  tmp=tmp/((float) len);
  return tmp;
}

float vec_var(const float *vec, int len) {
  float mean=0.0;
  float var=0.0;
  int i;
  for(i=0;i<len;i++) mean+=vec[i];
  mean=mean/((float) len);
  for(i=0;i<len;i++) var+=pow((vec[i]-mean),2);
  var/=((float) (len-1));
  var=sqrt(var);
  return var;
}

float vec_med(const float *vec, int len)
{
  double new_vec[len];
  double med;
  int i, pk;
  for(i=0;i<len;i++) {
    new_vec[i]= ((double) vec[i]); 
  }
  if(len==1) {
    med= vec[0];
    return ((float) med);
  }
  else if(len%2==0) {
    gsl_sort(new_vec,1,len);
    pk=(len-2)/2;
    med=(new_vec[pk]+new_vec[pk+1])/2.0;
    return ((float) med);
  }
  else {
    gsl_sort(new_vec,1,len);
    pk=(len-1)/2;
    med=new_vec[pk];
    return ((float) med);
  }
}

float vec_mad(const float *vec, int len)
{
  float new_vec[len];
  float med, mad;
  int i;

  med = vec_med(vec,len);
  for(i=0;i<len;i++) {
    new_vec[i]=fabs(vec[i] - med); 
  }
  mad = vec_med(new_vec,len);
  return mad;
}

int ranMultinom(const gsl_rng *r, float *p, int K) {
  int i, rr;
  float coin, sum; 
  sum = vec_sum(p,K);
  /*for(i=0;i<K-1;i++) fprintf(stderr, "%.2f\t", p[i]);
  fprintf(stderr, "%.2f\n", p[K-1]); */
  if(sum != 1.0) {
    for(i=0;i<K;i++) p[i] /= sum;
  }
  coin = gsl_ran_flat(r,0.0,1.0);
  sum = p[0];
  rr = 0;
  while(coin > sum) {
    rr++;
    sum += p[rr];
  }
  if(rr >= K) rr = K-1;
  return rr;
}


float geometric_mean(float *x, int n) {
  int i;
  float res = 0.0;
  for(i=0;i<n;i++) res += log(x[i]) / ((float) n);
  res = exp(res);
  return res;
}


