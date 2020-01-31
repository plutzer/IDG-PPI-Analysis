#include "saint.h"

double vec_sum(const double *vec, int len) {
  int i;
  double res;
  res=vec[0];
  for(i=1;i<len;i++) res+=vec[i];
  return res;
}

double vec_max(const double *vec, int len) {
  int i;
  double res;
  res=vec[0];
  for(i=1;i<len;i++) {
    if(res<vec[i]) res=vec[i];
  }
  return res;
}

double vec_min(const double *vec, int len) {
  int i;
  double res;
  res=vec[0];
  for(i=1;i<len;i++) {
    if(res>vec[i]) res=vec[i];
  }
  return res;
}

double vec_mean(const double *vec, int len) {
  double tmp=0.0;
  int i;
  for(i=0;i<len;i++) tmp+=vec[i];
  tmp=tmp/((double) len);
  return tmp;
}

double vec_var(const double *vec, int len) {
  double mean=0.0;
  double var=0.0;
  int i;
  for(i=0;i<len;i++) mean+=vec[i];
  mean=mean/((double) len);
  for(i=0;i<len;i++) var+=pow((vec[i]-mean),2);
  var/=((double) (len-1));
  var=sqrt(var);
  return var;
}

double vec_med(const double *vec, int len)
{
  double new_vec[len];
  double med;
  int i, pk;
  for(i=0;i<len;i++) {
    new_vec[i]=vec[i]; 
  }
  if(len==1) {
    med=vec[0];
    return med;
  }
  else if(len%2==0) {
    gsl_sort(new_vec,1,len);
    pk=(len-2)/2;
    med=(new_vec[pk]+new_vec[pk+1])/2.0;
    return med;
  }
  else {
    gsl_sort(new_vec,1,len);
    pk=(len-1)/2;
    med=new_vec[pk];
    return med;
  }
}

double vec_mad(const double *vec, int len)
{
  double new_vec[len];
  double med, mad;
  int i;

  med = vec_med(vec,len);
  for(i=0;i<len;i++) {
    new_vec[i]=fabs(vec[i] - med); 
  }
  mad = vec_med(new_vec,len);
  return mad;
}

double poisson_unscaled_pdf(const double c, const double lambda) {
  double res;
  res = c * log(lambda) - lambda;
  res = exp(res);
  return res;
}





