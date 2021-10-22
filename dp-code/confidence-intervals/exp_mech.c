/* This file contains the main steps of the Algorithm 5 EXPQ performed in c language to reduce computation time*/
#include "math.h"
void priv_median_exp(double* db, int* n, double* e, int* qi, double* probs, double* r, int* priv_qi) {
  double e_ = *e;
  int m = *qi;
  double r_ = *r;
  for(int i = 0; i < m; i += 1) {
    double utility = (i + 1) - m;
    probs[i] = (db[i + 1] - db[i])*exp((e_)*utility/2);
  }
  for(int i = m; i <= *n; i += 1) {
    double utility = m - i;
    probs[i] = (db[i + 1] - db[i])*exp((e_)*utility/2);
  }
  
  double sum = 0;
  for(int i = 0; i <= *n; i += 1) sum += probs[i];
  r_ *= sum;
  
  for(int i = 0; i <= *n; i += 1) {
    r_ -= probs[i];
    if(r_ < 0) {
      *priv_qi = i + 1;
      break;
    }
  }
}