#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;

//[[Rcpp::export]]
NumericVector first(NumericVector za, NumericVector zb, NumericVector stdv,
                    NumericVector nints,
                    double delta, LogicalVector bs) {
  double hh = (zb[0]-za[0])/nints[0];
  NumericVector last(nints[0]+1);

  for(int i = 0; i<nints[0]+1; ++i){
    last[i] = R::dnorm(za[0] + hh*i, 0, 1, 0);
    if(bs[0]){
      last[i] = R::dnorm(za[0] + hh*i, delta*stdv[0], 1, 0);
    }
  }
  return last;
}


// [[Rcpp::export]]
NumericVector init_int(NumericVector wj, NumericVector zj, double delta, NumericVector stdv) {
  int n = zj.size();
  NumericVector last(n);

  for(int i = 0; i<n; ++i){
    last[i] = wj[i]*(R::dnorm(zj[i], delta*stdv[0], 1, 0));
  }

  return last;
}

// [[Rcpp::export]]
double trap(NumericVector f, int n, double h) {
  double sum1 = f[0];

  for(int i = 1; i < n; ++i){
    sum1 += f[i]*2;
  }

  sum1 += f[n];
  return (h/2)*sum1;
}


// [[Rcpp::export]]
double fcab(NumericVector last, int nint, double zam1,
            double h, double x, NumericMatrix stdv, double delta,
            int k, LogicalVector bs) {
  NumericVector f(nint+1);
  NumericVector grid2(nint+1);

  for(int i = 0; i < nint+1; ++i){
    grid2[i] = (zam1+h*(i))*stdv(k-2,1);
    if(bs){
      f[i] = last[i]*stdv(k-1,1)/stdv(k-1,0)*R::dnorm((x-grid2[i])/stdv(k-1,0),
                          delta*stdv(k-1,0),1,0);
    } else {
      f[i] = last[i]*stdv(k-1,1)/stdv(k-1,0)*R::dnorm((grid2[i]-x)/stdv(k-1,0),
                          delta*stdv(k-1,0),1,0);
    }
  }
  return trap(f, nint, h);
}

// [[Rcpp::export]]
NumericVector other(NumericVector za, NumericVector zb, int k, NumericMatrix stdv,
                    NumericVector last, NumericVector nints,
                    double delta, LogicalVector bs) {
  NumericVector fn(nints[k-1]+1);
  double hh = (zb[k-1]-za[k-1])/nints[k-1];
  double hlast = (zb[k-2]-za[k-2])/nints[k-2];
  NumericVector grid1(500);

  for(int i = 0; i < nints[k-1]+1; ++i){
    grid1[i] = (za[k-1]+hh*i)*stdv(k-1,1);
    fn[i] = fcab(last, nints[k-2], za[k-2], hlast, grid1[i], stdv, delta,
                 k, bs);
  }
  return fn;
}

// [[Rcpp::export]]
NumericVector recur_int(int k,
                         NumericMatrix stdv, NumericVector zj,
                         NumericVector last, NumericVector zj_up,
                         NumericVector wj_up, double delta, bool bs) {
  NumericVector last_up = zj_up.size();

  for(int i = 0; i < zj_up.size(); ++i){
  for(int j = 0; j < last.size(); ++j){
      if(bs){
        last_up[i] += last[j]*stdv(k-1,1)/stdv(k-1,0)*R::dnorm((zj[j]*stdv(k-2,1)-zj_up[i]*stdv(k-1,1))/stdv(k-1,0),
                            delta*stdv(k-1,0),1,0);
      } else if(false == bs){
        last_up[i] += last[j]*stdv(k-1,1)/stdv(k-1,0)*R::dnorm((zj_up[i]*stdv(k-1,1)-zj[j]*stdv(k-2,1))/stdv(k-1,0),
                            delta*stdv(k-1,0),1,0);
      }
    }
  last_up[i] *= wj_up[i];
  }
  return last_up;
}

// [[Rcpp::export]]
double qpos(double xq, NumericVector last, int nint, int k, double zam1,
            double zbm1, NumericMatrix stdv, bool bs, double delta){
  NumericVector f1(nint+1);
  double hlast = (zbm1 - zam1)/nint;
  double grid3;

  for(int i = 0; i < nint + 1; ++i){
    grid3 = (zam1 + hlast*i)*stdv(k-2,1);
    if(false == bs && delta != 0){
      f1[i] = last[i]*(R::pnorm((grid3-xq)/stdv(k-1,0),-delta*stdv(k-1,0), 1, 1, 0));
    } else if(bs){
      f1[i] = last[i]*(R::pnorm((xq-grid3)/stdv(k-1,0),delta*stdv(k-1,0), 1, 1, 0));
    } else{
      f1[i] = last[i]*(R::pnorm((grid3-xq)/stdv(k-1,0),delta*stdv(k-1,0), 1, 1, 0));
    }
  }
  return(trap(f1,nint,hlast));
}

// [[Rcpp::export]]
double prob(double xq, NumericVector last, NumericVector zj, int k, NumericMatrix stdv,
                 bool bs, double delta){
  double p_out = 0;

  for(int i = 0; i < zj.size(); ++i){
    if(false == bs && delta != 0){
      p_out += last[i]*(R::pnorm((zj[i]*stdv(k-2,1)-xq)/stdv(k-1,0),-delta*stdv(k-1,0), 1, 1, 0));
    } else if(bs){
      p_out += last[i]*(R::pnorm((xq-zj[i]*stdv(k-2,1))/stdv(k-1,0),delta*stdv(k-1,0), 1, 1, 0));
    } else{
      p_out += last[i]*(R::pnorm((zj[i]*stdv(k-2,1)-xq)/stdv(k-1,0),delta*stdv(k-1,0), 1, 1, 0));
    }
  }
  return p_out;
}
