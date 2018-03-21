#ifndef _METRIC_TERMS_H_
#define _METRIC_TERMS_H_

struct MetricTerms
{
  MetricTerms(double a, double rg, double theta);

  double alpha(double r);
  double gammarr(double r);
  double gamma_p(double r, double ur);
  double dr_alpha(double r);
  double dr_gammarr(double r);
  double sqrt_gamma(double r);

  // parameters
  double a2_;
  double rg_;
  double theta_;

  // Memoized values
  double r_;
  double ur_;

  double cos2_;
  double sin2_;
  double r2a2_;
  double rho2_;
  double alpha_;
  double Sigma_;
  double gammarr_;
};


#endif  // _METRIC_TERMS_H_

