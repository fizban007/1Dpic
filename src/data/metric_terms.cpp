#include "data/metric_terms.h"

namespace Aperture {

const double eps = 1.0e-10;

MetricTerms::MetricTerms(double a, double rg, double theta) :
    a2_(a*a), rg_(rg), theta_(theta) {
  sin2_ = std::sin(theta) * std::sin(theta);
  cos2_ = std::cos(theta) * std::cos(theta);
}

double
MetricTerms::alpha(double r) {
  if (std::abs(r_ - r) < eps) {
    return alpha_;
  } else {
    r_ = r;
    r2a2_ = r*r + a2_;
    Sigma_ = (r2a2_) * (r2a2_) - a2_ * (r2a2_ - rg_*r) * sin2_;
    rho2_ = r*r + a2_ * cos2_;
    alpha_ = std::sqrt(rho2_ * (r2a2_ - rg_*r) / Sigma_);
    return alpha_;
  }
}

double
MetricTerms::gammarr(double r) {
  if (std::abs(r_ - r) < eps) {
    return gammarr_;
  } else {
    r_ = r;
    r2a2_ = r*r + a2_;
    rho2_ = r*r + a2_ * cos2_;
    gammarr_ = (r2a2_ - rg_*r) / rho2_;
    return gammarr_;
  }
}

double
MetricTerms::gamma_p(double r, double ur) {
  alpha(r);
  gammarr(r);
  return std::sqrt(1.0 + gammarr_ * ur * ur) / alpha_;
}

double
MetricTerms::dr_alpha(double r) {
  // TODO: Finish this thing
  alpha(r);
  double r3 = r*r*r;
  double result = -rg_ * (a2_*a2_*a2_ - 2.0*r3*r3 + a2_*r3*(-3.0*r + 2.0*rg_) + a2_*(r2a2_*r2a2_ - 2.0*r3*rg_)*(cos2_ - sin2_)) / (4.0 * alpha_ * Sigma_ * Sigma_);
  return result;
}

double
MetricTerms::dr_gammarr(double r) {
  gammarr(r);
  return (2.0 * r - rg_ - 2.0 * r * gammarr_) / rho2_;
}

double
MetricTerms::sqrt_gamma(double r) {
  alpha(r);
  return Sigma_ * sqrt(sin2_) / alpha_;
}

}
