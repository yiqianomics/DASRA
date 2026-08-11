#include <Rcpp.h>

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <limits>
#include <type_traits>
#include <vector>

// [[Rcpp::plugins(cpp11)]]

// The R reference evaluates vectorized arithmetic in separate, rounded steps.
// Prevent the compiler from contracting those steps into fused operations.
#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize ("fp-contract=off")
#endif

namespace {

constexpr int kNodes = 15;

// R's colSums() uses extended accumulation when long double has additional
// precision. Match that platform-dependent behavior and round to double before
// the subsequent panel arithmetic.
using Accumulator = typename std::conditional<
  (LDBL_MANT_DIG > DBL_MANT_DIG), long double, double
>::type;

constexpr std::array<double, kNodes> nodes = {{
  -0.9914553711208126, -0.9491079123427585, -0.8648644233597691,
  -0.7415311855993945, -0.5860872354676911, -0.4058451513773972,
  -0.2077849550078985, 0.0,
   0.2077849550078985,  0.4058451513773972,  0.5860872354676911,
   0.7415311855993945,  0.8648644233597691,  0.9491079123427585,
   0.9914553711208126
}};

constexpr std::array<double, kNodes> kronrod_weights = {{
  0.02293532201052922, 0.06309209262997855, 0.1047900103222502,
  0.1406532597155259, 0.1690047266392679, 0.1903505780647854,
  0.2044329400752989, 0.2094821410847278,
  0.2044329400752989, 0.1903505780647854, 0.1690047266392679,
  0.1406532597155259, 0.1047900103222502, 0.06309209262997855,
  0.02293532201052922
}};

constexpr std::array<double, 7> gauss_weights = {{
  0.1294849661688697, 0.2797053914892767, 0.3818300505051189,
  0.4179591836734694, 0.3818300505051189, 0.2797053914892767,
  0.1294849661688697
}};

struct Panel {
  double a;
  double b;
  std::array<double, 3> value;
  std::array<double, 3> error;
};

inline double softplus(const double x) {
  return x > 0.0 ? x + std::log1p(std::exp(-x)) : std::log1p(std::exp(x));
}

inline double expm1_minus_x(const double x) {
  if (std::isfinite(x) && std::abs(x) <= 0.5) {
    if (x == 0.0) return 0.0;
    double term = x * x / 2.0;
    double total = term;
    for (int order = 3; order <= 100; ++order) {
      term = term * x / static_cast<double>(order);
      const double updated = total + term;
      if (std::abs(term) <= std::numeric_limits<double>::epsilon() *
          std::max(std::abs(updated), std::numeric_limits<double>::min())) {
        total = updated;
        break;
      }
      total = updated;
    }
    return total;
  }
  return std::expm1(x) - x;
}

inline double log1p_minus_x(const double x) {
  if (!std::isfinite(x) || x <= -1.0) {
    Rcpp::stop("log1p remainder received an input outside its domain");
  }
  if (std::abs(x) <= 0.5) {
    if (x == 0.0) return 0.0;
    double term = -x * x / 2.0;
    double total = term;
    for (int order = 3; order <= 100; ++order) {
      term = term * (-x) * static_cast<double>(order - 1) /
        static_cast<double>(order);
      const double updated = total + term;
      if (std::abs(term) <= std::numeric_limits<double>::epsilon() *
          std::max(std::abs(updated), std::numeric_limits<double>::min())) {
        total = updated;
        break;
      }
      total = updated;
    }
    return total;
  }
  return std::log1p(x) - x;
}

inline double logistic(const double x) {
  return R::plogis(x, 0.0, 1.0, true, false);
}

inline double softplus_bregman(const double h_mode,
                               const double displacement,
                               const double p_mode) {
  double base;
  double stable_displacement;
  double derivative;
  if (h_mode > 0.0) {
    base = -h_mode;
    stable_displacement = -displacement;
    derivative = logistic(base);
  } else {
    base = h_mode;
    stable_displacement = displacement;
    derivative = logistic(base);
  }

  const double shifted_softplus = softplus(base + stable_displacement);
  const double base_softplus = softplus(base);
  const double linear_part = derivative * stable_displacement;
  const double softplus_difference = shifted_softplus - base_softplus;
  double out = softplus_difference - linear_part;
  if (std::abs(stable_displacement) <= 0.5) {
    const double exponential_remainder = expm1_minus_x(stable_displacement);
    const double x = derivative * std::expm1(stable_displacement);
    const double log_remainder = log1p_minus_x(x);
    const double weighted_exponential_remainder =
      derivative * exponential_remainder;
    out = log_remainder + weighted_exponential_remainder;
  }

  const double rounding_limit = 100.0 * std::numeric_limits<double>::epsilon() *
    std::max({1.0, std::abs(softplus(base + stable_displacement)),
              std::abs(softplus(base)),
              std::abs(derivative * stable_displacement)});
  if (!std::isfinite(out) || out < -rounding_limit) {
    Rcpp::stop("softplus Bregman remainder violated convexity");
  }
  return std::max(out, 0.0);
}

inline Panel make_panel(const double a,
                        const double b,
                        const double y,
                        const double N,
                        const double eta,
                        const double sigma,
                        const double mode,
                        const double p_mode,
                        const double mode_residual) {
  const double midpoint = (a + b) / 2.0;
  const double halfwidth = (b - a) / 2.0;
  const double h_mode = eta + sigma * mode;
  std::array<Accumulator, 3> kronrod = {{0.0, 0.0, 0.0}};
  std::array<Accumulator, 3> gauss = {{0.0, 0.0, 0.0}};
  int gauss_index = 0;

  for (int j = 0; j < kNodes; ++j) {
    const double node_offset = halfwidth * nodes[j];
    const double u = midpoint + node_offset;
    const double displacement = sigma * u;
    const double divergence = softplus_bregman(h_mode, displacement, p_mode);
    const double linear_term = mode_residual * u;
    const double u_squared = u * u;
    const double quadratic_term = u_squared / 2.0;
    const double upper_bound = linear_term - quadratic_term;
    const double scaled_divergence = N * divergence;
    const double log_ratio = upper_bound - scaled_divergence;
    if (!std::isfinite(log_ratio) ||
        log_ratio > upper_bound + 1e-10 * std::max(1.0, std::abs(upper_bound))) {
      Rcpp::stop("anchored posterior ratio violated strong log concavity");
    }
    const double base = std::exp(log_ratio);
    const std::array<double, 3> value = {{base, base * u, base * u_squared}};
    for (int k = 0; k < 3; ++k) {
      const double product = value[k] * kronrod_weights[j];
      kronrod[k] += static_cast<Accumulator>(product);
    }
    // R selects rows 2,4,...,14 (one-based), i.e. odd zero-based indices.
    if ((j % 2) == 1) {
      for (int k = 0; k < 3; ++k) {
        const double product = value[k] * gauss_weights[gauss_index];
        gauss[k] += static_cast<Accumulator>(product);
      }
      ++gauss_index;
    }
  }

  Panel panel;
  panel.a = a;
  panel.b = b;
  for (int k = 0; k < 3; ++k) {
    const double kronrod_sum = static_cast<double>(kronrod[k]);
    const double gauss_sum = static_cast<double>(gauss[k]);
    panel.value[k] = halfwidth * kronrod_sum;
    panel.error[k] = std::abs(panel.value[k] - halfwidth * gauss_sum);
  }
  return panel;
}

inline Rcpp::NumericVector array_to_vector(const std::array<double, 3>& x) {
  return Rcpp::NumericVector::create(x[0], x[1], x[2]);
}

} // namespace

// One-panel export is useful for validating the compiled numerical primitive.
// [[Rcpp::export(rng = false)]]
Rcpp::List doram_gk15_panel_cpp(const double a,
                                const double b,
                                const double y,
                                const double N,
                                const double eta,
                                const double sigma,
                                const double mode,
                                const double p_mode,
                                const double mode_residual) {
  const Panel panel = make_panel(
    a, b, y, N, eta, sigma, mode, p_mode, mode_residual
  );
  return Rcpp::List::create(
    Rcpp::_["a"] = panel.a,
    Rcpp::_["b"] = panel.b,
    Rcpp::_["value"] = array_to_vector(panel.value),
    Rcpp::_["error"] = array_to_vector(panel.error)
  );
}

// This kernel deliberately starts after posterior-mode computation and ends
// before likelihood anchoring. It only replaces the adaptive vector quadrature
// loop; model algebra, optimizer constraints, and inferential calculations stay
// in R.
// [[Rcpp::export(rng = false)]]
Rcpp::List doram_integrate_moments_cpp(const double y,
                                       const double N,
                                       const double eta,
                                       const double sigma,
                                       const double mode,
                                       const double p_mode,
                                       const double mode_residual,
                                       const double local_scale,
                                       const double relative_tolerance,
                                       const double radius,
                                       const int maximum_panels,
                                       Rcpp::NumericVector tail_error) {
  if (tail_error.size() != 3 || maximum_panels < 1 || local_scale <= 0.0 ||
      relative_tolerance <= 0.0 || radius <= 0.0) {
    Rcpp::stop("invalid adaptive quadrature controls");
  }

  std::vector<double> positive;
  positive.push_back(local_scale);
  while (positive.back() < radius / 2.0) {
    positive.push_back(2.0 * positive.back());
  }
  std::vector<double> kept;
  for (double x : positive) {
    const double clipped = std::min(x, radius);
    if (clipped > 0.0 && clipped < radius &&
        (kept.empty() || clipped != kept.back())) {
      kept.push_back(clipped);
    }
  }

  std::vector<double> breaks;
  breaks.reserve(2 * kept.size() + 3);
  breaks.push_back(-radius);
  for (auto it = kept.rbegin(); it != kept.rend(); ++it) {
    breaks.push_back(-*it);
  }
  breaks.push_back(0.0);
  for (double x : kept) breaks.push_back(x);
  breaks.push_back(radius);

  std::vector<Panel> panels;
  panels.reserve(std::min(maximum_panels, 256));
  for (std::size_t i = 0; i + 1 < breaks.size(); ++i) {
    panels.push_back(make_panel(
      breaks[i], breaks[i + 1], y, N, eta, sigma, mode, p_mode,
      mode_residual
    ));
  }
  const int initial_panels = static_cast<int>(panels.size());

  std::array<double, 3> value = {{0.0, 0.0, 0.0}};
  std::array<double, 3> error = {{0.0, 0.0, 0.0}};
  std::array<double, 3> transformed = {{R_PosInf, R_PosInf, R_PosInf}};
  bool converged = false;

  while (true) {
    value = {{0.0, 0.0, 0.0}};
    error = {{0.0, 0.0, 0.0}};
    for (const Panel& panel : panels) {
      for (int k = 0; k < 3; ++k) {
        value[k] += panel.value[k];
        error[k] += panel.error[k];
      }
    }
    // Match Reduce(`+`, panel_errors) + tail_error in the R implementation.
    for (int k = 0; k < 3; ++k) error[k] += tail_error[k];

    const double denominator = value[0];
    if (std::isfinite(denominator) && denominator > error[0]) {
      const double mean_u_now = value[1] / denominator;
      const double second_u_now = value[2] / denominator;
      const double log_error_bound = -std::log1p(-error[0] / denominator);
      const double mean_u_error_bound =
        (error[1] + std::abs(mean_u_now) * error[0]) /
        (denominator - error[0]);
      const double second_u_error_bound =
        (error[2] + std::abs(second_u_now) * error[0]) /
        (denominator - error[0]);
      transformed = {{
        log_error_bound,
        mean_u_error_bound / sigma,
        2.0 * std::abs(mode) * mean_u_error_bound + second_u_error_bound
      }};
      if (std::isfinite(transformed[0]) && std::isfinite(transformed[1]) &&
          std::isfinite(transformed[2]) &&
          transformed[0] <= relative_tolerance &&
          transformed[1] <= relative_tolerance &&
          transformed[2] <= relative_tolerance) {
        converged = true;
        break;
      }
    } else {
      transformed = {{R_PosInf, R_PosInf, R_PosInf}};
    }

    if (static_cast<int>(panels.size()) >= maximum_panels) break;

    std::size_t split_index = 0;
    double best_ratio = -R_PosInf;
    for (std::size_t i = 0; i < panels.size(); ++i) {
      double ratio;
      if (!std::isfinite(denominator) || denominator <= 0.0) {
        ratio = std::max({panels[i].error[0], panels[i].error[1],
                          panels[i].error[2]});
      } else {
        const double local_mean =
          (panels[i].error[1] + std::abs(value[1] / denominator) *
           panels[i].error[0]) / denominator;
        const double local_second =
          (panels[i].error[2] + std::abs(value[2] / denominator) *
           panels[i].error[0]) / denominator;
        ratio = std::max({
          panels[i].error[0] / denominator,
          local_mean / sigma,
          2.0 * std::abs(mode) * local_mean + local_second
        });
      }
      // which.max() selects the first maximum.
      if (ratio > best_ratio) {
        best_ratio = ratio;
        split_index = i;
      }
    }

    const Panel old = panels[split_index];
    const double midpoint = (old.a + old.b) / 2.0;
    const Panel left = make_panel(
      old.a, midpoint, y, N, eta, sigma, mode, p_mode, mode_residual
    );
    const Panel right = make_panel(
      midpoint, old.b, y, N, eta, sigma, mode, p_mode, mode_residual
    );
    panels.erase(panels.begin() + split_index);
    // Match append(panels[-split_index], children) in the R implementation.
    panels.push_back(left);
    panels.push_back(right);
  }

  if (!converged) {
    Rcpp::stop("vector Gauss--Kronrod integration did not meet its error contract");
  }

  return Rcpp::List::create(
    Rcpp::_["value"] = array_to_vector(value),
    Rcpp::_["error"] = array_to_vector(error),
    Rcpp::_["transformed_error"] = array_to_vector(transformed),
    Rcpp::_["panels"] = static_cast<int>(panels.size()),
    Rcpp::_["evaluations"] = 15 * initial_panels +
      30 * (static_cast<int>(panels.size()) - initial_panels)
  );
}
