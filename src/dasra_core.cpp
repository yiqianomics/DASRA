#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

inline double clamp_scalar(const double x, const double lo, const double hi) {
    return std::min(std::max(x, lo), hi);
}

inline double softplus_scalar(const double x) {
    return std::max(x, 0.0) + std::log1p(std::exp(-std::abs(x)));
}

inline bool integerish(const double x) {
    return std::abs(x - std::round(x)) <= 1e-8;
}

struct ModeResult {
    double mode;
    double curvature;
    double score;
    double bracket_width;
};

ModeResult latent_mode_one(const double y,
                           const double N,
                           const double eta,
                           const double sigma,
                           const int iterations,
                           const double score_tolerance,
                           const double bracket_tolerance) {
    const double sigma2 = sigma * sigma;
    const double prior_precision = 1.0 / sigma2;
    const double pseudo_p = clamp_scalar(
        (y + 0.5) / (N + 1.0), 1e-12, 1.0 - 1e-12
    );
    double x = (
        R::qlogis(pseudo_p, 0.0, 1.0, true, false) +
        prior_precision * eta
    ) / (1.0 + prior_precision);
    double lower = eta - sigma2 * (N - y);
    double upper = eta + sigma2 * y;
    x = std::min(std::max(x, lower), upper);

    for (int iter = 0; iter < iterations; ++iter) {
        const double p = R::plogis(x, 0.0, 1.0, true, false);
        const double score = y - N * p - (x - eta) * prior_precision;
        const double curvature = N * p * (1.0 - p) + prior_precision;
        if (score > 0.0) {
            lower = x;
        } else {
            upper = x;
        }
        const double newton = x + score / curvature;
        const double midpoint = lower + 0.5 * (upper - lower);
        const double bracket_width = upper - lower;
        const bool use_newton = R_FINITE(newton) && newton > lower &&
            newton < upper &&
            std::abs(newton - x) <= 0.5 * bracket_width;
        double next_x = use_newton ? newton : midpoint;
        const double newton_step_size = std::abs(score / curvature);
        const bool converged =
            newton_step_size <= score_tolerance * (1.0 + std::abs(x)) ||
            std::abs(upper - lower) <=
                bracket_tolerance * (1.0 + std::abs(x));
        if (converged) {
            break;
        }
        x = next_x;
    }

    const double p = R::plogis(x, 0.0, 1.0, true, false);
    const double score = y - N * p - (x - eta) * prior_precision;
    const double curvature = N * p * (1.0 - p) + prior_precision;
    const double newton_step_size = std::abs(score / curvature);
    if (!R_FINITE(x) || !R_FINITE(curvature) || curvature <= 0.0 ||
        newton_step_size > 1e-8 * (1.0 + std::abs(x))) {
        Rcpp::stop("Safeguarded latent-mode solver failed to converge.");
    }
    ModeResult result;
    result.mode = x;
    result.curvature = curvature;
    result.score = score;
    result.bracket_width = upper - lower;
    return result;
}

void validate_inputs(const Rcpp::NumericVector& y,
                     const Rcpp::NumericVector& N,
                     const Rcpp::NumericVector& eta,
                     const double sigma) {
    const R_xlen_t n = y.size();
    if (N.size() != n || eta.size() != n) {
        Rcpp::stop("y, N, and eta must have the same length.");
    }
    if (!R_FINITE(sigma) || sigma <= 0.0) {
        Rcpp::stop("sigma must be finite and positive.");
    }
    for (R_xlen_t i = 0; i < n; ++i) {
        const double yi = y[i];
        const double Ni = N[i];
        const double etai = eta[i];
        if (!R_FINITE(yi) || !R_FINITE(Ni) || !R_FINITE(etai) ||
            yi < 0.0 || Ni < 0.0 || yi > Ni ||
            !integerish(yi) || !integerish(Ni)) {
            Rcpp::stop(
                "y and N must be finite integer counts satisfying "
                "0 <= y <= N, and eta must be finite."
            );
        }
    }
}

void validate_solver_controls(const int iterations,
                              const double score_tolerance,
                              const double bracket_tolerance) {
    if (iterations < 1 || !R_FINITE(score_tolerance) ||
        !R_FINITE(bracket_tolerance) || score_tolerance <= 0.0 ||
        bracket_tolerance <= 0.0) {
        Rcpp::stop("Invalid latent-mode solver controls.");
    }
}

}  // namespace

// [[Rcpp::export]]
Rcpp::NumericVector dasra_gh_log_weights_cpp(
        const Rcpp::NumericVector& node,
        const int order) {
    if (order < 3 || node.size() != order) {
        Rcpp::stop("The quadrature order must match at least three nodes.");
    }

    const long double pi =
        3.141592653589793238462643383279502884L;
    const long double p0 = std::pow(pi, -0.25L);
    const long double sqrt_two = std::sqrt(2.0L);
    Rcpp::NumericVector log_weight(order);

    // At a root x_q of the physicists' Hermite polynomial H_order,
    // w_q = 1 / {order * p_{order-1}(x_q)^2}, where p_k is orthonormal
    // for exp(-x^2).  The recurrence is rescaled at every step so the
    // logarithm remains representable even when w_q itself underflows.
    for (int q = 0; q < order; ++q) {
        const long double x = static_cast<long double>(node[q]);
        if (!std::isfinite(static_cast<double>(x))) {
            Rcpp::stop("Quadrature nodes must be finite.");
        }

        long double previous = p0;
        long double current = sqrt_two * x * p0;
        long double log_scale = 0.0L;
        long double scale = std::max(std::abs(previous), std::abs(current));
        if (!(scale > 0.0L) || !std::isfinite(static_cast<double>(scale))) {
            Rcpp::stop("Hermite recurrence initialization failed.");
        }
        previous /= scale;
        current /= scale;
        log_scale += std::log(scale);

        for (int k = 1; k <= order - 2; ++k) {
            const long double next =
                std::sqrt(2.0L / static_cast<long double>(k + 1)) *
                    x * current -
                std::sqrt(static_cast<long double>(k) /
                          static_cast<long double>(k + 1)) * previous;
            scale = std::max(std::abs(current), std::abs(next));
            if (!(scale > 0.0L) ||
                !std::isfinite(static_cast<double>(scale))) {
                Rcpp::stop("Hermite recurrence failed.");
            }
            previous = current / scale;
            current = next / scale;
            log_scale += std::log(scale);
        }

        const long double log_polynomial =
            std::log(std::abs(current)) + log_scale;
        const long double value =
            -std::log(static_cast<long double>(order)) -
            2.0L * log_polynomial;
        if (!std::isfinite(static_cast<double>(value))) {
            Rcpp::stop("A Gauss-Hermite log weight is non-finite.");
        }
        log_weight[q] = static_cast<double>(value);
    }

    // Remove the negligible common roundoff drift while retaining the
    // standard Hermite normalization sum_q w_q = sqrt(pi).
    double maximum = R_NegInf;
    for (int q = 0; q < order; ++q) {
        maximum = std::max(maximum, log_weight[q]);
    }
    long double scaled_sum = 0.0L;
    for (int q = 0; q < order; ++q) {
        scaled_sum += std::exp(
            static_cast<long double>(log_weight[q] - maximum)
        );
    }
    const double log_sum = maximum +
        std::log(static_cast<double>(scaled_sum));
    const double correction = 0.5 * std::log(static_cast<double>(pi)) -
        log_sum;
    for (int q = 0; q < order; ++q) {
        log_weight[q] += correction;
    }
    return log_weight;
}

// [[Rcpp::export]]
Rcpp::NumericVector dasra_count_log_hy_adaptive_cpp(
        const Rcpp::NumericVector& y,
        const Rcpp::NumericVector& N,
        const Rcpp::NumericVector& eta,
        const double sigma,
        const Rcpp::NumericVector& node,
        const Rcpp::NumericVector& log_raw_weight,
        const int iterations = 80,
        const double score_tolerance = 1e-12,
        const double bracket_tolerance = 1e-12) {
    validate_inputs(y, N, eta, sigma);
    validate_solver_controls(
        iterations, score_tolerance, bracket_tolerance
    );
    const R_xlen_t n = y.size();
    const R_xlen_t Q = node.size();
    if (Q < 1 || log_raw_weight.size() != Q) {
        Rcpp::stop(
            "node and log_raw_weight must have the same positive length."
        );
    }
    for (R_xlen_t q = 0; q < Q; ++q) {
        if (!R_FINITE(node[q]) ||
            Rcpp::NumericVector::is_na(log_raw_weight[q]) ||
            log_raw_weight[q] == R_PosInf) {
            Rcpp::stop("Quadrature nodes and log weights are invalid.");
        }
    }

    const double sigma2 = sigma * sigma;
    const double pi = 3.141592653589793238462643383279502884;
    const double log_normalizer = std::log(sigma) + 0.5 * std::log(2.0 * pi);
    Rcpp::NumericVector answer(n);

    for (R_xlen_t i = 0; i < n; ++i) {
        const ModeResult mode_result = latent_mode_one(
            y[i], N[i], eta[i], sigma, iterations,
            score_tolerance, bracket_tolerance
        );
        const double mode = mode_result.mode;
        const double curvature = mode_result.curvature;
        const double scale = std::sqrt(2.0 / curvature);

        double maximum = -std::numeric_limits<double>::infinity();
        long double scaled_sum = 0.0L;
        for (R_xlen_t q = 0; q < Q; ++q) {
            const double x = mode + scale * node[q];
            const double centered = x - eta[i];
            const double log_kernel =
                y[i] * x - N[i] * softplus_scalar(x) -
                centered * centered / (2.0 * sigma2) - log_normalizer;
            const double term = log_kernel + node[q] * node[q] +
                log_raw_weight[q];
            if (term > maximum) {
                if (R_FINITE(maximum)) {
                    scaled_sum = scaled_sum * std::exp(maximum - term) +
                        1.0L;
                } else {
                    scaled_sum = 1.0L;
                }
                maximum = term;
            } else if (R_FINITE(term)) {
                scaled_sum += std::exp(term - maximum);
            }
        }
        const double log_sum_exp = R_FINITE(maximum) && scaled_sum > 0.0L ?
            maximum + std::log(static_cast<double>(scaled_sum)) : maximum;

        const double log_choose =
            R::lgammafn(N[i] + 1.0) - R::lgammafn(y[i] + 1.0) -
            R::lgammafn(N[i] - y[i] + 1.0);
        answer[i] = log_choose + 0.5 * std::log(2.0 / curvature) +
            log_sum_exp;
    }
    return answer;
}

// [[Rcpp::export]]
Rcpp::List dasra_count_moments_adaptive_cpp(
        const Rcpp::NumericVector& y,
        const Rcpp::NumericVector& N,
        const Rcpp::NumericVector& eta,
        const double sigma,
        const Rcpp::NumericVector& node,
        const Rcpp::NumericVector& log_raw_weight,
        const bool need_moments = true,
        const int iterations = 80,
        const double score_tolerance = 1e-12,
        const double bracket_tolerance = 1e-12) {
    validate_inputs(y, N, eta, sigma);
    validate_solver_controls(
        iterations, score_tolerance, bracket_tolerance
    );
    const R_xlen_t n = y.size();
    const R_xlen_t Q = node.size();
    if (Q < 1 || log_raw_weight.size() != Q) {
        Rcpp::stop(
            "node and log_raw_weight must have the same positive length."
        );
    }
    for (R_xlen_t q = 0; q < Q; ++q) {
        if (!R_FINITE(node[q]) ||
            Rcpp::NumericVector::is_na(log_raw_weight[q]) ||
            log_raw_weight[q] == R_PosInf) {
            Rcpp::stop("Quadrature nodes and log weights are invalid.");
        }
    }

    const double sigma2 = sigma * sigma;
    const double pi = 3.141592653589793238462643383279502884;
    const double log_normalizer =
        std::log(sigma) + 0.5 * std::log(2.0 * pi);
    Rcpp::NumericVector log_h(n);
    Rcpp::NumericVector posterior_m1(n);
    Rcpp::NumericVector posterior_m2(n);
    Rcpp::NumericVector score_eta(n);
    Rcpp::NumericVector score_omega(n);

    for (R_xlen_t i = 0; i < n; ++i) {
        const ModeResult mode_result = latent_mode_one(
            y[i], N[i], eta[i], sigma, iterations,
            score_tolerance, bracket_tolerance
        );
        const double mode = mode_result.mode;
        const double curvature = mode_result.curvature;
        const double scale = std::sqrt(2.0 / curvature);

        double maximum = -std::numeric_limits<double>::infinity();
        long double scaled_sum = 0.0L;
        long double scaled_m1 = 0.0L;
        long double scaled_m2 = 0.0L;
        for (R_xlen_t q = 0; q < Q; ++q) {
            const double x = mode + scale * node[q];
            const double centered = x - eta[i];
            const double log_kernel =
                y[i] * x - N[i] * softplus_scalar(x) -
                centered * centered / (2.0 * sigma2) - log_normalizer;
            const double term = log_kernel + node[q] * node[q] +
                log_raw_weight[q];
            if (term > maximum) {
                const long double rescale = R_FINITE(maximum) ?
                    std::exp(static_cast<long double>(maximum - term)) :
                    0.0L;
                scaled_sum = scaled_sum * rescale + 1.0L;
                scaled_m1 = scaled_m1 * rescale + centered;
                scaled_m2 = scaled_m2 * rescale + centered * centered;
                maximum = term;
            } else if (R_FINITE(term)) {
                const long double weight = std::exp(
                    static_cast<long double>(term - maximum)
                );
                scaled_sum += weight;
                scaled_m1 += weight * centered;
                scaled_m2 += weight * centered * centered;
            }
        }
        if (!R_FINITE(maximum) || !(scaled_sum > 0.0L)) {
            Rcpp::stop("Adaptive quadrature produced a non-finite marginal.");
        }
        const double log_sum_exp = maximum +
            std::log(static_cast<double>(scaled_sum));
        const double log_choose =
            R::lgammafn(N[i] + 1.0) - R::lgammafn(y[i] + 1.0) -
            R::lgammafn(N[i] - y[i] + 1.0);
        log_h[i] = log_choose + 0.5 * std::log(2.0 / curvature) +
            log_sum_exp;
        if (need_moments) {
            posterior_m1[i] = static_cast<double>(scaled_m1 / scaled_sum);
            posterior_m2[i] = static_cast<double>(scaled_m2 / scaled_sum);
            score_eta[i] = posterior_m1[i] / sigma2;
            score_omega[i] = posterior_m2[i] / sigma2 - 1.0;
        }
    }

    if (!need_moments) {
        return Rcpp::List::create(Rcpp::Named("log_h") = log_h);
    }
    return Rcpp::List::create(
        Rcpp::Named("log_h") = log_h,
        Rcpp::Named("posterior_m1") = posterior_m1,
        Rcpp::Named("posterior_m2") = posterior_m2,
        Rcpp::Named("score_eta") = score_eta,
        Rcpp::Named("score_omega") = score_omega
    );
}

// [[Rcpp::export]]
Rcpp::NumericVector dasra_mean_log_relative_cpp(
        const Rcpp::NumericVector& location,
        const double sigma,
        const Rcpp::NumericVector& z,
        const Rcpp::NumericVector& weight) {
    if (!R_FINITE(sigma) || sigma <= 0.0) {
        Rcpp::stop("sigma must be finite and positive.");
    }
    if (z.size() < 1 || weight.size() != z.size()) {
        Rcpp::stop("z and weight must have the same positive length.");
    }
    long double weight_sum = 0.0L;
    for (R_xlen_t q = 0; q < z.size(); ++q) {
        if (!R_FINITE(z[q]) || !R_FINITE(weight[q]) || weight[q] < 0.0) {
            Rcpp::stop("Quadrature nodes and weights must be finite.");
        }
        weight_sum += weight[q];
    }
    if (!(weight_sum > 0.0L)) {
        Rcpp::stop("Quadrature weights must have a positive sum.");
    }

    Rcpp::NumericVector answer(location.size());
    for (R_xlen_t i = 0; i < location.size(); ++i) {
        if (!R_FINITE(location[i])) {
            Rcpp::stop("location must contain only finite values.");
        }
        long double total = 0.0L;
        for (R_xlen_t q = 0; q < z.size(); ++q) {
            const double x = location[i] + sigma * z[q];
            const double log_relative = -softplus_scalar(-x);
            total += weight[q] * log_relative;
        }
        answer[i] = static_cast<double>(total / weight_sum);
    }
    return answer;
}

