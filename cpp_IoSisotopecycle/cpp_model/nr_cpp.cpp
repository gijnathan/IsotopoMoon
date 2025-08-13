#include <iostream>
#include <chrono>
#include <cmath>
#include <limits>

using namespace std;

/*
 * newton_raphson():
 *   Generic Newton–Raphson solver for a scalar function f(x) = 0.
 *
 * Parameters:
 *   x0   – initial guess for the root
 *   constants – extra parameters passed to eqs() and deriv() on each call
 *   tol  – stopping tolerance: stop when |f(x)| < tol
 *   step – step scaling factor (normally 1.0)
 *   eqs  – callable: double eqs(double x, const Constants& c)
 *           → returns the value of f(x)
 *   deriv – callable: double deriv(double x, const Constants& c)
 *           → returns the value of f'(x)
 *   maxiter – maximum allowed iterations before giving up
 *   per_iter_out – optional output stream for logging each iteration
 *       nullptr = quiet (no per-iteration output)
 *       &cout   = print per-iteration data to the terminal
 *   print_header – if true, print CSV header before per-iteration output
 *
 * Usage example for eqs/deriv:
 *   struct Params { double a; };
 *   Params constants{2.0};
 *   auto eqs   = [](double x, const Params& c) noexcept { return x*x - c.a; };
 *   auto deriv = [](double x, const Params&)   noexcept { return 2.0*x;     };
 */
template <typename Eq, typename Deriv, typename Constants>
double newton_raphson(double x0,
                      const Constants& constants,
                      double tol,
                      double step,
                      Eq eqs,
                      Deriv deriv,
                      int maxiter = 1000,
                      ostream* per_iter_out = nullptr,
                      bool print_header = true)
{
    auto dx = [&](double x) noexcept {
        const double f = eqs(x, constants);
        return (f >= 0.0) ? f : -f;
    };

    if (per_iter_out && print_header) {
        per_iter_out->setf(ios::fixed);
        per_iter_out->precision(12);
        *per_iter_out << "guessx,diff,step,f,df,f/df\n";
    }

    double delta = dx(x0);
    int n = 0;

    // Optional first row like the Python version
    if (per_iter_out) {
        *per_iter_out << x0 << ',' << delta << ',' << step << ",,,\n";
    }

    while (delta > tol && n < maxiter) {
        const double f  = eqs(x0, constants);
        const double df = deriv(x0, constants);
        if (!isfinite(df) || fabs(df) < numeric_limits<double>::min()) break;

        const double ratio = f / df;
        x0 -= step * ratio;

        delta = dx(x0);
        ++n;

        if (per_iter_out) {
            *per_iter_out << x0 << ',' << delta << ',' << step << ','
                          << f  << ',' << df << ',' << ratio << '\n';
        }
    }
    return x0;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);

    // Toggle this to enable/disable per-iteration CSV output
    const bool VERBOSE = false;

    // Example: solve x^2 - a = 0, where a = 2.0  (root = sqrt(2))
    struct Params { double a; } constants{2.0};

    // eqs(x, constants) → f(x)
    auto eqs   = [](double x, const Params& c) noexcept { return x*x - c.a; };

    // deriv(x, constants) → f'(x)
    auto deriv = [](double x, const Params&)   noexcept { return 2.0*x;     };

    const double x0   = 1.0;      // initial guess
    const double tol  = 1e-12;    // tolerance
    const double step = 1.0;      // NR step scaling

    using clock = chrono::steady_clock;
    auto t0 = clock::now();
    double root = newton_raphson(x0, constants, tol, step, eqs, deriv,
                                 1000, VERBOSE ? &cout : nullptr, true);
    auto t1 = clock::now();

    cout.setf(ios::fixed);
    cout.precision(12);
    cout << "root = " << root << '\n';
    cout << "elapsed_seconds = " 
         << chrono::duration<double>(t1 - t0).count() << '\n';

    return 0;
}
