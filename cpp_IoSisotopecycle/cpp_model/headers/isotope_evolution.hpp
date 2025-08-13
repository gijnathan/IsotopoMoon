#ifndef ISOTOPE_EVOLUTION_HPP
#define ISOTOPE_EVOLUTION_HPP
#include <unordered_map>
#include <cmath>
#include <map>
#include <string>
#include <type_traits>
#include <tuple>      // <-- needed by tuple returns
#include <iostream>   // <-- needed for cerr used in R2d
using namespace std;

// ###############
// ### OPTIONS ###
// ###############

// # options for calculating fractionation factors
// nofrac = "no" # if "yes", all fractionation factors = 1 (i.e., no fractionation); if "no", fractionation factors are pre-determined values
// MAF = "no" # if "yes", mass-anomalous (i.e., with powers not equal to canonical values) is possible; if "no", fractionation factors follow canonical power laws

// # options for calculating cap-delta values
// # "O2017_eq5" = equation (5) in Ono (2017): D3XS = d3XS - e*d34S 
// # "O2017_eq6" = equation (6) in Ono (2017):  D3X = d3X - ((d34+1.)^e3X - 1.)
// # "FW2003_pg3" = equations on pg. 3 of Farquar & Wing (2003): D3XS = d3XS - 1000*((1+d34S/1000)^e - 1)
// MIF_eq = "FW2003_pg3" 

// #####################
// ### SYSTEM VALUES ###
// #####################

// # reservoir mass bounds
const double M_mass_i_min = 2.95E+21;  // mantle S mol
const double M_mass_i_max = 2.24E+22; 
const double M_mass_i_min_sil = 7.26E+22;  // mantle silicate kg
const double M_mass_i_max_sil = 7.63E+22;  // kg
const double F_mass_i_min = 0.;  // frost S mol (crust)
const double F_mass_i_max = 1.92367E+20; 
const double SS_mass_i_min = 0.;  // silicates/sulfates S mol (crust)
const double SS_mass_i_max = 7.9E+19; 

// canonical mass-dependent fractionation powers
const double e33 = 0.515;
const double e36 = 1.9;

// Function declarations and definitions

// mass-dependent fractionation factors
inline double a_3X(double a34, int X) {
    double e = 0.0;
    if (X == 33) {
        e = e33;
    } else if (X == 36) {
        e = e36;
    }
    return pow(a34, e);
}

// obtain alphas (the fractionation factors) for different processes and isotopes
inline double alphas(const string& process, int X, const string& nofrac, const string& MAF) {
    double a = 1.0;

    if (nofrac == "yes") {
        return 1.0;
    }

    double a34 = 1.0;

    if (process == "gr") {
        a34 = 0.916820406;
    } else if (process == "th") {
        a34 = 0.9998;
    } else if (process == "ed") {
        a34 = 0.994;
    } else if (process == "pd") {
        a34 = 1.008725;
    } else if (process == "hg") {
        a34 = 0.998;
    } else if (process == "sq") {
        a34 = 0.9985;
    } else if (process == "cf") {
        a34 = 0.9998;
    } else if (
        process == "pu" || process == "ac" || process == "ei" || process == "rc" ||
        process == "ec" || process == "pi" || process == "mm" || process == "dg" ||
        process == "xt" || process == "fr" || process == "rm" || process == "vp" ||
        process == "dm"
    ) {
        a34 = 1.0;
    }

    if (X == 34) {
        a = a34;
    } else if (X == 33 || X == 36) {
        double a33 = 1.0, a36 = 1.0;

        if (MAF == "no") {
            a33 = a_3X(a34, 33);
            a36 = a_3X(a34, 36);
        } else if (MAF == "yes") {
            if (process == "gr") {
                a33 = 0.957445318;
                a36 = 0.840434662;
            } else if (process == "ed") {
                a33 = 0.996;
                a36 = 0.989;
            } else if (
                process == "pd" || process == "th" || process == "pu" || process == "ac" ||
                process == "ei" || process == "rc" || process == "ec" || process == "pi" ||
                process == "mm" || process == "dg" || process == "xt" || process == "fr" ||
                process == "hg" || process == "rm" || process == "vp" || process == "sq" ||
                process == "dm" || process == "cf"
            ) {
                a33 = a_3X(a34, 33);
                a36 = a_3X(a34, 36);
            }
        }

        if (X == 33) {
            a = a33;
        } else if (X == 36) {
            a = a36;
        }
    }

    return a;
}


// ----------------------
// Rate bounds structure
// ----------------------
struct RateBounds {
    double min;
    double max;
};

// ----------------------
// Surface constants
// ----------------------
const double Io_SA = 4.174e13;     // Io's surface area in m^2
const double sil_mag_rho = 3000.;  // silicate melt density in kg/m^3

// ----------------------
// Process SO2 release rates (kg/s)
// ----------------------
const unordered_map<string, RateBounds> rate_kgs = {
    {"pd", {3028.0, 3533.0}},
    {"pi", {205.0, 326.0}},
    {"ei", {105.0, 105.0}},
    {"ed", {1500.0, 1500.0}},
    {"ac", {30.0, 30.0}},
    {"rc", {1600.0, 1600.0}},
    {"ec", {110.0, 110.0}},
    {"pu", {1000.0, 1000.0}}  // only min provided, so max = min
};

// ----------------------
// VCDT reference isotope ratios
// ----------------------
const map<string, double> VCDT = {
    {"34", 1. / 22.6436},
    {"33", 1. / 126.948},
    {"36", 1. / 6515.}
};

// ratio R -> delta value (‰)
inline double R2d(double R, int n, const map<string, double>& VCDT) {
    double std = 1.0;

    if (n == 34) {
        std = VCDT.at("34");
    } else if (n == 33) {
        std = VCDT.at("33");
    } else if (n == 36) {
        std = VCDT.at("36");
    } else {
        cerr << "Invalid isotope number: " << n << endl;
        return NAN;
    }

    return ((R - std) / std) * 1000.0;
}

// delta value (‰) -> ratio R
inline double d2R(double d, int n, const map<string, double>& VCDT) {
    double std = 1.0;
    if (n == 34)      std = VCDT.at("34");
    else if (n == 33) std = VCDT.at("33");
    else if (n == 36) std = VCDT.at("36");
    else              return NAN; // invalid isotope
    return ((d / 1000.0) * std) + std;
}

// ratio-over-total: a / (a + b + c + d)
// from ratio (e.g., 33S/32S = a) to over total (e.g., 33S/ST = X), where c, b, d are the other ratios (i.e., 32S/32S, 34S/32S, 36S/32S)
inline double ROT(double a, double c, double b, double d) {
    return a / (a + b + c + d);
}


// kg/s SO2 -> mol S / s
inline double kgs2mSs(double kgs) {
    return (kgs * 1000.0) / 64.0;
}

// MIF: compute D33, D36 from ratios and chosen equation
// Cap-Delta notation
inline pair<double, double> MIF(double R33, double R34, double R36,
                                const map<string, double>& VCDT,
                                const string& MIF_eq) {
    double d33 = R2d(R33, 33, VCDT);
    double d34 = R2d(R34, 34, VCDT);
    double d36 = R2d(R36, 36, VCDT);

    auto calc = [&](double d34v, double d3X, double e3X) -> double {
        if (MIF_eq == "O2017_eq5") {
            return d3X - e3X * d34v;
        } else if (MIF_eq == "O2017_eq6") {
            return d3X - (pow(d34v + 1.0, e3X) - 1.0);
        } else { // "FW2003_pg3" default
            return d3X - 1000.0 * (pow(1.0 + d34v / 1000.0, e3X) - 1.0);
        }
    };

    double D33 = calc(d34, d33, e33);
    double D36 = calc(d34, d36, e36);
    return {D33, D36};
}

// calc_rate overloads to handle "N" or numeric f
// One function that accepts strings or numbers, Python-style
template <typename T>
inline double calc_rate(T&& f, double r_min, double r_max) {
    auto compute = [&](double fv) -> double {
        if (fv <= 1.0) return (fv * (r_max - r_min)) + r_min;
        return fv * ((r_max + r_min) / 2.0);
    };

    using U = std::decay_t<T>;
    if constexpr (std::is_arithmetic_v<U>) {
        // numeric input
        return compute(static_cast<double>(f));
    } else {
        // string-like input ("N", "0.5", "2", etc.)
        std::string s = std::string(f);
        if (s == "N" || s == "n") return 0.0;
        try {
            double fv = std::stod(s);
            return compute(fv);
        } catch (...) {
            // If parse fails, treat as 0
            return 0.0;
        }
    }
}


// #################
// ### FUNCTIONS ###
// #################

// initial mass
inline double initial_mass(double minimum, double maximum, double fraction) {
    return minimum + fraction * (maximum - minimum);
}

// total across reservoirs
inline double reservoir_totals(double M, double DM, double F, double SS, double S) {
    return M + DM + F + SS + S;
}

// mass balance as a fraction of initial
inline double mass_balance(double initial, double now) {
    return (initial - now) / initial;  
}

// Sulfur concentration from total sulfur in moles and total mass of reservoir in kg
inline double S_conc(double ST, double mass) {
    return ((ST * 32.)/1000.)/ mass;  
}

// mole S/s from mantle melting
inline double molSs_mm(double mag_kgs, double mag_S, double f) {
    return ((f * mag_kgs * mag_S) * 1000.) / 32.0;  // mm
}

// reservoir composition
inline tuple<double,double,double,double,double,double,double,double>
reservoir_R_S(double S32, double S33, double S34, double S36) {
    if (S32 == 0.0) {
        return {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    } else {
        double R33 = S33 / S32; // 33S/32S
        double R34 = S34 / S32; // 34S/32S
        double R36 = S36 / S32; // 36S/32S
        double ST  = S32 + S33 + S34 + S36;
        double T32 = S32 / ST;  // 32S/ST
        double T33 = S33 / ST;  // 33S/ST
        double T34 = S34 / ST;  // 34S/ST
        double T36 = S36 / ST;  // 36S/ST
        return {R33, R34, R36, T32, T33, T34, T36, ST};
    }
}

// delta and cap-delta composition for reservoir
inline tuple<double,double,double,double,double,double,double,double,double,double>
reservoir_isotope(double S32, double S33, double S34, double S36,
                  const map<string,double>& VCDT,
                  const string& MIF_eq)
{
    if (S32 == 0.0) {
        // R33, R34, R36, d33, d34, d36, D33, D36, ST, log_ST
        return {0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0, 0.0,0.0};
    }

    auto [R33, R34, R36, T32, T33, T34, T36, ST] = reservoir_R_S(S32, S33, S34, S36);
    double log_ST = (ST > 0.0) ? log10(ST) : 0.0;

    double d33 = R2d(R33, 33, VCDT);
    double d34 = R2d(R34, 34, VCDT);
    double d36 = R2d(R36, 36, VCDT);

    auto [D33, D36] = MIF(R33, R34, R36, VCDT, MIF_eq);

    // R33, R34, R36, d33, d34, d36, D33, D36, ST, log_ST
    return {R33, R34, R36, d33, d34, d36, D33, D36, ST, log_ST};
}


// calculate initial isotopic composition of a reservoir
// Returns:
//  log_ST, R33, R34, R36, T32, T33, T34, T36, D33, D36, S32, S33, S34, S36
inline tuple<double,double,double,double,double,double,double,double,
             double,double,double,double,double,double>
initial_reservoir(double ST, double d33, double d34, double d36,
                  const map<string,double>& VCDT,
                  const string& MIF_eq)
{
    if (ST > 0.0) {
        double log_ST = log10(ST);

        // delta -> ratios
        double R33 = d2R(d33, 33, VCDT);
        double R34 = d2R(d34, 34, VCDT);
        double R36 = d2R(d36, 36, VCDT);

        // ratios -> isotope fractions over total (Ti = iS / ST)
        double T32 = ROT(1.0,  R33, R34, R36);
        double T33 = ROT(R33,  R34, R36, 1.0);
        double T34 = ROT(R34,  R36, 1.0, R33);
        double T36 = ROT(R36,  1.0, R33, R34);

        // mass-independent anomalies (cap-deltas)
        auto [D33, D36] = MIF(R33, R34, R36, VCDT, MIF_eq);

        // absolute isotope amounts (moles) from fractions
        double S32 = ST * T32;
        double S33 = ST * T33;
        double S34 = ST * T34;
        double S36 = ST * T36;

        return {log_ST, R33, R34, R36, T32, T33, T34, T36, D33, D36, S32, S33, S34, S36};
    } else {
        // Python returns empty strings + zeros here; in C++ we return zeros
        return {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    }
}

// composition due to a process - mo, rf, sl, sq, hg
// Returns: R33, R34, R36, S32, S33, S34, S36
inline tuple<double,double,double,double,double,double,double>
process_R_S(double F,
            double res_33R, double res_34R, double res_36R,
            double a33,     double a34,     double a36)
{
    if (F == 0.0) {
        return {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    } else {
        double R33 = res_33R * a33;
        double R34 = res_34R * a34;
        double R36 = res_36R * a36;

        double S32 = F * ROT(1.0,  R33, R34, R36);
        double S33 = F * ROT(R33,  R34, R36, 1.0);
        double S34 = F * ROT(R34,  R36, 1.0, R33);
        double S36 = F * ROT(R36,  1.0, R33, R34);

        return {R33, R34, R36, S32, S33, S34, S36};
    }
}

// mass balance
inline double mantle(double M, double mo, double pl, double rs, double dm) 
{
    return M - mo - pl + rs + dm; // mantle outgassing, plutons, return silicate+sulfate, deep mantle mantle
}

// deep mantle mixing
inline double deep_mantle(double DM, double dm) {
    return DM - dm; // deep mantle mixing
}
// crustal frost
inline double frost(double F, double fr, double bu, double hg) {
    return F - fr + bu + hg; // frost remobilisation, burial
}

// silicate and sulfate sequestration
inline double silsulf(double SS, double rs, double sq, double pl) {
    return SS - rs + sq + pl; // return silicate+sulfate, sequestration, pluton
}

// space loss
inline double space(double S, double pu) {
    return S + pu; // space loss
}

// burial
inline double a_bu(double sn_F, double pd_F, double a_vp, double a_pd) {
    double bu_F = sn_F + pd_F; 
    return a_vp * (sn_F / bu_F) + a_pd * (pd_F / bu_F);
}

// plasma reactions
inline double a_pr(const unordered_map<string, double>& alphas,
                   const unordered_map<string, double>& rates) {
    // Expect required keys to exist; use .at() so missing keys throw and are easy to catch.
    const double pu = alphas.at("pu");
    const double gr = alphas.at("gr");
    const double th = alphas.at("th");
    const double ed_a = alphas.at("ed");
    const double ac_a = alphas.at("ac");
    const double ei_a = alphas.at("ei");
    const double rc_a = alphas.at("rc");
    const double ec_a = alphas.at("ec");
    const double pi_a = alphas.at("pi");

    const double pr_r = rates.at("pr");
    if (pr_r == 0.0) return 0.0;            // avoid divide-by-zero; define as 0

    const double ed_r = rates.at("ed") / pr_r;
    const double ac_r = rates.at("ac") / pr_r;
    const double ei_r = rates.at("ei") / pr_r;
    const double rc_r = rates.at("rc") / pr_r;
    const double ec_r = rates.at("ec") / pr_r;
    const double pi_r = rates.at("pi") / pr_r;

    return pu * ( gr * th * ( ed_a * ed_r + ac_a * ac_r + ei_a * ei_r + rc_a * rc_r + ec_a * ec_r )
                + pi_a * pi_r );
}



// ==================== Newton–Raphson (fast) ==========================

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
 *   per_iter_out – optional stream (&cout for CSV each iteration; nullptr for quiet)
 *   print_header – if true, print CSV header before per-iteration output
 */
template <typename Eq, typename Deriv, typename Constants>
inline double newton_raphson(double x0,
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

    if (per_iter_out) { // optional initial row like your Python version
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


// --------------------------
// fractionate() — C++ port
// --------------------------
// partition material between space and burial
struct FractionateConsts {
    // Matches Python tuple: (STT, STE, STC, S32T, S33T, S34T, S36T, a33, a34, a36)
    double STT, STE, STC;
    double S32T, S33T, S34T, S36T;
    double a33, a34, a36;
};

// Returns: S32C, S33C, S34C, S36C, S32E, S33E, S34E, S36E
inline tuple<double,double,double,double,double,double,double,double>
fractionate(const FractionateConsts& c,
            double nr_step,
            double nr_tol,
            double guessx,
            ostream* per_iter_out = nullptr,
            bool print_header = true)
{
    auto crust_comp = [&](double S32C) {
        // denominators appear multiple times; compute once
        const double den33 = S32C + c.S32T * c.a33 - S32C * c.a33;
        const double den34 = S32C + c.S32T * c.a34 - S32C * c.a34;
        const double den36 = S32C + c.S32T * c.a36 - S32C * c.a36;

        // minimal guard to avoid UB if a denominator hits 0 exactly
        // (real data shouldn’t land here; NR should avoid)
        const double eps = 1e-30;
        const double S33C = (c.S33T * S32C) / (fabs(den33) < eps ? copysign(eps, den33) : den33);
        const double S34C = (c.S34T * S32C) / (fabs(den34) < eps ? copysign(eps, den34) : den34);
        const double S36C = (c.S36T * S32C) / (fabs(den36) < eps ? copysign(eps, den36) : den36);
        return tuple<double,double,double>(S33C, S34C, S36C);
    };

    auto f = [&](double S32C, const FractionateConsts&) {
        auto [S33C, S34C, S36C] = crust_comp(S32C);
        const double mb = c.STC - S32C - S33C - S34C - S36C;
        return mb;
    };

    // auto df = [&](double S32C, const FractionateConsts&) {
    //     // Python:
    //     // -S32C*S33T*(a33 - 1)/(-S32C*a33 + S32C + S32T*a33)^2
    //     // -S32C*S34T*(a34 - 1)/(-S32C*a34 + S32C + S32T*a34)^2
    //     // -S32C*S36T*(a36 - 1)/(-S32C*a36 + S32C + S32T*a36)^2
    //     // -S33T/(-S32C*a33 + S32C + S32T*a33)
    //     // -S34T/(-S32C*a34 + S32C + S32T*a34)
    //     // -S36T/(-S32C*a36 + S32C + S32T*a36)
    //     // -1
    //     const double d33 = (-S32C * c.a33 + S32C + c.S32T * c.a33);
    //     const double d34 = (-S32C * c.a34 + S32C + c.S32T * c.a34);
    //     const double d36 = (-S32C * c.a36 + S32C + c.S32T * c.a36);

    //     const double eps = 1e-30;
    //     const double D33 = fabs(d33) < eps ? copysign(eps, d33) : d33;
    //     const double D34 = fabs(d34) < eps ? copysign(eps, d34) : d34;
    //     const double D36 = fabs(d36) < eps ? copysign(eps, d36) : d36;

    //     const double term33 = -S32C * c.S33T * (c.a33 - 1.0) / (D33 * D33) - c.S33T / D33;
    //     const double term34 = -S32C * c.S34T * (c.a34 - 1.0) / (D34 * D34) - c.S34T / D34;
    //     const double term36 = -S32C * c.S36T * (c.a36 - 1.0) / (D36 * D36) - c.S36T / D36;

    //     return term33 + term34 + term36 - 1.0;
    // };
    // df(S32C) — derivative, copy-paste checkable (matches Python term-by-term)
    auto df = [&](double S32C, const FractionateConsts&) {
        return
            // -S32C*S33T*(a33 - 1.)/(-S32C*a33 + S32C + S32T*a33)**2.
            - S32C * c.S33T * (c.a33 - 1.0)
                / pow(-S32C * c.a33 + S32C + c.S32T * c.a33, 2)

            // -S32C*S34T*(a34 - 1.)/(-S32C*a34 + S32C + S32T*a34)**2.
            - S32C * c.S34T * (c.a34 - 1.0)
                / pow(-S32C * c.a34 + S32C + c.S32T * c.a34, 2)

            // -S32C*S36T*(a36 - 1.)/(-S32C*a36 + S32C + S32T*a36)**2.
            - S32C * c.S36T * (c.a36 - 1.0)
                / pow(-S32C * c.a36 + S32C + c.S32T * c.a36, 2)

            // -S33T/(-S32C*a33 + S32C + S32T*a33)
            - c.S33T
                / (-S32C * c.a33 + S32C + c.S32T * c.a33)

            // -S34T/(-S32C*a34 + S32C + S32T*a34)
            - c.S34T
                / (-S32C * c.a34 + S32C + c.S32T * c.a34)

            // -S36T/(-S32C*a36 + S32C + S32T*a36)
            - c.S36T
                / (-S32C * c.a36 + S32C + c.S32T * c.a36)

            // -1
            - 1.0;
    };


    // Solve for S32C
    const double S32C = newton_raphson(guessx, c, nr_tol, nr_step, f, df,
                                       /*maxiter*/1000, per_iter_out, print_header);

    // Build reservoirs from S32C
    auto [S33C, S34C, S36C] = crust_comp(S32C);
    const double S32E = c.S32T - S32C;
    const double S33E = c.S33T - S33C;
    const double S34E = c.S34T - S34C;
    const double S36E = c.S36T - S36C;

    return {S32C, S33C, S34C, S36C, S32E, S33E, S34E, S36E};
}


////////////////////////////////////
//// CALCULATIONS  /////////////////
////////////////////////////////////

// FUNCTION: main function to calculate isotope evolution
// iso_evo() - main function to calculate isotope evolution

///////////////////////////////////////
////// RAYLEIGH BUT BOX MODEL //////
/////////////////////////////////////////

// FUNCTION: simplified isotope evolution function
// iso_evo_simple() - simplified isotope evolution function 

/////////////////////////////////////////
///// ATMOSPHERIC PROFILE /////////////
/////////////////////////////////////////

// FUNCTION: calculate atmospheric profile
// calc_atm() - calculate atmospheric profile

///////////
// FUNCTION: range of fractionatinon factors at a single step:
// range_frac_fac_ss()

// FUNCTION: range of fractionation factors for various proportions of space to burial:
// range_frac_fac_s2b()

////////////////////////////////////
// core mantle calculations ////////
////////////////////////////////////

// FUNCTION:
// ratio2ototal(ratio)

// FUNCTION: 34S/32S of mantle and core given mantle-core ratio
// ratio_core2mantle() 

// FUNCTION: isotopic composition for variable S fractions between core and mantle
// calc_iso_variable()








//////////////////////////////////////////////////////
#endif

