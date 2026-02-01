#ifndef ISOTOPE_EVOLUTION_HPP
#define ISOTOPE_EVOLUTION_HPP
#include <unordered_map>
#include <cmath>
#include <map>
#include <string>
#include <type_traits>
#include <tuple>      // <-- needed by tuple returns
#include <iostream>   // <-- needed for cerr used in R2d
#include "alphaset_types.hpp"  // <- this is where AlphaSet lives
#include "alphaset.hpp"        // <- this is where compute_all_alphas() lives, needed to get all alphas at once
#include "helpers.hpp"     // <- this is where MaybeDouble lives to handle "N" or numeric f inputs
#pragma once
#include <vector>
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
const double seconds_per_year = 365.0 * 24.0 * 60.0 * 60.0; // seconds in a year

// canonical mass-dependent fractionation powers
const double e33 = 0.515;
const double e36 = 1.9;


// ===========================================================================
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
    // cout << "Inside frost(): F = " << F << ", fr = " << fr << ", bu = " << bu << ", hg = " << hg << endl;
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
void iso_evo(
    double t_step_Myr,
    double end_time_Myr,
    double nr_step,
    double nr_tol,

    // deep mantle variables
    MaybeDouble DM_mass_f,
    MaybeDouble DM_mass_f_sil,
    MaybeDouble DM_ST,
    double DM_33d, double DM_34d, double DM_36d,

    // mantle variables
    MaybeDouble M_mass_f,
    MaybeDouble M_mass_f_sil,
    MaybeDouble M_ST,
    double M_33d, double M_34d, double M_36d,

    // crustal frost variables
    MaybeDouble F_mass_f,
    MaybeDouble F_ST,
    double F_33d, double F_34d, double F_36d,

    // crustal sulfate variables
    MaybeDouble SS_mass_f,
    MaybeDouble SS_ST,
    double SS_33d, double SS_34d, double SS_36d,

    MaybeDouble S_ST,
    double S_33d, double S_34d, double S_36d,

    const unordered_map<string, double>& rate_f,
    const string& oscillate,

    double resurf_cm_yr,
    double sil_mag_S,
    double thick_C,

    double f_S2,
    double f_pl2mo,
    double f_sq,
    double f_deep,
    double f_remobilised,
    double f_pu,

    // # options for calculating fractionation factors
    // nofrac = "no" # if "yes", all fractionation factors = 1 (i.e., no fractionation); if "no", fractionation factors are pre-determined values
    // MAF = "no" # if "yes", mass-anomalous (i.e., with powers not equal to canonical values) is possible; if "no", fractionation factors follow canonical power laws
    const string& nofrac, const string& MAF, 
    const string& MIF_eq
    
)
{
    // === set up vector for results ===
    vector<SulfurState> results;  

    // === TIMESTAMP ===
    auto now = chrono::system_clock::to_time_t(chrono::system_clock::now());
    cout << endl << "======================================================" << endl;
    cout <<         "=== Starting sulfur isotope evolution calculation ====" << endl;
    cout << "======================================================" << endl;
    cout << "Start time: " << std::put_time(std::localtime(&now), "%Y-%m-%d %H:%M:%S") << endl;


    // === TIME SETUP ===
    double time_step = t_step_Myr * 60.0 * 60.0 * 24.0 * 365.0 * 1.0e6;  // seconds per timestep
    int end = static_cast<int>((end_time_Myr + t_step_Myr) / t_step_Myr);  // number of steps

    // === RATE CONSTANTS (mol S / s) plasma and photodissociation ===
    double pd_r = kgs2mSs(calc_rate(rate_f.at("pd"),  rate_kgs.at("pd").min,  rate_kgs.at("pd").max));
    double pi_r =         calc_rate(rate_f.at("pi"),  rate_kgs.at("pi").min,  rate_kgs.at("pi").max);
    double ei_r =         calc_rate(rate_f.at("ei"),  rate_kgs.at("ei").min,  rate_kgs.at("ei").max);
    double ed_r =         calc_rate(rate_f.at("ed"),  rate_kgs.at("ed").min,  rate_kgs.at("ed").max);
    double ac_r =         calc_rate(rate_f.at("ac"),  rate_kgs.at("ac").min,  rate_kgs.at("ac").max);
    double rc_r =         calc_rate(rate_f.at("rc"),  rate_kgs.at("rc").min,  rate_kgs.at("rc").max);
    double ec_r =         calc_rate(rate_f.at("ec"),  rate_kgs.at("ec").min,  rate_kgs.at("ec").max);

    // === PLASMA REACTIONS & FRACTIONATION ===
    double pr_r = pi_r + ei_r + ed_r + ac_r + rc_r + ec_r; // total plasma reaction rate is sum of pi_r, ei_r, ed_r, ac_r, rc_r, ec_r

    unordered_map<string, double> rates_plasma = { //this maps "pr" and "pi" etc to their rates (plasma rates)
        {"pr", pr_r},
        {"pi", pi_r},
        {"ei", ei_r},
        {"ed", ed_r},
        {"ac", ac_r},
        {"rc", rc_r},
        {"ec", ec_r}
    };

    // get Set of all Alpha vals struct AlphaSet from header function, useful for the alphas called next
    AlphaSet alpha = compute_all_alphas(nofrac, MAF);

    unordered_map<string, double> a33_plasma = { // this maps "pu" and "gr" etc to their fractionation factors for a33)
        {"pu", alpha.a33["pu"]}, {"gr", alpha.a33["gr"]}, {"th", alpha.a33["th"]},
        {"ed", alpha.a33["ed"]}, {"ac", alpha.a33["ac"]}, {"ei", alpha.a33["ei"]},
        {"rc", alpha.a33["rc"]}, {"ec", alpha.a33["ec"]}, {"pi", alpha.a33["pi"]}
    };

    unordered_map<string, double> a34_plasma = { // this maps "pu" and "gr" etc to their fractionation factors for a34
        {"pu", alpha.a34["pu"]}, {"gr", alpha.a34["gr"]}, {"th", alpha.a34["th"]},
        {"ed", alpha.a34["ed"]}, {"ac", alpha.a34["ac"]}, {"ei", alpha.a34["ei"]},
        {"rc", alpha.a34["rc"]}, {"ec", alpha.a34["ec"]}, {"pi", alpha.a34["pi"]}
    };

    unordered_map<string, double> a36_plasma = { // this maps "pu" and "gr" etc to their fractionation factors for a36
        {"pu", alpha.a36["pu"]}, {"gr", alpha.a36["gr"]}, {"th", alpha.a36["th"]},
        {"ed", alpha.a36["ed"]}, {"ac", alpha.a36["ac"]}, {"ei", alpha.a36["ei"]},
        {"rc", alpha.a36["rc"]}, {"ec", alpha.a36["ec"]}, {"pi", alpha.a36["pi"]}
    };

    double a33_pr = a_pr(a33_plasma, rates_plasma); // weighted average fractionation factor for 33S in plasma reactions
    double a34_pr = a_pr(a34_plasma, rates_plasma); // weighted average fractionation factor for 34S in plasma reactions
    double a36_pr = a_pr(a36_plasma, rates_plasma); // weighted average fractionation factor for 36S in plasma reactions

    double pu_r = calc_rate(f_pu, rate_kgs.at("pu").min, pr_r); // pu rate in mol S / s
    pu_r = kgs2mSs(pu_r); // convert pu rate to mol S / s
    pr_r = kgs2mSs(pr_r); // convert pr rate to mol S / s

    double pu_frac = (pr_r != 0.0) ? pu_r / pr_r : 0.0; // calculate pu fraction of total plasma reactions

    // -------------------------
    // Mass of silicate crust and fraction of frost remobilised
    // -------------------------
    double radius = 1822.6 * 1.e3;  // IO RADIUS (in meters
    double mass_sil_C = (4.0 / 3.0) * M_PI * 
        (pow(radius, 3) - pow(radius - thick_C, 3)) * sil_mag_rho;
    
    // -------------------------
    // MANTLE INITIAL RESERVOIR: INITIAL RESERVOIR NUMBERS
    // -------------------------

    // Resolve M_ST if still flagged as "N"
    if (M_ST.flag == "N") {
        M_ST.value = M_mass_i_min + M_mass_f.value * (M_mass_i_max - M_mass_i_min);
        M_ST.flag = "Y";  // Optional: mark as resolved
    }

    // Compute silicate mantle mass (subtracting crust mass)
    double M_mass = (M_mass_i_min_sil + M_mass_f_sil.value * (M_mass_i_max_sil - M_mass_i_min_sil)) - mass_sil_C;

    // Sulfur concentration in mantle (mol/kg)
    double S_conc_M = S_conc(M_ST.value, M_mass);

    // Call initial_reservoir for mantle
    double logM_ST, M_33R, M_34R, M_36R;
    double M_32T_i, M_33T_i, M_34T_i, M_36T_i;
    double M_33D, M_36D;
    double M_32S, M_33S, M_34S, M_36S;

    tie(logM_ST, M_33R, M_34R, M_36R,
        M_32T_i, M_33T_i, M_34T_i, M_36T_i,
        M_33D, M_36D,
        M_32S, M_33S, M_34S, M_36S) = initial_reservoir(M_ST.value, M_33d, M_34d, M_36d, VCDT, MIF_eq);


    // -------------------------
    // DEEP MANTLE INITIAL RESERVOIR
    // -------------------------

    // Resolve DM_ST if flagged as "N"
    if (DM_ST.flag == "N") {
        if (DM_mass_f.flag != "N") {
            DM_ST.value = M_mass_i_min + DM_mass_f.value * (M_mass_i_max - M_mass_i_min);
        } else {
            DM_ST.value = 0.0;
        }
        DM_ST.flag = "Y";
    }

    // Resolve DM_mass using silicate-only bounds
    double DM_mass;
    if (DM_mass_f_sil.flag != "N") {
        DM_mass = M_mass_i_min_sil + DM_mass_f_sil.value * (M_mass_i_max_sil - M_mass_i_min_sil);
    } else {
        DM_mass = 0.0;
    }

    // Call initial_reservoir for deep mantle
    double logDM_ST, DM_33R, DM_34R, DM_36R;
    double DM_32T_i, DM_33T_i, DM_34T_i, DM_36T_i;
    double DM_33D, DM_36D;
    double DM_32S, DM_33S, DM_34S, DM_36S;

    tie(logDM_ST, DM_33R, DM_34R, DM_36R,
        DM_32T_i, DM_33T_i, DM_34T_i, DM_36T_i,
        DM_33D, DM_36D,
        DM_32S, DM_33S, DM_34S, DM_36S) = initial_reservoir(DM_ST.value, DM_33d, DM_34d, DM_36d, VCDT, MIF_eq);
    
    // -------------------------
    // CRUSTAL FROSTS INITIAL RESERVOIR
    // -------------------------
    if (F_ST.flag == "N") {
        F_ST.value = F_mass_i_min + F_mass_f.value * (F_mass_i_max - F_mass_i_min);
        F_ST.flag = "Y";
    }
    
    double logF_ST, F_33R, F_34R, F_36R;
    double F_32T_i, F_33T_i, F_34T_i, F_36T_i;
    double F_33D, F_36D;
    double F_32S, F_33S, F_34S, F_36S;

    tie(logF_ST, F_33R, F_34R, F_36R,
        F_32T_i, F_33T_i, F_34T_i, F_36T_i,
        F_33D, F_36D,
        F_32S, F_33S, F_34S, F_36S) = initial_reservoir(F_ST.value, F_33d, F_34d, F_36d, VCDT, MIF_eq);

    // -------------------------
    // CRUSTAL SULFATES INITIAL RESERVOIR
    // -------------------------
    if (SS_ST.flag == "N") {
        SS_ST.value = SS_mass_i_min + SS_mass_f.value * (SS_mass_i_max - SS_mass_i_min);
        SS_ST.flag = "Y";
    }

    double logSS_ST, SS_33R, SS_34R, SS_36R;
    double SS_32T_i, SS_33T_i, SS_34T_i, SS_36T_i;
    double SS_33D, SS_36D;
    double SS_32S, SS_33S, SS_34S, SS_36S;

    tie(logSS_ST, SS_33R, SS_34R, SS_36R,
        SS_32T_i, SS_33T_i, SS_34T_i, SS_36T_i,
        SS_33D, SS_36D,
        SS_32S, SS_33S, SS_34S, SS_36S) = initial_reservoir(SS_ST.value, SS_33d, SS_34d, SS_36d, VCDT, MIF_eq);


    // -------------------------
    // TOTAL SURFACE SULFUR INITIAL RESERVOIR
    // -------------------------
    if (S_ST.flag == "N") {
        S_ST.value = 0.0;
        S_ST.flag = "Y";
    }

    double logS_ST, S_33R, S_34R, S_36R;
    double S_32T_i, S_33T_i, S_34T_i, S_36T_i;
    double S_33D, S_36D;
    double S_32S, S_33S, S_34S, S_36S;

    tie(logS_ST, S_33R, S_34R, S_36R,
        S_32T_i, S_33T_i, S_34T_i, S_36T_i,
        S_33D, S_36D,
        S_32S, S_33S, S_34S, S_36S) = initial_reservoir(S_ST.value, S_33d, S_34d, S_36d, VCDT, MIF_eq);

    // -------------------------
    // MASS BALANCE CALCULATIONS - INITIAL STATE
    // -------------------------
    double Io_ST_initial = reservoir_totals(M_ST.value, DM_ST.value, F_ST.value, SS_ST.value, S_ST.value);
    double Io_32S_initial = reservoir_totals(M_32S, DM_32S, F_32S, SS_32S, S_32S);
    double Io_33S_initial = reservoir_totals(M_33S, DM_33S, F_33S, SS_33S, S_33S);
    double Io_34S_initial = reservoir_totals(M_34S, DM_34S, F_34S, SS_34S, S_34S);
    double Io_36S_initial = reservoir_totals(M_36S, DM_36S, F_36S, SS_36S, S_36S);

    // Mass balance checks (all using initial totals against themselves at t = 0)
    double mbST   = mass_balance(Io_ST_initial, Io_ST_initial);
    double mb32S  = mass_balance(Io_32S_initial, Io_32S_initial);
    double mb33S  = mass_balance(Io_33S_initial, Io_33S_initial);
    double mb34S  = mass_balance(Io_34S_initial, Io_34S_initial);
    double mb36S  = mass_balance(Io_36S_initial, Io_36S_initial);


    // // SETUP OUTPUT FILE: 
    ofstream out("time_evolution_cpp.txt");
    // if (!out) throw runtime_error("Could not open output file.");

    // // Build the Python-like header (determines total columns)
    vector<string> header = sulfur_python_like_header();

    // // Write the first 3 rows (rate-factor labels, values, then main header)
    write_python_like_preamble(out, header,
        rate_f.at("pd"), rate_f.at("pi"), rate_f.at("ei"), rate_f.at("ed"), rate_f.at("ac"),
        rate_f.at("rc"), rate_f.at("ec"),
        resurf_cm_yr, sil_mag_S, thick_C,
        f_remobilised, f_S2, f_pl2mo, f_sq, f_deep,
        M_mass
    );




    // -------------------------
    // INITIAL TIME STEP
    // -------------------------

    int n = 0;
    double t_s = static_cast<double>(n) * t_step_Myr;     // time in Myr
    double t_s_Myr = t_s;                                  // (duplicate for clarity with Python)
    
    // // --- write initial timestep row (n=0) like Python ---
    vector<string> init_row = sulfur_initial_row_python_like(
        0.0, t_s_Myr,
        S_conc_M,
        M_34d, M_33d, M_36d,
        S_34d, S_33d, S_36d,
        M_33D, M_36D,
        S_33D, S_36D,
        M_ST.value, S_ST.value,
        logM_ST, logS_ST,
        mbST
    );

    write_csv_row_padded(out, init_row, header.size());


    // IN PYTHON THE RESULTS ARE STORED IN A DATAFRAME LIKE THIS: 
    // # set up results table
    // results = pd.DataFrame([["rate factor pd","rate factor pi","rate factor ei","rate factor ed","rate factor ac",
    //                      "rate factor rc", "rate factor ec","resurfacing rate (cm/yr)","sil mag S (wf)", "crustal thickness (m)",
    //                      "fraction of f_crust_returned that is remolised","S2 to SO2","plutons from mantle melting","sulfate sequestration","deep mantle - fraction of mantle melting","mass silicate mantle"]])
    // results2 = pd.DataFrame([[rate_f["pd"],rate_f["pi"],rate_f["ei"],rate_f["ed"],rate_f["ac"],rate_f["rc"], 
    //                     rate_f["ec"],resurf_cm_yr,sil_mag_S, thick_C,f_remobilised,f_S2,f_pl2mo,f_sq,f_deep,M_mass]])
    // results = pd.concat([results, results2], ignore_index=True)
    // results2 = pd.DataFrame([["time step","time (Myr)","[S] mantle",
    //                      "d34S M","d34S F","d34S SS","d34S S","d34S DM","d34S og","d34S iS",
    //                      "d33S M","d33S F","d33S SS", "d33S S","d33S DM","d33S og","d33S iS",
    //                      "d36S M","d36S F","d36S SS","d36S S","d36S DM","d36S og","d36S iS",
    //                      "D33S M","D33S F","D33S SS","D33S S","D33S DM","D33S og","D33S iS",
    //                      "D36S M","D36S F","D36S SS","D36S S","D36S DM","D36S og","D36S iS",
    //                      "ST M","ST F", "ST SS","ST S","ST DM",
    //                      "log ST M","log ST F", "log ST SS","log ST S","log ST DM",
    //                      "F mo", "F fr", "F pr","F pu", "F sn", "F pd", "F hg", "F bu", "F pl","F rs","F sq","F ao","F dm",
    //                      "mo/mo+ca", "bu/bu+pl", 
    //                      "mbST","mb32","mb33","mb34","mb36","",
    //                      "R33 M", "R33 F", "R33 SS","R33 S","R33 DM","R33 og","R33 iS",
    //                      "R34 M", "R34 F", "R34 SS","R34 S","R34 DM","R34 og","R34 iS",
    //                      "R36 M", "R36 F", "R36 SS","R36 S","R36 DM","R36 og","R36 iS",
    //                      "a33_bu", "a34_bu", "a36_bu",
    //                      "a33_escape-burial", "a34_escape-burial", "a36_escape-burial"]])
    // results = pd.concat([results, results2], ignore_index=True)
    // results2 = pd.DataFrame([[n,t_s,S_conc_M,
    //                      M_34d,F_34d,SS_34d,S_34d,DM_34d,"","",
    //                      M_33d,F_33d,SS_33d,S_33d,DM_33d,"","",
    //                      M_36d,F_36d,SS_36d,S_36d,DM_36d,"","",
    //                      M_33D,F_33D,SS_33D,S_33D,DM_33D,"","",
    //                      M_36D,F_36D,SS_36D,S_36D,DM_36D,"","",
    //                      M_ST,F_ST,SS_ST,S_ST,DM_ST,
    //                      logM_ST,logF_ST,logSS_ST,logS_ST,logDM_ST,
    //                      "","","","","","","","","","","","","",
    //                      "","",
    //                      mbST,mb32S,mb33S,mb34S,mb36S,
    //                      "",
    //                      M_33R, F_33R, SS_33R, S_33R, DM_33R,"","",
    //                      M_34R, F_34R, SS_34R, S_34R, DM_36R,"","",
    //                      M_36R, F_36R, SS_36R, S_36R, DM_36R,"","",
    //                      "","","","","","","","",""]])         
    // results = pd.concat([results, results2], ignore_index=True)
    // results.to_csv('time_evolution.csv', index=False, header=False)
    
    

    // NOW THE BIG LOOP:
    // Purpose of loop: evolve sulfur reservoirs over time
    // vector<SulfurState> results; // to hold results at each timestep
    // Loop sections: 
    // 1. timesteps
    // 2. oscillating resufacing rate
    // 3. mantle and deep mantle melting (rates)
    // 4. crustal and atmospheric rates
    // 5. pl and rs fluxes
    // 6. flux for snow
    // 7. rates conversion to fluxes -- converted to flux (mol per timestep)
    // 8. "fractions", giving ratios of volcanic S from outgassing and S from atmosphere vs total
    // 9. burial and escape (effects on 33,34,36 fractionation)
    // 10. mo, rf, dm isotope inputs. 
    // 11. total input calculation 
    // 12. homogenous gas
    // 13. sulfate sequestration
    // 14. snow if negative
    // 15. atmospheric outgassing
    // 16. reservoirs fractionate
    // 17. instatntaneous isotope ratios of space reservoir
    // 18. plutonic and silicate/sulfate return
    // 19. update reservoirs (M_XYS, DM_XYS, F_XYS, SS_XYS, S_XYS)
    // 20. mass balance checks
    // 21. store results in struct and push to vector
    

    // ============================================================================
    // ============================================================================
    // Sulfur model main time loop (C++ version)
    // Here we go...
    // ============================================================================
    // ============================================================================
    // cout << "HERE we start time looping" << endl;
    int total_steps = static_cast<int>(end_time_Myr / t_step_Myr);
    for (int n = 1; n <= total_steps; ++n) { //updated to start at n=1 to match python version
    // for (int n = 1; n < end; ++n) {
        // --------------------
        // 1. TIMESTEPS
        // --------------------
        double t_s = n * time_step;
        double t_s_Myr = n * t_step_Myr;
        
        // --------------------
        // debugging
        // if (n < 5) {
        //     ostringstream oss_time, oss_F34d;
        //     oss_time << fixed << setprecision(17) << t_s_Myr;
        //     oss_F34d << fixed << setprecision(17) << F_34d;
            
        //     cout << endl;
        //     cout << "Time: " << oss_time.str() << " Myr"
        //         << " | F_34d: " << oss_F34d.str()
        //         << " | Time: " << current_datetime_string() << endl;
            // DEBUG_PRINT("time step (Myr): ", t_s_Myr);
            // DEBUG_PRINT("F_34d: ", F_34d);
            // DEBUG_PRINT("M_ST: ", M_ST.value);
            // DEBUG_PRINT("DM_ST: ", DM_ST.value);
            // DEBUG_PRINT("F_ST: ", F_ST.value);
            // DEBUG_PRINT("SS_ST: ", SS_ST.value);
            // DEBUG_PRINT("S_ST: ", S_ST.value);
        // }
       
        
        // --------------------

        // // periodic output:
        // if (n % 1000 == 0) {
        //     cout << "Time: " << t_s_Myr << " Myr"
        //          << " | F_34d: " << F_34d
        //          << " | Time: " << current_datetime_string() << endl; 
        // }
        // --------------------
        // 2. OSCILLATING RESURFACING RATE
        // --------------------
        if (oscillate == "Y") {
            if (t_s_Myr >= 4400.0) {
                resurf_cm_yr = 5.0;
                }
            else if (
                (t_s_Myr >= 100.0 && t_s_Myr < 200.0) || (t_s_Myr >= 300.0 && t_s_Myr < 400.0) ||
                (t_s_Myr >= 500.0 && t_s_Myr < 600.0) || (t_s_Myr >= 700.0 && t_s_Myr < 800.0) ||
                (t_s_Myr >= 900.0 && t_s_Myr < 1000.0) || (t_s_Myr >= 1100.0 && t_s_Myr < 1200.0) ||
                (t_s_Myr >= 1300.0 && t_s_Myr < 1400.0) || (t_s_Myr >= 1500.0 && t_s_Myr < 1600.0) ||
                (t_s_Myr >= 1700.0 && t_s_Myr < 1800.0) || (t_s_Myr >= 1900.0 && t_s_Myr < 2000.0) ||
                (t_s_Myr >= 2100.0 && t_s_Myr < 2200.0) || (t_s_Myr >= 2300.0 && t_s_Myr < 2400.0) ||
                (t_s_Myr >= 2500.0 && t_s_Myr < 2600.0) || (t_s_Myr >= 2700.0 && t_s_Myr < 2800.0) ||
                (t_s_Myr >= 2900.0 && t_s_Myr < 3000.0) || (t_s_Myr >= 3100.0 && t_s_Myr < 3200.0) ||
                (t_s_Myr >= 3300.0 && t_s_Myr < 3400.0) || (t_s_Myr >= 3500.0 && t_s_Myr < 3600.0) ||
                (t_s_Myr >= 3700.0 && t_s_Myr < 3800.0) || (t_s_Myr >= 3900.0 && t_s_Myr < 4000.0) ||
                (t_s_Myr >= 4100.0 && t_s_Myr < 4200.0) || (t_s_Myr >= 4300.0 && t_s_Myr < 4400.0)
                ) {
                    resurf_cm_yr = 1.0;
                }
            else {
                resurf_cm_yr = 9.0;
                }
            }
        // DEBUG_PRINT("resurf_cm_yr: ", resurf_cm_yr);
        // --------------------
        // 3. MANTLE AND DEEP MANTLE MELTING (RATES)
        // --------------------
        // Silicate magmatism in m/yr (adjusting for plutonic retention)
        double sil_mag_m_yr = (resurf_cm_yr / (1.0 - f_pl2mo)) * 0.01;  // silicate magmatism in cm/yr → m/yr
        double sil_mag_m_s = sil_mag_m_yr / seconds_per_year;          // silicate magmatism in m/s
        double sil_mag_kg_s = sil_mag_m_s * Io_SA * sil_mag_rho;       // silicate magmatism mass in kg/s volume flux * density

        // Fraction of crust returned in one time step
        double f_crust_return = (sil_mag_m_s * time_step) / thick_C; // volume fraction crust returned in a time step
        double f_remob = f_remobilised * f_crust_return; // fraction of crustal frosts that are remobilised
        
        // mantle melting rate
        double S_conc_M = S_conc(M_ST.value, M_mass);  // in wt%
        double mm_r;
        
        if (S_conc_M > sil_mag_S) {
            mm_r = molSs_mm(sil_mag_kg_s, sil_mag_S, 1.0);
        } else {
            mm_r = molSs_mm(sil_mag_kg_s, S_conc_M, 1.0);
        }

        // deep mantle melting rate
        double dm_r = 0.0;

        if (DM_ST.flag == "Y" && DM_ST.value > 0.0) {
            double S_conc_DM = S_conc(DM_ST.value, DM_mass);
            if (S_conc_DM > sil_mag_S) {
                dm_r = molSs_mm(sil_mag_kg_s, sil_mag_S, f_deep);
            } else {
                dm_r = molSs_mm(sil_mag_kg_s, S_conc_DM, f_deep);
            }
        } else {
            dm_r = 0.0; // Deep mantle inactive or no sulfur present
        }
        // DEBUG_PRINT("mm_r: ", mm_r);
        // DEBUG_PRINT("dm_r: ", dm_r);

        // 4. CRUSTAL AND ATMOSPHERIC RATES
        // -------------------------------
        // DEBUG_PRINT("A fr_r: ", fr_r);
        // cout << "--------------------------------" << "\n";
        // DEBUG_PRINT("A f_remob: ", f_remob);
        // DEBUG_PRINT("A F_ST: ", F_ST.value);
        // DEBUG_PRINT("A timestep: ", time_step);

        double mo_r = (1.0 - f_pl2mo) * mm_r;  // mantle outgassing rate
        double fr_r = (f_remob * F_ST.value) / time_step;  // remobilized frost release
        double hg_r = (mo_r + fr_r) * f_S2;  // photochemical loss to space (Hg)
        double sq_r = (mo_r + fr_r) * f_sq;  // surface sequestration (Sq)
        double ao_r = (mo_r + fr_r) - (sq_r + hg_r);  // atmospheric outgassing
        double sn_r = ao_r - (pu_r + pd_r);  // net sulfur to atmosphere

        // DEBUG_PRINT("mo_r: ", mo_r);
        // DEBUG_PRINT("B fr_r: ", fr_r);
        // DEBUG_PRINT("B f_remob: ", f_remob);
        // DEBUG_PRINT("B F_ST: ", F_ST.value);
        // DEBUG_PRINT("B timestep: ", time_step);
        // DEBUG_PRINT("hg_r: ", hg_r);
        // DEBUG_PRINT("sq_r: ", sq_r);
        // DEBUG_PRINT("ao_r: ", ao_r);
        // DEBUG_PRINT("sn_r: ", sn_r);

        // 5. PLUTONIC AND REGASSING FLUXES
        // -------------------------------
        double pl_r = mm_r * f_pl2mo;  // plutonic retention rate
        double rs_r = (f_crust_return * SS_ST.value) / time_step;  // subduction (regassing) rate
        // DEBUG_PRINT("pl_r: ", pl_r);
        // DEBUG_PRINT("rs_r: ", rs_r);
        // 6. FLUX FOR SNOW
        // ----------------
        double bu_r;

        if (sn_r < 0.0) {
            bu_r = pd_r;
        } else {
            bu_r = pd_r + sn_r;
        }
        // DEBUG_PRINT("bu_r: ", bu_r);
        // 7. RATES → FLUXES (mol per timestep) .... rates converted to flux (mol per timestep)
        // -------------------------------------
        double dm_F = dm_r * time_step;
        double mm_F = mm_r * time_step;
        double mo_F = mo_r * time_step;
        double fr_F = fr_r * time_step;
        double pd_F = pd_r * time_step;
        double hg_F = hg_r * time_step;
        double sq_F = sq_r * time_step;
        double ao_F = ao_r * time_step;
        double sn_F = sn_r * time_step;
        double bu_F = bu_r * time_step;
        double pl_F = pl_r * time_step;
        double pr_F = pr_r * time_step;
        double pu_F = pu_r * time_step;
        double rs_F = rs_r * time_step;

        // 8. FRACTIONS (Volcanic vs Frost Remobilization, Burial vs Escape)
        // ------------------------------------------------------------------
        // Avoid divide-by-zero with conditional expressions
        double mo_mo_ca = (mo_r + fr_r > 0.0) ? (mo_r / (mo_r + fr_r)) : 0.0; // amount of volcanic S from mantle outgassing vs. total (inc. frost remobilisation)
        double bu_bu_pl = (bu_r + pl_r + hg_r + sq_r > 0.0) ? (bu_r / (bu_r + pl_r + hg_r + sq_r)) : 0.0; // amount of S lost from atmosphere to surface vs total (inc. space loss)

        // 9. BURIAL AND ESCAPE FRACTIONATION FACTORS - effects on 33,34,36 fractionation
        // -------------------------------------------
        // Burial fractionation (depends on presence of snow flux)
        double a33_bu = a_bu((sn_r > 0.0 ? sn_F : 0.0), pd_F, alpha.a33["vp"], alpha.a33["pd"]);
        double a34_bu = a_bu((sn_r > 0.0 ? sn_F : 0.0), pd_F, alpha.a34["vp"], alpha.a34["pd"]);
        double a36_bu = a_bu((sn_r > 0.0 ? sn_F : 0.0), pd_F, alpha.a36["vp"], alpha.a36["pd"]);

        // Escape-burial fractionation
        double a33_eb = (a33_bu != 0.0) ? (a33_pr / a33_bu) : 0.0;
        double a34_eb = (a34_bu != 0.0) ? (a34_pr / a34_bu) : 0.0;
        double a36_eb = (a36_bu != 0.0) ? (a36_pr / a36_bu) : 0.0;

        // 10. mo, rf, dm isotope inputs
        // ------------------------------
        // Mantle outgassing
        auto [mo_33R, mo_34R, mo_36R, mo_32S, mo_33S, mo_34S, mo_36S] =
            process_R_S(mo_F, M_33R, M_34R, M_36R, alpha.a33["mo"], alpha.a34["mo"], alpha.a36["mo"]);

        // Crustal remobilisation
        auto [fr_33R, fr_34R, fr_36R, fr_32S, fr_33S, fr_34S, fr_36S] =
            process_R_S(fr_F, F_33R, F_34R, F_36R, alpha.a33["fr"], alpha.a34["fr"], alpha.a36["fr"]);

        // Deep mantle mixing
        auto [dm_33R, dm_34R, dm_36R, dm_32S, dm_33S, dm_34S, dm_36S] =
            process_R_S(dm_F, DM_33R, DM_34R, DM_36R, alpha.a33["dm"], alpha.a34["dm"], alpha.a36["dm"]);


        // 11. TOTAL INPUT CALCULATION
        // ------------------------------
        // ti = total input
        double ti_32S = mo_32S + fr_32S;
        double ti_33S = mo_33S + fr_33S;
        double ti_34S = mo_34S + fr_34S;
        double ti_36S = mo_36S + fr_36S;

        double ti_TS  = ti_32S + ti_33S + ti_34S + ti_36S;

        double ti_34R = (ti_32S != 0.0) ? (ti_34S / ti_32S) : 0.0;
        double ti_33R = (ti_32S != 0.0) ? (ti_33S / ti_32S) : 0.0;
        double ti_36R = (ti_32S != 0.0) ? (ti_36S / ti_32S) : 0.0;

        // 12. Homogeneous gas
        // --------------------
        auto [hg_33R, hg_34R, hg_36R, hg_32S, hg_33S, hg_34S, hg_36S] =
            process_R_S(hg_F, ti_33R, ti_34R, ti_36R, alpha.a33["hg"], alpha.a34["hg"], alpha.a36["hg"]);

        // 13. Sulfate sequestration
        // --------------------------
        auto [sq_33R, sq_34R, sq_36R, sq_32S, sq_33S, sq_34S, sq_36S] =
            process_R_S(sq_F, ti_33R, ti_34R, ti_36R, alpha.a33["sq"], alpha.a34["sq"], alpha.a36["sq"]);


        // 14. Snow (negative only — adds frost to atmosphere if sn_r < 0)
        // ---------------------------------------------------------------
        double sn_32S, sn_33S, sn_34S, sn_36S;

        if (sn_r < 0.0) {
            auto [sn_33R, sn_34R, sn_36R, temp_sn_32S, temp_sn_33S, temp_sn_34S, temp_sn_36S] =
                process_R_S(-1.0 * sn_F, F_33R, F_34R, F_36R,
                            alpha.a33["vp"], alpha.a34["vp"], alpha.a36["vp"]);
            sn_32S = temp_sn_32S;
            sn_33S = temp_sn_33S;
            sn_34S = temp_sn_34S;
            sn_36S = temp_sn_36S;
        } else {
            sn_32S = 0.0;
            sn_33S = 0.0;
            sn_34S = 0.0;
            sn_36S = 0.0;
        }

        // 15. ATMOSPHERIC OUTGASSING ao = atmospheric outgassing
        // --------------------------
        double ao_32S = ti_32S - (hg_32S + sq_32S) + sn_32S;
        double ao_33S = ti_33S - (hg_33S + sq_33S) + sn_33S;
        double ao_34S = ti_34S - (hg_34S + sq_34S) + sn_34S;
        double ao_36S = ti_36S - (hg_36S + sq_36S) + sn_36S;

        auto [ao_33R, ao_36R, ao_34R, ao_33d, ao_34d, ao_36d, ao_33D, ao_36D, ao_TS, logao_ST] =
            reservoir_isotope(ao_32S, ao_33S, ao_34S, ao_36S, VCDT, MIF_eq);

        // 16. RESERVOIRS FRACTIONATE
        // --------------------------
        // assuming space loss material and burial material have correct alpha between them
        double S32C, S33C, S34C, S36C;
        double S32E, S33E, S34E, S36E;
        // Set up constants for fractionation step
        FractionateConsts constants = {
            ao_TS, pu_F, bu_F,
            ao_32S, ao_33S, ao_34S, ao_36S,
            a33_eb, a34_eb, a36_eb
        };

        // Choose initial guess value for Newton-Raphson
        double S32_input;
        if (n == 1) {
            S32_input = ao_32S;
        } else {
            S32_input = S32C;
        }

        // Perform fractionation
        tie(S32C, S33C, S34C, S36C,
            S32E, S33E, S34E, S36E) =
            fractionate(constants, nr_step, nr_tol, S32_input);

        // comments from python version: 
        //  # fractionating space loss and photo-dissociation, and remainder going to frost... ONLY WORKS IF PD = 0
        // # pu_33R, pu_34R, pu_36R, pu_32S, pu_33S, pu_34S, pu_36S = process_R_S(pu_F,ao_33R,ao_34R,ao_36R,a33_pr,a34_pr,a36_pr)
        // # sn_32S, sn_33S, sn_34S, sn_36S = (ao_32S - pu_32S), (ao_33S - pu_33S), (ao_34S - pu_34S), (ao_36S - pu_36S) 
        // # S32E, S33E, S34E, S36E = pu_32S, pu_33S, pu_34S, pu_36S
        // # S32C, S33C, S34C, S36C = sn_32S, sn_33S, sn_34S, sn_36S
        // # E_33R, E_34R, E_36R, E_33d, E_34d, E_36d, E_33D, E_36D, E_ST, logE_ST = reservoir_isotope(S32E, S33E, S34E, S36E)
        // # C_33R, C_34R, C_36R, C_33d, C_34d, C_36d, C_33D, C_36D, C_ST, logC_ST  = reservoir_isotope(S32C, S33C, S34C, S36C)
        // #undo_32S, undo_33S, undo_34S, undo_36S = (S32E + S32C), (S33E + S33C), (S34E + S34C), (S36E + S36C)
        // #undo_ST,undo_34R, undo_33R, undo_36R, logundo_ST, undo_33d, undo_34d, undo_36d, undo_33D, undo_36D = reservoir_isotope(undo_32S, undo_33S, undo_34S, undo_36S)
        // #print(ao_33D, ao_36D,E_33D, E_36D, C_33D, C_36D,undo_33D, undo_36D)

        // ------------------------------
        // 17. Instantaneous isotope ratios of space reservoir
        // ------------------------------
        double iS_33R, iS_34R, iS_36R;
        double iS_33d, iS_34d, iS_36d;
        double iS_33D, iS_36D;
        double iS_TS, logiS_ST;
        tie(iS_33R, iS_34R, iS_36R,
            iS_33d, iS_34d, iS_36d,
            iS_33D, iS_36D,
            iS_TS, logiS_ST) =
            reservoir_isotope(S32E, S33E, S34E, S36E, VCDT, MIF_eq);

        // 18. plutonic and silicate/sulfate return
        // -----------------------------------------
        // pl = pluton
        auto [pl_33R, pl_34R, pl_36R,
            pl_32S, pl_33S, pl_34S, pl_36S] =
            process_R_S(pl_F, M_33R, M_34R, M_36R,
                        alpha.a33["pl"], alpha.a34["pl"], alpha.a36["pl"]);

        // rs = return silicate/sulfate
        auto [rs_33R, rs_34R, rs_36R,
            rs_32S, rs_33S, rs_34S, rs_36S] =
            process_R_S(rs_F, SS_33R, SS_34R, SS_36R,
                        alpha.a33["rm"], alpha.a34["rm"], alpha.a36["rm"]);

        // ------------------------------
        // 19. UPDATE RESERVOIRS (M_XYS, DM_XYS, F_XYS, SS_XYS, S_XYS)
        // ------------------------------
        // Mantle
        M_32S = mantle(M_32S, mo_32S, pl_32S, rs_32S, dm_32S);
        M_33S = mantle(M_33S, mo_33S, pl_33S, rs_33S, dm_33S);
        M_34S = mantle(M_34S, mo_34S, pl_34S, rs_34S, dm_34S);
        M_36S = mantle(M_36S, mo_36S, pl_36S, rs_36S, dm_36S);

        tie(M_33R, M_34R, M_36R,
            M_33d, M_34d, M_36d,
            M_33D, M_36D,
            M_ST.value, logM_ST) =
            reservoir_isotope(M_32S, M_33S, M_34S, M_36S,
                            VCDT, MIF_eq);

        // Deep Mantle
        DM_32S = deep_mantle(DM_32S, dm_32S);
        DM_33S = deep_mantle(DM_33S, dm_33S);
        DM_34S = deep_mantle(DM_34S, dm_34S);
        DM_36S = deep_mantle(DM_36S, dm_36S);

        tie(DM_33R, DM_34R, DM_36R,
            DM_33d, DM_34d, DM_36d,
            DM_33D, DM_36D,
            DM_ST.value, logDM_ST) =
            reservoir_isotope(DM_32S, DM_33S, DM_34S, DM_36S,
                            VCDT, MIF_eq);
        
        // cout << "Site A variables: " << endl;
        // cout << "----> fr_32S: " << fr_32S << endl;
        // cout << "----> fr_33S: " << fr_33S << endl;
        // cout << "----> fr_34S: " << fr_34S << endl;
        // cout << "----> fr_36S: " << fr_36S << endl;
        
        // cout << "---------> sn_32S: " << sn_32S << endl;
        // cout << "---------> sn_33S: " << sn_33S << endl;
        // cout << "---------> sn_34S: " << sn_34S << endl;
        // cout << "---------> sn_36S: " << sn_36S << endl;

        // Frost
        if (sn_F < 0.0) {
            fr_32S += sn_32S;
            fr_33S += sn_33S;
            fr_34S += sn_34S;
            fr_36S += sn_36S;
        }
        // cout << "Frost variables BEFORE: " << endl;
        // cout << "====>  F_32S: " <<  scientific << setprecision(16) << F_32S << endl;
        // cout << "====>  fr_32S: " <<  scientific << setprecision(16) << fr_32S << endl;
        // cout << "====>  S32C: " <<  scientific << setprecision(16) << S32C << endl;
        // cout << "====>  hg_32S: " <<  scientific << setprecision(16) << hg_32S << endl;
        // // cout << "====>  F_33S: " <<  F_33S << endl;
        // // cout << "====>  F_34S: " <<  F_34S << endl;
        // cout << "====>  F_36S: " <<  F_36S << endl;
        
        F_32S = frost(F_32S, fr_32S, S32C, hg_32S);
        F_33S = frost(F_33S, fr_33S, S33C, hg_33S);
        F_34S = frost(F_34S, fr_34S, S34C, hg_34S);
        F_36S = frost(F_36S, fr_36S, S36C, hg_36S); // ERROR: HAD BEEN USING F_36S with a fr_34S as an input! that's wrong
        
        // cout << "Frost variables AFTER: " << endl;
        // cout << "====>  F_32S: " <<  scientific << setprecision(16) << F_32S << endl;
        // cout << "====>  F_33S: " <<  scientific << setprecision(16) << F_33S << endl;
        // cout << "====>  F_34S: " <<  scientific << setprecision(16) << F_34S << endl;
        // cout << "====>  F_36S: " <<  scientific << setprecision(16) << F_36S << endl;

        // cout << "BEFORE F_ST.value: " << F_ST.value << endl;
        // cout << "The reservoir_isotope inputs are: " << endl;
        // cout << "-----> F_32S: " << scientific << setprecision(16) << F_32S << endl;
        // cout << "-----> F_33S: " << scientific << setprecision(16) << F_33S << endl;
        // cout << "-----> F_34S: " << scientific << setprecision(16) << F_34S << endl;
        // cout << "-----> F_36S: " << scientific << setprecision(16) << F_36S << endl;
        // cout << "-----> VCDT contents:" << endl;
        // for (const auto& [key, val] : VCDT) {
        //     cout << "        " << key << ": " << val << endl;
        // }
        // cout << "-----> MIF_eq: " << MIF_eq << endl;

        tie(F_33R, F_34R, F_36R,
            F_33d, F_34d, F_36d,
            F_33D, F_36D,
            F_ST.value, logF_ST) =
            reservoir_isotope(F_32S, F_33S, F_34S, F_36S,
                            VCDT, MIF_eq);
        
        // cout << "AFTER F_ST.value: " << F_ST.value << endl;

        // Silicate + Sulfate (SS)
        SS_32S = silsulf(SS_32S, rs_32S, sq_32S, pl_32S);
        SS_33S = silsulf(SS_33S, rs_33S, sq_33S, pl_33S);
        SS_34S = silsulf(SS_34S, rs_34S, sq_34S, pl_34S);
        SS_36S = silsulf(SS_36S, rs_36S, sq_36S, pl_36S);

        tie(SS_33R, SS_34R, SS_36R,
            SS_33d, SS_34d, SS_36d,
            SS_33D, SS_36D,
            SS_ST.value, logSS_ST) =
            reservoir_isotope(SS_32S, SS_33S, SS_34S, SS_36S,
                            VCDT, MIF_eq);

        // Space (S)
        S_32S = space(S_32S, S32E);
        S_33S = space(S_33S, S33E);
        S_34S = space(S_34S, S34E);
        S_36S = space(S_36S, S36E);

        tie(S_33R, S_34R, S_36R,
            S_33d, S_34d, S_36d,
            S_33D, S_36D,
            S_ST.value, logS_ST) =
            reservoir_isotope(S_32S, S_33S, S_34S, S_36S,
                            VCDT, MIF_eq);

        // ------------------------------
        // 20. MASS BALANCE CHECKS
        // ------------------------------

        // Total system mass and isotopes across all reservoirs
        double Io_ST  = reservoir_totals(M_ST.value,  DM_ST.value,  F_ST.value,  SS_ST.value,  S_ST.value);
        double Io_32S = reservoir_totals(M_32S,       DM_32S,       F_32S,       SS_32S,       S_32S);
        double Io_33S = reservoir_totals(M_33S,       DM_33S,       F_33S,       SS_33S,       S_33S);
        double Io_34S = reservoir_totals(M_34S,       DM_34S,       F_34S,       SS_34S,       S_34S);
        double Io_36S = reservoir_totals(M_36S,       DM_36S,       F_36S,       SS_36S,       S_36S);

        // Mass balance comparisons against initial totals
        mbST   = mass_balance(Io_ST_initial,  Io_ST);
        mb32S  = mass_balance(Io_32S_initial, Io_32S);
        mb33S  = mass_balance(Io_33S_initial, Io_33S);
        mb34S  = mass_balance(Io_34S_initial, Io_34S);
        mb36S  = mass_balance(Io_36S_initial, Io_36S);



        // ------------------------------
        // 21. STORE RESULTS IN STRUCT AND PUSH TO VECTOR 
        // ------------------------------
        // (print every 1000 steps or at end) this is a dataframe in the python version
        // ==========================================
        // === End of timestep: update full state ===
        SulfurState state;
        update_sulfur_state(state,
            static_cast<double>(n),              // t_step
            t_s_Myr,                             // t_step_Myr
            S_conc_M,                            // S_conc_M

            // δ values
            M_34d, F_34d, SS_34d, S_34d, DM_34d,
            M_33d, F_33d, SS_33d, S_33d, DM_33d,
            M_36d, F_36d, SS_36d, S_36d, DM_36d,

            // Δ values
            M_33D, F_33D, SS_33D, S_33D, DM_33D,
            M_36D, F_36D, SS_36D, S_36D, DM_36D,

            // ST values
            M_ST.value, F_ST.value, SS_ST.value, S_ST.value, DM_ST.value,

            // log(ST)
            logM_ST, logF_ST, logSS_ST, logS_ST, logDM_ST,

            // mass balance
            mbST, mb32S, mb33S, mb34S, mb36S,

            // Reference isotope ratios
            M_33R, F_33R, SS_33R, S_33R, DM_33R,
            M_34R, F_34R, SS_34R, S_34R, DM_34R,
            M_36R, F_36R, SS_36R, S_36R, DM_36R,

            // Escape-burial δ and Δ
            ao_34d, iS_34d,
            ao_33d, iS_33d,
            ao_36d, iS_36d,
            ao_33D, iS_33D,
            ao_36D, iS_36D,

            // Fluxes
            mo_F, fr_F, pr_F, pu_F, sn_F,
            pd_F, hg_F, bu_F, pl_F, rs_F,
            sq_F, ao_F, dm_F,

            // Ratios
            mo_mo_ca, bu_bu_pl,

            // Ref ratios (initial-state and escape-burial)
            ao_33R, iS_33R,
            ao_34R, iS_34R,
            ao_36R, iS_36R,

            // Escape-burial α values
            a33_bu, a34_bu, a36_bu,
            a33_eb, a34_eb, a36_eb
        );

        write_csv_row_padded(out, sulfur_state_to_python_like_row(state), header.size());

        if (n % 1000 == 0) out.flush();  // optional but nice


        // Push the updated state to store results in RAM. commment out if you don't want to store full results in RAM during calculation
        // results.push_back(state);

        // if (n % 1000 == 0 || n == static_cast<int>(end_time_Myr / t_step_Myr) - 1) {
        //     cout << "Last Step: " << n << ", Time = " << t_s_Myr << " Myr, F_34d = " << setprecision(16) << F_34d << endl;
        // }
        // if (n == total_steps) {
        //     cout << "Last step: " << n
        //         << " | Time = " << fixed << setprecision(1) << t_s_Myr << " Myr"
        //         << " | F_34d = " << setprecision(16) << F_34d << endl;
        // }
        if (n % 1000 == 0) {
            cout << "Step: " << n
                << " | Time = " << fixed << setprecision(1) << t_s_Myr << " Myr"
                << " | F_34d = " << setprecision(16) << F_34d << endl;
        }
        else if (n == total_steps) {
            cout << "Last step: " << n
                << " | Time = " << fixed << setprecision(1) << t_s_Myr << " Myr"
                << " | F_34d = " << setprecision(16) << F_34d << endl;

        }
        



    }
}
//     return; 
// }






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

