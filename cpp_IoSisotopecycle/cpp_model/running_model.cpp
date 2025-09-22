#include <iostream>
#include <iomanip>
#include <chrono>
#include <sstream>   // <-- needed for ostringstream
#include <cmath>     // <-- needed for round/fabs
#include "headers/isotope_evolution.hpp"
using namespace std;



// GETTING ALPHAS used to be done like this (below) but now it is in a header file that is then referenced in iso_evo()
// // getting alphas for a specific process and isotope
// // Example inputs
// string nofrac = "no";
// string MAF = "no";

// // Alpha 34
// double a34_gr = alphas("gr", 34, nofrac, MAF);
// double a34_th = alphas("th", 34, nofrac, MAF);
// double a34_pu = alphas("pu", 34, nofrac, MAF);
// double a34_ed = alphas("ed", 34, nofrac, MAF);
// double a34_ac = alphas("ac", 34, nofrac, MAF);
// double a34_ei = alphas("ei", 34, nofrac, MAF);
// double a34_rc = alphas("rc", 34, nofrac, MAF);
// double a34_ec = alphas("ec", 34, nofrac, MAF);
// double a34_pi = alphas("pi", 34, nofrac, MAF);
// double a34_mm = alphas("mm", 34, nofrac, MAF);
// double a34_dg = alphas("dg", 34, nofrac, MAF);
// double a34_xt = alphas("xt", 34, nofrac, MAF);
// double a34_fr = alphas("fr", 34, nofrac, MAF);
// double a34_pd = alphas("pd", 34, nofrac, MAF);
// double a34_hg = alphas("hg", 34, nofrac, MAF);
// double a34_rm = alphas("rm", 34, nofrac, MAF);
// double a34_vp = alphas("vp", 34, nofrac, MAF);
// double a34_sq = alphas("sq", 34, nofrac, MAF);
// double a34_dm = alphas("dm", 34, nofrac, MAF);
// double a34_cf = alphas("cf", 34, nofrac, MAF);
// double a34_mo = a34_mm * a34_dg;
// double a34_pl = a34_mm * a34_xt;

// // Alpha 33
// double a33_gr = alphas("gr", 33, nofrac, MAF);
// double a33_th = alphas("th", 33, nofrac, MAF);
// double a33_pu = alphas("pu", 33, nofrac, MAF);
// double a33_ed = alphas("ed", 33, nofrac, MAF);
// double a33_ac = alphas("ac", 33, nofrac, MAF);
// double a33_ei = alphas("ei", 33, nofrac, MAF);
// double a33_rc = alphas("rc", 33, nofrac, MAF);
// double a33_ec = alphas("ec", 33, nofrac, MAF);
// double a33_pi = alphas("pi", 33, nofrac, MAF);
// double a33_mm = alphas("mm", 33, nofrac, MAF);
// double a33_dg = alphas("dg", 33, nofrac, MAF);
// double a33_xt = alphas("xt", 33, nofrac, MAF);
// double a33_fr = alphas("fr", 33, nofrac, MAF);
// double a33_pd = alphas("pd", 33, nofrac, MAF);
// double a33_hg = alphas("hg", 33, nofrac, MAF);
// double a33_rm = alphas("rm", 33, nofrac, MAF);
// double a33_vp = alphas("vp", 33, nofrac, MAF);
// double a33_sq = alphas("sq", 33, nofrac, MAF);
// double a33_dm = alphas("dm", 33, nofrac, MAF);
// double a33_cf = alphas("cf", 33, nofrac, MAF);
// double a33_mo = a33_mm * a33_dg;
// double a33_pl = a33_mm * a33_xt;

// // Alpha 36
// double a36_gr = alphas("gr", 36, nofrac, MAF);
// double a36_th = alphas("th", 36, nofrac, MAF);
// double a36_pu = alphas("pu", 36, nofrac, MAF);
// double a36_ed = alphas("ed", 36, nofrac, MAF);
// double a36_ac = alphas("ac", 36, nofrac, MAF);
// double a36_ei = alphas("ei", 36, nofrac, MAF);
// double a36_rc = alphas("rc", 36, nofrac, MAF);
// double a36_ec = alphas("ec", 36, nofrac, MAF);
// double a36_pi = alphas("pi", 36, nofrac, MAF);
// double a36_mm = alphas("mm", 36, nofrac, MAF);
// double a36_dg = alphas("dg", 36, nofrac, MAF);
// double a36_xt = alphas("xt", 36, nofrac, MAF);
// double a36_fr = alphas("fr", 36, nofrac, MAF);
// double a36_pd = alphas("pd", 36, nofrac, MAF);
// double a36_hg = alphas("hg", 36, nofrac, MAF);
// double a36_rm = alphas("rm", 36, nofrac, MAF);
// double a36_vp = alphas("vp", 36, nofrac, MAF);
// double a36_sq = alphas("sq", 36, nofrac, MAF);
// double a36_dm = alphas("dm", 36, nofrac, MAF);
// double a36_cf = alphas("cf", 36, nofrac, MAF);
// double a36_mo = a36_mm * a36_dg;
// double a36_pl = a36_mm * a36_xt;




int main() {
    auto start = chrono::high_resolution_clock::now();

    auto fmt = [](double v) {
    std::ostringstream oss;
    // choose 6 vs 17 if you want; here I’ll just use 17 then strip
    oss << std::fixed << std::setprecision(17) << v;

    std::string s = oss.str();
    // strip trailing zeros and any trailing '.'
    if (s.find('.') != std::string::npos) {
        while (!s.empty() && s.back() == '0') s.pop_back();
        if (!s.empty() && s.back() == '.') s.pop_back();
    }
    // ensure .0 for integers
    if (s.find('.') == std::string::npos) s += ".0";
    return s;
};
//     // ==============================================================
//     // === TESTING INDIVIDUAL FUNCTIONS SECTION === //
//     // ==============================================================

//     // testing function-by-function comparison between cpp and python
//     cout << "Testing isotope evolution functions..." << endl;
//     cout << fixed << setprecision(17);
//     // testing a_3X()
//     cout << "Testing a_3X(): " << endl;
//     double a34 = 0.9998;
//     int X = 33;
//     double a = a_3X(a34, X);
//     cout << "\ta_" << X << " = " << fmt(a) << endl;

//     // testing alphas()
//     cout << "Testing alphas(): " << endl;
//     string process = "gr";
//     string nofrac_1 = "no";
//     string MAF_1 = "no";
//     double result = alphas(process, X, nofrac_1, MAF_1);
//     cout << "\talpha_" << X << " = " << fmt(result) << endl;

//     // --- selected alphas (defined before printing) ---
//     cout << "Testing selected alphas()... " << endl;
//     double a34_ec = alphas("ec", 34, nofrac_1, MAF_1);
//     double a33_pi = alphas("pi", 33, nofrac_1, MAF_1);
//     double a36_fr = alphas("fr", 36, nofrac_1, MAF_1);
//     cout << "\ta34_ec = " << fmt(a34_ec) << endl;
//     cout << "\ta33_pi = " << fmt(a33_pi) << endl;
//     cout << "\ta36_fr = " << fmt(a36_fr) << endl;

//     // --- rate_kgs ---
//     cout << "=== Test rate_kgs ===" << endl;
//     string proc = "pi";
//     if (rate_kgs.count(proc)) {
//         RateBounds r = rate_kgs.at(proc);
//         cout << "\t" << proc << ": min = " << fmt(r.min) << ", max = " << fmt(r.max) << endl;
//     }
    
//     // --- VCDT ---
//     cout << "\n=== Test VCDT ===" << endl;
//     for (const auto& kv : VCDT) {
//         cout << "\tdelta " << kv.first << " ratio = " << fmt(kv.second) << endl;
//     }

//     // --- R2d ---
//     double R = 0.045; // example ratio
//     int n = 34;
//     double delta = R2d(R, n, VCDT);
//     cout << "\n=== Test R2d() ===" << endl;
//     // cout << "R: " << R << endl;
//     cout << "\tInput ratio R = " << fmt(R) << endl;
//     cout << "\tIsotope n = " << n << endl;
//     cout << "\tDelta_" << n << " = " << fmt(delta) << endl;

//     // --- d2R ---
//     cout << "=== Test d2R ===" << endl;
//     double d = 10.0;  // +10‰
//     double Rback = d2R(d, n, VCDT);
//     cout << "\td = " << fmt(d) << " ‰, n = " << n << " -> R = " << fmt(Rback) << endl;

//     // --- ROT ---
//     cout << "\n=== Test ROT ===" << endl;
//     // example ratios (33/32 = a, 32/32=c, 34/32=b, 36/32=d)
//     a = 0.0079;
//     double c = 1.0;
//     double b = 0.0442;
//     double dd = 0.000153;
//     double Xrot = ROT(a, c, b, dd);
//     cout << "\tX (33S/ST) = " << fmt(Xrot) << "\n";

//     // --- kgs2mSs ---
//     cout << "\n=== Test kgs2mSs ===" << endl;
//     double kgs = 205.0; // kg/s SO2
//     double molSps = kgs2mSs(kgs);
//     cout << "\t" << fmt(kgs) << " kg/s SO2 -> " << fmt(molSps) << " mol S/s\n";

//     // --- MIF tests ---
//     cout << "\n=== Test MIF ===\n";
//     string MIF_eq_1 = "FW2003_pg3";

//     // 1) Sanity: using standard ratios should give ~0,0
//     double R33s = VCDT.at("33"), R34s = VCDT.at("34"), R36s = VCDT.at("36");
//     auto Dstd = MIF(R33s, R34s, R36s, VCDT, MIF_eq_1);
//     cout << "\tAt standards: D33 = " << fmt(Dstd.first) << ", D36 = " << fmt(Dstd.second) << "\n";

//     // 2) Non-trivial: construct ratios from chosen deltas
//     double d34 = 10.0, d33 = 5.0, d36 = 20.0;
//     double R34 = d2R(d34, 34, VCDT);
//     double R33 = d2R(d33, 33, VCDT);
//     double R36 = d2R(d36, 36, VCDT);
//     auto D = MIF(R33, R34, R36, VCDT, MIF_eq_1);
//     cout << "\tFor d34=" << fmt(d34) << "‰, d33=" << fmt(d33) << "‰, d36=" << fmt(d36) << "‰  -> D33 = " << fmt(D.first) << ", D36 = " << fmt(D.second) << "\n";

//     // --- calc_rate tests ---
//     cout << "\n=== Test calc_rate ===\n";
//     auto pd = rate_kgs.at("pd");
//     cout << "\tf=\"N\" -> " << fmt(calc_rate("N",   pd.min, pd.max)) << " kg/s\n";
//     cout << "\tf=0.5 -> " << fmt(calc_rate(0.5,   pd.min, pd.max)) << " kg/s\n";
//     cout << "\tf=\"0.5\" -> " << fmt(calc_rate("0.5", pd.min, pd.max)) << " kg/s\n";
//     cout << "\tf=2 -> " << fmt(calc_rate(2,     pd.min, pd.max)) << " kg/s\n";

//     // === Test initial_mass ===
//     double min_mass = 2.95e21, max_mass = 2.24e22;
//     cout << "=== Test initial_mass ===\n";
//     cout << "fraction=0.0 -> " << initial_mass(min_mass, max_mass, 0.0) << "\n";
//     cout << "fraction=1.0 -> " << initial_mass(min_mass, max_mass, 1.0) << "\n";
//     cout << "fraction=0.5 -> " << initial_mass(min_mass, max_mass, 0.5) << "\n";

//     // === Test reservoir_totals ===
//     double M = 1.0, DM = 2.0, F = 3.0, SS = 4.0, S = 5.0;
//     cout << "\n=== Test reservoir_totals ===\n";
//     cout << "Sum = " << reservoir_totals(M, DM, F, SS, S) << "\n";

//     // === Test mass_balance ===
//     double init = 100.0, now = 80.0;
//     cout << "\n=== Test mass_balance ===\n";
//     cout << "Mass balance fraction = " << mass_balance(init, now) << "\n";

//     // === Test S_conc ===
//     double ST = 1.0e6; // moles S
//     double mass = 5.0e9; // kg
//     cout << "\n=== Test S_conc ===\n";
//     cout << "Concentration = " << S_conc(ST, mass) << " (fraction by mass)\n";

//     // === Test molSs_mm ===
//     double mag_kgs = 200.0, mag_S = 0.01, f = 0.5;
//     cout << "\n=== Test molSs_mm ===\n";
//     cout << "mol S/s = " << molSs_mm(mag_kgs, mag_S, f) << "\n";

//     // === Test reservoir_R_S ===
//     cout << "\n=== Test reservoir_R_S ===\n";

//     // sample inputs
//     double test_S32 = 100.0, test_S33 = 1.0, test_S34 = 4.4, test_S36 = 0.15;

//     // unpack the tuple
//     auto [test_R33, test_R34, test_R36, test_T32, test_T33, test_T34, test_T36, test_ST] = reservoir_R_S(test_S32, test_S33, test_S34, test_S36);

//     // print
//     cout << "Inputs: S32=" << test_S32 << ", S33=" << test_S33 << ", S34=" << test_S34 << ", S36=" << test_S36 << "\n";
//     cout << "R33 = " << test_R33 << ", R34 = " << test_R34 << ", R36 = " << test_R36 << "\n";
//     cout << "T32 = " << test_T32 << ", T33 = " << test_T33 << ", T34 = " << test_T34 << ", T36 = " << test_T36 << "\n";
//     cout << "ST  = " << test_ST  << "\n";

//     // edge case: S32 == 0
//     auto [r33z, r34z, r36z, t32z, t33z, t34z, t36z, stz] = reservoir_R_S(0.0, 1.0, 1.0, 1.0);
//     cout << "\nS32 == 0 edge case:\n";
//     cout << "R33=" << r33z << ", R34=" << r34z << ", R36=" << r36z
//         << ", T32=" << t32z << ", T33=" << t33z << ", T34=" << t34z
//         << ", T36=" << t36z << ", ST=" << stz << "\n";


//     cout << "\n=== Test reservoir_isotope ===\n";
//     {
//     // sample reservoir
//     double S32 = 100.0, S33 = 1.0, S34 = 4.4, S36 = 0.15;
//     string MIF_eq_2 = "FW2003_pg3";

//     auto [R33, R34, R36, d33, d34, d36, D33, D36, ST, log_ST] =
//         reservoir_isotope(S32, S33, S34, S36, VCDT, MIF_eq_2);

//     cout << "Inputs: S32=" << S32 << ", S33=" << S33
//          << ", S34=" << S34 << ", S36=" << S36 << "\n";
//     cout << "R33=" << fmt(R33) << ", R34=" << fmt(R34) << ", R36=" << fmt(R36) << "\n";
//     cout << "d33=" << fmt(d33) << "‰, d34=" << fmt(d34) << "‰, d36=" << fmt(d36) << "‰\n";
//     cout << "D33=" << fmt(D33) << "‰, D36=" << fmt(D36) << "‰\n";
//     cout << "ST="  << fmt(ST)  << ", log_ST=" << fmt(log_ST) << "\n";
//     }

//     // edge case: S32 == 0
//     {
//         auto [R33z, R34z, R36z, d33z, d34z, d36z, D33z, D36z, STz, log_STz] =
//             reservoir_isotope(0.0, 1.0, 1.0, 1.0, VCDT, string("FW2003_pg3"));
//         cout << "\nS32 == 0 edge case:\n";
//         cout << "R33=" << R33z << ", R34=" << R34z << ", R36=" << R36z
//             << ", d33=" << d33z << ", d34=" << d34z << ", d36=" << d36z
//             << ", D33=" << D33z << ", D36=" << D36z
//             << ", ST=" << STz << ", log_ST=" << log_STz << "\n";
//     }

//     // === Test initial_reservoir ===
// cout << "\n=== Test initial_reservoir ===\n";
// {
//     // inputs
//     double ST_in = 1.0e9;      // total sulfur (mol)
//     double d33_in = 0.0;       // per mil
//     double d34_in = 10.0;      // per mil
//     double d36_in = 20.0;      // per mil
//     string MIF_eq_3 = "FW2003_pg3";

//     auto [log_ST_ir, R33_ir, R34_ir, R36_ir,
//           T32_ir, T33_ir, T34_ir, T36_ir,
//           D33_ir, D36_ir,
//           S32_ir, S33_ir, S34_ir, S36_ir] =
//         initial_reservoir(ST_in, d33_in, d34_in, d36_in, VCDT, MIF_eq_3);

//     cout << "Inputs: ST=" << fmt(ST_in)
//          << ", d33=" << fmt(d33_in) << "‰"
//          << ", d34=" << fmt(d34_in) << "‰"
//          << ", d36=" << fmt(d36_in) << "‰\n";

//     cout << "R33=" << fmt(R33_ir) << ", R34=" << fmt(R34_ir) << ", R36=" << fmt(R36_ir) << "\n";
//     cout << "T32=" << fmt(T32_ir) << ", T33=" << fmt(T33_ir)
//          << ", T34=" << fmt(T34_ir) << ", T36=" << fmt(T36_ir) << "\n";
//     cout << "D33=" << fmt(D33_ir) << "‰, D36=" << fmt(D36_ir) << "‰\n";
//     cout << "S32=" << fmt(S32_ir) << ", S33=" << fmt(S33_ir)
//          << ", S34=" << fmt(S34_ir) << ", S36=" << fmt(S36_ir) << "\n";
//     cout << "log_ST=" << fmt(log_ST_ir) << "\n";
// }

// // edge case: ST <= 0
// {
//     auto [log_ST0, R33_0, R34_0, R36_0,
//           T32_0, T33_0, T34_0, T36_0,
//           D33_0, D36_0,
//           S32_0, S33_0, S34_0, S36_0] =
//         initial_reservoir(0.0, 0.0, 0.0, 0.0, VCDT, string("FW2003_pg3"));

//     cout << "\nST <= 0 edge case:\n";
//     cout << "R33=" << R33_0 << ", R34=" << R34_0 << ", R36=" << R36_0
//          << ", T32=" << T32_0 << ", T33=" << T33_0
//          << ", T34=" << T34_0 << ", T36=" << T36_0 << "\n";
//     cout << "D33=" << D33_0 << ", D36=" << D36_0
//          << ", S32=" << S32_0 << ", S33=" << S33_0
//          << ", S34=" << S34_0 << ", S36=" << S36_0
//          << ", log_ST=" << log_ST0 << "\n";
// }


// cout << "\n=== Test process_R_S ===" << endl;

// double a_F = 1.5;
// double a_res_33R = 0.0079;
// double a_res_34R = 0.0442;
// double a_res_36R = 0.000153;
// double a_a33 = 0.999;
// double a_a34 = 0.998;
// double a_a36 = 1.001;

// // Unpack tuple directly
// auto [a_R33, a_R34, a_R36, a_S32, a_S33, a_S34, a_S36] = process_R_S(a_F, a_res_33R, a_res_34R, a_res_36R, a_a33, a_a34, a_a36);

// cout << "Inputs:\n";
// cout << "  F=" << a_F
//      << ", res_33R=" << a_res_33R
//      << ", res_34R=" << a_res_34R
//      << ", res_36R=" << a_res_36R << "\n";
// cout << "  a33=" << a_a33 << ", a34=" << a_a34 << ", a36=" << a_a36 << "\n";

// cout << "Outputs:\n";
// cout << "  R33=" << a_R33
//      << ", R34=" << a_R34
//      << ", R36=" << a_R36 << "\n";
// cout << "  S32=" << a_S32
//      << ", S33=" << a_S33
//      << ", S34=" << a_S34
//      << ", S36=" << a_S36 << "\n";

// // Edge case: F = 0.0
// double a_F_zero = 0.0;
// auto [a_ZR33, a_ZR34, a_ZR36, a_ZS32, a_ZS33, a_ZS34, a_ZS36] =
//     process_R_S(a_F_zero, a_res_33R, a_res_34R, a_res_36R, a_a33, a_a34, a_a36);

// cout << "\nEdge case: F=0.0\n";
// cout << "  Output: (" << a_ZR33 << ", " << a_ZR34 << ", " << a_ZR36
//      << ", " << a_ZS32 << ", " << a_ZS33 << ", " << a_ZS34 << ", " << a_ZS36 << ")\n";

// // test mantle, deep_mantle, frost, silsulf, space
// cout << "\n=== Test mantle, deep_mantle, frost, silsulf, space ===\n";

// double b_M = 1.0e22; // kg
// double b_mo = 1.0e21; // kg
// double b_pl = 5.0e20; // kg
// double b_rs = 2.0e20; // kg
// double b_dm = 3.0e20; // kg
// cout << "mantle: " << mantle(b_M, b_mo, b_pl, b_rs, b_dm) << endl;
// double b_DM = 4.0e20; // kg
// cout << " deep mantle: " << deep_mantle(b_DM, b_dm)  << endl;
// double b_F = 1.0e20; // kg
// double b_fr = 2.0e19; // kg
// double b_bu = 1.0e19; // kg
// double b_hg = 5.0e18; // kg
// cout << " frost: " << frost(b_F, b_fr, b_bu, b_hg) << endl;
// double b_SS = 6.0e19; // kg
// double b_sq = 3.0e19; // kg
// cout << " silsulf: " << silsulf(b_SS, b_rs, b_sq, b_pl) << endl;
// double b_S = 7.0e19; // kg
// double b_pu = 2.0e19; // kg
// cout << " space: " << space(b_S, b_pu)  << endl;



// std::cout << "\n=== Test a_bu and a_pr ===\n";

// // ---- a_bu test ----
// double sn_F = 2.0e19;
// double pd_F = 3.0e19;
// double a_vp = 0.995;
// double a_pd = 1.002;

// double a_bu_val = a_bu(sn_F, pd_F, a_vp, a_pd);
// std::cout << "a_bu: " << a_bu_val << "\n";

// // ---- a_pr test ----
// std::unordered_map<std::string, double> alphas = {
//     {"pu", 0.98}, {"gr", 1.01}, {"th", 0.97}, {"ed", 1.003},
//     {"ac", 0.999}, {"ei", 1.002}, {"rc", 1.001}, {"ec", 0.998},
//     {"pi", 1.004}
// };
// std::unordered_map<std::string, double> rates = {
//     {"pr", 10.0}, {"ed", 2.0}, {"ac", 1.5}, {"ei", 3.0},
//     {"rc", 1.0}, {"ec", 0.5}, {"pi", 2.0}
// };

// double a_pr_val = a_pr(alphas, rates);
// std::cout << "a_pr: " << a_pr_val << "\n";


// // === Test newton_raphson_like_python ===
// cout << "\n=== Test newton_raphson_like_python ===" << endl;
// {
//     // Solve x^2 - a = 0 for a = 2 (root = sqrt(2))
//     struct Params { double a; } constants{2.0};

//     auto eqs   = [](double x, const Params& c) noexcept { return x*x - c.a; };
//     auto deriv = [](double x, const Params&)   noexcept { return 2.0*x;     };

//     const double x0   = 1.0;     // initial guess
//     const double tol  = 1e-12;   // tolerance
//     const double step = 1.0;     // standard Newton step

//     using clock = chrono::steady_clock;
//     const auto t0 = clock::now();

//     // Quiet run (no per-iteration output). For CSV per-iter, pass &cout instead of nullptr.
//     double root = newton_raphson(x0, constants, tol, step, eqs, deriv,
//                                  /*maxiter*/1000, /*per_iter_out*/nullptr, /*print_header*/true);
//     //if you want per-iteration output do this: double root = newton_raphson(x0, constants, tol, step, eqs, deriv, 1000, &cout, true);
//     const auto t1 = clock::now();
//     const chrono::duration<double> elapsed = t1 - t0;

//     cout << "root = " << fmt(root) << '\n';
//     cout << "elapsed_seconds = " << fmt(elapsed.count()) << '\n';

    

// }


// cout << "\n=== Test fractionate (C++ port) ===\n";

// // Example numbers (use your real case values here)
// FractionateConsts fc{
//     /*STT*/ 0.0, /*STE*/ 0.0, /*STC*/ 10.0,   // only STC is used inside f/df
//     /*S32T*/ 100.0, /*S33T*/ 1.0, /*S34T*/ 4.4, /*S36T*/ 0.15,
//     /*a33*/ 0.999, /*a34*/ 0.998, /*a36*/ 1.001
// };

// double nr_step_1 = 1.0;
// double nr_tol_1  = 1e-12;
// double guessx  = 50.0;

// using clock = chrono::steady_clock;
// auto t0 = clock::now();

// // Quiet run: pass nullptr; for per-iteration CSV, pass &cout and set print_header=true
// auto [S32C, S33C, S34C, S36C, S32E, S33E, S34E, S36E] =
//     fractionate(fc, nr_step_1, nr_tol_1, guessx, /*per_iter_out*/nullptr, /*print_header*/true);
//     // if you want per iteration output do this:
//     //fractionate(fc, nr_step, nr_tol, guessx, &cout, true);


// auto t1 = clock::now();
// chrono::duration<double> frac_elapsed = t1 - t0;

// cout << "Crust:\n";
// cout << "  S32C=" << fmt(S32C) << ", S33C=" << fmt(S33C)
//      << ", S34C=" << fmt(S34C) << ", S36C=" << fmt(S36C) << "\n";
// cout << "Ejecta:\n";
// cout << "  S32E=" << fmt(S32E) << ", S33E=" << fmt(S33E)
//      << ", S34E=" << fmt(S34E) << ", S36E=" << fmt(S36E) << "\n";
// cout << "elapsed_seconds = " << fmt(frac_elapsed.count()) << "\n";


// // =============
// // END TESTING OF INDIVIDUAL FUNCTIONS
// // =============


auto print_inputs = [](double t_step_Myr, double end_time_Myr, double nr_step, double nr_tol,
                      const MaybeDouble& DM_mass_f, const MaybeDouble& DM_mass_f_sil, const MaybeDouble& DM_ST,
                      double DM_33d, double DM_34d, double DM_36d,
                      const MaybeDouble& M_mass_f, const MaybeDouble& M_mass_f_sil, const MaybeDouble& M_ST,
                      double M_33d, double M_34d, double M_36d,
                      const MaybeDouble& F_mass_f, const MaybeDouble& F_ST,
                      double F_33d, double F_34d, double F_36d,
                      const MaybeDouble& SS_mass_f, const MaybeDouble& SS_ST,
                      double SS_33d, double SS_34d, double SS_36d,
                      const MaybeDouble& S_ST, double S_33d, double S_34d, double S_36d,
                      const unordered_map<string,double>& rate_f,
                      const string& oscillate, double resurf_cm_yr, double sil_mag_S, double thick_C,
                      double f_S2, double f_pl2mo, double f_sq, double f_deep, double f_remobilised, double f_pu) {
    cout << "\n=== C++ MODEL INPUTS ===" << endl;
    cout << "t_step_Myr=" << t_step_Myr << ", end_time_Myr=" << end_time_Myr
         << ", nr_step=" << nr_step << ", nr_tol=" << nr_tol << endl;
    cout << "DM_mass_f=(" << DM_mass_f.flag << "," << DM_mass_f.value << "), "
         << "DM_mass_f_sil=(" << DM_mass_f_sil.flag << "," << DM_mass_f_sil.value << "), "
         << "DM_ST=(" << DM_ST.flag << "," << DM_ST.value << "), "
         << "DM_33d=" << DM_33d << ", DM_34d=" << DM_34d << ", DM_36d=" << DM_36d << endl;
    cout << "M_mass_f=(" << M_mass_f.flag << "," << M_mass_f.value << "), "
         << "M_mass_f_sil=(" << M_mass_f_sil.flag << "," << M_mass_f_sil.value << "), "
         << "M_ST=(" << M_ST.flag << "," << M_ST.value << "), "
         << "M_33d=" << M_33d << ", M_34d=" << M_34d << ", M_36d=" << M_36d << endl;
    cout << "F_mass_f=(" << F_mass_f.flag << "," << F_mass_f.value << "), "
         << "F_ST=(" << F_ST.flag << "," << F_ST.value << "), "
         << "F_33d=" << F_33d << ", F_34d=" << F_34d << ", F_36d=" << F_36d << endl;
    cout << "SS_mass_f=(" << SS_mass_f.flag << "," << SS_mass_f.value << "), "
         << "SS_ST=(" << SS_ST.flag << "," << SS_ST.value << "), "
         << "SS_33d=" << SS_33d << ", SS_34d=" << SS_34d << ", SS_36d=" << SS_36d << endl;
    cout << "S_ST=(" << S_ST.flag << "," << S_ST.value << "), "
         << "S_33d=" << S_33d << ", S_34d=" << S_34d << ", S_36d=" << S_36d << endl;

    cout << "rate_f: ";
    for (const auto& kv : rate_f) cout << kv.first << "=" << kv.second << " ";
    cout << endl;

    cout << "oscillate=" << oscillate
         << ", resurf_cm_yr=" << resurf_cm_yr
         << ", sil_mag_S=" << sil_mag_S
         << ", thick_C=" << thick_C << endl;
    cout << "f_S2=" << f_S2 << ", f_pl2mo=" << f_pl2mo
         << ", f_sq=" << f_sq << ", f_deep=" << f_deep
         << ", f_remobilised=" << f_remobilised << ", f_pu=" << f_pu << endl;
    cout << "===========================\n" << endl;
};



/////////////////////////////////////////////////////////////////
// iso_evo() function test
/////////////////////////////////////////////////////////////////

// Basic example inputs (replace with real values as needed)
double t_step_Myr = 1e-1; // timestep in Myr
double end_time_Myr = 200.1; // 4570 // end time in Myr
double nr_step = 1.;
double nr_tol = 1e5; //THIS SEEMS LIKE A WRONG VALUE BUT IT IS WHAT IS IN THE PYTHON CODE SO I NEED IT LIKE THIS TO MATCH THE CODES

// === RESERVOIRS === //
// MaybeDoubles turning on and off reservoirs: 
// deep mantle variables
MaybeDouble DM_mass_f{"N", 0.0}; // {"Y", value} to turn on, {"N", 0.0} to turn off
MaybeDouble DM_mass_f_sil{"N", 0.0};
MaybeDouble DM_ST{"N", 0.0};
double DM_33d = 0.0; double DM_34d = 0.0; double DM_36d = 0.0;

// mantle variables
MaybeDouble M_mass_f{"Y", 0.1257}; // mantle mass S fractionation of min and max [0.1257]
MaybeDouble M_mass_f_sil{"Y", M_mass_f.value};
MaybeDouble M_ST{"N", 0.0};
double M_33d = 0.0; double M_34d = 0.0; double M_36d = 0.0;

// crustal frost variables
MaybeDouble F_mass_f{"N", 0.0};
MaybeDouble F_ST{"N", 0.0};
double F_33d = 0.0; double F_34d = 0.0; double F_36d = 0.0;

// crustal sulfate variables
MaybeDouble SS_mass_f{"N", 0.0};
MaybeDouble SS_ST{"N", 0.0};
double SS_33d = 0.0; double SS_34d = 0.0; double SS_36d = 0.0;

MaybeDouble S_ST{"N", 0.0};
double S_33d = 0.0; double S_34d = 0.0; double S_36d= 0.0;

// rate_f: 
unordered_map<string, double> rate_f = {
    {"pd", 1.0}, {"pi", 1.0}, {"ei", 1.0},
    {"ed", 1.0}, {"ac", 1.0}, {"rc", 1.0}, {"ec", 1.0}//,
    // {"pu", 1.0}
};

// options for calculating cap-delta values
// "O2017_eq5" = equation (5) in Ono (2017): D3XS = d3XS - e*d34S 
// "O2017_eq6" = equation (6) in Ono (2017):  D3X = d3X - ((d34+1.)^e3X - 1.)
// "FW2003_pg3" = equations on pg. 3 of Farquar & Wing (2003): D3XS = d3XS - 1000*((1+d34S/1000)^e - 1)
string MIF_eq = "FW2003_pg3";
string nofrac = "no";
string MAF = "no";

//options
string oscillate = "N";
double resurf_cm_yr = 1.0; // cm/yr (today rate = 1.0 cm/yr)
double sil_mag_S = 0.001;   // fraction of S in mantle melt
double thick_C = 40000.0;  // crustal thickness m - min = 30000, max = 50000
double f_S2 = 0.2;       // fraction of S for SO2 due to homogeneous gas equilibria [0.2]
double f_pl2mo = 0.8; // fraction of mantle melting that goes to plutons [0.8]
double f_sq = 0.25; // fraction of sulfate sequestration from gas [0.25]
double f_deep = 0.0; // fraction of mantle melting rate that is added from the deep mantle if present [0.]
double f_remobilised = 1.0; // fraction of f_crust_returned that is remobilised [1.]
double f_pu = 0.5; // fraction of material that could be lost to space that is [0.5]


print_inputs(t_step_Myr, end_time_Myr, nr_step, nr_tol,
             DM_mass_f, DM_mass_f_sil, DM_ST, DM_33d, DM_34d, DM_36d,
             M_mass_f, M_mass_f_sil, M_ST, M_33d, M_34d, M_36d,
             F_mass_f, F_ST, F_33d, F_34d, F_36d,
             SS_mass_f, SS_ST, SS_33d, SS_34d, SS_36d,
             S_ST, S_33d, S_34d, S_36d,
             rate_f, oscillate, resurf_cm_yr, sil_mag_S, thick_C,
             f_S2, f_pl2mo, f_sq, f_deep, f_remobilised, f_pu);




// Call the main function
iso_evo(
    t_step_Myr, end_time_Myr, nr_step, nr_tol,   // time and Newton-Raphson params

    DM_mass_f, DM_mass_f_sil, DM_ST, DM_33d, DM_34d, DM_36d,  // DM
    M_mass_f, M_mass_f_sil, M_ST, M_33d, M_34d, M_36d,  // M
    F_mass_f, F_ST, F_33d, F_34d, F_36d,       // F
    SS_mass_f, SS_ST, SS_33d, SS_34d, SS_36d,       // SS
    S_ST, S_33d, S_34d, S_36d,            // S

    rate_f, oscillate, resurf_cm_yr, sil_mag_S, thick_C, //other parameters (rates, resurfacing, silicate melt S content, crustal thickness)
    
    f_S2, f_pl2mo, f_sq, f_deep, f_remobilised, f_pu, // fractions (S2/SO2, plutonic/magmatic, sulfate sequestration, deep mantle, remobilised crust, space loss

    nofrac, MAF, MIF_eq // no fractionation, MAF, MIF equation options
);





/////////////////////////////////////////////////////////////////
// Execution time measurement
/////////////////////////////////////////////////////////////////
auto end = chrono::high_resolution_clock::now();
chrono::duration<double> elapsed = end - start;
cout << "\nExecution time: " << elapsed.count() << " seconds" << endl;

return 0;
}

