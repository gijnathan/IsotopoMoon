#include <iostream>
#include <iomanip>
#include <chrono>
#include <sstream>
#include <cmath>

#include "headers/isotope_evolution.hpp"
#include "headers/config.hpp"
#include "headers/config_io.hpp"

using namespace std;

int main(int argc, char** argv) {
    auto start = chrono::high_resolution_clock::now();

    auto fmt = [](double v) {
        ostringstream oss;
        oss << fixed << setprecision(17) << v;
        string s = oss.str();
        if (s.find('.') != string::npos) {
            while (!s.empty() && s.back() == '0') s.pop_back();
            if (!s.empty() && s.back() == '.') s.pop_back();
        }
        if (s.find('.') == string::npos) s += ".0";
        return s;
    };

    // DEFINE print_inputs BEFORE YOU USE IT
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
        cout << "\n=== C++ MODEL INPUTS ===\n";
        cout << "t_step_Myr=" << t_step_Myr << ", end_time_Myr=" << end_time_Myr
             << ", nr_step=" << nr_step << ", nr_tol=" << nr_tol << "\n";
        cout << "===========================\n\n";
    };

    // LOAD CONFIG
    string cfg_path = (argc >= 2) ? argv[1] : "config.txt";
    ModelConfig cfg = load_config(cfg_path);

    print_inputs(cfg.t_step_Myr, cfg.end_time_Myr, cfg.nr_step, cfg.nr_tol,
                 cfg.DM_mass_f, cfg.DM_mass_f_sil, cfg.DM_ST, cfg.DM_33d, cfg.DM_34d, cfg.DM_36d,
                 cfg.M_mass_f, cfg.M_mass_f_sil, cfg.M_ST, cfg.M_33d, cfg.M_34d, cfg.M_36d,
                 cfg.F_mass_f, cfg.F_ST, cfg.F_33d, cfg.F_34d, cfg.F_36d,
                 cfg.SS_mass_f, cfg.SS_ST, cfg.SS_33d, cfg.SS_34d, cfg.SS_36d,
                 cfg.S_ST, cfg.S_33d, cfg.S_34d, cfg.S_36d,
                 cfg.rate_f, cfg.oscillate, cfg.resurf_cm_yr, cfg.sil_mag_S, cfg.thick_C,
                 cfg.f_S2, cfg.f_pl2mo, cfg.f_sq, cfg.f_deep, cfg.f_remobilised, cfg.f_pu);

    iso_evo(
        cfg.t_step_Myr, cfg.end_time_Myr, cfg.nr_step, cfg.nr_tol,
        cfg.DM_mass_f, cfg.DM_mass_f_sil, cfg.DM_ST, cfg.DM_33d, cfg.DM_34d, cfg.DM_36d,
        cfg.M_mass_f, cfg.M_mass_f_sil, cfg.M_ST, cfg.M_33d, cfg.M_34d, cfg.M_36d,
        cfg.F_mass_f, cfg.F_ST, cfg.F_33d, cfg.F_34d, cfg.F_36d,
        cfg.SS_mass_f, cfg.SS_ST, cfg.SS_33d, cfg.SS_34d, cfg.SS_36d,
        cfg.S_ST, cfg.S_33d, cfg.S_34d, cfg.S_36d,
        cfg.rate_f, cfg.oscillate, cfg.resurf_cm_yr, cfg.sil_mag_S, cfg.thick_C,
        cfg.f_S2, cfg.f_pl2mo, cfg.f_sq, cfg.f_deep, cfg.f_remobilised, cfg.f_pu,
        cfg.nofrac, cfg.MAF, cfg.MIF_eq
    );

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    cout << "\nExecution time: " << elapsed.count() << " seconds\n";
    return 0;
}
