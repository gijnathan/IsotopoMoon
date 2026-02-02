#pragma once
#include <string>
#include <unordered_map>
#include <stdexcept>
#include "helpers.hpp"
using namespace std;


inline bool is_on(const MaybeDouble& m) {
    return (m.flag == "Y" || m.flag == "y");
}

struct ModelConfig {
    // time + NR
    double t_step_Myr = 1e-1;
    double end_time_Myr = 4570.1;
    double nr_step = 1.0;
    double nr_tol  = 1e5;

    // reservoirs
    MaybeDouble DM_mass_f{"N", 0.0};
    MaybeDouble DM_mass_f_sil{"N", 0.0};
    MaybeDouble DM_ST{"N", 0.0};
    double DM_33d = 0.0, DM_34d = 0.0, DM_36d = 0.0;

    MaybeDouble M_mass_f{"Y", 0.1257};
    MaybeDouble M_mass_f_sil{"Y", 0.1257};
    MaybeDouble M_ST{"N", 0.0};
    double M_33d = 0.0, M_34d = 0.0, M_36d = 0.0;

    MaybeDouble F_mass_f{"N", 0.0};
    MaybeDouble F_ST{"N", 0.0};
    double F_33d = 0.0, F_34d = 0.0, F_36d = 0.0;

    MaybeDouble SS_mass_f{"N", 0.0};
    MaybeDouble SS_ST{"N", 0.0};
    double SS_33d = 0.0, SS_34d = 0.0, SS_36d = 0.0;

    MaybeDouble S_ST{"N", 0.0};
    double S_33d = 0.0, S_34d = 0.0, S_36d = 0.0;

    // rates (fractions)
    unordered_map<string,double> rate_f = {
        {"pd", 1.0}, {"pi", 1.0}, {"ei", 1.0},
        {"ed", 1.0}, {"ac", 1.0}, {"rc", 1.0}, {"ec", 1.0}
        // {"pu", 1.0}
    };

    // options
    string MIF_eq = "FW2003_pg3";
    string nofrac = "no";
    string MAF    = "no";

    string oscillate = "N";
    double resurf_cm_yr = 1.0;
    double sil_mag_S = 0.001;
    double thick_C = 40000.0;

    double f_S2 = 0.2;
    double f_pl2mo = 0.8;
    double f_sq = 0.25;
    double f_deep = 0.0;
    double f_remobilised = 1.0;
    double f_pu = 0.5;

    void validate() const {
        auto yn = [](const string& s){
            return s=="Y"||s=="y"||s=="N"||s=="n";
        };

        if (!yn(DM_mass_f.flag) || !yn(DM_mass_f_sil.flag) || !yn(DM_ST.flag) ||
            !yn(M_mass_f.flag)  || !yn(M_mass_f_sil.flag)  || !yn(M_ST.flag)  ||
            !yn(F_mass_f.flag)  || !yn(F_ST.flag)          ||
            !yn(SS_mass_f.flag) || !yn(SS_ST.flag)         ||
            !yn(S_ST.flag))
            throw runtime_error("One or more MaybeDouble flags are not Y/N.");

        if (t_step_Myr <= 0) throw runtime_error("t_step_Myr must be > 0.");
        if (end_time_Myr < 0) throw runtime_error("end_time_Myr must be >= 0.");
        if (nr_step <= 0) throw runtime_error("nr_step must be > 0.");
        if (thick_C <= 0) throw runtime_error("thick_C must be > 0.");

        auto in01 = [](double x){ return x >= 0.0 && x <= 1.0; };
        if (!in01(f_S2) || !in01(f_pl2mo) || !in01(f_sq) || !in01(f_deep) ||
            !in01(f_remobilised) || !in01(f_pu))
            throw runtime_error("One or more fractions are outside [0,1].");
    }
};
