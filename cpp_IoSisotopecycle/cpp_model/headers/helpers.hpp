#ifndef HELPERS_HPP
#define HELPERS_HPP

#include <string>
#include <chrono>
#include <ctime>
using namespace std;

inline string current_datetime_string() {
    auto now = chrono::system_clock::now();
    time_t now_time = chrono::system_clock::to_time_t(now);
    string result = ctime(&now_time);
    result.pop_back(); // Remove trailing newline
    return result;
}

// Struct for conditional double values
struct MaybeDouble {
    string flag;  // "N" means not yet initialized, "Y" means ready
    double value;
};

// Helper to quickly create MaybeDouble objects
inline MaybeDouble make_maybe(const string& flag, double val = 0.0) {
    return {flag, val};
}

// IN C++ WE WILL MAKE A STRUCT TO HOLD ALL THE VARIABLES AND THEN OUTPUT TO A CSV AT THE END

struct SulfurState {
    // Time tracking
    double t_step = 0.0;       // Integer step
    double t_step_Myr = 0.0;   // Time in Myr

    // Mantle sulfur concentration (mol/kg)
    double S_conc_M = 0.0;

    // δ values (‰) - delta values by reservoir
    double M_34d = 0.0, F_34d = 0.0, SS_34d = 0.0, S_34d = 0.0, DM_34d = 0.0;
    double M_33d = 0.0, F_33d = 0.0, SS_33d = 0.0, S_33d = 0.0, DM_33d = 0.0;
    double M_36d = 0.0, F_36d = 0.0, SS_36d = 0.0, S_36d = 0.0, DM_36d = 0.0;

    // Δ values (‰) - capital D
    double M_33D = 0.0, F_33D = 0.0, SS_33D = 0.0, S_33D = 0.0, DM_33D = 0.0;
    double M_36D = 0.0, F_36D = 0.0, SS_36D = 0.0, S_36D = 0.0, DM_36D = 0.0;

    // ST (total sulfur) per reservoir
    double M_ST = 0.0, F_ST = 0.0, SS_ST = 0.0, S_ST = 0.0, DM_ST = 0.0;

    // log(ST) per reservoir
    double logM_ST = 0.0, logF_ST = 0.0, logSS_ST = 0.0, logS_ST = 0.0, logDM_ST = 0.0;

    // Mass balance terms
    double mbST = 0.0, mb32S = 0.0, mb33S = 0.0, mb34S = 0.0, mb36S = 0.0;

    // Reference isotope ratios
    double M_33R = 0.0, F_33R = 0.0, SS_33R = 0.0, S_33R = 0.0, DM_33R = 0.0;
    double M_34R = 0.0, F_34R = 0.0, SS_34R = 0.0, S_34R = 0.0, DM_34R = 0.0;
    double M_36R = 0.0, F_36R = 0.0, SS_36R = 0.0, S_36R = 0.0, DM_36R = 0.0;

    // Optional exists to include the following to mirror the extra Python columns:
    // - Escape-burial fractionation factors (a33, a34, a36)
    double ao_34d = 0.0, iS_34d = 0.0;
    double ao_33d = 0.0, iS_33d = 0.0;
    double ao_36d = 0.0, iS_36d = 0.0;
    double ao_33D = 0.0, iS_33D = 0.0;
    double ao_36D = 0.0, iS_36D = 0.0;

    // - Fluxes (e.g., F_mo, F_pd, etc.)
    double mo_F = 0.0, fr_F = 0.0, pr_F = 0.0, pu_F = 0.0, sn_F = 0.0;
    double pd_F = 0.0, hg_F = 0.0, bu_F = 0.0, pl_F = 0.0, rs_F = 0.0;
    double sq_F = 0.0, ao_F = 0.0, dm_F = 0.0;

    // - Derived ratios like mo/(mo+ca), bu/(bu+pl)
    double mo_mo_ca = 0.0, bu_bu_pl = 0.0;

    // reference ratios: 
    double ao_33R = 0.0, iS_33R = 0.0;
    double ao_34R = 0.0, iS_34R = 0.0;
    double ao_36R = 0.0, iS_36R = 0.0;

    // escape burial alpha values: 
    double a33_bu = 0.0, a34_bu = 0.0, a36_bu = 0.0;
    double a33_eb = 0.0, a34_eb = 0.0, a36_eb = 0.0;

};

void update_sulfur_state(SulfurState& state,
                         double t_step, double t_step_Myr, double S_conc_M,

                         // δ values
                         double M_34d, double F_34d, double SS_34d, double S_34d, double DM_34d,
                         double M_33d, double F_33d, double SS_33d, double S_33d, double DM_33d,
                         double M_36d, double F_36d, double SS_36d, double S_36d, double DM_36d,

                         // Δ values
                         double M_33D, double F_33D, double SS_33D, double S_33D, double DM_33D,
                         double M_36D, double F_36D, double SS_36D, double S_36D, double DM_36D,

                         // ST per reservoir
                         double M_ST, double F_ST, double SS_ST, double S_ST, double DM_ST,

                         // log(ST) per reservoir
                         double logM_ST, double logF_ST, double logSS_ST, double logS_ST, double logDM_ST,

                         // Mass balance
                         double mbST, double mb32S, double mb33S, double mb34S, double mb36S,

                         // Reference isotope ratios
                         double M_33R, double F_33R, double SS_33R, double S_33R, double DM_33R,
                         double M_34R, double F_34R, double SS_34R, double S_34R, double DM_34R,
                         double M_36R, double F_36R, double SS_36R, double S_36R, double DM_36R,

                         // Escape-burial fractionation δ
                         double ao_34d, double iS_34d,
                         double ao_33d, double iS_33d,
                         double ao_36d, double iS_36d,
                         double ao_33D, double iS_33D,
                         double ao_36D, double iS_36D,

                         // Fluxes
                         double mo_F, double fr_F, double pr_F, double pu_F, double sn_F,
                         double pd_F, double hg_F, double bu_F, double pl_F, double rs_F,
                         double sq_F, double ao_F, double dm_F,

                         // Derived ratios
                         double mo_mo_ca, double bu_bu_pl,

                         // More reference ratios
                         double ao_33R, double iS_33R,
                         double ao_34R, double iS_34R,
                         double ao_36R, double iS_36R,

                         // Escape-burial alpha
                         double a33_bu, double a34_bu, double a36_bu,
                         double a33_eb, double a34_eb, double a36_eb) 
{
    state.t_step = t_step;
    state.t_step_Myr = t_step_Myr;
    state.S_conc_M = S_conc_M;

    state.M_34d = M_34d; state.F_34d = F_34d; state.SS_34d = SS_34d; state.S_34d = S_34d; state.DM_34d = DM_34d;
    state.M_33d = M_33d; state.F_33d = F_33d; state.SS_33d = SS_33d; state.S_33d = S_33d; state.DM_33d = DM_33d;
    state.M_36d = M_36d; state.F_36d = F_36d; state.SS_36d = SS_36d; state.S_36d = S_36d; state.DM_36d = DM_36d;

    state.M_33D = M_33D; state.F_33D = F_33D; state.SS_33D = SS_33D; state.S_33D = S_33D; state.DM_33D = DM_33D;
    state.M_36D = M_36D; state.F_36D = F_36D; state.SS_36D = SS_36D; state.S_36D = S_36D; state.DM_36D = DM_36D;

    state.M_ST = M_ST; state.F_ST = F_ST; state.SS_ST = SS_ST; state.S_ST = S_ST; state.DM_ST = DM_ST;
    state.logM_ST = logM_ST; state.logF_ST = logF_ST; state.logSS_ST = logSS_ST; state.logS_ST = logS_ST; state.logDM_ST = logDM_ST;

    state.mbST = mbST; state.mb32S = mb32S; state.mb33S = mb33S; state.mb34S = mb34S; state.mb36S = mb36S;

    state.M_33R = M_33R; state.F_33R = F_33R; state.SS_33R = SS_33R; state.S_33R = S_33R; state.DM_33R = DM_33R;
    state.M_34R = M_34R; state.F_34R = F_34R; state.SS_34R = SS_34R; state.S_34R = S_34R; state.DM_34R = DM_34R;
    state.M_36R = M_36R; state.F_36R = F_36R; state.SS_36R = SS_36R; state.S_36R = S_36R; state.DM_36R = DM_36R;

    state.ao_34d = ao_34d; state.iS_34d = iS_34d;
    state.ao_33d = ao_33d; state.iS_33d = iS_33d;
    state.ao_36d = ao_36d; state.iS_36d = iS_36d;
    state.ao_33D = ao_33D; state.iS_33D = iS_33D;
    state.ao_36D = ao_36D; state.iS_36D = iS_36D;

    state.mo_F = mo_F; state.fr_F = fr_F; state.pr_F = pr_F; state.pu_F = pu_F; state.sn_F = sn_F;
    state.pd_F = pd_F; state.hg_F = hg_F; state.bu_F = bu_F; state.pl_F = pl_F; state.rs_F = rs_F;
    state.sq_F = sq_F; state.ao_F = ao_F; state.dm_F = dm_F;

    state.mo_mo_ca = mo_mo_ca;
    state.bu_bu_pl = bu_bu_pl;

    state.ao_33R = ao_33R; state.iS_33R = iS_33R;
    state.ao_34R = ao_34R; state.iS_34R = iS_34R;
    state.ao_36R = ao_36R; state.iS_36R = iS_36R;

    state.a33_bu = a33_bu; state.a34_bu = a34_bu; state.a36_bu = a36_bu;
    state.a33_eb = a33_eb; state.a34_eb = a34_eb; state.a36_eb = a36_eb;
};

#define DEBUG_PRINT(label, value) \
    cout << "     " << label << " = " << setprecision(17) << value << endl; 


// // struct for reservoir properties
// struct Reservoir {
//     double mass;
//     double ST;
//     double d33S;
//     double d34S;
//     double d36S;
// };

// SulfurState update_sulfur_state(
//     double t_step, double t_step_Myr,
//     const Reservoir& M, const Reservoir& F, const Reservoir& SS, const Reservoir& S, const Reservoir& DM,
//     const string& MIF_eq
// ) {
//     SulfurState state;

//     // Time
//     state.t_step = t_step;
//     state.t_step_Myr = t_step_Myr;

//     // Mantle sulfur concentration
//     state.S_conc_M = S_conc(M.ST, M.mass);  // mol/kg

//     // δ-values (‰)
//     state.M_33d = M.d33S; state.F_33d = F.d33S; state.SS_33d = SS.d33S; state.S_33d = S.d33S; state.DM_33d = DM.d33S;
//     state.M_34d = M.d34S; state.F_34d = F.d34S; state.SS_34d = SS.d34S; state.S_34d = S.d34S; state.DM_34d = DM.d34S;
//     state.M_36d = M.d36S; state.F_36d = F.d36S; state.SS_36d = SS.d36S; state.S_36d = S.d36S; state.DM_36d = DM.d36S;

//     // Δ-values (‰)
//     state.M_33D = calc_capD33(M.d33S, M.d34S, MIF_eq);
//     state.F_33D = calc_capD33(F.d33S, F.d34S, MIF_eq);
//     state.SS_33D = calc_capD33(SS.d33S, SS.d34S, MIF_eq);
//     state.S_33D = calc_capD33(S.d33S, S.d34S, MIF_eq);
//     state.DM_33D = calc_capD33(DM.d33S, DM.d34S, MIF_eq);

//     state.M_36D = calc_capD36(M.d36S, M.d34S, MIF_eq);
//     state.F_36D = calc_capD36(F.d36S, F.d34S, MIF_eq);
//     state.SS_36D = calc_capD36(SS.d36S, SS.d34S, MIF_eq);
//     state.S_36D = calc_capD36(S.d36S, S.d34S, MIF_eq);
//     state.DM_36D = calc_capD36(DM.d36S, DM.d34S, MIF_eq);

//     // ST (mol S)
//     state.M_ST = M.ST; state.F_ST = F.ST; state.SS_ST = SS.ST; state.S_ST = S.ST; state.DM_ST = DM.ST;

//     // log(ST) — avoid log(0)
//     state.logM_ST  = (M.ST  > 0.0) ? log10(M.ST)  : -99.0;
//     state.logF_ST  = (F.ST  > 0.0) ? log10(F.ST)  : -99.0;
//     state.logSS_ST = (SS.ST > 0.0) ? log10(SS.ST) : -99.0;
//     state.logS_ST  = (S.ST  > 0.0) ? log10(S.ST)  : -99.0;
//     state.logDM_ST = (DM.ST > 0.0) ? log10(DM.ST) : -99.0;

//     // Mass balance checks (these can be updated if you track isotope masses explicitly)
//     state.mbST   = M.ST + F.ST + SS.ST + S.ST + DM.ST;

//     // Reference isotope ratios (‰ -> fractional R = δ/1000 + 1 * R_std)
//     state.M_33R = d2R(M.d33S, 33, VCDT);
//     state.F_33R = d2R(F.d33S, 33, VCDT);
//     state.SS_33R = d2R(SS.d33S, 33, VCDT);
//     state.S_33R = d2R(S.d33S, 33, VCDT);
//     state.DM_33R = d2R(DM.d33S, 33, VCDT);

//     state.M_34R = d2R(M.d34S, 34, VCDT);
//     state.F_34R = d2R(F.d34S, 34, VCDT);
//     state.SS_34R = d2R(SS.d34S, 34, VCDT);
//     state.S_34R = d2R(S.d34S, 34, VCDT);
//     state.DM_34R = d2R(DM.d34S, 34, VCDT);

//     state.M_36R = d2R(M.d36S, 36, VCDT);
//     state.F_36R = d2R(F.d36S, 36, VCDT);
//     state.SS_36R = d2R(SS.d36S, 36, VCDT);
//     state.S_36R = d2R(S.d36S, 36, VCDT);
//     state.DM_36R = d2R(DM.d36S, 36, VCDT);

//     return state;
// }



#endif // HELPERS_HPP
