#ifndef HELPERS_HPP
#define HELPERS_HPP

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <stdexcept>

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

// Function to update SulfurState struct

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


// helpers to write csv rows with padding
inline void write_csv_row_padded(ofstream& out,
                                const vector<string>& fields,
                                size_t total_cols)
{
    // write given fields
    for (size_t i = 0; i < fields.size(); ++i) {
        if (i) out << ",";
        out << fields[i];
    }
    // pad with empty fields to reach total_cols
    for (size_t i = fields.size(); i < total_cols; ++i) {
        out << ",";
        // empty cell => nothing between commas
    }
    out << "\n";
}

inline string csv_num(double v, int prec = 16) {
    ostringstream oss;
    oss.setf(ios::fixed);
    oss << setprecision(prec) << v;
    return oss.str();
}



/// make a Python-like header for sulfur isotope CSV output for what is being output each timestep
inline vector<string> sulfur_python_like_header() {
    return {
        "time step","time (Myr)","[S] mantle",
        "d34S M","d34S F","d34S SS","d34S S","d34S DM","d34S og","d34S iS",
        "d33S M","d33S F","d33S SS","d33S S","d33S DM","d33S og","d33S iS",
        "d36S M","d36S F","d36S SS","d36S S","d36S DM","d36S og","d36S iS",
        "D33S M","D33S F","D33S SS","D33S S","D33S DM","D33S og","D33S iS",
        "D36S M","D36S F","D36S SS","D36S S","D36S DM","D36S og","D36S iS",
        "ST M","ST F","ST SS","ST S","ST DM",
        "log ST M","log ST F","log ST SS","log ST S","log ST DM",
        "F mo","F fr","F pr","F pu","F sn","F pd","F hg","F bu","F pl","F rs","F sq","F ao","F dm",
        "mo/mo+ca","bu/bu+pl",
        "mbST","mb32","mb33","mb34","mb36","",
        "R33 M","R33 F","R33 SS","R33 S","R33 DM","R33 og","R33 iS",
        "R34 M","R34 F","R34 SS","R34 S","R34 DM","R34 og","R34 iS",
        "R36 M","R36 F","R36 SS","R36 S","R36 DM","R36 og","R36 iS",
        "a33_bu","a34_bu","a36_bu",
        "a33_escape-burial","a34_escape-burial","a36_escape-burial",
        "","",""  // <- your Python header shows trailing ",,,"
    };
}


/// preamble that has the input parameters at the top of the CSV
inline void write_python_like_preamble(ofstream& out,
                                      const vector<string>& main_header,
                                      // your 16 preamble values:
                                      double rate_pd, double rate_pi, double rate_ei, double rate_ed, double rate_ac,
                                      double rate_rc, double rate_ec,
                                      double resurf_cm_yr, double sil_mag_S, double thick_C,
                                      double f_remobilised, double f_S2, double f_pl2mo, double f_sq, double f_deep,
                                      double M_mass)
{
    size_t total_cols = main_header.size();

    vector<string> labels = {
        "rate factor pd","rate factor pi","rate factor ei","rate factor ed","rate factor ac",
        "rate factor rc","rate factor ec","resurfacing rate (cm/yr)","sil mag S (wf)",
        "crustal thickness (m)","fraction of f_crust_returned that is remolised","S2 to SO2",
        "plutons from mantle melting","sulfate sequestration","deep mantle - fraction of mantle melting",
        "mass silicate mantle"
    };

    vector<string> values = {
        csv_num(rate_pd), csv_num(rate_pi), csv_num(rate_ei), csv_num(rate_ed), csv_num(rate_ac),
        csv_num(rate_rc), csv_num(rate_ec),
        csv_num(resurf_cm_yr), csv_num(sil_mag_S), csv_num(thick_C),
        csv_num(f_remobilised), csv_num(f_S2), csv_num(f_pl2mo), csv_num(f_sq), csv_num(f_deep),
        // IMPORTANT: if you want scientific notation like Python sometimes produces,
        // you can special-case M_mass formatting. Otherwise fixed is fine.
        csv_num(M_mass)
    };

    // Row 0: labels + padding commas
    write_csv_row_padded(out, labels, total_cols);

    // Row 1: values + padding commas
    write_csv_row_padded(out, values, total_cols);

    // Row 2: main header (already full width)
    write_csv_row_padded(out, main_header, total_cols);
}

// INIITIAL STATE OUTPUT of SulfurState in Python-like format
inline vector<string> sulfur_initial_row_python_like(double t_step, double t_Myr,
                                                     double S_conc_M,
                                                     double M_34d, double M_33d, double M_36d,
                                                     double S_34d, double S_33d, double S_36d,
                                                     double M_33D, double M_36D,
                                                     double S_33D, double S_36D,
                                                     double M_ST, double S_ST,
                                                     double logM_ST, double logS_ST,
                                                     double mbST) {
    vector<string> row;

    // time + conc
    row.push_back(csv_num(t_step));
    row.push_back(csv_num(t_Myr));
    row.push_back(csv_num(S_conc_M));

    // d34S M, F, SS, S, DM, og, iS  (Python uses many blanks here)
    row.push_back(csv_num(M_34d));
    row.push_back(""); // F_34d
    row.push_back(""); // SS_34d
    row.push_back(csv_num(S_34d));
    row.push_back(""); // DM_34d
    row.push_back(""); // og
    row.push_back(""); // iS

    // d33S
    row.push_back(csv_num(M_33d));
    row.push_back("");
    row.push_back("");
    row.push_back(csv_num(S_33d));
    row.push_back("");
    row.push_back("");
    row.push_back("");

    // d36S
    row.push_back(csv_num(M_36d));
    row.push_back("");
    row.push_back("");
    row.push_back(csv_num(S_36d));
    row.push_back("");
    row.push_back("");
    row.push_back("");

    // D33S
    row.push_back(csv_num(M_33D));
    row.push_back("");
    row.push_back("");
    row.push_back(csv_num(S_33D));
    row.push_back("");
    row.push_back("");
    row.push_back("");

    // D36S
    row.push_back(csv_num(M_36D));
    row.push_back("");
    row.push_back("");
    row.push_back(csv_num(S_36D));
    row.push_back("");
    row.push_back("");
    row.push_back("");

    // ST M, F, SS, S, DM
    row.push_back(csv_num(M_ST));
    row.push_back("0.0");
    row.push_back("0.0");
    row.push_back(csv_num(S_ST));
    row.push_back("0.0");

    // log ST M, F, SS, S, DM
    row.push_back(csv_num(logM_ST));
    row.push_back("");
    row.push_back("");
    row.push_back(csv_num(logS_ST));
    row.push_back("");

    // Fluxes F mo..F dm (13 fields) -> blanks like python
    for (int i = 0; i < 13; i++) row.push_back("");

    // mo/mo+ca, bu/bu+pl
    row.push_back("");
    row.push_back("");

    // mbST, mb32, mb33, mb34, mb36, blank
    row.push_back(csv_num(mbST));
    row.push_back(""); row.push_back(""); row.push_back(""); row.push_back("");
    row.push_back("");

    // R33, R34, R36 blocks (21 fields total) -> blanks
    for (int i = 0; i < 21; i++) row.push_back("");

    // a33_bu,a34_bu,a36_bu,a33_eb,a34_eb,a36_eb
    for (int i = 0; i < 6; i++) row.push_back("");

    // trailing ",,,"
    row.push_back(""); row.push_back(""); row.push_back("");

    return row;
}


/// Convert one SulfurState to one "Python-like" CSV row (strings), matching the header ordering.
/// NOTE: We fill "og" columns with blanks "" (like your Python output),
/// and we add 3 trailing blanks to match the header's final ",,,".
inline vector<string> sulfur_state_to_python_like_row(const SulfurState& s) {
    vector<string> row;

    // time + conc
    row.push_back(csv_num(s.t_step));
    row.push_back(csv_num(s.t_step_Myr));
    row.push_back(csv_num(s.S_conc_M));

    // d34S: M,F,SS,S,DM,og(blank),iS
    row.push_back(csv_num(s.M_34d));
    row.push_back(csv_num(s.F_34d));
    row.push_back(csv_num(s.SS_34d));
    row.push_back(csv_num(s.S_34d));
    row.push_back(csv_num(s.DM_34d));
    row.push_back("");                 // d34S og (blank)
    row.push_back(csv_num(s.iS_34d));   // d34S iS

    // d33S
    row.push_back(csv_num(s.M_33d));
    row.push_back(csv_num(s.F_33d));
    row.push_back(csv_num(s.SS_33d));
    row.push_back(csv_num(s.S_33d));
    row.push_back(csv_num(s.DM_33d));
    row.push_back("");                 // d33S og
    row.push_back(csv_num(s.iS_33d));   // d33S iS

    // d36S
    row.push_back(csv_num(s.M_36d));
    row.push_back(csv_num(s.F_36d));
    row.push_back(csv_num(s.SS_36d));
    row.push_back(csv_num(s.S_36d));
    row.push_back(csv_num(s.DM_36d));
    row.push_back("");                 // d36S og
    row.push_back(csv_num(s.iS_36d));   // d36S iS

    // D33S
    row.push_back(csv_num(s.M_33D));
    row.push_back(csv_num(s.F_33D));
    row.push_back(csv_num(s.SS_33D));
    row.push_back(csv_num(s.S_33D));
    row.push_back(csv_num(s.DM_33D));
    row.push_back("");                 // D33S og
    row.push_back(csv_num(s.iS_33D));   // D33S iS

    // D36S
    row.push_back(csv_num(s.M_36D));
    row.push_back(csv_num(s.F_36D));
    row.push_back(csv_num(s.SS_36D));
    row.push_back(csv_num(s.S_36D));
    row.push_back(csv_num(s.DM_36D));
    row.push_back("");                 // D36S og
    row.push_back(csv_num(s.iS_36D));   // D36S iS

    // ST
    row.push_back(csv_num(s.M_ST));
    row.push_back(csv_num(s.F_ST));
    row.push_back(csv_num(s.SS_ST));
    row.push_back(csv_num(s.S_ST));
    row.push_back(csv_num(s.DM_ST));

    // log ST
    row.push_back(csv_num(s.logM_ST));
    row.push_back(csv_num(s.logF_ST));
    row.push_back(csv_num(s.logSS_ST));
    row.push_back(csv_num(s.logS_ST));
    row.push_back(csv_num(s.logDM_ST));

    // Fluxes
    row.push_back(csv_num(s.mo_F));
    row.push_back(csv_num(s.fr_F));
    row.push_back(csv_num(s.pr_F));
    row.push_back(csv_num(s.pu_F));
    row.push_back(csv_num(s.sn_F));

    row.push_back(csv_num(s.pd_F));
    row.push_back(csv_num(s.hg_F));
    row.push_back(csv_num(s.bu_F));
    row.push_back(csv_num(s.pl_F));
    row.push_back(csv_num(s.rs_F));

    row.push_back(csv_num(s.sq_F));
    row.push_back(csv_num(s.ao_F));
    row.push_back(csv_num(s.dm_F));

    // ratios
    row.push_back(csv_num(s.mo_mo_ca));
    row.push_back(csv_num(s.bu_bu_pl));

    // mass balance + blank column
    row.push_back(csv_num(s.mbST));
    row.push_back(csv_num(s.mb32S));
    row.push_back(csv_num(s.mb33S));
    row.push_back(csv_num(s.mb34S));
    row.push_back(csv_num(s.mb36S));
    row.push_back(""); // blank column after mb36

    // R33: M,F,SS,S,DM,og(blank),iS
    row.push_back(csv_num(s.M_33R));
    row.push_back(csv_num(s.F_33R));
    row.push_back(csv_num(s.SS_33R));
    row.push_back(csv_num(s.S_33R));
    row.push_back(csv_num(s.DM_33R));
    row.push_back(csv_num(s.ao_33R));  // you used ao_33R as "og" analog in your later Python
    row.push_back(csv_num(s.iS_33R));

    // R34
    row.push_back(csv_num(s.M_34R));
    row.push_back(csv_num(s.F_34R));
    row.push_back(csv_num(s.SS_34R));
    row.push_back(csv_num(s.S_34R));
    row.push_back(csv_num(s.DM_34R));
    row.push_back(csv_num(s.ao_34R));
    row.push_back(csv_num(s.iS_34R));

    // R36
    row.push_back(csv_num(s.M_36R));
    row.push_back(csv_num(s.F_36R));
    row.push_back(csv_num(s.SS_36R));
    row.push_back(csv_num(s.S_36R));
    row.push_back(csv_num(s.DM_36R));
    row.push_back(csv_num(s.ao_36R));
    row.push_back(csv_num(s.iS_36R));

    // alpha terms
    row.push_back(csv_num(s.a33_bu));
    row.push_back(csv_num(s.a34_bu));
    row.push_back(csv_num(s.a36_bu));

    row.push_back(csv_num(s.a33_eb));
    row.push_back(csv_num(s.a34_eb));
    row.push_back(csv_num(s.a36_eb));

    // match your header's trailing ",,,"
    row.push_back("");
    row.push_back("");
    row.push_back("");

    return row;
}





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
