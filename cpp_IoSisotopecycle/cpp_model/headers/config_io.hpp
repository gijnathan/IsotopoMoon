#pragma once
#include "config.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>

using namespace std;

inline string trim(string s) {
    auto notspace = [](unsigned char c){ return !isspace(c); };
    s.erase(s.begin(), find_if(s.begin(), s.end(), notspace));
    s.erase(find_if(s.rbegin(), s.rend(), notspace).base(), s.end());
    return s;
}

inline MaybeDouble parse_maybe_double(const string& v) {
    // expects "Y,0.1" or "N,0"
    auto comma = v.find(',');
    if (comma == string::npos)
        throw runtime_error("MaybeDouble must look like Y,0.1 or N,0");

    string flag = trim(v.substr(0, comma));
    string val  = trim(v.substr(comma + 1));
    double d = stod(val);
    return MaybeDouble{flag, d};
}

inline ModelConfig load_config(const string& path) {
    ifstream in(path);
    if (!in) throw runtime_error("Could not open config file: " + path);

    ModelConfig cfg;

    string line;
    int lineno = 0;

    while (getline(in, line)) {
        ++lineno;

        // strip comments (# ...)
        auto hash = line.find('#');
        if (hash != string::npos) line = line.substr(0, hash);

        line = trim(line);
        if (line.empty()) continue;

        auto eq = line.find('=');
        if (eq == string::npos)
            throw runtime_error("Config parse error at line " + to_string(lineno) + ": missing '='");

        string key = trim(line.substr(0, eq));
        string val = trim(line.substr(eq + 1));

        // rates: rate_f.xxx
        if (key.rfind("rate_f.", 0) == 0) {
            string proc = key.substr(string("rate_f.").size());
            cfg.rate_f[proc] = stod(val);
            continue;
        }

        auto set_double = [&](double& target){ target = stod(val); };
        auto set_string = [&](string& target){ target = val; };
        auto set_maybe  = [&](MaybeDouble& target){ target = parse_maybe_double(val); };

        // map keys to fields
        if      (key == "t_step_Myr") set_double(cfg.t_step_Myr);
        else if (key == "end_time_Myr") set_double(cfg.end_time_Myr);
        else if (key == "nr_step") set_double(cfg.nr_step);
        else if (key == "nr_tol") set_double(cfg.nr_tol);

        else if (key == "DM_mass_f") set_maybe(cfg.DM_mass_f);
        else if (key == "DM_mass_f_sil") set_maybe(cfg.DM_mass_f_sil);
        else if (key == "DM_ST") set_maybe(cfg.DM_ST);
        else if (key == "DM_33d") set_double(cfg.DM_33d);
        else if (key == "DM_34d") set_double(cfg.DM_34d);
        else if (key == "DM_36d") set_double(cfg.DM_36d);

        else if (key == "M_mass_f") set_maybe(cfg.M_mass_f);
        else if (key == "M_mass_f_sil") set_maybe(cfg.M_mass_f_sil);
        else if (key == "M_ST") set_maybe(cfg.M_ST);
        else if (key == "M_33d") set_double(cfg.M_33d);
        else if (key == "M_34d") set_double(cfg.M_34d);
        else if (key == "M_36d") set_double(cfg.M_36d);

        else if (key == "F_mass_f") set_maybe(cfg.F_mass_f);
        else if (key == "F_ST") set_maybe(cfg.F_ST);
        else if (key == "F_33d") set_double(cfg.F_33d);
        else if (key == "F_34d") set_double(cfg.F_34d);
        else if (key == "F_36d") set_double(cfg.F_36d);

        else if (key == "SS_mass_f") set_maybe(cfg.SS_mass_f);
        else if (key == "SS_ST") set_maybe(cfg.SS_ST);
        else if (key == "SS_33d") set_double(cfg.SS_33d);
        else if (key == "SS_34d") set_double(cfg.SS_34d);
        else if (key == "SS_36d") set_double(cfg.SS_36d);

        else if (key == "S_ST") set_maybe(cfg.S_ST);
        else if (key == "S_33d") set_double(cfg.S_33d);
        else if (key == "S_34d") set_double(cfg.S_34d);
        else if (key == "S_36d") set_double(cfg.S_36d);

        else if (key == "MIF_eq") set_string(cfg.MIF_eq);
        else if (key == "nofrac") set_string(cfg.nofrac);
        else if (key == "MAF") set_string(cfg.MAF);

        else if (key == "oscillate") set_string(cfg.oscillate);
        else if (key == "resurf_cm_yr") set_double(cfg.resurf_cm_yr);
        else if (key == "sil_mag_S") set_double(cfg.sil_mag_S);
        else if (key == "thick_C") set_double(cfg.thick_C);

        else if (key == "f_S2") set_double(cfg.f_S2);
        else if (key == "f_pl2mo") set_double(cfg.f_pl2mo);
        else if (key == "f_sq") set_double(cfg.f_sq);
        else if (key == "f_deep") set_double(cfg.f_deep);
        else if (key == "f_remobilised") set_double(cfg.f_remobilised);
        else if (key == "f_pu") set_double(cfg.f_pu);

        else {
            throw runtime_error("Unknown key at line " + to_string(lineno) + ": " + key);
        }
    }

    cfg.validate();
    return cfg;
}
