#ifndef ALPHASET_HPP
#define ALPHASET_HPP

#include <string>
#include <unordered_map>
#include <vector>
#include "alphaset_types.hpp"
using namespace std;

// Forward declaration of alphas() so compiler knows it exists
double alphas(const string& process, int X, const string& nofrac, const string& MAF);

inline AlphaSet compute_all_alphas(const string& nofrac, const string& MAF) {
    AlphaSet alpha;

    vector<string> processes = {
        "gr", "th", "pu", "ed", "ac", "ei", "rc", "ec", "pi",
        "mm", "dg", "xt", "fr", "pd", "hg", "rm", "vp", "sq", "dm", "cf"
    };

    for (const auto& p : processes) {
        alpha.a33[p] = alphas(p, 33, nofrac, MAF);
        alpha.a34[p] = alphas(p, 34, nofrac, MAF);
        alpha.a36[p] = alphas(p, 36, nofrac, MAF);
    }

    // Composite processes
    alpha.a33["pl"] = alpha.a33["mm"] * alpha.a33["xt"];
    alpha.a34["pl"] = alpha.a34["mm"] * alpha.a34["xt"];
    alpha.a36["pl"] = alpha.a36["mm"] * alpha.a36["xt"];

    alpha.a33["mo"] = alpha.a33["mm"] * alpha.a33["dg"];
    alpha.a34["mo"] = alpha.a34["mm"] * alpha.a34["dg"];
    alpha.a36["mo"] = alpha.a36["mm"] * alpha.a36["dg"];

    return alpha;
}

#endif
