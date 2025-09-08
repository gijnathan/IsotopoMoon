#ifndef ALPHASET_TYPES_HPP
#define ALPHASET_TYPES_HPP

#include <string>
#include <unordered_map>
using namespace std;

struct AlphaSet {
    unordered_map<string, double> a33;
    unordered_map<string, double> a34;
    unordered_map<string, double> a36;
};

#endif
