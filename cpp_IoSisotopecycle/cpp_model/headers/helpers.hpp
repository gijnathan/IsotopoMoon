#ifndef HELPERS_HPP
#define HELPERS_HPP

#include <string>
using namespace std;

// Struct for conditional double values
struct MaybeDouble {
    string flag;  // "N" means not yet initialized, "Y" means ready
    double value;
};

// Helper to quickly create MaybeDouble objects
inline MaybeDouble make_maybe(const string& flag, double val = 0.0) {
    return {flag, val};
}

#endif // HELPERS_HPP
