#ifndef ISOTOPE_EVOLUTION_HPP
#define ISOTOPE_EVOLUTION_HPP

#include <cmath>

// Function declarations and definitions

inline double compute_decay(double initial_amount, double decay_constant, double time) {
    return initial_amount * std::exp(-decay_constant * time);
}

// Add more functions as needed

#endif

