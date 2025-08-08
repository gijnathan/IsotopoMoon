#include <iostream>
#include "headers/isotope_evolution.hpp"
using namespace std;

int main() {
    double N0 = 100.0;
    double lambda = 0.01;
    double t = 50.0;

    double N = compute_decay(N0, lambda, t);
    cout << "Remaining amount after " << t << " units of time: " << N << endl;

    return 0;
}

