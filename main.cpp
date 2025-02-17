#include <iostream>

#include "include/amp.hpp"

// BENCHMARK KINEMATICS
static double po1 = 0.0957420393030742;
static double ph1 = 3.308298445140734;
static double pz1 = 0.0019461535323972434;
static double po2 = 0.07809375969212062;
static double ph2 = 4.233660395861847;
static double pz2 = 0.009562805013332953;
static double po3 = 1.;
static double ph3 = 0.8607512110484707;
static double pz3 = -7.03336;

int main() {
    const double T = 10;
    const double alphaS = 0.12;
    const double gs = sqrt(4 * M_PI * alphaS);
    const double mH = 125.09;
    const double mt = 20;
    const double mg = sqrt(2) * gs * T;

    CollisionInt C(mt, mH, mt, mg);
    std::cout << C(po1, ph1, pz1, po2, ph2, pz2, po3, ph3) << "\n";

    return 0;
}