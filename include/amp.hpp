#pragma once

#include <cmath>
#include <iostream>

#include "utils.hpp"

class CollisionInt {
   private:
    const double T = 10;
    const double alphaS = 0.12;
    const double el = 0.332;
    const double m1, m2, m3, m4;
    const double mW = 80.;
    const double gs = sqrt(4 * M_PI * alphaS);
    const double sW = sqrt(0.223);
    double pz31, pz32;
    double po1, ph1, pz1, po2, ph2, pz2, po3, ph3;
    double E1, E2, E3, E4;

   public:
    CollisionInt(const double &m11, const double &m22, const double &m33,
                 const double &m44)
        : m1(m11), m2(m22), m3(m33), m4(m44) {};
    double En(const double &po, const double &pz, const double &m);
    double lips(const double &po, const double &En);
    double feq(const double &x, const double &s);
    double Pf(const double &pz3);
    double DR(const double &pz3);
    void pz3sol();
    double amp(const double &pz3);
    double operator()(const double &po11, const double &ph11,
                      const double &pz11, const double &po22,
                      const double &ph22, const double &pz22,
                      const double &po33, const double &ph33);
    ~CollisionInt() {};
};
