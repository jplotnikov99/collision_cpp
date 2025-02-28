#pragma once

#include <cmath>
#include <iostream>

#include "constants.hpp"
#include "integration_methods.hpp"
#include "utils.hpp"

class CollisionIntM0 {
   private:
    const bool massive = false;
    const double mtinf = gs / sqrt(6.);
    double cos12;
    double cos13;
    double cos23;
    double pz31, pz32;
    double p1p2n, p1p3n, p2p3n;

   public:
    CollisionIntM0() {};
    double lips(const double &p,const double &sth);
    double feq(const double &x, const double &s);
    double Pf(const double &E1, const double &E2, const double &E3,
              const double &E4);
    double amp(const double &p1, const double &p2, const double &p3);
    double operator()(const double &p1, const double &ph1, const double &th1,
                      const double &p2, const double &ph2, const double &th2,
                      const double &ph3, const double &th3);
    ~CollisionIntM0() {};
};
