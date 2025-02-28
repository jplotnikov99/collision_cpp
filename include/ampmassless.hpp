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
    double po1, ph1, pz1, po2, ph2, pz2, po3, ph3, p4x, p4y, p4z;
    double E1, E2, E3, E4;

   public:
    CollisionIntM0() {};
    void set_phase_space(const double &po11, const double &ph11,
                         const double &pz11, const double &po22,
                         const double &ph22, const double &pz22,
                         const double &ph33);
    double En(const double &po, const double &pz);
    double lips(const double &po, const double &En);
    double feq(const double &x, const double &s);
    double Pf(const double &pz3);
    double DR(const double &pz3);
    void pz3sol();
    void po3sol(double &zero1, double &zero2);
    double amp(const double &pz3);
    double operator()(const double &po33);
    ~CollisionIntM0() {};
};

class IntMonteM0 {
   private:
    CollisionIntM0 Col;

   public:
    IntMonteM0() {};
    double operator()(const double &po11, const double &ph11,
                      const double &pz11, const double &po22,
                      const double &ph22, const double &pz22,
                      const double &ph33);
    ~IntMonteM0() {};
};