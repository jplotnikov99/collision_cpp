#pragma once

#include <cmath>
#include <complex>
#include <iostream>

#include "constants.hpp"
#include "integration_methods.hpp"
#include "utils.hpp"

double feq(const double &x, const double &s);

double ReL1(const double &w, const double &k, const double &p);

double ImL1(const double &w, const double &k, const double &p);

double ReL2(const double &w, const double &k, const double &p);

double ImL2(const double &w, const double &k, const double &p);

class T1 {
   private:
    const double w, k;
    bool isreal = true;

   public:
    T1(const double &ww, const double &kk) : w(ww), k(kk) {};
    void setre(const bool &re);
    double operator()(const double &p);
    ~T1() {};
};

class T2 {
   private:
    double w, k;
    bool isreal = true;

   public:
    T2(const double &ww, const double &kk) : w(ww), k(kk) {};
    void setre(const bool &re);
    double operator()(const double &p);
    ~T2() {};
};

class CollisionIntM0 {
   private:
    enum DenType { mass, htl, full };
    DenType dentype = full;
    const std::complex<double> I = {0., 1.};
    const double mtinf = gs / sqrt(6.);
    double p1p2n, p1p3n, p2p3n;

   public:
    CollisionIntM0() {};
    double lips(const double &p, const double &sth);
    double Pf(const double &E1, const double &E2, const double &E3,
              const double &E4);
    std::complex<double> den(const double &p1, const double &p3,
                             const double &t);
    double amp(const double &p1, const double &p2, const double &p3);
    double operator()(const double &p1, const double &ph1, const double &th1,
                      const double &p2, const double &ph2, const double &th2,
                      const double &ph3, const double &th3);
    ~CollisionIntM0() {};
};
