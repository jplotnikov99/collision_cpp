#pragma once

#include <cmath>
#include <complex>
#include <iostream>

#include "constants.hpp"
#include "integration_methods.hpp"
#include "utils.hpp"

double theta(const double& x);

double feq(const double& x, const double& s);

double L1(const double& w, const double& k, const double& p);

double L2(const double& w, const double& k, const double& p);

double ReL1(const double& w, const double& k, const double& p);

double ImL1(const double& w, const double& k, const double& p);

double ReL2(const double& w, const double& k, const double& p);

double ImL2(const double& w, const double& k, const double& p);

class T1int {
   private:
    const double w, k;
    bool isreal = true;
    bool imint = true;

   public:
    T1int(const double& ww, const double& kk) : w(ww), k(kk) {};
    void setre(const bool& re);
    void setImInt(const bool& type);
    double operator()(const double& p);
    ~T1int() {};
};

class T2int {
   private:
    const double w, k;
    bool isreal = true;

   public:
    T2int(const double& ww, const double& kk) : w(ww), k(kk) {};
    void setre(const bool& re);
    double operator()(const double& p);
    ~T2int() {};
};

class T1_v2 {
   private:
    bool isreal = true;

   public:
    T1_v2() {};
    void setre(const bool& re);
    double operator()(const double& w, const double& k);
    ~T1_v2() {};
};

class T2_v2 {
   private:
    bool isreal = true;

   public:
    T2_v2() {};
    void setre(const bool& re);
    double operator()(const double& w, const double& k);
    ~T2_v2() {};
};

class T1 {
   private:
    const double w, k;
    bool isreal = true;

   public:
    T1(const double& ww, const double& kk) : w(ww), k(kk) {};
    void setre(const bool& re);
    double operator()(const double& p);
    ~T1() {};
};

class T2 {
   private:
    double w, k;
    bool isreal = true;

   public:
    T2(const double& ww, const double& kk) : w(ww), k(kk) {};
    void setre(const bool& re);
    double operator()(const double& p);
    ~T2() {};
};

class CollisionIntM0 {
   private:
    enum DenType { mass, htl, full };
    DenType dentype = full;
    const std::complex<double> I = {0., 1.};
    const double mtinf = gs / sqrt(6.);
    double p1p2n, p1p3n, p2p3n, epse;

   public:
    CollisionIntM0() {};
    double lips(const double& p, const double& sth);
    double Pf(const double& E1, const double& E2, const double& E3,
              const double& E4);
    double den(const double& p1, const double& p2, const double& p3,
               const double& w, const double& k, const double& s,
               const double& t);
    double amp(const double& p1, const double& p2, const double& p3);
    double operator()(const double& p1, const double& ph1, const double& th1,
                      const double& p2, const double& ph2, const double& th2,
                      const double& ph3, const double& th3);
    ~CollisionIntM0() {};
};

class T1im {
   private:
    const double w, k;

   public:
    T1im(const double ww, const double kk) : w(ww), k(kk) {};
    double operator()(const double& x);
    ~T1im() {};
};