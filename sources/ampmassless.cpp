#include "../include/ampmassless.hpp"

double CollisionIntM0::lips(const double &p, const double &sth) {
    return 0.002015720902074968 * p * sth;
}

double CollisionIntM0::feq(const double &x, const double &s) {
    if (x > 40)
        return exp(-x);
    else
        return 1 / (exp(x) + s);
}

double CollisionIntM0::Pf(const double &E1, const double &E2, const double &E3,
                          const double &E4) {
    return feq(E1, 1.) * feq(E2, -1.) * (1. - feq(E3, 1.)) * (1 + feq(E4, -1.));
}

double CollisionIntM0::amp(const double &p1, const double &p2,
                           const double &p3) {
    const double s = 2. * p1 * p2 * (1. - p1p2n);
    const double t = -2. * p1 * p3 * (1. - p1p3n);
    const std::complex<double> D =
        (1. + a) * (1. + a) * t + 2. * (1. + a) * b * w + b * b;
    /* double res = (-2. * s * t * el * el * gs * gs * mt_pole * mt_pole) /
                 (SQR(t - mtinf * mtinf) * mW * mW * sW * sW); */
    std::complex<double> res = (-2. * s * t * el * el * gs * gs * mt_pole * mt_pole) /
                 (D * std::conj(D) * mW * mW * sW * sW);
    if (res.real() < 0) return 0.;
    return res.real();
}

double CollisionIntM0::operator()(const double &p1, const double &ph1,
                                  const double &th1, const double &p2,
                                  const double &ph2, const double &th2,
                                  const double &ph3, const double &th3) {
    double res = 0;
    const double s1 = sin(th1), s2 = sin(th2), s3 = sin(th3);
    p1p2n = cos(th1) * cos(th2) + cos(ph1 - ph2) * s1 * s2;
    p1p3n = cos(th1) * cos(th3) + cos(ph1 - ph3) * s1 * s3;
    p2p3n = cos(th2) * cos(th3) + cos(ph2 - ph3) * s2 * s3;
    const double F = p1 + p2 - (p1 * p1p3n + p2 * p2p3n);
    const double p3 = p1 * p2 * (1 - p1p2n) / F;
    w = std::abs(p3 - p1);
    k = std::sqrt(SQR(p1) + SQR(p3) - 2 * p1 * p3 * p1p3n);
    std::complex<double> logwk =
        std::log(std::abs((w + k) / (w - k))) + I * M_PI;
    a = SQR(mtinf) / SQR(k) * (1. - w / (2. * k) * logwk);
    b = SQR(mtinf) / k * (-w / k + (SQR(w) / SQR(k) - 1.) * logwk / 2.);
    if (p3 > p1 + p2) {
        std::cout << "test\n";
        exit(1);
    }
    const double p4 = p1 + p2 - p3;
    res = lips(p1, s1) * lips(p2, s2) * lips(p3, s3) * 2 * M_PI /
          std::abs(2. * F) * amp(p1, p2, p3) * Pf(p1, p2, p3, p4);
    return res;
}