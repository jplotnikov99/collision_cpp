#include "../include/ampmassless.hpp"

double theta(const double &x) { return x < 0 ? 0. : 1.; }

double feq(const double &x, const double &s) {
    if (x > 40)
        return exp(-x);
    else
        return 1 / (exp(x) + s);
}

double ReL1(const double &w, const double &k, const double &p) {
    return log(std::abs((SQR(k + 2. * p) - w * w) / (SQR(k - 2. * p) - w * w)));
}

double ImL1(const double &w, const double &k, const double &p) {
    double res = 0.;
    double wp = (w + k) / 2.;
    double wm = (w - k) / 2.;
    if ((p + wp) / (p + wm) < 0.) res += M_PI;
    if ((p - wp) / (p - wm) < 0.) res -= M_PI;
    return res;
}

double ReL2(const double &w, const double &k, const double &p) {
    return std::log(
        std::abs((SQR(k + w) - 4. * p * p) / (SQR(k - w) - 4. * p * p)));
}

double ImL2(const double &w, const double &k, const double &p) {
    double res = 0.;
    double wp = (w + k) / 2.;
    double wm = (w - k) / 2.;
    if ((p + wp) / (p + wm) < 0.) res += M_PI;
    if ((p - wp) / (p - wm) < 0.) res += M_PI;
    return res;
}

void T1::setre(const bool &re) { isreal = re; }

double T1::operator()(const double &u) {
    double res = 0.;
    const double p = (1. - u) / u;
    double l1;
    double l2;
    double logwk = 0.;
    if (isreal) {
        logwk = log(std::abs((w + k) / (w - k)));
        l1 = ReL1(w, k, p);
        l2 = ReL2(w, k, p);
    } else {
        if ((w - k) < 0.) logwk = M_PI;
        l1 = ImL1(w, k, p);
        l2 = ImL2(w, k, p);
        res += (theta(k + 2. * p - w) * theta(-(k - w) * (k - 2. * p + w)) *
                    feq(std::abs(p - w), -1.) -
                theta(2. * p + w - k) * theta((k - w) * (k + 2. * p + w)) *
                    feq(std::abs(p + w), -1.)) *
               2. * M_PI * p * feq(p, 1.);
    }
    res += (p * (2. * logwk - l2) * feq(p, 1.) +
            (2. * p * logwk - p * l2 - w * l1) * feq(p, -1.));
    return res / (u * u);
}

void T2::setre(const bool &re) { isreal = re; }

double T2::operator()(const double &u) {
    double res = 0.;
    const double p = (1. - u) / u;
    double l1;
    if (isreal) {
        l1 = ReL1(w, k, p);
        return ((4. * p + (w * w - k * k) / (2. * k) * l1) * feq(p, 1.) +
                (4. * p - (w * w - k * k) / (2. * k) * l1) * feq(p, -1.)) /
               (u * u);
    } else {
        l1 = ImL1(w, k, p);
        /* res += (theta(k + 2. * p - w) * theta(-(k - w) * (k - 2. * p + w)) *
                    feq(std::abs(p - w), -1.) +
                theta(2. * p + w - k) * theta((k - w) * (k + 2. * p + w)) *
                    feq(std::abs(p + w), -1.)) *
               M_PI * feq(p, 1.) * (w * w - k * k) / k; */
        res += (w * w - k * k) / (2. * k) * l1 * (feq(p, 1.) - feq(p, -1.));
        return res / (u * u);
    }
}

double CollisionIntM0::lips(const double &p, const double &sth) {
    return 0.002015720902074968 * p * sth;
}

double CollisionIntM0::Pf(const double &E1, const double &E2, const double &E3,
                          const double &E4) {
    return feq(E1, 1.) * feq(E2, -1.) * (1. - feq(E3, 1.)) * (1 + feq(E4, -1.));
}

double CollisionIntM0::den(const double &p1, const double &p2, const double &p3,
                           const double &w, const double &k, const double &s,
                           const double &t) {
    std::complex<double> a, b, den, num;
    switch (dentype) {
        case mass:
            return s * t / (t - mtinf * mtinf);
            break;
        case htl: {
            std::complex<double> logwk = std::log(std::abs((w + k) / (w - k)));
            if (w - k < 0) logwk += I * M_PI;
            a = SQR(mtinf) / SQR(k) * (1. - w / (2. * k) * logwk);
            b = SQR(mtinf) / k * (-w / k + (SQR(w) / SQR(k) - 1.) * logwk / 2.);
            den = ((1. + a) * (1. + a) * t + 2. * (1. + a) * b * w + b * b);
            num = (1. + a) * std::conj(1. + a) * s * t +
                  b * std::conj(b) * (s + t - 4. * p2 * p3) +
                  2. * (a * std::conj(b) + b).real() *
                      (p1 * (s + t) + p2 * t - p3 * s) -
                  4. * (std::conj(a) * b + b).imag() * epse;
            return num.real() / (SQR(den.real()) + SQR(den.imag()));
        } break;
        case full: {
            std::complex<double> tu, tk;
            T1 t1(w, k);
            T2 t2(w, k);
            tu = adap_gauss_kronrod_15(t1, 0., 1., 1e-6);
            t1.setre(false);
            tu += I * adap_gauss_kronrod_15(t1, 0., 1., 1e-6);
            tu *= SQR(mtinf) / (k * M_PI * M_PI);
            tk = adap_gauss_kronrod_15(t2, 0., 1., 1e-6);
            t2.setre(false);
            tk += I * adap_gauss_kronrod_15(t2, 0., 1., 1e-6);
            tk *= SQR(mtinf) / (M_PI * M_PI);
            a = (tk - w * tu) / (k * k);
            b = (w * w / (k * k) - 1.) * tu - w / (k * k) * tk;
            den = ((1. + a) * (1. + a) * t + 2. * (1. + a) * b * w + b * b);
            num = (1. + a) * std::conj(1. + a) * s * t +
                  b * std::conj(b) * (s + t - 4. * p2 * p3) +
                  2. * (a * std::conj(b) + b).real() *
                      (p1 * (s + t) + p2 * t - p3 * s) -
                  4. * (std::conj(a) * b + b).imag() * epse;
            return num.real() / (SQR(den.real()) + SQR(den.imag()));
        } break;

        default:
            std::cout << "Wrong dentype in CollisionIntM0\n";
            exit(1);
            break;
    }
}

double CollisionIntM0::amp(const double &p1, const double &p2,
                           const double &p3) {
    const double s = 2. * p1 * p2 * (1. - p1p2n);
    const double t = -2. * p1 * p3 * (1. - p1p3n);
    const double w = std::abs(p1 - p3);
    const double k = std::sqrt(SQR(p1) + SQR(p3) - 2 * p1 * p3 * p1p3n);
    const double kins = den(p1, p2, p3, w, k, s, t);
    // const std::complex<double> Ds = den(p1, p2, p1p2n, s);
    // const double Dtsq = SQR(Dt.real()) + SQR(Dt.imag());
    // const double Dssq = 1. / (SQR(Ds.real()) + SQR(Ds.imag()));
    // const double Dts = (1. / (Dt * std::conj(Ds))).real();
    double res = (-2. * el * el * gs * gs * mt_pole * mt_pole) /
                 (mW * mW * sW * sW) * kins;
    if (res < 0) return 0.;
    return res;
}

double CollisionIntM0::operator()(const double &p1, const double &ph1,
                                  const double &th1, const double &p2,
                                  const double &ph2, const double &th2,
                                  const double &ph3, const double &th3) {
    double res = 0;
    const double s1 = sin(th1), s2 = sin(th2), s3 = sin(th3);
    const double c1 = cos(th1), c2 = cos(th2), c3 = cos(th3);
    p1p2n = c1 * c2 + cos(ph1 - ph2) * s1 * s2;
    p1p3n = c1 * c3 + cos(ph1 - ph3) * s1 * s3;
    p2p3n = c2 * c3 + cos(ph2 - ph3) * s2 * s3;
    epse = s3 * (c1 * s2 * sin(ph2 - ph3) - s1 * c2 * sin(ph1 - ph3)) +
           s1 * s2 * c3 * sin(ph1 - ph2);
    const double F = p1 + p2 - (p1 * p1p3n + p2 * p2p3n);
    const double p3 = p1 * p2 * (1 - p1p2n) / F;
    epse = epse * p1 * p2 * p3;
    if (p3 > p1 + p2) {
        std::cout << "test\n";
        exit(1);
    }
    const double p4 = p1 + p2 - p3;
    res = lips(p1, s1) * lips(p2, s2) * lips(p3, s3) * 2 * M_PI /
          std::abs(2. * F) * amp(p1, p2, p3) * Pf(p1, p2, p3, p4);
    return res;
}

double T1im::operator()(const double &x) {
    return (w + x) * feq(std::abs(x), -1) * feq(std::abs(w + x), 1) /
               feq(std::abs(w + x), -1) +
           x * feq(std::abs(x), 1);
}