#include "amp.hpp"

double CollisionInt::En(const double &po, const double &pz, const double &m) {
    return sqrt(po * po + pz * pz + m * m);
}

double CollisionInt::lips(const double &po, const double &En) {
    return 0.002015720902074968 * po / En;
}

double CollisionInt::feq(const double &x, const double &s) {
    if (x > 30)
        return exp(-x);
    else
        return 1 / (exp(x) + s);
}

double CollisionInt::Pf(const double &pz3) {
    return feq(E1 / T, 1.) * feq(E2 / T, -1.) * (1. - feq(E3 / T, 1.)) *
           (1 + feq(E4 / T, -1.));
}

double CollisionInt::DR(const double &pz3) {
    return (pz1 + pz2 - pz3) /
               sqrt((m4 * m4) + (pz1 + pz2 - pz3) * (pz1 + pz2 - pz3) +
                    SQR(po1 * cos(ph1) + po2 * cos(ph2) - po3 * cos(ph3)) +
                    SQR(po1 * sin(ph1) + po2 * sin(ph2) - po3 * sin(ph3))) -
           pz3 / sqrt((m3 * m3) + (po3 * po3) + (pz3 * pz3));
}

void CollisionInt::pz3sol() {
    const double cos12 = cos(ph1 - ph2);
    const double cos13 = cos(ph1 - ph3);
    const double cos23 = cos(ph2 - ph3);
    const double den = 2. * (E1 + E2 - pz1 - pz2) * (E1 + E2 + pz1 + pz2);
    const double num1 =
        2. * (pz1 + pz2) *
            (po1 * po2 * cos12 - po3 * (po1 * cos13 + po2 * cos23)) +
        (pz1 + pz2) * (-((E1 + E2) * (E1 + E2)) - (m3 * m3) + (m4 * m4) +
                       (po1 * po1) + (po2 * po2) + ((pz1 + pz2) * (pz1 + pz2)));
    const double num2 = sqrt(
        ((E1 + E2) * (E1 + E2)) *
        (4 * E2 * (E1 * E1 * E1) + (E1 * E1 * E1 * E1) + (E2 * E2 * E2 * E2) +
         (m3 * m3 * m3 * m3) - 2 * (m3 * m3) * (m4 * m4) + (m4 * m4 * m4 * m4) -
         2 * (m3 * m3) * (po1 * po1) + 2 * (m4 * m4) * (po1 * po1) +
         (po1 * po1 * po1 * po1) - 2 * (m3 * m3) * (po2 * po2) +
         2 * (m4 * m4) * (po2 * po2) + 4 * (po1 * po1) * (po2 * po2) +
         (po2 * po2 * po2 * po2) + 2 * (po1 * po1) * (po3 * po3) +
         2 * (po2 * po2) * (po3 * po3) + 2 * (m3 * m3) * (pz1 * pz1) +
         2 * (m4 * m4) * (pz1 * pz1) + 2 * (po1 * po1) * (pz1 * pz1) +
         2 * (po2 * po2) * (pz1 * pz1) + 4 * (po3 * po3) * (pz1 * pz1) +
         4 * pz1 * pz2 *
             ((m3 * m3) + (m4 * m4) + (po1 * po1) + (po2 * po2) +
              2 * (po3 * po3) + (pz1 * pz1)) +
         (pz1 * pz1 * pz1 * pz1) +
         2 *
             ((m3 * m3) + (m4 * m4) + (po1 * po1) + (po2 * po2) +
              2 * (po3 * po3) + 3 * (pz1 * pz1)) *
             (pz2 * pz2) +
         4 * pz1 * (pz2 * pz2 * pz2) + (pz2 * pz2 * pz2 * pz2) +
         4 * E1 * E2 *
             ((E2 * E2) - (m3 * m3) - (m4 * m4) - (po1 * po1) - (po2 * po2) -
              2 * (po3 * po3) - ((pz1 + pz2) * (pz1 + pz2))) -
         2 * (E2 * E2) *
             ((m3 * m3) + (m4 * m4) + (po1 * po1) + (po2 * po2) +
              2 * (po3 * po3) + ((pz1 + pz2) * (pz1 + pz2))) +
         (E1 * E1) * (6 * (E2 * E2) -
                      2 * ((m3 * m3) + (m4 * m4) + (po1 * po1) + (po2 * po2) +
                           2 * (po3 * po3) + ((pz1 + pz2) * (pz1 + pz2)))) +
         2 * (cos(2 * (ph1 - ph2)) * (po1 * po1) * (po2 * po2) +
              2 * po1 * po2 * cos12 *
                  (-2 * po3 * (po1 * cos13 + po2 * cos23) -
                   ((E1 + E2) * (E1 + E2)) - (m3 * m3) + (m4 * m4) +
                   (po1 * po1) + (po2 * po2) + (po3 * po3) +
                   ((pz1 + pz2) * (pz1 + pz2))) +
              po3 *
                  (po3 * cos(2 * ph3) *
                       (2 * po1 * po2 * cos(ph1 + ph2) +
                        cos(2 * ph1) * (po1 * po1) +
                        cos(2 * ph2) * (po2 * po2)) -
                   2 * (po1 * cos(ph1) + po2 * cos(ph2)) * cos(ph3) *
                       (-((E1 + E2) * (E1 + E2)) - (m3 * m3) + (m4 * m4) +
                        (po1 * po1) + (po2 * po2) +
                        ((pz1 + pz2) * (pz1 + pz2))) -
                   2 * (-((E1 + E2) * (E1 + E2)) - (m3 * m3) + (m4 * m4) + (po1 * po1) + (po2 * po2) + ((pz1 + pz2) * (pz1 + pz2))) *
                       (po1 * sin(ph1) + po2 * sin(ph2)) * sin(ph3) +
                   2 * po3 * (po1 * cos(ph1) + po2 * cos(ph2)) *
                       (po1 * sin(ph1) + po2 * sin(ph2)) * sin(2 * ph3)))));

    pz31 = -(num1 + num2) / den;
    pz32 = -(num1 - num2) / den;
}

double CollisionInt::amp(const double &pz3) {
    const double s = 2. * (E1 * E2 - po1 * po2 * cos(ph1 - ph2) - pz1 * pz2) +
                     m1 * m1 + m2 * m2;
    const double t = -2. * (E1 * E3 - po1 * po3 * cos(ph1 - ph3) - pz1 * pz3) +
                     m1 * m1 + m1 * m1;
    return (
        (-4 * (el * el) * (gs * gs) * (m1 * m1) *
         (2 * (m2 * m2 * m2 * m2) * (s - (m1 * m1)) * (t - (m1 * m1)) +
          (s + t) * (s * (m4 * m4 * m4 * m4) +
                     t * (((m4 * m4) + s) * ((m4 * m4) + s)) + s * (t * t)) -
          (m2 * m2) * (2 * s * t * (s + t) +
                       2 * (m1 * m1) *
                           (2 * s * t - 3 * (s + t) * (m4 * m4) +
                            (m1 * m1) * (-5 * (s + t) + 3 * (m4 * m4) +
                                         6 * (m1 * m1))) +
                       (m4 * m4) * (4 * s * t + (s * s) + (t * t))) +
          (m1 * m1) *
              (-4 * (s + t) * (m4 * m4 * m4 * m4) -
               (s + t) * ((s - t) * (s - t)) +
               4 * (m4 * m4) * (3 * s * t + (s * s) + (t * t)) +
               (m1 * m1) * (-2 * s * t - 26 * (s + t) * (m4 * m4) +
                            4 * (m4 * m4 * m4 * m4) +
                            4 * (m1 * m1) *
                                (-4 * (s + t) + 7 * (m4 * m4) + 5 * (m1 * m1)) +
                            5 * (s * s) + 5 * (t * t))))) /
        (3. * (mW * mW) * ((-(m1 * m1) + s) * (-(m1 * m1) + s)) * (sW * sW) *
         ((-(m1 * m1) + t) * (-(m1 * m1) + t))));
}

double CollisionInt::operator()(const double &po11, const double &ph11,
                                const double &pz11, const double &po22,
                                const double &ph22, const double &pz22,
                                const double &po33, const double &ph33) {
    double res = 0;
    po1 = po11;
    ph1 = ph11;
    pz1 = pz11;
    po2 = po22;
    ph2 = ph22;
    pz2 = pz22;
    po3 = po33;
    ph3 = ph33;
    E1 = En(po1, pz1, m1);
    E2 = En(po2, pz2, m2);
    pz3sol();
    E3 = En(po3, pz31, m1);
    E4 = E1 + E2 - E3;
    res += lips(po1, E1) * lips(po2, E2) * lips(po3, E3) * lips(1., E4) *
           SQR(SQR(2 * M_PI)) / std::abs(DR(pz31)) * amp(pz31) * Pf(pz31);
    E3 = En(po3, pz32, m1);
    E4 = E1 + E2 - E3;
    res += lips(po1, E1) * lips(po2, E2) * lips(po3, E3) * lips(1., E4) *
           SQR(SQR(2 * M_PI)) / std::abs(DR(pz32)) * amp(pz32) * Pf(pz32);
    return res;
}