#include "../include/amp.hpp"

void CollisionInt::set_phase_space(const double &po11, const double &ph11,
                                   const double &pz11, const double &po22,
                                   const double &ph22, const double &pz22,
                                   const double &ph33) {
    po1 = po11;
    ph1 = ph11;
    pz1 = pz11;
    po2 = po22;
    ph2 = ph22;
    pz2 = pz22;
    ph3 = ph33;
    E1 = En(po1, pz1, m1);
    E2 = En(po2, pz2, m2);
    cos12 = cos(ph1 - ph2);
    cos13 = cos(ph1 - ph3);
    cos23 = cos(ph2 - ph3);
}

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
    return feq(E1, 1.) * feq(E2, -1.) * (1. - feq(E3, 1.)) * (1 + feq(E4, -1.));
}

double CollisionInt::DR(const double &pz3) {
    return (pz1 + pz2 - pz3) /
               sqrt((m4 * m4) + (pz1 + pz2 - pz3) * (pz1 + pz2 - pz3) +
                    p4x * p4x + p4y * p4y) -
           pz3 / sqrt((m3 * m3) + (po3 * po3) + (pz3 * pz3));
}

void CollisionInt::pz3sol() {
    const double a = m3 * m3 + po3 * po3;
    const double b = m4 * m4 + p4x * p4x + p4y * p4y;
    const double c = pz1 + pz2;
    const double d = E1 + E2;
    const double den = 2. * (c * c - d * d);
    const double num1 = -c * (a - b - c * c + d * d);
    const double num2 = d * sqrt(SQR(a - b + c * c - d * d) + 2. * b * den);
    pz31 = (num1 + num2) / den;
    pz32 = (num1 - num2) / den;
}

void CollisionInt::po3sol(double &zero1, double &zero2) {
    const double a = m3 * m3 + m4 * m4 + po1 * po1 + po2 * po2;
    const double b = (pz1 + pz2) * (pz1 + pz2);
    const double po12 = po1 * po1;
    const double po22 = po2 * po2;
    const double pz12 = pz1 * pz1;
    double den = (4 * po1 * po2 * cos13 * cos23 - 2 * ((E1 + E2) * (E1 + E2)) +
                  po12 + cos(2 * (ph1 - ph3)) * po12 + po22 +
                  cos(2 * (ph2 - ph3)) * po22 + 2 * b);
    double num1 = ((po1 * cos13 + po2 * cos23) *
                   (2 * po1 * po2 * cos12 - ((E1 + E2) * (E1 + E2)) -
                    2. * (m3 * m3) + a + b));
    double num2 =
        (E1 + E2 - pz1 - pz2) * (E1 + E2 + pz1 + pz2) *
        (4 * E2 * (E1 * E1 * E1) + (E1 * E1 * E1 * E1) + (E2 * E2 * E2 * E2) +
         (m3 * m3 * m3 * m3) - 2 * (m3 * m3) * (m4 * m4) + (m4 * m4 * m4 * m4) +
         2 * (m4 * m4) * po12 + po12 * po12 + 2 * (m4 * m4) * po22 +
         4 * po12 * po22 + 2 * cos(2 * (ph1 - ph2)) * po12 * po22 +
         2 * (m3 * m3) *
             (2 * po1 * po2 * cos(ph1 + ph2 - 2 * ph3) +
              cos(2 * (ph1 - ph3)) * po12 + cos(2 * (ph2 - ph3)) * po22) +
         po22 * po22 + 2 * (m3 * m3) * pz12 + 2 * (m4 * m4) * pz12 +
         2 * po12 * pz12 + 2 * po22 * pz12 + 4 * pz1 * pz2 * (a + pz12) +
         pz12 * pz12 + 2 * (a + 3 * pz12) * (pz2 * pz2) +
         4 * pz1 * (pz2 * pz2 * pz2) + (pz2 * pz2 * pz2 * pz2) +
         4 * E1 * E2 * ((E2 * E2) - a - b) +
         4 * po1 * po2 * cos12 * (-((E1 + E2) * (E1 + E2)) + a - m3 * m3 + b) -
         2 * (E2 * E2) * (a + b) + (E1 * E1) * (6 * (E2 * E2) - 2 * (a + b)));
    if (num2 < 0) {
        zero1 = 0.;
        zero2 = 0.;
    } else {
        zero1 = (num1 - sqrt(num2)) / den;
        zero2 = (num1 + sqrt(num2)) / den;
    }
}

double CollisionInt::amp(const double &pz3) {
    if (massive) {
        const double s =
            2. * (E1 * E2 - po1 * po2 * cos12 - pz1 * pz2) + m1 * m1 + m2 * m2;
        const double t =
            -2. * (E1 * E3 - po1 * po3 * cos13 - pz1 * pz3) + m1 * m1 + m1 * m1;
        return 4 * el * el * gs * gs * mt_pole * mt_pole *
               (2 * (s - m1 * m1) * (-t + m1 * m1) * m2 * m2 * m2 * m2 +
                m2 * m2 *
                    (2 * s * t * (s + t) +
                     2 * m1 * m1 *
                         (2 * s * t - (s + t) * m4 * m4 +
                          m1 * m1 * (-5 * (s + t) + 6 * m1 * m1 + m4 * m4)) +
                     m4 * m4 * (s * s + t * t)) -
                (s + t - 2 * m1 * m1) *
                    (s * t * (s + t) - 10 * m1 * m1 * m1 * m1 * m1 * m1 +
                     m1 * m1 * m1 * m1 * (3 * (s + t) - 8 * m4 * m4) -
                     m1 * m1 *
                         (-4 * s * t - 4 * (s + t) * m4 * m4 + s * s +
                          t * t))) /
               (3. * mW * mW * (-m1 * m1 + s) * (-m1 * m1 + s) * sW * sW *
                (-m1 * m1 + t) * (-m1 * m1 + t));
    } else {
        const double s = 2. * (E1 * E2 - po1 * po2 * cos12 - pz1 * pz2);
        const double t = -2. * (E1 * E3 - po1 * po3 * cos13 - pz1 * pz3);
        return -4 * el * el * gs * gs * mt_pole * mt_pole *
               (-(s + t) * (s + 2 * t) * mtinf * mtinf +
                t * mtinf * mtinf * mtinf * mtinf + t * (s + t) * (s + t)) /
               (3. * s * mW * mW * sW * sW *
                ((mtinf * mtinf - t) * (mtinf * mtinf - t)));
    }
}

double CollisionInt::operator()(const double &po33) {
    double res = 0;

    po3 = po33;
    p4x = po1 * cos(ph1) + po2 * cos(ph2) - po3 * cos(ph3);
    p4y = po1 * sin(ph1) + po2 * sin(ph2) - po3 * sin(ph3);
    pz3sol();
    E3 = En(po3, pz31, m1);
    E4 = sqrt(m4 * m4 + p4x * p4x + p4y * p4y + SQR(pz1 + pz2 - pz31));
    if (std::abs(E1 + E2 - E3 - E4) < E1 * 1e-8)
        res += lips(po1, E1) * lips(po2, E2) * lips(po3, E3) * lips(1., E4) *
               SQR(SQR(2 * M_PI)) / std::abs(DR(pz31)) * amp(pz31) * Pf(pz31);
    if (res < 0) return 0.;
    E3 = En(po3, pz32, m1);
    E4 = sqrt(m4 * m4 + p4x * p4x + p4y * p4y + SQR(pz1 + pz2 - pz32));
    if (std::abs(E1 + E2 - E3 - E4) < E1 * 1e-8)
        res += lips(po1, E1) * lips(po2, E2) * lips(po3, E3) * lips(1., E4) *
               SQR(SQR(2 * M_PI)) / std::abs(DR(pz32)) * amp(pz32) * Pf(pz32);
    if (res < 0) return 0.;

    return res;
}

double IntMonte::operator()(const double &po11, const double &ph11,
                            const double &pz11, const double &po22,
                            const double &ph22, const double &pz22,
                            const double &ph33) {
    double minpo3, maxpo3;
    double res;
    Col.set_phase_space(po11, ph11, pz11, po22, ph22, pz22, ph33);
    Col.po3sol(maxpo3, minpo3);
    if (maxpo3 == minpo3) return 0.;
    if (minpo3 < 0) minpo3 = 0.;
    res = adap_gauss_kronrod_15(Col, minpo3, maxpo3, 1e-6);
    // std::cout << po11 << " " << ph11 << " " << pz11 << " " << po22 << " "
    //           << ph22 << " " << pz22 << " " << ph33 << " " << res << "\n";
    return res;
}