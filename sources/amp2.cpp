#include "amp.hpp"

double amp(const double &po1, const double &ph1, const double &pz1,
           const double &po2, const double &ph2, const double &pz2,
           const double &po3, const double &ph3) {
    return ((4 * (el * el) * (gs * gs) * (mt_pole * mt_pole) *
             (2 * (s - (m1 * m1)) * (-t + (m1 * m1)) * (m2 * m2 * m2 * m2) +
              (m2 * m2) * (2 * s * t * (s + t) +
                           2 * (m1 * m1) *
                               (2 * s * t - (s + t) * (m4 * m4) +
                                (m1 * m1) * (-5 * (s + t) + 6 * (m1 * m1) +
                                             (m4 * m4))) +
                           (m4 * m4) * ((s * s) + (t * t))) -
              (s + t - 2 * (m1 * m1)) *
                  (s * t * (s + t) - 10 * (m1 * m1 * m1 * m1 * m1 * m1) +
                   (m1 * m1 * m1 * m1) * (3 * (s + t) - 8 * (m4 * m4)) -
                   (m1 * m1) * (-4 * s * t - 4 * (s + t) * (m4 * m4) + (s * s) +
                                (t * t))))) /
            (3. * (mW * mW) * ((-(m1 * m1) + s) * (-(m1 * m1) + s)) *
             (sW * sW) * ((-(m1 * m1) + t) * (-(m1 * m1) + t))));
}
