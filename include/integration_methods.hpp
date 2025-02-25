#pragma once
#include <iostream>

#include "abscissa.hpp"

template <class FUNC>
double adap_gauss_kronrod_15(FUNC &f, const double l, const double r,
                             const double err, const size_t depth = 0) {
    double I1, I2, y[15];
    double h = (r - l) / 2;
    for (int i = 0; i < 15; i++) {
        y[i] = f((kronx_15[i] + 1) * h + l);
    }
    I1 = h * (0.022935322010529224963732008058970 * (y[0] + y[14]) +
              0.063092092629978553290700663189204 * (y[1] + y[13]) +
              0.104790010322250183839876322541518 * (y[2] + y[12]) +
              0.140653259715525918745189590510238 * (y[3] + y[11]) +
              0.169004726639267902826583426598550 * (y[4] + y[10]) +
              0.190350578064785409913256402421014 * (y[5] + y[9]) +
              0.204432940075298892414161999234649 * (y[6] + y[8]) +
              0.209482141084727828012999174891714 * y[7]);
    if ((I1 == 0.) || (depth > 20)) {
        return I1;
    }
    I2 = h * (0.129484966168869693270611432679082 * (y[1] + y[13]) +
              0.279705391489276667901467771423780 * (y[3] + y[11]) +
              0.381830050505118944950369775488975 * (y[5] + y[9]) +
              0.417959183673469387755102040816327 * y[7]);

    if (std::abs((I2 / I1 - 1)) < err) {
        return I1;
    }
    double m = (l + r) / 2;
    return adap_gauss_kronrod_15(f, l, m, err, depth + 1) +
           adap_gauss_kronrod_15(f, m, r, err, depth + 1);
}