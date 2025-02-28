#include "amp.hpp"

double amp(const double &po1, const double &ph1, const double &pz1,
           const double &po2, const double &ph2, const double &pz2,
           const double &po3, const double &ph3) {
    return ((-4 * (el * el) * (gs * gs) * (mtpole * mtpole) *
             (-((s + t) * (s + 2 * t) * (mtinf * mtinf)) +
              t * (mtinf * mtinf * mtinf * mtinf) + t * ((s + t) * (s + t)))) /
            (3. * s * (mW * mW) * (sW * sW) *
             (((mtinf * mtinf) - t) * ((mtinf * mtinf) - t))));
}
