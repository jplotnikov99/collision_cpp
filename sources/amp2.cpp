#include "amp.hpp"

double amp(const double &po1, const double &ph1, const double &pz1,
           const double &po2, const double &ph2, const double &pz2,
           const double &po3, const double &ph3) {
    return ((-4 * s * (el * el) * (gs * gs) * (mtpole * mtpole)) /
            (3. * t * (mW * mW) * (sW * sW)));
}
