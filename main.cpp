#include <chrono>
#include <fstream>
#include <iostream>

#include "cuba.h"
#include "include/amp.hpp"
#include "include/constants.hpp"

#define NDIM 7
#define NCOMP 1
#define USERDATA NULL
#define NVEC 1
#define EPSREL 0
#define EPSABS 0
#define VERBOSE 2
#define LAST 4
#define SEED 0
#define MINEVAL 0
#define MAXEVAL 1e4

#define NSTART 1e4
#define NINCREASE 0
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL
#define SPIN NULL

#define NNEW 1000
#define NMIN 2
#define FLATNESS 25.

#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0

double T = 100;
double mH = sqrt(11. / 6.) * el / sW * T;
double mt = gs * T / sqrt(6);
double mg = sqrt(2) * gs * T;

static int Integrand(const int *ndim, const cubareal xx[], const int *ncomp,
                     cubareal ff[], void *userdata) {
    IntMonte Im(mt, mH, mt, mg, T);
    double po1 = (1 - xx[0]) / xx[0];
    double ph1 = xx[1] * 2 * M_PI;
    double pz1 = (1 - xx[2]) / xx[2];
    double po2 = (1 - xx[3]) / xx[3];
    double ph2 = xx[4] * 2 * M_PI;
    double pz2 = (1 - xx[5]) / xx[5];
    double ph3 = xx[6] * 2 * M_PI;
    double jacobian = SQR(2 * M_PI) * 2 * M_PI /
                      (SQR(xx[0]) * SQR(xx[2]) * SQR(xx[3]) * SQR(xx[5]));
    ff[0] = T * T * T * T *
            (Im(po1, ph1, pz1, po2, ph2, pz2, ph3) +
             Im(po1, ph1, pz1, po2, ph2, -pz2, ph3) +
             Im(po1, ph1, -pz1, po2, ph2, pz2, ph3) +
             Im(po1, ph1, -pz1, po2, ph2, -pz2, ph3)) *
            jacobian;

    return 0;
}

void warm_up_vegas(integrand_t integrand, int points, int iterations, int phase,
                   cubareal integral[], cubareal error[], cubareal prob[]) {
    int last_only = 1, smoothing = 1;
    int comp, nregions, neval, fail;
    int maxpoints = iterations * points;
    int cuba_flags = VERBOSE + last_only * 4 + smoothing * 8;
    // char *state = "my_state";
    Vegas(NDIM, NCOMP, Integrand, USERDATA, NVEC, EPSREL, EPSABS, cuba_flags,
          SEED, MINEVAL, maxpoints, points, NINCREASE, NBATCH, phase, STATEFILE,
          SPIN, &neval, &fail, integral, error, prob);
}

void gridded_vegas(integrand_t integrand, int points, int iterations, int phase,
                   cubareal integral[], cubareal error[], cubareal prob[]) {
    int last_only = 0, smoothing = 1;
    int comp, nregions, neval, fail;
    int maxpoints = iterations * points;
    int cuba_flags = VERBOSE + last_only * 4 + smoothing * 8;
    // char *state = "my_state";
    Vegas(NDIM, NCOMP, Integrand, USERDATA, NVEC, EPSREL, EPSABS, cuba_flags,
          SEED, MINEVAL, maxpoints, points, NINCREASE, NBATCH, phase, STATEFILE,
          SPIN, &neval, &fail, integral, error, prob);
}

int main() {
    using namespace std::chrono;

    int comp, nregions, neval, fail;
    cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];
    std::ofstream outfile("T_dependence_2.dat", std::ios::out | std::ios::app);
    auto start = high_resolution_clock::now();
    // for (double i = 0; i < 100; i++) {
    T = 100;
    mH = sqrt(11. / 6.) * el / sW;
    mt = gs / sqrt(6);
    mg = sqrt(2) * gs;
    warm_up_vegas(Integrand, 1e4, 20, -1, integral, error, prob);
    gridded_vegas(Integrand, 1e5, 5, 1, integral, error, prob);
    for (size_t i = 0; i < NCOMP; i++) {
        outfile << 3. * integral[i] / (2. * SQR(M_PI) * pow(T, 4)) << "\t"
                << error[i] << "\n";
    }
    //}

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    //outfile << duration.count() << "\n";
    std::cout << "Computation time:\n" << duration.count() << "ms" << std::endl;

    //}
    outfile.close();
    return 0;
}