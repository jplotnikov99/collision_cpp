#include "../include/utils.hpp"

double SQR(const double &x) { return x * x; }

double SIGN(const double a, const double b) {
    assert(std::abs(b) != 0);
    return std::abs(a / b) * b;
}

double generate_random(const double a, const double b) {
    double random = rand();
    random /= (double)RAND_MAX;
    return a + (b - a) * random;
}