#pragma once

#include <vector>
#include <cmath>
#include <assert.h>

namespace solv{
    std::vector<double> operator +(const std::vector<double> a,const  std::vector<double> b);
    std::vector<double> operator -(const std::vector<double> a,const  std::vector<double> b);
    std::vector<double> operator /(std::vector<double> a, std::vector<double> b);
    std::vector<double> operator /(std::vector<double> a, double b);
    std::vector<double> vabs(std::vector<double> a);
    double max(std::vector<double> a);

}