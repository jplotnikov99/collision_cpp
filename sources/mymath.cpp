#include "mymath.hpp"

std::vector<double> solv::operator+(std::vector<double> a, std::vector<double> b)
{
    std::vector<double> sum;
    int k = a.size();

    assert(k == b.size());

    // std::cout << "test\n";

    for (size_t i = 0; i < k; i++)
    {
        sum.push_back(a.at(i) + b.at(i));
    }
    return sum;
}

std::vector<double> solv::operator-(std::vector<double> a, std::vector<double> b)
{
    std::vector<double> diff;
    int k = a.size();

    assert(k == b.size());

    // std::cout << "test\n";

    for (size_t i = 0; i < k; i++)
    {
        diff.push_back(a.at(i) - b.at(i));
    }
    return diff;
}

std::vector<double> solv::operator/(std::vector<double> a, std::vector<double> b)
{
    std::vector<double> frac;
    int k = a.size();

    assert(k == b.size());

    // std::cout << "test\n";

    for (size_t i = 0; i < k; i++)
    {
        frac.push_back(a.at(i) / b.at(i));
    }
    return frac;
}

std::vector<double> solv::operator/(std::vector<double> a, double b)
{
    std::vector<double> res;
    
    assert(b != 0);

    for (auto it : a)
    {
        res.push_back(it / b);
    }

    return res;
}

std::vector<double> solv::vabs(std::vector<double> a)
{
    std::vector<double> abs;

    for (auto it : a)
    {
        abs.push_back(fabs(it));
    }
    return abs;
}

double solv::max(std::vector<double> a)
{
    double max = 0;
    double temp;
    for (auto it : a)
    {
        temp = it;
        if (temp > max)
            max = temp;
    }
    return max;
}