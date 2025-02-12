#include <iostream>

#include "include/framework.hpp"
int main(){

    Board b(12);
    findComb(b);
    std::cout << "The total number of solutions is: " << solutions << std::endl;

    return 0;
}