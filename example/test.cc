#include "ubsmear.h"

#include <iostream>

int main()
{
    const auto matrix = ubsmear::UBFileHelper::ReadMatrix("matrix.txt");
    matrix.Print();
    std::cout << std::endl;

    const auto matrix2 = ubsmear::UBFileHelper::ReadMatrix("matrix2.txt");
    matrix2.Print();
    std::cout << std::endl;

    const auto newMatrix = matrix * matrix2;
    newMatrix.Print();

    return 0;
}
