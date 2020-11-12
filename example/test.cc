#include "ubsmear.h"

#include <iostream>

int main()
{
    const auto matrix = ubsmear::UBFileHelper::ReadMatrix("matrix.txt");
    matrix.Print();

    const auto colVect = ubsmear::UBFileHelper::ReadColumnVector("colVect.txt");
    colVect.Print();

    return 0;
}
