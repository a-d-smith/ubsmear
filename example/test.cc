#include "ubsmear.h"

#include <iostream>

int main()
{
    // Load the input matrix
    const auto matrix = ubsmear::UBFileHelper::ReadMatrix("matrix2.txt");

    // Get the eigenvalues & eigenvectors
    const auto [eigenvalues, eigenvectorMatrix] = ubsmear::UBMatrixHelper::GetEigenvectorMatrix(matrix, 1e-6);

    std::cout << "Input matrix" << std::endl;
    matrix.Print();

    std::cout << "Eigenvalues" << std::endl;
    eigenvalues.Print();

    std::cout << "Eigenvector matrix" << std::endl;
    eigenvectorMatrix.Print();

    // Reconstruct the original matrix from its eigenvalues and eigenvectors
    auto eigenvalueMatrix = ubsmear::UBMatrixHelper::GetUnitMatrix(matrix.GetRows());
    for (unsigned i = 0; i < matrix.GetRows(); ++i)
        eigenvalueMatrix.SetElement(i, i, eigenvalues.At(i, 0));

    const auto reconstructedMatrix = eigenvectorMatrix * eigenvalueMatrix * eigenvectorMatrix.GetTranspose();

    std::cout << "Reconstructed matrix" << std::endl;
    reconstructedMatrix.Print();

    // Get the difference between the input and the reconstructed matrix
    const auto differenceMatrix = reconstructedMatrix - matrix;
    std::cout << "Difference matrix" << std::endl;
    differenceMatrix.Print();

    return 0;
}
