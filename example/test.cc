#include "ubsmear.h"

#include <iostream>

int main()
{
    // Load the input matrix
    const auto matrix = ubsmear::UBFileHelper::ReadMatrix("matrix.txt");

    // Get the eigenvalues & eigenvectors
    const auto [eigenvalues, eigenvectorMatrix] = ubsmear::UBMatrixHelper::GetEigenDecomposition(matrix, 1e-6);

    // Reconstruct the original matrix from its eigenvalues and eigenvectors
    const auto reconstructedMatrix = ubsmear::UBMatrixHelper::GetMatrixFromEigenDecomposition(eigenvalues, eigenvectorMatrix);

    // Get the difference between the input and the reconstructed matrix
    const auto differenceMatrix = reconstructedMatrix - matrix;

    std::cout << "Input matrix" << std::endl;
    matrix.Print();

    std::cout << "Eigenvalues" << std::endl;
    eigenvalues.Print();

    std::cout << "Eigenvector matrix" << std::endl;
    eigenvectorMatrix.Print();

    std::cout << "Reconstructed matrix" << std::endl;
    reconstructedMatrix.Print();

    std::cout << "Difference matrix" << std::endl;
    differenceMatrix.Print();

    return 0;
}
