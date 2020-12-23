#ifndef UBSMEAR_UBSMEAR_HELPERS_UBMATRIXHELPER
#define UBSMEAR_UBSMEAR_HELPERS_UBMATRIXHELPER

#include "ubsmear/Objects/UBMatrix.h"

#include <stdexcept>
#include <cmath>
#include <algorithm>

namespace ubsmear
{

/**
* @brief A helper class for common tasks in linear algebra
*/
class UBMatrixHelper
{
    public:
        /**
        * @brief Get a matrix of the specified dimensions with zero for every elememt
        *
        * @param nRows the number of rows
        * @param nCols the number of columns
        *
        * @return the zero matrix
        */
        static UBMatrix GetZeroMatrix(const size_t nRows, const size_t nCols);

        /**
        * @brief Get a unit square matrix of the specified size
        *
        * @param size the size of the matrix (both rows and columns)
        *
        * @return the unit matrix
        */
        static UBMatrix GetUnitMatrix(const size_t size);

        /**
        * @brief Get a diagonal matrix (all off-diagonals are zero) with diagonal entries given by the elements of the input column vector
        *
        * @param columnVector the input column vector which defines the diagonals of the output matrix
        *
        * @return the diagonal matrix
        */
        static UBMatrix GetDiagonalMatrix(const UBMatrix &columnVector);

        /**
        * @brief Get a (square) Givens rotation matrix of the specified size, which acts to rotate by a specified angle in the plane
                 corresponding to the specified row and column indices
        *
        * @param size the size of the matrix (both rows and columns)
        * @param rowIndex the row index (first axis of the plane of rotation)
        * @param columnIndex the column index (second axis of the plane of rotation)
        * @param angle the angle through which to rotate (radians)
        *
        * @return the rotation matrix
        */
        static UBMatrix GetGivensRotationMatrix(const size_t size, const size_t rowIndex, const size_t columnIndex, const float angle);

        /**
        * @brief Get the matrix of whose columns are the eigenvectors of the input matrix. The eigenvectors are normalised to unity.
                 This method in an implementation of the Jacobi eigenvalue algorithm for real and symmetric matrices.
        *
        * @param matrix the input real symmetric matrix
        * @param precision the maximum allowed RMS of the off-diagonal elements of the matrix when transformed by the matrix of eigenvector
        *
        * @return a pair, first is the column vector of eigenvalues, second is the matrix of eigenvectors
        */
        static std::pair<UBMatrix, UBMatrix> GetEigenDecomposition(const UBMatrix &matrix, const float precision);

        /**
        * @brief Reconstruct a real symmetric matrix from it's eigenvectors and eigenvalues.
        *
        * @param eigenvalues the input column vector of eigenvalues
        * @param eigenvectorMatrix the input matrix whose columns are the eigenvectors
        *
        * @return the real symmetric matrix which has the supplied eigenvalues and eigenvectors
        */
        static UBMatrix GetMatrixFromEigenDecomposition(const UBMatrix &eigenvalues, const UBMatrix &eigenvectorMatrix);

        /**
        * @brief Reconstruct a real symmetric matrix from a partial set of it's eigenvectors and eigenvalues.
        *
        * @param eigenvalues the input column vectors of L eigenvalues
        * @param eigenvectorMatrix the input (N x L) matrix of eigenvectors, whose columns correspond the eigenvectors
        *
        * @return the (N x N) real symmetric matrix which has the supplied eigenvalues and eigenvectors
        */
        static UBMatrix GetMatrixFromPartialEigenDecomposition(const UBMatrix &eigenvalues, const UBMatrix &eigenvectorMatrix);

        /**
        * @brief Check if the input matrix is square
        *
        * @param matrix the input matrix
        *
        * @return boolean, true if square
        */
        static bool IsSquare(const UBMatrix &matrix);

        /**
        * @brief Check if the input matrix is symmetric
        *
        * @param matrix the input matrix
        *
        * @return boolean, true if symmetric
        */
        static bool IsSymmetric(const UBMatrix &matrix);

    private:

        /**
        * @brief Get the row and column of the off-diagonal element of the input matrix with the largest absolute value
        *
        * @param matrix the input matrix
        *
        * @return a pair (first = row index, second = column index), giving the index of the largest off-diagonal matrix element
        */
        static std::pair<size_t, size_t> GetLargestOffDiagonalElement(const UBMatrix &matrix);

        /**
        * @brief Get the sum of the squares of the off-diagonal terms in the input matrix
        *
        * @param matrix the input matrix
        *
        * @return the off-diagonal squared sum
        */
        static float GetOffDiagonalSquaredSum(const UBMatrix &matrix);
};

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBMatrix UBMatrixHelper::GetZeroMatrix(const size_t nRows, const size_t nCols)
{
    // Get a vector of zeros for the elements of the matrix
    std::vector<float> elements(nRows*nCols, 0.f);

    // Construct the matrix
    return UBMatrix(elements, nRows, nCols);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBMatrix UBMatrixHelper::GetUnitMatrix(const size_t size)
{
    // Get a zero matrix
    auto matrix = UBMatrixHelper::GetZeroMatrix(size, size);

    // Set the diagonals to be one
    for (size_t i = 0; i < size; ++i)
        matrix.SetElement(i, i, 1.f);

    return matrix;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBMatrix UBMatrixHelper::GetDiagonalMatrix(const UBMatrix &columnVector)
{
    // Check that the input matrix is really a column vector
    if (columnVector.GetColumns() != 1u)
        throw std::invalid_argument("UBMatrixHelper::GetDiagonalMatrix - Input matrix has more than one column, it should be a column vector");

    // Get the size of the output matrix
    const auto size = columnVector.GetRows();

    // Get a zero matrix
    auto matrix = UBMatrixHelper::GetZeroMatrix(size, size);

    // Set the diagonals
    for (size_t i = 0; i < size; ++i)
        matrix.SetElement(i, i, columnVector.At(i, 0));

    return matrix;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBMatrix UBMatrixHelper::GetGivensRotationMatrix(const size_t size, const size_t rowIndex, const size_t columnIndex, const float angle)
{
    // Check that the supplied row & column indices are not out of bounds
    if (rowIndex >= size)
        throw std::out_of_range("UBMatrixHelper::GetGivensRotationMatrix - supplied row index is larger than the supplied matrix size");

    if (columnIndex >= size)
        throw std::out_of_range("UBMatrixHelper::GetGivensRotationMatrix - supplied column index is larger than the supplied matrix size");

    if (rowIndex == columnIndex)
        throw std::invalid_argument("UBMatrixHelper::GetGivensRotationMatrix - supplied row and column indices are equal");

    // Get a unit matrix
    auto matrix = UBMatrixHelper::GetUnitMatrix(size);

    // Set the diagonals in the specified plane to be the cosine of the angle
    const auto cosTheta = std::cos(angle);
    matrix.SetElement(rowIndex, rowIndex, cosTheta);
    matrix.SetElement(columnIndex, columnIndex, cosTheta);

    // Set the off-diagonals in the specified plane to plus and minus the sine of the angle
    const auto sinTheta = std::sin(angle);
    matrix.SetElement(rowIndex, columnIndex, -sinTheta);
    matrix.SetElement(columnIndex, rowIndex, +sinTheta);

    return matrix;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline std::pair<UBMatrix, UBMatrix> UBMatrixHelper::GetEigenDecomposition(const UBMatrix &matrix, const float precision)
{
    // Here we implement the Jacobi eigenvalue algorithm. A summary of this algorithm is described below.
    //
    // The intention of this algorithm perfom eigenvalue decomposition, i.e. to express the input matrix, M in diagonal form: M = Q D Q^-1
    //  - D    = diagonal matrix whose elements are the eigenvalues of M
    //  - Q    = matrix whose columns are the eigenvectors of M
    //  - Q^-1 = the inverse matrix of Q
    //
    // If M is a real symmetric matrix (which we insist to be true), then the matrix of eigenvectors, Q, is an orthogonal matrix: Q^-1 = Q^T
    //  - Q^T  = the transform matrix of Q (switching rows for columns)
    //
    // The Jacobi eigenvalue algorithm is iterative. Each iteration (labelled with n), results in a new matrix D_n. This matrix is formed by
    // applying a rotation matrix, R_n, to the outcome of the previous iteration: D_n = R_n  D_n-1  R_n^T
    //  - D_n   = the outcome of the nth iteration of the algorithm
    //  - D_n-1 = the outcome of the (n-1)th iteration of the algorithm
    //  - R_n   = the rotation matrix used in the nth iteration
    //  - R_n^T = the transpose of the rotation matrix, R_n
    //  - D_0   = M, the algorithm is seeded with the input matrix
    //
    // Note that rotation matrices are orthogonal, and the matrix product of two orthogonal matrices is itself orthogonal. Hence, the matrix
    // product: (R_1 R_2 R_3 ... R_N) is an orthogonal matrix. The aim of the algorithm is to apply successive rotations such that D_n
    // converges toward a diagonal matrix. After N iterations, if the sequence has converged, then we can make the approximation:
    //  - D   = D_N
    //  - Q^T = R_N ... R_3 R_2 R_1
    //
    // In particular, R_n is chosen to be a Givens rotation matrix (i.e. a rotation in a plane defined by two coordinate axes). The plane
    // chosen to correspond to the off-diagonal element of D_n-1 with the largest absolute value (known as the pivot). The angle of rotation
    // is chosen such that this pivot element becomes zero after the rotation. The angle is given by:
    //
    // tan(2 theta) = 2 [ij] / ([jj] - [ii])
    //  - theta = angle of rotation
    //  - [ij] = the value of D_n-1 in the ith row and jth column. i and j are the row and column of the pivot.
    //
    // To determine if D_n has converged to a diagonal matrix, we consider the sum of the squares of its off-diagonal entries, X_n^2. If the
    // matrix has size s, then there are s(s-1) off-diagonal elements. We define: S = s(s-1)/2, where series of S rotations is known as a
    // "sweep". One can show that after S rotations, we have:
    //
    // X_{n+S} <= r X_n
    //  - X_n = square-root of the sum of the squares of the off-diagonal elements of D_n
    //  - X_{n+S} = the above, after a further S iterations of the algorithm
    //  - r = (1 - 1/S)^(S/2) = the worst-case-scenario linear convergence rate after one sweep
    //    - For large matrices, r tends to e^(1/2), which is approximately 0.61
    //  - S = half the number of the number of off-diagonal elements

    // Insist that the input matrix is symmetric
    if (!UBMatrixHelper::IsSymmetric(matrix))
        throw std::invalid_argument("UBMatrixHelper::GetEigenDecomposition - Input matrix isn't symmetric");

    // Get the size of the matrix (it is square)
    const auto size = matrix.GetRows();

    // If the size is one (i.e. the matrix is just a scalar, then the answer is trivial)
    if (size == 1)
    {
        return std::pair<UBMatrix, UBMatrix> (
            {{ matrix.At(0, 0) }, 1, 1},   // Eigenvalue is the entry of the matrix
            {{ 1.f }, 1, 1}                // Normalised eigenvector is just the number 1.0
        );
    }

    // Determine the maximum number of iterations we should use - start by assuming one
    size_t maxIterations = 1;
    float squaredSumTarget = -std::numeric_limits<float>::max(); // The target square of sum of diagonals to reach the desired precision

    // ATTN for size 1, the answer is trivial. For size 2, only one iteration is required to get the answer precicely
    if (size > 2)
    {
        // Get the sweep length
        const size_t sweep = size * (size - 1) / 2;

        // Get the initial squared sum of the off-diagonals
        const auto squaredSumInitial = UBMatrixHelper::GetOffDiagonalSquaredSum(matrix);

        // Get the off-diagonal squared sum we would expect at the target precision
        squaredSumTarget = 2 * sweep * std::pow(precision, 2.f);

        // ATTN if we are already within the desired precision, then just run one iteration
        if (squaredSumInitial > squaredSumTarget && squaredSumInitial > std::numeric_limits<float>::epsilon())
        {
            // Get the worst-case-scenario convergence rate for a single sweep
            const auto convergenceRate = std::pow((1.f - (1.f / static_cast<float>(sweep))), static_cast<float>(sweep) / 2.f);

            // Work out how many iterations we would need to reach the desired precision in the worst case scenario
            maxIterations = std::ceil(sweep * std::log(squaredSumTarget/squaredSumInitial) / std::log(convergenceRate));
        }
    }

    // Make a copy of the input matrix which will converge to be diagonal
    auto diagMatrix = matrix;

    // Make a matrix to hold the total transformation we applied to the matrix to make it diagonal (start with a unit matrix)
    auto transformationMatrix = UBMatrixHelper::GetUnitMatrix(size);

    // Repeat for each iteration of the algorithm
    for (size_t i = 0; i < maxIterations; ++i)
    {
        // Identify the pivot
        const auto &[pivotRow, pivotCol] = UBMatrixHelper::GetLargestOffDiagonalElement(diagMatrix);

        // Get the angle through which we should rotate
        const auto valueRowCol = diagMatrix.At(pivotRow, pivotCol);
        const auto valueRowRow = diagMatrix.At(pivotRow, pivotRow);
        const auto valueColCol = diagMatrix.At(pivotCol, pivotCol);
        const auto theta = std::atan2(2 * valueRowCol, valueColCol - valueRowRow) / 2;

        // Get the rotation matrix and it's transpose
        const auto rotationMatrix = UBMatrixHelper::GetGivensRotationMatrix(size, pivotRow, pivotCol, theta);
        const auto rotationMatrixTranspose = rotationMatrix.GetTranspose();

        // Apply the rotation
        diagMatrix = rotationMatrix * diagMatrix * rotationMatrixTranspose;
        transformationMatrix = rotationMatrix * transformationMatrix;

        // Get the squared sum of the off-diagonals after this iteration
        const auto squaredSum = UBMatrixHelper::GetOffDiagonalSquaredSum(diagMatrix);

        // If we have reached the desired precision, then stop
        if (squaredSum < squaredSumTarget)
            break;
    }

    // ATTN the Jacobi eigenvalue algorithm is now complete!
    // - The diagonals of diagMatrix are the eigenvalues
    // - The transpose of transformationMatrix is the matrix whose columns are the eigenvectors
    //
    // Below we put this information into a reproducible form, i.e. we sort the eigenvalues by size, and normalize the eigenvectors

    // Extract the eigenvalues from the diagonalised matrix along with their indices
    std::vector< std::pair<size_t, float> > eigenvalueIndexPairs;
    for (size_t i = 0; i < size; ++i)
        eigenvalueIndexPairs.emplace_back(i, diagMatrix.At(i, i));

    // Sort the eigenvalues by their magnitude (largest first)
    std::sort(eigenvalueIndexPairs.begin(), eigenvalueIndexPairs.end(), [] (const auto &a, const auto &b) { return a.second > b.second; });

    // Get the elements of the matrix whose rows are the eigenvectors of the input matrix - sorted by eigenvalue and normalised to unity
    std::vector<float> eigenvectorMatrixElements;
    std::vector<float> eigenvalues;
    for (const auto &[iRow, eigenvalue] : eigenvalueIndexPairs)
    {
        // ATTN we have a choice of the sign of the normalisation
        // Here we choose the sign such that the largest element of the eigenvector is positive
        float maxValue = -std::numeric_limits<float>::max();
        int sign = 1;
        float norm2 = 0.f;
        for (size_t iCol = 0; iCol < size; ++iCol)
        {
            // Sum the squares of the elements of the eigenvector
            const auto value = transformationMatrix.At(iRow, iCol);
            norm2 += std::pow(value, 2.f);

            // Get the sign of the largest element of the eigenvector
            if (std::abs(value) >= maxValue)
            {
                maxValue = std::abs(value);
                sign = (value > 0.f) ? 1 : -1;
            }
        }

        if (norm2 <= std::numeric_limits<float>::epsilon())
            throw std::logic_error("UBMatrixHelper::GetEigenDecomposition - Found an eigenvector with invalid norm");

        const auto norm = std::pow(norm2, 0.5f) * sign;

        // Normalise the elements of the eigenvector and store them
        for (size_t iCol = 0; iCol < size; ++iCol)
            eigenvectorMatrixElements.push_back(transformationMatrix.At(iRow, iCol) / norm);

        // Store the current eigenvalue
        eigenvalues.push_back(eigenvalue);
    }

    // Return the output matrices as a pair
    return std::pair<UBMatrix, UBMatrix> (
        { eigenvalues, size, 1 },
        UBMatrix(eigenvectorMatrixElements, size, size).GetTranspose()
    );
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBMatrix UBMatrixHelper::GetMatrixFromEigenDecomposition(const UBMatrix &eigenvalues, const UBMatrix &eigenvectorMatrix)
{
    // Check that the input matrices have consistent dimensions
    if (eigenvalues.GetColumns() != 1u)
        throw std::invalid_argument("UBMatrixHelper::GetMatrixFromEigenDecomposition - input eigenvalues is not a column vector");

    const auto size = eigenvalues.GetRows();
    if (eigenvectorMatrix.GetRows() != size)
        throw std::invalid_argument("UBMatrixHelper::GetMatrixFromEigenDecomposition - the number of rows of the input eigenvector matrix, doesn't match the number of supplied eigenvalues");

    if (eigenvectorMatrix.GetColumns() != size)
        throw std::invalid_argument("UBMatrixHelper::GetMatrixFromEigenDecomposition - the number of columns of the input eigenvector matrix, doesn't match the number of supplied eigenvalues");

    return UBMatrixHelper::GetMatrixFromPartialEigenDecomposition(eigenvalues, eigenvectorMatrix);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBMatrix UBMatrixHelper::GetMatrixFromPartialEigenDecomposition(const UBMatrix &eigenvalues, const UBMatrix &eigenvectorMatrix)
{
    // Check that the input matrices have consistent dimensions
    if (eigenvalues.GetColumns() != 1u)
        throw std::invalid_argument("UBMatrixHelper::GetMatrixFromPartialEigenDecomposition - input eigenvalues is not a column vector");

    const auto nEigenvalues = eigenvalues.GetRows();

    if (eigenvectorMatrix.GetColumns() != nEigenvalues)
        throw std::invalid_argument("UBMatrixHelper::GetMatrixFromPartialEigenDecomposition - the number of columns of the input eigenvector matrix, doesn't match the number of supplied eigenvalues");

    const auto size = eigenvectorMatrix.GetRows();
    if (size < nEigenvalues)
        throw std::invalid_argument("UBMatrixHelper::GetMatrixFromPartialEigenDecomposition - the number of rows of the input eigenvector matrix is smaller than the number of columns");

    // Get the diagonal matrix with eigenvalues along the diagonal
    const auto eigenvalueMatrix = ubsmear::UBMatrixHelper::GetDiagonalMatrix(eigenvalues);

    // Perform the transformation
    // ATTN here we are assuming that the output matrix is real and symmetric, and hence the inverse of the matrix of eigenvectors is given
    // by it's transpose. This is not true for a general matrix.
    return (eigenvectorMatrix * eigenvalueMatrix * eigenvectorMatrix.GetTranspose());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline std::pair<size_t, size_t> UBMatrixHelper::GetLargestOffDiagonalElement(const UBMatrix &matrix)
{
    float maxElementValue = -std::numeric_limits<float>::max();
    std::pair<size_t, size_t> maxElementIndex;
    bool foundElement = false;

    for (size_t iRow = 0; iRow < matrix.GetRows(); ++iRow)
    {
        for (size_t iCol = 0; iCol < matrix.GetColumns(); ++iCol)
        {
            // Skip diagonal entries
            if (iRow == iCol)
                continue;

            // Skip any values that aren't bigger then the current largest value found
            const auto value = std::abs(matrix.At(iRow, iCol));
            if (value < maxElementValue)
                continue;

            // Store this element
            foundElement = true;
            maxElementValue = value;
            maxElementIndex.first = iRow;
            maxElementIndex.second = iCol;
        }
    }

    if (!foundElement)
        throw std::logic_error("UBMatrixHelper::GetLargestOffDiagonalElement - Largest off-diagonal element doesn't exist");

    return maxElementIndex;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline float UBMatrixHelper::GetOffDiagonalSquaredSum(const UBMatrix &matrix)
{
    float sum = 0.f;

    for (size_t iRow = 0; iRow < matrix.GetRows(); ++iRow)
    {
        for (size_t iCol = 0; iCol < matrix.GetColumns(); ++iCol)
        {
            // Skip diagonal entries
            if (iRow == iCol)
                continue;

            sum += std::pow(matrix.At(iRow, iCol), 2.f);
        }
    }

    return sum;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline bool UBMatrixHelper::IsSquare(const UBMatrix &matrix)
{
    return (matrix.GetRows() == matrix.GetColumns());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline bool UBMatrixHelper::IsSymmetric(const UBMatrix &matrix)
{
    // A symmetric matrix is also square
    if (!UBMatrixHelper::IsSquare(matrix))
        return false;

    for (size_t iRow = 0; iRow < matrix.GetRows(); ++iRow)
    {
        for (size_t iCol = 0; iCol < iRow; ++iCol)
        {
            // Check if the matrix elements are identical under the exchange of rows and columns
            const auto diff = std::abs(matrix.At(iRow, iCol) - matrix.At(iCol, iRow));
            if (diff >= std::numeric_limits<float>::epsilon())
                return false;
        }
    }

    return true;
}

} // namespace ubsmear

#endif
