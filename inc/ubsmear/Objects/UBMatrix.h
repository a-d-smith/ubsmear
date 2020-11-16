#ifndef UBSMEAR_UBSMEAR_OBJECTS_UBMATRIX
#define UBSMEAR_UBSMEAR_OBJECTS_UBMATRIX

#include <stdexcept>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <functional>

namespace ubsmear
{

/*
* @brief A matrix class
*/
class UBMatrix
{
    public:

        /**
        * @brief Constructor using explicitly provided matrix elements, flattened into a vector
        *
        * @param matrixElements the flattened vector of matrix elements with index [row*nCols + col]
        * @param nRows the number of rows
        * @param nColumns the number of columns
        */
        UBMatrix(const std::vector<float> &matrixElements, const size_t nRows, const size_t nColumns);

        /**
        * @brief Constructor using explicitly provided matrix elements
        *
        * @param matrixElements the input matrix elements with indices [row][column]
        */
        UBMatrix(const std::vector< std::vector<float> > &matrixElements);

        /**
        * @brief Get the matrix element at the supplied row & column index
        *
        * @param rowIndex the desired rox index (counting from zero)
        * @param columnIndex the desired column index (counting from zero)
        *
        * @return the matrix element
        */
        float At(const size_t rowIndex, const size_t columnIndex) const;

        /**
        * @brief Get the number of rows
        *
        * @return the number of rows
        */
        size_t GetRows() const;

        /**
        * @brief Get the number of columns
        *
        * @return the number of columns
        */
        size_t GetColumns() const;

        /**
        * @brief Print the matrix to the screen
        */
        void Print() const;

        /**
        * @brief Get the transpose of this matrix (switching rows for columns)
        *
        * @return the transpose matrix
        */
        UBMatrix GetTranspose() const;

        /**
        * @brief Get a copy of this matrix where the rows are normalized to one
        *
        * @return the row normalized matrix
        */
        UBMatrix GetRowNormalized() const;

        /**
        * @brief Get a copy of this matrix where the columns are normalized to one
        *
        * @return the column-normalized matrix
        */
        UBMatrix GetColumnNormalized() const;

        /**
        * @brief Set the matrix element at the supplied row and column
        *
        * @param rowIndex the input row index
        * @param columnIndex the input column index
        * @param value the value to use
        */
        void SetElement(const size_t rowIndex, const size_t columnIndex, const float value);

        /**
        * @brief Overloaded + operator for matrix addition
        *
        * @param matrixL the first matrix to add
        * @param matrixR the second matrix to add
        *
        * @return the sum matrixL + matrixR
        */
        friend UBMatrix operator+(const UBMatrix &matrixL, const UBMatrix &matrixR);

        /**
        * @brief Overloaded - operator for matrix subtraction
        *
        * @param matrixL the first matrix to add
        * @param matrixR the second matrix to add
        *
        * @return the sum matrixL - matrixR
        */
        friend UBMatrix operator-(const UBMatrix &matrixL, const UBMatrix &matrixR);

        /**
        * @brief Overloaded * operator for matrix multiplication
        *
        * @param matrixL the left matrix
        * @param matrixR the right matrix
        *
        * @return the result of the matrix multiplcation matrixL * matrixR
        */
        friend UBMatrix operator*(const UBMatrix &matrixL, const UBMatrix &matrixR);

        /**
        * @brief Overloaded * operator for right multiplication by a scalar
        *
        * @param matrix the matrix
        * @param scalar the scalar
        *
        * @return the result of the multiplcation matrix * scalar
        */
        friend UBMatrix operator*(const UBMatrix &matrix, const float scalar);

        /**
        * @brief Overloaded * operator for left multiplication by a scalar
        *
        * @param scalar the scalar
        * @param matrix the matrix
        *
        * @return the result of the multiplcation scalar * matrix
        */
        friend UBMatrix operator*(const float scalar, const UBMatrix &matrix);

    private:

        /**
        * @brief Apply an element-wise operation between two matrices
        *
        * @param matrixL the left matrix
        * @param matrixR the right matrix
        * @param operation the element wise operation to apply
        *
        * @return the result of the element wise operation
        */
        friend UBMatrix ElementWiseOperation(const UBMatrix &matrixL, const UBMatrix &matrixR, const std::function<float(const float, const float)> &operation);


        size_t             m_nRows;    ///< The number of rows
        size_t             m_nCols;    ///< The number of columns
        std::vector<float> m_elements; ///< The matrix elements stored in a flattened vector
};

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBMatrix::UBMatrix(const std::vector<float> &matrixElements, const size_t nRows, const size_t nColumns) :
    m_nRows(nRows),
    m_nCols(nColumns),
    m_elements(matrixElements)
{
    // Validate the inputs
    if (nRows == 0)
        throw std::invalid_argument("UBMatrix::UBMatrix - Supplied number of rows can't be zero");

    if (nColumns == 0)
        throw std::invalid_argument("UBMatrix::UBMatrix - Supplied number of columns can't be zero");

    const auto nElements = matrixElements.size();
    if (nElements != nRows*nColumns)
    {
        throw std::invalid_argument("UBMatrix::UBMatrix - " + std::to_string(nElements) + " matrix elements were supplied. Expected " +
            std::to_string(nRows*nColumns) + " for a " + std::to_string(nRows) + "x" + std::to_string(nColumns) + " matrix");
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBMatrix::UBMatrix(const std::vector< std::vector<float> > &matrixElements)
{
    if (matrixElements.empty())
        throw std::invalid_argument("UBMatrix::UBMatrix - vector of matrix elements is empty");

    // Store the dimensions of the matrix
    m_nRows = matrixElements.size();
    m_nCols = matrixElements.front().size();

    // Store the matrix elements
    for (const auto &row : matrixElements)
    {
        // Check that this row has the same number of columns as the first
        if (row.size() != m_nCols)
            throw std::invalid_argument("UBMatrix::UBMatrix - the number of entries in each row is inconsistent");

        // Store the matrix elements in a flattened vector
        m_elements.insert(m_elements.end(), row.begin(), row.end());
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline float UBMatrix::At(const size_t rowIndex, const size_t columnIndex) const
{
    // Check that the requested index is not out of bounds
    const auto maxRowIndex = m_nRows - 1u;
    if (rowIndex > maxRowIndex)
        throw std::out_of_range("UBMatrix::At - the supplied row index, " + std::to_string(rowIndex) + " is out of bounds: 0 -> " + std::to_string(maxRowIndex));

    const auto maxColIndex = m_nCols - 1u;
    if (columnIndex > maxColIndex)
        throw std::out_of_range("UBMatrix::At - the supplied column index, " + std::to_string(columnIndex) + " is out of bounds: 0 -> " + std::to_string(maxColIndex));

    const auto flattenedIndex = (rowIndex * m_nCols) + columnIndex;
    return m_elements.at(flattenedIndex);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline size_t UBMatrix::GetRows() const
{
    return m_nRows;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline size_t UBMatrix::GetColumns() const
{
    return m_nCols;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline void UBMatrix::Print() const
{
    // Convert all of the matrix elements to strings, and keep track of the maximum width of each column
    std::vector<size_t> colWidths(m_nCols, 0);
    for (size_t iRow = 0; iRow < m_nRows; ++iRow)
    {
        for (size_t iCol = 0; iCol < m_nCols; ++iCol)
        {
            const auto valueString = std::to_string(this->At(iRow, iCol));
            colWidths.at(iCol) = std::max(colWidths.at(iCol), valueString.length());
        }
    }

    // Now actually print the matrix
    for (size_t iRow = 0; iRow < m_nRows; ++iRow)
    {
        std::cout << "(";

        for (size_t iCol = 0; iCol < m_nCols; ++iCol)
        {
            const auto valueString = std::to_string(this->At(iRow, iCol));
            std::cout << " " << std::setw(colWidths.at(iCol)) << valueString << " ";
        }

        std::cout << ")" << std::endl;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBMatrix UBMatrix::GetTranspose() const
{
    // Get the elements of the transpose matrix
    std::vector<float> elementsT;
    for (size_t i = 0; i < m_elements.size(); ++i)
    {
        // Get the row & column index corresponding to i for the *transpose* matrix
        const auto iRowT = i / m_nRows;
        const auto iColT = i % m_nRows;

        // Get the corresponding element switching rows for columns
        elementsT.push_back(this->At(iColT, iRowT));
    }

    // Construct the transpose matrix
    return UBMatrix(elementsT, m_nCols, m_nRows);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBMatrix UBMatrix::GetRowNormalized() const
{
    std::vector<float> elements;

    for (size_t iRow = 0; iRow < m_nRows; ++iRow)
    {
        // Get the sum of the elements in this row
        float rowSum = 0.f;
        for (size_t iCol = 0; iCol < m_nCols; ++iCol)
        {
            rowSum += this->At(iRow, iCol);
        }

        if (std::abs(rowSum) <= std::numeric_limits<float>::epsilon())
            throw std::logic_error("UBMatrix::GetRowNormalized - Found a row which sums to zero, can't normalize it");

        // Normalize the elements in the row
        const auto norm = 1.f / rowSum;
        for (size_t iCol = 0; iCol < m_nCols; ++iCol)
        {
            elements.push_back(this->At(iRow, iCol) * norm);
        }
    }

    return UBMatrix(elements, m_nRows, m_nCols);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBMatrix UBMatrix::GetColumnNormalized() const
{
    // Transpose the matrix, normalize it by row, then transpose it back
    return this->GetTranspose().GetRowNormalized().GetTranspose();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline void UBMatrix::SetElement(const size_t rowIndex, const size_t columnIndex, const float value)
{
    // Check that the requested index is not out of bounds
    const auto maxRowIndex = m_nRows - 1u;
    if (rowIndex > maxRowIndex)
        throw std::out_of_range("UBMatrix::SetElement - the supplied row index, " + std::to_string(rowIndex) + " is out of bounds: 0 -> " + std::to_string(maxRowIndex));

    const auto maxColIndex = m_nCols - 1u;
    if (columnIndex > maxColIndex)
        throw std::out_of_range("UBMatrix::SetElement - the supplied column index, " + std::to_string(columnIndex) + " is out of bounds: 0 -> " + std::to_string(maxColIndex));

    // Set the value of the matrix element at this index
    const auto flattenedIndex = (rowIndex * m_nCols) + columnIndex;
    m_elements.at(flattenedIndex) = value;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBMatrix ElementWiseOperation(const UBMatrix &matrixL, const UBMatrix &matrixR, const std::function<float(const float, const float)> &operation)
{
    // Check the input matrices have the same dimensions
    if (matrixL.m_nRows != matrixR.m_nRows)
        throw std::invalid_argument("UBMatrix::ElementWiseOperation - Input matrices have different numbers of rows");

    if (matrixL.m_nCols != matrixR.m_nCols)
        throw std::invalid_argument("UBMatrix::ElementWiseOperation - Input matrices have different numbers of columns");

    if (matrixL.m_elements.size() != matrixR.m_elements.size())
        throw std::invalid_argument("UBMatrix::ElementWiseOperation - Sanity check failed! Input matrices have same dimensions, but different number of elements.");

    const auto nElements = matrixL.m_elements.size();
    const auto nRows = matrixL.m_nRows;
    const auto nCols = matrixL.m_nCols;

    // Make a new vector to hold the elements
    std::vector<float> newElements(nElements);
    for (size_t i = 0; i < nElements; ++i)
    {
        // Do the element-wise operation
        newElements.at(i) = operation(matrixL.m_elements.at(i), matrixR.m_elements.at(i));
    }

    // Construct a new matrix with the added elements
    return UBMatrix(newElements, nRows, nCols);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBMatrix operator+(const UBMatrix &matrixL, const UBMatrix &matrixR)
{
    return ElementWiseOperation(matrixL, matrixR, [] (const auto &l, const auto &r) { return l + r; } );
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBMatrix operator-(const UBMatrix &matrixL, const UBMatrix &matrixR)
{
    return ElementWiseOperation(matrixL, matrixR, [] (const auto &l, const auto &r) { return l - r; } );
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBMatrix operator*(const UBMatrix &matrixL, const UBMatrix &matrixR)
{
    // Check the input matrices have compatible dimensions
    if (matrixL.m_nCols != matrixR.m_nRows)
    {
        throw std::invalid_argument("UBMatrix - Can't multiple matrices with dimensions: (" + std::to_string(matrixL.m_nRows) + ", " +
            std::to_string(matrixL.m_nCols) + ")*(" + std::to_string(matrixR.m_nRows) + ", " + std::to_string(matrixR.m_nCols) + ")");
    }

    const auto nRows = matrixL.m_nRows;
    const auto nCols = matrixR.m_nCols;
    const auto nDummy = matrixL.m_nCols;

    // Make a new vector to hold the elements
    std::vector< std::vector<float> > newElements;
    for (size_t iRow = 0; iRow < nRows; ++iRow)
    {
        std::vector<float> row(nCols, 0.f);

        for (size_t iCol = 0; iCol < nCols; ++iCol)
        {
            for (size_t iDummy = 0; iDummy < nDummy; ++iDummy)
            {
                row.at(iCol) += matrixL.At(iRow, iDummy) * matrixR.At(iDummy, iCol);
            }
        }

        newElements.push_back(row);
    }

    // Construct a new matrix with the added elements
    return UBMatrix(newElements);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBMatrix operator*(const UBMatrix &matrix, const float scalar)
{
    // Make a new vector to hold the elements
    std::vector<float> newElements;
    for (size_t i = 0; i < matrix.m_elements.size(); ++i)
    {
        newElements.push_back(matrix.m_elements.at(i) * scalar);
    }

    // Construct a new matrix with the added elements
    return UBMatrix(newElements, matrix.m_nRows, matrix.m_nCols);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBMatrix operator*(const float scalar, const UBMatrix &matrix)
{
    return matrix * scalar;
}

} // namespace ubsmear

#endif
