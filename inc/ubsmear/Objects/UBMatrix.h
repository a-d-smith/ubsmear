#ifndef UBSMEAR_UBSMEAR_OBJECTS_UBMATRIX
#define UBSMEAR_UBSMEAR_OBJECTS_UBMATRIX

#include <stdexcept>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

namespace ubsmear
{

/*
* @brief A simple matrix class
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

    private:

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
        throw std::invalid_argument("UBMatrix::At - the supplied row index, " + std::to_string(rowIndex) + " is out of bounds: 0 -> " + std::to_string(maxRowIndex));

    const auto maxColIndex = m_nCols - 1u;
    if (columnIndex > maxColIndex)
        throw std::invalid_argument("UBMatrix::At - the supplied column index, " + std::to_string(columnIndex) + " is out of bounds: 0 -> " + std::to_string(maxColIndex));

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

} // namespace ubsmear

#endif
