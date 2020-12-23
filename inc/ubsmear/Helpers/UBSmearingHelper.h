#ifndef UBSMEAR_UBSMEAR_HELPERS_UBSMEARINGHELPER
#define UBSMEAR_UBSMEAR_HELPERS_UBSMEARINGHELPER

#include "ubsmear/Objects/UBMatrix.h"
#include "ubsmear/Objects/UBXSecMeta.h"

#include "ubsmear/Helpers/UBMatrixHelper.h"

#include <stdexcept>
#include <random>

namespace ubsmear
{

/**
* @brief A helper class for smearing a predicted cross section so it can be compared to a data measurement
*/
class UBSmearingHelper
{
    public:

        /**
        * @brief Flatten an input (N x N) matrix into an (N^2 x 1) column vector.
        *        The original 2D matrix index [i, j] is flattend to a 1D index [N*i + j]
        *
        * @param matrix the input matrix
        *
        * @return the output column vector containing the flattened elements of the input matrix
        */
        static UBMatrix Flatten(const UBMatrix &matrix);

        /**
        * @brief Unflatten an input (N^2 x 1) column vector back to the (N x N) matrix from which it was producted
        *        The input 1D column vector index [i] is unflattened to the 2D index [i/N, i%N]
        *
        * @param flattenedMatrix the input flattened matrix, supplied as a column vector
        *
        * @return the output unflattened matrix
        */
        static UBMatrix Unflatten(const UBMatrix &flattenedMatrix);

        /**
        * @brief Get the vector of standard deviations from the vector of eigenvalues of a covariance matrix
        *
        * @param eigenvalues the input vector of eigenvalues
        *
        * @return the ouput vector of standard deviations
        */
        static UBMatrix GetStandardDeviationVector(const UBMatrix &eigenvalues);

        /**
        * @brief Generate a random vector drawn from a multivariate Gaussian distribution (defined by a covariance matrix and a mean vector)
        *
        * @param stdVector the column vector of standard deviations (i.e. square roots of the eigenvalues of the covariance matrix)
        * @param axisMatrix the matrix whose columns are the normalised eigenvectors of the covariance matrix
        * @param meanVector the mean vector
        *
        * @return a random vector draw from the multivariate Gaussian distribution supplied
        */
        static UBMatrix GetGaussianThrow(const UBMatrix &stdVector, const UBMatrix &axisMatrix, const UBMatrix &meanVector);

        /**
        * @brief Smear the input predicted cross-section by the smearing matrix
        *
        * @param metadata the input metadata which contains the information about the binning
        * @param prediction the input predicted differential cross-section in truth-space, supplied as a column vector
        * @param smearingMatrix the input smearing matrix
        *
        * @return the smeared prediction in reco space
        */
        static UBMatrix Smear(const UBXSecMeta &metadata, const UBMatrix &prediction, const UBMatrix &smearingMatrix);

        /**
        * @brief Smear the input predicted cross-section by the smearing matrix, and construct a covariance matrix for the result using the
        *        covariance matrix for the smearing matrix.
        *
        * @param metadata the input metadata which contains the information about the binning
        * @param prediction the input predicted differential cross-section in truth space, supplied as a column vector
        * @param predictionCovarianceMatrix the covariance matrix to encode an uncertainty on the prediction
        * @param smearingMatrix the input smearing matrix
        * @param smearingCovarianceMatrix the covariance matrix for the flattened smearing matrix
        * @param nUniverses the number of universes to use when propagating the uncertainy on the smearing matrix to the smeared prediction
        * @param precision the smaller the more precise the eigenvalues and eigenvectors of the smearing matrix will be, used to generate the universes
        *
        * @return a pair, first is the smeared prediction in reco space, second is the covariance matrix for the smeared prediction
        */
        static std::pair<UBMatrix, UBMatrix> Smear(const UBXSecMeta &metadata, const UBMatrix &prediction, const UBMatrix &predictionCovarianceMatrix, const UBMatrix &smearingMatrix, const UBMatrix &smearingCovarianceMatrix, const size_t nUniverses, const float precision);

        /**
        * @brief Trim underflow and overflow bins from the input matrix
        *
        * @param matrix the input matrix, can either be a square matrix or a column vector with the number of bins given in the metadata
        * @param metadata the input metadata, which contains information about the binning
        *
        * @return the trimmed matrix
        */
        static UBMatrix TrimUnderOverflowBins(const UBMatrix &matrix, const UBXSecMeta &metadata);

        /**
        * @brief Get a copy of the input vector with any negative values replace with zero
        *
        * @param vector the input vector
        *
        * @return the output positive semidefinite vector
        */
        static UBMatrix GetPositiveSemidefiniteVector(const UBMatrix &vector);

    private:

        static std::default_random_engine m_generator; ///< The random number generator
};

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

std::default_random_engine UBSmearingHelper::m_generator = std::default_random_engine();

// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBMatrix UBSmearingHelper::Flatten(const UBMatrix &matrix)
{
    // Check that the input matrix is square
    if (matrix.GetRows() != matrix.GetColumns())
        throw std::invalid_argument("UBSmearingHelper::Flatten - input matrix is not square");

    const auto size = matrix.GetRows();

    std::vector<float> flattenedMatrixElements;
    for (size_t iRow = 0; iRow < size; ++iRow)
    {
        for (size_t iCol = 0; iCol < size; ++iCol)
        {
            flattenedMatrixElements.push_back(matrix.At(iRow, iCol));
        }
    }

    return UBMatrix(flattenedMatrixElements, size*size, 1);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBMatrix UBSmearingHelper::Unflatten(const UBMatrix &flattenedMatrix)
{
    // Check that the input is a column vector
    if (flattenedMatrix.GetColumns() != 1)
        throw std::invalid_argument("UBSmearingHelper::Unflatten - the input flattened matrix is not a column vector");

    // The number of rows of the flattened matrix is the square of the size of the unflattened matrix
    const auto sizeSquared = flattenedMatrix.GetRows();

    // Get the integer closest to the square root of sizeSquared
    const auto size = static_cast<size_t>(std::round(std::pow(static_cast<float>(sizeSquared), 0.5f)));

    // Check that the input flattened matrix has a square number of entries
    if (size*size != sizeSquared)
    {
        throw std::invalid_argument("UBSmearingHelper::Unflatten - the input flattened matrix has " + std::to_string(sizeSquared)
                                    + " entries, this should be a square number");
    }

    // Read off the elements of the unflattened matrix
    std::vector< std::vector<float> > matrixElements;
    size_t flattenedIndex = 0;
    for (size_t iRow = 0; iRow < size; ++iRow)
    {
        // Add a new row
        matrixElements.emplace_back();
        auto &row = matrixElements.back();

        // Fill the entries of the row
        for (size_t iCol = 0; iCol < size; ++iCol)
            row.push_back(flattenedMatrix.At(flattenedIndex++, 0));
    }

    return UBMatrix(matrixElements);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBMatrix UBSmearingHelper::GetStandardDeviationVector(const UBMatrix &eigenvalues)
{
    // Check that the input object is a column vector
    if (eigenvalues.GetColumns() != 1u)
        throw std::invalid_argument("UBSmearingHelper::GetStandardDeviationVector - input eigenvalues is not a column vector");

    const auto size = eigenvalues.GetRows();

    std::vector<float> matrixElements;
    for (size_t iRow = 0; iRow < size; ++iRow)
    {
        const auto eigenvalue = eigenvalues.At(iRow, 0);

        // Handle eigenvalues of zero
        if (std::abs(eigenvalue) <= std::numeric_limits<float>::epsilon())
        {
            matrixElements.push_back(0.f);
            continue;
        }

        // Insist all eigenvalues are non-negative
        if (eigenvalue < 0.f)
            throw std::invalid_argument("UBSmearingHelper::GetStandardDeviationVector - found a negative input eigenvalue: " + std::to_string(eigenvalue));

        matrixElements.push_back(std::pow(eigenvalue, 0.5f));
    }

    return UBMatrix(matrixElements, size, 1);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBMatrix UBSmearingHelper::GetGaussianThrow(const UBMatrix &stdVector, const UBMatrix &axisMatrix, const UBMatrix &meanVector)
{
    // Check that the input matrices have consistent dimensions
    if (stdVector.GetColumns() != 1u)
        throw std::invalid_argument("UBSmearingHelper::GetGaussianThrow - input stdVector is not a column vector");

    const auto size = stdVector.GetRows();

    if (meanVector.GetColumns() != 1u)
        throw std::invalid_argument("UBSmearingHelper::GetGaussianThrow - input meanVector is not a column vector");

    if (meanVector.GetRows() != size)
        throw std::invalid_argument("UBSmearingHelper::GetGaussianThrow - input meanVector is a different size to the supplied vector of stdVector");

    if (axisMatrix.GetRows() != size)
        throw std::invalid_argument("UBSmearingHelper::GetGaussianThrow - the number of rows of the input eigenvector matrix, doesn't match the number of supplied stdVector");

    if (axisMatrix.GetColumns() != size)
        throw std::invalid_argument("UBSmearingHelper::GetGaussianThrow - the number of columns of the input eigenvector matrix, doesn't match the number of supplied stdVector");

    // Define a Gaussian distribution with mean zero, and standard deviation of one
    std::normal_distribution<float> gaussian(0.f, 1.f);

    // Working in the basis of eigenvectors of the covariance matrix, generate a vector of random Gaussian numbers using the stdVector
    std::vector<float> randomVectorElements;
    for (size_t i = 0; i < size; ++i)
        randomVectorElements.push_back(gaussian(m_generator) * stdVector.At(i, 0));

    UBMatrix randomVector(randomVectorElements, size, 1);

    // Use the matrix of eigenvectors to transform this random vector into the same basis as the input mean vector, and add on the mean
    return ( (axisMatrix * randomVector) + meanVector );
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBMatrix UBSmearingHelper::Smear(const UBXSecMeta &metadata, const UBMatrix &prediction, const UBMatrix &smearingMatrix)
{
    // Check that the input matrices have consistent dimensions
    if (prediction.GetColumns() != 1u)
        throw std::invalid_argument("UBSmearingHelper::Smear - input prediction is not a column vector");

    const auto size = metadata.GetNBins();

    if (prediction.GetRows() != size)
        throw std::invalid_argument("UBSmearingHelper::Smear - the number of rows of the input prediction doesn't match the number of bins in the metadata");

    if (smearingMatrix.GetRows() != size)
        throw std::invalid_argument("UBSmearingHelper::Smear - the number of rows of the input smearing matrix doesn't match the number of bins in the metadata");

    if (smearingMatrix.GetColumns() != size)
        throw std::invalid_argument("UBSmearingHelper::Smear - the number of columns of the input smearing matrix doesn't match the number of bins in the metadata");

    // Get the bin widths
    const auto binWidths = metadata.GetBinWidths();

    // Scale the prediction up by the bin widths (to get a quantity that's proportional to an event rate so it can be smeared)
    const auto scaledPrediction = ElementWiseOperation(prediction, binWidths, [](const auto &l, const auto& r) { return l * r; });

    // Smear the scaled prediction using the smearing matrix
    const auto smearedScaledPrediction = smearingMatrix * scaledPrediction;

    // Scale the smeared prediction back down by the bin widths
    return ElementWiseOperation(smearedScaledPrediction, binWidths, [](const auto &l, const auto& r)
    {
        if (r <= std::numeric_limits<float>::epsilon())
            throw std::invalid_argument("UBSmearingHelper::Smear - input metadata contains a bin with a zero or negative width");

        return l / r;
    });
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline std::pair<UBMatrix, UBMatrix> UBSmearingHelper::Smear(const UBXSecMeta &metadata, const UBMatrix &prediction, const UBMatrix &predictionCovarianceMatrix, const UBMatrix &smearingMatrix, const UBMatrix &smearingCovarianceMatrix, const size_t nUniverses, const float precision)
{
    // Check that the input matrices have consistent dimensions
    if (prediction.GetColumns() != 1u)
        throw std::invalid_argument("UBSmearingHelper::Smear - input prediction is not a column vector");

    const auto size = metadata.GetNBins();

    if (prediction.GetRows() != size)
        throw std::invalid_argument("UBSmearingHelper::Smear - the number of rows of the input prediction doesn't match the number of bins in the metadata");

    if (predictionCovarianceMatrix.GetRows() != size)
        throw std::invalid_argument("UBSmearingHelper::Smear - the number of rows of the input prediction covariance matrix doesn't match the number of bins in the metadata");

    if (predictionCovarianceMatrix.GetColumns() != size)
        throw std::invalid_argument("UBSmearingHelper::Smear - the number of columns of the input prediction covariance matrix doesn't match the number of bins in the metadata");

    if (smearingMatrix.GetRows() != size)
        throw std::invalid_argument("UBSmearingHelper::Smear - the number of rows of the input smearing matrix doesn't match the number of bins in the metadata");

    if (smearingMatrix.GetColumns() != size)
        throw std::invalid_argument("UBSmearingHelper::Smear - the number of columns of the input smearing matrix doesn't match the number of bins in the metadata");

    if (smearingCovarianceMatrix.GetRows() != size*size)
        throw std::invalid_argument("UBSmearingHelper::Smear - the number of rows of the input smearing covariance matrix doesn't match the number of bins in the metadata");

    if (smearingCovarianceMatrix.GetColumns() != size*size)
        throw std::invalid_argument("UBSmearingHelper::Smear - the number of columns of the input smearing covariance matrix doesn't match the number of bins in the metadata");

    // Check we have a valid number of universes
    if (nUniverses == 0)
        throw std::invalid_argument("UBSmearingHelper::Smear - the number of universes can't be zero");

    // Get the eigenvalues and eigenvectors of the prediction covariance matrix
    const auto [eigenvaluesPrediction, eigenvectorMatrixPrediction] = UBMatrixHelper::GetEigenDecomposition(predictionCovarianceMatrix, precision);

    // Get the eigenvalues and eigenvectors of the smearing covariance matrix
    const auto [eigenvaluesSmearing, eigenvectorMatrixSmearing] = UBMatrixHelper::GetEigenDecomposition(smearingCovarianceMatrix, precision);

    // ATTN Because the eigenvalues we find are not exact, it's possible that some may be found to be a very small negative number when the
    // true eigenvalue is zero. Here we manually clamp any negative values to be zero so that the eigenvalues are definitely non-negative
    const auto eigenvaluesPredictionClamped = UBSmearingHelper::GetPositiveSemidefiniteVector(eigenvaluesPrediction);
    const auto eigenvaluesSmearingClamped = UBSmearingHelper::GetPositiveSemidefiniteVector(eigenvaluesSmearing);

    // ATTN the eigenvalues of the covariance matrices give the variance of the distribution along the eigenvectors
    // Here we take their square root to obtain the standard devation
    const auto stdVectorPrediction = UBSmearingHelper::GetStandardDeviationVector(eigenvaluesPredictionClamped);
    const auto stdVectorSmearing = UBSmearingHelper::GetStandardDeviationVector(eigenvaluesSmearingClamped);

    // Get the flattened smearing matrix
    const auto flattenedSmearingMatrix = UBSmearingHelper::Flatten(smearingMatrix);

    // Get the smeared prediction
    const auto smearedPrediction = UBSmearingHelper::Smear(metadata, prediction, smearingMatrix);

    // Setup an empty matrix covariance matrix for the smearde prediction
    auto smearedPredictionCovarianceMatrix = UBMatrixHelper::GetZeroMatrix(size, size);

    // Compute the normalistaion factor for speed
    const auto norm = 1.f / static_cast<float>(nUniverses);

    // Repeat for each universe
    for (size_t iUni = 0; iUni < nUniverses; ++iUni)
    {
        // Generate a variation of the prediction according to its covariance matrix
        const auto predictionUni = UBSmearingHelper::GetGaussianThrow(stdVectorPrediction, eigenvectorMatrixPrediction, prediction);

        // Generate a variation of the smearing matrix according to its covariance matrix
        const auto smearingMatrixUni = UBSmearingHelper::Unflatten(UBSmearingHelper::GetGaussianThrow(stdVectorSmearing, eigenvectorMatrixSmearing, flattenedSmearingMatrix));

        // Get the smeared prediction in this universe
        const auto smearedPredictionUni = UBSmearingHelper::Smear(metadata, predictionUni, smearingMatrixUni);

        // Get the difference between the smeared prediction in this universe and in the central-value universe
        const auto diff = smearedPredictionUni - smearedPrediction;

        // Populate the prediction covariance matrix with this universe
        for (unsigned int iRow = 0; iRow < size; ++iRow)
        {
            for (unsigned int iCol = 0; iCol <= iRow; ++iCol)
            {
                // Get the contribution of this universe to the output covariance matrix
                const auto contribution = diff.At(iRow, 0) * diff.At(iCol, 0) * norm;
                const auto currentValue = smearedPredictionCovarianceMatrix.At(iRow, iCol);
                const auto newValue = currentValue + contribution;

                // ATTN The covariance matrix is symmetric by construction
                smearedPredictionCovarianceMatrix.SetElement(iRow, iCol, newValue);
                smearedPredictionCovarianceMatrix.SetElement(iCol, iRow, newValue);
            }
        }
    }

    return {smearedPrediction, smearedPredictionCovarianceMatrix};
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBMatrix UBSmearingHelper::TrimUnderOverflowBins(const UBMatrix &matrix, const UBXSecMeta &metadata)
{
    const auto nBins = metadata.GetNBins();

    // Check the input matrix has the correct dimensions
    if (matrix.GetRows() != nBins)
        throw std::invalid_argument("UBSmearingHelper::TrimUnderOverflowBins - the number of rows of the matrix doesn't match the supplied number of bins");

    // Check if the input matrix is a column vector, if not then insist that it is square
    const auto isColumnVector = (matrix.GetColumns() == 1u);

    if (!isColumnVector && matrix.GetColumns() != nBins)
        throw std::invalid_argument("UBSmearingHelper::TrimUnderOverflowBins - the number of columns of the matrix doesn't match the supplied number of bins");

    // If the metadata doesn't prescribe an underflow or overflow bin, then no trimming is required
    if (!metadata.HasUnderflow() && !metadata.HasOverflow())
        return matrix;

    // Otherwise, extract the relevant elements from the input matrix
    std::vector< std::vector<float> > elements;
    for (size_t iRow = 0; iRow < nBins; ++iRow)
    {
        // Skip underflow and overflow rows
        if (metadata.IsUnderOverflowBin(iRow))
            continue;

        // Add a new row to the output matrix elements
        elements.emplace_back();
        auto &row = elements.back();

        // Handle case where input matrix is a column vector
        if (isColumnVector)
        {
            row.emplace_back(matrix.At(iRow, 0));
            continue;
        }

        // Handle case where input matrix is square
        for (size_t iCol = 0; iCol < nBins; ++iCol)
        {
            // Skip underflow and overflow columns
            if (metadata.IsUnderOverflowBin(iCol))
                continue;

            row.emplace_back(matrix.At(iRow, iCol));
        }
    }

    return UBMatrix(elements);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBMatrix UBSmearingHelper::GetPositiveSemidefiniteVector(const UBMatrix &vector)
{
    if (vector.GetColumns() != 1u)
        throw std::invalid_argument("UBSmearingHelper::GetPositiveSemidefiniteVector input object isn't a column vector");

    const auto size = vector.GetRows();

    std::vector<float> elements;
    for (size_t iRow = 0; iRow < size; ++iRow)
    {
        const auto element = vector.At(iRow, 0);

        // Only use positive elements, set all other elements to zero
        elements.push_back(element > 0.f ? element : 0.f);
    }

    return UBMatrix(elements, size, 1);
}

} // namespace ubsmear

#endif
