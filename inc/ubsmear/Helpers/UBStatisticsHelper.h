#ifndef UBSMEAR_UBSMEAR_HELPERS_UBSTATISTICSHELPER
#define UBSMEAR_UBSMEAR_HELPERS_UBSTATISTICSHELPER

#include "ubsmear/Objects/UBMatrix.h"

#include "ubsmear/Helpers/UBMatrixHelper.h"

#include <stdexcept>

namespace ubsmear
{

/*
* @brief A helper class for statistical analysis of cross-section data
*/
class UBStatisticsHelper
{
    public:

        /**
        * @brief Get the chi2 test statistic when comparing a predicted (smeared) cross-section to data.
        *
        * @param smearedPrediction the input smeared prediction for the cross-section, supplied as a column vector
        * @param data the input data cross-section (in reco space), supplied as a column vector
        * @param covarianceMatrix the total covariance matrix, with contributions from the data & smeared prediction
        * @param precision the precision to use when finding the inverse of the covariance matrix, smaller values are more precise
        *
        * @return a pair, first is the chi2, second is the number of degrees of freedom
        */
        static std::pair<float, size_t> GetChi2(const UBMatrix &smearedPrediction, const UBMatrix &data, const UBMatrix &covarianceMatrix, const float precision);
};

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

inline std::pair<float, size_t> UBStatisticsHelper::GetChi2(const UBMatrix &smearedPrediction, const UBMatrix &data, const UBMatrix &covarianceMatrix, const float precision)
{
    // Check the input objects have consistent dimensions
    if (smearedPrediction.GetColumns() != 1u)
        throw std::invalid_argument("UBStatisticsHelper::GetChi2 - input smeared prediction is not a column vector");

    const auto nBins = smearedPrediction.GetRows();

    if (data.GetColumns() != 1u)
        throw std::invalid_argument("UBStatisticsHelper::GetChi2 - input data is not a column vector");

    if (data.GetRows() != nBins)
        throw std::invalid_argument("UBStatisticsHelper::GetChi2 - input data has a different number of bins to the input smeared prediction");

    if (covarianceMatrix.GetRows() != nBins || covarianceMatrix.GetColumns() != nBins)
        throw std::invalid_argument("UBStatisticsHelper::GetChi2 - dimensions of the input covariance matrix, doesn't matrch the input data or smeared prediction");

    // Insist that the covariance matrix is symmetric
    if (!UBMatrixHelper::IsSymmetric(covarianceMatrix))
        throw std::invalid_argument("UBStatisticsHelper::GetChi2 - the input covariance matrix is not symmetric");

    // Get the eigenvalues and eigenvectors of the covariance matrix
    const auto &[eigenvalues, eigenvectorMatrix] = UBMatrixHelper::GetEigenDecomposition(covarianceMatrix, precision);

    // Get the matrix which transforms an input vector into the basis of eigenvectors
    const auto transformationMatrix = eigenvectorMatrix.GetTranspose();

    // Transform the smeared prediction and input data into the basis of eigenvectors
    const auto smearedPredictionEigen = transformationMatrix * smearedPrediction;
    const auto dataEigen = transformationMatrix * data;

    // Get the difference between the smeared prediction and the data in the eigenbasis
    const auto diffEigen = dataEigen - smearedPredictionEigen;

    // Calculate the test statistic
    float chi2 = 0.f;
    size_t dof = 0u;

    for (size_t iBin = 0; iBin < nBins; ++iBin)
    {
        // The eigenvalues of the covariance matrix give the variance along each eigenvector
        const auto variance = eigenvalues.At(iBin, 0);

        // Check if the eigenvalue is zero.
        // This represents a combination of bins in which there is no variance
        if (std::abs(variance) <= std::numeric_limits<float>::epsilon())
        {
            // ATTN here we are not allowing the input covariance matrix to be non-invertible
            // we could instead just skip this bin, and calculated the test-statistic in the subspace where the covariance matrix is invertible
            throw std::logic_error("UBStatisticsHelper::GetChi2 - input covariance matrix has a zero eigenvalue");
        }

        // Check that we don't have a negative eigenvalue - this should never happen for a true covariance matrix
        if (variance < 0.f)
            throw std::logic_error("UBStatisticsHelper::GetChi2 - input covariance matrix has negative eigenvalue");

        // Add to the total chi2 and increase the number of degrees of freedom
        const auto diff = diffEigen.At(iBin, 0);
        chi2 += (diff * diff) / variance;
        dof += 1u;
    }

    return {chi2, dof};
}

} // namespace ubsmear

#endif
