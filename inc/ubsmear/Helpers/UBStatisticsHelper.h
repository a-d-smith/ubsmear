#ifndef UBSMEAR_UBSMEAR_HELPERS_UBSTATISTICSHELPER
#define UBSMEAR_UBSMEAR_HELPERS_UBSTATISTICSHELPER

#include "ubsmear/Objects/UBMatrix.h"

#include "ubsmear/Helpers/UBMatrixHelper.h"

#include <stdexcept>
#include <cmath>

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

        /**
        * @brief Get the p-value from an input chi2 test statistic. The p-value is the probability of observing a test statistic at least as unlikely as the one supplied
        *
        * @param chi2 the input test statistic that follows a chi2 distribution
        * @param degreesOfFreedom the number of degrees of freedom of the chi2 distribution
        * @param precision the fractional error allowed on the returned p-value
        *
        * @return the p-value, or -float max if the p-value couldn't be found
        */
        static float GetPValue(const float chi2, const size_t degreesOfFreedom, const float precision);
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
        throw std::invalid_argument("UBStatisticsHelper::GetChi2 - dimensions of the input covariance matrix, doesn't match the input data or smeared prediction");

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
            throw std::logic_error("UBStatisticsHelper::GetChi2 - input covariance matrix has negative eigenvalue: " + std::to_string(variance));

        // Add to the total chi2 and increase the number of degrees of freedom
        const auto diff = diffEigen.At(iBin, 0);
        chi2 += (diff * diff) / variance;
        dof += 1u;
    }

    return {chi2, dof};
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline float UBStatisticsHelper::GetPValue(const float chi2, const size_t degreesOfFreedom, const float precision)
{
    if (precision < std::numeric_limits<float>::epsilon())
        throw std::invalid_argument("UBStatisticsHelper::GetPValue - input precision is <= 0, it should be 0 -> 1 not inclusive");

    if (precision >= 1.f)
        throw std::invalid_argument("UBStatisticsHelper::GetPValue - input precision is >= 1, it should be 0 -> 1 not inclusive");

    // The cumulative distribution function of the chi2 distribution with k degrees of freedom is given by:
    //
    // F(x; k) = g(k/2, x/2) / G(k/2)
    //  - g is the lower incomplete gamma function
    //  - G is the gamma function
    //
    // This gives the probability of a test statistic of x or less, to get the p-value we need 1 - F(x; k).

    // Without any extra dependencies, C++ includes the gamma function from cmath. However, the lower incomplete gamma function isn't
    // supplied. Here we use a power-series expansion of F(x; k):
    //
    // F(x; k) = (x/2)^(k/2) * exp(-x/2) * sum_r^inf [ (x/2)^r / G(k/2 + r + 1) ]
    //
    // This equation is simplified by the substitution y = x/2, s = k/2

    const auto y = chi2 / 2.f;
    const auto s = static_cast<float>(degreesOfFreedom) / 2.f;

    // Get the multiplictative factor that's outside of the sum
    const auto factor = std::pow(y, s) * std::exp(-y);

    // Define the function for the term that gets summed
    const auto getTerm = [&](const size_t r)
    {
        // ATTN gamma function has no zeros
        return (std::pow(y, static_cast<float>(r)) / std::tgamma(s + r + 1));
    };

    // Do the sum until the value is equal to the last iteration within the precision specified
    float sum = 0.f;
    float lastSum = 0.f;
    size_t r = 0;

    do
    {
        // If we have reached the maximum possible index, then there's no more that we can do
        if (r == std::numeric_limits<size_t>::max())
        {
            std::cerr << "UBStatisticsHelper::GetPValue - WARNING, Unable to calculate p-value, reached maximum number of terms without converging to desired precision" << std::endl;
            return -std::numeric_limits<float>::max();
        }

        lastSum = sum;
        sum += getTerm(r++);

        // Check we haven't reached the maximum possible number
        if (!std::isnormal(sum))
        {
            std::cerr << "UBStatisticsHelper::GetPValue - WARNING, Unable to calculate p-value, sum is not converging fast enough" << std::endl;
            return -std::numeric_limits<float>::max();
        }

    }
    while ( std::abs(sum - lastSum) > precision * sum );

    // Calculate the p-value
    const auto pValue = 1.f - (factor * sum);

    return pValue;
}

} // namespace ubsmear

#endif
