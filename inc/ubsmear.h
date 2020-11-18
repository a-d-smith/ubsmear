#ifndef UBSMEAR_UBSMEAR
#define UBSMEAR_UBSMEAR

#include "ubsmear/Objects/UBMatrix.h"
#include "ubsmear/Objects/UBXSecMeta.h"

#include "ubsmear/Helpers/UBFileHelper.h"
#include "ubsmear/Helpers/UBMatrixHelper.h"
#include "ubsmear/Helpers/UBSmearingHelper.h"
#include "ubsmear/Helpers/UBStatisticsHelper.h"

/**
* @brief The ubsmear namespace holds the functionality required to compare a cross-section prediction to MicroBooNE data.
*
*        The functions exposed at the top level of this namespace should provide all the functionality required to use MicroBooNE
*        forward-folded differential cross-section data. The user should provide the input data as matrices using the UBMatrix class, and
*        information about the binning using the UBXSecMeta class (both of which are defined in ubsmear/Objects). The matrices can either be
*        constructed from an std::vector of their element, or read in directly from a txt file. These top-level functions make use of
*        lower-level machinary that's defined in ubsmear/Helpers.
*/
namespace ubsmear
{

/**
* @brief Apply the forward-folding methodology to smear an input differential cross-section prediction (defined in truth space) to the
*        corresponding prediction in reco space so that it can be compared directly to data.
*
* @param metadata an object containing information about the number of bins, and if an underflow or overflow bin is used
* @param prediction the input predicted differential cross-section in truth space
* @param predictionCovarianceMatrix the covariance matrix encoding any uncertainties on the input prediction
* @param smearingMatrix the input smearing matrix that defines the transformation from truth to reco space
* @param smearingCovarianceMatrix the covariance matrix encoding any uncertainties on the smearing matrix
* @param nUniverses the number of universes to use when propagating the smearing matrix uncertainties to the smeared prediction
* @param precision the precision to use when calculating eigenvalues and eigenvectors of a covariance matrix. For more see UBMatrixHelper::GetEigenDecomposition()
*
* @return a pair of matrices, the first is the smeared prediction (in reco space), the second is the covariance matrix encoding the total error on the smeared prediction
*/
std::pair<UBMatrix, UBMatrix> ForwardFold(const UBXSecMeta & metadata,
                                          const UBMatrix   & prediction,
                                          const UBMatrix   & predictionCovarianceMatrix,
                                          const UBMatrix   & smearingMatrix,
                                          const UBMatrix   & smearingCovarianceMatrix,
                                          const size_t       nUniverses,
                                          const float        precision);

/**
* @brief Compare a forward-folded differential cross-section prediction to measured data, and calculate a test statistic that is distributed
*        according to a chi2 distribution
*
* @param smearedPrediction the input smeared prediction (obtained using ubsmear::ForwardFold)
* @param smearedPredictionCovariance the covariance matrix for the smeared prediction (obtained using ubsmear::ForwardFold)
* @param data the input cross-section data (in reco space)
* @param dataCovariance the input covariance matrix encoding the uncertainties on the data
* @param precision the precision to use when calculating eigenvalues and eigenvectors of a covariance matrix. For more see UBMatrixHelper::GetEigenDecomposition()
*
* @return a test statistic that is distributed according to a chi2, supplied as a pair: first is the statistic, second if the number of degrees of freedom of the chi2 distribution
*/
std::pair<float, size_t> GetChi2(const UBMatrix & smearedPrediction,
                                 const UBMatrix & smearedPredictionCovariance,
                                 const UBMatrix & data,
                                 const UBMatrix & dataCovariance,
                                 const float      precision);

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

std::pair<UBMatrix, UBMatrix> ForwardFold(const UBXSecMeta &metadata, const UBMatrix &prediction, const UBMatrix &predictionCovarianceMatrix, const UBMatrix &smearingMatrix, const UBMatrix &smearingCovarianceMatrix, const size_t nUniverses, const float precision)
{
    // Apply the smearing matrix on the input predicted cross-section - include any and all underflow/overflow bins while smearing
    const auto &[smearedPrediction, smearedPredictionCovariance] = UBSmearingHelper::Smear(prediction, predictionCovarianceMatrix, smearingMatrix, smearingCovarianceMatrix, nUniverses, precision);

    // Trim any underflow or overflow bins, so the smeared prediction is on the same footing as the data
    return {
        UBSmearingHelper::TrimUnderOverflowBins(smearedPrediction, metadata),
        UBSmearingHelper::TrimUnderOverflowBins(smearedPredictionCovariance, metadata)
    };
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::pair<float, size_t> GetChi2(const UBMatrix &smearedPrediction, const UBMatrix &smearedPredictionCovariance, const UBMatrix &data, const UBMatrix &dataCovariance, const float precision)
{
    // Get the chi2 test statistic using the sum of the prediction & data covariance matrices
    return UBStatisticsHelper::GetChi2(smearedPrediction, data, (smearedPredictionCovariance + dataCovariance), precision);
}

} // namespace ubsmear

#endif
