#include "ubsmear.h"

#include <iostream>

using namespace ubsmear;

int main()
{
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Define the binning of the cross-section
    // -------------------------------------------------------------------------------------------------------------------------------------

    // In this case our dummy input data has 5 bins + an overflow
    const auto nBinsTotal = 6u;
    const auto hasUnderflow = false;
    const auto hasOverflow = true;

    UBXSecMeta metadata(nBinsTotal, hasUnderflow, hasOverflow);

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Load the dummy input data from the supplied files
    // -------------------------------------------------------------------------------------------------------------------------------------

    // Define the path directory containing our input data
    const std::string dataFilePath = "dummyData/";

    // Load the input data
    const auto data = UBFileHelper::ReadColumnVector(dataFilePath + "data.txt");
    const auto dataCovarianceMatrix = UBFileHelper::ReadMatrix(dataFilePath + "dataCovarianceMatrix.txt");

    const auto smearingMatrix = UBFileHelper::ReadMatrix(dataFilePath + "smearingMatrix.txt");
    const auto smearingCovarianceMatrix = UBFileHelper::ReadMatrix(dataFilePath + "smearingCovarianceMatrix.txt");

    const auto prediction = UBFileHelper::ReadColumnVector(dataFilePath + "prediction.txt");
    const auto predictionCovarianceMatrix = UBFileHelper::ReadMatrix(dataFilePath + "predictionCovarianceMatrix.txt");

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Print out the input data
    // -------------------------------------------------------------------------------------------------------------------------------------

    std::cout << "Cross-section data (reco-space):" << std::endl;
    data.Print();

    std::cout << "Cross-section data covariance matrix:" << std::endl;
    dataCovarianceMatrix.Print();

    std::cout << "Smearing matrix:" << std::endl;
    smearingMatrix.Print();

    std::cout << "Smearing covariance matrix:" << std::endl;
    smearingCovarianceMatrix.Print();

    std::cout << "Cross-section prediction (truth-space):" << std::endl;
    prediction.Print();

    std::cout << "Cross-section prediction covariance matrix:" << std::endl;
    predictionCovarianceMatrix.Print();

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Apply the smearing matrix to the predicted cross-section
    // -------------------------------------------------------------------------------------------------------------------------------------

    // Define the number of universe to use when propagating the uncertainties on the smearing matrix to the smeared prediction
    const auto nUniverses = 10000u;

    // Define the precision to use when calculating eigenvalues and eigenvectors of a covariance matrix
    const auto precision = std::numeric_limits<float>::epsilon();

    // Run the forward folding
    const auto &[smearedPrediction, smearedPredictionCovarianceMatrix] = ForwardFold(metadata,                                  // The binning information
                                                                                     prediction, predictionCovarianceMatrix,    // The input prediction
                                                                                     smearingMatrix, smearingCovarianceMatrix,  // The smearing matrix
                                                                                     nUniverses, precision);                    // Configuration options

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Print the smeared prediction
    // -------------------------------------------------------------------------------------------------------------------------------------

    std::cout << "Smeared prediction (reco-space):" << std::endl;
    smearedPrediction.Print();

    std::cout << "Smeared prediction covariance matrix" << std::endl;
    smearedPredictionCovarianceMatrix.Print();

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Calculate the chi2 test statistic between the smeared prediction and the data
    // -------------------------------------------------------------------------------------------------------------------------------------

    const auto &[chi2, degreesOfFreedom] = GetChi2(smearedPrediction, smearedPredictionCovarianceMatrix,  // The smeared (forward-folded) prediction
                                                   data, dataCovarianceMatrix,                            // The data
                                                   precision);                                            // Configuration options

    // For interest get the p-value from the chi2 - although better to just use a lookup table
    const auto pValue = UBStatisticsHelper::GetPValue(chi2, degreesOfFreedom, precision);

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Print the test statistic
    // -------------------------------------------------------------------------------------------------------------------------------------

    std::cout << "Chi2 / dof = " << chi2 << " / " << degreesOfFreedom << " = " << (chi2 / degreesOfFreedom) << std::endl;
    std::cout << "p-value = " << pValue << std::endl;

    return 0;
}
