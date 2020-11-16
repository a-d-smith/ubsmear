#include "ubsmear.h"

#include <iostream>

using namespace ubsmear;

int main()
{
    // Define some dummy data
    // The numbers here are really just make up to test the functionality of the code - but they are chosen to be semi-realistic

    // Define the number of bins
    const auto nBins = 3u;

    // Define a smearing matrix
    // Here we multiply two matricies, one whose columns normalise to one which describes the smearing of selected signal events, and a
    // second diagonal matirx which describes the selection efficiency in each true bin.
    const auto smearingMatrix = (

        // Smearing of selected signal events
        UBMatrix({
            {7, 1, 1},
            {3, 6, 1},
            {1, 2, 4}
        }).GetColumnNormalized() *

        // Selection efficiency
        UBMatrixHelper::GetDiagonalMatrix(UBMatrix( {0.6, 0.5, 0.4}, nBins, 1))
    );

    std::cout << "Smearing matrix" << std::endl;
    smearingMatrix.Print();

    // Define the covariance matrix for the flattened smearing matrix
    // To this, let's first define a fractional covariance matrix, by some supplying partial eigenvalues and eigenvectors
    const auto fractionalSmearingCovarianceMatrix = UBMatrixHelper::GetMatrixFromPartialEigenDecomposition(

        // The eigenvalues of the fractional smearing-covariance matrix
        UBMatrix( {0.2, 0.15, 0.1 }, 3, 1 ),

        // The eigenvectors of
        UBMatrix({
            {12, 1, 0,   3, 0,  0,   2, 0, 0 },
            {0,  1, 0,   1, 11, 0,   0, 1, 0 },
            {0,  0, 2,   0, 0,  3,   1, 0, 10}
        }).GetRowNormalized().GetTranspose()
    );

    // Now build the full covariance matrix from the fractional one
    const auto flattenedSmearingMatrix = UBSmearingHelper::Flatten(smearingMatrix);
    std::vector<float> smearingCovarianceMatrixElements;
    for (size_t iRow = 0; iRow < nBins*nBins; ++iRow)
    {
        for (size_t iCol = 0; iCol < nBins*nBins; ++iCol)
        {
            const auto entry = flattenedSmearingMatrix.At(iRow, 0) * flattenedSmearingMatrix.At(iCol, 0) * fractionalSmearingCovarianceMatrix.At(iRow, iCol);
            smearingCovarianceMatrixElements.push_back(entry);
        }
    }

    const UBMatrix smearingCovarianceMatrix(smearingCovarianceMatrixElements, nBins*nBins, nBins*nBins);

    std::cout << "Smearing covariance matrix" << std::endl;
    smearingCovarianceMatrix.Print();

    // Let's define the true underlying cross-section as a column vector
    const auto truth = UBMatrix( {1, 5, 2}, nBins, 1);
    std::cout << "True cross-section" << std::endl;
    truth.Print();

    // Now smear the truth to get our dummy "data"
    const auto data = UBSmearingHelper::Smear(truth, smearingMatrix);
    std::cout << "Measured data" << std::endl;
    data.Print();

    // Now make a prediction for the true underlying cross-section
    const auto prediction = UBMatrix( {1.3, 4.6, 1.5}, nBins, 1);
    std::cout << "Predicted cross-section" << std::endl;
    prediction.Print();

    // And put some uncertainty on the prediction - here we just use a diagonal matrix, with the 10% of the prediction itself as the variance
    const auto predictionCovarianceMatrix = UBMatrixHelper::GetDiagonalMatrix(0.1f * prediction);
    std::cout << "Predicted cross-section covariance matrix" << std::endl;
    predictionCovarianceMatrix.Print();

    // Smear the prediction to be able to compare it to the data and get it's covariance matrix
    const auto nUniverses = 100000u;
    const auto precision = 1e-7f;
    const auto &[smearedPrediction, smearedPredictionCovariance] = UBSmearingHelper::Smear(prediction, predictionCovarianceMatrix, smearingMatrix, smearingCovarianceMatrix, nUniverses, precision);

    std::cout << "Smeared prediction" << std::endl;
    smearedPrediction.Print();

    std::cout << "Smeared prediction covariance matrix" << std::endl;
    smearedPredictionCovariance.Print();


    return 0;
}
