#ifndef UBSMEAR_UBSMEAR_HELPERS_UBFILEHELPER
#define UBSMEAR_UBSMEAR_HELPERS_UBFILEHELPER

#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>

namespace ubsmear
{

/*
* @brief A helper class for reading csv files
*/
class UBFileHelper
{
    public:
        /**
        * @brief Read a matrix from the supplied whitespace delimited file of floats
        *
        * @param fileName the input file name
        *
        * @return the matrix
        */
        static UBMatrix ReadMatrix(const std::string &fileName);

        /**
        * @brief Read a column vector from the suplied whitespace delimited file of float. The file should contain a single line.
        *
        * @param fileName the input file name
        *
        * @return the matrix
        */
        static UBMatrix ReadColumnVector(const std::string &fileName);

    private:

        /**
        * @brief Read a matrix from the supplied whitespace delimited file of floats
        *
        * @param fileName the input file name
        * @param data the output matrix data
        */
        static void ReadMatrixData(const std::string &fileName, std::vector< std::vector<float> > &data);
};

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBMatrix UBFileHelper::ReadMatrix(const std::string &fileName)
{
    // Read the matrix data from the file
    std::vector< std::vector<float> > data;
    UBFileHelper::ReadMatrixData(fileName, data);

    // Make a matrix from the data
    try
    {
        return UBMatrix(data);
    }
    catch (const std::exception &error)
    {
        throw std::runtime_error("UBFileHelper::ReadMatrix - Problem constructing matrix from file \"" + fileName + "\". " +
            "See error below: \n\n" + error.what());
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline UBMatrix UBFileHelper::ReadColumnVector(const std::string &fileName)
{
    // Read the matrix data from the file
    std::vector< std::vector<float> > data;
    UBFileHelper::ReadMatrixData(fileName, data);

    // Ensure we only have one line
    if (data.size() != 1)
    {
        throw std::runtime_error("UBFileHelper::ReadColumnVector - Problem constructing matrix from file \"" + fileName + "\". Found " +
            std::to_string(data.size()) + " lines of data, but expected 1");
    }

    // Make a matrix from the data
    try
    {
        // Make an Nx1 matrix from the data (i.e. a column vector)
        return UBMatrix(data.front(), data.front().size(), 1);
    }
    catch (const std::exception &error)
    {
        throw std::runtime_error("UBFileHelper::ReadColumnVector - Problem constructing matrix from file \"" + fileName + "\". " +
            "See error below: \n\n" + error.what());
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

inline void UBFileHelper::ReadMatrixData(const std::string &fileName, std::vector< std::vector<float> > &data)
{
    if (!data.empty())
        throw std::invalid_argument("UBFileHelper::ReadMatrixData - Vector supplied to accept matrix data isn't empty");

    // Open the input file
    std::ifstream file(fileName);
    if (!file.is_open())
        throw std::runtime_error("UBFileHelper::ReadMatrixData - Failed to open file \"" + fileName + "\"");

    // Read the file line-by-line
    std::string line;
    while (std::getline(file, line))
    {
        // Add a new row to the values
        data.emplace_back();
        auto &row = data.back();

        // Extract each value on the line
        float value;
        std::stringstream stream(line);
        while (!stream.eof())
        {
            stream >> value;

            if (stream.fail() || stream.bad())
            {
                throw std::runtime_error("UBFileHelper::ReadMatrixData - Problem reading line " + std::to_string(data.size()) + " of \"" +
                    fileName + "\". It must contain only floating point numbers and whitespace.");
            }

            row.push_back(value);
        }
    }
}

} // namespace ubsmear

#endif
