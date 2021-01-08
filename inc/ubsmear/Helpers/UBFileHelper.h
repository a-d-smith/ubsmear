#ifndef UBSMEAR_UBSMEAR_HELPERS_UBFILEHELPER
#define UBSMEAR_UBSMEAR_HELPERS_UBFILEHELPER

#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cctype>

namespace ubsmear
{

/**
* @brief A helper class for reading input data files
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

        static std::string m_whitespace; ///< The characters to consider as whitespace
};

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

inline std::string UBFileHelper::m_whitespace = " \f\n\r\t\v"; // Space, form feed (f), line feed (n), carriage return (r), horizontal tab (t), vertical tab (v)

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
    size_t lineNumber = 0;
    while (std::getline(file, line))
    {
        lineNumber++;

        // Trim the line of any preceeding or trailing whitespace
        line.erase(0, line.find_first_not_of(m_whitespace));
        line.erase(line.find_last_not_of(m_whitespace) + 1);

        // Skip empty lines
        if (line.empty())
            continue;

        // Skip comment lines (starting with #)
        if (line.front() == '#')
            continue;

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
                // Clear any data that's currently been read
                data.clear();

                throw std::runtime_error("UBFileHelper::ReadMatrixData - Problem reading line " + std::to_string(lineNumber) + " of \"" +
                    fileName + "\". It must contain only floating point numbers and whitespace. "
                    + "Comments lines begin with # as the first character.");
            }

            row.push_back(value);
        }
    }
}

} // namespace ubsmear

#endif
